/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include <cstring>
#include <iomanip>
#include <sstream>

#include <ext.h>
#include <ext_buffer.h>
#include <ext_obex.h>
#include <ext_path.h>
#include <ext_systhread.h>
#include <z_dsp.h>

#include <OpenScofo.hpp>

static t_class *openscofo_class = nullptr;
#ifdef OPENSCOFO_LUA
int luaopen_max(lua_State *L);
#endif

struct MaxAction {
    double time;
    bool isLua;
    std::string Receiver;
    std::string LuaCode;
    t_atom *MaxArgs;
    int MaxArgsSize;
};

// ─────────────────────────────────────
class MaxOpenScofo {
  public:
    t_pxobject MaxObject;
    t_systhread_mutex Mutex;
    t_sample Sample;
    std::string PatchDir;

    // Clock
    t_clock *ClockEvent;
    t_clock *ClockActions;
    t_clock *ClockInfo;

    spdlog::level::level_enum log;

    // Actions
    std::vector<MaxAction> Actions;

    // Mir
    std::vector<OpenScofo::Descriptors> RequestMIR;
    bool MirOutput = false;

    // OpenScofo
    OpenScofo::OpenScofo *OpenScofo;
    std::unique_ptr<OpenScofo::Description> Desc;
    int Event;
    int StateIndex;
    float Tempo;
    bool Following;

    // Audio
    std::vector<double> inBuffer;
    int FFTSize;
    int HopSize;
    int BlockSize;
    int Sr;
    int BlockIndex;
    bool JustDescription;

    // Outlet
    void *EventOut;
    void *TempoOut;
    void *DescOut;
};

static void openscofo_error_callback(const spdlog::details::log_msg &log, void *data);
static void openscofo_assist(MaxOpenScofo *x, void *b, long m, long a, char *s);

static void openscofo_clear_actions(MaxOpenScofo *x) {
    if (x->ClockActions)
        clock_unset(x->ClockActions);
    if (x->ClockEvent)
        clock_unset(x->ClockEvent);
    if (x->ClockInfo)
        clock_unset(x->ClockInfo);

    for (MaxAction &action : x->Actions) {
        delete[] action.MaxArgs;
        action.MaxArgs = nullptr;
    }
    x->Actions.clear();
}

// ─────────────────────────────────────
static void openscofo_score(MaxOpenScofo *x, t_symbol *s) {
    if (!s || !s->s_name || !s->s_name[0]) {
        object_error((t_object *)x, "No score file provided");
        return;
    }

    char filename[MAX_PATH_CHARS];
    strncpy_zero(filename, s->s_name, MAX_PATH_CHARS);

    const char *ext = strrchr(filename, '.');
    if (!ext || strchr(ext, '/')) {
        snprintf(filename, MAX_PATH_CHARS, "%s.scofo", s->s_name);
    }

    short pathId = 0;
    t_fourcc fileType = 0;
    if (locatefile_extended(filename, &pathId, &fileType, nullptr, 0) != 0) {
        object_error((t_object *)x, "Can't find score file %s", s->s_name);
        return;
    }

    char scorePath[MAX_PATH_CHARS];
    path_toabsolutesystempath(pathId, filename, scorePath);

    systhread_mutex_lock(x->Mutex);
    x->Following = false;
    x->JustDescription = false;
    systhread_mutex_unlock(x->Mutex);
    openscofo_clear_actions(x);

    systhread_mutex_lock(x->Mutex);
    bool ok = x->OpenScofo->LoadScore(scorePath);
    if (ok) {
        x->OpenScofo->SetCurrentEvent(0);
        x->OpenScofo->ResetRuntimeState();
        x->FFTSize = x->OpenScofo->GetFFTSize();
        x->HopSize = x->OpenScofo->GetHopSize();
        x->inBuffer.assign(static_cast<size_t>(x->FFTSize), 0.0);
    }
    systhread_mutex_unlock(x->Mutex);

    if (ok) {
        object_post((t_object *)x, "Score loaded");
    } else {
        object_error((t_object *)x, "Score has errors");
        return;
    }

    x->Event = 0;
    x->StateIndex = 0;
    x->Tempo = x->OpenScofo->GetCurrentBPM();
    outlet_float(x->TempoOut, x->Tempo);
    outlet_float(x->EventOut, x->Event);

#ifdef OPENSCOFO_LUA
    std::string LuaCode = x->OpenScofo->GetLuaCode();
    bool result = x->OpenScofo->LuaExecute(LuaCode.c_str());

    if (!result) {
        std::string error = x->OpenScofo->LuaGetError();
        object_error((t_object *)x, "Lua error");
        object_error((t_object *)x, "%s", error.c_str());
    }
#endif

    systhread_mutex_lock(x->Mutex);
    x->Following = false;
    systhread_mutex_unlock(x->Mutex);
}

// ─────────────────────────────────────
static void openscofo_start(MaxOpenScofo *x) {
    if (!x->OpenScofo->ScoreIsLoaded()) {
        object_error((t_object *)x, "Score not loaded");
        return;
    }
    systhread_mutex_lock(x->Mutex);
    x->Following = false;
    x->JustDescription = false;
    systhread_mutex_unlock(x->Mutex);
    openscofo_clear_actions(x);

    systhread_mutex_lock(x->Mutex);
    x->OpenScofo->SetCurrentEvent(0);
    x->OpenScofo->ResetRuntimeState();
    x->Event = 0;
    x->StateIndex = 0;
    x->Tempo = x->OpenScofo->GetCurrentBPM();
    systhread_mutex_unlock(x->Mutex);

    outlet_float(x->TempoOut, x->Tempo);
    outlet_float(x->EventOut, x->Event);

    systhread_mutex_lock(x->Mutex);
    x->Following = true;
    systhread_mutex_unlock(x->Mutex);
    object_post((t_object *)x, "Start following");
}

// ─────────────────────────────────────
static void openscofo_output_descriptors(MaxOpenScofo *x, OpenScofo::Description &Desc) {
    for (auto it = x->RequestMIR.rbegin(); it != x->RequestMIR.rend(); ++it) {
        OpenScofo::Descriptors d = *it;
        switch (d) {
        case OpenScofo::Descriptors::INVALID:
            object_error((t_object *)x, "[openscofo~] Invalid descriptors");
            break;
        case OpenScofo::Descriptors::MFCC:
        case OpenScofo::Descriptors::CHROMA:
        case OpenScofo::Descriptors::LOGMEL:
        case OpenScofo::Descriptors::MAGNITUDE: {
            std::vector<double> DescArray = x->OpenScofo->GetDescriptionArray(Desc, d);
            int DescSize = DescArray.size();
            std::vector<t_atom> DescAtoms(DescSize);
            for (int i = 0; i < DescSize; ++i) {
                atom_setfloat(&DescAtoms[i], (t_float)DescArray[i]);
            }
            outlet_anything(x->DescOut, gensym(x->OpenScofo->GetDescriptionId(d)), DescSize, DescAtoms.data());
            break;
        }
        case OpenScofo::Descriptors::ONNX: {
            for (const auto &ONNXDesc : Desc.ONNX) {
                std::vector<t_atom> onnxAtoms(2);
                atom_setsym(&onnxAtoms[0], gensym(ONNXDesc.first.c_str()));
                atom_setfloat(&onnxAtoms[1], ONNXDesc.second);
                outlet_anything(x->DescOut, gensym("onnx"), 2, onnxAtoms.data());
            }
            break;
        }
        default:
            double DescValue = x->OpenScofo->GetDescriptionFloat(Desc, d);
            std::vector<t_atom> DescAtoms(1);
            atom_setfloat(&DescAtoms[0], (t_float)DescValue);
            outlet_anything(x->DescOut, gensym(x->OpenScofo->GetDescriptionId(d)), 1, DescAtoms.data());
        }
    }
}

// ─────────────────────────────────────
static std::vector<OpenScofo::Descriptors> openscofo_get_descriptors(MaxOpenScofo *x, int start, int argc, t_atom *argv) {
    std::vector<OpenScofo::Descriptors> Descriptors;
    Descriptors.reserve(argc - start);
    for (int i = start; i < argc; i++) {
        t_symbol *sym = atom_getsym(argv + i);
        Descriptors.push_back(x->OpenScofo->GetDescriptorsEnum(sym->s_name));
    }
    return Descriptors;
}

// ─────────────────────────────────────
static void openscofo_get(MaxOpenScofo *x, t_symbol *s, long argc, t_atom *argv) {
    (void)s;
    if (argc < 1) {
        object_error((t_object *)x, "Wrong number of arguments");
        return;
    }

    if (argv[0].a_type != A_SYM) {
        object_error((t_object *)x, "First argument of set method must be a symbol");
        return;
    }

    std::string method = atom_getsym(argv)->s_name;
    if (method == "description") {
        if (argc < 2 || argv[1].a_type != A_SYM) {
            object_error((t_object *)x, "description method requires <buffer~name>");
            return;
        }

        t_symbol *bufferName = atom_getsym(argv + 1);
        t_buffer_ref *bufferRef = buffer_ref_new((t_object *)x, bufferName);
        if (bufferRef == nullptr) {
            object_error((t_object *)x, "failed to create buffer reference");
            return;
        }

        t_buffer_obj *bufferObj = buffer_ref_getobject(bufferRef);
        if (bufferObj == nullptr) {
            object_error((t_object *)x, "buffer %s not found", bufferName->s_name);
            object_free(bufferRef);
            return;
        }

        long frameCount = buffer_getframecount(bufferObj);
        if (frameCount <= 0) {
            object_error((t_object *)x, "buffer %s is empty", bufferName->s_name);
            object_free(bufferRef);
            return;
        }

        long start = 0;
        if (argc == 2) {
            start = 0;
        } else if (argc == 3) {
            start = atom_getlong(argv + 2);
            if (start < 0 || start >= frameCount) {
                object_error((t_object *)x, "invalid start index %ld for buffer size %ld and fftsize %d", start,
                             frameCount, x->OpenScofo->GetFFTSize());
                object_free(bufferRef);
                return;
            }
        } else {
            object_error((t_object *)x, "Wrong arguments");
            object_free(bufferRef);
            return;
        }

        float *samples = buffer_locksamples(bufferObj);
        if (samples == nullptr) {
            object_error((t_object *)x, "failed to read buffer samples from %s", bufferName->s_name);
            object_free(bufferRef);
            return;
        }

        int fftsize = x->OpenScofo->GetFFTSize();
        std::vector<double> audioBuffer(fftsize, 0.0);
        for (int i = 0; i < fftsize; i++) {
            const long src = start + i;
            if (src >= frameCount) {
                break;
            }
            audioBuffer[i] = static_cast<double>(samples[src]);
        }

        buffer_unlocksamples(bufferObj);
        object_free(bufferRef);

        x->OpenScofo->ProcessBlock(audioBuffer.data(), audioBuffer.size());
        OpenScofo::Description Desc = x->OpenScofo->GetDescription();
        openscofo_output_descriptors(x, Desc);
    }
}

// ─────────────────────────────────────
static void openscofo_set(MaxOpenScofo *x, t_symbol *s, long argc, t_atom *argv) {
    (void)s;

    if (argc < 1) {
        object_error((t_object *)x, "Wrong number of arguments");
        return;
    }

    if (argv[0].a_type != A_SYM) {
        object_error((t_object *)x, "First argument of set method must be a symbol");
        return;
    }

    std::string method = atom_getsym(argv)->s_name;
    if (method == "event") {
        if (argc < 2) {
            object_error((t_object *)x, "Wrong number of arguments");
            return;
        }

        long f = atom_getlong(argv + 1);
        x->Event = (int)f;
        x->OpenScofo->SetCurrentEvent((int)f);
    } else if (method == "onnxmodel") {
        if (argc < 3 || argv[1].a_type != A_SYM) {
            object_error((t_object *)x, "onnxmodel requires: set onnxmodel <path> <descriptors...>");
            return;
        }

        std::vector<OpenScofo::Descriptors> desc = openscofo_get_descriptors(x, 2, argc, argv);

        if (desc.size() == 0) {
            object_error((t_object *)x, "ONNX models require the descriptor order used in training");
            return;
        }

        char filename[MAX_PATH_CHARS];
        strncpy_zero(filename, atom_getsym(argv + 1)->s_name, MAX_PATH_CHARS);

        short pathId = 0;
        t_fourcc fileType = 0;
        if (locatefile_extended(filename, &pathId, &fileType, nullptr, 0) != 0) {
            object_error((t_object *)x, "Can't find ONNX model %s", atom_getsym(argv + 1)->s_name);
            return;
        }

        char modelPath[MAX_PATH_CHARS];
        path_toabsolutesystempath(pathId, filename, modelPath);
        x->OpenScofo->LoadONNXModel(modelPath, desc);
    } else if (method == "verbosity") {
        if (argc < 2) {
            object_error((t_object *)x, "Wrong number of arguments");
            return;
        }

        long f = atom_getlong(argv + 1);
        switch (f) {
        case 0:
            x->log = spdlog::level::warn;
            break;
        case 1:
            x->log = spdlog::level::info;
            break;
        case 2:
            x->log = spdlog::level::debug;
            break;
        case 3:
            x->log = spdlog::level::trace;
            break;
        default:
            object_error((t_object *)x, "Invalid verbosity value %ld", f);
            return;
        }
        x->OpenScofo->SetLogLevel(x->log);
    } else if (method == "section") {
        if (argc != 2) {
            object_error((t_object *)x, "section requires one name");
            return;
        }

        std::string section;
        if (atom_gettype(argv + 1) == A_SYM) {
            section = atom_getsym(argv + 1)->s_name;
        } else if (atom_gettype(argv + 1) == A_LONG) {
            section = std::to_string(atom_getlong(argv + 1));
        } else if (atom_gettype(argv + 1) == A_FLOAT) {
            std::ostringstream stream;
            stream << std::setprecision(15) << atom_getfloat(argv + 1);
            section = stream.str();
        } else {
            object_error((t_object *)x, "section name must be a symbol or number");
            return;
        }

        if (x->OpenScofo->SetCurrentSection(section)) {
            x->StateIndex = x->OpenScofo->GetCurrentStateIndex();
            outlet_float(x->TempoOut, x->OpenScofo->GetCurrentBPM());
        }
    } else if (method == "justdescription" || method == "description") {
        if (argc < 2) {
            object_error((t_object *)x, "Wrong number of arguments");
            return;
        }
        long f = atom_getlong(argv + 1);
        x->JustDescription = f != 0;
    } else {
        object_error((t_object *)x, "Unknown method");
    }
}

// ─────────────────────────────────────
static void openscofo_section(MaxOpenScofo *x, t_symbol *s, long argc, t_atom *argv) {
    t_atom args[2];
    atom_setsym(args, gensym("section"));
    if (argc == 1) {
        args[1] = argv[0];
    }
    openscofo_set(x, s, argc + 1, args);
}

// ─────────────────────────────────────
static void openscofo_following(MaxOpenScofo *x, long f) {
    if (!x->OpenScofo->ScoreIsLoaded()) {
        object_error((t_object *)x, "Score not loaded");
        return;
    }
    if (f == 1) {
        x->Following = true;
    } else {
        x->Following = false;
    }
}

// ─────────────────────────────────────
static void openscofo_luaexecute(MaxOpenScofo *x, std::string code) {
#ifdef OPENSCOFO_LUA
    if (!x->OpenScofo->LuaExecute(code)) {
        std::string error = x->OpenScofo->LuaGetError();
        object_error((t_object *)x, "Lua error");
        object_error((t_object *)x, "%s", error.c_str());
    }
#endif
}

// ─────────────────────────────────────
static void openscofo_maxsend(MaxOpenScofo *x, std::string r, int argc, t_atom *argv) {
    t_object *receiver = (t_object *)gensym(r.c_str())->s_thing;
    if (!receiver) {
        object_error((t_object *)x, "Receiver %s not found", r.c_str());
        return;
    }

    if (argc == 0) {
        object_method_typed(receiver, gensym("bang"), 0, nullptr, nullptr);
    } else {
        object_method_typed(receiver, gensym("list"), argc, argv, nullptr);
    }
}

// ─────────────────────────────────────
static t_atom *openscofo_convertargs(OpenScofo::ScoreAction &action) {
    int size = action.Args.size();
    t_atom *maxArgs = new t_atom[size];

    for (int i = 0; i < size; i++) {
        std::variant<float, int, std::string> arg = action.Args[i];
        if (std::holds_alternative<float>(arg)) {
            atom_setfloat(&maxArgs[i], std::get<float>(arg));
        } else if (std::holds_alternative<int>(arg)) {
            atom_setfloat(&maxArgs[i], std::get<int>(arg));
        } else if (std::holds_alternative<std::string>(arg)) {
            atom_setsym(&maxArgs[i], gensym(std::get<std::string>(arg).c_str()));
        }
    }
    return maxArgs;
}

// ─────────────────────────────────────
static void openscofo_tickactions(void *xv) {
    MaxOpenScofo *x = (MaxOpenScofo *)xv;
    const double currentTime = gettime();
    const double nextBlock = 1000.0 / x->Sr * x->BlockSize;
    const double nextTime = currentTime + nextBlock;

    std::vector<MaxAction>::iterator it = x->Actions.begin();
    while (it != x->Actions.end()) {
        MaxAction &curAction = *it;
        if (curAction.time <= nextTime) {
            if (curAction.isLua) {
                openscofo_luaexecute(x, curAction.LuaCode);
            } else {
                openscofo_maxsend(x, curAction.Receiver, curAction.MaxArgsSize, curAction.MaxArgs);
            }
            delete[] curAction.MaxArgs;
            it = x->Actions.erase(it);
        } else {
            ++it;
        }
    }
}

// ─────────────────────────────────────
static void openscofo_tickinfo(MaxOpenScofo *xv) {
    MaxOpenScofo *x = (MaxOpenScofo *)xv;
    if (x->MirOutput) {
        OpenScofo::Description Desc = x->OpenScofo->GetDescription();
        openscofo_output_descriptors(x, Desc);
    }
}

// ─────────────────────────────────────
static void openscofo_ticknewevent(void *xv) {
    MaxOpenScofo *x = (MaxOpenScofo *)xv;

    int prevStateIndex = x->StateIndex;

    x->Event = x->OpenScofo->GetCurrentScorePosition();
    x->StateIndex = x->OpenScofo->GetCurrentStateIndex();

    OpenScofo::EventActions audioStateActions = x->OpenScofo->GetAudioStateChangeActions();
    for (OpenScofo::ScoreAction &action : audioStateActions) {
        t_atom *args = openscofo_convertargs(action);
        openscofo_maxsend(x, action.Receiver, action.Args.size(), args);
        delete[] args;
    }

    if (prevStateIndex == x->StateIndex) {
        return;
    }

    outlet_float(x->TempoOut, x->OpenScofo->GetCurrentBPM());
    outlet_int(x->EventOut, static_cast<int>(x->Event));

    OpenScofo::EventActions actions = x->OpenScofo->GetCurrentEventActions();

    for (OpenScofo::ScoreAction &act : actions) {
        double time = act.Time;

        if (!act.AbsoluteTime) {
            act.Time = 60.0 / x->OpenScofo->GetCurrentBPM() * act.Time * 1000;
            time = act.Time;
        }

        if (time == 0) {
            if (act.isLua) {
                openscofo_luaexecute(x, act.Lua);
            } else {
                t_atom *maxArgs = openscofo_convertargs(act);
                openscofo_maxsend(x, act.Receiver, act.Args.size(), maxArgs);
                delete[] maxArgs;
            }
        } else {
            // If Max has an equivalent to Pd's clock_getsystimeafter(), prefer it.
            // Otherwise retain gettime() but note the semantic difference.
            double sysTime = gettime() + time;

            int size = act.Args.size();
            std::string receiver = act.Receiver;
            t_atom *maxArgs = openscofo_convertargs(act);

            MaxAction action = {sysTime, act.isLua, receiver, act.Lua, maxArgs, size};
            x->Actions.push_back(action);
        }
    }
}

// ─────────────────────────────────────
static void openscofo_perform_score(MaxOpenScofo *x, double **ins, long sampleframes) {
    if (systhread_mutex_trylock(x->Mutex) != 0) {
        return;
    }

    if (!x->Following && !x->JustDescription) {
        systhread_mutex_unlock(x->Mutex);
        return;
    }

    x->BlockIndex = 0;
    bool ok = x->OpenScofo->ProcessBlock(ins[0], sampleframes);
    if (!ok) {
        systhread_mutex_unlock(x->Mutex);
        return;
    }

    clock_delay(x->ClockActions, 0);
    clock_delay(x->ClockEvent, 0);
    clock_delay(x->ClockInfo, 0);
    systhread_mutex_unlock(x->Mutex);
}

// ─────────────────────────────────────
static void openscofo_perform64(t_object *obj, t_object *dsp64, double **ins, long numins, double **outs, long numouts,
                             long sampleframes, long flags, void *userparam) {
    auto *x = (MaxOpenScofo *)obj;
    (void)dsp64;
    (void)numins;
    (void)outs;
    (void)numouts;
    (void)flags;
    (void)userparam;
    openscofo_perform_score(x, ins, sampleframes);
}

// ─────────────────────────────────────
static void openscofo_dsp64(MaxOpenScofo *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize,
                         long flags) {
    (void)count;
    (void)flags;
    x->BlockSize = maxvectorsize;
    x->BlockIndex = 0;
    x->Sr = (int)samplerate;

    systhread_mutex_lock(x->Mutex);
    OpenScofo::Configuration config = x->OpenScofo->GetConfiguration();
    if ((int)config.SR != x->Sr) {
        config.SR = (float)x->Sr;
        x->OpenScofo->SetConfiguration(config);
    }
    systhread_mutex_unlock(x->Mutex);

    x->inBuffer.resize((size_t)x->FFTSize, 0.0f);
    dsp_add64(dsp64, (t_object *)x, openscofo_perform64, 0, nullptr);
}

// ─────────────────────────────────────
static void openscofo_error_callback(const spdlog::details::log_msg &log, void *data) {
    MaxOpenScofo *x = static_cast<MaxOpenScofo *>(data);
    spdlog::level::level_enum maxlevel = x->log;
    if (log.level < maxlevel) {
        return;
    }

    std::string text(log.payload.data(), log.payload.size());
    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
        object_error((t_object *)x, "%s", text.c_str());
        break;
    case spdlog::level::info:
    case spdlog::level::warn:
        object_post((t_object *)x, "%s", text.c_str());
        break;
    case spdlog::level::debug:
    case spdlog::level::trace:
        object_post((t_object *)x, "%s", text.c_str());
        break;
    default:
        break;
    }
}

// ─────────────────────────────────────
static void openscofo_assist(MaxOpenScofo *x, void *b, long m, long a, char *s) {
    (void)x;
    (void)b;
    if (m == ASSIST_OUTLET) {
        switch (a) {
        case 0:
            snprintf(s, 256, "Score Event Index");
            break;
        case 1:
            snprintf(s, 256, "Tempo in BPM of the current performance");
            break;
        case 2:
            snprintf(s, 256, "Descriptor output");
            break;
        }
    } else {
        switch (a) {
        case 0:
            snprintf(s, 256, "Signal Input");
            break;
        }
    }
}

// ─────────────────────────────────────
static void *openscofo_new(t_symbol *s, long argc, t_atom *argv) {
    MaxOpenScofo *x = (MaxOpenScofo *)object_alloc(openscofo_class);
    (void)s;

    if (!x) {
        return nullptr;
    }

    x->Mutex = nullptr;
    systhread_mutex_new(&x->Mutex, 0);
    dsp_setup((t_pxobject *)x, 1);

    x->FFTSize = 2048;
    x->HopSize = 512;
    x->Sr = 48000;
    x->Following = false;
    x->Event = -1;
    x->StateIndex = 0;
    x->BlockIndex = 0;
    x->JustDescription = false;
    x->log = spdlog::level::warn;

    x->OpenScofo = new OpenScofo::OpenScofo((float)x->Sr, (float)x->FFTSize, (float)x->HopSize);

    x->RequestMIR = openscofo_get_descriptors(x, 0, argc, argv);
    if (x->RequestMIR.size() > 0) {
        x->DescOut = outlet_new(x, nullptr);
        x->MirOutput = true;
    } else {
        x->MirOutput = false;
    }

    x->TempoOut = outlet_new(x, "float");
    x->EventOut = outlet_new(x, "float");

    x->ClockEvent = clock_new(x, (method)openscofo_ticknewevent);
    x->ClockActions = clock_new(x, (method)openscofo_tickactions);
    x->ClockInfo = clock_new(x, (method)openscofo_tickinfo);

    char patchPath[MAX_PATH_CHARS];
    short pathId = path_getdefault();
    path_toabsolutesystempath(pathId, nullptr, patchPath);
    x->PatchDir = patchPath;

    x->OpenScofo->SetErrorCallback(openscofo_error_callback, static_cast<void *>(x));
    x->OpenScofo->SetLogLevel(x->log);
    x->OpenScofo->SetRequestedDescriptors(x->RequestMIR);

#ifdef OPENSCOFO_LUA
    x->OpenScofo->LuaAddModule("max", luaopen_max);
    x->OpenScofo->LuaAddPath(x->PatchDir);
    x->OpenScofo->LuaAddPointer(x, "_maxobj");
#endif

    return (void *)x;
}

// ─────────────────────────────────────
static void openscofo_free(MaxOpenScofo *x) {
    dsp_free((t_pxobject *)x);
    openscofo_clear_actions(x);

    if (x->ClockEvent) {
        object_free(x->ClockEvent);
    }
    if (x->ClockActions) {
        object_free(x->ClockActions);
    }
    if (x->ClockInfo) {
        object_free(x->ClockInfo);
    }

    systhread_mutex_lock(x->Mutex);
    delete x->OpenScofo;
    x->OpenScofo = nullptr;
    systhread_mutex_unlock(x->Mutex);
    systhread_mutex_free(x->Mutex);
}

// ─────────────────────────────────────
int C74_EXPORT main() {
    t_class *c =
        class_new("openscofo~", (method)openscofo_new, (method)openscofo_free, (long)sizeof(MaxOpenScofo), 0L, A_GIMME, 0);

    object_post(nullptr, "[openscofo~] version %s (%s), by Charles K. Neimog", OPENSCOFO_VERSION, OPENSCOFO_BUILD_TIME);

    class_addmethod(c, (method)openscofo_score, "score", A_SYM, 0);
    class_addmethod(c, (method)openscofo_start, "start", A_NOTHING, 0);
    class_addmethod(c, (method)openscofo_following, "follow", A_LONG, 0);
    class_addmethod(c, (method)openscofo_section, "section", A_GIMME, 0);
    class_addmethod(c, (method)openscofo_set, "set", A_GIMME, 0);
    class_addmethod(c, (method)openscofo_get, "get", A_GIMME, 0);

    class_addmethod(c, (method)stdinletinfo, "inletinfo", A_CANT, 0);
    class_addmethod(c, (method)openscofo_assist, "assist", A_CANT, 0);

    class_addmethod(c, (method)openscofo_dsp64, "dsp64", A_CANT, 0);

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    openscofo_class = c;
    return 0;
}
