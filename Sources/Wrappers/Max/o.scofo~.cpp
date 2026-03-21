#include <filesystem>
#include <cstring>

#include <ext.h>
#include <ext_buffer.h>
#include <ext_obex.h>
#include <z_dsp.h>

#include <OpenScofo.hpp>

static t_class *oscofo_class = nullptr;
#ifdef OSCOFO_LUA
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

static void oscofo_error_callback(const spdlog::details::log_msg &log, void *data);
static void oscofo_assist(MaxOpenScofo *x, void *b, long m, long a, char *s);

// ─────────────────────────────────────
static void oscofo_score(MaxOpenScofo *x, t_symbol *s) {
    if (!s) {
        object_error((t_object *)x, "No score file provided");
        return;
    }

    bool ok;
    std::string scorePath = s->s_name;
    if (!std::filesystem::exists(s->s_name)) {
        scorePath = x->PatchDir + "/" + s->s_name;
    }

    ok = x->OpenScofo->ParseScore(scorePath);
    if (ok) {
        object_post((t_object *)x, "Score loaded");
    } else {
        object_error((t_object *)x, "Score has errors");
        return;
    }
    x->OpenScofo->SetCurrentEvent(0);

    x->Event = 0;
    outlet_float(x->TempoOut, x->OpenScofo->GetLiveBPM());
    outlet_float(x->EventOut, 0);

    x->FFTSize = x->OpenScofo->GetFFTSize();
    x->HopSize = x->OpenScofo->GetHopSize();
    x->inBuffer.resize(x->FFTSize, 0.0f);

#ifdef OSCOFO_LUA
    std::string LuaCode = x->OpenScofo->GetLuaCode();
    bool result = x->OpenScofo->LuaExecute(LuaCode.c_str());

    if (!result) {
        std::string error = x->OpenScofo->LuaGetError();
        object_error((t_object *)x, "Lua error");
        object_error((t_object *)x, "%s", error.c_str());
    }
#endif
}

// ─────────────────────────────────────
static void oscofo_start(MaxOpenScofo *x) {
    if (!x->OpenScofo->ScoreIsLoaded()) {
        object_error((t_object *)x, "Score not loaded");
        return;
    }
    x->Actions.clear();
    x->OpenScofo->SetCurrentEvent(0);
    x->Event = 0;

    outlet_float(x->TempoOut, x->OpenScofo->GetLiveBPM());
    outlet_float(x->EventOut, x->Event);

    x->Following = true;
    object_post((t_object *)x, "Start following");
}

// ─────────────────────────────────────
static void oscofo_output_descriptors(MaxOpenScofo *x, OpenScofo::Description &Desc) {
    for (auto it = x->RequestMIR.rbegin(); it != x->RequestMIR.rend(); ++it) {
        OpenScofo::Descriptors d = *it;
        switch (d) {
        case OpenScofo::Descriptors::INVALID:
            object_error((t_object *)x, "[o.scofo~] Invalid descriptors");
            break;
        case OpenScofo::Descriptors::MFCC:
        case OpenScofo::Descriptors::CHROMA:
        case OpenScofo::Descriptors::MELOGRAM:
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
                SETSYMBOL(&onnxAtoms[0], gensym(ONNXDesc.first.c_str()));
                SETFLOAT(&onnxAtoms[1], ONNXDesc.second);
                outlet_anything(x->DescOut, gensym("onnx"), 2, onnxAtoms.data());
            }
            break;
        }
        default:
            double DescValue = x->OpenScofo->GetDescriptionFloat(Desc, d);
            std::vector<t_atom> DescAtoms(1);
            atom_setfloat(&DescAtoms[0], (t_float)DescValue);
            outlet_anything(x->DescOut, gensym(x->OpenScofo->GetDescriptionId(d)), DescValue, DescAtoms.data());
        }
    }
}

// ─────────────────────────────────────
static std::vector<OpenScofo::Descriptors> oscofo_get_descriptors(MaxOpenScofo *x, int start, int argc, t_atom *argv) {
    std::vector<OpenScofo::Descriptors> Descriptors;
    Descriptors.reserve(argc - start);
    for (int i = start; i < argc; i++) {
        t_symbol *sym = atom_getsym(argv + i);
        Descriptors.push_back(x->OpenScofo->GetDescriptorsEnum(sym->s_name));
    }
    return Descriptors;
}

// ─────────────────────────────────────
static void oscofo_get(MaxOpenScofo *x, t_symbol *s, long argc, t_atom *argv) {
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
    if (method == "descriptors") {
        if (argc < 2 || argv[1].a_type != A_SYM) {
            object_error((t_object *)x, "descriptors method requires <buffer~name>");
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

        OpenScofo::Description Desc = x->OpenScofo->GetAudioDescription(audioBuffer);
        oscofo_output_descriptors(x, Desc);
    }
}

// ─────────────────────────────────────
static void oscofo_set(MaxOpenScofo *x, t_symbol *s, long argc, t_atom *argv) {
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

        const char *path = atom_getsym(argv + 1)->s_name;
        std::vector<OpenScofo::Descriptors> desc = oscofo_get_descriptors(x, 2, argc, argv);

        if (desc.size() == 0) {
            object_error((t_object *)x, "ONNX models require the descriptor order used in training");
            return;
        }

        x->OpenScofo->LoadONNXModel(path, desc);
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
        object_error((t_object *)x, "Section method not implemented");
    } else if (method == "justdescription") {
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
static void oscofo_following(MaxOpenScofo *x, long f) {
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
static void oscofo_luaexecute(MaxOpenScofo *x, std::string code) {
#ifdef OSCOFO_LUA
    if (!x->OpenScofo->LuaExecute(code)) {
        std::string error = x->OpenScofo->LuaGetError();
        object_error((t_object *)x, "Lua error");
        object_error((t_object *)x, "%s", error.c_str());
    }
#endif
}

// ─────────────────────────────────────
static void oscofo_maxsend(MaxOpenScofo *x, std::string r, int argc, t_atom *argv) {
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
static t_atom *oscofo_convertargs(OpenScofo::ScoreAction &action) {
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
static void oscofo_tickactions(void *xv) {
    MaxOpenScofo *x = (MaxOpenScofo *)xv;
    const double currentTime = gettime();
    const double nextBlock = 1000.0 / x->Sr * x->BlockSize;
    const double nextTime = currentTime + nextBlock;

    std::vector<MaxAction>::iterator it = x->Actions.begin();
    while (it != x->Actions.end()) {
        MaxAction &curAction = *it;
        if (currentTime <= curAction.time && curAction.time <= nextTime) {
            if (curAction.isLua) {
                oscofo_luaexecute(x, curAction.LuaCode);
            } else {
                oscofo_maxsend(x, curAction.Receiver, curAction.MaxArgsSize, curAction.MaxArgs);
                delete[] curAction.MaxArgs;
            }
            it = x->Actions.erase(it);
        } else {
            ++it;
        }
    }
}

// ─────────────────────────────────────
static void oscofo_tickinfo(MaxOpenScofo *xv) {
    MaxOpenScofo *x = (MaxOpenScofo *)xv;
    if (x->MirOutput) {
        OpenScofo::Description Desc;
        if (x->JustDescription && x->Desc) {
            Desc = *x->Desc;
        } else {
            Desc = x->OpenScofo->GetDescription();
        }
        oscofo_output_descriptors(x, Desc);
    }
}

// ─────────────────────────────────────
static void oscofo_ticknewevent(void *xv) {
    MaxOpenScofo *x = (MaxOpenScofo *)xv;

    int prevEvent = x->Event;
    x->Event = x->OpenScofo->GetEventIndex();
    if (prevEvent == x->Event || x->Event == 0) {
        return;
    }

    outlet_float(x->TempoOut, x->OpenScofo->GetLiveBPM());
    outlet_float(x->EventOut, x->OpenScofo->GetEventIndex());
    OpenScofo::EventActions actions = x->OpenScofo->GetEventActions(x->Event);

    for (OpenScofo::ScoreAction &act : actions) {
        double time = act.Time;
        if (!act.AbsoluteTime) {
            act.Time = 60.0 / x->OpenScofo->GetLiveBPM() * act.Time * 1000;
            time = act.Time;
        }

        if (time == 0) {
            if (act.isLua) {
                oscofo_luaexecute(x, act.Lua);
            } else {
                t_atom *maxArgs = oscofo_convertargs(act);
                oscofo_maxsend(x, act.Receiver, act.Args.size(), maxArgs);
                delete[] maxArgs;
            }
        } else {
            double sysTime = gettime() + time;
            int size = act.Args.size();
            std::string receiver = act.Receiver;
            t_atom *maxArgs = oscofo_convertargs(act);
            MaxAction action = {sysTime, act.isLua, receiver, act.Lua, maxArgs, size};
            x->Actions.push_back(action);
        }
    }
}

// ─────────────────────────────────────
static void oscofo_perform_descriptors(MaxOpenScofo *x, double **ins, long sampleframes) {
    x->BlockIndex += sampleframes;
    std::copy(x->inBuffer.begin() + sampleframes, x->inBuffer.end(), x->inBuffer.begin());
    std::copy(ins[0], ins[0] + sampleframes, x->inBuffer.end() - sampleframes);

    if (x->BlockIndex != x->HopSize) {
        clock_delay(x->ClockActions, 0);
        return;
    }

    x->BlockIndex = 0;
    x->Desc = std::make_unique<OpenScofo::Description>(x->OpenScofo->GetAudioDescription(x->inBuffer));
    clock_delay(x->ClockInfo, 0);
}

// ─────────────────────────────────────
static void oscofo_perform_score(MaxOpenScofo *x, double **ins, long sampleframes) {
    x->BlockIndex += sampleframes;
    std::copy(x->inBuffer.begin() + sampleframes, x->inBuffer.end(), x->inBuffer.begin());
    std::copy(ins[0], ins[0] + sampleframes, x->inBuffer.end() - sampleframes);

    if (x->BlockIndex != x->HopSize) {
        clock_delay(x->ClockActions, 0);
        return;
    }

    if (!x->OpenScofo->ScoreIsLoaded() || !x->Following) {
        return;
    }

    x->BlockIndex = 0;
    bool ok = x->OpenScofo->ProcessBlock(x->inBuffer);
    if (!ok) {
        return;
    }

    clock_delay(x->ClockActions, 0);
    clock_delay(x->ClockEvent, 0);
    clock_delay(x->ClockInfo, 0);
}

// ─────────────────────────────────────
static void oscofo_perform64(t_object *obj, t_object *dsp64, double **ins, long numins, double **outs, long numouts,
                             long sampleframes, long flags, void *userparam) {
    auto *x = (MaxOpenScofo *)obj;
    (void)dsp64;
    (void)numins;
    (void)outs;
    (void)numouts;
    (void)flags;
    (void)userparam;

    if (x->JustDescription) {
        oscofo_perform_descriptors(x, ins, sampleframes);
    } else {
        oscofo_perform_score(x, ins, sampleframes);
    }
}

// ─────────────────────────────────────
static void oscofo_dsp64(MaxOpenScofo *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize,
                         long flags) {
    (void)count;
    (void)flags;
    x->BlockSize = maxvectorsize;
    x->BlockIndex = 0;
    x->Sr = (int)samplerate;
    x->inBuffer.resize((size_t)x->FFTSize, 0.0f);
    dsp_add64(dsp64, (t_object *)x, oscofo_perform64, 0, nullptr);
}

// ─────────────────────────────────────
static void oscofo_error_callback(const spdlog::details::log_msg &log, void *data) {
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
static void oscofo_assist(MaxOpenScofo *x, void *b, long m, long a, char *s) {
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
static void *oscofo_new(t_symbol *s, long argc, t_atom *argv) {
    MaxOpenScofo *x = (MaxOpenScofo *)object_alloc(oscofo_class);
    (void)s;

    if (!x) {
        return nullptr;
    }

    dsp_setup((t_pxobject *)x, 1);

    x->EventOut = outlet_new(x, "int");
    x->TempoOut = outlet_new(x, "float");

    x->RequestMIR = oscofo_get_descriptors(x, 0, argc, argv);
    if (x->RequestMIR.size() > 0) {
        x->DescOut = outlet_new(x, nullptr);
        x->MirOutput = true;
    } else {
        x->DescOut = outlet_new(x, nullptr);
        x->MirOutput = false;
    }

    x->ClockEvent = clock_new(x, (method)oscofo_ticknewevent);
    x->ClockActions = clock_new(x, (method)oscofo_tickactions);
    x->ClockInfo = clock_new(x, (method)oscofo_tickinfo);

    x->FFTSize = 2048;
    x->HopSize = 512;
    x->Sr = (int)sys_getsr();
    x->Following = false;
    x->Event = -1;
    x->BlockIndex = 0;
    x->JustDescription = false;
    x->log = spdlog::level::warn;

    char patchPath[MAX_PATH_CHARS];
    short pathId = path_getdefault();
    path_toabsolutesystempath(pathId, nullptr, patchPath);
    x->PatchDir = patchPath;

    x->OpenScofo = new OpenScofo::OpenScofo((float)x->Sr, (float)x->FFTSize, (float)x->HopSize);
    x->OpenScofo->SetErrorCallback(oscofo_error_callback, static_cast<void *>(x));
    x->OpenScofo->SetLogLevel(x->log);

#ifdef OSCOFO_LUA
    x->OpenScofo->LuaAddModule("max", luaopen_max);
    x->OpenScofo->LuaAddPath(x->PatchDir);
    x->OpenScofo->LuaAddPointer(x, "_maxobj");
#endif

    return (void *)x;
}

// ─────────────────────────────────────
static void oscofo_free(MaxOpenScofo *x) {
    for (MaxAction &action : x->Actions) {
        if (!action.isLua && action.MaxArgs) {
            delete[] action.MaxArgs;
        }
    }
    x->Actions.clear();

    if (x->ClockEvent) {
        object_free(x->ClockEvent);
    }
    if (x->ClockActions) {
        object_free(x->ClockActions);
    }
    if (x->ClockInfo) {
        object_free(x->ClockInfo);
    }

    delete x->OpenScofo;
    dsp_free((t_pxobject *)x);
}

// ─────────────────────────────────────
void ext_main(void) {
    t_class *c =
        class_new("o.scofo~", (method)oscofo_new, (method)oscofo_free, (long)sizeof(MaxOpenScofo), 0L, A_GIMME, 0);

    object_post(nullptr, "[o.scofo~] version %d.%d.%d (%s), by Charles K. Neimog", OSCOFO_VERSION_MAJOR,
                OSCOFO_VERSION_MINOR, OSCOFO_VERSION_PATCH, OSCOFO_BUILD_TIME);

    class_addmethod(c, (method)oscofo_score, "score", A_SYM, 0);
    class_addmethod(c, (method)oscofo_start, "start", A_NOTHING, 0);
    class_addmethod(c, (method)oscofo_following, "follow", A_LONG, 0);
    class_addmethod(c, (method)oscofo_set, "set", A_GIMME, 0);
    class_addmethod(c, (method)oscofo_get, "get", A_GIMME, 0);

    class_addmethod(c, (method)stdinletinfo, "inletinfo", A_CANT, 0);
    class_addmethod(c, (method)oscofo_assist, "assist", A_CANT, 0);

    class_addmethod(c, (method)oscofo_dsp64, "dsp64", A_CANT, 0);

    class_dspinit(c);
    class_register(CLASS_BOX, c);
    oscofo_class = c;
}
