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

// ─────────────────────────────────────
struct Action {
    double time;
    bool isLua;
    std::string Receiver;
    std::string LuaCode;
    t_atom *Args;
    int MaxArgsSize;
};

// ─────────────────────────────────────
class MaxOpenScofo {
  public:
    t_pxobject MaxObject;
    t_sample Sample;
    std::string PatchDir;

    enum MIR {
        MFCC = 0,
        LOUDNESS,
        RMS,
        POWER,
        SILENCE,
        CHROMA,
        CQT,
    };

    // Clock
    t_clock *ClockEvent;
    t_clock *ClockInfo;
    t_clock *ClockActions;

    spdlog::level::level_enum log;

    // Actions
    std::vector<Action> Actions;

    // Mir
    std::vector<MIR> RequestMIR;
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
    void *InfoOut;
};

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
static void oscofo_output_descriptiors(MaxOpenScofo *x, OpenScofo::Description &Desc) {
    for (MaxOpenScofo::MIR v : x->RequestMIR) {
        if (v == MaxOpenScofo::MIR::MFCC) {
            size_t mfccSize = Desc.MFCC.size();
            std::vector<t_atom> mfccAtoms(mfccSize);
            for (size_t i = 0; i < mfccSize; ++i) {
                atom_setfloat(&mfccAtoms[i], (float)Desc.MFCC[i]);
            }
            outlet_anything(x->InfoOut, gensym("mfcc"), mfccSize, mfccAtoms.data());
        } else if (v == MaxOpenScofo::MIR::CHROMA) {
            size_t chromaSize = Desc.Chroma.size();
            std::vector<t_atom> chromaAtoms(chromaSize);
            for (size_t i = 0; i < chromaSize; ++i) {
                atom_setfloat(&chromaAtoms[i], (float)Desc.Chroma[i]);
            }
            outlet_anything(x->InfoOut, gensym("chroma"), chromaSize, chromaAtoms.data());
        } else if (v == MaxOpenScofo::MIR::CQT) {
            size_t cqtSize = Desc.PseudoCQT.size();
            std::vector<t_atom> cqtAtoms(cqtSize);
            for (size_t i = 0; i < cqtSize; ++i) {
                atom_setfloat(&cqtAtoms[i], (float)Desc.PseudoCQT[i]);
            }
            outlet_anything(x->InfoOut, gensym("cqt"), cqtSize, cqtAtoms.data());
        } else if (v == MaxOpenScofo::MIR::POWER) {
            size_t powerSize = Desc.Power.size();
            std::vector<t_atom> powerAtoms(powerSize);
            for (size_t i = 0; i < powerSize; ++i) {
                atom_setfloat(&powerAtoms[i], (float)Desc.Power[i]);
            }
            outlet_anything(x->InfoOut, gensym("power"), powerSize, powerAtoms.data());
        } else if (v == MaxOpenScofo::MIR::LOUDNESS) {
            std::vector<t_atom> loudnessAtoms(1);
            atom_setfloat(&loudnessAtoms[0], (float)Desc.Loudness);
            outlet_anything(x->InfoOut, gensym("loudness"), 1, loudnessAtoms.data());
        } else if (v == MaxOpenScofo::MIR::RMS) {
            std::vector<t_atom> rmsAtoms(1);
            atom_setfloat(&rmsAtoms[0], (float)Desc.RMS);
            outlet_anything(x->InfoOut, gensym("rms"), 1, rmsAtoms.data());
        } else if (v == MaxOpenScofo::MIR::SILENCE) {
            std::vector<t_atom> silenceAtoms(1);
            atom_setfloat(&silenceAtoms[0], (float)Desc.SilenceProb);
            outlet_anything(x->InfoOut, gensym("silence"), 1, silenceAtoms.data());
        }
    }
}

// ─────────────────────────────────────
static void oscofo_score(MaxOpenScofo *x, t_symbol *s) {
    // check if file exists
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
    x->Event = -1;
    outlet_float(x->TempoOut, x->OpenScofo->GetLiveBPM());
    outlet_float(x->EventOut, 0);

    // Update Audio
    x->FFTSize = x->OpenScofo->GetFFTSize();
    x->HopSize = x->OpenScofo->GetHopSize();
    x->inBuffer.resize(x->FFTSize, 0.0f);

    // Get Lua Code

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
    x->OpenScofo->SetCurrentEvent(0);
    x->Event = -1;

    // clear actions
    x->Actions.clear();

    outlet_float(x->TempoOut, x->OpenScofo->GetLiveBPM());
    outlet_float(x->EventOut, 0);
    x->Following = true;
    object_post((t_object *)x, "Start following");
}

// ─────────────────────────────────────
static void oscofo_set(MaxOpenScofo *x, t_symbol *s, long argc, t_atom *argv) {
    (void)s;

    if (argc < 1) {
        object_error((t_object *)x, "Wrong number of arguments");
        return;
    }

    if (argv[0].a_type != A_SYM) {
        object_error((t_object *)x, "First argument must be a symbol");
        return;
    }

    std::string method = atom_getsym(argv)->s_name;
    if (method == "event") {
        if (argc < 2) {
            object_error((t_object *)x, "set event requires one argument");
            return;
        }
        long f = atom_getlong(argv + 1);
        x->Event = f;
        x->OpenScofo->SetCurrentEvent(f);
        object_post((t_object *)x, "Event set to %d", (int)f);
    } else if (method == "verbosity") {
        if (argc < 2) {
            object_error((t_object *)x, "set verbosity requires one argument");
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
            object_error((t_object *)x, "set justdescription requires one argument");
            return;
        }
        long f = atom_getlong(argv + 1);
        x->JustDescription = f != 0;
    } else {
        object_error((t_object *)x, "Unknown method");
    }
}

// ─────────────────────────────────────
static void oscofo_get(MaxOpenScofo *x, t_symbol *s, long argc, t_atom *argv) {
    (void)s;

    if (argc < 1) {
        object_error((t_object *)x, "Wrong number of arguments");
        return;
    }

    if (argv[0].a_type != A_SYM) {
        object_error((t_object *)x, "First argument of get must be a symbol");
        return;
    }

    std::string method = atom_getsym(argv)->s_name;
    if (method != "descriptors") {
        object_error((t_object *)x, "Unknown get method");
        return;
    }

    if (argc < 2 || argv[1].a_type != A_SYM) {
        object_error((t_object *)x, "get descriptors requires a buffer~ name");
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
    if (argc >= 3) {
        start = atom_getlong(argv + 2);
        if (start < 0 || start >= frameCount) {
            object_error((t_object *)x, "invalid start index %ld for buffer size %ld", start, frameCount);
            object_free(bufferRef);
            return;
        }
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
        long src = start + i;
        if (src >= frameCount) {
            break;
        }
        audioBuffer[i] = static_cast<double>(samples[src]);
    }

    buffer_unlocksamples(bufferObj);
    object_free(bufferRef);

    OpenScofo::Description Desc = x->OpenScofo->GetAudioDescription(audioBuffer);
    oscofo_output_descriptiors(x, Desc);
}

// ─────────────────────────────────────
static void oscofo_following(MaxOpenScofo *x, long f) {
    if (!x->OpenScofo->ScoreIsLoaded()) {
        object_error((t_object *)x, "Score not loaded");
        return;
    }
    if (f == 1) {
        x->Following = true;
        object_post((t_object *)x, "Following!");
    } else {
        object_post((t_object *)x, "Not Following!");
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
    t_symbol *sym = gensym(r.c_str());
    t_object *receiver = sym->s_thing;
    if (receiver == nullptr) {
        object_error((t_object *)x, "Receiver '%s' not found", r.c_str());
        return;
    }

    if (argc == 0) {
        object_method_typed(receiver, gensym("bang"), 0, nullptr, nullptr);
    } else {
        object_method_typed(receiver, gensym("list"), argc, argv, nullptr);
    }
}
// ─────────────────────────────────────
static t_atom *oscofo_convertargs(MaxOpenScofo *x, OpenScofo::Action &action) {
    (void)x;
    int size = action.Args.size();
    t_atom *MaxArgs = new t_atom[size];

    for (int i = 0; i < size; i++) {
        std::variant<float, int, std::string> arg = action.Args[i];
        if (std::holds_alternative<float>(arg)) {
            atom_setfloat(&MaxArgs[i], std::get<float>(arg));
        } else if (std::holds_alternative<int>(arg)) {
            atom_setlong(&MaxArgs[i], std::get<int>(arg));
        } else if (std::holds_alternative<std::string>(arg)) {
            atom_setsym(&MaxArgs[i], gensym(std::get<std::string>(arg).c_str()));
        }
    }
    return MaxArgs;
}

// ─────────────────────────────────────
static void oscofo_tickactions(MaxOpenScofo *x) {
    const double CurrentTime = gettime();
    const double nextBlock = 1000.0 / x->Sr * x->BlockSize;
    const double NextTime = CurrentTime + nextBlock;

    std::vector<Action>::iterator it = x->Actions.begin();
    while (it != x->Actions.end()) {
        Action &CurAction = *it;
        if ((CurrentTime <= CurAction.time && CurAction.time <= NextTime) || CurAction.time < CurrentTime) {
            if (CurAction.isLua) {
                oscofo_luaexecute(x, CurAction.LuaCode);
            } else {
                oscofo_maxsend(x, CurAction.Receiver, CurAction.MaxArgsSize, CurAction.Args);
                delete[] CurAction.Args;
            }
            it = x->Actions.erase(it);
        } else {
            ++it;
        }
    }
}

// ─────────────────────────────────────
static void oscofo_tickinfo(MaxOpenScofo *x) {
    if (x->MirOutput) {
        OpenScofo::Description Desc;
        if (x->JustDescription && x->Desc) {
            Desc = *x->Desc;
        } else {
            Desc = x->OpenScofo->GetDescription();
        }
        oscofo_output_descriptiors(x, Desc);
    }
}

// ─────────────────────────────────────
static void oscofo_ticknewevent(MaxOpenScofo *x) {
    int PrevEvent = x->Event;
    x->Event = x->OpenScofo->GetEventIndex();
    if (PrevEvent == x->Event || x->Event == 0) {
        return;
    }

    outlet_float(x->TempoOut, x->OpenScofo->GetLiveBPM());
    outlet_float(x->EventOut, x->OpenScofo->GetEventIndex());
    OpenScofo::ActionVec Actions = x->OpenScofo->GetEventActions(x->Event);

    for (OpenScofo::Action &Act : Actions) {
        double time = Act.Time;
        if (!Act.AbsoluteTime) {
            Act.Time = 60.0 / x->OpenScofo->GetLiveBPM() * Act.Time * 1000;
            time = Act.Time;
        }
        if (time == 0) {
            if (Act.isLua) {
                oscofo_luaexecute(x, Act.Lua);
            } else {
                t_atom *MaxArgs = oscofo_convertargs(x, Act);
                oscofo_maxsend(x, Act.Receiver, Act.Args.size(), MaxArgs);
                delete[] MaxArgs;
            }
        } else {
            double actionTime = gettime() + time;
            int size = Act.Args.size();
            std::string receiver = Act.Receiver;
            t_atom *PdArgs = oscofo_convertargs(x, Act);
            Action action = {actionTime, Act.isLua, receiver, Act.Lua, PdArgs, size};
            x->Actions.push_back(action);
        }
    }
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

    x->BlockIndex += sampleframes;
    std::copy(x->inBuffer.begin() + sampleframes, x->inBuffer.end(), x->inBuffer.begin());
    std::copy(ins[0], ins[0] + sampleframes, x->inBuffer.end() - sampleframes);
    if (x->BlockIndex != x->HopSize) {
        clock_delay(x->ClockActions, 0);
        return;
    }

    if (x->JustDescription) {
        x->BlockIndex = 0;
        x->Desc = std::make_unique<OpenScofo::Description>(x->OpenScofo->GetAudioDescription(x->inBuffer));
        clock_delay(x->ClockInfo, 0);
        return;
    }

    if (!x->OpenScofo->ScoreIsLoaded() || !x->Following) {
        x->BlockIndex = 0;
        return;
    }

    x->BlockIndex = 0;
    bool ok = x->OpenScofo->ProcessBlock(x->inBuffer);
    if (!ok) {
        x->OpenScofo->ClearErrors();
        return;
    }
    clock_delay(x->ClockActions, 0);
    clock_delay(x->ClockEvent, 0);
    clock_delay(x->ClockInfo, 0);
    return;
}

// ─────────────────────────────────────
static void oscofo_dsp64(MaxOpenScofo *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize,
                         long flags) {
    (void)count;
    (void)flags;
    x->BlockSize = maxvectorsize;
    x->BlockIndex = 0;
    x->Sr = samplerate;
    x->inBuffer.resize(x->FFTSize, 0.0f);
    dsp_add64(dsp64, (t_object *)x, oscofo_perform64, 0, nullptr);
}

// ─────────────────────────────────────
static void *oscofo_new(t_symbol *s, long argc, t_atom *argv) {
    MaxOpenScofo *x = (MaxOpenScofo *)object_alloc(oscofo_class);
    (void)s;
    if (!x) {
        object_error((t_object *)x, "Error creating object");
        return nullptr;
    }

    dsp_setup((t_pxobject *)x, 1);

    x->InfoOut = outlet_new(x, nullptr);
    x->TempoOut = outlet_new(x, "float");
    x->EventOut = outlet_new(x, "int");
    x->ClockEvent = clock_new(x, (method)oscofo_ticknewevent);
    x->ClockActions = clock_new(x, (method)oscofo_tickactions);
    x->ClockInfo = clock_new(x, (method)oscofo_tickinfo);
    x->FFTSize = 4096.0f;
    x->HopSize = 1024.0f;
    x->Sr = sys_getsr();
    x->Following = false;
    x->Event = -1;
    x->log = spdlog::level::warn;
    x->JustDescription = false;

    char PatchPath[MAX_PATH_CHARS];
    short PathId = path_getdefault();
    path_toabsolutesystempath(PathId, NULL, PatchPath);
    x->PatchDir = PatchPath;

    x->OpenScofo = new OpenScofo::OpenScofo(x->Sr, x->FFTSize, x->HopSize);
    x->OpenScofo->SetErrorCallback(oscofo_error_callback, static_cast<void *>(x));
    x->OpenScofo->SetLogLevel(x->log);

#ifdef OSCOFO_LUA
    x->OpenScofo->LuaAddModule("max", luaopen_max);
    x->OpenScofo->LuaAddPath(x->PatchDir);
    x->OpenScofo->LuaAddPointer(x, "_maxobj");
#endif

    // Args
    while (argc) {
        if ((argv)->a_type == A_SYM) {
            t_symbol *sym = atom_getsym(argv);
            if (strcmp(sym->s_name, "mfcc") == 0) {
                x->MirOutput = true;
                x->RequestMIR.push_back(MaxOpenScofo::MIR::MFCC);
            } else if (strcmp(sym->s_name, "rms") == 0) {
                x->MirOutput = true;
                x->RequestMIR.push_back(MaxOpenScofo::MIR::RMS);
            } else if (strcmp(sym->s_name, "loudness") == 0) {
                x->MirOutput = true;
                x->RequestMIR.push_back(MaxOpenScofo::MIR::LOUDNESS);
            } else if (strcmp(sym->s_name, "silence") == 0) {
                x->MirOutput = true;
                x->RequestMIR.push_back(MaxOpenScofo::MIR::SILENCE);
            } else if (strcmp(sym->s_name, "cqt") == 0) {
                x->MirOutput = true;
                x->RequestMIR.push_back(MaxOpenScofo::MIR::CQT);
            } else if (strcmp(sym->s_name, "chroma") == 0) {
                x->MirOutput = true;
                x->RequestMIR.push_back(MaxOpenScofo::MIR::CHROMA);
            } else if (strcmp(sym->s_name, "power") == 0) {
                x->MirOutput = true;
                x->RequestMIR.push_back(MaxOpenScofo::MIR::POWER);
            }
        }
        argc--, argv++;
    }

    return (x);
}

// ─────────────────────────────────────
static void oscofo_free(MaxOpenScofo *x) {
    for (Action &action : x->Actions) {
        if (!action.isLua && action.Args) {
            delete[] action.Args;
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
    object_post(nullptr, "[oscofo~] version %d.%d.%d, by Charles K. Neimog", OSCOFO_VERSION_MAJOR, OSCOFO_VERSION_MINOR,
                OSCOFO_VERSION_PATCH);
    // message methods
    class_addmethod(c, (method)oscofo_set, "set", A_GIMME, 0);
    class_addmethod(c, (method)oscofo_get, "get", A_GIMME, 0);
    class_addmethod(c, (method)oscofo_score, "score", A_SYM, 0);
    class_addmethod(c, (method)oscofo_following, "follow", A_LONG, 0);
    class_addmethod(c, (method)oscofo_start, "start", A_NOTHING, 0);

    // user methods
    class_addmethod(c, (method)stdinletinfo, "inletinfo", A_CANT, 0);
    class_addmethod(c, (method)oscofo_assist, "assist", A_CANT, 0);

    // dsp methods
    class_addmethod(c, (method)oscofo_dsp64, "dsp64", A_CANT, 0);

    // register the class
    class_dspinit(c);
    class_register(CLASS_BOX, c);
    oscofo_class = c;
}
