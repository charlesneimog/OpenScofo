#include <filesystem>

#include <m_pd.h>
#include <OpenScofo.hpp>

static t_class *OpenScofoObj;

#ifdef OSCOFO_LUA
int luaopen_pd(lua_State *L);
#endif

// ─────────────────────────────────────
struct Action {
    double time;
    bool isLua;
    std::string Receiver;
    std::string LuaCode;
    t_atom *MaxArgs;
    int PdArgsSize;
};

// ─────────────────────────────────────
class PdOpenScofo {
  public:
    t_object PdObject;
    t_sample Sample;
    t_canvas *Canvas;
    std::string PatchDir;

    enum MIR {
        MFCC = 0,
        LOUDNESS,
        RMS,
        POWER,
        SILENCE,
        CQT,
    };

    // Clock
    t_clock *ClockEvent;
    t_clock *ClockActions;
    t_clock *ClockInfo;

    spdlog::level::level_enum log;

    // Actions
    std::vector<Action> Actions;

    // Mir
    std::vector<MIR> RequestMIR;
    bool MirOutput = false;

    // OpenScofo
    OpenScofo::OpenScofo *OpenScofo;
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

    // Outlet
    t_outlet *EventOut;
    t_outlet *TempoOut;
    t_outlet *DescOut;
};

// ─────────────────────────────────────
static void oscofo_score(PdOpenScofo *x, t_symbol *s) {
    // check if file exists
    if (!s) {
        pd_error(x, "[o.scofo~] No score file provided");
        return;
    }

    bool ok;
    std::string scorePath = s->s_name;
    if (!std::filesystem::exists(s->s_name)) {
        scorePath = x->PatchDir + "/" + s->s_name;
    }

    ok = x->OpenScofo->ParseScore(scorePath);
    if (ok) {
        logpost(x, 2, "[o.scofo~] Score loaded");
    } else {
        return;
    }
    x->OpenScofo->SetCurrentEvent(0);

    x->Event = 0;
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
        pd_error(x, "[o.scofo~] Lua error");
        pd_error(x, "[o.scofo~] %s", error.c_str());
        logpost(x, 1, "");
    }
#endif
}

// ─────────────────────────────────────
static void oscofo_start(PdOpenScofo *x) {
    if (!x->OpenScofo->ScoreIsLoaded()) {
        pd_error(x, "[o.scofo~] Score not loaded");
        return;
    }
    x->Actions.clear();

    outlet_float(x->TempoOut, x->OpenScofo->GetLiveBPM());
    x->OpenScofo->SetCurrentEvent(0);
    x->Event = 0;

    x->Following = true;
    logpost(x, 2, "[o.scofo~] Start following");
}

// ─────────────────────────────────────
static void oscofo_set(PdOpenScofo *x, t_symbol *s, int argc, t_atom *argv) {
    if (argv[0].a_type != A_SYMBOL) {
        pd_error(x, "[o.scofo~] First argument of set method must be a symbol");
        return;
    }
    std::string method = atom_getsymbol(argv)->s_name;
    if (method == "event") {
        int f = atom_getint(argv + 1);
        x->Event = f;
        x->OpenScofo->SetCurrentEvent(f);
    } else if (method == "verbosity") {
        int f = atom_getint(argv + 1);
        switch (f) {
        case 0: {
            x->log = spdlog::level::warn;
            break;
        }
        case 1:
            x->log = spdlog::level::info;
            break;
        case 2:
            x->log = spdlog::level::debug;
            break;
        case 3:
            x->log = spdlog::level::trace;
            break;
        }

    } else if (method == "section") {
        pd_error(x, "[o.scofo~] Section method not implemented");
    } else {
        pd_error(x, "[o.scofo~] Unknown method");
    }
}

// ─────────────────────────────────────
static void oscofo_following(PdOpenScofo *x, t_float f) {
    if (!x->OpenScofo->ScoreIsLoaded()) {
        pd_error(x, "[o.scofo~] Score not loaded");
        return;
    }
    if (f == 1) {
        x->Following = true;
    } else {
        x->Following = false;
    }
}
// ─────────────────────────────────────
static void oscofo_luaexecute(PdOpenScofo *x, std::string code) {
#if OSCOFO_LUA
    if (!x->OpenScofo->LuaExecute(code)) {
        std::string error = x->OpenScofo->LuaGetError();
        pd_error(x, "[o.scofo~] Lua error");
        pd_error(x, "[o.scofo~] %s", error.c_str());
    }
#endif
}

// ─────────────────────────────────────
static void oscofo_pdsend(PdOpenScofo *x, std::string r, int argc, t_atom *argv) {
    t_pd *receiver = gensym(r.c_str())->s_thing;
    if (!receiver) {
        pd_error(x, "[o.scofo~] Receiver %s not found", r.c_str());
        return;
    }

    if (argc == 0) {
        pd_bang(receiver);
    } else {
        pd_list(receiver, &s_list, argc, argv);
    }
}

// ─────────────────────────────────────
static t_atom *oscofo_convertargs(PdOpenScofo *x, OpenScofo::Action &action) {
    int size = action.Args.size();
    t_atom *PdArgs = new t_atom[size];

    for (int i = 0; i < size; i++) {
        std::variant<float, int, std::string> arg = action.Args[i];
        if (std::holds_alternative<float>(arg)) {
            SETFLOAT(&PdArgs[i], std::get<float>(arg));
        } else if (std::holds_alternative<int>(arg)) {
            SETFLOAT(&PdArgs[i], std::get<int>(arg));
        } else if (std::holds_alternative<std::string>(arg)) {
            SETSYMBOL(&PdArgs[i], gensym(std::get<std::string>(arg).c_str()));
        }
    }
    return PdArgs;
}

// ─────────────────────────────────────
static void oscofo_tickactions(PdOpenScofo *x) {
    const double CurrentTime = clock_getlogicaltime();
    const double nextBlock = 1000.0 / x->Sr * x->BlockSize;
    const double NextTime = clock_getsystimeafter(nextBlock);

    std::vector<Action>::iterator it = x->Actions.begin();
    while (it != x->Actions.end()) {
        Action &CurAction = *it;
        if (CurrentTime <= CurAction.time && CurAction.time <= NextTime) {
            if (CurAction.isLua) {
                oscofo_luaexecute(x, CurAction.LuaCode);
            } else {
                oscofo_pdsend(x, CurAction.Receiver, CurAction.PdArgsSize, CurAction.MaxArgs);
                delete[] CurAction.MaxArgs;
            }
            it = x->Actions.erase(it);
        } else {
            ++it;
        }
    }
}

// ─────────────────────────────────────
static void oscofo_tickinfo(PdOpenScofo *x) {
    if (x->MirOutput) {
        OpenScofo::Description Desc = x->OpenScofo->GetDescription();
        for (PdOpenScofo::MIR v : x->RequestMIR) {
            if (v == PdOpenScofo::MIR::MFCC) {
                size_t mfccSize = Desc.MFCC.size();
                std::vector<t_atom> mfccAtoms(mfccSize);
                for (size_t i = 0; i < mfccSize; ++i) {
                    SETFLOAT(&mfccAtoms[i], (t_float)Desc.MFCC[i]);
                }
                outlet_anything(x->DescOut, gensym("mfcc"), mfccSize, mfccAtoms.data());
            } else if (v == PdOpenScofo::MIR::CQT) {
                size_t mfccSize = Desc.PseudoCQT.size();
                std::vector<t_atom> Atoms(mfccSize);
                for (size_t i = 0; i < mfccSize; ++i) {
                    SETFLOAT(&Atoms[i], (t_float)Desc.PseudoCQT[i]);
                }
                outlet_anything(x->DescOut, gensym("cqt"), mfccSize, Atoms.data());
            } else if (v == PdOpenScofo::MIR::POWER) {
                size_t mfccSize = Desc.Power.size();
                std::vector<t_atom> mfccAtoms(mfccSize);
                for (size_t i = 0; i < mfccSize; ++i) {
                    SETFLOAT(&mfccAtoms[i], (t_float)Desc.MFCC[i]);
                }
                outlet_anything(x->DescOut, gensym("power"), mfccSize, mfccAtoms.data());
            } else if (v == PdOpenScofo::MIR::LOUDNESS) {
                std::vector<t_atom> mfccAtoms(1);
                SETFLOAT(&mfccAtoms[0], (t_float)Desc.Loudness);
                outlet_anything(x->DescOut, gensym("loudness"), 1, mfccAtoms.data());
            } else if (v == PdOpenScofo::MIR::RMS) {
                std::vector<t_atom> mfccAtoms(1);
                SETFLOAT(&mfccAtoms[0], (t_float)Desc.RMS);
                outlet_anything(x->DescOut, gensym("rms"), 1, mfccAtoms.data());
            } else if (v == PdOpenScofo::MIR::SILENCE) {
                std::vector<t_atom> mfccAtoms(1);
                SETFLOAT(&mfccAtoms[0], (t_float)Desc.SilenceProb);
                outlet_anything(x->DescOut, gensym("silence"), 1, mfccAtoms.data());
            }
        }
    }
}

// ─────────────────────────────────────
static void oscofo_ticknewevent(PdOpenScofo *x) {
    int PrevEvent = x->Event;
    x->Event = x->OpenScofo->GetEventIndex();
    if (PrevEvent == x->Event || x->Event == 0) {
        return;
    }

    outlet_float(x->TempoOut, x->OpenScofo->GetLiveBPM());
    outlet_float(x->EventOut, x->OpenScofo->GetEventIndex());
    OpenScofo::ActionVec Actions = x->OpenScofo->GetEventActions(x->Event - 1);

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
                t_atom *PdArgs = oscofo_convertargs(x, Act);
                oscofo_pdsend(x, Act.Receiver, Act.Args.size(), PdArgs);
                delete[] PdArgs;
            }
        } else {
            double sysTime = clock_getsystimeafter(time);
            int size = Act.Args.size();
            std::string receiver = Act.Receiver;
            t_atom *PdArgs = oscofo_convertargs(x, Act);
            Action action = {sysTime, Act.isLua, receiver, Act.Lua, PdArgs, size};
            x->Actions.push_back(action);
        }
    }
}

// ─────────────────────────────────────
static t_int *oscofo_perform(t_int *w) {
    PdOpenScofo *x = (PdOpenScofo *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    int n = static_cast<int>(w[3]);

    if (!x->OpenScofo->ScoreIsLoaded() || !x->Following) {
        return (w + 4);
    }

    x->BlockIndex += n;
    std::copy(x->inBuffer.begin() + n, x->inBuffer.end(), x->inBuffer.begin());
    std::copy(in, in + n, x->inBuffer.end() - n);

    if (x->BlockIndex != x->HopSize) {
        clock_delay(x->ClockActions, 0);
        return (w + 4);
    }

    x->BlockIndex = 0;
    bool ok = x->OpenScofo->ProcessBlock(x->inBuffer);
    if (!ok) {
        return (w + 4);
    }

    clock_delay(x->ClockActions, 0);
    clock_delay(x->ClockEvent, 0);
    clock_delay(x->ClockInfo, 0);
    return (w + 4);
}

// ─────────────────────────────────────
static void oscofo_adddsp(PdOpenScofo *x, t_signal **sp) {
    x->BlockSize = sp[0]->s_n;
    x->BlockIndex = 0;
    x->inBuffer.resize((size_t)x->FFTSize, (double)0.0f);
    dsp_add(oscofo_perform, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

// ─────────────────────────────────────
static void oscofo_callback(const spdlog::details::log_msg &log, void *data) {
    PdOpenScofo *x = static_cast<PdOpenScofo *>(data);
    spdlog::level::level_enum pdlevel = x->log;
    if (log.level < pdlevel) {
        return;
    }
    std::string text(log.payload.data(), log.payload.size());
    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
        logpost(x, 1, "[o.scofo~] %s", text.c_str());
        break;
    case spdlog::level::info:
    case spdlog::level::warn:
        logpost(x, 2, "[o.scofo~] %s", text.c_str());
        break;
    case spdlog::level::debug:
    case spdlog::level::trace:
        logpost(x, 3, "[o.scofo~] %s", text.c_str());
        break;
    default:
        break;
    }
}

// ─────────────────────────────────────
static void *oscofo_new(t_symbol *s, int argc, t_atom *argv) {
    PdOpenScofo *x = (PdOpenScofo *)pd_new(OpenScofoObj);
    if (!x) {
        pd_error(x, "[o.scofo~] Error creating object");
        return nullptr;
    }

    // default parameters
    x->FFTSize = 4096;
    x->HopSize = 1024;
    x->Sr = (int)sys_getsr();
    x->Following = false;
    x->Event = -1;

    // Outlets
    x->EventOut = outlet_new(&x->PdObject, &s_float);
    x->TempoOut = outlet_new(&x->PdObject, &s_float);

    // Args
    bool DescOut = false;
    while (argc) {
        if ((argv)->a_type == A_SYMBOL) {
            t_symbol *sym = atom_getsymbol(argv);
            if (strcmp(sym->s_name, "mfcc") == 0) {
                DescOut = true;
                x->RequestMIR.push_back(PdOpenScofo::MIR::MFCC);
            } else if (strcmp(sym->s_name, "rms") == 0) {
                DescOut = true;
                x->RequestMIR.push_back(PdOpenScofo::MIR::RMS);
            } else if (strcmp(sym->s_name, "loudness") == 0) {
                DescOut = true;
                x->RequestMIR.push_back(PdOpenScofo::MIR::LOUDNESS);
            } else if (strcmp(sym->s_name, "silence") == 0) {
                DescOut = true;
                x->RequestMIR.push_back(PdOpenScofo::MIR::SILENCE);
            } else if (strcmp(sym->s_name, "cqt") == 0) {
                DescOut = true;
                x->RequestMIR.push_back(PdOpenScofo::MIR::CQT);
            }
        }
        argc--, argv++;
    }

    if (DescOut) {
        x->DescOut = outlet_new(&x->PdObject, &s_list);
        x->MirOutput = true;
    }

    // Schedule
    x->ClockEvent = clock_new(x, (t_method)oscofo_ticknewevent);
    x->ClockActions = clock_new(x, (t_method)oscofo_tickactions);
    x->ClockInfo = clock_new(x, (t_method)oscofo_tickinfo);

    // Current Dir
    x->Canvas = canvas_getcurrent();
    x->PatchDir = canvas_getdir(x->Canvas)->s_name;

    // OpenScofo Library
    x->OpenScofo = new OpenScofo::OpenScofo((float)x->Sr, (float)x->FFTSize, (float)x->HopSize);
    x->OpenScofo->SetErrorCallback(oscofo_callback, static_cast<void *>(x));

    x->log = (spdlog::level::warn);

#ifdef NDEBUG
#else
    // x->OpenScofo->SetLogLevel(spdlog::level::info);
#endif

#ifdef OSCOFO_LUA
    x->OpenScofo->LuaAddModule("pd", luaopen_pd);
    x->OpenScofo->LuaAddPath(x->PatchDir);
    x->OpenScofo->LuaAddPointer(x, "_pdobj");
#endif
    return (void *)x;
}

// ─────────────────────────────────────
static void oscofo_free(PdOpenScofo *x) {
    delete x->OpenScofo;
}

// ─────────────────────────────────────
extern "C" void setup_o0x2escofo_tilde(void) {
    OpenScofoObj = class_new(gensym("o.scofo~"), (t_newmethod)oscofo_new, (t_method)oscofo_free, sizeof(PdOpenScofo),
                             CLASS_DEFAULT, A_GIMME, A_NULL);

    post("[o.scofo~] version %d.%d.%d (%s %s), by Charles K. Neimog", OSCOFO_VERSION_MAJOR, OSCOFO_VERSION_MINOR,
         OSCOFO_VERSION_PATCH, __DATE__, __TIME__);

    // message methods
    class_addmethod(OpenScofoObj, (t_method)oscofo_score, gensym("score"), A_SYMBOL, 0);
    class_addmethod(OpenScofoObj, (t_method)oscofo_start, gensym("start"), A_NULL, 0);
    class_addmethod(OpenScofoObj, (t_method)oscofo_following, gensym("follow"), A_FLOAT, 0);
    class_addmethod(OpenScofoObj, (t_method)oscofo_set, gensym("set"), A_GIMME, 0);

    // dsp
    CLASS_MAINSIGNALIN(OpenScofoObj, PdOpenScofo, Sample);
    class_addmethod(OpenScofoObj, (t_method)oscofo_adddsp, gensym("dsp"), A_CANT, 0);
}
