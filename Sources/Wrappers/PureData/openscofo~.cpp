#include <filesystem>

#include <m_pd.h>
#include <g_canvas.h>
#include <m_imp.h>
#include <s_stuff.h>

#include <OpenScofo.hpp>

static t_class *OpenScofoObj;

#ifdef OPENSCOFO_LUA
int luaopen_pd(lua_State *L);
#endif

// ─────────────────────────────────────
struct PdAction {
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

    // Clock
    t_clock *ClockEvent;
    t_clock *ClockActions;
    t_clock *ClockInfo;

    spdlog::level::level_enum log;

    // Actions
    std::vector<PdAction> Actions;

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
    int FFTSize;
    int HopSize;
    int BlockSize;
    int Sr;
    int BlockIndex;
    bool JustDescription;

    // Outlet
    t_outlet *EventOut;
    t_outlet *TempoOut;
    t_outlet *DescOut;
};

// ─────────────────────────────────────
static void oscofo_score(PdOpenScofo *x, t_symbol *s) {
    // check if file exists
    if (!s) {
        pd_error(x, "[openscofo~] No score file provided");
        return;
    }

    bool ok;
    std::string scorePath = s->s_name;
    if (!std::filesystem::exists(s->s_name)) {
        scorePath = x->PatchDir + "/" + s->s_name;
    }

    ok = x->OpenScofo->LoadScore(scorePath);
    if (ok) {
        logpost(x, 2, "[openscofo~] Score loaded");
    } else {
        logpost(x, 1, "[openscofo~] Score has errors");
        return;
    }
    x->OpenScofo->SetCurrentEvent(0);
    x->JustDescription = false;

    x->Event = 0;
    outlet_float(x->TempoOut, x->OpenScofo->GetCurrentBPM());
    outlet_float(x->EventOut, 0);

    // Update Audio
    x->FFTSize = x->OpenScofo->GetFFTSize();
    x->HopSize = x->OpenScofo->GetHopSize();

    // Get Lua Code

#ifdef OPENSCOFO_LUA
    std::string LuaCode = x->OpenScofo->GetLuaCode();
    bool result = x->OpenScofo->LuaExecute(LuaCode.c_str());

    if (!result) {
        std::string error = x->OpenScofo->LuaGetError();
        pd_error(x, "[openscofo~] Lua error");
        pd_error(x, "[openscofo~] %s", error.c_str());
    }
#endif
}

// ─────────────────────────────────────
static void oscofo_start(PdOpenScofo *x) {
    if (!x->OpenScofo->ScoreIsLoaded()) {
        pd_error(x, "[openscofo~] Score not loaded");
        return;
    }
    x->Actions.clear();
    x->OpenScofo->SetCurrentEvent(0);
    x->Event = 0;

    outlet_float(x->TempoOut, x->OpenScofo->GetCurrentBPM());
    outlet_float(x->EventOut, x->Event);

    x->Following = true;
    logpost(x, 2, "[openscofo~] Start following");
}

// ─────────────────────────────────────
static void oscofo_output_descriptors(PdOpenScofo *x, OpenScofo::Description &Desc) {
    for (auto it = x->RequestMIR.rbegin(); it != x->RequestMIR.rend(); ++it) {
        OpenScofo::Descriptors d = *it;
        switch (d) {
        case OpenScofo::Descriptors::INVALID:
            pd_error(x, "[openscofo~] Invalid descriptors");
            break;
        case OpenScofo::Descriptors::MFCC:
        case OpenScofo::Descriptors::CHROMA:
        case OpenScofo::Descriptors::MELOGRAM:
        case OpenScofo::Descriptors::MAGNITUDE: {
            std::vector<double> DescArray = x->OpenScofo->GetDescriptionArray(Desc, d);
            int DescSize = DescArray.size();
            std::vector<t_atom> DescAtoms(DescSize);
            for (int i = 0; i < DescSize; ++i) {
                SETFLOAT(&DescAtoms[i], (t_float)DescArray[i]);
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
            SETFLOAT(&DescAtoms[0], static_cast<t_float>(DescValue));
            outlet_anything(x->DescOut, gensym(x->OpenScofo->GetDescriptionId(d)), 1, DescAtoms.data());
        }
    }
}

// ─────────────────────────────────────
static std::vector<OpenScofo::Descriptors> oscofo_get_descriptors(PdOpenScofo *x, int start, int argc, t_atom *argv) {
    std::vector<OpenScofo::Descriptors> Descriptors;
    Descriptors.reserve(argc - start);
    for (int i = start; i < argc; i++) {
        t_symbol *sym = atom_getsymbol(argv + i);
        Descriptors.push_back(x->OpenScofo->GetDescriptorsEnum(sym->s_name));
    }
    return Descriptors;
}

// ─────────────────────────────────────
static void oscofo_get(PdOpenScofo *x, t_symbol *s, int argc, t_atom *argv) {
    (void)s;
    if (argc < 1) {
        pd_error(x, "[openscofo~] Wrong number of arguments");
        return;
    }

    if (argv[0].a_type != A_SYMBOL) {
        pd_error(x, "[openscofo~] First argument of set method must be a symbol");
        return;
    }

    std::string method = atom_getsymbol(argv)->s_name;
    if (method == "descriptors") {
        if (argc < 2) {
            pd_error(x, "[openscofo~] descriptors method require <arrayname>");
            return;
        }

        t_garray *pdarray;
        const char *arrayname = atom_getsymbol(argv + 1)->s_name;
        t_symbol *pd_symbol = gensym(arrayname);
        if (!(pdarray = (t_garray *)pd_findbyclass(pd_symbol, garray_class))) {
            pd_error(x, "[openscofo~] array %s not found", arrayname);
            return;
        } else {
            int vecsize;
            t_word *vec;
            if (!garray_getfloatwords(pdarray, &vecsize, &vec) || vec == nullptr) {
                pd_error(x, "[openscofo~] failed to read array %s", arrayname);
                return;
            }

            int fftsize = x->OpenScofo->GetFFTSize();
            if (vecsize <= 0) {
                pd_error(x, "[openscofo~] array %s is empty", arrayname);
                return;
            }

            int start = 0;
            if (argc == 2) {
                start = 0;
            } else if (argc == 3) {
                start = atom_getint(argv + 2);
                if (start < 0 || start >= vecsize) {
                    pd_error(x, "[openscofo~] invalid start index %d for array size %d and fftsize %d", start, vecsize,
                             fftsize);
                    return;
                }
            } else {
                pd_error(x, "[openscofo~] Wrong arguments");
                return;
            }

            std::vector<double> AudioBuffer(fftsize, 0.0);
            for (int i = 0; i < fftsize; i++) {
                const int src = start + i;
                if (src >= vecsize) {
                    break;
                }
                AudioBuffer[i] = static_cast<double>(vec[src].w_float);
            }

            x->OpenScofo->ProcessBlock(AudioBuffer.data(), AudioBuffer.size());
            OpenScofo::Description Desc = x->OpenScofo->GetDescription();
            oscofo_output_descriptors(x, Desc);
        }
    }
}

// ─────────────────────────────────────
static void oscofo_set(PdOpenScofo *x, t_symbol *s, int argc, t_atom *argv) {
    (void)s;

    if (argv[0].a_type != A_SYMBOL) {
        pd_error(x, "[openscofo~] First argument of set method must be a symbol");
        return;
    }

    std::string method = atom_getsymbol(argv)->s_name;
    if (method == "event") {
        if (argc > 0) {
            pd_error(x, "[openscofo~] Wrong number of arguments");
            return;
        }

        int f = atom_getint(argv + 1);
        x->Event = f;
        x->OpenScofo->SetCurrentEvent(f);
    } else if (method == "onnxmodel") {
        std::vector<OpenScofo::Descriptors> Desc = oscofo_get_descriptors(x, 2, argc, argv);

        if (Desc.size() == 0) {
            pd_error(x, "[openscofo~] ONNX models required the Descriptors order used on train");
            return;
        }

        char dirbuf[MAXPDSTRING], *nameptr;
        char fullpath[MAXPDSTRING];
        int fd = canvas_open(x->Canvas, atom_getsymbol(argv + 1)->s_name, "", dirbuf, &nameptr, MAXPDSTRING, 1);
        sys_close(fd);
        snprintf(fullpath, MAXPDSTRING, "%s/%s", dirbuf, nameptr);
        x->OpenScofo->LoadONNXModel(fullpath, Desc);
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

    } else if (method == "mfcc") {
    } else if (method == "section") {
        pd_error(x, "[openscofo~] Section method not implemented");
    } else if (method == "justdescription") {
        int f = atom_getint(argv + 1);
        x->JustDescription = f != 0;
        canvas_update_dsp();
    } else {
        pd_error(x, "[openscofo~] Unknown method");
    }
}

// ─────────────────────────────────────
static void oscofo_following(PdOpenScofo *x, t_float f) {
    if (!x->OpenScofo->ScoreIsLoaded()) {
        pd_error(x, "[openscofo~] Score not loaded");
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
#if OPENSCOFO_LUA
    if (!x->OpenScofo->LuaExecute(code)) {
        std::string error = x->OpenScofo->LuaGetError();
        pd_error(x, "[openscofo~] Lua error");
        pd_error(x, "[openscofo~] %s", error.c_str());
    }
#endif
}

// ─────────────────────────────────────
static void oscofo_pdsend(PdOpenScofo *x, std::string r, int argc, t_atom *argv) {
    t_pd *receiver = gensym(r.c_str())->s_thing;
    if (!receiver) {
        pd_error(x, "[openscofo~] Receiver %s not found", r.c_str());
        return;
    }

    if (argc == 0) {
        pd_bang(receiver);
    } else {
        pd_list(receiver, &s_list, argc, argv);
    }
}

// ─────────────────────────────────────
static t_atom *oscofo_convertargs(OpenScofo::ScoreAction &action) {
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

    std::vector<PdAction>::iterator it = x->Actions.begin();
    while (it != x->Actions.end()) {
        PdAction &CurAction = *it;
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
        oscofo_output_descriptors(x, Desc);
    }
}

// ─────────────────────────────────────
static void oscofo_ticknewevent(PdOpenScofo *x) {
    int PrevStateIndex = x->StateIndex;
    x->Event = x->OpenScofo->GetCurrentScorePosition();
    x->StateIndex = x->OpenScofo->GetCurrentStateIndex();
    if (PrevStateIndex == x->StateIndex) {
        return;
    }

    outlet_float(x->TempoOut, x->OpenScofo->GetCurrentBPM());
    outlet_float(x->EventOut, x->OpenScofo->GetCurrentScorePosition());
    OpenScofo::EventActions Actions = x->OpenScofo->GetCurrentEventActions();

    for (OpenScofo::ScoreAction &Act : Actions) {
        double time = Act.Time;
        if (!Act.AbsoluteTime) {
            Act.Time = 60.0 / x->OpenScofo->GetCurrentBPM() * Act.Time * 1000;
            time = Act.Time;
        }

        if (time == 0) {
            if (Act.isLua) {
                oscofo_luaexecute(x, Act.Lua);
            } else {
                t_atom *PdArgs = oscofo_convertargs(Act);
                oscofo_pdsend(x, Act.Receiver, Act.Args.size(), PdArgs);
                delete[] PdArgs;
            }
        } else {
            double sysTime = clock_getsystimeafter(time);
            int size = Act.Args.size();
            std::string receiver = Act.Receiver;
            t_atom *PdArgs = oscofo_convertargs(Act);
            PdAction action = {sysTime, Act.isLua, receiver, Act.Lua, PdArgs, size};
            x->Actions.push_back(action);
        }
    }
}

// ─────────────────────────────────────
static t_int *oscofo_perform_score(t_int *w) {
    PdOpenScofo *x = (PdOpenScofo *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    int n = static_cast<int>(w[3]);
    if (!x->Following && !x->JustDescription) {
        return (w + 4);
    }

    bool ok = x->OpenScofo->ProcessBlock(in, n);
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
    dsp_add(oscofo_perform_score, 3, x, sp[0]->s_vec, sp[0]->s_n);
}

// ─────────────────────────────────────
static void oscofo_error_callback(const spdlog::details::log_msg &log, void *data) {
    PdOpenScofo *x = static_cast<PdOpenScofo *>(data);
    spdlog::level::level_enum pdlevel = x->log;
    if (log.level < pdlevel) {
        return;
    }

    std::string text(log.payload.data(), log.payload.size());
    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
        logpost(x, 1, "[openscofo~] %s", text.c_str());
        break;
    case spdlog::level::info:
    case spdlog::level::warn:
        logpost(x, 2, "[openscofo~] %s", text.c_str());
        break;
    case spdlog::level::debug:
    case spdlog::level::trace:
        logpost(x, 3, "[openscofo~] %s", text.c_str());
        break;
    default:
        break;
    }
}

// ─────────────────────────────────────
static void *oscofo_new(t_symbol *s, int argc, t_atom *argv) {
    PdOpenScofo *x = (PdOpenScofo *)pd_new(OpenScofoObj);
    (void)s;

    if (!x) {
        pd_error(x, "[openscofo~] Error creating object");
        return nullptr;
    }

    // default parameters
    x->FFTSize = 2048;
    x->HopSize = 512;
    x->Sr = (int)sys_getsr();
    x->Following = false;
    x->Event = -1;
    x->StateIndex = 0;

    // Outlets
    x->EventOut = outlet_new(&x->PdObject, &s_float);
    x->TempoOut = outlet_new(&x->PdObject, &s_float);

    // Args
    x->RequestMIR = oscofo_get_descriptors(x, 0, argc, argv);
    if (x->RequestMIR.size() > 0) {
        x->DescOut = outlet_new(&x->PdObject, &s_list);
        x->MirOutput = true;
        x->JustDescription = true;
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
    x->OpenScofo->SetErrorCallback(oscofo_error_callback, static_cast<void *>(x));

    x->log = (spdlog::level::warn);

#ifdef NDEBUG
#else
    // x->OpenScofo->SetLogLevel(spdlog::level::info);
#endif

#ifdef OPENSCOFO_LUA
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
extern "C" void openscofo_tilde_setup(void) {
    OpenScofoObj = class_new(gensym("openscofo~"), (t_newmethod)oscofo_new, (t_method)oscofo_free, sizeof(PdOpenScofo),
                             CLASS_DEFAULT, A_GIMME, A_NULL);

    post("[openscofo~] version %s (%s) (openscofo~ %s %s), by Charles K. Neimog\n\n", OPENSCOFO_VERSION,
         OPENSCOFO_BUILD_TIME, __DATE__, __TIME__);

    // message methods
    class_addmethod(OpenScofoObj, (t_method)oscofo_score, gensym("score"), A_SYMBOL, 0);
    class_addmethod(OpenScofoObj, (t_method)oscofo_start, gensym("start"), A_NULL, 0);
    class_addmethod(OpenScofoObj, (t_method)oscofo_following, gensym("follow"), A_FLOAT, 0);
    class_addmethod(OpenScofoObj, (t_method)oscofo_set, gensym("set"), A_GIMME, 0);
    class_addmethod(OpenScofoObj, (t_method)oscofo_get, gensym("get"), A_GIMME, 0);

    // dsp
    CLASS_MAINSIGNALIN(OpenScofoObj, PdOpenScofo, Sample);
    class_addmethod(OpenScofoObj, (t_method)oscofo_adddsp, gensym("dsp"), A_CANT, 0);
}
