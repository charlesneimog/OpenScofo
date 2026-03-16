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
static void oscofo_output_descriptiors(MaxOpenScofo *x, OpenScofo::Description &Desc) {
    for (auto it = x->RequestMIR.rbegin(); it != x->RequestMIR.rend(); ++it) {
        OpenScofo::Descriptors v = *it;
        if (v == OpenScofo::Descriptors::MFCC) {
            size_t mfccSize = Desc.MFCC.size();
            std::vector<t_atom> mfccAtoms(mfccSize);
            for (size_t i = 0; i < mfccSize; ++i) {
                atom_setfloat(&mfccAtoms[i], (float)Desc.MFCC[i]);
            }
            outlet_anything(x->DescOut, gensym("mfcc"), mfccSize, mfccAtoms.data());
        } else if (v == OpenScofo::Descriptors::CHROMA) {
            size_t chromaSize = Desc.Chroma.size();
            std::vector<t_atom> chromaAtoms(chromaSize);
            for (size_t i = 0; i < chromaSize; ++i) {
                atom_setfloat(&chromaAtoms[i], (float)Desc.Chroma[i]);
            }
            outlet_anything(x->DescOut, gensym("chroma"), chromaSize, chromaAtoms.data());
        } else if (v == OpenScofo::Descriptors::POWER) {
            size_t powerSize = Desc.Magnitude.size();
            std::vector<t_atom> powerAtoms(powerSize);
            for (size_t i = 0; i < powerSize; ++i) {
                atom_setfloat(&powerAtoms[i], (float)Desc.Magnitude[i]);
            }
            outlet_anything(x->DescOut, gensym("power"), powerSize, powerAtoms.data());
        } else if (v == OpenScofo::Descriptors::LOUDNESS) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.Loudness);
            outlet_anything(x->DescOut, gensym("loudness"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::SILENCEPROB) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.SilenceProb);
            outlet_anything(x->DescOut, gensym("silence"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::CENTROID) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.SpectralCentroid);
            outlet_anything(x->DescOut, gensym("centroid"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::SPREAD) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.SpectralSpread);
            outlet_anything(x->DescOut, gensym("spread"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::FLATNESS) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.SpectralFlatness);
            outlet_anything(x->DescOut, gensym("flatness"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::FLUX) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.SpectralFlux);
            outlet_anything(x->DescOut, gensym("flux"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::IRREGULARITY) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.SpectralIrregularity);
            outlet_anything(x->DescOut, gensym("irregularity"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::YIN) {
            std::vector<t_atom> atoms(2);
            atom_setfloat(&atoms[0], (float)Desc.Pitch);
            atom_setfloat(&atoms[1], (float)Desc.PitchConfidence);
            outlet_anything(x->DescOut, gensym("yin"), 2, atoms.data());
        } else if (v == OpenScofo::Descriptors::HARMONICITY) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.Harmonicity);
            outlet_anything(x->DescOut, gensym("harmonicity"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::PERCUSSIVEPROB) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.PercussiveTechProb);
            outlet_anything(x->DescOut, gensym("perc"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::EXTENDEDPROB) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.NoiseTechProb);
            outlet_anything(x->DescOut, gensym("ext"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::ONSET) {
            if (Desc.Onset) {
                std::vector<t_atom> atoms(1);
                atom_setfloat(&atoms[0], (float)Desc.Onset);
                outlet_anything(x->DescOut, gensym("onset"), 1, atoms.data());
            }
        } else if (v == OpenScofo::Descriptors::RMS) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.RMS);
            outlet_anything(x->DescOut, gensym("rms"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::CHROMA) {
        } else if (v == OpenScofo::Descriptors::HFR) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.HighFreqRatio);
            outlet_anything(x->DescOut, gensym("hfr"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::ZCR) {
            std::vector<t_atom> atoms(1);
            atom_setfloat(&atoms[0], (float)Desc.ZeroCrossingRate);
            outlet_anything(x->DescOut, gensym("zcr"), 1, atoms.data());
        } else if (v == OpenScofo::Descriptors::ONNX) {
            for (const auto &onnxDesc : Desc.ONNX) {
                std::vector<t_atom> onnxAtoms(2);
                atom_setsym(&onnxAtoms[0], gensym(onnxDesc.first.c_str()));
                atom_setfloat(&onnxAtoms[1], (float)onnxDesc.second);
                outlet_anything(x->DescOut, gensym("onnx"), 2, onnxAtoms.data());
            }
        }
    }
}

// ─────────────────────────────────────
static std::vector<OpenScofo::Descriptors> oscofo_get_descriptors(MaxOpenScofo *x, int start, long argc, t_atom *argv) {
    std::vector<OpenScofo::Descriptors> descriptors;
    for (long i = start; i < argc; i++) {
        if ((argv + i)->a_type != A_SYM) {
            object_error((t_object *)x, "Invalid argument type at index %ld", i);
            continue;
        }

        t_symbol *sym = atom_getsym(argv + i);
        if (strcmp(sym->s_name, "mfcc") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::MFCC);
        } else if (strcmp(sym->s_name, "rms") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::RMS);
        } else if (strcmp(sym->s_name, "loudness") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::LOUDNESS);
        } else if (strcmp(sym->s_name, "chroma") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::CHROMA);
        } else if (strcmp(sym->s_name, "silence") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::SILENCEPROB);
        } else if (strcmp(sym->s_name, "centroid") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::CENTROID);
        } else if (strcmp(sym->s_name, "zcr") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::ZCR);
        } else if (strcmp(sym->s_name, "hfr") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::HFR);
        } else if (strcmp(sym->s_name, "spread") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::SPREAD);
        } else if (strcmp(sym->s_name, "flatness") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::FLATNESS);
        } else if (strcmp(sym->s_name, "flux") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::FLUX);
        } else if (strcmp(sym->s_name, "irregularity") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::IRREGULARITY);
        } else if (strcmp(sym->s_name, "harmonicity") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::HARMONICITY);
        } else if (strcmp(sym->s_name, "ext") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::EXTENDEDPROB);
        } else if (strcmp(sym->s_name, "perc") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::PERCUSSIVEPROB);
        } else if (strcmp(sym->s_name, "onset") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::ONSET);
        } else if (strcmp(sym->s_name, "yin") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::YIN);
        } else if (strcmp(sym->s_name, "onnx") == 0) {
            descriptors.push_back(OpenScofo::Descriptors::ONNX);
        } else {
            object_error((t_object *)x, "Invalid argument: %s", sym->s_name);
        }
    }
    return descriptors;
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
        oscofo_output_descriptiors(x, Desc);
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
static void oscofo_tickactions(MaxOpenScofo *x) {
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
