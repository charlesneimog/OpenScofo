#include <OpenScofo.hpp>

// ╭─────────────────────────────────────╮
// │     Construstor and Destructor      │
// ╰─────────────────────────────────────╯
namespace OpenScofo {

#if defined(OPENSCOFO_LUA)
int luaopen_OpenScofo(lua_State *L);
#endif

//  ─────────────────────────────────────
OpenScofo::OpenScofo(float Sr, float FftSize, float HopSize)
    : m_MDP(Sr, FftSize, HopSize), m_MIR(Sr, FftSize, HopSize), m_Score(FftSize, HopSize) {

    m_States = States();
    m_Desc = Description();
    m_Sr = Sr;
    m_FFTSize = FftSize;
    m_HopSize = HopSize;
    m_InBuffer.reserve(FftSize);
    m_BlockIndex = 0;

#if defined(OPENSCOFO_LUA)
    InitLuaModule();
#endif

#if defined(NDEBUG)
    spdlog::set_level(spdlog::level::info);
#else
    spdlog::set_level(spdlog::level::debug);
    spdlog::enable_backtrace(32);
#endif

    // --- Create OpenScofoLog sink ---
    m_Log = std::make_shared<OpenScofoLog<std::mutex>>();
    m_Log->SetCallback(nullptr, nullptr, &m_HasErrors); // ensures error flag updates
    std::vector<spdlog::sink_ptr> sinks{m_Log};

#ifndef NDEBUG
    auto consoleSink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    sinks.push_back(consoleSink);
    m_Log->set_pattern("%v"); // keep only message for callback sink
#endif

    auto logger = std::make_shared<spdlog::logger>("OpenScofo", sinks.begin(), sinks.end());
    spdlog::set_default_logger(logger);
    SetNewAudioParameters(Sr, FftSize, HopSize);
}

//  ─────────────────────────────────────
void OpenScofo::SetNewAudioParameters(float Sr, float FFTSize, float HopSize) {
    size_t NHalf = FFTSize / 2 + 1;
    if (m_FFTSize == FFTSize && m_HopSize == HopSize && m_Sr == Sr && m_Desc.Magnitude.size() == NHalf) {
        spdlog::debug("Everything allocated for FFTSize {}, NHalf {}", FFTSize, NHalf);
        return;
    }
    m_Sr = Sr;
    m_FFTSize = FFTSize;
    m_HopSize = HopSize;
    m_MDP.UpdateAudioParameters(Sr, FFTSize, HopSize);
    m_MIR.UpdateAudioParameters(Sr, FFTSize, HopSize);
    m_Score.UpdateAudioParameters(FFTSize, HopSize);
    m_InBuffer.resize(FFTSize);
    std::fill(m_InBuffer.begin(), m_InBuffer.end(), 0.0);
    m_BlockIndex = 0;

    spdlog::debug("Allocated Description size for Window Size {}, NHalf {}", FFTSize, NHalf);

    if (NHalf != m_Desc.Magnitude.size()) {
        m_Desc.Magnitude.resize(NHalf);
        m_Desc.Power.resize(NHalf);
        m_Desc.SpectralMagnitudeNorm.resize(NHalf);
        m_Desc.SpectralMagnitudeFrameNorm.resize(NHalf);
        m_Desc.ReverbSpectralPower.resize(NHalf);
    }
}

// ╭─────────────────────────────────────╮
// │               Errors                │
// ╰─────────────────────────────────────╯
void OpenScofo::SetErrorCallback(std::function<void(const spdlog::details::log_msg &, void *data)> cb, void *data) {
    if (m_Log) {
        m_Log->SetCallback(cb, data, &m_HasErrors);
#if defined(NDEBUG)
        spdlog::set_level(spdlog::level::info);
#else
        spdlog::set_level(spdlog::level::debug);
#endif
    } else {
        std::cerr << "Not possible to create Log" << std::endl;
    }
}

// ─────────────────────────────────────
void OpenScofo::SetLogLevel(spdlog::level::level_enum level) {
    auto logger = spdlog::default_logger();
    spdlog::set_level(level);
}

// ─────────────────────────────────────
void OpenScofo::ClearErrors() {
    if (m_HasErrors == spdlog::level::critical) {
        spdlog::error(
            "Critical error encountered. Recovery is not possible. Please restart OpenScofo or report the error.");
        return;
    } else {
        m_HasErrors = spdlog::level::off;
    }
}

// ╭─────────────────────────────────────╮
// │                ONNX                 │
// ╰─────────────────────────────────────╯
void OpenScofo::LoadONNXModel(fs::path Model, std::vector<Descriptors> Descriptors) {
    if (Model.extension() != ".onnx") {
        spdlog::error("OpenScofo just work with onnx models trained with the object py.train from pd-xlab library");
        return;
    }

    // m_MIR.ONNXInit(Model, Descriptors);
    std::thread([ModelCopy = std::move(Model), DescriptorsCopy = Descriptors, this]() mutable {
        m_LoadingONNX = true;
        spdlog::warn("Loading ONNX model, wait...");
        m_MIR.ONNXInit(ModelCopy, DescriptorsCopy);
        m_LoadingONNX = false;
        spdlog::warn("ONNX Model ready");
    }).detach();
}

// ╭─────────────────────────────────────╮
// │                 Lua                 │
// ╰─────────────────────────────────────╯
#if defined(OPENSCOFO_LUA)
void OpenScofo::InitLuaModule() {
    m_LuaState = luaL_newstate();
    luaL_openlibs(m_LuaState); // NOTE: Rethink if I load all functions
    lua_newtable(m_LuaState);
    lua_pushlightuserdata(m_LuaState, this);
    lua_setfield(m_LuaState, -2, "pointer");
    lua_setglobal(m_LuaState, "_OpenScofo");
    luaL_requiref(m_LuaState, "OpenScofo", luaopen_OpenScofo, 1);
}

// ─────────────────────────────────────
bool OpenScofo::LuaAddModule(std::string name, lua_CFunction func) {
    if (m_LuaState == nullptr) {
        return false;
    }
    luaL_requiref(m_LuaState, name.c_str(), func, 1);
    if (lua_isnil(m_LuaState, -1)) {
        return false;
    }
    return true;
}

// ─────────────────────────────────────
bool OpenScofo::LuaExecute(std::string code) {
    if (m_LuaState == nullptr) {
        return false;
    }
    int status = luaL_loadstring(m_LuaState, code.c_str());
    if (status == LUA_OK) {
        status = lua_pcall(m_LuaState, 0, LUA_MULTRET, 0);
        if (status != LUA_OK) {
            return false;
        }
        return true;
    } else {
        return false;
    }
}

// ─────────────────────────────────────
bool OpenScofo::LuaAddPointer(void *pointer, const char *name) {
    if (m_LuaState == nullptr) {
        return false;
    }
    lua_pushlightuserdata(m_LuaState, pointer);
    lua_setglobal(m_LuaState, name);
    return true;
}

// ─────────────────────────────────────
void OpenScofo::LuaAddPath(std::string path) {
    if (m_LuaState == nullptr) {
        return;
    }

    lua_getglobal(m_LuaState, "package");
    lua_getfield(m_LuaState, -1, "path");
    const char *current_path = lua_tostring(m_LuaState, -1);
    if (path.back() != '/') {
        lua_pushfstring(m_LuaState, "%s;%s/?.lua", current_path, path.c_str());
    } else {
        lua_pushfstring(m_LuaState, "%s;%s?.lua", current_path, path.c_str());
    }

    lua_setfield(m_LuaState, -3, "path");
    lua_pop(m_LuaState, 1);
}

// ─────────────────────────────────────
std::string OpenScofo::LuaGetError() {
    if (m_LuaState == nullptr) {
        return "m_LuaState is null";
    }
    if (lua_isstring(m_LuaState, -1)) {
        std::string errorMsg = lua_tostring(m_LuaState, -1);
        lua_pop(m_LuaState, 1);
        return errorMsg;
    }
    return "Unknown error";
}
#endif

// ╭─────────────────────────────────────╮
// │            Set Functions            │
// ╰─────────────────────────────────────╯
void OpenScofo::SetPitchTemplateSigma(double Sigma) {
    m_MDP.SetPitchTemplateSigma(Sigma);
    m_MDP.UpdateAudioTemplate();
}

// ─────────────────────────────────────
void OpenScofo::SetAmplitudeDecay(double decay) {
    m_MDP.SetAmplitudeDecay(decay);
}

// ─────────────────────────────────────
void OpenScofo::SetHarmonics(int Harmonics) {
    m_MDP.SetHarmonics(Harmonics);
    m_MDP.UpdateAudioTemplate();
}

// ─────────────────────────────────────
double OpenScofo::GetdBValue() {
    return 0;
}

// ─────────────────────────────────────
void OpenScofo::SetdBTreshold(double dB) {
    m_MDP.SetdBTreshold(dB);
    m_MIR.SetdBTreshold(dB);
}

// ─────────────────────────────────────
void OpenScofo::SetTunning(double Tunning) {
    m_Score.SetTunning(Tunning);
}

// ─────────────────────────────────────
void OpenScofo::SetCurrentEvent(int Event) {
    m_MDP.SetCurrentEvent(Event);

    if (Event == 0) {
        m_CurrentScorePosition = 0;
        return;
    }

    if (Event > 0 && static_cast<size_t>(Event) < m_States.size()) {
        m_CurrentScorePosition = m_States[Event].ScorePos;
        return;
    }

    m_CurrentScorePosition = Event;
}

// ╭─────────────────────────────────────╮
// │            Get Functions            │
// ╰─────────────────────────────────────╯
int OpenScofo::GetCurrentScorePosition() {
    return m_CurrentScorePosition;
}

// ─────────────────────────────────────
double OpenScofo::GetLiveBPM() {
    return m_MDP.GetLiveBPM();
}

// ─────────────────────────────────────
EventActions OpenScofo::GetEventActions(int Index) {
    return m_MDP.GetEventActions(Index);
}

// ─────────────────────────────────────
double OpenScofo::GetKappa() {
    return m_MDP.GetKappa();
}

// ─────────────────────────────────────
double OpenScofo::GetPitchProb(double f) {
    return m_MDP.GetPitchProbability(f);
}

// ─────────────────────────────────────
std::string OpenScofo::GetLuaCode() {
    return m_Score.GetLuaCode();
}

// ╭─────────────────────────────────────╮
// │          Helpers Functions          │
// ╰─────────────────────────────────────╯
bool OpenScofo::ScoreIsLoaded() {
    return m_Score.ScoreIsLoaded();
}

// ─────────────────────────────────────
const char *OpenScofo::GetDescriptionId(Descriptors d) {
    switch (d) {
    case Descriptors::MFCC:
        return "mfcc";
    case Descriptors::MELOGRAM:
        return "logmel";
    case Descriptors::RMS:
        return "rms";
    case Descriptors::LOUDNESS:
        return "loudness";
    case Descriptors::STDDEV:
        return "stddev";
    case Descriptors::CHROMA:
        return "chroma";
    case Descriptors::SILENCEPROB:
        return "silence";
    case Descriptors::CENTROID:
        return "centroid";
    case Descriptors::ZCR:
        return "zcr";
    case Descriptors::HFR:
        return "hfr";
    case Descriptors::SPREADHZ:
        return "spread";
    case Descriptors::FLATNESS:
        return "flatness";
    case Descriptors::ENTROPY:
        return "entropy";
    case Descriptors::ROLLOFF:
        return "rolloff";
    case Descriptors::FLUX:
        return "flux";
    case Descriptors::SKEWNESS:
        return "skewness";
    case Descriptors::SLOPE:
        return "slope";
    case Descriptors::KURTOSIS:
        return "kurtosis";
    case Descriptors::EXTENDEDPROB:
        return "ext";
    case Descriptors::ONSET:
        return "onset";
    case Descriptors::YIN:
        return "yin";
    case Descriptors::ONNX:
        return "onnx";
    default:
        return "unknown";
    }
}

// ─────────────────────────────────────
Descriptors OpenScofo::GetDescriptorsEnum(const char *s) {
    if (strcmp(s, "mfcc") == 0) {
        return Descriptors::MFCC;
    } else if (strcmp(s, "logmel") == 0) {
        return Descriptors::MELOGRAM;
    } else if (strcmp(s, "rms") == 0) {
        return Descriptors::RMS;
    } else if (strcmp(s, "loudness") == 0) {
        return Descriptors::LOUDNESS;
    } else if (strcmp(s, "stddev") == 0) {
        return Descriptors::STDDEV;
    } else if (strcmp(s, "chroma") == 0) {
        return Descriptors::CHROMA;
    } else if (strcmp(s, "silence") == 0) {
        return Descriptors::SILENCEPROB;
    } else if (strcmp(s, "harmonicity") == 0) {
        return Descriptors::HARMONICITY;
    } else if (strcmp(s, "centroid") == 0) {
        return Descriptors::CENTROID;
    } else if (strcmp(s, "zcr") == 0) {
        return Descriptors::ZCR;
    } else if (strcmp(s, "hfr") == 0) {
        return Descriptors::HFR;
    } else if (strcmp(s, "spread") == 0) {
        return Descriptors::SPREADHZ;
    } else if (strcmp(s, "flatness") == 0) {
        return Descriptors::FLATNESS;
    } else if (strcmp(s, "entropy") == 0) {
        return Descriptors::ENTROPY;
    } else if (strcmp(s, "rolloff") == 0) {
        return Descriptors::ROLLOFF;
    } else if (strcmp(s, "flux") == 0) {
        return Descriptors::FLUX;
    } else if (strcmp(s, "skewness") == 0) {
        return Descriptors::SKEWNESS;
    } else if (strcmp(s, "slope") == 0) {
        return Descriptors::SLOPE;
    } else if (strcmp(s, "kurtosis") == 0) {
        return Descriptors::KURTOSIS;
    } else if (strcmp(s, "ext") == 0) {
        return Descriptors::EXTENDEDPROB;
    } else if (strcmp(s, "onset") == 0) {
        return Descriptors::ONSET;
    } else if (strcmp(s, "yin") == 0) {
        return Descriptors::YIN;
    } else if (strcmp(s, "onnx") == 0) {
        return Descriptors::ONNX;
    } else {
        spdlog::error("Invalid descriptors argument: {}", s);
        return Descriptors::INVALID;
    }
}

// ─────────────────────────────────────
double OpenScofo::GetDescriptionFloat(Description &Desc, Descriptors d) {
    switch (d) {
    // Scalar descriptors
    case Descriptors::ONSET:
        return Desc.Onset;
    case Descriptors::SILENCEPROB:
        return Desc.SilenceProb;
    case Descriptors::EXTENDEDPROB:
        return Desc.ExtendedTechProb;
    case Descriptors::DB:
        return Desc.dB;
    case Descriptors::RMS:
        return Desc.RMS;
    case Descriptors::MAXAMP:
        return Desc.MaxAmp;
    case Descriptors::LOUDNESS:
        return Desc.Loudness;
    case Descriptors::HARMONICITY:
        return Desc.Harmonicity;
    case Descriptors::FLATNESS:
        return Desc.SpectralFlatness;
    case Descriptors::ENTROPY:
        return Desc.SpectralEntropy;
    case Descriptors::ROLLOFF:
        return Desc.SpectralRolloff;
    case Descriptors::FLUX:
        return Desc.SpectralFlux;
    case Descriptors::IRREGULARITY:
        return Desc.SpectralIrregularity;
    case Descriptors::CREST:
        return Desc.SpectralCrest;
    case Descriptors::CENTROID:
        return Desc.SpectralCentroid;
    case Descriptors::CENTROIDVEL:
        return Desc.CentroidVelocity;
    case Descriptors::SPREADHZ:
        return Desc.SpectralSpreadHz;
    case Descriptors::SPREADVARIANCE:
        return Desc.SpectralSpreadVariance;
    case Descriptors::SKEWNESS:
        return Desc.SpectralSkewness;
    case Descriptors::SLOPE:
        return Desc.SpectralSlope;
    case Descriptors::KURTOSIS:
        return Desc.SpectralKurtosis;
    case Descriptors::HFR:
        return Desc.HighFreqRatio;
    case Descriptors::ZCR:
        return Desc.ZeroCrossingRate;
    case Descriptors::STDDEV:
        return Desc.StdDev;
    case Descriptors::YIN:
        return Desc.Pitch;
    case Descriptors::YINCONFIDENCE:
        return Desc.PitchConfidence;

    // Vector descriptors: cannot return single value
    case Descriptors::MFCC:
    case Descriptors::MELOGRAM:
    case Descriptors::MAGNITUDE:
    case Descriptors::POWERARRAY:
    case Descriptors::CHROMA:
    case Descriptors::ONNX:
        spdlog::error("Descriptor '{}' is vector-valued; cannot return a single double", GetDescriptionId(d));
        return -1.0;

    default:
        spdlog::error("Invalid descriptor '{}'", static_cast<int>(d));
        return -1.0;
    }
}

// ─────────────────────────────────────
std::vector<double> &OpenScofo::GetDescriptionArray(Description &Desc, Descriptors d) {
    switch (d) {
    case Descriptors::MFCC:
        return Desc.MFCC;
    case Descriptors::CHROMA:
        return Desc.Chroma;
    case Descriptors::MELOGRAM:
        return Desc.LogMelSpectrum;
    case Descriptors::POWERARRAY:
        return Desc.Power;
    case Descriptors::MAGNITUDE:
        return Desc.Magnitude;
    default:
        spdlog::error("Descriptor '{}' is not an array/vector type", GetDescriptionId(d));
        throw std::runtime_error("Descriptor is not a vector");
    }
}

// ╭─────────────────────────────────────╮
// │ Python Research and Test Functions  │
// ╰─────────────────────────────────────╯
States &OpenScofo::GetStates() {
    return m_MDP.GetStates();
}

// ─────────────────────────────────────
std::vector<double> OpenScofo::GetPitchTemplate(double Freq) {
    return m_MDP.GetPitchTemplate(Freq);
}

// ─────────────────────────────────────
std::vector<double> OpenScofo::GetSpectrumPower() const {
    return m_Desc.SpectralMagnitudeFrameNorm;
}

// ─────────────────────────────────────
double OpenScofo::GetSr() {
    return m_Sr;
}

// ─────────────────────────────────────
double OpenScofo::GetFFTSize() {
    return m_FFTSize;
}

// ─────────────────────────────────────
double OpenScofo::GetHopSize() {
    return m_HopSize;
}

// ─────────────────────────────────────
double OpenScofo::GetBlockDuration() {
    return m_MDP.GetBlockDuration();
}

// ╭─────────────────────────────────────╮
// │           Main Functions            │
// ╰─────────────────────────────────────╯
bool OpenScofo::ParseScore(fs::path ScorePath) {
    ClearErrors();

    m_CurrentScorePosition = 0;

    m_Score.UpdateAudioParameters(m_FFTSize, m_HopSize);
    m_States.clear();
    m_States = m_Score.Parse(ScorePath);

    // Timbre detection
    if (m_Score.HasTimbreModel()) {
        auto model = m_Score.GetTimbreModel();
        std::vector<std::string> descriptors = m_Score.GetTimbreModelDescriptors();
        spdlog::warn("Loading ONNX model, wait...");
        std::vector<Descriptors> DescEnum;
        for (auto d : descriptors) {
            DescEnum.push_back(GetDescriptorsEnum(d.c_str()));
        }
        m_MIR.ONNXInit(model, DescEnum);
        spdlog::warn("ONNX Model ready");
    }

    m_FFTSize = m_Score.GetFFTSize();
    m_HopSize = m_Score.GetHopSize();
    SetNewAudioParameters(m_Sr, m_FFTSize, m_HopSize);

    // Parse Config
    m_MDP.SetPitchTemplateSigma(m_Score.GetPitchTemplateSigma());

    // Add States
    m_MDP.SetScoreStates(m_States);
    m_Mode = SCOREFOLLOWER;

    if (m_HasErrors != spdlog::level::err && m_HasErrors != spdlog::level::critical) {
        return true;
    } else {
        return false;
    }
}

// ─────────────────────────────────────
Description OpenScofo::GetDescription() {
    return m_Desc;
}

// ─────────────────────────────────────
int OpenScofo::GetCurrentBufferIndex() {
    return m_MDP.GetCurrentBufferIndex();
}

// ─────────────────────────────────────
template <OpenScofoPrecision T> bool OpenScofo::ProcessBlock(const T *AudioBuffer, size_t n) {
    m_BlockIndex += n;

    std::copy(m_InBuffer.begin() + n, m_InBuffer.end(), m_InBuffer.begin());
    std::transform(AudioBuffer, AudioBuffer + n, m_InBuffer.end() - n, [](T x) { return static_cast<double>(x); });

    if (m_BlockIndex != m_HopSize) {
        return true;
    }

    m_BlockIndex = 0;

    switch (m_Mode) {
    case SCOREFOLLOWER:
        m_MIR.GetDescription(m_InBuffer, m_Desc);
        m_CurrentScorePosition = m_MDP.GetEvent(m_Desc);
        m_MIR.AddReverb(m_Desc, 0.01);
        break;

    case DESCRIPTORS:
        m_MIR.GetDescription(m_InBuffer, m_Desc);
        break;
    }

    return true;
}

// ─────────────────────────────────────
template bool OpenScofo::ProcessBlock<float>(const float *, size_t);
template bool OpenScofo::ProcessBlock<double>(const double *, size_t);

} // namespace OpenScofo
