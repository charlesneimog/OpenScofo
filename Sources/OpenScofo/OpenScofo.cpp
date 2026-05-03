#include <OpenScofo.hpp>

// ╭─────────────────────────────────────╮
// │     Construstor and Destructor      │
// ╰─────────────────────────────────────╯
namespace OpenScofo {

#if defined(OPENSCOFO_LUA)
int luaopen_OpenScofo(lua_State *L);
#endif

//  ─────────────────────────────────────
/**
 * @brief Initialize OpenScofo processing pipeline.
 *
 * @param Sr Sampling rate (Hz)
 * @param FftSize FFT window size
 * @param HopSize Hop size (samples)
 *
 * @note Initializes Forward model, MIR extractor, and score handler.
 * @note Configures global spdlog logger (overwrites default).
 */
OpenScofo::OpenScofo(float Sr, float FftSize, float HopSize)
    : m_Forward(Sr, FftSize, HopSize), m_MIR(Sr, FftSize, HopSize), m_Score(FftSize, HopSize) {

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

    // #if defined(NDEBUG)
    //     spdlog::set_level(spdlog::level::info);
    // #else
    spdlog::set_level(spdlog::level::debug);
    spdlog::enable_backtrace(32);
    // #endif

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
/**
 * @brief Update audio processing parameters.
 *
 * @param Sr Sampling rate (Hz)
 * @param FFTSize FFT window size
 * @param HopSize Hop size (samples)
 *
 * @note Propagates parameters to forward model, MIR extractor, and score.
 * @note Reallocates internal buffers and descriptor arrays if sizes change.
 * @note Resets input buffer and processing state.
 *
 * @warning Not thread-safe. Must not be called during audio processing.
 */
void OpenScofo::SetNewAudioParameters(float Sr, float FFTSize, float HopSize) {
    size_t NHalf = FFTSize / 2 + 1;
    if (m_FFTSize == FFTSize && m_HopSize == HopSize && m_Sr == Sr && m_Desc.Magnitude.size() == NHalf) {
        spdlog::debug("Everything allocated for FFTSize {}, NHalf {}", FFTSize, NHalf);
        return;
    }
    m_Sr = Sr;
    m_FFTSize = FFTSize;
    m_HopSize = HopSize;
    m_Forward.UpdateAudioParameters(Sr, FFTSize, HopSize);
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
/**
 * @brief Set callback for log/error messages.
 *
 * @param cb Callback invoked on each log message
 * @param data User-defined pointer passed to the callback
 *
 * @note The callback is triggered by the internal logging sink.
 * @note Updates the internal error flag (m_HasErrors) automatically.
 * @note Log level is set to debug (debug builds) or info (release builds).
 *
 * @warning Overwrites any previously registered callback.
 */
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
/**
 * @brief Set logging verbosity level.
 *
 * @param level spdlog log level
 *
 * @note Affects the global default spdlog logger.
 */
void OpenScofo::SetLogLevel(spdlog::level::level_enum level) {
    auto logger = spdlog::default_logger();
    spdlog::set_level(level);
}

// ─────────────────────────────────────
/**
 * @brief Reset internal error state.
 *
 * @note Clears m_HasErrors unless a critical error was previously set.
 *
 * @warning If a critical error occurred, the state is not reset and recovery
 * is not possible without reinitializing the instance.
 */
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
/**
 * @brief Load an ONNX model for descriptor inference.
 *
 * @param Model Path to .onnx model file
 * @param Descriptors List of descriptors expected by the model
 *
 * @note Only .onnx models are supported.
 * @note Delegates initialization to the MIR module.
 *
 * @warning Invalid file extension or incompatible descriptors will result in an error log.
 */
void OpenScofo::LoadONNXModel(fs::path Model, std::vector<Descriptors> Descriptors) {
    if (Model.extension() != ".onnx") {
        spdlog::error("OpenScofo just work with onnx models trained with the object py.train from pd-xlab library");
        return;
    }

    m_MIR.ONNXInit(Model, Descriptors);
}

// ╭─────────────────────────────────────╮
// │                 Lua                 │
// ╰─────────────────────────────────────╯
#if defined(OPENSCOFO_LUA)
/**
 * @brief Initialize embedded Lua runtime and OpenScofo bindings.
 *
 * @note Creates a new Lua state and opens standard libraries.
 * @note Exposes a global `_OpenScofo` table with a lightuserdata pointer to this instance.
 * @note Registers the OpenScofo Lua module via `luaL_requiref`.
 */
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
/**
 * @brief Register a Lua module into the current Lua state.
 *
 * @param name Module name exposed to Lua
 * @param func Lua C function used to initialize the module
 *
 * @return true if module was successfully registered, false otherwise
 *
 * @note Requires a valid Lua state (m_LuaState != nullptr).
 * @note Uses luaL_requiref, so the module may be cached by Lua.
 */
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
/**
 * @brief Execute a Lua code string in the current Lua state.
 *
 * @param code Lua source code to execute
 *
 * @return true if execution succeeded, false on load/runtime error
 */
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
/**
 * @brief Expose a raw pointer to Lua as a global lightuserdata.
 *
 * @param pointer C/C++ pointer to expose
 * @param name Global variable name in Lua
 *
 * @return true if Lua state is valid and pointer was set, false otherwise
 *
 * @note Stored as lightuserdata (no ownership or lifetime management).
 * @warning Lua code can access this pointer without safety checks.
 */
bool OpenScofo::LuaAddPointer(void *pointer, const char *name) {
    if (m_LuaState == nullptr) {
        return false;
    }
    lua_pushlightuserdata(m_LuaState, pointer);
    lua_setglobal(m_LuaState, name);
    return true;
}

// ─────────────────────────────────────
/**
 * @brief Extend Lua module search path.
 *
 * @param path Directory to add to package.path
 *
 * @note Appends a ".lua" search pattern for the given directory.
 * @note Modifies global Lua `package.path`.
 */
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
/**
 * @brief Retrieve and pop the last Lua error message.
 *
 * @return Error string if present, otherwise a default message
 *
 * @note Reads error from the top of the Lua stack and removes it.
 * @note If no valid string is found, returns a generic error message.
 */
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
/**
 * @brief Set pitch template variance parameter.
 *
 * @param Sigma Smoothing / spread parameter for pitch template
 *
 * @note Updates forward model and rebuilds its internal audio template.
 */
void OpenScofo::SetPitchTemplateSigma(double Sigma) {
    m_Forward.SetPitchTemplateSigma(Sigma);
    m_Forward.UpdateAudioTemplate();
}

// ─────────────────────────────────────
/**
 * @brief Set amplitude decay factor for the forward model.
 *
 * @param decay Decay coefficient applied to amplitude tracking
 *
 * @note Directly updates forward model parameter (no recomputation triggered here).
 */
void OpenScofo::SetAmplitudeDecay(double decay) {
    m_Forward.SetAmplitudeDecay(decay);
}

// ─────────────────────────────────────
/**
 * @brief Set number of harmonics used in pitch template (default: 10).
 *
 * @param Harmonics Number of harmonic components to model
 *
 * @note Updates forward model and rebuilds its audio template.
 */
void OpenScofo::SetHarmonics(int Harmonics) {
    m_Forward.SetHarmonics(Harmonics);
    m_Forward.UpdateAudioTemplate();
}

// ─────────────────────────────────────
/**
 * @brief Set amplitude threshold in dB (default: -60).
 *
 * @param dB Threshold value in decibels
 *
 * @note Applies to both forward model and MIR feature extraction.
 */
void OpenScofo::SetdBTreshold(double dB) {
    m_MIR.SetdBTreshold(dB);
}

// ─────────────────────────────────────
/**
 * @brief Set reference tuning frequency for A4 (default: 440).
 *
 * @param Tunning Base tuning reference (e.g., 440 Hz standard A4)
 *
 * @note Affects pitch-related interpretation in the scoring system.
 */
void OpenScofo::SetTunning(double Tunning) {
    m_Score.SetTunning(Tunning);
}

// ─────────────────────────────────────
/**
 * @brief Set active score event and reset decoding state.
 *
 * @param Event Event index in the loaded score (0 = reset)
 *
 * @note Resets forward model state, buffers, and descriptors.
 * @note Updates current score position based on event mapping if valid.
 */
void OpenScofo::SetCurrentEvent(int Event) {
    m_CurrentScorePosition = 0;
    m_BlockIndex = 0;
    m_Forward.ResetDecoding();
    m_Forward.SetCurrentEvent(Event);
    std::fill(m_InBuffer.begin(), m_InBuffer.end(), 0.0);

    std::fill(m_Desc.Magnitude.begin(), m_Desc.Magnitude.end(), 0.0);
    std::fill(m_Desc.Power.begin(), m_Desc.Power.end(), 0.0);
    std::fill(m_Desc.SpectralMagnitudeNorm.begin(), m_Desc.SpectralMagnitudeNorm.end(), 0.0);
    std::fill(m_Desc.SpectralMagnitudeFrameNorm.begin(), m_Desc.SpectralMagnitudeFrameNorm.end(), 0.0);
    std::fill(m_Desc.ReverbSpectralPower.begin(), m_Desc.ReverbSpectralPower.end(), 0.0);

    if (Event == 0) {
        m_CurrentScorePosition = 0;
        return;
    }

    if (Event > 0 && static_cast<size_t>(Event) < m_States.size()) {
        m_CurrentScorePosition = m_States[Event].ScorePos;
        return;
    }
}

// ╭─────────────────────────────────────╮
// │            Get Functions            │
// ╰─────────────────────────────────────╯
/**
 * @brief Get current position in the score (following Antescofo, Rest does not count for this).
 *
 * @return Score position index computed by the forward model.
 */
int OpenScofo::GetCurrentScorePosition() {
    return m_CurrentScorePosition;
}

// ─────────────────────────────────────
/**
 * @brief Get current state index from the forward model.
 *
 * @return Index of the active internal state.
 */
int OpenScofo::GetCurrentStateIndex() {
    return m_Forward.GetCurrentStateIndex();
}

// ─────────────────────────────────────
/**
 * @brief Get estimated current tempo.
 *
 * @return Current BPM estimate from the forward model.
 */
double OpenScofo::GetCurrentBPM() {
    return m_Forward.GetCurrentBPM();
}

// ─────────────────────────────────────
/**
 * @brief Get actions associated with the current score event.
 *
 * @return Event action list (empty if no score is loaded).
 *
 * @note Delegates to the forward model when a score is active.
 */
EventActions OpenScofo::GetCurrentEventActions() {
    if (ScoreIsLoaded()) {
        return m_Forward.GetCurrentEventActions();
    } else {
        return {};
    }
}

// ─────────────────────────────────────
/**
 * @brief Compute pitch probability for a given frequency.
 *
 * @param Freq Frequency in Hz
 *
 * @return Probability score from forward model
 *
 * @note Uses current description frame as input to the model.
 */
double OpenScofo::GetPitchProb(double Freq) {
    m_Forward.SetDescription(m_Desc);
    return m_Forward.GetPitchProbability(Freq);
}

// ─────────────────────────────────────
/**
 * @brief Return Lua code string defined in global events using LUA {}.
 *
 * @return Lua script as a string
 */
std::string OpenScofo::GetLuaCode() {
    return m_Score.GetLuaCode();
}

// ╭─────────────────────────────────────╮
// │          Helpers Functions          │
// ╰─────────────────────────────────────╯
/**
 * @brief Check if a score is currently loaded.
 *
 * @return true if score data is available, false otherwise
 */
bool OpenScofo::ScoreIsLoaded() {
    return m_Score.ScoreIsLoaded();
}

// ─────────────────────────────────────
/**
 * @brief Convert descriptor enum to string identifier.
 *
 * @param d Descriptor enum value
 *
 * @return Human-readable identifier string (e.g. "mfcc", "chroma")
 */
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
    case Descriptors::EXTENDEDTECHNIQUE:
        return "ext";
    case Descriptors::ODSONSET:
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
/**
 * @brief Convert string identifier to descriptor enum.
 *
 * @param s Descriptor name (e.g. "mfcc", "chroma")
 *
 * @return Corresponding Descriptors enum value, or INVALID on failure
 *
 * @note Logs an error if the string is not recognized.
 */
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
        return Descriptors::EXTENDEDTECHNIQUE;
    } else if (strcmp(s, "onset") == 0) {
        return Descriptors::ODSONSET;
    } else if (strcmp(s, "irregularity") == 0) {
        return Descriptors::IRREGULARITY;
    } else if (strcmp(s, "kurtosis") == 0) {
        return Descriptors::KURTOSIS;
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
/**
 * @brief Extract scalar descriptor value from a Description.
 *
 * @param Desc Audio description container
 * @param d Descriptor type
 *
 * @return Scalar value for the requested descriptor, or -1.0 on error / invalid type
 *
 * @note Logs an error if the descriptor is vector-valued or invalid.
 */
double OpenScofo::GetDescriptionFloat(Description &Desc, Descriptors d) {
    switch (d) {
    // Scalar descriptors
    case Descriptors::ODSONSET:
        return Desc.Onset;
    case Descriptors::SILENCEPROB:
        return Desc.SilenceProb;
    case Descriptors::EXTENDEDTECHNIQUE:
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
/**
 * @brief Access vector-valued descriptor data.
 *
 * @param Desc Audio description container
 * @param d Descriptor type (must be vector-valued)
 *
 * @return Reference to internal descriptor array
 *
 * @throws std::runtime_error if descriptor is not vector-valued
 *
 * @note Logs an error before throwing for invalid descriptor types.
 */
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
        spdlog::critical("Descriptor '{}' is not an array/vector type, returning Magnitude", GetDescriptionId(d));
        return Desc.Magnitude;
    }
}

// ╭─────────────────────────────────────╮
// │ Python Research and Test Functions  │
// ╰─────────────────────────────────────╯
/**
 * @brief Access internal score state machine.
 *
 * @return Reference to forward model states container.
 *
 * @note Exposes internal mutable state (no copy is made).
 */
States &OpenScofo::GetStates() {
    return m_Forward.GetStates();
}

// ─────────────────────────────────────
/**
 * @brief Generate pitch template for a given frequency.
 *
 * @param Freq Target frequency in Hz
 *
 * @return Pitch template vector computed by the forward model
 */
std::vector<double> OpenScofo::GetPitchTemplate(double Freq) {
    return m_Forward.GetPitchTemplate(Freq);
}

// ─────────────────────────────────────
/**
 * @brief Get current sampling rate.
 *
 * @return Sampling rate in Hz
 */
double OpenScofo::GetSr() {
    return m_Sr;
}

// ─────────────────────────────────────
/**
 * @brief Get FFT window size.
 *
 * @return FFT size in samples
 */
double OpenScofo::GetFFTSize() {
    return m_FFTSize;
}

// ─────────────────────────────────────
/**
 * @brief Get hop size used for frame processing.
 *
 * @return Hop size in samples
 */
double OpenScofo::GetHopSize() {
    return m_HopSize;
}

// ─────────────────────────────────────
/**
 * @brief Get processing block duration in seconds.
 *
 * @return Block duration in seconds (derived from hop size and sampling rate)
 */
double OpenScofo::GetBlockDuration() {
    return m_Forward.GetBlockDuration();
}

// ╭─────────────────────────────────────╮
// │           Main Functions            │
// ╰─────────────────────────────────────╯
/**
 * @brief Get current audio description frame.
 *
 * @return Copy of the current descriptor structure
 *
 * @note Returns by value (snapshot, not live reference).
 */
Description OpenScofo::GetDescription() {
    return m_Desc;
}

// ─────────────────────────────────────
/**
 * @brief Get current processing buffer index.
 *
 * @return Current index within the analysis buffer (forward model state)
 */
int OpenScofo::GetCurrentBufferIndex() {
    return m_Forward.GetCurrentBufferIndex();
}

// ─────────────────────────────────────
/**
 * @brief Load and initialize a score from file.
 *
 * @param ScorePath Path to score file
 *
 * @return true if loading and initialization succeeded, false otherwise
 *
 * @note Resets internal state and reinitializes processing pipeline.
 * @note May load and validate an ONNX model if present in the score.
 */
bool OpenScofo::LoadScore(fs::path ScorePath) {
    ClearErrors();

    m_CurrentScorePosition = 0;

    m_Score.UpdateAudioParameters(m_FFTSize, m_HopSize);
    m_States.clear();
    m_States = m_Score.Parse(ScorePath);

    // Timbre/Extended Tech detection
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
    m_Forward.SetPitchTemplateSigma(m_Score.GetPitchTemplateSigma());

    // Add States
    m_Forward.SetScoreStates(m_States);
    m_Mode = SCOREFOLLOWER;

    // verify the states.
    const std::vector<std::string> &ONNXLabels = m_MIR.GetONNXLabels();
    for (auto &state : m_States) {
        for (auto &audioState : state.AudioStates) {
            if (audioState.Type == LABEL) {
                const auto &label = audioState.Label;
                auto it = std::find(ONNXLabels.begin(), ONNXLabels.end(), label);
                if (it == ONNXLabels.end()) {
                    std::string validLabels = "[";
                    for (size_t i = 0; i < ONNXLabels.size(); ++i) {
                        validLabels += ONNXLabels[i];
                        if (i + 1 < ONNXLabels.size()) {
                            validLabels += ", ";
                        }
                    }
                    validLabels += "]";
                    spdlog::error("Extended Technique Label '{}' is not valid. Valid labels: {}", label, validLabels);
                    return false;
                }
            }
        }
    }

    if (m_HasErrors != spdlog::level::err && m_HasErrors != spdlog::level::critical) {
        return true;
    } else {
        return false;
    }
}

// ─────────────────────────────────────
/**
 * @brief Process an incoming audio block.
 *
 * @tparam T Audio sample precision type (float/double)
 * @param AudioBuffer Input audio buffer
 * @param n Number of samples in buffer
 *
 * @return true if processing succeeded
 *
 * @note Maintains internal circular buffer state.
 * @note Triggers analysis every hop size.
 * @note Updates descriptors and score position depending on mode.
 */
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
        m_CurrentScorePosition = m_Forward.GetEvent(m_Desc);
        m_MIR.AddReverb(m_Desc, 0.01);
        break;

    case DESCRIPTORS:
        m_MIR.GetDescription(m_InBuffer, m_Desc);
        m_Forward.SetDescription(m_Desc);
        break;
    }

    return true;
}

// ─────────────────────────────────────
template bool OpenScofo::ProcessBlock<float>(const float *, size_t);
template bool OpenScofo::ProcessBlock<double>(const double *, size_t);

} // namespace OpenScofo
