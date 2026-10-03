/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file OpenScofo.cpp
 * @brief Public score-following API and audio processing pipeline.
 *
 * @note Coordinates score parsing, descriptor extraction, forward inference, and optional Lua execution.
 * @warning Serialize configuration changes, score loading, and audio processing on each instance.
 */

#include <OpenScofo.hpp>
#include <algorithm>

// ╭─────────────────────────────────────╮
// │     Construstor and Destructor      │
// ╰─────────────────────────────────────╯
namespace OpenScofo {

#if defined(OPENSCOFO_LUA)
int luaopen_OpenScofo(lua_State *L);
#endif

//  ─────────────────────────────────────
/**
 * @brief Initialize the audio analysis and score-following pipeline.
 *
 * @param Sr Sampling rate in Hz.
 * @param FftSize Analysis window size in samples.
 * @param HopSize Analysis hop size in samples.
 *
 * @note Initializes optional Lua bindings, installs a logging sink, and applies the initial audio configuration.
 * @warning Replaces the global default spdlog logger. Use a positive sample rate and supported FFT and hop sizes.
 */
OpenScofo::OpenScofo(float Sr, float FftSize, float HopSize) : m_Forward(), m_MIR(), m_Score() {
    m_States = States();
    m_Desc = Description();
    m_Config = Configuration();
    m_Config.SR = Sr;
    m_Config.FFTSize = FftSize;
    m_Config.HOPSize = HopSize;

    m_InBuffer.reserve(m_Config.FFTSize);
    m_BlockIndex = 0;

#if defined(OPENSCOFO_LUA)
    InitLuaModule();
#endif

    spdlog::set_level(spdlog::level::debug);
    spdlog::enable_backtrace(32);

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
    UpdateConfiguration(m_Config);
}

//  ─────────────────────────────────────
/**
 * @brief Apply audio configuration and resize descriptor buffers.
 *
 * @param Config Audio and score-following configuration to apply.
 *
 * @note Updates the forward model and MIR extractor, clears the input window, and resets the hop counter.
 * @warning Reallocates processing buffers; do not call concurrently with ProcessBlock().
 */
void OpenScofo::UpdateConfiguration(Configuration &Config) {
    m_Config = Config;

    size_t NHalf = static_cast<size_t>(m_Config.FFTSize / 2 + 1);
    m_Forward.UpdateConfiguration(m_Config);
    m_MIR.UpdateConfiguration(m_Config);

    m_InBuffer.resize(static_cast<size_t>(m_Config.FFTSize));
    std::fill(m_InBuffer.begin(), m_InBuffer.end(), 0.0);
    m_BlockIndex = 0;

    spdlog::debug("Allocated Description size for Window Size {}, NHalf {}", Config.FFTSize, NHalf);

    m_Desc.Magnitude.resize(NHalf);
    m_Desc.Power.resize(NHalf);
    m_Desc.SpectralMagnitudeNorm.resize(NHalf);
    m_Desc.SpectralMagnitudeFrameNorm.resize(NHalf);
    m_Desc.SpectralPowerFrameNorm.resize(NHalf);
    m_Desc.ReverbSpectralPower.resize(NHalf);
    m_Desc.LogMelSpectrum.resize(m_Config.MFCCMels);
    m_Desc.MFCC.resize(m_Config.MFCCCount);
    m_Desc.Chroma.resize(m_Config.ChromaSize);
}

// ╭─────────────────────────────────────╮
// │               Errors                │
// ╰─────────────────────────────────────╯
/**
 * @brief Register the callback used to deliver log and error messages.
 *
 * @param cb Callback invoked by the logging sink for each delivered log message.
 * @param data Caller-owned context passed to the callback.
 *
 * @note The logging sink also updates the error flag; release builds use info level and debug builds use debug
 * level.
 * @warning Replaces the previous callback. Keep its context valid while the callback is registered.
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
 * @brief Set the minimum logging level.
 *
 * @param level Minimum spdlog message level to emit.
 *
 * @note Changes the global spdlog level, affecting other users of its default logger.
 */
void OpenScofo::SetLogLevel(spdlog::level::level_enum level) {
    auto logger = spdlog::default_logger();
    spdlog::set_level(level);
}

// ─────────────────────────────────────
/**
 * @brief Apply a new processing configuration.
 *
 * @param Config Audio and score-following configuration to apply.
 *
 * @note Delegates to UpdateConfiguration(), resetting the input window and rebuilding MIR resources.
 * @warning Do not change configuration concurrently with audio processing.
 */
void OpenScofo::SetConfiguration(Configuration &Config) {
    UpdateConfiguration(Config);
}

// ─────────────────────────────────────
/**
 * @brief Replace the set of enabled audio descriptors.
 *
 * @param Descriptors Descriptor set to enable.
 *
 * @note Removes INVALID entries and duplicates, sorts the set, and rebuilds MIR resources only when the set
 * changes.
 * @warning Changing the set can allocate memory and reset analysis history; serialize with audio processing.
 */
void OpenScofo::SetRequestedDescriptors(std::vector<Descriptors> Descriptors) {
    Descriptors.erase(std::remove(Descriptors.begin(), Descriptors.end(), INVALID), Descriptors.end());
    std::sort(Descriptors.begin(), Descriptors.end());
    Descriptors.erase(std::unique(Descriptors.begin(), Descriptors.end()), Descriptors.end());

    if (m_Config.RequestedDescriptors == Descriptors) {
        return;
    }

    m_Config.RequestedDescriptors = std::move(Descriptors);
    m_MIR.UpdateConfiguration(m_Config);
}

// ─────────────────────────────────────
/**
 * @brief Enable every supported audio descriptor.
 *
 * @note Calls SetRequestedDescriptors() with the complete descriptor set, including ONNX.
 * @warning ONNX inference still requires a loaded model. Enabling descriptors can rebuild analysis resources.
 */
void OpenScofo::ActivateAllDescriptors() {
    SetRequestedDescriptors({ODSONSET,
                             LOUDNESS,
                             DB,
                             MAXAMP,
                             RMS,
                             STDDEV,
                             MAGNITUDE,
                             POWERARRAY,
                             SILENCEPROB,
                             MFCC,
                             CHROMA,
                             LOGMEL,
                             ZCR,
                             HFR,
                             CENTROID,
                             SPREADHZ,
                             SPREADVARIANCE,
                             CREST,
                             FLATNESS,
                             ENTROPY,
                             ROLLOFF,
                             CENTROIDVEL,
                             FLUX,
                             SKEWNESS,
                             SLOPE,
                             KURTOSIS,
                             IRREGULARITY,
                             HARMONICITY,
                             YIN,
                             YINCONFIDENCE,
                             EXTENDEDTECHNIQUE,
                             ONNX});
}

// ─────────────────────────────────────
/**
 * @brief Enable a descriptor and refresh the current audio description.
 *
 * @param Descriptor Descriptor to enable.
 *
 * @note Ignores INVALID, avoids duplicate requests, and analyzes the current input window when it is available.
 * @warning Can allocate memory, reset MIR history, and run analysis; serialize with audio processing.
 */
void OpenScofo::RequestDescriptor(Descriptors Descriptor) {
    if (Descriptor == INVALID) {
        return;
    }

    if (std::find(m_Config.RequestedDescriptors.begin(), m_Config.RequestedDescriptors.end(), Descriptor) ==
        m_Config.RequestedDescriptors.end()) {
        m_Config.RequestedDescriptors.push_back(Descriptor);
        m_MIR.UpdateConfiguration(m_Config);
    }

    if (!m_InBuffer.empty()) {
        m_MIR.GetDescription(m_InBuffer, m_Desc);
    }
}

// ─────────────────────────────────────
/**
 * @brief Clear a recoverable error status.
 *
 * @note Sets the stored status to info unless a critical error has already been recorded.
 * @warning A critical status remains set; this method does not repair failed processing resources.
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
 * @brief Load an ONNX classifier for audio descriptor inference.
 *
 * @param Model Path to an ONNX model file.
 * @param Descriptors Descriptors to compute, or model inputs when loading an ONNX model.
 *
 * @note Requires the .onnx extension and delegates model initialization to the MIR extractor.
 * @warning Invalid paths or incompatible models are reported through logging; serialize loading with processing.
 */
void OpenScofo::LoadONNXModel(fs::path Model, std::vector<Descriptors> Descriptors) {
    if (Model.extension() != ".onnx") {
        spdlog::error("OpenScofo just work with onnx models. Model {} is not valid", Model.string());
        return;
    }

    m_MIR.ONNXInit(Model, Descriptors);
}

// ╭─────────────────────────────────────╮
// │                 Lua                 │
// ╰─────────────────────────────────────╯
#if defined(OPENSCOFO_LUA)
/**
 * @brief Create the embedded Lua runtime and register OpenScofo bindings.
 *
 * @note Opens standard libraries and stores this instance as lightuserdata in the global _OpenScofo table.
 * @warning Closes any existing Lua runtime and releases its pending timers before creating a new one.
 */
void OpenScofo::InitLuaModule() {
    CloseLuaModule();
    m_LuaCurrentSample = 0;
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
 * @brief Register a module in the embedded Lua runtime.
 *
 * @param name Name exposed in the Lua runtime.
 * @param func Lua C function that opens the module.
 *
 * @return True if a non-nil module was registered; false if the runtime is absent or the result is nil.
 *
 * @note Uses luaL_requiref(), which caches loaded modules, and leaves the module value on the Lua stack.
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
 * @brief Compile and execute Lua source in the embedded runtime.
 *
 * @param code Lua source code to execute.
 *
 * @return True on successful execution; false if no runtime exists or compilation or execution fails.
 *
 * @note Retains Lua results or the error object on the stack; LuaGetError() retrieves and pops an error.
 * @warning Runs synchronously and may execute arbitrary registered callbacks; serialize access to the Lua state.
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
 * @brief Expose a caller-owned pointer as a Lua global.
 *
 * @param pointer Caller-owned pointer exposed as lightuserdata.
 * @param name Name exposed in the Lua runtime.
 *
 * @return True if the Lua runtime exists and the global was assigned; false otherwise.
 *
 * @note Stores lightuserdata without taking ownership or managing the pointed-to object.
 * @warning Keep the pointed-to object alive for as long as Lua code can access it.
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
 * @brief Append a directory to the Lua module search path.
 *
 * @param path Directory to append to the Lua module search path.
 *
 * @note Adds a directory/?.lua pattern to package.path when the Lua runtime exists.
 * @warning The directory path must not be empty.
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
 * @brief Retrieve and remove the top Lua stack value as an error message.
 *
 * @return The error text or a fallback diagnostic.
 *
 * @note Returns a fallback message if no string value is available, or if the Lua runtime does not exist.
 * @warning Pops the top stack value only if it can be converted to a string; other values remain on the stack.
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
 * @brief Select an internal score state and reset input analysis buffers.
 *
 * @param Event Zero-based internal score state index; zero selects the initial state.
 *
 * @note Resets forward decoding for a valid index and maps that state to the public score position.
 * @warning Invalid indices are logged by the forward model, but input buffers and the public position still reset.
 */
void OpenScofo::SetCurrentEvent(int Event) {
    m_CurrentScorePosition = 0;
    m_BlockIndex = 0;
    m_Forward.SetCurrentEvent(Event);
    std::fill(m_InBuffer.begin(), m_InBuffer.end(), 0.0);

    std::fill(m_Desc.Magnitude.begin(), m_Desc.Magnitude.end(), 0.0);
    std::fill(m_Desc.Power.begin(), m_Desc.Power.end(), 0.0);
    std::fill(m_Desc.SpectralMagnitudeNorm.begin(), m_Desc.SpectralMagnitudeNorm.end(), 0.0);
    std::fill(m_Desc.SpectralMagnitudeFrameNorm.begin(), m_Desc.SpectralMagnitudeFrameNorm.end(), 0.0);
    std::fill(m_Desc.SpectralPowerFrameNorm.begin(), m_Desc.SpectralPowerFrameNorm.end(), 0.0);
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

// ─────────────────────────────────────
/**
 * @brief Restart decoding at the first state of a named section.
 *
 * @param Section Section name without surrounding quotes.
 *
 * @return True if the section was found and selected; false otherwise.
 *
 * @note A successful selection clears input and spectral buffers and refreshes the public score position.
 * @warning Selecting a section resets decoding history; serialize with audio processing.
 */
bool OpenScofo::SetCurrentSection(const std::string &Section) {
    if (!m_Forward.SetCurrentSection(Section)) {
        return false;
    }

    m_BlockIndex = 0;
    std::fill(m_InBuffer.begin(), m_InBuffer.end(), 0.0);
    std::fill(m_Desc.Magnitude.begin(), m_Desc.Magnitude.end(), 0.0);
    std::fill(m_Desc.Power.begin(), m_Desc.Power.end(), 0.0);
    std::fill(m_Desc.SpectralMagnitudeNorm.begin(), m_Desc.SpectralMagnitudeNorm.end(), 0.0);
    std::fill(m_Desc.SpectralMagnitudeFrameNorm.begin(), m_Desc.SpectralMagnitudeFrameNorm.end(), 0.0);
    std::fill(m_Desc.SpectralPowerFrameNorm.begin(), m_Desc.SpectralPowerFrameNorm.end(), 0.0);
    std::fill(m_Desc.ReverbSpectralPower.begin(), m_Desc.ReverbSpectralPower.end(), 0.0);

    const int StateIndex = m_Forward.GetCurrentStateIndex();
    m_CurrentScorePosition = m_States[static_cast<size_t>(StateIndex)].ScorePos;
    return true;
}

// ╭─────────────────────────────────────╮
// │            Get Functions            │
// ╰─────────────────────────────────────╯
/**
 * @brief Read the current public score position.
 *
 * @return Current public score position.
 *
 * @note Returns the position last recorded by this API; rests do not advance the public event numbering.
 */
int OpenScofo::GetCurrentScorePosition() {
    return m_CurrentScorePosition;
}

// ─────────────────────────────────────
/**
 * @brief Read the active internal forward-model state index.
 *
 * @return Zero-based index of the active internal state.
 *
 * @note Internal indices include boundary and silence states and can differ from public score positions.
 */
int OpenScofo::GetCurrentStateIndex() {
    return m_Forward.GetCurrentStateIndex();
}

// ─────────────────────────────────────
/**
 * @brief Read the current tempo estimate.
 *
 * @return Estimated tempo in beats per minute.
 *
 * @note Delegates to the forward model; the initial estimate comes from the selected score state.
 */
double OpenScofo::GetCurrentBPM() {
    return m_Forward.GetCurrentBPM();
}

// ─────────────────────────────────────
/**
 * @brief Copy the actions associated with the current score event.
 *
 * @return Actions for the active event, or an empty list.
 *
 * @note Returns an empty list when no score has been loaded; reading does not consume event actions.
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
 * @brief Retrieve pending audio-state-change actions.
 *
 * @return Pending audio-state-change actions.
 *
 * @note Returns an empty list when no score is loaded; otherwise drains the forward-model action queue.
 * @warning Each queued action is returned only once.
 */
EventActions OpenScofo::GetAudioStateChangeActions() {
    if (ScoreIsLoaded()) {
        return m_Forward.GetAudioStateChangeActions();
    }
    return {};
}

// ─────────────────────────────────────
/**
 * @brief Evaluate pitch evidence against the current audio frame.
 *
 * @param Freq Target pitch frequency in Hz.
 *
 * @return Pitch evidence score from the forward model.
 *
 * @note Copies the current description to the forward model before evaluating its spectral pitch template.
 */
double OpenScofo::GetPitchProb(double Freq) {
    m_Forward.SetDescription(m_Desc);
    return m_Forward.GetPitchProbability(Freq);
}

// ─────────────────────────────────────
/**
 * @brief Copy the global Lua source collected from the score.
 *
 * @return Collected Lua source code.
 *
 * @note Returns source from score-level LUA blocks without executing it.
 */
std::string OpenScofo::GetLuaCode() {
    return m_Score.GetLuaCode();
}

// ╭─────────────────────────────────────╮
// │          Helpers Functions          │
// ╰─────────────────────────────────────╯
/**
 * @brief Read the score parser loaded flag.
 *
 * @return True when the parser reports a loaded score; false otherwise.
 *
 * @note Reflects the parser status rather than independently validating decoder or ONNX readiness.
 */
bool OpenScofo::ScoreIsLoaded() {
    return m_Score.ScoreIsLoaded();
}

// ─────────────────────────────────────
/**
 * @brief Convert a descriptor enum to its public identifier.
 *
 * @param d Descriptor enum value.
 *
 * @return Null-terminated descriptor identifier.
 *
 * @note Returns a string literal; unknown enum values map to "unknown".
 */
const char *OpenScofo::GetDescriptionId(Descriptors d) {
    switch (d) {
    case Descriptors::ODSONSET:
        return "onset";
    case Descriptors::LOUDNESS:
        return "loudness";
    case Descriptors::DB:
        return "db";
    case Descriptors::MAXAMP:
        return "maxamp";
    case Descriptors::RMS:
        return "rms";
    case Descriptors::STDDEV:
        return "stddev";
    case Descriptors::MAGNITUDE:
        return "magnitude";
    case Descriptors::POWERARRAY:
        return "power";
    case Descriptors::SILENCEPROB:
        return "silence";
    case Descriptors::MFCC:
        return "mfcc";
    case Descriptors::CHROMA:
        return "chroma";
    case Descriptors::LOGMEL:
        return "logmel";
    case Descriptors::ZCR:
        return "zcr";
    case Descriptors::HFR:
        return "hfr";
    case Descriptors::CENTROID:
        return "centroid";
    case Descriptors::SPREADHZ:
        return "spread";
    case Descriptors::SPREADVARIANCE:
        return "spread_variance";
    case Descriptors::CREST:
        return "crest";
    case Descriptors::FLATNESS:
        return "flatness";
    case Descriptors::ENTROPY:
        return "entropy";
    case Descriptors::ROLLOFF:
        return "rolloff";
    case Descriptors::CENTROIDVEL:
        return "centroid_velocity";
    case Descriptors::FLUX:
        return "flux";
    case Descriptors::SKEWNESS:
        return "skewness";
    case Descriptors::SLOPE:
        return "slope";
    case Descriptors::KURTOSIS:
        return "kurtosis";
    case Descriptors::IRREGULARITY:
        return "irregularity";
    case Descriptors::HARMONICITY:
        return "harmonicity";
    case Descriptors::YIN:
        return "yin";
    case Descriptors::YINCONFIDENCE:
        return "yin_confidence";
    case Descriptors::EXTENDEDTECHNIQUE:
        return "ext";
    case Descriptors::ONNX:
        return "onnx";
    default:
        return "unknown";
    }
}

// ─────────────────────────────────────
/**
 * @brief Convert a public descriptor identifier to its enum.
 *
 * @param s Null-terminated descriptor identifier.
 *
 * @return Matching enum value, or INVALID for an unknown name.
 *
 * @note Accepts supported aliases and logs an error for unknown identifiers.
 * @warning The identifier pointer must not be null.
 */
Descriptors OpenScofo::GetDescriptorsEnum(const char *s) {
    if (strcmp(s, "mfcc") == 0) {
        return Descriptors::MFCC;
    } else if (strcmp(s, "logmel") == 0) {
        return Descriptors::LOGMEL;
    } else if (strcmp(s, "rms") == 0) {
        return Descriptors::RMS;
    } else if (strcmp(s, "loudness") == 0) {
        return Descriptors::LOUDNESS;
    } else if (strcmp(s, "db") == 0) {
        return Descriptors::DB;
    } else if (strcmp(s, "maxamp") == 0 || strcmp(s, "max_amp") == 0) {
        return Descriptors::MAXAMP;
    } else if (strcmp(s, "magnitude") == 0) {
        return Descriptors::MAGNITUDE;
    } else if (strcmp(s, "power") == 0 || strcmp(s, "powerarray") == 0) {
        return Descriptors::POWERARRAY;
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
    } else if (strcmp(s, "spread_variance") == 0) {
        return Descriptors::SPREADVARIANCE;
    } else if (strcmp(s, "crest") == 0) {
        return Descriptors::CREST;
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
    } else if (strcmp(s, "yin") == 0) {
        return Descriptors::YIN;
    } else if (strcmp(s, "yin_confidence") == 0 || strcmp(s, "pitch_confidence") == 0) {
        return Descriptors::YINCONFIDENCE;
    } else if (strcmp(s, "onnx") == 0) {
        return Descriptors::ONNX;
    } else {
        spdlog::error("Invalid descriptors argument: {}", s);
        return Descriptors::INVALID;
    }
}

// ─────────────────────────────────────
/**
 * @brief Read a scalar value from an audio description.
 *
 * @param Desc Description containing the requested descriptor.
 * @param d Descriptor enum value.
 *
 * @return Requested scalar value, or -1.0 for an unsupported descriptor.
 *
 * @note Vector-valued and unsupported descriptors produce an error log and return -1.0.
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
    case Descriptors::LOGMEL:
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
 * @brief Access an array-valued descriptor in an audio description.
 *
 * @param Desc Description containing the requested descriptor.
 * @param d Descriptor enum value.
 *
 * @return Mutable reference to the requested array, or the magnitude array as a fallback.
 *
 * @note Supports MFCC, chroma, log-mel, power, and magnitude arrays.
 * @warning Unsupported descriptors log a critical message and return Magnitude; no exception is thrown.
 */
std::vector<double> &OpenScofo::GetDescriptionArray(Description &Desc, Descriptors d) {
    switch (d) {
    case Descriptors::MFCC:
        return Desc.MFCC;
    case Descriptors::CHROMA:
        return Desc.Chroma;
    case Descriptors::LOGMEL:
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
 * @brief Access the mutable forward-model score states.
 *
 * @return Reference to the forward-model states.
 *
 * @note Returns the actual state container rather than a copy.
 * @warning Mutations can invalidate decoder assumptions or references; do not modify states during processing.
 */
States &OpenScofo::GetStates() {
    return m_Forward.GetStates();
}

// ─────────────────────────────────────
/**
 * @brief Build or retrieve a spectral pitch template.
 *
 * @param Freq Target pitch frequency in Hz.
 *
 * @return Pitch template bins; an unsupported frequency can produce an empty template.
 *
 * @note Delegates to the forward model and returns a copy of its cached template.
 */
std::vector<double> OpenScofo::GetPitchTemplate(double Freq) {
    return m_Forward.GetPitchTemplate(Freq);
}

// ─────────────────────────────────────
/**
 * @brief Read the configured sampling rate.
 *
 * @return Sampling rate in Hz.
 *
 * @note Reads the value stored in the current configuration.
 */
double OpenScofo::GetSr() {
    return m_Config.SR;
}

// ─────────────────────────────────────
/**
 * @brief Read the configured analysis window size.
 *
 * @return FFT window size in samples.
 *
 * @note Reads the value stored in the current configuration.
 */
double OpenScofo::GetFFTSize() {
    return m_Config.FFTSize;
}

// ─────────────────────────────────────
/**
 * @brief Read the configured analysis hop size.
 *
 * @return Hop size in samples.
 *
 * @note Reads the value stored in the current configuration.
 */
double OpenScofo::GetHopSize() {
    return m_Config.HOPSize;
}

// ─────────────────────────────────────
/**
 * @brief Read the duration of one analysis hop.
 *
 * @return Hop duration in seconds.
 *
 * @note Delegates to the forward model, which computes hop size divided by sampling rate.
 */
double OpenScofo::GetBlockDuration() {
    return m_Forward.GetBlockDuration();
}

// ╭─────────────────────────────────────╮
// │           Main Functions            │
// ╰─────────────────────────────────────╯
/**
 * @brief Copy the most recent audio description.
 *
 * @return Copy of the current description.
 *
 * @note The returned value is a snapshot and does not track subsequent frames.
 */
Description OpenScofo::GetDescription() {
    return m_Desc;
}

// ─────────────────────────────────────
/**
 * @brief Copy the current processing configuration.
 *
 * @return Copy of the current configuration.
 *
 * @note Changing the returned copy does not apply it; use SetConfiguration() to apply changes.
 */
Configuration OpenScofo::GetConfiguration() {
    return m_Config;
}

// ─────────────────────────────────────
/**
 * @brief Read the forward-model circular history index.
 *
 * @return Current slot in the forward-model history buffer.
 *
 * @note This index tracks inference frames, rather than the number of input samples buffered.
 */
int OpenScofo::GetCurrentBufferIndex() {
    return m_Forward.GetCurrentBufferIndex();
}

// ─────────────────────────────────────
/**
 * @brief Parse a score and initialize its analysis and decoding configuration.
 *
 * @param ScorePath Path to the score file.
 *
 * @return True if initialization completes without an error or critical log status; false otherwise.
 *
 * @note Preserves requested descriptors, enables score-required features, and validates technique labels against
 * ONNX.
 * @warning Replaces score and processing state even when loading fails; serialize with audio processing.
 */
bool OpenScofo::LoadScore(fs::path ScorePath) {
    ClearErrors();

    m_CurrentScorePosition = 0;
    const std::vector<Descriptors> requestedDescriptors = m_Config.RequestedDescriptors;
    int SR = m_Config.SR;
    auto [newConfig, newStates] = m_Score.Parse(ScorePath);
    newConfig.RequestedDescriptors = requestedDescriptors;
    m_Config = newConfig;
    m_States = newStates;

    if (m_Config.SR != SR) {
        spdlog::error("Sample rate mismatch: OpenScofo is running at {} Hz, but the score file requires {} Hz.",
                      m_Config.SR, SR);
    }

    auto requestScoreDescriptor = [&](Descriptors Descriptor) {
        if (std::find(newConfig.RequestedDescriptors.begin(), newConfig.RequestedDescriptors.end(), Descriptor) ==
            newConfig.RequestedDescriptors.end()) {
            newConfig.RequestedDescriptors.push_back(Descriptor);
        }
    };

    for (const ScoreState &state : m_States) {
        for (const Observation &audioState : state.Observations) {
            if (audioState.Type == LABEL) {
                requestScoreDescriptor(ONNX);
                requestScoreDescriptor(EXTENDEDTECHNIQUE);
            } else if (audioState.Type == ONSET) {
                requestScoreDescriptor(ODSONSET);
            }
        }
        for (const MarkovMicroState &microState : state.MicroStates) {
            for (const Observation &audioState : microState.Observations) {
                if (audioState.Type == LABEL) {
                    requestScoreDescriptor(ONNX);
                    requestScoreDescriptor(EXTENDEDTECHNIQUE);
                } else if (audioState.Type == ONSET) {
                    requestScoreDescriptor(ODSONSET);
                }
            }
        }
    }

    UpdateConfiguration(newConfig);

    // Timbre/Extended Tech detection
    if (fs::exists(newConfig.TimbreONNXModel)) {
        std::vector<std::string> descriptors = newConfig.ONNXDescriptors;
        spdlog::info("Loading ONNX model, wait...");
        std::vector<Descriptors> DescEnum;
        for (auto d : descriptors) {
            DescEnum.push_back(GetDescriptorsEnum(d.c_str()));
        }
        m_MIR.ONNXInit(newConfig.TimbreONNXModel, DescEnum);
        spdlog::info("ONNX Model ready");
    }

    // Add States
    m_Forward.SetScoreStates(m_States);
    if (m_Config.SectionRestrict && !m_States.empty()) {
        const int StateIndex = m_Forward.GetCurrentStateIndex();
        m_CurrentScorePosition = m_States[static_cast<size_t>(StateIndex)].ScorePos;
    }
    m_Mode = SCOREFOLLOWER;

    // verify the states.
    const std::vector<std::string> &ONNXLabels = m_MIR.GetONNXLabels();
    for (auto &state : m_States) {
        for (const Observation &audioState : state.Observations) {
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
                    spdlog::error("Extended Technique Label '{}' is not valid on line {}. Valid labels: {}", label,
                                  state.Line, validLabels);
                    return false;
                }
            }
        }
        for (const MarkovMicroState &microState : state.MicroStates) {
            for (const Observation &audioState : microState.Observations) {
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
                        spdlog::error("Extended Technique Label '{}' is not valid on line {}. Valid labels: {}", label,
                                      state.Line, validLabels);
                        return false;
                    }
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
 * @brief Append input samples and process an audio analysis frame when a hop is due.
 *
 * @tparam T Audio sample type, constrained to float or double.
 * @param AudioBuffer Readable input samples in float or double precision.
 * @param n Number of input samples; must not exceed the configured FFT window size.
 *
 * @return True after buffering or processing the block.
 *
 * @note Shifts the input window and performs at most one analysis per call. Lua timers run on the calling thread.
 * @warning AudioBuffer must cover n samples and n must fit the FFT window. Calls exceeding one hop do not process
 * every hop.
 */
template <OpenScofoPrecision T> bool OpenScofo::ProcessBlock(const T *AudioBuffer, size_t n) {
#if defined(OPENSCOFO_LUA)
    m_LuaCurrentSample += n;
    ProcessLuaTimers();
#endif
    m_BlockIndex += n;

    std::copy(m_InBuffer.begin() + n, m_InBuffer.end(), m_InBuffer.begin());
    std::transform(AudioBuffer, AudioBuffer + n, m_InBuffer.end() - n, [](T x) { return static_cast<double>(x); });

    if (m_BlockIndex < m_Config.HOPSize) {
        return true;
    }

    m_BlockIndex -= m_Config.HOPSize;

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
