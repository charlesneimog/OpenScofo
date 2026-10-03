/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file lua.cpp
 * @brief Embedded Lua bindings, table serialization, and sample-clock timers.
 *
 * @note Compiled only when OPENSCOFO_LUA is enabled; timer callbacks run during audio block processing.
 * @warning Lua state and timer operations must be serialized with processing and runtime destruction.
 */

#include <OpenScofo.hpp>

#if defined(OPENSCOFO_LUA)

#include <limits>

namespace OpenScofo {

// ─────────────────────────────────────
/**
 * @brief Release the embedded Lua runtime when the instance is destroyed.
 *
 * @note Calls CloseLuaModule() to release pending callback references before closing Lua.
 * @warning Destroy the instance only after audio processing and Lua access have stopped.
 */
OpenScofo::~OpenScofo() {
    CloseLuaModule();
}

// ─────────────────────────────────────
/**
 * @brief Release pending timers and close the embedded Lua state.
 *
 * @note Unreferences callback functions and data; repeated calls are harmless after the state is closed.
 * @warning Invalidates all external references to this Lua state and cancels its pending callbacks.
 */
void OpenScofo::CloseLuaModule() {
    if (m_LuaState == nullptr) {
        return;
    }
    for (const auto &[Key, Timer] : m_LuaTimers) {
        luaL_unref(m_LuaState, LUA_REGISTRYINDEX, Timer.CallbackRef);
        luaL_unref(m_LuaState, LUA_REGISTRYINDEX, Timer.DataRef);
    }
    m_LuaTimers.clear();
    lua_close(m_LuaState);
    m_LuaState = nullptr;
}

// ─────────────────────────────────────
/**
 * @brief Schedule a callback against the processed-audio sample clock.
 *
 * @param DelayMs Finite nonnegative delay in milliseconds, rounded up to an audio sample.
 * @param CallbackRef Lua registry reference to the callback function.
 * @param DataRef Lua registry reference to callback data, or LUA_REFNIL.
 *
 * @return Nonzero timer identifier, or zero when the delay, sample deadline, or identifier is out of range.
 *
 * @note A successful insertion takes ownership of the registry references until execution, cancellation, or runtime
 * closure.
 * @warning References must belong to this instance Lua state. On failure, the caller remains responsible for
 * releasing them.
 */
uint64_t OpenScofo::ScheduleLuaCallback(double DelayMs, int CallbackRef, int DataRef) {
    const double DelaySamples = std::ceil(DelayMs * 0.001 * m_Config.SR);
    if (!std::isfinite(DelaySamples) || DelaySamples < 0 ||
        DelaySamples >= static_cast<double>(std::numeric_limits<uint64_t>::max() - m_LuaCurrentSample) ||
        m_LuaNextTimerId > static_cast<uint64_t>(LUA_MAXINTEGER)) {
        return 0;
    }
    const uint64_t Id = m_LuaNextTimerId++;
    m_LuaTimers.emplace(LuaTimerKey{m_LuaCurrentSample + static_cast<uint64_t>(DelaySamples), Id},
                        LuaTimer{CallbackRef, DataRef});
    return Id;
}

// ─────────────────────────────────────
/**
 * @brief Cancel a pending callback by timer identifier.
 *
 * @param Id Identifier returned by ScheduleLuaCallback().
 *
 * @return True if a pending timer was removed; false if the identifier was not found.
 *
 * @note Releases the callback and data registry references when a matching timer is removed.
 */
bool OpenScofo::CancelLuaCallback(uint64_t Id) {
    for (auto It = m_LuaTimers.begin(); It != m_LuaTimers.end(); ++It) {
        if (It->first.second == Id) {
            luaL_unref(m_LuaState, LUA_REGISTRYINDEX, It->second.CallbackRef);
            luaL_unref(m_LuaState, LUA_REGISTRYINDEX, It->second.DataRef);
            m_LuaTimers.erase(It);
            return true;
        }
    }
    return false;
}

// ─────────────────────────────────────
/**
 * @brief Execute callbacks whose sample-clock deadlines have passed.
 *
 * @note Removes each timer before invocation so callbacks can alter other timers; logs Lua errors and releases
 * references.
 * @warning Runs callbacks synchronously on the ProcessBlock() caller thread; deadlines are checked at block
 * boundaries.
 */
void OpenScofo::ProcessLuaTimers() {
    while (!m_LuaTimers.empty()) {
        auto It = m_LuaTimers.begin();
        if (It->first.first > m_LuaCurrentSample) {
            break;
        }
        const LuaTimer Timer = It->second;
        m_LuaTimers.erase(It); // Callbacks may schedule or cancel other timers.
        lua_rawgeti(m_LuaState, LUA_REGISTRYINDEX, Timer.CallbackRef);
        if (Timer.DataRef == LUA_REFNIL) {
            lua_pushnil(m_LuaState);
        } else {
            lua_rawgeti(m_LuaState, LUA_REGISTRYINDEX, Timer.DataRef);
        }
        if (lua_pcall(m_LuaState, 1, 0, 0) != LUA_OK) {
            const char *Error = lua_tostring(m_LuaState, -1);
            spdlog::error("Lua timer callback: {}", Error ? Error : "Unknown error");
            lua_pop(m_LuaState, 1);
        }
        luaL_unref(m_LuaState, LUA_REGISTRYINDEX, Timer.CallbackRef);
        luaL_unref(m_LuaState, LUA_REGISTRYINDEX, Timer.DataRef);
    }
}

// ─────────────────────────────────────
/**
 * @brief Resolve the instance pointer stored in the global Lua binding table.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return Instance pointer, or null when the binding table or pointer is unavailable.
 *
 * @note Restores the stack after inspecting _OpenScofo.pointer.
 * @warning The table stores non-owning lightuserdata; the pointed-to instance must remain alive.
 */
static OpenScofo *GetCurrentOpenScofo(lua_State *L) {
    lua_getglobal(L, "_OpenScofo");
    if (!lua_istable(L, -1)) {
        lua_pop(L, 1);
        return nullptr;
    }

    lua_getfield(L, -1, "pointer");
    void *pointer = lua_touserdata(L, -1);
    lua_pop(L, 2);
    return static_cast<OpenScofo *>(pointer);
}

// ─────────────────────────────────────
/**
 * @brief Push numeric values as a Lua sequence table.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 * @param values Numeric values to copy into a Lua table.
 *
 * @note Adds one table to the stack and uses one-based Lua array indices.
 */
static void PushNumberVector(lua_State *L, const std::vector<double> &values) {
    lua_createtable(L, static_cast<int>(values.size()), 0);
    for (size_t i = 0; i < values.size(); ++i) {
        lua_pushnumber(L, values[i]);
        lua_rawseti(L, -2, static_cast<int>(i + 1));
    }
}

// ─────────────────────────────────────
/**
 * @brief Push an audio observation as a Lua table.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 * @param state Score or audio state to serialize into a Lua table.
 *
 * @note Adds one table containing type, frequency, MIDI pitch, observation index, and label.
 */
static void PushAudioState(lua_State *L, const Observation &state) {
    lua_createtable(L, 0, 5);
    lua_pushinteger(L, state.Type);
    lua_setfield(L, -2, "type");
    lua_pushnumber(L, state.Freq);
    lua_setfield(L, -2, "freq");
    lua_pushnumber(L, state.Midi);
    lua_setfield(L, -2, "midi");
    lua_pushinteger(L, static_cast<lua_Integer>(state.Index));
    lua_setfield(L, -2, "index");
    lua_pushlstring(L, state.Label.data(), state.Label.size());
    lua_setfield(L, -2, "label");
}

/**
 * @brief Push audio observations as a Lua sequence table.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 * @param observations Audio observations to serialize in their original order.
 *
 * @note Adds one table of serialized observations using one-based Lua indices.
 */
static void PushObservations(lua_State *L, const std::vector<Observation> &observations) {
    lua_createtable(L, static_cast<int>(observations.size()), 0);
    for (size_t i = 0; i < observations.size(); ++i) {
        PushAudioState(L, observations[i]);
        lua_rawseti(L, -2, static_cast<int>(i + 1));
    }
}

// // ─────────────────────────────────────
/**
 * @brief Push an audio descriptor snapshot as a Lua table.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 * @param desc Audio description snapshot to serialize.
 *
 * @note Copies selected scalar descriptors and magnitude, MFCC, and chroma arrays; adds one table to the stack.
 */
static void PushDescription(lua_State *L, const Description &desc) {
    lua_createtable(L, 0, 20);

    lua_pushboolean(L, desc.Onset);
    lua_setfield(L, -2, "onset");
    lua_pushnumber(L, desc.SilenceProb);
    lua_setfield(L, -2, "silence");

    lua_pushnumber(L, desc.dB);
    lua_setfield(L, -2, "db");
    lua_pushnumber(L, desc.RMS);
    lua_setfield(L, -2, "rms");
    lua_pushnumber(L, desc.MaxAmp);
    lua_setfield(L, -2, "max_amp");
    lua_pushnumber(L, desc.Loudness);
    lua_setfield(L, -2, "loudness");

    lua_pushnumber(L, desc.Harmonicity);
    lua_setfield(L, -2, "harmonicity");
    lua_pushnumber(L, desc.SpectralFlatness);
    lua_setfield(L, -2, "spectral_flatness");
    lua_pushnumber(L, desc.SpectralEntropy);
    lua_setfield(L, -2, "spectral_entropy");
    lua_pushnumber(L, desc.SpectralRolloff);
    lua_setfield(L, -2, "spectral_rolloff");
    lua_pushnumber(L, desc.SpectralFlux);
    lua_setfield(L, -2, "spectral_flux");
    lua_pushnumber(L, desc.SpectralIrregularity);
    lua_setfield(L, -2, "spectral_irregularity");
    lua_pushnumber(L, desc.SpectralIrregularityJensen);
    lua_setfield(L, -2, "spectral_irregularity_jensen");
    lua_pushnumber(L, desc.SpectralIrregularityKrimphoff);
    lua_setfield(L, -2, "spectral_irregularity_krimphoff");
    lua_pushnumber(L, desc.SpectralCrest);
    lua_setfield(L, -2, "spectral_crest");
    lua_pushnumber(L, desc.SpectralCentroid);
    lua_setfield(L, -2, "spectral_centroid");
    lua_pushnumber(L, desc.CentroidVelocity);
    lua_setfield(L, -2, "centroid_velocity");
    lua_pushnumber(L, desc.SpectralSpreadHz);
    lua_setfield(L, -2, "spectral_spread_hz");
    lua_pushnumber(L, desc.HighFreqRatio);
    lua_setfield(L, -2, "high_freq_ratio");

    PushNumberVector(L, desc.SpectralMagnitudeNorm);
    lua_setfield(L, -2, "spectral_magnitude");
    PushNumberVector(L, desc.MFCC);
    lua_setfield(L, -2, "mfcc");
    PushNumberVector(L, desc.Chroma);
    lua_setfield(L, -2, "chroma");
}

// ─────────────────────────────────────
/**
 * @brief Push a score state and its runtime data as a Lua table.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 * @param state Score or audio state to serialize into a Lua table.
 *
 * @note Includes timing, forward probabilities, audio observations, topology, and microstate data.
 */
static void PushMarkovState(lua_State *L, const ScoreState &state) {
    lua_createtable(L, 0, 15);
    lua_pushinteger(L, state.ScorePos);
    lua_setfield(L, -2, "position");
    lua_pushlstring(L, state.Section.data(), state.Section.size());
    lua_setfield(L, -2, "section");
    lua_pushinteger(L, state.Type);
    lua_setfield(L, -2, "type");
    lua_pushinteger(L, state.HSMMType);
    lua_setfield(L, -2, "markov");

    PushNumberVector(L, state.Forward);
    lua_setfield(L, -2, "forward");
    lua_pushnumber(L, state.BPMExpected);
    lua_setfield(L, -2, "bpm_expected");
    lua_pushnumber(L, state.BPMObserved);
    lua_setfield(L, -2, "bpm_observed");
    lua_pushnumber(L, state.OnsetExpected);
    lua_setfield(L, -2, "onset_expected");
    lua_pushnumber(L, state.OnsetObserved);
    lua_setfield(L, -2, "onset_observed");
    lua_pushnumber(L, state.PhaseExpected);
    lua_setfield(L, -2, "phase_expected");
    lua_pushnumber(L, state.PhaseObserved);
    lua_setfield(L, -2, "phase_observed");
    lua_pushnumber(L, state.IOIPhiN);
    lua_setfield(L, -2, "ioi_phi_n");
    lua_pushnumber(L, state.IOIHatPhiN);
    lua_setfield(L, -2, "ioi_hat_phi_n");
    lua_pushnumber(L, state.Duration);
    lua_setfield(L, -2, "duration");
    lua_pushinteger(L, state.Line);
    lua_setfield(L, -2, "line");

    PushObservations(L, state.Observations);
    lua_setfield(L, -2, "audiostates");

    lua_pushinteger(L, state.MicroTopologyType);
    lua_setfield(L, -2, "micro_topology");
    lua_pushboolean(L, state.IsInterEventSilence);
    lua_setfield(L, -2, "inter_event_silence");
    lua_pushinteger(L, state.BestMicroStateIndex);
    lua_setfield(L, -2, "best_microstate_index");
    lua_createtable(L, static_cast<int>(state.MicroStates.size()), 0);
    for (size_t i = 0; i < state.MicroStates.size(); ++i) {
        lua_createtable(L, 0, 2);
        PushObservations(L, state.MicroStates[i].Observations);
        lua_setfield(L, -2, "observations");
        lua_pushnumber(L, state.MicroStates[i].DurationWeight);
        lua_setfield(L, -2, "duration_weight");
        lua_rawseti(L, -2, static_cast<int>(i + 1));
    }
    lua_setfield(L, -2, "microstates");
}

// ─────────────────────────────────────
/**
 * @brief Bind internal event selection to Lua.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return Zero Lua return values.
 *
 * @note Reads an integer from argument 1 and resets the selected instance state.
 * @warning Raises a Lua error for a missing instance or invalid argument type; selection resets decoding history.
 */
static int OpenScofoSetCurrentEvent(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    self->SetCurrentEvent(static_cast<int>(luaL_checkinteger(L, 1)));
    return 0;
}

// ─────────────────────────────────────
/**
 * @brief Bind section selection to Lua.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return One Lua return value: a boolean.
 *
 * @note Reads a section name from argument 1 and pushes the selection success flag.
 * @warning Raises a Lua error for a missing instance or invalid argument type; selection resets decoding history.
 */
static int OpenScofoSetCurrentSection(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    lua_pushboolean(L, self->SetCurrentSection(luaL_checkstring(L, 1)));
    return 1;
}

// ─────────────────────────────────────
/**
 * @brief Expose the current tempo estimate to Lua.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return One Lua return value: tempo in BPM.
 *
 * @note Resolves the bound instance and pushes its BPM estimate.
 * @warning Raises a Lua error when no instance pointer is available.
 */
static int OpenScofoGetLiveBPM(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    lua_pushnumber(L, self->GetCurrentBPM());
    return 1;
}

// ─────────────────────────────────────
/**
 * @brief Expose the current public score position to Lua.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return One Lua return value: the public score position.
 *
 * @note Pushes the public score position rather than the internal decoder state index.
 * @warning Raises a Lua error when no instance pointer is available.
 */
static int OpenScofoGetEventIndex(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    lua_pushinteger(L, self->GetCurrentScorePosition());
    return 1;
}

// ─────────────────────────────────────
/**
 * @brief Expose score state snapshots as a Lua sequence.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return One Lua return value: the score state table.
 *
 * @note Copies the current states and serializes each with one-based Lua table indices.
 * @warning Raises a Lua error when no instance pointer is available.
 */
static int OpenScofoGetStates(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");

    States states = self->GetStates();
    lua_createtable(L, static_cast<int>(states.size()), 0);
    for (size_t i = 0; i < states.size(); ++i) {
        PushMarkovState(L, states[i]);
        lua_rawseti(L, -2, static_cast<int>(i + 1));
    }
    return 1;
}

// ─────────────────────────────────────
/**
 * @brief Expose the latest audio description snapshot to Lua.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return One Lua return value: the descriptor table.
 *
 * @note Copies the description and serializes its supported fields.
 * @warning Raises a Lua error when no instance pointer is available.
 */
static int OpenScofoGetCurrentDescription(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");

    Description desc = self->GetDescription();

    lua_newtable(L);
    PushDescription(L, desc);
    return 1;
}

// ─────────────────────────────────────
/**
 * @brief Bind activation of all descriptors to Lua.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return Zero Lua return values.
 *
 * @note Enables optional stages, including ONNX when a model is loaded.
 * @warning Raises a Lua error for a missing instance; activation can rebuild and allocate analysis resources.
 */
static int OpenScofoActivateAllDescriptors(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    self->ActivateAllDescriptors();
    return 0;
}

// ─────────────────────────────────────
/**
 * @brief Bind sample-clock callback scheduling to Lua.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return One Lua return value: the nonzero timer identifier.
 *
 * @note Reads delay, callback, and optional data from arguments 1 through 3 and retains them in the Lua registry.
 * @warning Invalid delays, callbacks, or exhausted timer ranges raise Lua errors; callbacks run during audio
 * processing.
 */
static int OpenScofoSchedule(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    luaL_checktype(L, 1, LUA_TNUMBER);
    const double DelayMs = lua_tonumber(L, 1);
    luaL_argcheck(L, std::isfinite(DelayMs) && DelayMs >= 0, 1, "delay must be finite and non-negative");
    luaL_checktype(L, 2, LUA_TFUNCTION);
    lua_settop(L, 3); // Missing data becomes nil.
    lua_pushvalue(L, 2);
    const int CallbackRef = luaL_ref(L, LUA_REGISTRYINDEX);
    lua_pushvalue(L, 3);
    const int DataRef = luaL_ref(L, LUA_REGISTRYINDEX);
    const uint64_t Id = self->ScheduleLuaCallback(DelayMs, CallbackRef, DataRef);
    if (Id == 0) {
        luaL_unref(L, LUA_REGISTRYINDEX, CallbackRef);
        luaL_unref(L, LUA_REGISTRYINDEX, DataRef);
        return luaL_error(L, "Lua timer delay or ID is out of range");
    }
    lua_pushinteger(L, static_cast<lua_Integer>(Id));
    return 1;
}

// ─────────────────────────────────────
/**
 * @brief Bind timer cancellation to Lua.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return One Lua return value: whether a timer was cancelled.
 *
 * @note Reads an integer identifier and pushes false for nonpositive or unknown identifiers.
 * @warning Raises a Lua error for a missing instance or invalid argument type.
 */
static int OpenScofoCancel(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    const lua_Integer Id = luaL_checkinteger(L, 1);
    lua_pushboolean(L, Id > 0 && self->CancelLuaCallback(static_cast<uint64_t>(Id)));
    return 1;
}

// ─────────────────────────────────────
static const luaL_Reg oscofo_funcs[] = {
    {"schedule", OpenScofoSchedule},
    {"cancel", OpenScofoCancel},
    {"activate_all_descriptors", OpenScofoActivateAllDescriptors},
    {"set_current_event", OpenScofoSetCurrentEvent},
    {"set_current_section", OpenScofoSetCurrentSection},
    {"get_live_bpm", OpenScofoGetLiveBPM},
    {"get_event_index", OpenScofoGetEventIndex},
    {"get_states", OpenScofoGetStates},
    {"get_audio_description", OpenScofoGetCurrentDescription},

    {NULL, NULL},
};

// ─────────────────────────────────────
/**
 * @brief Register the OpenScofo functions and return the Lua module table.
 *
 * @param L Lua state whose stack supplies arguments or receives results.
 *
 * @return One Lua return value: the OpenScofo module table.
 *
 * @note Reuses the _OpenScofo table when present and registers it in package.loaded when available.
 */
int luaopen_OpenScofo(lua_State *L) {
    lua_getglobal(L, "_OpenScofo");
    if (!lua_istable(L, -1)) {
        lua_pop(L, 1);
        lua_newtable(L);
        lua_pushvalue(L, -1);
        lua_setglobal(L, "_OpenScofo");
    }

    const int moduleIndex = lua_absindex(L, -1);
    luaL_setfuncs(L, oscofo_funcs, 0);

    lua_getglobal(L, "package");
    if (lua_istable(L, -1)) {
        lua_getfield(L, -1, "loaded");
        if (lua_istable(L, -1)) {
            lua_pushvalue(L, moduleIndex);
            lua_setfield(L, -2, "OpenScofo");
        }
        lua_pop(L, 1);
    }
    lua_pop(L, 1);

    return 1;
}

} // namespace OpenScofo

#endif
