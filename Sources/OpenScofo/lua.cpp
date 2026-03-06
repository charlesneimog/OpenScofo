#include <OpenScofo.hpp>

#if defined(OSCOFO_LUA)

namespace OpenScofo {

// ─────────────────────────────────────
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
static void PushNumberVector(lua_State *L, const std::vector<double> &values) {
    lua_createtable(L, static_cast<int>(values.size()), 0);
    for (size_t i = 0; i < values.size(); ++i) {
        lua_pushnumber(L, values[i]);
        lua_rawseti(L, -2, static_cast<int>(i + 1));
    }
}

// ─────────────────────────────────────
static void PushAudioState(lua_State *L, const AudioState &state) {
    lua_createtable(L, 0, 4);
    lua_pushinteger(L, state.Type);
    lua_setfield(L, -2, "type");
    lua_pushnumber(L, state.Freq);
    lua_setfield(L, -2, "freq");
    lua_pushnumber(L, state.Midi);
    lua_setfield(L, -2, "midi");
    lua_pushinteger(L, static_cast<lua_Integer>(state.Index));
    lua_setfield(L, -2, "index");
}

// // ─────────────────────────────────────
// static void PushDescription(lua_State *L, const Description &desc) {
//     lua_createtable(L, 0, 16);
//
//     lua_pushboolean(L, desc.Silence);
//     lua_setfield(L, -2, "silence");
//     lua_pushboolean(L, desc.Onset);
//     lua_setfield(L, -2, "onset");
//     lua_pushnumber(L, desc.SilenceProb);
//     lua_setfield(L, -2, "silence_prob");
//
//     lua_pushnumber(L, desc.dB);
//     lua_setfield(L, -2, "db");
//     lua_pushnumber(L, desc.RMS);
//     lua_setfield(L, -2, "rms");
//     lua_pushnumber(L, desc.MaxAmp);
//     lua_setfield(L, -2, "max_amp");
//     lua_pushnumber(L, desc.Loudness);
//     lua_setfield(L, -2, "loudness");
//
//     lua_pushnumber(L, desc.Harmonicity);
//     lua_setfield(L, -2, "harmonicity");
//     lua_pushnumber(L, desc.SpectralFlatness);
//     lua_setfield(L, -2, "spectral_flatness");
//     lua_pushnumber(L, desc.SpectralFlux);
//     lua_setfield(L, -2, "spectral_flux");
//     lua_pushnumber(L, desc.StdDev);
//     lua_setfield(L, -2, "std_dev");
//
//     PushNumberVector(L, desc.Power);
//     lua_setfield(L, -2, "power");
//     PushNumberVector(L, desc.SpectralPower);
//     lua_setfield(L, -2, "spectral_power");
//     PushNumberVector(L, desc.NormSpectralPower);
//     lua_setfield(L, -2, "norm_spectral_power");
//     PushNumberVector(L, desc.ReverbSpectralPower);
//     lua_setfield(L, -2, "reverb_spectral_power");
//     PushNumberVector(L, desc.PseudoCQT);
//     lua_setfield(L, -2, "pseudo_cqt");
//     PushNumberVector(L, desc.MFCC);
//     lua_setfield(L, -2, "mfcc");
//     PushNumberVector(L, desc.Chroma);
//     lua_setfield(L, -2, "chroma");
// }

// ─────────────────────────────────────
static void PushMarkovState(lua_State *L, const MarkovState &state) {
    lua_createtable(L, 0, 14);
    lua_pushinteger(L, state.ScorePos);
    lua_setfield(L, -2, "position");
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

    lua_createtable(L, static_cast<int>(state.AudioStates.size()), 0);
    for (size_t i = 0; i < state.AudioStates.size(); ++i) {
        PushAudioState(L, state.AudioStates[i]);
        lua_rawseti(L, -2, static_cast<int>(i + 1));
    }
    lua_setfield(L, -2, "audiostates");
}

// ─────────────────────────────────────
static int OpenScofoSetDbThreshold(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    self->SetdBTreshold(luaL_checknumber(L, 1));
    return 0;
}

// ─────────────────────────────────────
static int OpenScofoSetTuning(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    self->SetTunning(luaL_checknumber(L, 1));
    return 0;
}

// ─────────────────────────────────────
static int OpenScofoSetCurrentEvent(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    self->SetCurrentEvent(static_cast<int>(luaL_checkinteger(L, 1)));
    return 0;
}

// ─────────────────────────────────────
static int OpenScofoSetAmplitudeDecay(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    self->SetAmplitudeDecay(luaL_checknumber(L, 1));
    return 0;
}

// ─────────────────────────────────────
static int OpenScofoSetHarmonics(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    self->SetHarmonics(static_cast<int>(luaL_checkinteger(L, 1)));
    return 0;
}

// ─────────────────────────────────────
static int OpenScofoSetPitchTemplateSigma(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    self->SetPitchTemplateSigma(luaL_checknumber(L, 1));
    return 0;
}

// ─────────────────────────────────────
static int OpenScofoGetLiveBPM(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    lua_pushnumber(L, self->GetLiveBPM());
    return 1;
}

// ─────────────────────────────────────
static int OpenScofoGetEventIndex(lua_State *L) {
    OpenScofo *self = GetCurrentOpenScofo(L);
    if (self == nullptr)
        return luaL_error(L, "OpenScofo pointer is null");
    lua_pushinteger(L, self->GetEventIndex());
    return 1;
}

// ─────────────────────────────────────
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
static const luaL_Reg oscofo_funcs[] = {
    {"set_db_threshold", OpenScofoSetDbThreshold},
    {"set_tuning", OpenScofoSetTuning},
    {"set_current_event", OpenScofoSetCurrentEvent},
    {"set_amplitude_decay", OpenScofoSetAmplitudeDecay},
    {"set_harmonics", OpenScofoSetHarmonics},
    {"set_pitch_template_sigma", OpenScofoSetPitchTemplateSigma},
    {"get_live_bpm", OpenScofoGetLiveBPM},
    {"get_event_index", OpenScofoGetEventIndex},
    {"get_states", OpenScofoGetStates},
    {NULL, NULL},
};

// ─────────────────────────────────────
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
