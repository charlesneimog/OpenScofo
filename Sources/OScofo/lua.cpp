#include <OScofo.hpp>
#include <sol/sol.hpp>

#if defined(OSCOFO_LUA)

namespace OScofo {

int luaopen_oscofo(lua_State *L) {
    sol::state_view lua(L);
    sol::table m = lua.create_table();
    sol::usertype<OScofo> oscofo_type = lua.new_usertype<OScofo>("OScofo", sol::constructors<OScofo(float, float, float)>());

    // Error callback example
    oscofo_type.set(sol::call_constructor, [](float sr, float fft_size, float hop) {
        auto obj = new OScofo(sr, fft_size, hop);
        obj->SetErrorCallBack([](const std::string &msg) { throw std::runtime_error(msg); });
        return obj;
    });

    // Methods (same as pybind11)
    oscofo_type["set_db_threshold"] = &OScofo::SetdBTreshold;
    oscofo_type["set_tuning"] = &OScofo::SetTunning;
    oscofo_type["set_current_event"] = &OScofo::SetCurrentEvent;
    oscofo_type["set_amplitude_decay"] = &OScofo::SetAmplitudeDecay;
    oscofo_type["set_harmonics"] = &OScofo::SetHarmonics;
    oscofo_type["set_pitch_template_sigma"] = &OScofo::SetPitchTemplateSigma;
    oscofo_type["get_live_bpm"] = &OScofo::GetLiveBPM;
    oscofo_type["get_event_index"] = &OScofo::GetEventIndex;
    oscofo_type["get_error"] = &OScofo::GetErrorMessage;
    oscofo_type["get_states"] = &OScofo::GetStates;
    oscofo_type["get_pitch_template"] = &OScofo::GetPitchTemplate;
    oscofo_type["get_cqt_template"] = &OScofo::GetCQTTemplate;
    oscofo_type["get_audio_description"] = &OScofo::GetAudioDescription;
    oscofo_type["get_time_coherence_template"] = sol::overload(
        [](OScofo &self, int pos) { return self.GetTimeCoherenceTemplate(pos, 0); },
        [](OScofo &self, int pos, int time_in_event) { return self.GetTimeCoherenceTemplate(pos, time_in_event); }
    );
    oscofo_type["get_time_coherence_confiability"] = &OScofo::GetTimeCoherenceConfiability;

    m["OScofo"] = oscofo_type;

    // ─── Description class ───
    sol::usertype<Description> desc_type = lua.new_usertype<Description>("Description", sol::constructors<Description()>());
    desc_type["mfcc"] = &Description::MFCC;
    desc_type["onset"] = &Description::Onset;
    desc_type["silence_prob"] = &Description::SilenceProb;
    desc_type["spectral_power"] = &Description::SpectralPower;
    desc_type["norm_spectral_power"] = &Description::NormSpectralPower;
    desc_type["pseudo_cqt"] = &Description::PseudoCQT;
    desc_type["loudness"] = &Description::Loudness;
    desc_type["spectral_flux"] = &Description::SpectralFlux;
    desc_type["spectral_flatness"] = &Description::SpectralFlatness;
    desc_type["harmonicity"] = &Description::Harmonicity;
    desc_type["db"] = &Description::dB;
    desc_type["rms"] = &Description::RMS;
    desc_type["power"] = &Description::Power;
    m["Description"] = desc_type;

    // ─── MacroState / State class ───
    sol::usertype<MacroState> state_type = lua.new_usertype<MacroState>("State", sol::constructors<MacroState()>());
    state_type["position"] = &MacroState::ScorePos;
    state_type["type"] = &MacroState::Type;
    state_type["markov"] = &MacroState::Markov;
    state_type["forward"] = &MacroState::Forward;
    state_type["bpm_expected"] = &MacroState::BPMExpected;
    state_type["bpm_observed"] = &MacroState::BPMObserved;
    state_type["onset_expected"] = &MacroState::OnsetExpected;
    state_type["onset_observed"] = &MacroState::OnsetObserved;
    state_type["phase_expected"] = &MacroState::PhaseExpected;
    state_type["phase_observed"] = &MacroState::PhaseObserved;
    state_type["ioi_phi_n"] = &MacroState::IOIPhiN;
    state_type["ioi_hat_phi_n"] = &MacroState::IOIHatPhiN;
    state_type["audiostates"] = &MacroState::AudioStates;
    state_type["duration"] = &MacroState::Duration;
    state_type["line"] = &MacroState::Line;
    m["State"] = state_type;

    // ─── AudioState class ───
    sol::usertype<AudioState> audio_type = lua.new_usertype<AudioState>("AudioState", sol::constructors<AudioState()>());
    audio_type["freq"] = &AudioState::Freq;
    audio_type["index"] = &AudioState::Index;
    m["AudioState"] = audio_type;

    // ─── Push the module table ───
    m.push();
    return 1;
}

} // namespace OScofo

#endif
