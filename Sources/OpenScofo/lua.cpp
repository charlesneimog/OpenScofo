#include <OpenScofo.hpp>
#include <sol/sol.hpp>

#if defined(OSCOFO_LUA)

namespace OpenScofo {

int luaopen_oscofo(lua_State *L) {
    sol::state_view lua(L);
    sol::table m = lua.create_table();
    sol::usertype<OpenScofo> oscofo_type = lua.new_usertype<OpenScofo>("OpenScofo", sol::constructors<OpenScofo(float, float, float)>());

    // Error callback example
    oscofo_type.set(sol::call_constructor, [](float sr, float fft_size, float hop) {
        auto obj = new OpenScofo(sr, fft_size, hop);
        obj->SetErrorCallBack([](const std::string &msg) { throw std::runtime_error(msg); });
        return obj;
    });

    // Methods (same as pybind11)
    oscofo_type["set_db_threshold"] = &OpenScofo::SetdBTreshold;
    oscofo_type["set_tuning"] = &OpenScofo::SetTunning;
    oscofo_type["set_current_event"] = &OpenScofo::SetCurrentEvent;
    oscofo_type["set_amplitude_decay"] = &OpenScofo::SetAmplitudeDecay;
    oscofo_type["set_harmonics"] = &OpenScofo::SetHarmonics;
    oscofo_type["set_pitch_template_sigma"] = &OpenScofo::SetPitchTemplateSigma;
    oscofo_type["get_live_bpm"] = &OpenScofo::GetLiveBPM;
    oscofo_type["get_event_index"] = &OpenScofo::GetEventIndex;
    oscofo_type["get_error"] = &OpenScofo::GetErrorMessage;
    oscofo_type["get_states"] = &OpenScofo::GetStates;
    oscofo_type["get_pitch_template"] = &OpenScofo::GetPitchTemplate;
    oscofo_type["get_cqt_template"] = &OpenScofo::GetCQTTemplate;
    oscofo_type["get_audio_description"] = &OpenScofo::GetAudioDescription;
    oscofo_type["get_time_coherence_template"] = sol::overload(
        [](OpenScofo &self, int pos) { return self.GetTimeCoherenceTemplate(pos, 0); },
        [](OpenScofo &self, int pos, int time_in_event) { return self.GetTimeCoherenceTemplate(pos, time_in_event); }
    );
    oscofo_type["get_time_coherence_confiability"] = &OpenScofo::GetTimeCoherenceConfiability;

    m["OpenScofo"] = oscofo_type;

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

    // ─── MarkovState / State class ───
    sol::usertype<MarkovState> state_type = lua.new_usertype<MarkovState>("State", sol::constructors<MarkovState()>());
    state_type["position"] = &MarkovState::ScorePos;
    state_type["type"] = &MarkovState::Type;
    state_type["markov"] = &MarkovState::HSMMType;
    state_type["forward"] = &MarkovState::Forward;
    state_type["bpm_expected"] = &MarkovState::BPMExpected;
    state_type["bpm_observed"] = &MarkovState::BPMObserved;
    state_type["onset_expected"] = &MarkovState::OnsetExpected;
    state_type["onset_observed"] = &MarkovState::OnsetObserved;
    state_type["phase_expected"] = &MarkovState::PhaseExpected;
    state_type["phase_observed"] = &MarkovState::PhaseObserved;
    state_type["ioi_phi_n"] = &MarkovState::IOIPhiN;
    state_type["ioi_hat_phi_n"] = &MarkovState::IOIHatPhiN;
    state_type["audiostates"] = &MarkovState::AudioStates;
    state_type["duration"] = &MarkovState::Duration;
    state_type["line"] = &MarkovState::Line;
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

} // namespace OpenScofo

#endif
