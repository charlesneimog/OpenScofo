#include <emscripten/bind.h>
#include <OpenScofo.hpp>
#include <vector>
#include <string>

using namespace emscripten;

// Error callback wrapper for JS console
void js_error_callback(const spdlog::details::log_msg &log, void *data) {
    (void)data;
    std::string text(log.payload.data(), log.payload.size());
    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
        // In JS, throw an exception
        throw std::runtime_error(text);
    case spdlog::level::info:
    case spdlog::level::debug:
        printf("%s\n", text.c_str());
        break;
    default:
        break;
    }
}

EMSCRIPTEN_BINDINGS(OpenScofo_module) {

    // Bind std::vector<double> explicitly
    register_vector<double>("VectorDouble");

    // Enums
    enum_<OpenScofo::EventType>("EventType")
        .value("REST", OpenScofo::REST)
        .value("NOTE", OpenScofo::NOTE)
        .value("CHORD", OpenScofo::CHORD)
        .value("TRILL", OpenScofo::TRILL)
        .value("MULTI", OpenScofo::MULTI);

    enum_<OpenScofo::HMMType>("HMMType").value("SEMIMARKOV", OpenScofo::SEMIMARKOV).value("MARKOV", OpenScofo::MARKOV);

    // AudioState
    class_<OpenScofo::AudioState>("AudioState")
        .constructor<>()
        .property("frequency", &OpenScofo::AudioState::Freq)
        .property("index", &OpenScofo::AudioState::Index);

    // Description
    class_<OpenScofo::Description>("Description")
        .constructor<>()
        .property("mfcc", &OpenScofo::Description::MFCC)
        .property("logmelspectrum", &OpenScofo::Description::LogMelSpectrum)
        .property("chroma", &OpenScofo::Description::Chroma)
        .property("onset", &OpenScofo::Description::Onset)
        .property("silence", &OpenScofo::Description::SilenceProb)
        .property("ext", &OpenScofo::Description::ExtendedTechProb)
        .property("max_amp", &OpenScofo::Description::MaxAmp)
        .property("loudness", &OpenScofo::Description::Loudness)
        .property("flux", &OpenScofo::Description::SpectralFlux)
        .property("irregularity", &OpenScofo::Description::SpectralIrregularity)
        .property("crest", &OpenScofo::Description::SpectralCrest)
        .property("centroid", &OpenScofo::Description::SpectralCentroid)
        .property("magnitude", &OpenScofo::Description::SpectralMagnitudeNorm)
        .property("spreadhz", &OpenScofo::Description::SpectralSpreadHz)
        .property("spread_variance", &OpenScofo::Description::SpectralSpreadVariance)
        .property("flatness", &OpenScofo::Description::SpectralFlatness)
        .property("skewness", &OpenScofo::Description::SpectralSkewness)
        .property("slope", &OpenScofo::Description::SpectralSlope)
        .property("kurtosis", &OpenScofo::Description::SpectralKurtosis)
        .property("centroid_velocity", &OpenScofo::Description::CentroidVelocity)
        .property("hfr", &OpenScofo::Description::HighFreqRatio)
        .property("harmonicity", &OpenScofo::Description::Harmonicity)
        .property("zcr", &OpenScofo::Description::ZeroCrossingRate)
        .property("stddev", &OpenScofo::Description::StdDev)
        .property("pitch", &OpenScofo::Description::Pitch)
        .property("pitch_confidence", &OpenScofo::Description::PitchConfidence)
        .property("db", &OpenScofo::Description::dB)
        .property("rms", &OpenScofo::Description::RMS)
        .property("magnitude", &OpenScofo::Description::Magnitude);

    // MarkovState
    class_<OpenScofo::MarkovState>("MarkovState")
        .constructor<>()
        .property("position", &OpenScofo::MarkovState::ScorePos)
        .property("type", &OpenScofo::MarkovState::Type)
        .property("markov", &OpenScofo::MarkovState::HSMMType)
        .property("forward", &OpenScofo::MarkovState::Forward)
        .property("bpm_expected", &OpenScofo::MarkovState::BPMExpected)
        .property("bpm_observed", &OpenScofo::MarkovState::BPMObserved)
        .property("onset_expected", &OpenScofo::MarkovState::OnsetExpected)
        .property("onset_observed", &OpenScofo::MarkovState::OnsetObserved)
        .property("phase_expected", &OpenScofo::MarkovState::PhaseExpected)
        .property("phase_observed", &OpenScofo::MarkovState::PhaseObserved)
        .property("ioi_phi_n", &OpenScofo::MarkovState::IOIPhiN)
        .property("ioi_hat_phi_n", &OpenScofo::MarkovState::IOIHatPhiN)
        .property("audio_states", &OpenScofo::MarkovState::AudioStates)
        .property("duration", &OpenScofo::MarkovState::Duration)
        .property("line", &OpenScofo::MarkovState::Line);

    // OpenScofo class
    class_<OpenScofo::OpenScofo>("OpenScofo")
        .constructor<float, float, float>()
        .function("parse_score", &OpenScofo::OpenScofo::ParseScore)
        .function("set_db_threshold", &OpenScofo::OpenScofo::SetdBTreshold)
        .function("set_tuning", &OpenScofo::OpenScofo::SetTunning)
        .function("set_current_event", &OpenScofo::OpenScofo::SetCurrentEvent)
        .function("load_onnx_model", &OpenScofo::OpenScofo::LoadONNXModel)
        .function("set_amplitude_decay", &OpenScofo::OpenScofo::SetAmplitudeDecay)
        .function("set_harmonics", &OpenScofo::OpenScofo::SetHarmonics)
        .function("set_pitch_template_sigma", &OpenScofo::OpenScofo::SetPitchTemplateSigma)
        .function("get_live_bpm", &OpenScofo::OpenScofo::GetLiveBPM)
        .function("get_event_index", &OpenScofo::OpenScofo::GetEventIndex)
        .function("get_states", &OpenScofo::OpenScofo::GetStates)
        .function("get_pitch_template", &OpenScofo::OpenScofo::GetPitchTemplate)
        .function("get_block_duration", &OpenScofo::OpenScofo::GetBlockDuration)
        .function("process_block", &OpenScofo::OpenScofo::ProcessBlock);
}
