/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

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
    enum_<OpenScofo::AudioDescType>("AudioDescType")
        .value("PITCH", OpenScofo::PITCH)
        .value("SILENCE", OpenScofo::SILENCE)
        .value("LABEL", OpenScofo::LABEL)
        .value("ONSET", OpenScofo::ONSET);
    class_<OpenScofo::Observation>("AudioState")
        .constructor<>()
        .property("type", &OpenScofo::Observation::Type)
        .property("frequency", &OpenScofo::Observation::Freq)
        .property("midi", &OpenScofo::Observation::Midi)
        .property("label", &OpenScofo::Observation::Label)
        .property("index", &OpenScofo::Observation::Index);

    enum_<OpenScofo::MicroTopology>("MicroTopology")
        .value("NO_MICROSTATES", OpenScofo::NO_MICROSTATES)
        .value("UNORDERED", OpenScofo::UNORDERED)
        .value("LEFT_RIGHT", OpenScofo::LEFT_RIGHT);

    register_vector<OpenScofo::Observation>("VectorObservation");
    class_<OpenScofo::MarkovMicroState>("MarkovMicroState")
        .constructor<>()
        .property("observations", &OpenScofo::MarkovMicroState::Observations)
        .property("duration_weight", &OpenScofo::MarkovMicroState::DurationWeight);
    register_vector<OpenScofo::MarkovMicroState>("VectorMarkovMicroState");

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
        .property("irregularity_jensen", &OpenScofo::Description::SpectralIrregularityJensen)
        .property("irregularity_krimphoff", &OpenScofo::Description::SpectralIrregularityKrimphoff)
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
    class_<OpenScofo::ScoreState>("MarkovState")
        .constructor<>()
        .property("position", &OpenScofo::ScoreState::ScorePos)
        .property("section", &OpenScofo::ScoreState::Section)
        .property("type", &OpenScofo::ScoreState::Type)
        .property("markov", &OpenScofo::ScoreState::HSMMType)
        .property("forward", &OpenScofo::ScoreState::Forward)
        .property("bpm_expected", &OpenScofo::ScoreState::BPMExpected)
        .property("bpm_observed", &OpenScofo::ScoreState::BPMObserved)
        .property("onset_expected", &OpenScofo::ScoreState::OnsetExpected)
        .property("onset_observed", &OpenScofo::ScoreState::OnsetObserved)
        .property("phase_expected", &OpenScofo::ScoreState::PhaseExpected)
        .property("phase_observed", &OpenScofo::ScoreState::PhaseObserved)
        .property("ioi_phi_n", &OpenScofo::ScoreState::IOIPhiN)
        .property("ioi_hat_phi_n", &OpenScofo::ScoreState::IOIHatPhiN)
        .property("audio_states", &OpenScofo::ScoreState::Observations)
        .property("microstates", &OpenScofo::ScoreState::MicroStates)
        .property("micro_topology", &OpenScofo::ScoreState::MicroTopologyType)
        .property("best_microstate_index", &OpenScofo::ScoreState::BestMicroStateIndex)
        .property("duration", &OpenScofo::ScoreState::Duration)
        .property("line", &OpenScofo::ScoreState::Line);

    // OpenScofo class
    class_<OpenScofo::OpenScofo>("OpenScofo")
        .constructor<float, float, float>()
        .function("load_score", &OpenScofo::OpenScofo::LoadScore)
        .function("activate_all_descriptors", &OpenScofo::OpenScofo::ActivateAllDescriptors)
        .function("set_current_section", &OpenScofo::OpenScofo::SetCurrentSection)
        .function("get_current_bpm", &OpenScofo::OpenScofo::GetCurrentBPM)
        .function("get_current_score_position", &OpenScofo::OpenScofo::GetCurrentScorePosition)
        .function("get_states", &OpenScofo::OpenScofo::GetStates)
        .function(
            "process_block", +[](OpenScofo::OpenScofo &self, emscripten::val input) {
                std::vector<double> vec = emscripten::vecFromJSArray<double>(input);
                self.ProcessBlock(vec.data(), vec.size());
            });
}
