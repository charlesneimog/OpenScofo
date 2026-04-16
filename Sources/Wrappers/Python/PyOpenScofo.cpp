#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/filesystem.h>

#include <OpenScofo.hpp>

namespace nb = nanobind;

template <typename T> static bool process_python_block(OpenScofo::OpenScofo &self, const T *data, size_t size) {
    if (size == 0) {
        return true;
    }

    const size_t hop = static_cast<size_t>(self.GetHopSize());
    if (hop == 0) {
        throw nb::value_error("OpenScofo hop size cannot be zero");
    }

    bool ok = true;
    size_t offset = 0;

    // Accept both hop-sized and FFT-sized buffers from Python by internally
    // chunking larger inputs into hop-sized processing steps.
    while (offset < size) {
        const size_t chunk = std::min(hop, size - offset);
        ok = self.ProcessBlock<T>(data + offset, chunk) && ok;
        offset += chunk;
    }

    return ok;
}

static void python_error_callback(const spdlog::details::log_msg &log, void *data) {
    (void)data;
    std::string text(log.payload.data(), log.payload.size());
    switch (log.level) {
    case spdlog::level::critical:
        throw nb::value_error(text.c_str());
    case spdlog::level::err:
        throw nb::value_error(text.c_str());
        break;
    case spdlog::level::info:
        nb::print(nb::str(text.c_str()));
        break;
    case spdlog::level::debug:
        nb::print(nb::str(("\033[90m" + text + "\033[0m").c_str()));
        break;
    default:
        break;
    }
}

NB_MODULE(_OpenScofo, m) {

    m.doc() = "OpenScofo Python bindings";
    m.attr("__version__") = OPENSCOFO_VERSION;

    nb::enum_<OpenScofo::EventType>(m, "EventType")
        .value("REST", OpenScofo::REST)
        .value("NOTE", OpenScofo::NOTE)
        .value("CHORD", OpenScofo::CHORD)
        .value("TRILL", OpenScofo::TRILL)
        .value("MULTI", OpenScofo::MULTI)
        .export_values();

    nb::enum_<OpenScofo::HMMType>(m, "HMMType")
        .value("SEMIMARKOV", OpenScofo::SEMIMARKOV)
        .value("MARKOV", OpenScofo::MARKOV)
        .export_values();

    nb::class_<OpenScofo::AudioState>(m, "AudioState")
        .def(nb::init<>())
        .def_rw("type", &OpenScofo::AudioState::Type)
        .def_rw("frequency", &OpenScofo::AudioState::Freq)
        .def_rw("midi", &OpenScofo::AudioState::Midi)
        .def_rw("label", &OpenScofo::AudioState::Label)
        .def_rw("index", &OpenScofo::AudioState::Index);

    // Description Class
    nb::class_<OpenScofo::Description>(m, "Description")
        .def(nb::init<>())
        // Onset
        .def_rw("onset", &OpenScofo::Description::Onset)

        // Probability
        .def_rw("silence", &OpenScofo::Description::SilenceProb)  // SILENCEPROB
        .def_rw("ext", &OpenScofo::Description::ExtendedTechProb) // EXTENDEDPROB
        // .def_rw("window_last_onset", &OpenScofo::Description::WindowLastOnset)

        // Amplitude
        .def_rw("db", &OpenScofo::Description::dB)
        .def_rw("rms", &OpenScofo::Description::RMS) // RMS
        .def_rw("max_amp", &OpenScofo::Description::MaxAmp)
        .def_rw("loudness", &OpenScofo::Description::Loudness) // LOUDNESS

        // Spectral (following your mapping)
        .def_rw("flatness", &OpenScofo::Description::SpectralFlatness) // FLATNESS
        .def_rw("entropy", &OpenScofo::Description::SpectralEntropy)   // ENTROPY
        .def_rw("rolloff", &OpenScofo::Description::SpectralRolloff)   // ROLLOFF
        .def_rw("flux", &OpenScofo::Description::SpectralFlux)         // FLUX
        .def_rw("skewness", &OpenScofo::Description::SpectralSkewness) // SKEWNESS
        .def_rw("slope", &OpenScofo::Description::SpectralSlope)       // SLOPE
        .def_rw("kurtosis", &OpenScofo::Description::SpectralKurtosis) // KURTOSIS
        .def_rw("centroid", &OpenScofo::Description::SpectralCentroid) // CENTROID
        .def_rw("spreadhz", &OpenScofo::Description::SpectralSpreadHz) // SPREADHZ
        .def_rw("zcr", &OpenScofo::Description::ZeroCrossingRate)      // ZCR
        .def_rw("hfr", &OpenScofo::Description::HighFreqRatio)         // HFR

        // Other spectral-related (no enum mapping provided → keep concise names)
        .def_rw("harmonicity", &OpenScofo::Description::Harmonicity)
        .def_rw("irregularity", &OpenScofo::Description::SpectralIrregularity)
        .def_rw("crest", &OpenScofo::Description::SpectralCrest)
        .def_rw("centroid_velocity", &OpenScofo::Description::CentroidVelocity)
        .def_rw("spread_variance", &OpenScofo::Description::SpectralSpreadVariance)

        // Pitch
        .def_rw("yin", &OpenScofo::Description::Pitch) // YIN
        .def_rw("pitch", &OpenScofo::Description::Pitch)
        .def_rw("pitch_confidence", &OpenScofo::Description::PitchConfidence)

        // Stats
        .def_rw("stddev", &OpenScofo::Description::StdDev) // STDDEV

        // Arrays
        .def_rw("power", &OpenScofo::Description::Power)
        .def_rw("magnitude", &OpenScofo::Description::Magnitude)
        // .def_rw("spectral_magnitude_norm", &OpenScofo::Description::SpectralMagnitudeNorm)
        // .def_rw("spectral_magnitude_frame_norm", &OpenScofo::Description::SpectralMagnitudeFrameNorm)

        .def_rw("logmelspectrum", &OpenScofo::Description::LogMelSpectrum) // MELOGRAM
        // .def_rw("reverb_spectral_power", &OpenScofo::Description::ReverbSpectralPower)

        .def_rw("mfcc", &OpenScofo::Description::MFCC)     // MFCC
        .def_rw("chroma", &OpenScofo::Description::Chroma) // CHROMA

        // ONNX
        .def_rw("onnx", &OpenScofo::Description::ONNX); // ONNX

    // State Class

    nb::class_<OpenScofo::MarkovState>(m, "MarkovState")
        .def(nb::init<>())

        // Core
        .def_rw("index", &OpenScofo::MarkovState::Index)
        .def_rw("score_pos", &OpenScofo::MarkovState::ScorePos)
        .def_rw("markov_index", &OpenScofo::MarkovState::MarkovIndex)
        .def_rw("audio_states", &OpenScofo::MarkovState::AudioStates)

        // State actions
        .def_rw("hsmm_type", &OpenScofo::MarkovState::HSMMType)
        .def_rw("type", &OpenScofo::MarkovState::Type)
        .def_rw("actions", &OpenScofo::MarkovState::Actions)

        // Inference
        .def_rw("init_prob", &OpenScofo::MarkovState::InitProb)
        .def_rw("forward", &OpenScofo::MarkovState::Forward)
        .def_rw("exit_prob", &OpenScofo::MarkovState::ExitProb)
        .def_rw("best_obs", &OpenScofo::MarkovState::BestObs)

        // Time
        .def_rw("upper_bound", &OpenScofo::MarkovState::UpperBound)
        .def_rw("bpm_expected", &OpenScofo::MarkovState::BPMExpected)
        .def_rw("bpm_observed", &OpenScofo::MarkovState::BPMObserved)
        .def_rw("onset_expected", &OpenScofo::MarkovState::OnsetExpected)
        .def_rw("onset_observed", &OpenScofo::MarkovState::OnsetObserved)
        .def_rw("phase_expected", &OpenScofo::MarkovState::PhaseExpected)
        .def_rw("phase_observed", &OpenScofo::MarkovState::PhaseObserved)
        .def_rw("ioi_phi_n", &OpenScofo::MarkovState::IOIPhiN)
        .def_rw("ioi_hat_phi_n", &OpenScofo::MarkovState::IOIHatPhiN)
        .def_rw("duration", &OpenScofo::MarkovState::Duration)

        // Configuration
        .def_rw("phase_coupling", &OpenScofo::MarkovState::PhaseCoupling)
        .def_rw("sync_strength", &OpenScofo::MarkovState::SyncStrength)

        // Error
        .def_rw("line", &OpenScofo::MarkovState::Line);

    nb::class_<OpenScofo::OpenScofo>(m, "OpenScofo")
        .def("__init__",
             [](OpenScofo::OpenScofo *self, float sr, float fft_size, float hop) {
                 new (self) OpenScofo::OpenScofo(sr, fft_size, hop);
                 self->SetErrorCallback(python_error_callback);
             })

        .def("__repr__",
             [](OpenScofo::OpenScofo &self) {
                 return std::format("<OpenScofo(sr={}, fft_size={}, hop={})>", self.GetSr(), self.GetFFTSize(),
                                    self.GetHopSize());
             })

        // Score
        .def("parse_score", &OpenScofo::OpenScofo::ParseScore)
        .def("score_is_loaded", &OpenScofo::OpenScofo::ScoreIsLoaded)

        // Config
        .def("set_db_threshold", &OpenScofo::OpenScofo::SetdBTreshold)
        .def("set_tuning", &OpenScofo::OpenScofo::SetTunning)
        .def("set_current_event", &OpenScofo::OpenScofo::SetCurrentEvent)

        // ONNX
        .def("load_onnx_model", &OpenScofo::OpenScofo::LoadONNXModel)

        // Pitch template
        .def("set_amplitude_decay", &OpenScofo::OpenScofo::SetAmplitudeDecay)
        .def("set_harmonics", &OpenScofo::OpenScofo::SetHarmonics)
        .def("set_pitch_template_sigma", &OpenScofo::OpenScofo::SetPitchTemplateSigma)

        // Getters
        .def("get_live_bpm", &OpenScofo::OpenScofo::GetLiveBPM)
        .def("get_current_score_position", &OpenScofo::OpenScofo::GetCurrentScorePosition)
        .def("get_kappa", &OpenScofo::OpenScofo::GetKappa)
        .def("get_db_value", &OpenScofo::OpenScofo::GetdBValue)
        .def("get_event_actions", &OpenScofo::OpenScofo::GetEventActions)
        .def("get_lua_code", &OpenScofo::OpenScofo::GetLuaCode)
        .def("get_pitch_prob", &OpenScofo::OpenScofo::GetPitchProb)
        .def("get_states", &OpenScofo::OpenScofo::GetStates)
        .def("get_pitch_template", &OpenScofo::OpenScofo::GetPitchTemplate)
        .def("get_spectrum_power", &OpenScofo::OpenScofo::GetSpectrumPower)
        .def("get_sr", &OpenScofo::OpenScofo::GetSr)
        .def("get_fft_size", &OpenScofo::OpenScofo::GetFFTSize)
        .def("get_hop_size", &OpenScofo::OpenScofo::GetHopSize)
        .def("get_block_duration", &OpenScofo::OpenScofo::GetBlockDuration)
        .def("get_description", &OpenScofo::OpenScofo::GetDescription)
        .def("get_current_buffer_index", &OpenScofo::OpenScofo::GetCurrentBufferIndex)

        // Process (template wrapper)
        .def("process_block",
             [](OpenScofo::OpenScofo &self, nb::ndarray<const float, nb::shape<-1>, nb::c_contig> audio) {
                 return process_python_block<float>(self, audio.data(), audio.size());
             })
        .def("process_block",
             [](OpenScofo::OpenScofo &self, nb::ndarray<const double, nb::shape<-1>, nb::c_contig> audio) {
                 return process_python_block<double>(self, audio.data(), audio.size());
             })

        // Logging
        .def("set_log_level", &OpenScofo::OpenScofo::SetLogLevel);
}
