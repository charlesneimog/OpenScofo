/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unordered_map.h>
#include <nanobind/stl/vector.h>

#include <OpenScofo.hpp>

namespace nb = nanobind;

template <typename F> static auto checked_call(F func) {
    try {
        auto result = func();
        if (PyErr_Occurred()) throw nb::python_error();
        return result;
    } catch (...) { if (PyErr_Occurred()) throw nb::python_error(); throw; }
}
template <typename T> static bool process_python_block(OpenScofo::OpenScofo &self, const T *data, size_t size) {
    if (size == 0) {
        return true;
    }

    const size_t hop = static_cast<size_t>(self.GetHopSize());
    assert(hop != 0);

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
    case spdlog::level::critical: {
        nb::gil_scoped_acquire gil;
        PyErr_SetString(PyExc_RuntimeError, text.c_str());
        break;
    }
    case spdlog::level::err: {
        nb::gil_scoped_acquire gil;
        PyErr_SetString(PyExc_RuntimeError, text.c_str());
        break;
    }
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
        .value("PTECH", OpenScofo::PTECH)
        .value("UTECH", OpenScofo::UTECH)
        .export_values();

    nb::enum_<OpenScofo::Descriptors>(m, "Descriptors")
        .value("INVALID", OpenScofo::INVALID)

        .value("ODSONSET", OpenScofo::ODSONSET)

        // Amplitude
        .value("LOUDNESS", OpenScofo::LOUDNESS)
        .value("DB", OpenScofo::DB)
        .value("MAXAMP", OpenScofo::MAXAMP)
        .value("RMS", OpenScofo::RMS)
        .value("STDDEV", OpenScofo::STDDEV)
        .value("MAGNITUDE", OpenScofo::MAGNITUDE)
        .value("POWERARRAY", OpenScofo::POWERARRAY)
        .value("SILENCEPROB", OpenScofo::SILENCEPROB)

        // Spectral Arrays
        .value("MFCC", OpenScofo::MFCC)
        .value("CHROMA", OpenScofo::CHROMA)
        .value("LOGMEL", OpenScofo::LOGMEL)

        // Spectral
        .value("ZCR", OpenScofo::ZCR)
        .value("HFR", OpenScofo::HFR)
        .value("CENTROID", OpenScofo::CENTROID)
        .value("SPREADHZ", OpenScofo::SPREADHZ)
        .value("SPREADVARIANCE", OpenScofo::SPREADVARIANCE)
        .value("CREST", OpenScofo::CREST)
        .value("FLATNESS", OpenScofo::FLATNESS)
        .value("ENTROPY", OpenScofo::ENTROPY)
        .value("ROLLOFF", OpenScofo::ROLLOFF)
        .value("CENTROIDVEL", OpenScofo::CENTROIDVEL)
        .value("FLUX", OpenScofo::FLUX)
        .value("SKEWNESS", OpenScofo::SKEWNESS)
        .value("SLOPE", OpenScofo::SLOPE)
        .value("KURTOSIS", OpenScofo::KURTOSIS)
        .value("IRREGULARITY", OpenScofo::IRREGULARITY)
        .value("HARMONICITY", OpenScofo::HARMONICITY)

        // Pitch
        .value("YIN", OpenScofo::YIN)
        .value("YINCONFIDENCE", OpenScofo::YINCONFIDENCE)

        // Percussive
        .value("EXTENDEDTECHNIQUE", OpenScofo::EXTENDEDTECHNIQUE)

        // AI
        .value("ONNX", OpenScofo::ONNX)
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
        .def_rw("irregularity_jensen", &OpenScofo::Description::SpectralIrregularityJensen)
        .def_rw("irregularity_krimphoff", &OpenScofo::Description::SpectralIrregularityKrimphoff)
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

        .def_rw("logmel", &OpenScofo::Description::LogMelSpectrum) // MELOGRAM

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

    nb::class_<OpenScofo::Configuration>(m, "Configuration")
        .def(nb::init<>())

        // Audio parameters
        .def_rw("sr", &OpenScofo::Configuration::SR)
        .def_rw("fft_size", &OpenScofo::Configuration::FFTSize)
        .def_rw("hop_size", &OpenScofo::Configuration::HOPSize)

        // Tuning
        .def_rw("tuning_a4", &OpenScofo::Configuration::TuningA4)

        // Pitch template
        .def_rw("pitch_template_sigma", &OpenScofo::Configuration::PitchTemplateSigma)
        .def_rw("pitch_template_harmonics", &OpenScofo::Configuration::PitchTemplateHarmonics)

        // MFCC
        .def_rw("mfcc_mels", &OpenScofo::Configuration::MFCCMels)
        .def_rw("mfcc_count", &OpenScofo::Configuration::MFCCCount)

        // Onset
        .def_rw("onset_type", &OpenScofo::Configuration::OnsetType)
        .def_rw("med_span", &OpenScofo::Configuration::MedSpan)

        // Silence
        .def_rw("db_threshold", &OpenScofo::Configuration::dBTreshold)

        // YIN / spectral
        .def_rw("spectral_rolloff_cutoff", &OpenScofo::Configuration::SpectralRolloffCutoff)
        .def_rw("yin_threshold", &OpenScofo::Configuration::YINThreshold)
        .def_rw("yin_min_freq", &OpenScofo::Configuration::YINMinFrequency)
        .def_rw("yin_max_freq", &OpenScofo::Configuration::YINMaxFrequency)

        // Chroma
        .def_rw("chroma_size", &OpenScofo::Configuration::ChromaSize)
        .def_rw("chroma_center_octave", &OpenScofo::Configuration::ChromaCenterOctave)
        .def_rw("chroma_octave_width", &OpenScofo::Configuration::ChromaOctaveWidth)

        // ZCR
        .def_rw("zcr_center", &OpenScofo::Configuration::ZCRCenter)
        .def_rw("zcr_pad", &OpenScofo::Configuration::ZCRPad)
        .def_rw("zcr_zero_pos", &OpenScofo::Configuration::ZCRZeroPos)
        .def_rw("zcr_threshold", &OpenScofo::Configuration::ZCRThreshold)

        // Temporal model
        .def_rw("sync_strength", &OpenScofo::Configuration::SyncStrength)
        .def_rw("phase_coupling", &OpenScofo::Configuration::PhaseCoupling)

        // ONNX
        .def_rw("audio_state_change_receiver", &OpenScofo::Configuration::AudioStateChangeReceiver)
        .def_rw("timbre_onnx_model", &OpenScofo::Configuration::TimbreONNXModel)
        .def_rw("onnx_descriptors", &OpenScofo::Configuration::ONNXDescriptors)
        .def_rw("requested_descriptors", &OpenScofo::Configuration::RequestedDescriptors);

    nb::class_<OpenScofo::OpenScofo>(m, "OpenScofo")
        .def("__init__",
             [](OpenScofo::OpenScofo *self, float sr, float fft_size, float hop) {
                 new (self) OpenScofo::OpenScofo(sr, fft_size, hop);
                 self->SetErrorCallback(python_error_callback);
             })

        .def("__repr__",
             [](OpenScofo::OpenScofo &self) {
                 std::ostringstream oss;
                 oss << "<OpenScofo(sr=" << self.GetSr() << ", fft_size=" << self.GetFFTSize()
                     << ", hop=" << self.GetHopSize() << ")>";

                 return oss.str();
             })

        // Score
        .def("load_score", [](OpenScofo::OpenScofo &self, const std::filesystem::path &path) { return checked_call([&] { return self.LoadScore(path); }); })
        .def("score_is_loaded", &OpenScofo::OpenScofo::ScoreIsLoaded)

        // Config
        .def("set_current_event", &OpenScofo::OpenScofo::SetCurrentEvent)

        // ONNX
        .def("load_onnx_model", &OpenScofo::OpenScofo::LoadONNXModel)

        // Setters
        .def("set_configuration", &OpenScofo::OpenScofo::SetConfiguration)
        .def("set_requested_descriptors", &OpenScofo::OpenScofo::SetRequestedDescriptors)
        .def("request_descriptor", &OpenScofo::OpenScofo::RequestDescriptor)

        // Getters
        .def("get_configuration", &OpenScofo::OpenScofo::GetConfiguration)
        .def("get_current_bpm", &OpenScofo::OpenScofo::GetCurrentBPM)
        .def("get_current_score_position", &OpenScofo::OpenScofo::GetCurrentScorePosition)
        .def("get_current_event_actions", &OpenScofo::OpenScofo::GetCurrentEventActions)
        .def("get_audio_state_change_actions", &OpenScofo::OpenScofo::GetAudioStateChangeActions)
        .def("get_lua_code", &OpenScofo::OpenScofo::GetLuaCode)
        .def("get_pitch_prob", &OpenScofo::OpenScofo::GetPitchProb)
        .def("get_states", &OpenScofo::OpenScofo::GetStates)
        .def("get_pitch_template", &OpenScofo::OpenScofo::GetPitchTemplate)
        .def("get_sr", &OpenScofo::OpenScofo::GetSr)
        .def("get_fft_size", &OpenScofo::OpenScofo::GetFFTSize)
        .def("get_hop_size", &OpenScofo::OpenScofo::GetHopSize)
        .def("get_block_duration", &OpenScofo::OpenScofo::GetBlockDuration)
        .def("get_description", &OpenScofo::OpenScofo::GetDescription)
        .def("get_current_buffer_index", &OpenScofo::OpenScofo::GetCurrentBufferIndex)

        // Process (template wrapper)
        .def("process_block",
             [](OpenScofo::OpenScofo &self, nb::ndarray<const float, nb::shape<-1>, nb::c_contig> audio) {
                 return checked_call([&] { return process_python_block<float>(self, audio.data(), audio.size()); });
             })
        .def("process_block",
             [](OpenScofo::OpenScofo &self, nb::ndarray<const double, nb::shape<-1>, nb::c_contig> audio) {
                 return checked_call([&] { return process_python_block<double>(self, audio.data(), audio.size()); });
             })

        // Logging
        .def("set_log_level", &OpenScofo::OpenScofo::SetLogLevel);
}
