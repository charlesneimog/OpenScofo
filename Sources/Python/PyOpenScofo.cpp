#include <OpenScofo.hpp>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

static void python_error_callback(const spdlog::details::log_msg &log, void *data) {
    (void)data;
    std::string text(log.payload.data(), log.payload.size());
    switch (log.level) {
    case spdlog::level::critical:
        throw py::value_error(text);
    case spdlog::level::err:
        throw py::value_error(text);
        break;
    case spdlog::level::info:
        py::print(text);
        break;
    case spdlog::level::debug:
        py::print("\033[90m" + text + "\033[0m");
        break;
    default:
        break;
    }
}

PYBIND11_MODULE(_OpenScofo, m) {

    py::class_<OpenScofo::OpenScofo>(m, "OpenScofo")
        .def(py::init([](float sr, float fft_size, float hop) {
            auto obj = new OpenScofo::OpenScofo(sr, fft_size, hop);
            obj->SetErrorCallback(python_error_callback);
            return obj;
        }))

        // Python
        .def("__repr__",
             [](OpenScofo::OpenScofo &self) {
                 return std::format("<OpenScofo(sr={}, fft_size={}, hop={})>", self.GetSr(), self.GetFFTSize(),
                                    self.GetHopSize());
             })

        // Score
        .def("parse_score", &OpenScofo::OpenScofo::ParseScore)

        // Config
        .def("set_db_threshold", &OpenScofo::OpenScofo::SetdBTreshold)
        .def("set_tuning", &OpenScofo::OpenScofo::SetTunning)
        .def("set_current_event", &OpenScofo::OpenScofo::SetCurrentEvent)

        // pitch template
        .def("set_amplitude_decay", &OpenScofo::OpenScofo::SetAmplitudeDecay)
        .def("set_harmonics", &OpenScofo::OpenScofo::SetHarmonics)
        .def("set_pitch_template_sigma", &OpenScofo::OpenScofo::SetPitchTemplateSigma)

        // Get Info
        .def("get_live_bpm", &OpenScofo::OpenScofo::GetLiveBPM)
        .def("get_event_index", &OpenScofo::OpenScofo::GetEventIndex)
        .def("get_states", &OpenScofo::OpenScofo::GetStates)
        .def("get_pitch_template", &OpenScofo::OpenScofo::GetPitchTemplate)
        .def("get_cqt_template", &OpenScofo::OpenScofo::GetCQTTemplate)
        .def("get_block_duration", &OpenScofo::OpenScofo::GetBlockDuration)

        // Help & Test Functions
        .def("get_audio_description", &OpenScofo::OpenScofo::GetAudioDescription)

        // Process
        .def("process_block", [](OpenScofo::OpenScofo &self, py::array_t<double> audio) {
            py::buffer_info bufInfo = audio.request();
            if (bufInfo.ndim != 1) {
                throw std::runtime_error("Input array must be 1-dimensional");
            }
            std::vector<double> cpp_audio(static_cast<double *>(bufInfo.ptr),
                                          static_cast<double *>(bufInfo.ptr) + bufInfo.shape[0]);
            return self.ProcessBlock(cpp_audio);
        });

    // Description Class
    py::class_<OpenScofo::Description>(m, "Description")
        .def(py::init<>())
        .def_readwrite("mfcc", &OpenScofo::Description::MFCC)
        .def_readwrite("onset", &OpenScofo::Description::Onset)
        .def_readwrite("silence_prob", &OpenScofo::Description::SilenceProb)
        .def_readwrite("spectral_power", &OpenScofo::Description::SpectralPower)
        .def_readwrite("norm_spectral_power", &OpenScofo::Description::NormSpectralPower)
        .def_readwrite("pseudo_cqt", &OpenScofo::Description::PseudoCQT)

        .def_readwrite("loudness", &OpenScofo::Description::Loudness)
        .def_readwrite("spectral_flux", &OpenScofo::Description::SpectralFlux)
        .def_readwrite("spectral_flatness", &OpenScofo::Description::SpectralFlatness)
        .def_readwrite("harmonicity", &OpenScofo::Description::Harmonicity)

        .def_readwrite("db", &OpenScofo::Description::dB)
        .def_readwrite("rms", &OpenScofo::Description::RMS)

        .def_readwrite("power", &OpenScofo::Description::Power);

    // State Class
    py::class_<OpenScofo::MarkovState>(m, "State")
        .def(py::init<>())
        .def_readwrite("position", &OpenScofo::MarkovState::ScorePos)
        .def_readwrite("type", &OpenScofo::MarkovState::Type)
        .def_readwrite("markov", &OpenScofo::MarkovState::HSMMType)
        .def_readwrite("forward", &OpenScofo::MarkovState::Forward)
        .def_readwrite("bpm_expected", &OpenScofo::MarkovState::BPMExpected)
        .def_readwrite("bpm_observed", &OpenScofo::MarkovState::BPMObserved)
        .def_readwrite("onset_expected", &OpenScofo::MarkovState::OnsetExpected)
        .def_readwrite("onset_observed", &OpenScofo::MarkovState::OnsetObserved)
        .def_readwrite("phase_expected", &OpenScofo::MarkovState::PhaseExpected)
        .def_readwrite("phase_observed", &OpenScofo::MarkovState::PhaseObserved)
        .def_readwrite("ioi_phi_n", &OpenScofo::MarkovState::IOIPhiN)
        .def_readwrite("ioi_hat_phi_n", &OpenScofo::MarkovState::IOIHatPhiN)
        .def_readwrite("audio_states", &OpenScofo::MarkovState::AudioStates)
        .def_readwrite("duration", &OpenScofo::MarkovState::Duration)
        .def_readwrite("line", &OpenScofo::MarkovState::Line);

    py::class_<OpenScofo::AudioState>(m, "AudioState")
        .def(py::init<>())
        .def_readwrite("frequency", &OpenScofo::AudioState::Freq)
        .def_readwrite("index", &OpenScofo::AudioState::Index);

    py::enum_<OpenScofo::EventType>(m, "EventType")
        .value("REST", OpenScofo::REST)
        .value("NOTE", OpenScofo::NOTE)
        .value("CHORD", OpenScofo::CHORD)
        .value("TRILL", OpenScofo::TRILL)
        .value("MULTI", OpenScofo::MULTI)
        .export_values(); // optional
}
