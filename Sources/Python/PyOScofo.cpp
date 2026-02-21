#include <OpenScofo.hpp>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_OpenScofo, m) {

    py::class_<OpenScofo::OpenScofo>(m, "OpenScofo")
        .def(py::init([](float sr, float fft_size, float hop) {
            auto obj = new OpenScofo::OpenScofo(sr, fft_size, hop);
            obj->SetErrorCallBack([](const std::string &msg) {
                py::gil_scoped_acquire acquire;
                throw py::value_error(msg);
            });
            return obj;
        }))

        // Python
        .def("__repr__",
             [](OpenScofo::OpenScofo &self) {
                 return std::format("<OpenScofo(sr={}, fft_size={}, hop={})>", self.GetSr(), self.GetFFTSize(), self.GetHopSize());
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
        .def("get_error", &OpenScofo::OpenScofo::GetErrorMessage)
        .def("get_states", &OpenScofo::OpenScofo::GetStates)
        .def("get_pitch_template", &OpenScofo::OpenScofo::GetPitchTemplate)
        .def("get_cqt_template", &OpenScofo::OpenScofo::GetCQTTemplate)

        // Help & Test Functions
        .def("get_audio_description", &OpenScofo::OpenScofo::GetAudioDescription)

        // Time Template
           .def("get_time_coherence_template",
               &OpenScofo::OpenScofo::GetTimeCoherenceTemplate,
               py::arg("pos"),
               py::arg("time_in_event") = 0)
           .def("get_time_coherence_confiability",
               &OpenScofo::OpenScofo::GetTimeCoherenceConfiability,
               py::arg("event_values"))

        // Process
        .def("process_block", [](OpenScofo::OpenScofo &self, py::array_t<double> audio) {
            py::buffer_info bufInfo = audio.request();
            if (bufInfo.ndim != 1) {
                throw std::runtime_error("Input array must be 1-dimensional");
            }
            std::vector<double> cpp_audio(static_cast<double *>(bufInfo.ptr), static_cast<double *>(bufInfo.ptr) + bufInfo.shape[0]);
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
    py::class_<OpenScofo::MacroState>(m, "State")
        .def(py::init<>())
        .def_readwrite("position", &OpenScofo::MacroState::ScorePos)
        .def_readwrite("type", &OpenScofo::MacroState::Type)
        .def_readwrite("markov", &OpenScofo::MacroState::Markov)
        .def_readwrite("forward", &OpenScofo::MacroState::Forward)
        .def_readwrite("bpm_expected", &OpenScofo::MacroState::BPMExpected)
        .def_readwrite("bpm_observed", &OpenScofo::MacroState::BPMObserved)
        .def_readwrite("onset_expected", &OpenScofo::MacroState::OnsetExpected)
        .def_readwrite("onset_observed", &OpenScofo::MacroState::OnsetObserved)
        .def_readwrite("phase_expected", &OpenScofo::MacroState::PhaseExpected)
        .def_readwrite("phase_observed", &OpenScofo::MacroState::PhaseObserved)
        .def_readwrite("ioi_phi_n", &OpenScofo::MacroState::IOIPhiN)
        .def_readwrite("ioi_hat_phi_n", &OpenScofo::MacroState::IOIHatPhiN)
        .def_readwrite("audiostates", &OpenScofo::MacroState::AudioStates)
        .def_readwrite("duration", &OpenScofo::MacroState::Duration)
        .def_readwrite("line", &OpenScofo::MacroState::Line)
        .def("__repr__", &OpenScofo::MacroState::__repr__)
        .def("__str__", &OpenScofo::MacroState::__repr__);

    py::class_<OpenScofo::AudioState>(m, "AudioState")
        .def(py::init<>())
        .def_readwrite("freq", &OpenScofo::AudioState::Freq)
        .def_readwrite("index", &OpenScofo::AudioState::Index);
}
