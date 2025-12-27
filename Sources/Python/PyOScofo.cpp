#include <OScofo.hpp>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_OScofo, m) {

    py::class_<OScofo::OScofo>(m, "OScofo")
        .def(py::init([](float sr, float fft_size, float hop) {
            auto obj = new OScofo::OScofo(sr, fft_size, hop);
            obj->SetErrorCallBack([](const std::string &msg) {
                py::gil_scoped_acquire acquire;
                throw py::value_error(msg);
            });
            return obj;
        }))

        // Score
        .def("parse_score", &OScofo::OScofo::ParseScore)

        // Config
        .def("set_pitch_template_sigma", &OScofo::OScofo::SetPitchTemplateSigma)
        .def("set_harmonics", &OScofo::OScofo::SetHarmonics)
        .def("set_db_threshold", &OScofo::OScofo::SetdBTreshold)
        .def("set_tuning", &OScofo::OScofo::SetTunning)
        .def("set_current_event", &OScofo::OScofo::SetCurrentEvent)

        // Get Info
        .def("get_live_bpm", &OScofo::OScofo::GetLiveBPM)
        .def("get_event_index", &OScofo::OScofo::GetEventIndex)
        .def("get_error", &OScofo::OScofo::GetErrorMessage)
        .def("get_states", &OScofo::OScofo::GetStates)
        .def("get_pitch_template", &OScofo::OScofo::GetPitchTemplate)

        // Help & Test Functions
        .def("get_audio_description", &OScofo::OScofo::GetAudioDescription)

        // Process
        .def("process_block", [](OScofo::OScofo &self, py::array_t<double> audio) {
            py::buffer_info bufInfo = audio.request();
            if (bufInfo.ndim != 1) {
                throw std::runtime_error("Input array must be 1-dimensional");
            }
            std::vector<double> cpp_audio(static_cast<double *>(bufInfo.ptr), static_cast<double *>(bufInfo.ptr) + bufInfo.shape[0]);
            return self.ProcessBlock(cpp_audio);
        });

    // Description Class
    py::class_<OScofo::Description>(m, "Description")
        .def(py::init<>())
        .def_readwrite("mfcc", &OScofo::Description::MFCC)
        .def_readwrite("onset", &OScofo::Description::Onset)
        .def_readwrite("silence_prob", &OScofo::Description::SilenceProb)
        .def_readwrite("spectral_power", &OScofo::Description::SpectralPower)
        .def_readwrite("norm_spectral_power", &OScofo::Description::NormSpectralPower)

        .def_readwrite("loudness", &OScofo::Description::Loudness)
        .def_readwrite("spectral_flux", &OScofo::Description::SpectralFlux)

        .def_readwrite("db", &OScofo::Description::dB)
        .def_readwrite("rms", &OScofo::Description::RMS)

        .def_readwrite("power", &OScofo::Description::Power);

    // State Class
    py::class_<OScofo::MacroState>(m, "State")
        .def(py::init<>())
        .def_readwrite("index", &OScofo::MacroState::Index)
        .def_readwrite("position", &OScofo::MacroState::ScorePos)
        .def_readwrite("type", &OScofo::MacroState::Type)
        .def_readwrite("markov", &OScofo::MacroState::Markov)
        .def_readwrite("freqs", &OScofo::MacroState::Freqs)
        .def_readwrite("kl_div", &OScofo::MacroState::Obs)
        .def_readwrite("forward", &OScofo::MacroState::Forward)
        .def_readwrite("bpm_expected", &OScofo::MacroState::BPMExpected)
        .def_readwrite("bpm_observed", &OScofo::MacroState::BPMObserved)
        .def_readwrite("onset_expected", &OScofo::MacroState::OnsetExpected)
        .def_readwrite("onset_observed", &OScofo::MacroState::OnsetObserved)
        .def_readwrite("phase_expected", &OScofo::MacroState::PhaseExpected)
        .def_readwrite("phase_observed", &OScofo::MacroState::PhaseObserved)
        .def_readwrite("ioi_phi_n", &OScofo::MacroState::IOIPhiN)
        .def_readwrite("ioi_hat_phi_n", &OScofo::MacroState::IOIHatPhiN)
        .def_readwrite("duration", &OScofo::MacroState::Duration)
        .def_readwrite("line", &OScofo::MacroState::Line)
        .def("__repr__", &OScofo::MacroState::__repr__)
        .def("__str__", &OScofo::MacroState::__repr__);
}
