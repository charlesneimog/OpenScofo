#pragma once

#include <math.h>
#include <vector>
#include <array>
#include <filesystem>
#include <unordered_map>

#include <fstream>

#include <fftw3.h>
#include <onnx.h>
#include <onsetsds.h>

#define CURRENT_ONNX_OPSET 24
#include "log.hpp"
#include "states.hpp"
namespace OScofo {

#ifndef TWO_PI
#define TWO_PI (2 * M_PI)
#endif

namespace fs = std::filesystem;

using Spectrum = std::vector<std::complex<double>>;
using Matrix = std::vector<std::vector<double>>; // row-major: rows x cols

// ╭─────────────────────────────────────╮
// │     Music Information Retrieval     │
// ╰─────────────────────────────────────╯
class MIR {
  public:
    MIR(float Sr, float WindowSize, float HopSize);
    ~MIR();

    void SetdBTreshold(double dB);
    void GetDescription(std::vector<double> &In, Description &Desc);
    double GetdB();

    // AI
    void LoadONNXModel(fs::path path);

    // Time coherence
    void BuildTimeCoherenceTemplate(States &States);

    // Error handling
    bool HasErrors();
    std::vector<std::string> GetErrorMessage();
    void SetError(const std::string &message);
    void ClearError();

  private:
    // Helpers
    std::vector<double> m_WindowingFunc;
    double Mtof(double Note, double Tunning);
    double Ftom(double Freq, double Tunning);
    double Freq2Bin(double freq, double n, double Sr);

    // FFT
    double *m_FFTIn;
    fftw_complex *m_FFTOut;
    fftw_plan m_FFTPlan;
    void GetFFTDescriptions(std::vector<double> &In, Description &Desc);

    // MFCC
    void MFCCInit();
    void MFCCExec(Description &Desc);
    int m_MFCCMels = 40;
    int m_MFCC = 13;
    std::vector<std::vector<double>> m_MFCCFilter;
    std::vector<std::vector<double>> m_DCTBasis;
    std::vector<double> m_MFCCEnergy;

    // Time coherence
    void BuildSingleEventPdf(MacroState &ev, double dt);

    // Onsets
    bool m_OnsetInit = false;
    void OnsetInit();
    void OnsetExec(Description &Desc);
    OnsetsDS *m_ODS = nullptr;
    float *m_ODSData = nullptr;
    int m_OnsetFFTSize = 512;
    int m_MedSpan = 20;
    int m_Accum = 0;

    // Machine Learning
    bool m_ONNXModelLoaded = false;
    struct onnx_context_t *m_OnnxCTX;
    int m_TensorCount;
    std::unordered_map<std::string, float> m_ONNXResults;

    // Env
    double m_dBTreshold = -50;
    void GetSignalPower(std::vector<double> &In, Description &Desc);
    void GetSpectralFlux(Description &Desc);

    // Audio
    float m_FFTSize;
    float m_BlockSize;
    float m_HopSize;
    float m_Sr;
    double m_dB;
    std::vector<double> m_PreviousSpectralPower;

    // Time
    double m_EventTimeElapsed = 0.0; // ms

    // Errors
    bool m_HasErrors = false;
    std::vector<std::string> m_Errors;
};
} // namespace OScofo
