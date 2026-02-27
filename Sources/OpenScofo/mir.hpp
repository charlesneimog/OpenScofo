#pragma once

#include <math.h>
#include <vector>
#include <array>
#include <filesystem>
#include <unordered_map>
#include <cstdint>

#include <fstream>

#include <fftw3.h>
#include <onnx.h>
#include <onsetsds.h>

#define CURRENT_ONNX_OPSET 24
#include "log.hpp"
#include "states.hpp"
namespace OpenScofo {

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
    void GetDescription(std::vector<double> &In, Description &Desc, States &States);
    void AddReverb(Description &Desc, double decay);

    double GetdB();
    void LoadONNXModel(fs::path path);
    void BuildTimeCoherenceTemplate(States &States);
    bool HasErrors();
    std::vector<std::string> GetErrorMessage();
    void SetError(const std::string &message);
    void ClearError();

    // Tests
    std::vector<std::pair<int, int>> GetCQT();

  private:
    double Mtof(double Note, double Tunning);
    double Ftom(double Freq, double Tunning);
    double Freq2Bin(double freq, double n, double Sr);
    void GetFFTDescriptions(std::vector<double> &In, Description &Desc);
    // MFCC
    void MFCCInit();
    void MFCCExec(Description &Desc);
    // Time coherence
    void BuildSingleEventPdf(MarkovState &ev, double dt);

    const std::vector<float> &GetTimeCoherenceGaussianKernel(double sigmaSeconds, double dt, int templateMax) const;
    // Onset
    void OnsetInit();
    void OnsetExec(Description &Desc);
    // Spectral Flux
    void SpectralFluxInit();
    void SpectralFluxExec(Description &Desc);
    // Spectral Flatness
    void SpectralFlatnessInit();
    void SpectralFlatnessExec(Description &Desc);
    // Harmonicity
    void SpectralHarmonicityInit();
    void SpectralHarmonicityExec(Description &Desc);
    // Chroma
    void SpectralChromaInit();
    void SpectralChromaExec(Description &Desc);

    // CQT
    void CQTInit();

    // Get Signal
    void GetSignalPower(std::vector<double> &In, Description &Desc);
    void GetSpectralFlux(Description &Desc);

  private:
    // FFT
    double *m_FFTIn;
    fftw_complex *m_FFTOut;
    fftw_plan m_FFTPlan;
    std::vector<std::pair<int, int>> m_FakeCQT;
    std::vector<double> m_WindowingFunc;

    // Onsets
    bool m_OnsetInit = false;
    OnsetsDS *m_ODS = nullptr;
    float *m_ODSData = nullptr;
    int m_OnsetFFTSize = 512;
    int m_MedSpan = 20;
    int m_Accum = 0;

    // MFCC
    int m_MFCCMels = 40;
    int m_MFCC = 13;
    std::vector<std::vector<float>> m_MFCCFilter;
    std::vector<std::vector<float>> m_DCTBasis;
    std::vector<float> m_MFCCEnergy;

    // Chroma
    std::vector<int> m_ChromaBinMap;
    size_t m_ChromaSize = 12;

    // Machine Learning
    bool m_ONNXModelLoaded = false;
    struct onnx_context_t *m_OnnxCTX;
    int m_TensorCount;
    std::unordered_map<std::string, float> m_ONNXResults;

    // Env
    double m_dBTreshold = -50;

    // Audio
    float m_FFTSize;
    float m_BlockSize;
    float m_HopSize;
    float m_Sr;
    double m_dB;
    std::vector<double> m_PreviousSpectralPower;

    // Time
    double m_EventTimeElapsed = 0.0; // ms
    unsigned m_TimeCoherenceTemplateSize = 1024;

    // Errors
    bool m_HasErrors = false;
    std::vector<std::string> m_Errors;

    // Time coherence peak kernels (keyed by quantized sigma+dt)
    mutable std::unordered_map<std::uint64_t, std::vector<float>> m_TimeCoherenceKernelCache;
};
} // namespace OpenScofo
