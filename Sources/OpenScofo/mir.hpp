#pragma once

#include <math.h>
#include <vector>
#include <array>
#include <filesystem>
#include <unordered_map>
#include <complex>
#include <cstdint>

#include <fstream>

#include <fftw3.h>
#include <onnx.h>
#include <onsetsds.h>

#define CURRENT_ONNX_OPSET 24
#include "log.hpp"
#include "states.hpp"
namespace OpenScofo {

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
    void UpdateAudioParameters(float Sr, float WindowSize, float HopSize);

    void SetdBTreshold(double dB);
    void GetDescription(std::vector<double> &In, Description &Desc);
    void AddReverb(Description &Desc, double decay);

    double GetdB();
    void LoadONNXModel(fs::path path);

  private:
    double Mtof(double Note, double Tunning);
    double Ftom(double Freq, double Tunning);
    double Freq2Bin(double freq, double n, double Sr);
    void GetSpectralDescriptions(Description &Desc);
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
    // Zero Crossing Rate
    void ZeroCrossingRateInit();
    void ZeroCrossingRateExec(std::vector<double> &In, Description &Desc);

    // FFTW
    void FFTWInit();

    // CQT
    void CQTInit();

    // Get Signal
    void GetSignalPower(std::vector<double> &In, Description &Desc);
    void GetSpectralFlux(Description &Desc);
    void YINInit();
    void YINExec(std::vector<double> &In, Description &Desc);

  private:
    // FFT
    double *m_FFTIn = nullptr;
    fftw_complex *m_FFTOut = nullptr;
    fftw_plan m_FFTPlan = nullptr;
    std::vector<std::pair<int, int>> m_FakeCQT;
    std::vector<double> m_WindowingFunc;
    double m_PrevPercussiveProb;
    double m_PrevRMS;
    double m_PeakFlux;
    double m_PeakDeltaRMS;
    double m_PeakFlatness;

    // Onsets
    bool m_OnsetInit = false;
    OnsetsDS *m_ODS = nullptr;
    float *m_ODSData = nullptr;
    int m_OnsetFFTSize = 512;
    int m_MedSpan = 50;
    int m_Accum = 0;

    // MFCC
    int m_MFCCMels = 40;
    int m_MFCC = 13;
    std::vector<std::vector<double>> m_MFCCFilter;
    std::vector<std::vector<double>> m_DCTBasis;
    std::vector<double> m_MFCCEnergy;
    std::vector<std::pair<int, int>> m_MFCCActiveBins;

    // Chroma
    std::vector<int> m_ChromaBinMap;
    size_t m_ChromaSize = 12;

    // Zero-crossing rate (librosa-like defaults)
    int m_ZCRFrameLength = 2048;
    int m_ZCRHopLength = 512;
    bool m_ZCRCenter = true;
    bool m_ZCRPad = false;
    bool m_ZCRZeroPos = true;
    double m_ZCRThreshold = 1e-10;
    std::vector<double> m_ZCRScratch;

    // Machine Learning
    bool m_ONNXModelLoaded = false;
    struct onnx_context_t *m_OnnxCTX = nullptr;
    int m_TensorCount = 0;
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
    std::vector<double> m_SpectralPrefix;
    std::vector<double> m_YINDifference;
    std::vector<double> m_YINCMNDF;
    double m_PrevCentroid = 0.0;
    double m_YINThreshold = 0.15;
    double m_YINMinFrequency = 50.0;
    double m_YINMaxFrequency = 2000.0;

    // Time
    double m_EventTimeElapsed = 0.0; // ms
};
} // namespace OpenScofo
