#pragma once

#include <math.h>
#include <vector>
#include <array>
#include <filesystem>
#include <unordered_map>
#include <complex>
#include <cstdint>
#include <numbers>

#include <fstream>

#include <fftw3.h>
#include <onnx.h>

#define CURRENT_ONNX_OPSET 24
#include "log.hpp"
#include "states.hpp"
namespace OpenScofo {

#include <onsetsds.h>

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
    void GetDescription(const std::vector<double> &In, Description &Desc);
    void AddReverb(Description &Desc, double decay);

    double GetdB();

    // AI
    void ONNXInit(fs::path path, std::vector<Descriptors> Descriptors);
    std::vector<std::string> GetONNXLabels();

  private:
    double HzToOcts(double frequency, double tuning, int binsPerOctave) const;
    double PositiveRemainder(double value, double modulus) const;
    void GetSpectralDescriptions(Description &Desc);

    // MFCC
    void MFCCInit();
    void MFCCExec(Description &Desc);

    // Onset
    void OnsetInit();
    void OnsetExec(Description &Desc);
    // Extended Technique
    void ExtendedTechExec(Description &Desc);

    // Chroma
    void SpectralChromaInit();
    void SpectralChromaExec(Description &Desc);

    // Zero Crossing Rate
    void ZeroCrossingRateInit();
    void ZeroCrossingRateExec(const std::vector<double> &In, Description &Desc);

    // FFTW
    void FFTWInit();

    void ONNXExec(Description &Desc);

    // Get Signal
    void InitITURFilters(void);
    void GetSignalPower(const std::vector<double> &In, Description &Desc);
    void GetSpectralFlux(Description &Desc);
    void YINInit();
    void YINExec(const std::vector<double> &In, Description &Desc);

  private:
    // FFT
    double *m_FullFFTIn = nullptr;
    fftw_complex *m_FullFFTOut = nullptr;
    fftw_plan m_FullFFTPlan = nullptr;
    std::vector<double> m_FullWindowingFunc;

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
    std::vector<float> m_OnsetFFTFrame;

    // MFCC
    int m_MFCCMels = 40;
    int m_MFCCCount = 13;
    std::vector<std::vector<double>> m_MFCCFilter;
    std::vector<std::vector<double>> m_DCTBasis;
    std::vector<double> m_MFCCEnergy;
    std::vector<std::pair<int, int>> m_MFCCActiveBins;

    // Chroma
    Matrix m_ChromaFilter;
    int m_ChromaSize = 12;
    double m_ChromaA440 = 440.0;
    double m_ChromaTuning = 0.0;
    double m_ChromaCenterOctave = 5.0;
    double m_ChromaOctaveWidth = 2.0;

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
    struct onnx_context_t *m_ONNXContext = nullptr;
    std::vector<std::string> m_ONNXLabels;
    std::vector<Descriptors> m_ONNXDescriptors;
    std::vector<float> m_ONNXDescriptorsArray;
    std::vector<std::function<void(const Description &, float *&)>> m_Writers;
    struct onnx_tensor_t *m_InputTensor;
    struct onnx_tensor_t *m_OutputTensor;
    int m_ONNXDescriptorsSize = 0;

    // Env
    double m_dBTreshold = -50;
    const std::array<double, 3> m_48kB1 = {1.53512485958697, -2.69169618940638, 1.19839281085285};
    const std::array<double, 3> m_48kA1 = {1.0, -1.69065929318241, 0.73248077421585};
    const std::array<double, 3> m_48kB2 = {1.0, -2.0, 1.0};
    const std::array<double, 3> m_48kA2 = {1.0, -1.99004745483398, 0.99007225036621};

    std::array<double, 3> m_B1;
    std::array<double, 3> m_A1;
    std::array<double, 3> m_B2;
    std::array<double, 3> m_A2;

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
    double m_SpectralRolloffCutoff = 0.85;
    double m_YINThreshold = 0.15;
    double m_YINMinFrequency = 50.0;
    double m_YINMaxFrequency = 2000.0;

    // Time
    double m_EventTimeElapsed = 0.0; // ms
};
} // namespace OpenScofo
