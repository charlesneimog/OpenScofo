/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

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

#include <pffft/pffft.h>

#include "log.hpp"
#include "onnx.hpp"
#include "states.hpp"
namespace OpenScofo {

#include <onsetsds.h>

namespace fs = std::filesystem;

// ╭─────────────────────────────────────╮
// │     Music Information Retrieval     │
// ╰─────────────────────────────────────╯
class MIR {
  public:
    MIR() = default;
    ~MIR();
    MIR(const MIR &) = delete;
    MIR &operator=(const MIR &) = delete;
    void UpdateConfiguration(const Configuration &Config);

    void SetdBTreshold(double dB);
    void GetDescription(const std::vector<double> &In, Description &Desc);
    void AddReverb(Description &Desc, double decay);

    double GetdB();

    // AI
    void ONNXInit(fs::path path, std::vector<Descriptors> Descriptors);
    std::vector<std::string> GetONNXLabels();

  private:
    Configuration m_Config;

    double HzToOcts(double frequency, double tuning, int binsPerOctave) const;
    double PositiveRemainder(double value, double modulus) const;
    bool DescriptorRequested(Descriptors Descriptor) const;
    void UpdateDescriptorFlags();
    void ComputeScalarFeatures(Description &Desc, const SpectralAccumulators &acc, size_t NHalf);
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

    // FFT
    void FFTInit();
    // Get Signal
    void InitITURFilters(void);
    void GetSignalPower(const std::vector<double> &In, Description &Desc);
    void GetSpectralFlux(Description &Desc);
    void YINInit();
    void YINExec(const std::vector<double> &In, Description &Desc);

  private:
    // FFT
    float *m_FullFFTIn = nullptr;
    float *m_FullFFTOut = nullptr;
    float *m_FullFFTWork = nullptr;
    PFFFT_Setup *m_FullFFTSetup = nullptr;
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
    std::vector<float> m_OnsetFFTFrame;

    // MFCC
    std::vector<std::vector<double>> m_MFCCFilter;
    std::vector<std::vector<double>> m_DCTBasis;
    std::vector<double> m_MFCCEnergy;
    std::vector<std::pair<int, int>> m_MFCCActiveBins;

    // Chroma
    std::vector<std::vector<double>> m_ChromaFilter;
    std::vector<double> m_ZCRScratch;

    // Machine Learning
    ONNXModel m_ONNXModel;

    // Env
    const std::array<double, 3> m_48kB1 = {1.53512485958697, -2.69169618940638, 1.19839281085285};
    const std::array<double, 3> m_48kA1 = {1.0, -1.69065929318241, 0.73248077421585};
    const std::array<double, 3> m_48kB2 = {1.0, -2.0, 1.0};
    const std::array<double, 3> m_48kA2 = {1.0, -1.99004745483398, 0.99007225036621};

    std::array<double, 3> m_B1;
    std::array<double, 3> m_A1;
    std::array<double, 3> m_B2;
    std::array<double, 3> m_A2;

    // Audio
    std::vector<double> m_PreviousSpectralPower;
    std::vector<double> m_SpectralPrefix;
    std::vector<double> m_YINDifference;
    std::vector<double> m_YINCMNDF;
    double m_PrevCentroid = 0.0;

    bool m_NeedYIN = false;
    bool m_NeedMFCC = false;
    bool m_NeedChroma = false;
    bool m_NeedZCR = false;
    bool m_NeedOnset = false;
    bool m_NeedExtendedTech = false;
    bool m_NeedONNX = false;

    // Time
    double m_EventTimeElapsed = 0.0; // ms
};
} // namespace OpenScofo
