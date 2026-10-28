/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file features.cpp
 * @brief Temporal, mel, chroma, onset, power, and loudness descriptors.
 */

#include "mir.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │           Onset Detector            │
// ╰─────────────────────────────────────╯
/**
 * @brief Allocate and initialize the onset detector.
 *
 * @note Replaces existing detector storage and prepares an interleaved complex FFT frame for OnsetsDS.
 */
void MIR::OnsetInit() {
    m_OnsetInit = false;

    const size_t nbytes = onsetsds_memneeded(m_Config.OnsetType, m_Config.FFTSize, m_Config.MedSpan);
    delete[] m_ODSData;
    m_ODSData = nullptr;
    delete m_ODS;
    m_ODS = nullptr;

    m_ODSData = new float[nbytes / sizeof(float)];
    m_ODS = new OnsetsDS();
    if (!m_ODS || !m_ODSData) {
        spdlog::critical("Not possible to initialize the onset detector");
        return;
    }

    onsetsds_init(m_ODS, m_ODSData, ODS_FFT_FFTW3_R2C, m_Config.OnsetType, m_Config.FFTSize, m_Config.MedSpan,
                  m_Config.SR);
    m_OnsetFFTFrame.assign(static_cast<size_t>(2 * (m_Config.FFTSize / 2 + 1)), 0.0f);
    m_OnsetInit = true;
}

// ─────────────────────────────────────
/**
 * @brief Update onset evidence from the current FFT output.
 *
 * @param Desc Audio description to populate or update in place.
 *
 * @note Converts PFFFT output to the OnsetsDS format and stores its postprocessed detection-function value.
 * @warning Requires a current FFT frame; returns without updating Desc if the detector is not initialized.
 */
void MIR::OnsetExec(Description &Desc) {
    if (!m_OnsetInit)
        return;

    const size_t nBins = static_cast<size_t>(m_Config.FFTSize / 2 + 1);
    for (size_t i = 0; i < nBins; ++i) {
        double re = 0.0;
        double im = 0.0;
        const size_t half = static_cast<size_t>(m_Config.FFTSize) / 2;
        if (i == 0) {
            re = static_cast<double>(m_FullFFTOut[0]);
            im = 0.0;
        } else if (i == half) {
            re = static_cast<double>(m_FullFFTOut[1]);
            im = 0.0;
        } else {
            const size_t idx = 2 * i;
            re = static_cast<double>(m_FullFFTOut[idx]);
            im = static_cast<double>(m_FullFFTOut[idx + 1]);
        }
        m_OnsetFFTFrame[2 * i] = static_cast<float>(re);
        m_OnsetFFTFrame[2 * i + 1] = static_cast<float>(im);
    }

    (void)onsetsds_process(m_ODS, m_OnsetFFTFrame.data());
    Desc.Onset = m_ODS->odfvalpost;
}

// ╭─────────────────────────────────────╮
// │         Power and Amplitude         │
// ╰─────────────────────────────────────╯
// Check https://github.com/klangfreund/LUFSMeter (use MIT)
/**
 * @brief Adapt the two loudness-weighting filter stages to the sample rate.
 *
 * @note Derives shelving and high-pass coefficients from the stored 48 kHz reference coefficients.
 * @warning Requires a positive configured sample rate.
 */
void MIR::InitITURFilters() {
    // Stage 1: shelving filter
    double KoverQ1 = (2.0 - 2.0 * m_48kA1[2]) / (m_48kA1[2] - m_48kA1[1] + 1.0);
    double K1 = std::sqrt((m_48kA1[1] + m_48kA1[2] + 1.0) / (m_48kA1[2] - m_48kA1[1] + 1.0));
    double Q1 = K1 / KoverQ1;
    double arctanK1 = std::atan(K1);
    double VB1 = (m_48kB1[0] - m_48kB1[2]) / (1.0 - m_48kA1[2]);
    double VH1 = (m_48kB1[0] - m_48kB1[1] + m_48kB1[2]) / (m_48kA1[2] - m_48kA1[1] + 1.0);
    double VL1 = (m_48kB1[0] + m_48kB1[1] + m_48kB1[2]) / (m_48kA1[1] + m_48kA1[2] + 1.0);

    double Knew1 = std::tan(arctanK1 * 48000.0 / m_Config.SR);
    double commonFactor1 = 1.0 / (1.0 + Knew1 / Q1 + Knew1 * Knew1);

    m_B1[0] = (VH1 + VB1 * Knew1 / Q1 + VL1 * Knew1 * Knew1) * commonFactor1;
    m_B1[1] = 2.0 * (VL1 * Knew1 * Knew1 - VH1) * commonFactor1;
    m_B1[2] = (VH1 - VB1 * Knew1 / Q1 + VL1 * Knew1 * Knew1) * commonFactor1;
    m_A1[0] = 1.0;
    m_A1[1] = 2.0 * (Knew1 * Knew1 - 1.0) * commonFactor1;
    m_A1[2] = (1.0 - Knew1 / Q1 + Knew1 * Knew1) * commonFactor1;

    // Stage 2: high-pass filter
    double KoverQ2 = (2.0 - 2.0 * m_48kA2[2]) / (m_48kA2[2] - m_48kA2[1] + 1.0);
    double K2 = std::sqrt((m_48kA2[1] + m_48kA2[2] + 1.0) / (m_48kA2[2] - m_48kA2[1] + 1.0));
    double Q2 = K2 / KoverQ2;
    double arctanK2 = std::atan(K2);
    double VB2 = (m_48kB2[0] - m_48kB2[2]) / (1.0 - m_48kA2[2]);
    double VH2 = (m_48kB2[0] - m_48kB2[1] + m_48kB2[2]) / (m_48kA2[2] - m_48kA2[1] + 1.0);
    double VL2 = (m_48kB2[0] + m_48kB2[1] + m_48kB2[2]) / (m_48kA2[1] + m_48kA2[2] + 1.0);

    double Knew2 = std::tan(arctanK2 * 48000.0 / m_Config.SR);
    double commonFactor2 = 1.0 / (1.0 + Knew2 / Q2 + Knew2 * Knew2);

    m_B2[0] = (VH2 + VB2 * Knew2 / Q2 + VL2 * Knew2 * Knew2) * commonFactor2;
    m_B2[1] = 2.0 * (VL2 * Knew2 * Knew2 - VH2) * commonFactor2;
    m_B2[2] = (VH2 - VB2 * Knew2 / Q2 + VL2 * Knew2 * Knew2) * commonFactor2;
    m_A2[0] = 1.0;
    m_A2[1] = 2.0 * (Knew2 * Knew2 - 1.0) * commonFactor2;
    m_A2[2] = (1.0 - Knew2 / Q2 + Knew2 * Knew2) * commonFactor2;
}

// ─────────────────────────────────────
/**
 * @brief Compute RMS, decibel level, weighted loudness, and silence evidence.
 *
 * @param In Input audio frame containing the configured FFT window samples.
 * @param Desc Audio description to populate or update in place.
 *
 * @note Resets filter delay values for each frame and floors silent decibel and loudness values to -100.
 * @warning In must be nonempty and loudness filter coefficients must be initialized.
 */
void MIR::GetSignalPower(const std::vector<double> &In, Description &Desc) {
    double x1_1 = 0.0, x2_1 = 0.0;
    double y1_1 = 0.0, y2_1 = 0.0;
    double x1_2 = 0.0, x2_2 = 0.0;
    double y1_2 = 0.0, y2_2 = 0.0;

    double z = 0.0;
    double z_loudness = 0.0;
    for (double sample : In) {
        double s1 = m_B1[0] * sample + m_B1[1] * x1_1 + m_B1[2] * x2_1 - m_A1[1] * y1_1 - m_A1[2] * y2_1;
        x2_1 = x1_1;
        x1_1 = sample;
        y2_1 = y1_1;
        y1_1 = s1;
        double s2 = m_B2[0] * s1 + m_B2[1] * x1_2 + m_B2[2] * x2_2 - m_A2[1] * y1_2 - m_A2[2] * y2_2;
        x2_2 = x1_2;
        x1_2 = s1;
        y2_2 = y1_2;
        y1_2 = s2;

        z_loudness += s2 * s2;
        z += sample * sample;
    }

    // Compute RMS
    double rms = std::sqrt(z / In.size());
    Desc.RMS = rms;

    // Convert RMS to dB
    Desc.dB = 20.0 * std::log10(rms);
    if (std::isinf(Desc.dB)) {
        Desc.dB = -100; // handle silence
    }

    // Loudness (based on sum of squares)
    double meanSquare = z_loudness / In.size();
    if (meanSquare <= 0.0) {
        Desc.Loudness = -100.0; // substitui -inf
    } else {
        Desc.Loudness = -0.691 + 10.0 * std::log10(meanSquare);
    }

    // Compute silence probability
    // TODO: Add the capability to set this via Score
    const double L0 = -60.0;
    const double alpha = 0.25;
    Desc.SilenceProb = 1.0 / (1.0 + std::exp(alpha * (Desc.Loudness - L0)));
}

// ╭─────────────────────────────────────╮
// │                MFCC                 │
// ╰─────────────────────────────────────╯
/**
 * @brief Precompute the mel filterbank and orthonormal DCT-II basis.
 *
 * @note Uses Slaney mel spacing and area normalization, recording active FFT bin ranges for each filter.
 * @warning Requires positive sample rate, FFT size, mel count, and MFCC count.
 */
void MIR::MFCCInit() {
    const int FFTSize = m_Config.FFTSize;
    const int NumBins = FFTSize / 2 + 1;
    const int NumMels = m_Config.MFCCMels;
    const int NumMFCC = m_Config.MFCCCount;
    const double SR = static_cast<double>(m_Config.SR);

    // librosa.filters.mel(htk=False, norm="slaney")
    constexpr double FMin = 0.0;
    const double FMax = SR * 0.5;

    constexpr double FSp = 200.0 / 3.0;
    constexpr double MinLogHz = 1000.0;
    constexpr double MinLogMel = MinLogHz / FSp;
    const double LogStep = std::log(6.4) / 27.0;

    auto HzToMel = [&](double hz) noexcept {
        return (hz < MinLogHz) ? hz / FSp : MinLogMel + std::log(hz / MinLogHz) / LogStep;
    };

    auto MelToHz = [&](double mel) noexcept {
        return (mel < MinLogMel) ? mel * FSp : MinLogHz * std::exp((mel - MinLogMel) * LogStep);
    };

    const double MelMin = HzToMel(FMin);
    const double MelMax = HzToMel(FMax);
    std::vector<double> HzPts(NumMels + 2);
    {
        const double MelStep = (MelMax - MelMin) / (NumMels + 1);
        for (int i = 0; i < NumMels + 2; ++i) {
            HzPts[i] = MelToHz(MelMin + MelStep * i);
        }
    }

    // FFT bin frequencies (exact np.fft.rfftfreq behavior)
    std::vector<double> FFTFreqs(NumBins);
    {
        const double BinHz = SR / static_cast<double>(FFTSize);
        for (int k = 0; k < NumBins; ++k) {
            FFTFreqs[k] = static_cast<double>(k) * BinHz;
        }
    }

    // Mel filterbank
    // Exact librosa.filters.mel(..., norm="slaney")
    m_MFCCFilter.assign(NumMels, std::vector<double>(NumBins, 0.0));
    m_MFCCActiveBins.assign(NumMels, {0, 0});

    for (int m = 0; m < NumMels; ++m) {
        const double Left = HzPts[m];
        const double Center = HzPts[m + 1];
        const double Right = HzPts[m + 2];

        const double InvLeftWidth = 1.0 / (Center - Left);
        const double InvRightWidth = 1.0 / (Right - Center);

        // Slaney area normalization
        const double ENorm = 2.0 / (Right - Left);
        int First = NumBins;
        int Last = -1;
        double *Filter = m_MFCCFilter[m].data();
        for (int k = 0; k < NumBins; ++k) {
            const double f = FFTFreqs[k];
            double w;
            if (f >= Left && f <= Center) {
                w = (f - Left) * InvLeftWidth;
            } else if (f > Center && f <= Right) {
                w = (Right - f) * InvRightWidth;
            } else {
                continue;
            }
            const double v = w * ENorm;
            Filter[k] = v;
            First = std::min(First, k);
            Last = k;
        }
        if (Last < 0) {
            First = 0;
            Last = 0;
        }
        m_MFCCActiveBins[m] = {First, Last};
    }

    // DCT-II basis
    // Exact scipy.fftpack.dct(..., type=2, norm="ortho")
    m_DCTBasis.assign(NumMFCC, std::vector<double>(NumMels));
    const double Scale0 = std::sqrt(1.0 / NumMels);
    const double Scale = std::sqrt(2.0 / NumMels);
    const double Factor = std::numbers::pi / NumMels;
    for (int k = 0; k < NumMFCC; ++k) {
        const double Norm = (k == 0) ? Scale0 : Scale;
        double *Basis = m_DCTBasis[k].data();
        for (int n = 0; n < NumMels; ++n) {
            Basis[n] = Norm * std::cos(Factor * (n + 0.5) * k);
        }
    }
    m_MFCCEnergy.resize(NumMels);
}

// ─────────────────────────────────────────────────────────────────────────────

/**
 * @brief Compute log-mel energies and MFCC coefficients from the power spectrum.
 *
 * @param Desc Audio description to populate or update in place.
 *
 * @note Floors mel power at 1e-10, clips log-mel values to an 80 dB range, and applies the precomputed DCT-II.
 * @warning Initialize the filterbank and size Power, LogMelSpectrum, and MFCC arrays before calling.
 */
void MIR::MFCCExec(Description &Desc) {
    constexpr double kAmin = 1e-10;
    constexpr double kTopDb = 80.0;

    const int NumMels = m_Config.MFCCMels;
    const int NumMFCC = m_Config.MFCCCount;

    // Mel projection + power_to_db(ref=1.0, top_db=80)
    double MaxLog = -std::numeric_limits<double>::infinity();
    for (int m = 0; m < NumMels; ++m) {
        const auto &[First, Last] = m_MFCCActiveBins[m];
        const double *Filter = m_MFCCFilter[m].data();
        const double *Power = Desc.Power.data();
        double MelEnergy = 0.0;
        // Sparse accumulation using active range only
        for (int k = First; k <= Last; ++k) {
            MelEnergy += Filter[k] * Power[k];
        }
        m_MFCCEnergy[m] = MelEnergy;
        const double LogMel = 10.0 * std::log10(std::max(kAmin, MelEnergy));
        Desc.LogMelSpectrum[m] = LogMel;
        MaxLog = std::max(MaxLog, LogMel);
    }

    // librosa power_to_db(..., top_db=80)
    const double Floor = MaxLog - kTopDb;
    for (int m = 0; m < NumMels; ++m) {
        Desc.LogMelSpectrum[m] = std::max(Desc.LogMelSpectrum[m], Floor);
    }

    // MFCC = DCT-II(log-mel)
    const double *LogMel = Desc.LogMelSpectrum.data();
    for (int k = 0; k < NumMFCC; ++k) {
        const double *Basis = m_DCTBasis[k].data();
        double Sum = 0.0;
        for (int n = 0; n < NumMels; ++n) {
            Sum += Basis[n] * LogMel[n];
        }
        Desc.MFCC[k] = Sum;
    }
}

// ╭─────────────────────────────────────╮
// │               Chroma                │
// ╰─────────────────────────────────────╯
/**
 * @brief Convert frequency to octaves relative to the adjusted A4 reference.
 *
 * @param frequency Positive frequency in Hz.
 * @param tuning Tuning offset measured in chroma bins.
 * @param binsPerOctave Positive number of chroma bins per octave.
 *
 * @return Octave coordinate relative to adjusted A4 / 16.
 *
 * @note Uses configured A4 tuning and a chroma-bin tuning offset.
 * @warning Frequency, A4 tuning, and binsPerOctave must be positive.
 */
double MIR::HzToOcts(double frequency, double tuning, int binsPerOctave) const {
    const double a440 = m_Config.TuningA4 * std::pow(2.0, tuning / static_cast<double>(binsPerOctave));
    return std::log2(frequency / (a440 / 16.0));
}

// ─────────────────────────────────────
/**
 * @brief Wrap a value to a nonnegative remainder.
 *
 * @param value Value to wrap.
 * @param modulus Positive period of the remainder.
 *
 * @return Remainder in [0, modulus).
 *
 * @note Adds the modulus when std::fmod() returns a negative value.
 * @warning The modulus must be positive.
 */
double MIR::PositiveRemainder(double value, double modulus) const {
    double result = std::fmod(value, modulus);
    if (result < 0.0) {
        result += modulus;
    }
    return result;
}

// ─────────────────────────────────────
/**
 * @brief Build the tuned spectral chroma filterbank.
 *
 * @note Normalizes frequency columns, applies octave weighting, and rotates the pitch-class origin toward C.
 * @warning Requires positive chroma dimensions and octave width and an FFT window with at least two samples.
 */
void MIR::SpectralChromaInit() {
    const size_t nHalf = m_Config.FFTSize / 2 + 1;
    m_ChromaFilter.assign(m_Config.ChromaSize, std::vector<double>(nHalf, 0.0));

    const double tuning = 0.0;
    std::vector<double> frqbins(m_Config.FFTSize, 0.0);
    if (m_Config.FFTSize > 0) {
        for (size_t k = 1; k < m_Config.FFTSize; ++k) {
            const double frequency =
                static_cast<double>(k) * static_cast<double>(m_Config.SR) / static_cast<double>(m_Config.FFTSize);
            frqbins[k] = static_cast<double>(m_Config.ChromaSize) *
                         HzToOcts(frequency, tuning, static_cast<int>(m_Config.ChromaSize));
        }
        frqbins[0] = frqbins[1] - 1.5 * static_cast<double>(m_Config.ChromaSize);
    } else {
        spdlog::critical("FFT is smaller than 1, this should not happen");
        return;
    }

    std::vector<double> binwidthbins(m_Config.FFTSize, 1.0);
    for (size_t k = 0; k + 1 < m_Config.FFTSize; ++k) {
        binwidthbins[k] = std::max(frqbins[k + 1] - frqbins[k], 1.0);
    }

    const double nChroma2 = std::round(static_cast<double>(m_Config.ChromaSize) / 2.0);
    for (size_t k = 0; k < nHalf; ++k) {
        double columnNorm = 0.0;
        for (int chroma = 0; chroma < m_Config.ChromaSize; ++chroma) {
            const double distance = PositiveRemainder(frqbins[k] - static_cast<double>(chroma) + nChroma2 +
                                                          10.0 * static_cast<double>(m_Config.ChromaSize),
                                                      static_cast<double>(m_Config.ChromaSize)) -
                                    nChroma2;
            const double weight = std::exp(-0.5 * std::pow(2.0 * distance / binwidthbins[k], 2.0));
            m_ChromaFilter[chroma][k] = weight;
            columnNorm += weight * weight;
        }

        if (columnNorm > 0.0) {
            const double invNorm = 1.0 / std::sqrt(columnNorm);
            for (int chroma = 0; chroma < m_Config.ChromaSize; ++chroma) {
                m_ChromaFilter[chroma][k] *= invNorm;
            }
        }

        const double octaveWeight = std::exp(
            -0.5 * std::pow((frqbins[k] / static_cast<double>(m_Config.ChromaSize) - m_Config.ChromaCenterOctave) /
                                m_Config.ChromaOctaveWidth,
                            2.0));
        for (int chroma = 0; chroma < m_Config.ChromaSize; ++chroma) {
            m_ChromaFilter[chroma][k] *= octaveWeight;
        }
    }

    const int chromaShift = 3 * (m_Config.ChromaSize / 12);
    if (chromaShift > 0 && chromaShift < m_Config.ChromaSize) {
        std::vector<std::vector<double>> rolled(m_Config.ChromaSize, std::vector<double>(nHalf, 0.0));
        for (int chroma = 0; chroma < m_Config.ChromaSize; ++chroma) {
            rolled[chroma] = m_ChromaFilter[(chroma + chromaShift) % m_Config.ChromaSize];
        }
        m_ChromaFilter.swap(rolled);
    }
}

// ─────────────────────────────────────
/**
 * @brief Project power-spectrum bins onto the chroma filterbank.
 *
 * @param Desc Audio description to populate or update in place.
 *
 * @note Clears the output and accumulates unnormalized energy for each chroma bin.
 * @warning Requires an initialized nonempty filterbank and a Chroma array sized to the configured chroma count.
 */
void MIR::SpectralChromaExec(Description &Desc) {
    std::fill(Desc.Chroma.begin(), Desc.Chroma.end(), 0.0);
    const size_t nHalf = std::min(Desc.Power.size(), m_ChromaFilter[0].size());
    for (int chroma = 0; chroma < m_Config.ChromaSize; ++chroma) {
        double energy = 0.0;
        const auto &filter = m_ChromaFilter[chroma];

        for (size_t k = 0; k < nHalf; ++k) {
            energy += filter[k] * Desc.Power[k];
        }
        Desc.Chroma[chroma] = energy;
    }
}

// ╭─────────────────────────────────────╮
// │         Zero Crossing Rate          │
// ╰─────────────────────────────────────╯
/**
 * @brief Allocate scratch storage for zero-crossing analysis.
 *
 * @note Reserves extra edge-padding space when centered analysis is enabled.
 */
void MIR::ZeroCrossingRateInit() {
    const size_t frame = static_cast<size_t>(std::max(1.0f, m_Config.FFTSize));
    const size_t pad = m_Config.ZCRCenter ? (frame / 2) : 0;
    m_ZCRScratch.resize(frame + (2 * pad));
}

// ─────────────────────────────────────
/**
 * @brief Compute the configured zero-crossing rate.
 *
 * @param In Input audio frame containing the configured FFT window samples.
 * @param Desc Audio description to populate or update in place.
 *
 * @note Applies threshold, zero-sign, padding, and centering options and divides the crossing count by FFTSize.
 * @warning Requires a nonempty full-sized input frame and initialized scratch storage when centering is enabled.
 */
void MIR::ZeroCrossingRateExec(const std::vector<double> &In, Description &Desc) {
    const double *yData = nullptr;
    if (m_Config.ZCRCenter) {
        const size_t pad = m_Config.FFTSize / 2;
        const size_t inSize = In.size();
        double *dst = m_ZCRScratch.data();
        const double edgeLeft = In.front();
        const double edgeRight = In.back();
        std::fill_n(dst, pad, edgeLeft);
        std::copy(In.begin(), In.end(), dst + pad);
        std::fill_n(dst + pad + inSize, pad, edgeRight);
        yData = dst;
    } else {
        yData = In.data();
    }

    size_t crossings = 0;
    if (m_Config.ZCRPad) {
        crossings += 1;
    }

    const double threshold = m_Config.ZCRThreshold;
    if (m_Config.ZCRZeroPos) {
        double prev = yData[0];
        if (std::abs(prev) <= threshold) {
            prev = 0.0;
        }

        for (int i = 1; i < m_Config.FFTSize; ++i) {
            double curr = yData[static_cast<size_t>(i)];
            if (std::abs(curr) <= threshold)
                curr = 0.0;

            crossings += static_cast<size_t>(std::signbit(prev) != std::signbit(curr));
            prev = curr;
        }
    } else {
        double prev = yData[0];
        if (std::abs(prev) <= threshold) {
            prev = 0.0;
        }

        int prevSign = (prev > 0.0) - (prev < 0.0);
        for (int i = 1; i < m_Config.FFTSize; ++i) {
            double curr = yData[static_cast<size_t>(i)];
            if (std::abs(curr) <= threshold)
                curr = 0.0;

            const int currSign = (curr > 0.0) - (curr < 0.0);
            crossings += static_cast<size_t>(prevSign != currSign);
            prevSign = currSign;
        }
    }

    Desc.ZeroCrossingRate = static_cast<double>(crossings) / static_cast<double>(m_Config.FFTSize);
}

} // namespace OpenScofo
