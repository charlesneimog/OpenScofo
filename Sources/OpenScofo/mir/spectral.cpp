/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file spectral.cpp
 * @brief FFT processing and scalar spectral descriptor extraction.
 */

#include "mir.hpp"
#include <algorithm>
#include <cmath>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │          Set|Get Functions          │
// ╰─────────────────────────────────────╯
/**
 * @brief Allocate a real FFT setup and its aligned working buffers.
 *
 * @note Precomputes a periodic Hann window and reports allocation failures as critical logs.
 * @warning Requires a supported configured FFT size and previously released FFT resources.
 */
void MIR::FFTInit() {
    const size_t fftSize = static_cast<size_t>(m_Config.FFTSize);
    m_FullFFTSetup = pffft_new_setup(static_cast<int>(fftSize), PFFFT_REAL);
    if (!m_FullFFTSetup) {
        spdlog::critical("pffft_new_setup failed");
        return;
    }

    m_FullFFTIn = static_cast<float *>(pffft_aligned_malloc(fftSize * sizeof(float)));
    m_FullFFTOut = static_cast<float *>(pffft_aligned_malloc(fftSize * sizeof(float)));
    m_FullFFTWork = static_cast<float *>(pffft_aligned_malloc(fftSize * sizeof(float)));
    if (!m_FullFFTIn || !m_FullFFTOut || !m_FullFFTWork) {
        if (m_FullFFTIn) {
            pffft_aligned_free(m_FullFFTIn);
            m_FullFFTIn = nullptr;
        }
        if (m_FullFFTOut) {
            pffft_aligned_free(m_FullFFTOut);
            m_FullFFTOut = nullptr;
        }
        if (m_FullFFTWork) {
            pffft_aligned_free(m_FullFFTWork);
            m_FullFFTWork = nullptr;
        }
        pffft_destroy_setup(m_FullFFTSetup);
        m_FullFFTSetup = nullptr;
        spdlog::critical("pffft_aligned_malloc failed");
        return;
    }

    // Match librosa/scipy get_window('hann', N, fftbins=True): periodic Hann.
    m_FullWindowingFunc.resize(m_Config.FFTSize);
    for (size_t i = 0; i < m_Config.FFTSize; i++) {
        m_FullWindowingFunc[i] = 0.5 * (1.0 - cos(2.0 * std::numbers::pi * i / m_Config.FFTSize));
    }
}

// ╭─────────────────────────────────────╮
// │              SPECTRAL               │
// ╰─────────────────────────────────────╯
/**
 * @brief Compute the spectrum and accumulate all first-pass spectral sums.
 *
 * @param Desc Audio description to populate or update in place.
 * @param Sums Zero-initialized spectral sums for the current frame.
 * @param BinCount Number of nonnegative-frequency bins, including DC and Nyquist.
 *
 * @note Keeps spectrum generation, moments, flux, irregularity, and related sums in one fused pass.
 */
void MIR::ComputeSpectrum(Description &Desc, SpectralAccumulators &Sums, int BinCount) {
    const double binWidth = static_cast<double>(m_Config.SR) / static_cast<double>(m_Config.FFTSize);
    const double fftScale = 1.0 / static_cast<double>(m_Config.FFTSize);
    const double minimumPower = 1e-10;
    const int highFrequencyStartBin = BinCount / 4;
    const size_t nyquistBin = static_cast<size_t>(m_Config.FFTSize) / 2;

    Desc.MaxAmp = 0.0;
    Desc.SpectralFlux = 0.0;

    double previousNormalizedMagnitude = 0.0;
    double magnitudeTwoBinsBack = 0.0;
    double previousMagnitude = 0.0;

    // Handle DC separately so the main loop only needs the Nyquist special case.
    {
        const double real = static_cast<double>(m_FullFFTOut[0]);
        const double imaginary = 0.0;
        // const double power = real * real + imaginary * imaginary;
        const float real32 = static_cast<float>(real);
        const float imaginary32 = static_cast<float>(imaginary);
        const float power32 = real32 * real32 + imaginary32 * imaginary32;
        const double power = static_cast<double>(power32);
        const double magnitude = std::sqrt(power);
        const double normalizedMagnitude = magnitude * fftScale;

        Desc.Power[0] = power;
        Desc.Magnitude[0] = magnitude;
        Desc.SpectralMagnitudeNorm[0] = normalizedMagnitude;
        Desc.MaxAmp = normalizedMagnitude;

        Sums.sumMagCrest += static_cast<float>(magnitude);
        Sums.maxMagCrest = static_cast<float>(magnitude);

        Sums.SumPower += normalizedMagnitude;
        Sums.irregularityDenominator += normalizedMagnitude * normalizedMagnitude;

        const double flatnessPower = std::max(1e-10, power);
        const double logPower = std::log(flatnessPower);

        Sums.logSumPower += std::log(flatnessPower);
        Sums.linSumPower += flatnessPower;
        Sums.spectralEnergySum += power;

        if (power > 0.0) {
            Sums.weightedLogPower += power * (power >= minimumPower ? logPower : std::log(power));
        }

        const double fluxDifference = magnitude - m_PreviousSpectralMagnitude[0];
        Desc.SpectralFlux += fluxDifference * fluxDifference;
        m_PreviousSpectralMagnitude[0] = magnitude;

        previousNormalizedMagnitude = normalizedMagnitude;
        previousMagnitude = magnitude;
    }

    for (int i = 1; i < BinCount; ++i) {
        double real = 0.0;
        double imaginary = 0.0;
        const size_t bin = static_cast<size_t>(i);
        if (bin == nyquistBin) {
            real = static_cast<double>(m_FullFFTOut[1]);
        } else {
            const size_t outputIndex = 2 * bin;
            real = static_cast<double>(m_FullFFTOut[outputIndex]);
            imaginary = static_cast<double>(m_FullFFTOut[outputIndex + 1]);
        }

        const double power = real * real + imaginary * imaginary;
        const double magnitude = std::sqrt(power);
        const double normalizedMagnitude = magnitude * fftScale;

        Desc.Power[i] = power;
        Desc.Magnitude[i] = magnitude;
        Desc.SpectralMagnitudeNorm[i] = normalizedMagnitude;
        Desc.MaxAmp = std::max(Desc.MaxAmp, normalizedMagnitude);

        const float crestMagnitude = static_cast<float>(magnitude);
        Sums.sumMagCrest += crestMagnitude;
        Sums.maxMagCrest = std::max(Sums.maxMagCrest, crestMagnitude);

        const double frequency = i * binWidth;
        const double frequencySquared = frequency * frequency;
        const double binIndex = static_cast<double>(i);

        Sums.SumPower += normalizedMagnitude;
        Sums.sumFreq += frequency;
        Sums.sumFreqSq += frequencySquared;
        Sums.sumFreq2 += frequencySquared * normalizedMagnitude;
        Sums.sumFreq3 += frequencySquared * frequency * normalizedMagnitude;
        Sums.sumFreq4 += frequencySquared * frequencySquared * normalizedMagnitude;
        Sums.WSumFreqs += frequency * normalizedMagnitude;
        Sums.sumIndex += binIndex * normalizedMagnitude;
        Sums.sumIndex2 += binIndex * binIndex * normalizedMagnitude;
        Sums.irregularityDenominator += normalizedMagnitude * normalizedMagnitude;

        const double flatnessPower = std::max(minimumPower, power);
        const double logPower = std::log(flatnessPower);

        Sums.logSumPower += logPower;
        Sums.linSumPower += flatnessPower;
        Sums.spectralEnergySum += power;

        if (power > 0.0) {
            Sums.weightedLogPower += power * (power >= minimumPower ? logPower : std::log(power));
        }

        const double fluxDifference = magnitude - m_PreviousSpectralMagnitude[i];
        Desc.SpectralFlux += fluxDifference * fluxDifference;
        m_PreviousSpectralMagnitude[i] = magnitude;

        // The comparison becomes 1.0 above the high-frequency boundary and 0.0 below it.
        Sums.highFreqEnergy += normalizedMagnitude * static_cast<double>(i >= highFrequencyStartBin);

        const double jensenDifference = previousNormalizedMagnitude - normalizedMagnitude;
        Sums.irregularityJensenNumerator += jensenDifference * jensenDifference;
        if (i > 1) {
            const double localMagnitudeMean = (magnitudeTwoBinsBack + previousMagnitude + magnitude) / 3.0;
            Sums.irregularityKrimphoffSum += std::abs(previousMagnitude - localMagnitudeMean);
        }

        Sums.harmonicitySum += normalizedMagnitude;
        Sums.harmonicityPeak = std::max(Sums.harmonicityPeak, normalizedMagnitude);

        previousNormalizedMagnitude = normalizedMagnitude;
        magnitudeTwoBinsBack = previousMagnitude;
        previousMagnitude = magnitude;
    }
}

// ─────────────────────────────────────
/**
 * @brief Finalize scalar spectral features from accumulated moments.
 *
 * @param Desc Audio description to populate or update in place.
 * @param Sums Spectral sums and moments accumulated over the current FFT bins.
 * @param BinCount Number of nonnegative-frequency bins, including DC and Nyquist.
 *
 * @note Updates centroid history, spread, shape, irregularity, crest, flatness, harmonicity, and high-frequency
 * ratio.
 * @warning BinCount must be positive and the sums must describe the current frame.
 */
void MIR::ComputeScalarFeatures(Description &Desc, const SpectralAccumulators &Sums, int BinCount) {
    Desc.SpectralFlux = std::sqrt(Desc.SpectralFlux);

    const double magnitudeSumWithEpsilon = Sums.SumPower + 1e-12;
    const double inverseMagnitudeSum = 1.0 / magnitudeSumWithEpsilon;

    Desc.SpectralCentroid = Sums.WSumFreqs * inverseMagnitudeSum;
    Desc.CentroidVelocity = std::abs(Desc.SpectralCentroid - m_PrevCentroid);
    m_PrevCentroid = Desc.SpectralCentroid;

    const double meanFrequency = Desc.SpectralCentroid;
    const double meanSquaredFrequency = Sums.sumFreq2 * inverseMagnitudeSum;
    const double meanCubedFrequency = Sums.sumFreq3 * inverseMagnitudeSum;
    const double meanFourthPowerFrequency = Sums.sumFreq4 * inverseMagnitudeSum;
    const double thirdCentralMoment = meanCubedFrequency - 3.0 * meanFrequency * meanSquaredFrequency +
                                      2.0 * meanFrequency * meanFrequency * meanFrequency;
    const double fourthCentralMoment = meanFourthPowerFrequency - 4.0 * meanFrequency * meanCubedFrequency +
                                       6.0 * meanFrequency * meanFrequency * meanSquaredFrequency -
                                       3.0 * meanFrequency * meanFrequency * meanFrequency * meanFrequency;
    const double frequencyVariance = std::max(0.0, meanSquaredFrequency - meanFrequency * meanFrequency);
    Desc.SpectralSpreadHz = Sums.SumPower > 0.0 ? std::sqrt(frequencyVariance) : 0.0;

    // Essentia DistributionShape spread: normalized second central moment over FFT bin indices.
    const double meanBinIndex = Sums.sumIndex * inverseMagnitudeSum;
    const double meanSquaredBinIndex = Sums.sumIndex2 * inverseMagnitudeSum;
    const double binIndexVariance = std::max(0.0, meanSquaredBinIndex - meanBinIndex * meanBinIndex);
    const double binIndexRange = static_cast<double>(std::max(1, BinCount - 1));
    const double inverseBinIndexRange = 1.0 / binIndexRange;
    Desc.SpectralSpreadVariance =
        Sums.SumPower > 0.0 ? (binIndexVariance * inverseBinIndexRange * inverseBinIndexRange) : 0.0;

    Desc.SpectralSkewness = 0.0;
    Desc.SpectralKurtosis = 0.0;

    if (Desc.SpectralSpreadHz > 1e-9 && Sums.SumPower > 1e-12) {
        const double spread = Desc.SpectralSpreadHz;
        const double spreadSquared = spread * spread;
        Desc.SpectralSkewness = thirdCentralMoment / (spreadSquared * spread);
        Desc.SpectralKurtosis = fourthCentralMoment / (spreadSquared * spreadSquared) - 3.0;
    }

    Desc.SpectralIrregularityJensen =
        Sums.irregularityDenominator > 0.0 ? (Sums.irregularityJensenNumerator / Sums.irregularityDenominator) : -1.0;
    Desc.SpectralIrregularityKrimphoff =
        Sums.irregularityKrimphoffSum > 0.0 ? std::log10(Sums.irregularityKrimphoffSum) : -1.0;
    Desc.SpectralIrregularity = Desc.SpectralIrregularityKrimphoff;

    if (Sums.maxMagCrest == 0.0f) {
        Desc.SpectralCrest = 0.0;
    } else {
        const float meanCrestMagnitude = Sums.sumMagCrest / static_cast<float>(BinCount);
        Desc.SpectralCrest =
            meanCrestMagnitude > 0.0f ? static_cast<double>(Sums.maxMagCrest / meanCrestMagnitude) : 0.0;
    }

    Desc.SpectralFlatness = std::exp(Sums.logSumPower / BinCount) / (Sums.linSumPower / BinCount);
    Desc.Harmonicity = Sums.harmonicityPeak / (Sums.harmonicitySum + 1e-12);
    Desc.HighFreqRatio = Sums.highFreqEnergy * inverseMagnitudeSum;

    const double binCount = static_cast<double>(BinCount);
    const double slopeDenominator = (binCount * Sums.sumFreqSq - Sums.sumFreq * Sums.sumFreq) + 1e-12;
    const double slope = (binCount * Sums.WSumFreqs - Sums.sumFreq * Sums.SumPower) / slopeDenominator;
    Desc.SpectralSlope = slope * inverseMagnitudeSum;
}

// ─────────────────────────────────────
/**
 * @brief Normalize the spectrum and compute all second-pass statistics.
 *
 * @param Desc Audio description to populate or update in place.
 * @param Sums Spectral sums accumulated during the first pass.
 * @param BinCount Number of nonnegative-frequency bins, including DC and Nyquist.
 *
 * @note Keeps frame normalization, variance, entropy, rolloff, and prefix sums in one fused pass.
 */
void MIR::NormalizeSpectrum(Description &Desc, const SpectralAccumulators &Sums, int BinCount) {
    const double magnitudeSumWithEpsilon = Sums.SumPower + 1e-12;
    const double inverseMagnitudeSum = 1.0 / magnitudeSumWithEpsilon;
    const double meanNormalizedMagnitude = 1.0 / static_cast<double>(BinCount);
    const double inverseSpectralEnergy = Sums.spectralEnergySum > 0.0 ? (1.0 / Sums.spectralEnergySum) : 0.0;
    const double rolloffCutoffPower = std::clamp(m_Config.SpectralRolloffCutoff, 0.0, 1.0) * Sums.spectralEnergySum;

    double normalizedMagnitudeVariance = 0.0;
    double cumulativePower = 0.0;
    size_t rolloffBin = 0;
    bool rolloffFound = false;

    Desc.SpectralEntropy = 0.0;
    m_SpectralPrefix[0] = 0.0;

    for (int i = 0; i < BinCount; ++i) {
        const double frameNormalizedMagnitude = (Desc.SpectralMagnitudeNorm[i] + 1e-12) * inverseMagnitudeSum;
        Desc.SpectralMagnitudeFrameNorm[i] = frameNormalizedMagnitude;
        Desc.SpectralPowerFrameNorm[i] = Desc.Power[i] * inverseSpectralEnergy;

        const double differenceFromMean = frameNormalizedMagnitude - meanNormalizedMagnitude;
        normalizedMagnitudeVariance += differenceFromMean * differenceFromMean;

        cumulativePower += Desc.Power[i];
        m_SpectralPrefix[i + 1] = cumulativePower;

        if (!rolloffFound && cumulativePower >= rolloffCutoffPower) {
            rolloffBin = i;
            rolloffFound = true;
        }
    }

    constexpr double invLn2 = 1.4426950408889634;

    Desc.SpectralEntropy = Sums.spectralEnergySum > 0.0 ? std::log2(Sums.spectralEnergySum) -
                                                              Sums.weightedLogPower * inverseSpectralEnergy * invLn2
                                                        : 0.0;

    Desc.StdDev = std::sqrt(normalizedMagnitudeVariance * meanNormalizedMagnitude);
    const double frequencyPerBin = (m_Config.SR * 0.5) / std::max(1, BinCount - 1);
    Desc.SpectralRolloff = static_cast<double>(rolloffBin) * frequencyPerBin;
}

// ─────────────────────────────────────
/**
 * @brief Transform the windowed input and compute spectral descriptors.
 *
 * @param Desc Audio description to populate or update in place.
 *
 * @note Executes one fused accumulation pass and one fused normalization pass after the FFT.
 * @warning Requires initialized FFT resources and descriptor arrays sized to FFTSize / 2 + 1.
 */
void MIR::GetSpectralDescriptions(Description &Desc) {
    const int binCount = m_Config.FFTSize / 2 + 1;
    pffft_transform_ordered(m_FullFFTSetup, m_FullFFTIn, m_FullFFTOut, m_FullFFTWork, PFFFT_FORWARD);

    SpectralAccumulators sums;
    ComputeSpectrum(Desc, sums, binCount);
    ComputeScalarFeatures(Desc, sums, binCount);
    NormalizeSpectrum(Desc, sums, binCount);
}

// ─────────────────────────────────────
/**
 * @brief Clear the reverberation spectrum used by pitch evidence.
 *
 * @param Desc Audio description to populate or update in place.
 * @param decay Reserved reverberation decay parameter; currently ignored.
 *
 * @note Reverberation accumulation is currently disabled; the decay argument is ignored.
 * @warning ReverbSpectralPower must cover every bin in SpectralMagnitudeFrameNorm.
 */
void MIR::AddReverb(Description &Desc, double decay) {
    (void)decay;
    for (size_t i = 0; i < Desc.SpectralMagnitudeFrameNorm.size(); i++) {
        Desc.ReverbSpectralPower[i] =
            0; //(Desc.ReverbSpectralPower[i] * decay) * (Desc.SpectralMagnitudeFrameNorm[i] * decay);
    }
}

} // namespace OpenScofo
