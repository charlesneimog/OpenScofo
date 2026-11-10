/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.

    The harmonic pitch-voting algorithm is adapted from Miller Puckette's
    BSD-licensed sigmund~ implementation (2005).
*/

/**
 * @file sigmund.cpp
 * @brief Sigmund-style pitch and confidence estimation over the shared FFT.
 */

#include "mir.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>

namespace OpenScofo {
namespace {

constexpr size_t kPeakCount = 20;
constexpr size_t kSalientPeakCount = 12;
constexpr size_t kHarmonicCount = 16;
constexpr int kStepsPerOctave = 48;
constexpr float kAmplitudePowerLaw = 0.5f;
constexpr float kQualityThreshold = 0.4f;
constexpr double kPeakThresholdFactor = 0.6;
constexpr float kPi = std::numbers::pi_v<float>;
constexpr float kLogTwo = std::numbers::ln2_v<float>;

struct SigmundPeak {
    float Frequency = 0.0f;
    float Amplitude = 0.0f;
    float Decibels = 0.0f;
    bool Selected = false;
};

struct SigmundEstimate {
    float Frequency = 0.0f;
    float Quality = 0.0f;
};

int IntegerLogTwo(size_t Value) {
    int Result = -1;
    while (Value != 0) {
        Value >>= 1U;
        ++Result;
    }
    return Result;
}

float SinX(float Theta, float SinTheta) {
    if (Theta > -0.003f && Theta < 0.003f) {
        return 1.0f;
    }
    return SinTheta / Theta;
}

float HannWindowMagnitude(float PiDetune, float SinPiDetune) {
    return 0.5f * SinX(PiDetune, SinPiDetune) -
           0.25f * (SinX(PiDetune + kPi, SinPiDetune) + SinX(PiDetune - kPi, SinPiDetune));
}

float PowerToDecibels(float Value) {
    if (Value <= 0.0f) {
        return 0.0f;
    }
    return std::max(0.0f, 100.0f + 10.0f * std::log10(Value));
}

void RemaskPeak(std::span<double> Mask, int MaxBin, int BestIndex, double MaximumPower) {
    // The reference mask measures distance on a half-bin-spaced spectrum.
    // Scale by two for the shared N-point FFT and cover the same bandwidth.
    const int FirstBin = std::max(1, BestIndex - 25);
    const int LastBin = std::min(MaxBin, BestIndex + 26);
    for (int Bin = FirstBin; Bin < LastBin; ++Bin) {
        const int BinDistance = Bin - BestIndex;
        const double ReferenceDistance = 2.0 * static_cast<double>(BinDistance);
        const double DistanceSquared = ReferenceDistance * ReferenceDistance;
        double NewMask = MaximumPower / (1.0 + DistanceSquared * DistanceSquared);
        if (std::abs(BinDistance) <= 1) {
            NewMask = 2.0 * MaximumPower;
        }
        Mask[static_cast<size_t>(Bin)] = std::max(Mask[static_cast<size_t>(Bin)], NewMask);
    }
}

size_t ExtractPeaks(size_t FrameSize, float SampleRate, std::span<const double> Power, std::span<double> Mask,
                    std::span<SigmundPeak, kPeakCount> Peaks) {
    std::fill(Mask.begin(), Mask.end(), 0.0);
    if (Power.size() < 3 || Mask.size() < Power.size() || FrameSize == 0 || !(SampleRate > 0.0f)) {
        return 0;
    }

    const int MaxBin = static_cast<int>(Power.size()) - 1; // Exclude Nyquist so both neighbors exist.
    const double BinWidth = static_cast<double>(SampleRate) / static_cast<double>(FrameSize);
    double TotalPower = 0.0;
    for (int Bin = 1; Bin < MaxBin; ++Bin) {
        const double BinPower = Power[static_cast<size_t>(Bin)];
        if (std::isfinite(BinPower) && BinPower > 0.0) {
            TotalPower += BinPower;
        }
    }

    size_t FoundCount = 0;
    while (FoundCount < Peaks.size()) {
        double MaximumPower = 0.0;
        int BestIndex = -1;
        for (int Bin = 1; Bin < MaxBin; ++Bin) {
            const size_t Index = static_cast<size_t>(Bin);
            const double BinPower = Power[Index];
            if (!(BinPower > MaximumPower) || !std::isfinite(BinPower) || !(BinPower > Mask[Index])) {
                continue;
            }

            const double LeftPower = std::max(0.0, Power[Index - 1]);
            const double RightPower = std::max(0.0, Power[Index + 1]);
            const bool IsLocalMaximum = BinPower >= LeftPower && BinPower > RightPower;
            if (IsLocalMaximum && BinPower > kPeakThresholdFactor * (LeftPower + RightPower)) {
                MaximumPower = BinPower;
                BestIndex = Bin;
            }
        }

        if (!(TotalPower > 0.0) || MaximumPower < 1.0e-10 * TotalPower || BestIndex < 0) {
            break;
        }

        RemaskPeak(Mask, MaxBin, BestIndex, MaximumPower);

        const size_t Index = static_cast<size_t>(BestIndex);
        const double LeftMagnitude = std::sqrt(std::max(0.0, Power[Index - 1]));
        const double CenterMagnitude = std::sqrt(MaximumPower);
        const double RightMagnitude = std::sqrt(std::max(0.0, Power[Index + 1]));
        const double Denominator = LeftMagnitude + 2.0 * CenterMagnitude + RightMagnitude;
        double Detune = 0.0;
        if (Denominator > 0.0 && std::isfinite(Denominator)) {
            Detune = 2.0 * (RightMagnitude - LeftMagnitude) / Denominator;
            Detune = std::clamp(Detune, -0.5, 0.5);
        }

        const float PiDetune = kPi * static_cast<float>(Detune);
        const float Gain = HannWindowMagnitude(PiDetune, std::sin(PiDetune));
        if (!(Gain > 0.0f) || !std::isfinite(Gain)) {
            continue;
        }

        const double Amplitude = CenterMagnitude / (static_cast<double>(FrameSize) * Gain);
        const double Frequency = (static_cast<double>(BestIndex) + Detune) * BinWidth;
        if (!std::isfinite(Amplitude) || !(Amplitude > 0.0) || !std::isfinite(Frequency) || !(Frequency > 0.0)) {
            continue;
        }

        SigmundPeak &Peak = Peaks[FoundCount++];
        Peak.Frequency = static_cast<float>(Frequency);
        Peak.Amplitude = static_cast<float>(Amplitude);
        Peak.Decibels = PowerToDecibels(Peak.Amplitude);
    }
    return FoundCount;
}

SigmundEstimate EstimatePitch(size_t FrameSize, float SampleRate, std::span<SigmundPeak> Peaks,
                              std::span<std::array<float, 2>> PitchPoints,
                              const std::array<float, kHarmonicCount> &HarmonicWeights) {
    SigmundEstimate Result;
    if (Peaks.empty() || PitchPoints.size() < 2 || FrameSize == 0 || !(SampleRate > 0.0f)) {
        return Result;
    }

    for (auto &Point : PitchPoints) {
        Point.fill(0.0f);
    }
    for (SigmundPeak &Peak : Peaks) {
        Peak.Selected = false;
    }

    std::array<SigmundPeak *, kSalientPeakCount> SalientPeaks{};
    size_t SalientCount = 0;
    while (SalientCount < SalientPeaks.size()) {
        SigmundPeak *BestPeak = nullptr;
        float BestSalience = -std::numeric_limits<float>::infinity();
        for (SigmundPeak &Peak : Peaks) {
            if (!Peak.Selected && Peak.Decibels > BestSalience) {
                BestSalience = Peak.Decibels;
                BestPeak = &Peak;
            }
        }
        if (BestPeak == nullptr) {
            break;
        }
        BestPeak->Selected = true;
        SalientPeaks[SalientCount++] = BestPeak;
    }

    const float BinWidth = SampleRate / static_cast<float>(FrameSize);
    float SumLoudness = 0.0f;
    for (size_t PeakIndex = 0; PeakIndex < SalientCount; ++PeakIndex) {
        const SigmundPeak &Peak = *SalientPeaks[PeakIndex];
        const float WeightIndex = static_cast<float>(kStepsPerOctave) / kLogTwo * std::log(Peak.Frequency / BinWidth);
        const float Loudness = std::pow(Peak.Amplitude, kAmplitudePowerLaw);
        for (size_t Harmonic = 0; Harmonic < HarmonicWeights.size(); ++Harmonic) {
            int LowIndex = static_cast<int>(
                WeightIndex -
                static_cast<float>(kStepsPerOctave) / kLogTwo * std::log(static_cast<float>(Harmonic + 1)) - 0.5f);
            const int HighIndex = LowIndex + 2;
            if (HighIndex < 0) {
                break;
            }
            if (HighIndex >= static_cast<int>(PitchPoints.size())) {
                continue;
            }
            LowIndex = std::max(0, LowIndex);
            for (int Index = LowIndex; Index <= HighIndex; ++Index) {
                auto &Point = PitchPoints[static_cast<size_t>(Index)];
                Point[0] += Loudness * HarmonicWeights[Harmonic];
                Point[1] += Loudness;
            }
        }
        SumLoudness += Loudness;
    }

    size_t BestBin = 0;
    float BestWeight = -std::numeric_limits<float>::infinity();
    for (size_t Index = 0; Index < PitchPoints.size(); ++Index) {
        if (PitchPoints[Index][0] > BestWeight) {
            BestWeight = PitchPoints[Index][0];
            BestBin = Index;
        }
    }

    if (BestBin + 1 >= PitchPoints.size() || !(SumLoudness > 0.0f) || !(PitchPoints[BestBin][1] > 0.0f) ||
        PitchPoints[BestBin][1] < SumLoudness * kQualityThreshold) {
        return Result;
    }

    Result.Quality = PitchPoints[BestBin][1] / SumLoudness;
    const size_t PreviousBin = BestBin == 0 ? 0 : BestBin - 1;
    const float InterpolationDenominator =
        PitchPoints[BestBin + 1][0] + PitchPoints[BestBin][0] + PitchPoints[PreviousBin][0];
    float RefinedBin = static_cast<float>(BestBin);
    if (InterpolationDenominator != 0.0f) {
        RefinedBin += (PitchPoints[BestBin + 1][0] - PitchPoints[PreviousBin][0]) / InterpolationDenominator;
    }

    float Frequency = BinWidth * std::exp(kLogTwo / static_cast<float>(kStepsPerOctave) * RefinedBin);
    float SumWeight = 0.0f;
    float SumFrequency = 0.0f;
    for (size_t PeakIndex = 0; PeakIndex < SalientCount; ++PeakIndex) {
        const SigmundPeak &Peak = *SalientPeaks[PeakIndex];
        const float Harmonic = Peak.Frequency / Frequency;
        const int RoundedHarmonic = static_cast<int>(std::lround(Harmonic));
        const float Inharmonicity = Harmonic - static_cast<float>(RoundedHarmonic);
        if (RoundedHarmonic >= 1 && RoundedHarmonic <= static_cast<int>(kHarmonicCount) &&
            std::abs(Inharmonicity) < 0.015f * static_cast<float>(RoundedHarmonic)) {
            const float Weight = Peak.Amplitude * static_cast<float>(RoundedHarmonic);
            SumWeight += Weight;
            SumFrequency += Weight * Peak.Frequency / static_cast<float>(RoundedHarmonic);
        }
    }
    if (SumWeight > 0.0f) {
        Frequency = SumFrequency / SumWeight;
    }

    if (std::isfinite(Frequency) && Frequency > 0.0f && std::isfinite(Result.Quality)) {
        Result.Frequency = Frequency;
        Result.Quality = std::clamp(Result.Quality, 0.0f, 1.0f);
    } else {
        Result = {};
    }
    return Result;
}

} // namespace

// ─────────────────────────────────────
/**
 * @brief Preallocate Sigmund's peak mask and harmonic-voting scratch storage.
 *
 * @note Sigmund shares the main N-point Hann spectrum and owns no FFT setup or transform buffers.
 */
void MIR::SigmundInit() {
    const size_t FrameSize = static_cast<size_t>(m_Config.FFTSize);
    if (FrameSize < 2) {
        m_SigmundMask.clear();
        m_SigmundPitchPoints.clear();
        return;
    }

    m_SigmundMask.assign(FrameSize / 2 + 1, 0.0);
    const int PitchOctaves = IntegerLogTwo(FrameSize);
    m_SigmundPitchPoints.assign(static_cast<size_t>(kStepsPerOctave) * static_cast<size_t>(PitchOctaves), {});
    for (size_t Harmonic = 0; Harmonic < m_SigmundHarmonicWeights.size(); ++Harmonic) {
        m_SigmundHarmonicWeights[Harmonic] = 6.0f / (6.0f + static_cast<float>(Harmonic));
    }
}

// ─────────────────────────────────────
/**
 * @brief Estimate pitch and confidence from the current shared FFT spectrum.
 *
 * @param Desc Audio description containing the current Hann-windowed power spectrum and receiving the result.
 *
 * @note Performs no FFT, zero-padding, peak output, or allocation after initialization.
 */
void MIR::SigmundExec(Description &Desc) {
    Desc.Pitch = 0.0;
    Desc.PitchConfidence = 0.0;

    const size_t FrameSize = static_cast<size_t>(m_Config.FFTSize);
    const size_t BinCount = FrameSize / 2 + 1;
    const size_t PitchPointCount =
        static_cast<size_t>(kStepsPerOctave) * static_cast<size_t>(std::max(0, IntegerLogTwo(FrameSize)));
    if (FrameSize < 2 || !(m_Config.SR > 0.0f) || Desc.Power.size() < BinCount || m_SigmundMask.size() != BinCount ||
        m_SigmundPitchPoints.size() != PitchPointCount || m_SigmundPitchPoints.empty()) {
        return;
    }

    std::array<SigmundPeak, kPeakCount> Peaks{};
    const size_t FoundCount = ExtractPeaks(FrameSize, m_Config.SR, std::span<const double>(Desc.Power.data(), BinCount),
                                           m_SigmundMask, Peaks);
    const SigmundEstimate Estimate =
        EstimatePitch(FrameSize, m_Config.SR, std::span<SigmundPeak>(Peaks.data(), FoundCount), m_SigmundPitchPoints,
                      m_SigmundHarmonicWeights);

    if (Estimate.Frequency < m_Config.YINMinFrequency || Estimate.Frequency > m_Config.YINMaxFrequency) {
        return;
    }
    Desc.Pitch = Estimate.Frequency;
    Desc.PitchConfidence = Estimate.Quality;
}

} // namespace OpenScofo
