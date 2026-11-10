/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file pitch.cpp
 * @brief YIN fundamental-frequency and pitch-confidence estimation.
 */

#include "mir.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │                Pitch                │
// ╰─────────────────────────────────────╯
/**
 * @brief Allocate YIN difference and normalization scratch arrays.
 *
 * @note Sizes the arrays from half the configured FFT window plus interpolation padding.
 */
void MIR::YINInit() {
    const size_t frameSize = static_cast<size_t>(std::max(2.0f, m_Config.FFTSize));
    const size_t half = frameSize / 2;
    const size_t allocSize = half + 2;
    m_YINDifference.assign(allocSize, 0.0);
    m_YINCMNDF.assign(allocSize, 1.0);
}

// ─────────────────────────────────────
/**
 * @brief Estimate fundamental frequency and confidence using YIN.
 *
 * @param In Input audio frame containing the configured FFT window samples.
 * @param Desc Audio description to populate or update in place.
 *
 * @note Uses a normalized difference function and parabolic refinement; unusable estimates produce zero pitch and
 * confidence.
 * @warning Initialize YIN scratch arrays and use a positive sample rate and valid positive frequency bounds.
 */
void MIR::YINExec(const std::vector<double> &In, Description &Desc) {
    const size_t frame = In.size();

    if (frame < 2) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    const size_t minTau = static_cast<size_t>(m_Config.SR / m_Config.YINMaxFrequency);
    const size_t maxTauByPitch = static_cast<size_t>(std::ceil(m_Config.SR / m_Config.YINMinFrequency));
    const size_t maxTau = std::min({frame / 2, m_YINDifference.size() - 1, maxTauByPitch});

    if (maxTau <= minTau) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    double *const diff = m_YINDifference.data();
    double *const cmndf = m_YINCMNDF.data();
    const double *const data = In.data();

    std::fill_n(diff, maxTau + 1, 0.0);

    // ─────────────────────────────────────
    // YIN difference function
    //
    // d(tau) = sum_i (x_i - x_{i+tau})²
    //
    // Keep tau as the inner loop because both diff[tau]
    // and data[i + tau] are then accessed sequentially.
    // TODO: LOOP
    for (size_t i = 0; i + minTau < frame; ++i) {
        const double x = data[i];

        const size_t limit = std::min(maxTau, frame - i - 1);

        double *d = diff + minTau;
        const double *lag = data + i + minTau;
        const double *const lagEnd = data + i + limit + 1;

        // TODO: LOOP
        while (lag < lagEnd) {
            const double delta = x - *lag++;
            *d++ += delta * delta;
        }
    }

    // ─────────────────────────────────────
    // Cumulative mean normalized difference
    //
    // diff[1 ... minTau-1] is zero in the existing
    // implementation, so those iterations can be skipped.
    cmndf[0] = 1.0;

    if (minTau > 1) {
        std::fill(cmndf + 1, cmndf + minTau, 1.0);
    }

    double cumulative = 0.0;
    // TODO: LOOP
    for (size_t tau = minTau; tau <= maxTau; ++tau) {
        cumulative += diff[tau];
        if (cumulative <= 0.0) {
            cmndf[tau] = 1.0;
        } else {
            cmndf[tau] = diff[tau] * static_cast<double>(tau) / cumulative;
        }
    }

    // ─────────────────────────────────────
    // First look for YIN's threshold crossing.
    //
    // If one exists, the original global-minimum scan is
    // unnecessary because it is overwritten afterward.
    size_t tauEstimate = 0;
    double bestValue = std::numeric_limits<double>::infinity();

    // TODO: LOOP
    for (size_t tau = minTau; tau <= maxTau; ++tau) {
        const double value = cmndf[tau];

        if (value < m_Config.YINThreshold) {
            tauEstimate = tau;
            bestValue = value;

            // TODO: LOOP
            while (tau + 1 <= maxTau && cmndf[tau + 1] < bestValue) {
                tau++;
                tauEstimate = tau;
                bestValue = cmndf[tau];
            }

            break;
        }
    }

    // No threshold crossing: fall back to the global minimum,
    // exactly as the previous implementation did.
    if (tauEstimate == 0) {
        // TODO: LOOP
        for (size_t tau = minTau; tau <= maxTau; ++tau) {
            const double value = cmndf[tau];

            if (value < bestValue) {
                bestValue = value;
                tauEstimate = tau;
            }
        }
    }

    if (tauEstimate == 0 || !std::isfinite(bestValue)) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    // ─────────────────────────────────────
    // Parabolic interpolation
    double RefinedTau = static_cast<double>(tauEstimate);
    if (tauEstimate > minTau && tauEstimate + 1 <= maxTau) {
        const double left = cmndf[tauEstimate - 1];
        const double center = cmndf[tauEstimate];
        const double right = cmndf[tauEstimate + 1];
        const double denominator = left - 2.0 * center + right;
        if (std::abs(denominator) > 1e-12) {
            const double offset = 0.5 * (left - right) / denominator;
            RefinedTau += std::clamp(offset, -1.0, 1.0);
        }
    }

    const double Confidence = std::clamp(1.0 - bestValue, 0.0, 1.0);
    if (RefinedTau <= 0.0 || Confidence <= 0.0) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    const double Pitch = m_Config.SR / RefinedTau;
    if (Pitch < m_Config.YINMinFrequency || Pitch > m_Config.YINMaxFrequency) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    Desc.Pitch = Pitch;
    Desc.PitchConfidence = Confidence;
}

} // namespace OpenScofo
