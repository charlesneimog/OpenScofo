/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file pitch.cpp
 * @brief Spectral pitch templates and pitch-probability evaluation.
 */
#include "OpenScofo.hpp"
#include <algorithm>

namespace OpenScofo {

// ─────────────────────────────────────
// GONG 2015 (adapted)
/**
 * @brief Cache a normalized harmonic spectral template for a pitch.
 *
 * @param Freq Pitch frequency in Hz.
 *
 * @note Uses stretched Gaussian harmonic peaks and a positive floor, keyed by the rounded fundamental FFT bin.
 * @warning Requires positive audio configuration values; frequencies outside (0, Nyquist] are ignored.
 */
void OnlineForward::BuildPitchTemplate(double Freq) {
    const double m_MinHarmonicDecay = 0.2;
    const double m_MaxHarmonicDecay = 1.8;
    const double binWidth = m_Sr / m_FFTSize;

    if (Freq > m_Sr / 2.0 || Freq <= 0)
        return;

    double m_MinF0 = 32.70;
    double m_MaxF0 = 4186.0;
    double f0Norm = std::log2(Freq / m_MinF0) / std::log2(m_MaxF0 / m_MinF0);
    f0Norm = std::clamp(f0Norm, 0.0, 1.0);

    double pitchDecayCurve = m_PitchTemplateAmplitudeDecay;
    double shaped = std::pow(1.0 - f0Norm, pitchDecayCurve);
    double beta = m_MaxHarmonicDecay - shaped * (m_MaxHarmonicDecay - m_MinHarmonicDecay);

    double rootBinFreq = std::round(Freq / binWidth);

    if (m_PitchTemplates.find(rootBinFreq) != m_PitchTemplates.end()) {
        return;
    }

    // Keep a tiny positive floor across the template. Besides avoiding
    // log(0), this is part of the original KL-divergence model: bins outside
    // the harmonic peaks must not be treated as unmatched observed energy.
    std::vector<double> templateBins(m_FFTSize / 2, 1e-24);

    const double sigmaLog = m_PitchTemplateSigma / 12.0;
    const double sigmaConst = std::pow(2.0, sigmaLog) - 1.0;

    const double B = 0.0001;
    for (int k = 1; k <= m_Harmonics; ++k) {
        double stretch = std::sqrt(1.0 + B * k * k);
        double harmonicFreqHz = Freq * k * stretch;
        if (harmonicFreqHz >= m_Sr / 2.0)
            break;

        double sigmaHz = harmonicFreqHz * sigmaConst;
        if (sigmaHz < binWidth * 0.75) {
            sigmaHz = binWidth * 0.75;
        }
        double envelope = std::exp(-beta * (k - 1));
        if (k > 1) {
            envelope *= 1.25;
        }
        if (envelope < 1e-5)
            break;

        double rangeHz = 4.0 * sigmaHz;
        int minBin = static_cast<int>(std::floor((harmonicFreqHz - rangeHz) / binWidth));
        int maxBin = static_cast<int>(std::ceil((harmonicFreqHz + rangeHz) / binWidth));

        minBin = std::max(0, minBin);
        maxBin = std::min(static_cast<int>(m_FFTSize / 2) - 1, maxBin);

        double twoSigmaSq = 2.0 * sigmaHz * sigmaHz;
        double normalizationFactor = 1.0 / (sigmaHz * std::sqrt(2.0 * std::numbers::pi));

        for (int i = minBin; i <= maxBin; ++i) {
            double binFreq = i * binWidth;
            double diff = binFreq - harmonicFreqHz;
            double exponent = -(diff * diff) / twoSigmaSq;
            double gaussian = normalizationFactor * std::exp(exponent);
            templateBins[i] += envelope * gaussian;
        }
    }

    double sum = std::accumulate(templateBins.begin(), templateBins.end(), 0.0);
    if (sum > 1e-12) {
        double invSum = 1.0 / sum;
        for (auto &val : templateBins) {
            val *= invSum;
        }
    }

    std::vector<PitchTemplateEntry> precomputedTemplate;
    precomputedTemplate.reserve(templateBins.size());
    for (size_t i = 0; i < templateBins.size(); ++i) {
        double P = templateBins[i];
        if (P > 0.0) {
            precomputedTemplate.push_back({i, P, P * std::log(P)});
        }
    }

    m_PitchTemplatesPrecomputed[rootBinFreq] = std::move(precomputedTemplate);
    m_PitchTemplates[rootBinFreq] = std::move(templateBins);
}

// ─────────────────────────────────────
// GONG 2015 (adapted)
/**
 * @brief Rebuild pitch templates for the installed score.
 *
 * @note Scans note observations and microstate observations after clearing both template caches.
 */
void OnlineForward::UpdatePitchTemplate() {
    int StateSize = (int)m_States.size();
    m_PitchTemplates.clear();
    m_PitchTemplatesPrecomputed.clear();

    for (int h = 0; h < StateSize; h++) {
        if (m_States[h].Type == NOTE || m_States[h].MicroTopologyType != NO_MICROSTATES) {
            for (const Observation &SubState : m_States[h].Observations) {
                if (SubState.Type == PITCH) {
                    BuildPitchTemplate(SubState.Freq);
                }
            }
            for (const MarkovMicroState &MicroState : m_States[h].MicroStates) {
                for (const Observation &SubState : MicroState.Observations) {
                    if (SubState.Type == PITCH) {
                        BuildPitchTemplate(SubState.Freq);
                    }
                }
            }
        }
    }
}

// ─────────────────────────────────────
// GONG 2015 (adapted)
/**
 * @brief Retrieve a cached template, building it if needed.
 *
 * @param Freq Pitch frequency in Hz.
 *
 * @return Copy of the spectral template, or an empty vector if no template was built.
 *
 * @note Templates sharing a rounded fundamental FFT bin share the same cache entry.
 */
PitchTemplateArray OnlineForward::GetPitchTemplate(double Freq) {
    BuildPitchTemplate(Freq);
    double rootBinFreq = std::round(Freq / (m_Sr / m_FFTSize));
    return m_PitchTemplates[rootBinFreq];
}

// ─────────────────────────────────────
// CONT (2010) section 3.1;
// CUVILLIER (2016) section 2.2.2;
// GONG (2015)
/**
 * @brief Evaluate spectral pitch evidence using a cached harmonic template.
 *
 * @param Freq Pitch frequency in Hz.
 *
 * @return Pitch evidence, or the smallest positive double when the pitch or available spectrum is unusable.
 *
 * @note Uses exponential KL-divergence scoring with reverberation and spectral-deviation terms.
 */
double OnlineForward::GetPitchProbability(double Freq) {
    constexpr double minProb = std::numeric_limits<double>::min();

    if (Freq <= 0.0 || m_FFTSize <= 0.0 || m_Sr <= 0.0) {
        return minProb;
    }

    const double binWidth = m_Sr / m_FFTSize;
    const double rootBinFreq = std::round(Freq / binWidth);

    auto it = m_PitchTemplates.find(rootBinFreq);
    if (it == m_PitchTemplates.end()) {
        BuildPitchTemplate(Freq);
        it = m_PitchTemplates.find(rootBinFreq);

        if (it == m_PitchTemplates.end()) {
            return minProb;
        }
    }

    const PitchTemplateArray &pitchTemplate = it->second;
    const auto &reverb = m_Desc.ReverbSpectralPower;
    const auto &spectrum = m_Desc.SpectralMagnitudeFrameNorm;

    const size_t bins = std::min({static_cast<size_t>(m_FFTSize / 2), pitchTemplate.size(), spectrum.size()});

    if (bins == 0) {
        return minProb;
    }

    const double *templateData = pitchTemplate.data();
    const double *spectrumData = spectrum.data();
    const double *reverbData = reverb.data();

    const size_t reverbBins = std::min(bins, reverb.size());

    double klDiv = 0.0;

    // Bins where reverb data exists.
    size_t i = 0;
    for (; i < reverbBins; ++i) {
        const double P = templateData[i] + reverbData[i];
        const double Q = spectrumData[i];

        if (P > 0.0 && Q > 0.0) {
            klDiv += P * std::log(P / Q);
        } else if (P == 0.0 && Q >= 0.0) {
            klDiv += Q;
        }
    }

    // Remaining bins: reverb = 0.
    for (; i < bins; ++i) {
        const double P = templateData[i];
        const double Q = spectrumData[i];

        if (P > 0.0 && Q > 0.0) {
            klDiv += P * std::log(P / Q);
        } else if (P == 0.0 && Q >= 0.0) {
            klDiv += Q;
        }
    }

    klDiv /= 1.0 + m_Desc.StdDev;
    return std::exp(-m_PitchScalingFactor * klDiv);
}

} // namespace OpenScofo
