/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file classification.cpp
 * @brief ONNX integration and extended-technique classification.
 */

#include "mir.hpp"
#include <cmath>
#include <utility>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │          Machine Learning           │
// ╰─────────────────────────────────────╯
/**
 * @brief Load a classifier and initialize its audio feature dependencies.
 *
 * @param path Path to an ONNX classifier file.
 * @param Descriptors Ordered model input descriptors.
 *
 * @note Returns early when loading fails; successful loading refreshes flags and required optional analysis
 * resources.
 * @warning Model loading and feature initialization can allocate memory; serialize with processing.
 */
void MIR::ONNXInit(fs::path path, std::vector<Descriptors> Descriptors) {
    if (!m_ONNXModel.Load(path, std::move(Descriptors), m_Config)) {
        return;
    }

    UpdateDescriptorFlags();
    if (m_NeedOnset) {
        OnsetInit();
    }
    if (m_NeedMFCC) {
        MFCCInit();
    }
    if (m_NeedChroma) {
        SpectralChromaInit();
    }
    if (m_NeedZCR) {
        ZeroCrossingRateInit();
    }
    if (m_NeedYIN) {
        SigmundInit();
    }
}

// ─────────────────────────────────────
/**
 * @brief Copy the class labels of the loaded ONNX model.
 *
 * @return Copy of the classifier labels.
 *
 * @note Delegates to the model label list; an unloaded model has no labels.
 */
std::vector<std::string> MIR::GetONNXLabels() {
    return m_ONNXModel.GetLabels();
}

// ╭─────────────────────────────────────╮
// │        Percussive Technique         │
// ╰─────────────────────────────────────╯
/**
 * @brief Estimate extended-technique evidence from flux and pitch confidence.
 *
 * @param Desc Audio description to populate or update in place.
 *
 * @note Applies a sigmoid to spectral flux weighted by the complement of pitch confidence.
 * @warning Compute spectral flux and YIN pitch confidence before calling this stage.
 */
void MIR::ExtendedTechExec(Description &Desc) {
    // Desc.ExtendedTechProb = (1.0f - Desc.Harmonicity);
    Desc.ExtendedTechProb = Desc.SpectralFlux;
    // Harmonic, confidently pitched frames should not be classified as an
    // extended technique. This is the confidence term used by the original
    // detector; ZCR is already represented by the spectral/noise features.
    Desc.ExtendedTechProb *= (1.0f - Desc.PitchConfidence);
    // Desc.ExtendedTechProb *= abs(m_ODS->odfvalpost);
    float steepness = 5.0f;
    Desc.ExtendedTechProb = 1.0f / (1.0f + std::exp(-steepness * (Desc.ExtendedTechProb - 0.5f)));
}

} // namespace OpenScofo
