/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file mir.cpp
 * @brief Audio descriptor orchestration, configuration, and resource management.
 *
 * @note Computes core frame features and enables optional analysis stages according to requested descriptors.
 * @warning Analysis uses mutable FFT buffers and feature history; serialize processing and reconfiguration.
 */

#include "mir.hpp"
#include <algorithm>
#include <cmath>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │Constructor and Destructor Functions │
// ╰─────────────────────────────────────╯
/**
 * @brief Release FFT and onset-detector resources.
 *
 * @note Destroys the FFT setup and frees its aligned buffers and onset storage.
 * @warning Do not destroy the extractor while another thread is using its resources.
 */
MIR::~MIR() {
    if (m_FullFFTSetup != nullptr) {
        pffft_destroy_setup(m_FullFFTSetup);
        m_FullFFTSetup = nullptr;
    }

    if (m_FullFFTIn != nullptr) {
        pffft_aligned_free(m_FullFFTIn);
        m_FullFFTIn = nullptr;
    }
    if (m_FullFFTOut != nullptr) {
        pffft_aligned_free(m_FullFFTOut);
        m_FullFFTOut = nullptr;
    }
    if (m_FullFFTWork != nullptr) {
        pffft_aligned_free(m_FullFFTWork);
        m_FullFFTWork = nullptr;
    }

    if (m_ODSData) {
        delete[] m_ODSData;
    }
    if (m_ODS) {
        delete m_ODS;
    }

    /// save values
}

// ─────────────────────────────────────
/**
 * @brief Rebuild audio analysis resources for a configuration.
 *
 * @param Config Configuration containing audio dimensions and requested descriptor settings.
 *
 * @note Resets spectral history and initializes requested optional descriptors and loudness filters.
 * @warning FFT sizes below 512 are rejected. Use a PFFFT-supported size and positive sample rate; serialize
 * processing.
 */
void MIR::UpdateConfiguration(const Configuration &Config) {
    m_Config = Config;
    UpdateDescriptorFlags();
    m_PrevCentroid = 0.0;
    m_PreviousSpectralMagnitude.assign(static_cast<size_t>(std::lround(m_Config.FFTSize / 2.0f)) + 1, 0.0);
    m_SpectralPrefix.resize(m_Config.FFTSize / 2 + 2);

    if (m_FullFFTSetup != nullptr) {
        pffft_destroy_setup(m_FullFFTSetup);
        m_FullFFTSetup = nullptr;
    }
    if (m_FullFFTIn != nullptr) {
        pffft_aligned_free(m_FullFFTIn);
        m_FullFFTIn = nullptr;
    }
    if (m_FullFFTOut != nullptr) {
        pffft_aligned_free(m_FullFFTOut);
        m_FullFFTOut = nullptr;
    }
    if (m_FullFFTWork != nullptr) {
        pffft_aligned_free(m_FullFFTWork);
        m_FullFFTWork = nullptr;
    }

    if (Config.FFTSize < 512) {
        spdlog::critical("OpenScofo requires FFTSize higher then 256");
        return;
    }

    FFTInit();

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
    InitITURFilters();

    spdlog::debug("Init MIR audio parameters using SR {}, FFTSize {}, HopSize {}", Config.SR, Config.FFTSize,
                  Config.HOPSize);
}

// ─────────────────────────────────────
/**
 * @brief Check whether a descriptor was explicitly requested.
 *
 * @param Descriptor Descriptor enum to search for in the request list.
 *
 * @return True if the descriptor is present in the request list; false otherwise.
 *
 * @note Searches the configured request list without expanding feature dependencies.
 */
bool MIR::DescriptorRequested(Descriptors Descriptor) const {
    return std::find(m_Config.RequestedDescriptors.begin(), m_Config.RequestedDescriptors.end(), Descriptor) !=
           m_Config.RequestedDescriptors.end();
}

// ─────────────────────────────────────
/**
 * @brief Resolve which optional analysis stages are required.
 *
 * @note Includes dependencies of extended-technique analysis and the loaded ONNX model when ONNX is requested.
 */
void MIR::UpdateDescriptorFlags() {
    m_NeedYIN =
        DescriptorRequested(YIN) || DescriptorRequested(YINCONFIDENCE) || DescriptorRequested(EXTENDEDTECHNIQUE);
    m_NeedMFCC = DescriptorRequested(MFCC) || DescriptorRequested(LOGMEL);
    m_NeedChroma = DescriptorRequested(CHROMA);
    m_NeedZCR = DescriptorRequested(ZCR) || DescriptorRequested(EXTENDEDTECHNIQUE);
    m_NeedExtendedTech = DescriptorRequested(EXTENDEDTECHNIQUE);
    m_NeedOnset = DescriptorRequested(ODSONSET) || m_NeedExtendedTech;
    m_NeedONNX = DescriptorRequested(ONNX);

    if (!m_NeedONNX) {
        return;
    }

    for (const Descriptors d : m_ONNXModel.GetDescriptors()) {
        switch (d) {
        case MFCC:
        case LOGMEL:
            m_NeedMFCC = true;
            break;
        case CHROMA:
            m_NeedChroma = true;
            break;
        case YIN:
        case YINCONFIDENCE:
            m_NeedYIN = true;
            break;
        case ZCR:
            m_NeedZCR = true;
            break;
        case ODSONSET:
            m_NeedOnset = true;
            break;
        case EXTENDEDTECHNIQUE:
            m_NeedYIN = true;
            m_NeedZCR = true;
            m_NeedOnset = true;
            m_NeedExtendedTech = true;
            break;
        default:
            break;
        }
    }
}

// ╭─────────────────────────────────────╮
// │            Main Function            │
// ╰─────────────────────────────────────╯
/**
 * @brief Analyze an audio frame and update its requested descriptors.
 *
 * @param In Input audio frame containing the configured FFT window samples.
 * @param Desc Audio description to populate or update in place.
 *
 * @note Always computes signal power and spectral features, then runs enabled optional stages and ONNX inference.
 * @warning Requires initialized processing resources, a full FFT window, and descriptor arrays matching
 * configuration.
 */
void MIR::GetDescription(const std::vector<double> &In, Description &Desc) {
    // 1. Temporal Domain
    GetSignalPower(In, Desc);
    if (m_NeedZCR) {
        ZeroCrossingRateExec(In, Desc);
    } else {
        Desc.ZeroCrossingRate = 0.0;
    }

    // 2. Frequency Domain (windowing + FFT)
    const double *x = In.data();
    const double *w = m_FullWindowingFunc.data();
    for (size_t i = 0; i < m_Config.FFTSize; ++i) {
        m_FullFFTIn[i] = static_cast<float>(x[i] * w[i]);
    }

    GetSpectralDescriptions(Desc);

    if (m_NeedYIN) {
        SigmundExec(Desc);
    } else {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
    }

    if (m_NeedOnset) {
        OnsetExec(Desc);
    } else {
        Desc.Onset = 0.0;
    }

    if (m_NeedExtendedTech) {
        ExtendedTechExec(Desc);
    } else {
        Desc.ExtendedTechProb = 0.0;
    }

    if (m_NeedMFCC) {
        MFCCExec(Desc);
    }

    if (m_NeedChroma) {
        SpectralChromaExec(Desc);
    }

    if (m_ONNXModel.IsLoaded() && m_NeedONNX) {
        m_ONNXModel.Execute(Desc);
    } else {
        Desc.ONNX.clear();
    }
}

} // namespace OpenScofo
