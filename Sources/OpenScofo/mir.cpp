/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file mir.cpp
 * @brief Audio descriptor extraction, spectral analysis, and classifier integration.
 *
 * @note Computes core frame features and enables optional analysis stages according to requested descriptors.
 * @warning Analysis uses mutable FFT buffers and feature history; serialize processing and reconfiguration.
 */

#include "mir.hpp"
#include <algorithm>
#include <limits>
#include <utility>

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
    m_PreviousSpectralPower.assign(static_cast<size_t>(std::lround(m_Config.FFTSize / 2.0f)) + 1, 0.0);
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
        YINInit();
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
        YINInit();
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
    for (size_t i = 0; i + minTau < frame; ++i) {
        const double x = data[i];

        const size_t limit = std::min(maxTau, frame - i - 1);

        double *d = diff + minTau;
        const double *lag = data + i + minTau;
        const double *const lagEnd = data + i + limit + 1;

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

    for (size_t tau = minTau; tau <= maxTau; ++tau) {
        const double value = cmndf[tau];

        if (value < m_Config.YINThreshold) {
            tauEstimate = tau;
            bestValue = value;

            while (tau + 1 <= maxTau && cmndf[tau + 1] < bestValue) {
                ++tau;
                tauEstimate = tau;
                bestValue = cmndf[tau];
            }

            break;
        }
    }

    // No threshold crossing: fall back to the global minimum,
    // exactly as the previous implementation did.
    if (tauEstimate == 0) {
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
    double refinedTau = static_cast<double>(tauEstimate);

    if (tauEstimate > minTau && tauEstimate + 1 <= maxTau) {

        const double left = cmndf[tauEstimate - 1];
        const double center = cmndf[tauEstimate];
        const double right = cmndf[tauEstimate + 1];

        const double denominator = left - 2.0 * center + right;

        if (std::abs(denominator) > 1e-12) {
            const double offset = 0.5 * (left - right) / denominator;

            refinedTau += std::clamp(offset, -1.0, 1.0);
        }
    }

    const double confidence = std::clamp(1.0 - bestValue, 0.0, 1.0);

    if (refinedTau <= 0.0 || confidence <= 0.0) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    const double pitch = m_Config.SR / refinedTau;

    if (pitch < m_Config.YINMinFrequency || pitch > m_Config.YINMaxFrequency) {

        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    Desc.Pitch = pitch;
    Desc.PitchConfidence = confidence;
}

// ╭─────────────────────────────────────╮
// │              SPECTRAL               │
// ╰─────────────────────────────────────╯
/**
 * @brief Finalize scalar spectral features from accumulated moments.
 *
 * @param Desc Audio description to populate or update in place.
 * @param acc Spectral sums and moments accumulated over the current FFT bins.
 * @param NHalf Number of nonnegative-frequency bins, including DC and Nyquist.
 *
 * @note Updates centroid history, spread, shape, irregularity, crest, flatness, harmonicity, and high-frequency
 * ratio.
 * @warning NHalf must be positive and the accumulators must describe the current frame.
 */
void MIR::ComputeScalarFeatures(Description &Desc, const SpectralAccumulators &acc, size_t NHalf) {
    Desc.SpectralFlux = std::sqrt(Desc.SpectralFlux);

    const double SumPowerEps = acc.SumPower + 1e-12;
    const double invSum = 1.0 / SumPowerEps;

    // Centroid related
    Desc.SpectralCentroid = acc.WSumFreqs * invSum;
    Desc.CentroidVelocity = std::abs(Desc.SpectralCentroid - m_PrevCentroid);
    m_PrevCentroid = Desc.SpectralCentroid;

    // Spectral Spread for Librosa
    const double E1 = Desc.SpectralCentroid;
    const double E2 = acc.sumFreq2 * invSum;
    const double E3 = acc.sumFreq3 * invSum;
    const double E4 = acc.sumFreq4 * invSum;
    const double mu3 = E3 - 3.0 * E1 * E2 + 2.0 * E1 * E1 * E1;
    const double mu4 = E4 - 4.0 * E1 * E3 + 6.0 * E1 * E1 * E2 - 3.0 * E1 * E1 * E1 * E1;
    const double spectralSpreadHzVariance = std::max(0.0, E2 - E1 * E1);
    Desc.SpectralSpreadHz = acc.SumPower > 0.0 ? std::sqrt(spectralSpreadHzVariance) : 0.0;

    // Essentia DistributionShape spread: normalized second central moment over FFT bin indices.
    const double EIndex = acc.sumIndex * invSum;
    const double EIndex2 = acc.sumIndex2 * invSum;
    const double rawIndexVariance = std::max(0.0, EIndex2 - EIndex * EIndex);
    const double indexRange = static_cast<double>(std::max<size_t>(1, NHalf - 1));
    const double invIndexRange = 1.0 / indexRange;
    Desc.SpectralSpreadVariance = acc.SumPower > 0.0 ? (rawIndexVariance * invIndexRange * invIndexRange) : 0.0;

    // Skewness, Kurtosis
    const double spreadEps = Desc.SpectralSpreadHz + 1e-12;
    Desc.SpectralSkewness = mu3 / (spreadEps * spreadEps * spreadEps);
    Desc.SpectralKurtosis = mu4 / (spreadEps * spreadEps * spreadEps * spreadEps) - 3.0;

    // Spectral Irregularity
    Desc.SpectralIrregularityJensen =
        acc.irregularityDenominator > 0.0 ? (acc.irregularityJensenNumerator / acc.irregularityDenominator) : -1.0;
    Desc.SpectralIrregularityKrimphoff =
        acc.irregularityKrimphoffSum > 0.0 ? std::log10(acc.irregularityKrimphoffSum) : -1.0;
    Desc.SpectralIrregularity = Desc.SpectralIrregularityKrimphoff;

    // Spectral Crest
    if (acc.maxMagCrest == 0.0f) {
        Desc.SpectralCrest = 0.0;
    } else {
        const float meanMagCrest = acc.sumMagCrest / static_cast<float>(NHalf);
        Desc.SpectralCrest = meanMagCrest > 0.0f ? static_cast<double>(acc.maxMagCrest / meanMagCrest) : 0.0;
    }

    Desc.SpectralFlatness = std::exp(acc.logSumPower / NHalf) / (acc.linSumPower / NHalf);
    Desc.Harmonicity = acc.harmonicityPeak / (acc.harmonicitySum + 1e-12);
    Desc.HighFreqRatio = acc.highFreqEnergy * invSum;

    const double K = static_cast<double>(NHalf);
    const double denom = (K * acc.sumFreqSq - acc.sumFreq * acc.sumFreq) + 1e-12;
    const double slope = (K * acc.WSumFreqs - acc.sumFreq * acc.SumPower) / denom;
    Desc.SpectralSlope = slope * invSum;
}

// ─────────────────────────────────────
/**
 * @brief Transform the windowed input and compute spectral descriptors.
 *
 * @param Desc Audio description to populate or update in place.
 *
 * @note Produces magnitude, power, normalized spectra, and scalar features while updating inter-frame flux history.
 * @warning Requires initialized FFT resources and descriptor arrays sized to FFTSize / 2 + 1.
 */
void MIR::GetSpectralDescriptions(Description &Desc) {
    const int NHalf = m_Config.FFTSize / 2 + 1;
    const double binWidth = static_cast<double>(m_Config.SR) / static_cast<double>(m_Config.FFTSize);
    const double invN = 1.0 / static_cast<double>(m_Config.FFTSize);
    constexpr double amin = 1e-10;
    const int hfStart = NHalf / 4;

    pffft_transform_ordered(m_FullFFTSetup, m_FullFFTIn, m_FullFFTOut, m_FullFFTWork, PFFFT_FORWARD);

    Desc.MaxAmp = 0.0;
    Desc.SpectralFlux = 0.0;

    SpectralAccumulators acc;
    double prevNorm = 0.0;
    double prevPrevMag = 0.0;
    double prevMag = 0.0;

    // Process first iteration (i = 0) explicitly to avoid "if (i > 0)" in the hot loop
    {
        double re = static_cast<double>(m_FullFFTOut[0]);
        double im = 0.0;
        const double p = re * re + im * im;
        const double mag = std::sqrt(p);
        const double norm = mag * invN;

        Desc.Power[0] = p;
        Desc.Magnitude[0] = mag;
        Desc.SpectralMagnitudeNorm[0] = norm;
        Desc.MaxAmp = norm;

        acc.maxMag = mag;
        acc.sumMagCrest += static_cast<float>(mag);
        acc.maxMagCrest = static_cast<float>(mag);

        acc.SumPower += norm;
        acc.irregularityDenominator += norm * norm;

        const double v = std::max(amin, p);
        acc.logSumPower += std::log(v);
        acc.linSumPower += v;
        acc.spectralEnergySum += p;

        const double diff = mag - m_PreviousSpectralPower[0];
        Desc.SpectralFlux += diff * diff;
        m_PreviousSpectralPower[0] = mag;

        // Note: i=0 is never >= hfStart (assuming FFT >= 8), so no highFreqEnergy addition here
        prevNorm = norm;
        prevMag = mag;
    }

    // Main Accumulation Loop (i > 0)
    for (int i = 1; i < NHalf; ++i) {
        double re = 0.0;
        double im = 0.0;
        const size_t half = static_cast<size_t>(m_Config.FFTSize) / 2;
        const size_t bin = static_cast<size_t>(i);
        if (bin == 0) {
            re = static_cast<double>(m_FullFFTOut[0]);
            im = 0.0;
        } else if (bin == half) {
            re = static_cast<double>(m_FullFFTOut[1]);
            im = 0.0;
        } else {
            const size_t idx = 2 * bin;
            re = static_cast<double>(m_FullFFTOut[idx]);
            im = static_cast<double>(m_FullFFTOut[idx + 1]);
        }
        const double p = re * re + im * im;
        const double mag = std::sqrt(p);
        const double norm = mag * invN;

        Desc.Power[i] = p;
        Desc.Magnitude[i] = mag;
        Desc.SpectralMagnitudeNorm[i] = norm;
        Desc.MaxAmp = std::max(Desc.MaxAmp, norm);

        acc.maxMag = std::max(acc.maxMag, mag);
        const float magCrest = static_cast<float>(mag);
        acc.sumMagCrest += magCrest;
        acc.maxMagCrest = std::max(acc.maxMagCrest, magCrest);

        const double freq = i * binWidth;
        const double freq2 = freq * freq;
        const double index = static_cast<double>(i);

        acc.SumPower += norm;
        acc.sumFreq += freq;
        acc.sumFreqSq += freq2;
        acc.sumFreq2 += freq2 * norm;
        acc.sumFreq3 += freq2 * freq * norm;
        acc.sumFreq4 += freq2 * freq2 * norm;
        acc.WSumFreqs += freq * norm;
        acc.sumIndex += index * norm;
        acc.sumIndex2 += index * index * norm;
        acc.irregularityDenominator += norm * norm;

        const double v = std::max(amin, p);
        acc.logSumPower += std::log(v);
        acc.linSumPower += v;
        acc.spectralEnergySum += p;

        const double diff = mag - m_PreviousSpectralPower[i];
        Desc.SpectralFlux += diff * diff;
        m_PreviousSpectralPower[i] = mag;

        // Branchless execution: static_cast resolves to 1.0 or 0.0
        acc.highFreqEnergy += norm * static_cast<double>(i >= hfStart);

        const double d = prevNorm - norm;
        acc.irregularityJensenNumerator += d * d;
        if (i > 1) {
            const double localAvg = (prevPrevMag + prevMag + mag) / 3.0;
            acc.irregularityKrimphoffSum += std::abs(prevMag - localAvg);
        }

        acc.harmonicitySum += norm;
        acc.harmonicityPeak = std::max(acc.harmonicityPeak, norm);

        prevNorm = norm;
        prevPrevMag = prevMag;
        prevMag = mag;
    }

    // SCALAR CALCULATIONS
    ComputeScalarFeatures(Desc, acc, NHalf);

    // PASS 2: NORMALIZE, ENTROPY, ROLLOFF & PREFIX FUSION
    const double sumPowerEps = acc.SumPower + 1e-12;
    const double invSum = 1.0 / sumPowerEps;
    const double Mean = 1.0 / static_cast<double>(NHalf);
    const double invSpectralEnergy = acc.spectralEnergySum > 0.0 ? (1.0 / acc.spectralEnergySum) : 0.0;
    const double rolloffCutoffEnergy = std::clamp(m_Config.SpectralRolloffCutoff, 0.0, 1.0) * acc.linSumPower;

    double Variance = 0.0;
    double cumulativeEnergy = 0.0;
    size_t rolloffBin = NHalf - 1;
    bool rolloffFound = false;

    Desc.SpectralEntropy = 0.0;
    m_SpectralPrefix[0] = 0.0;

    for (int i = 0; i < NHalf; ++i) {
        const double normSp = (Desc.SpectralMagnitudeNorm[i] + 1e-12) * invSum;
        Desc.SpectralMagnitudeFrameNorm[i] = normSp;
        Desc.SpectralPowerFrameNorm[i] = Desc.Power[i] * invSpectralEnergy;
        const double diffMean = normSp - Mean;
        Variance += diffMean * diffMean;
        const double currentPower = Desc.Power[i];

        if (invSpectralEnergy > 0.0) {
            const double prob = currentPower * invSpectralEnergy;
            if (prob > 0.0) {
                Desc.SpectralEntropy -= prob * std::log2(prob);
            }
        }
        cumulativeEnergy += currentPower;
        m_SpectralPrefix[i + 1] = cumulativeEnergy;
        if (!rolloffFound && cumulativeEnergy >= rolloffCutoffEnergy) {
            rolloffBin = i;
            rolloffFound = true;
        }
    }

    Desc.StdDev = std::sqrt(Variance * Mean);
    const double binToHz = (m_Config.SR * 0.5) / std::max(1, NHalf - 1);
    Desc.SpectralRolloff = static_cast<double>(rolloffBin) * binToHz;
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

    if (m_NeedYIN) {
        YINExec(In, Desc);
    } else {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
    }

    // 2. Frequency Domain (windowing + FFT)
    const double *x = In.data();
    const double *w = m_FullWindowingFunc.data();
    for (size_t i = 0; i < m_Config.FFTSize; ++i) {
        m_FullFFTIn[i] = static_cast<float>(x[i] * w[i]);
    }

    GetSpectralDescriptions(Desc);

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
