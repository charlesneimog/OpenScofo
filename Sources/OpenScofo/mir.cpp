/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "mir.hpp"
#include <algorithm>
#include <limits>
#include <utility>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │Constructor and Destructor Functions │
// ╰─────────────────────────────────────╯
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

    if (m_ONNXModelLoaded) {
        if (m_ONNXContext) {
            onnx_context_free(m_ONNXContext);
            m_ONNXContext = nullptr;
        }
        m_ONNXModelLoaded = false;
    }

    /// save values
}

// ─────────────────────────────────────
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
bool MIR::DescriptorRequested(Descriptors Descriptor) const {
    return std::find(m_Config.RequestedDescriptors.begin(), m_Config.RequestedDescriptors.end(), Descriptor) !=
           m_Config.RequestedDescriptors.end();
}

// ─────────────────────────────────────
void MIR::UpdateDescriptorFlags() {
    m_NeedYIN = DescriptorRequested(YIN) || DescriptorRequested(YINCONFIDENCE);
    m_NeedMFCC = DescriptorRequested(MFCC) || DescriptorRequested(LOGMEL);
    m_NeedChroma = DescriptorRequested(CHROMA);
    m_NeedZCR = DescriptorRequested(ZCR) || DescriptorRequested(EXTENDEDTECHNIQUE);
    m_NeedExtendedTech = DescriptorRequested(EXTENDEDTECHNIQUE);
    m_NeedOnset = DescriptorRequested(ODSONSET) || m_NeedExtendedTech;
    m_NeedONNX = DescriptorRequested(ONNX);

    if (!m_NeedONNX) {
        return;
    }

    for (const Descriptors d : m_ONNXDescriptors) {
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
void MIR::ONNXInit(fs::path path, std::vector<Descriptors> Descriptors) {
    auto u8 = path.u8string();
    std::string path_utf8(u8.begin(), u8.end());

    if (m_ONNXModelLoaded && m_ONNXModelPath == path) {
        spdlog::debug("ONNX model already loaded, reusing existing context");
        return;
    }
    m_ONNXModelPath = path;

    if (m_ONNXModelLoaded) {
        spdlog::debug("Loading different ONNX model, cleaning up old one");
        if (m_ONNXContext) {
            onnx_context_free(m_ONNXContext);
            m_ONNXContext = nullptr;
        }
        m_ONNXModelLoaded = false;
        m_ONNXLabels.clear();
        m_Writers.clear();
    }

    m_ONNXContext = onnx_context_alloc_from_file(path_utf8.c_str(), nullptr, 0);
    if (m_ONNXContext == nullptr) {
        spdlog::error("Failed to load ONNX model: {}.", path.string());
        return;
    }

    if (m_ONNXContext && m_ONNXContext->g != nullptr) {
        struct onnx_graph_t *g = m_ONNXContext->g;
        for (int i = 0; i < g->nlen; i++) {
            struct onnx_node_t *n = &g->nodes[i];
            if (n->opset > CURRENT_ONNX_OPSET) {
                spdlog::error("Unsupported opset => {} {}.", n->proto->op_type, n->opset);
                return;
            }
        }
    }

    // Get labels
    m_ONNXLabels.clear();
    bool TreeEnsembleClassifierFound = false;
    struct onnx_graph_t *g = m_ONNXContext->g;
    for (int i = 0; i < g->nlen; i++) {
        struct onnx_node_t *n = &g->nodes[i];
        if (!n || strcmp(n->proto->op_type, "TreeEnsembleClassifier") != 0)
            continue;

        TreeEnsembleClassifierFound = true;
        for (size_t k = 0; k < n->proto->n_attribute; k++) {
            Onnx__AttributeProto *attr = n->proto->attribute[k];
            if (strcmp(attr->name, "classlabels_strings") == 0) {
                for (size_t v = 0; v < attr->n_strings; v++) {
                    ProtobufCBinaryData *str = &attr->strings[v];
                    std::string label(reinterpret_cast<char *>(str->data), str->len);
                    m_ONNXLabels.push_back(label);
                }
            }
        }
    }

    if (!TreeEnsembleClassifierFound) {
        spdlog::error(
            "TreeEnsembleClassifier not found in model, please use PureData py.train to train models for OpenScofo");
        return;
    }

    m_ONNXModelLoaded = true;
    m_ONNXDescriptors = Descriptors;
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
    spdlog::debug("ONNX Model Loaded");
    // Count
    m_ONNXDescriptorsSize = 0;
    for (size_t i = 0; i < m_ONNXDescriptors.size(); i++) {
        switch (m_ONNXDescriptors[i]) {
        case MFCC:
            m_ONNXDescriptorsSize += m_Config.MFCCCount;
            break;
        case CHROMA:
            m_ONNXDescriptorsSize += m_Config.ChromaSize;
            break;
        case POWERARRAY:
            m_ONNXDescriptorsSize += (m_Config.FFTSize / 2) + 1;
            break;
        case MAGNITUDE:
            m_ONNXDescriptorsSize += (m_Config.FFTSize / 2) + 1;
            break;
        case LOGMEL:
            m_ONNXDescriptorsSize += m_Config.MFCCMels;
            break;
        default:
            m_ONNXDescriptorsSize++;
        }

        spdlog::debug("Descriptors ID: {}, {}", (int)m_ONNXDescriptors[i], m_ONNXDescriptorsSize);
    }

    spdlog::debug("ONNX Model array is {}", m_ONNXDescriptorsSize);
    m_ONNXDescriptorsArray.resize(m_ONNXDescriptorsSize);

    for (auto d : m_ONNXDescriptors) {
        switch (d) {
        case MFCC:
            m_Writers.push_back([this](const Description &desc, float *&out) {
                for (int i = 0; i < m_Config.MFCCCount; ++i)
                    *out++ = desc.MFCC[i];
            });
            break;
        case CHROMA:
            m_Writers.push_back([this](const Description &desc, float *&out) {
                for (int i = 0; i < m_Config.ChromaSize; ++i)
                    *out++ = desc.Chroma[i];
            });
            break;
        case POWERARRAY:
            m_Writers.push_back([](const Description &desc, float *&out) {
                for (int i = 0; i < (int)desc.Power.size(); ++i)
                    *out++ = desc.Power[i];
            });
            break;
        case LOGMEL:
            m_Writers.push_back([this](const Description &desc, float *&out) {
                for (int i = 0; i < m_Config.MFCCMels; ++i)
                    *out++ = desc.LogMelSpectrum[i];
            });
            break;
        case MAGNITUDE:
            m_Writers.push_back([](const Description &desc, float *&out) {
                for (double v : desc.Magnitude)
                    *out++ = v;
            });
            break;
        case LOUDNESS:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.Loudness; });
            break;
        case RMS:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.RMS; });
            break;
        case ZCR:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.ZeroCrossingRate; });
            break;
        case HFR:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.HighFreqRatio; });
            break;
        case CENTROID:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralCentroid; });
            break;
        case SPREADHZ:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralSpreadHz; });
            break;
        case FLATNESS:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralFlatness; });
            break;
        case ENTROPY:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralEntropy; });
            break;
        case ROLLOFF:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralRolloff; });
            break;
        case FLUX:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralFlux; });
            break;
        case IRREGULARITY:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralIrregularity; });
            break;
        case HARMONICITY:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.Harmonicity; });
            break;
        case YIN:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.Pitch; });
            break;
        case SILENCEPROB:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SilenceProb; });
            break;
        case EXTENDEDTECHNIQUE:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.ExtendedTechProb; });
            break;
        case ODSONSET:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.Onset ? 1.0f : 0.0f; });
            break;
        case STDDEV:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.StdDev; });
            break;
        case SKEWNESS:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralSkewness; });
            break;
        case SLOPE:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralSlope; });
            break;
        case DB:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.dB; });
            break;
        case MAXAMP:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.MaxAmp; });
            break;
        case SPREADVARIANCE:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralSpreadVariance; });
            break;
        case CREST:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralCrest; });
            break;
        case KURTOSIS:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralKurtosis; });
            break;
        case CENTROIDVEL:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.CentroidVelocity; });
            break;
        case YINCONFIDENCE:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.PitchConfidence; });
            break;
        case ONNX:
            break;
        case INVALID:
            spdlog::error("Invalid descriptor");
            break;
        }
    }

    m_InputTensor = onnx_tensor_search(m_ONNXContext, "features");
    if (m_InputTensor == nullptr) {
        spdlog::error("Tensor 'features', for input, not found");
        return;
    }

    // probabilities
    m_OutputTensor = onnx_tensor_search(m_ONNXContext, "probabilities");
    if (m_OutputTensor == nullptr) {
        spdlog::error("Tensor 'probabilities', for output, not found");
        return;
    }

    if ((size_t)m_ONNXDescriptorsSize != m_InputTensor->ndata) {
        spdlog::error("Tensor input expect {} values, selected descriptors give {}", m_InputTensor->ndata,
                      m_ONNXDescriptorsSize);
        return;
    }

    size_t LabelSize = m_ONNXLabels.size();
    if (LabelSize != m_OutputTensor->ndata) {
        spdlog::error("Tensor output expect {} values, but labels provided are {}", m_InputTensor->ndata,
                      m_ONNXDescriptorsSize);
        return;
    }
}

// ─────────────────────────────────────
void MIR::ONNXExec(Description &Desc) {

    float *out = m_ONNXDescriptorsArray.data();
    for (auto &w : m_Writers) {
        w(Desc, out);
    }

    onnx_tensor_apply(m_InputTensor, m_ONNXDescriptorsArray.data(), m_ONNXDescriptorsSize);
    onnx_run(m_ONNXContext);

    switch (m_OutputTensor->type) {
    case ONNX_TENSOR_TYPE_FLOAT32: {
        float *data = (float *)m_OutputTensor->datas;
        for (size_t i = 0; i < m_ONNXLabels.size(); i++) {
            Desc.ONNX[m_ONNXLabels[i]] = data[i];
        }
        break;
    }
    case ONNX_TENSOR_TYPE_FLOAT64: {
        double *data = (double *)m_OutputTensor->datas;
        for (size_t i = 0; i < m_ONNXLabels.size(); i++) {
            Desc.ONNX[m_ONNXLabels[i]] = (float)data[i];
        }
        break;
    }
    default:
        spdlog::error("Tensor output type not supported {}", (int)m_OutputTensor->type);
    }
}

// ─────────────────────────────────────
std::vector<std::string> MIR::GetONNXLabels() {
    return m_ONNXLabels;
}

// ╭─────────────────────────────────────╮
// │           Onset Detector            │
// ╰─────────────────────────────────────╯
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
inline void ReadPffftBin(const float *out, size_t fftSize, size_t bin, double &re, double &im) {
    const size_t half = fftSize / 2;
    if (bin == 0) {
        re = static_cast<double>(out[0]);
        im = 0.0;
        return;
    }
    if (bin == half) {
        re = static_cast<double>(out[1]);
        im = 0.0;
        return;
    }
    const size_t idx = 2 * bin;
    re = static_cast<double>(out[idx]);
    im = static_cast<double>(out[idx + 1]);
}

// ─────────────────────────────────────
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
void MIR::ExtendedTechExec(Description &Desc) {
    Desc.ExtendedTechProb = (1.0f - Desc.Harmonicity);
    Desc.ExtendedTechProb *= Desc.SpectralFlux;
    Desc.ExtendedTechProb *= (1.0f - Desc.ZeroCrossingRate);
    Desc.ExtendedTechProb *= abs(m_ODS->odfvalpost);
    float steepness = 5.0f;
    Desc.ExtendedTechProb = 1.0f / (1.0f + std::exp(-steepness * (Desc.ExtendedTechProb - 0.5f)));
}

// ╭─────────────────────────────────────╮
// │         Power and Amplitude         │
// ╰─────────────────────────────────────╯
// Check https://github.com/klangfreund/LUFSMeter (use MIT)
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
void MIR::YINInit() {
    const size_t frameSize = static_cast<size_t>(std::max(2.0f, m_Config.FFTSize));
    const size_t half = frameSize / 2;
    const size_t allocSize = half + 2;
    m_YINDifference.assign(allocSize, 0.0);
    m_YINCMNDF.assign(allocSize, 1.0);
}

// ─────────────────────────────────────
void MIR::YINExec(const std::vector<double> &In, Description &Desc) {
    const size_t frame = In.size();
    const size_t minTau = m_Config.SR / m_Config.YINMaxFrequency;
    const size_t maxTauByPitch = std::ceil(m_Config.SR / m_Config.YINMinFrequency);
    const size_t maxTau = std::min({frame / 2, m_YINDifference.size() - 1, maxTauByPitch});
    if (maxTau <= minTau) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    double *Diff = m_YINDifference.data();
    double *Cmnfg = m_YINCMNDF.data();
    std::fill_n(Diff, maxTau + 1, 0.0);
    std::fill_n(Cmnfg, maxTau + 1, 1.0);

    const double *data = In.data();
    const size_t frameMinus1 = frame - 1;

    // Walk the input once and update a contiguous lag window for each sample.
    for (size_t i = 0; i < frameMinus1; ++i) {
        const double x_i = data[i];
        const size_t remaining = frameMinus1 - i;
        const size_t limit = std::min(maxTau, remaining);
        const double *lagPtr = data + i + 1;
        for (size_t tau = minTau; tau <= limit; ++tau) {
            const double delta = x_i - lagPtr[tau - 1];
            Diff[tau] += delta * delta;
        }
    }

    double cumulative = 0.0;
    Cmnfg[0] = 1.0;
    for (size_t tau = 1; tau <= maxTau; ++tau) {
        cumulative += Diff[tau];
        if (cumulative <= 0.0) {
            Cmnfg[tau] = 1.0;
        } else {
            Cmnfg[tau] = (Diff[tau] * static_cast<double>(tau)) / cumulative;
        }
    }

    size_t tauEstimate = 0;
    double bestValue = std::numeric_limits<double>::infinity();
    for (size_t tau = minTau; tau <= maxTau; ++tau) {
        const double value = Cmnfg[tau];
        if (value < bestValue) {
            bestValue = value;
            tauEstimate = tau;
        }
    }

    for (size_t tau = minTau; tau <= maxTau; ++tau) {
        const double value = Cmnfg[tau];
        if (value < m_Config.YINThreshold) {
            double prevValue = value;
            tauEstimate = tau;
            while (tau + 1 <= maxTau && Cmnfg[tau + 1] < prevValue) {
                prevValue = Cmnfg[++tau];
                tauEstimate = tau;
            }
            bestValue = prevValue;
            break;
        }
    }

    if (tauEstimate == 0 || !std::isfinite(bestValue)) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    double refinedTau = static_cast<double>(tauEstimate);
    if (tauEstimate > minTau && tauEstimate + 1 <= maxTau) {
        const double left = Cmnfg[tauEstimate - 1];
        const double center = Cmnfg[tauEstimate];
        const double right = Cmnfg[tauEstimate + 1];
        const double denominator = left - (2.0 * center) + right;
        if (std::abs(denominator) > 1e-12) {
            const double offset = 0.5 * (left - right) / denominator;
            refinedTau += std::clamp(offset, -1.0, 1.0);
        }
    }

    const double Confidence = std::clamp(1.0 - bestValue, 0.0, 1.0);
    if (refinedTau <= 0.0 || Confidence <= 0.0) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    const double Pitch = m_Config.SR / refinedTau;
    if (Pitch < m_Config.YINMinFrequency || Pitch > m_Config.YINMaxFrequency) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    Desc.Pitch = Pitch;
    Desc.PitchConfidence = Confidence;
}

// ╭─────────────────────────────────────╮
// │              SPECTRAL               │
// ╰─────────────────────────────────────╯
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
    const double invSum = 1.0 / acc.SumPower;
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
        const double normSp = (Desc.SpectralMagnitudeNorm[i]) * invSum;
        Desc.SpectralMagnitudeFrameNorm[i] = normSp;
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
double MIR::HzToOcts(double frequency, double tuning, int binsPerOctave) const {
    const double a440 = m_Config.TunningA4 * std::pow(2.0, tuning / static_cast<double>(binsPerOctave));
    return std::log2(frequency / (a440 / 16.0));
}

// ─────────────────────────────────────
double MIR::PositiveRemainder(double value, double modulus) const {
    double result = std::fmod(value, modulus);
    if (result < 0.0) {
        result += modulus;
    }
    return result;
}

// ─────────────────────────────────────
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
void MIR::ZeroCrossingRateInit() {
    const size_t frame = static_cast<size_t>(std::max(1.0f, m_Config.FFTSize));
    const size_t pad = m_Config.ZCRCenter ? (frame / 2) : 0;
    m_ZCRScratch.resize(frame + (2 * pad));
}

// ─────────────────────────────────────
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

    if (m_ONNXModelLoaded && m_NeedONNX) {
        ONNXExec(Desc);
    } else {
        Desc.ONNX.clear();
    }
}

} // namespace OpenScofo
