#include "mir.hpp"
#include <algorithm>
#include <limits>
#include <utility>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │Constructor and Destructor Functions │
// ╰─────────────────────────────────────╯
MIR::MIR(float Sr, float FftSize, float HopSize) {
    UpdateAudioParameters(Sr, FftSize, HopSize);
}

// ─────────────────────────────────────
MIR::~MIR() {
    if (m_FullFFTPlan != nullptr) {
        fftw_destroy_plan(m_FullFFTPlan);
        m_FullFFTPlan = nullptr;
    }

    if (m_FullFFTIn != nullptr) {
        fftw_free(m_FullFFTIn);
        m_FullFFTIn = nullptr;
    }
    if (m_FullFFTOut != nullptr) {
        fftw_free(m_FullFFTOut);
        m_FullFFTOut = nullptr;
    }

    if (m_OnsetInit) {
        delete[] m_ODSData;
        delete m_ODS;
    }

    if (m_ONNXModelLoaded) {
        if (m_ONNXContext) {
            onnx_context_free(m_ONNXContext);
            m_ONNXContext = nullptr;
        }
        m_ONNXModelLoaded = false;
    }
}

// ╭─────────────────────────────────────╮
// │                Utils                │
// ╰─────────────────────────────────────╯
double MIR::Mtof(double note, double tunning) {
    return tunning * std::pow(2.0, (note - 69) / 12.0);
}

// ─────────────────────────────────────
double MIR::Ftom(double freq, double tunning) {
    return 69 + 12 * log2(freq / tunning);
}

// ─────────────────────────────────────
double MIR::Freq2Bin(double Freq, double n, double sr) {
    double bin;
    bin = Freq * n / sr;
    return round(bin);
}

// ─────────────────────────────────────
void MIR::UpdateAudioParameters(float Sr, float FftSize, float HopSize) {
    const bool sameAudioConfig = (m_HopSize == HopSize && m_FFTSize == FftSize && m_Sr == Sr);
    if (sameAudioConfig && m_FullFFTPlan != nullptr && m_FullFFTIn != nullptr && m_FullFFTOut != nullptr &&
        m_FullWindowingFunc.size() == static_cast<size_t>(FftSize)) {
        return;
    }

    m_HopSize = HopSize;
    m_FFTSize = FftSize;
    m_Sr = Sr;
    m_OnsetFFTSize = std::max(2, static_cast<int>(std::lround(FftSize)));
    m_BlockSize = HopSize;
    m_Accum = 0;
    m_PrevCentroid = 0.0;

    m_PreviousSpectralPower.assign(static_cast<size_t>(std::lround(FftSize / 2.0f)) + 1, 0.0);

    if (m_FullFFTPlan != nullptr) {
        fftw_destroy_plan(m_FullFFTPlan);
        m_FullFFTPlan = nullptr;
    }
    if (m_FullFFTIn != nullptr) {
        fftw_free(m_FullFFTIn);
        m_FullFFTIn = nullptr;
    }
    if (m_FullFFTOut != nullptr) {
        fftw_free(m_FullFFTOut);
        m_FullFFTOut = nullptr;
    }

    m_PreviousSpectralPower.assign(static_cast<size_t>(std::lround(m_FFTSize / 2.0f)) + 1, 0.0);

    FFTWInit();
    OnsetInit();
    MFCCInit();
    SpectralChromaInit();
    ZeroCrossingRateInit();
    YINInit();
    InitITURFilters();

    spdlog::debug("Init MIR audio parameters using SR {}, FFTSize {}, HopSize {}", Sr, FftSize, HopSize);
}

// ╭─────────────────────────────────────╮
// │          Set|Get Functions          │
// ╰─────────────────────────────────────╯
void MIR::SetdBTreshold(double dB) {
    m_dBTreshold = dB;
}

// ─────────────────────────────────────
void MIR::FFTWInit() {
    int WindowHalf = round(m_FFTSize / 2);
    m_FullFFTIn = fftw_alloc_real(m_FFTSize);
    if (!m_FullFFTIn) {
        spdlog::critical("fftw_alloc_real failed");
        return;
    }

    m_FullFFTOut = fftw_alloc_complex(WindowHalf + 1);
    if (!m_FullFFTOut) {
        fftw_free(m_FullFFTIn);
        spdlog::critical("fftw_alloc_complex failed");
        return;
    }

    m_FullFFTPlan = fftw_plan_dft_r2c_1d((int)m_FFTSize, m_FullFFTIn, m_FullFFTOut, FFTW_ESTIMATE);

    // Match librosa/scipy get_window('hann', N, fftbins=True): periodic Hann.
    m_FullWindowingFunc.resize(m_FFTSize);
    for (size_t i = 0; i < m_FFTSize; i++) {
        m_FullWindowingFunc[i] = 0.5 * (1.0 - cos(2.0 * std::numbers::pi * i / m_FFTSize));
    }
}

// ╭─────────────────────────────────────╮
// │          Machine Learning           │
// ╰─────────────────────────────────────╯
void MIR::ONNXInit(fs::path path, std::vector<Descriptors> Descriptors) {
    auto u8 = path.u8string();
    std::string path_utf8(u8.begin(), u8.end());
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
    spdlog::info("ONNX Model Loaded");
    // Count
    m_ONNXDescriptorsSize = 0;
    for (size_t i = 0; i < m_ONNXDescriptors.size(); i++) {
        switch (m_ONNXDescriptors[i]) {
        case MFCC:
            m_ONNXDescriptorsSize += m_MFCCCount;
            break;
        case CHROMA:
            m_ONNXDescriptorsSize += m_ChromaSize;
            break;
        case MELOGRAM:
            m_ONNXDescriptorsSize += m_MFCCMels;
            break;
        default:
            m_ONNXDescriptorsSize++;
        }

        spdlog::info("Descriptors ID: {}, {}", (int)m_ONNXDescriptors[i], m_ONNXDescriptorsSize);
    }

    spdlog::info("ONNX Model array is {}", m_ONNXDescriptorsSize);
    m_ONNXDescriptorsArray.resize(m_ONNXDescriptorsSize);

    for (auto d : m_ONNXDescriptors) {
        switch (d) {
        case MFCC:
            m_Writers.push_back([this](const Description &desc, float *&out) {
                for (int i = 0; i < m_MFCCCount; ++i)
                    *out++ = desc.MFCC[i];
            });
            break;
        case CHROMA:
            m_Writers.push_back([this](const Description &desc, float *&out) {
                for (int i = 0; i < m_ChromaSize; ++i)
                    *out++ = desc.Chroma[i];
            });
            break;
        case MELOGRAM:
            m_Writers.push_back([this](const Description &desc, float *&out) {
                for (int i = 0; i < m_MFCCMels; ++i)
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
        case EXTENDEDPROB:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.ExtendedTechProb; });
            break;
        case ONSET:
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
        case KURTOSIS:
            m_Writers.push_back([](const Description &desc, float *&out) { *out++ = desc.SpectralKurtosis; });
            break;
        case ONNX:
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
    if (!m_ONNXModelLoaded) {
        return;
    }

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
            Desc.ONNX[m_ONNXLabels[i]] = data[i] * Desc.ExtendedTechProb;
        }
        break;
    }
    case ONNX_TENSOR_TYPE_FLOAT64: {
        double *data = (double *)m_OutputTensor->datas;
        for (size_t i = 0; i < m_ONNXLabels.size(); i++) {
            Desc.ONNX[m_ONNXLabels[i]] = (float)data[i] * Desc.ExtendedTechProb;
        }
        break;
    }
    default:
        spdlog::error("Tensor output type not supported {}", (int)m_OutputTensor->type);
    }
}

// ╭─────────────────────────────────────╮
// │           Onset Detector            │
// ╰─────────────────────────────────────╯
void MIR::OnsetInit() {
    m_OnsetInit = false;

    const int onsetFFTSize = std::max(2, m_OnsetFFTSize);
    const size_t nbytes = onsetsds_memneeded(ODS_ODF_MKL, onsetFFTSize, m_MedSpan);
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

    onsetsds_init(m_ODS, m_ODSData, ODS_FFT_FFTW3_R2C, ODS_ODF_MKL, onsetFFTSize, m_MedSpan, m_Sr);
    m_OnsetFFTFrame.assign(static_cast<size_t>(2 * (onsetFFTSize / 2 + 1)), 0.0f);
    m_OnsetInit = true;
}

// ─────────────────────────────────────
void MIR::OnsetExec(Description &Desc) {
    if (!m_OnsetInit)
        return;

    const size_t nBins = static_cast<size_t>(m_OnsetFFTSize / 2 + 1);
    const size_t requiredSize = 2 * nBins;
    if (m_OnsetFFTFrame.size() != requiredSize) {
        m_OnsetFFTFrame.assign(requiredSize, 0.0f);
    }

    // TODO: Would be great without another copy
    for (size_t i = 0; i < nBins; ++i) {
        m_OnsetFFTFrame[2 * i] = static_cast<float>(m_FullFFTOut[i][0]);
        m_OnsetFFTFrame[2 * i + 1] = static_cast<float>(m_FullFFTOut[i][1]);
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
    Desc.ExtendedTechProb *= (1.0f - Desc.PitchConfidence);
    Desc.ExtendedTechProb *= abs(m_ODS->odfvalpost);
    float steepness = 5.0f;
    Desc.ExtendedTechProb = 1.0f / (1.0f + std::exp(-steepness * (Desc.ExtendedTechProb - 0.5f)));
}

// ╭─────────────────────────────────────╮
// │         Power and Amplitude         │
// ╰─────────────────────────────────────╯
// Check https://github.com/klangfreund/LUFSMeter
void MIR::InitITURFilters() {
    // Stage 1: shelving filter
    double KoverQ1 = (2.0 - 2.0 * m_48kA1[2]) / (m_48kA1[2] - m_48kA1[1] + 1.0);
    double K1 = std::sqrt((m_48kA1[1] + m_48kA1[2] + 1.0) / (m_48kA1[2] - m_48kA1[1] + 1.0));
    double Q1 = K1 / KoverQ1;
    double arctanK1 = std::atan(K1);
    double VB1 = (m_48kB1[0] - m_48kB1[2]) / (1.0 - m_48kA1[2]);
    double VH1 = (m_48kB1[0] - m_48kB1[1] + m_48kB1[2]) / (m_48kA1[2] - m_48kA1[1] + 1.0);
    double VL1 = (m_48kB1[0] + m_48kB1[1] + m_48kB1[2]) / (m_48kA1[1] + m_48kA1[2] + 1.0);

    double Knew1 = std::tan(arctanK1 * 48000.0 / m_Sr);
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

    double Knew2 = std::tan(arctanK2 * 48000.0 / m_Sr);
    double commonFactor2 = 1.0 / (1.0 + Knew2 / Q2 + Knew2 * Knew2);

    m_B2[0] = (VH2 + VB2 * Knew2 / Q2 + VL2 * Knew2 * Knew2) * commonFactor2;
    m_B2[1] = 2.0 * (VL2 * Knew2 * Knew2 - VH2) * commonFactor2;
    m_B2[2] = (VH2 - VB2 * Knew2 / Q2 + VL2 * Knew2 * Knew2) * commonFactor2;
    m_A2[0] = 1.0;
    m_A2[1] = 2.0 * (Knew2 * Knew2 - 1.0) * commonFactor2;
    m_A2[2] = (1.0 - Knew2 / Q2 + Knew2 * Knew2) * commonFactor2;
}

// ─────────────────────────────────────
void MIR::GetSignalPower(std::vector<double> &In, Description &Desc) {
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
    const size_t frameSize = static_cast<size_t>(std::max(2.0f, m_FFTSize));
    const size_t half = frameSize / 2;
    const size_t allocSize = half + 2;
    m_YINDifference.assign(allocSize, 0.0);
    m_YINCMNDF.assign(allocSize, 1.0);
}

// ─────────────────────────────────────
void MIR::YINExec(std::vector<double> &In, Description &Desc) {
    const size_t frame = In.size();
    const size_t minTau = m_Sr / m_YINMaxFrequency;
    const size_t maxTauByPitch = std::ceil(m_Sr / m_YINMinFrequency);
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
        if (value < m_YINThreshold) {
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

    const double Pitch = m_Sr / refinedTau;
    if (Pitch < m_YINMinFrequency || Pitch > m_YINMaxFrequency) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    Desc.Pitch = Pitch;
    Desc.PitchConfidence = Confidence;
}

// ─────────────────────────────────────
void MIR::GetSpectralDescriptions(Description &Desc) {
    const size_t NHalf = m_FFTSize / 2 + 1;
    const double binWidth = static_cast<double>(m_Sr) / static_cast<double>(m_FFTSize);
    const size_t hfStart = NHalf / 4;
    const double invN = 1.0 / static_cast<double>(m_FFTSize);
    constexpr double amin = 1e-10;

    fftw_execute(m_FullFFTPlan);

    Desc.MaxAmp = 0.0;
    Desc.SpectralFlux = 0.0;

    double SumPower = 0.0;
    double weightedSumFreqs = 0.0;
    double irregularityNumerator = 0.0;
    double irregularityDenominator = 0.0;
    double logSumPower = 0.0;
    double linSumPower = 0.0;
    double harmonicityPeak = 0.0;
    double harmonicitySum = 0.0;
    double highFreqEnergy = 0.0;

    double sumFreq2 = 0.0;
    double sumFreq3 = 0.0;
    double sumFreq4 = 0.0;
    double sumIndex = 0.0;
    double sumIndex2 = 0.0;

    // Spectral Slope
    double sumFreq = 0.0;
    double sumFreqSq = 0.0;

    // Crest
    double maxMag = 0.0;
    float sumMagCrest = 0.0f;
    float maxMagCrest = 0.0f;

    // Essentia Flux keeps a previous spectrum with the same number of bins.
    if (m_PreviousSpectralPower.size() != NHalf) {
        m_PreviousSpectralPower.assign(NHalf, 0.0);
    }

    for (size_t i = 0; i < NHalf; ++i) {
        const double re = m_FullFFTOut[i][0];
        const double im = m_FullFFTOut[i][1];

        const double p = re * re + im * im;
        const double mag = std::sqrt(p);
        const double norm = mag * invN;

        Desc.Power[i] = p;
        Desc.Magnitude[i] = mag;
        Desc.SpectralMagnitudeNorm[i] = norm;
        Desc.MaxAmp = std::max(Desc.MaxAmp, norm);

        maxMag = std::max(maxMag, mag);
        const float magCrest = static_cast<float>(mag);
        sumMagCrest += magCrest;
        maxMagCrest = std::max(maxMagCrest, magCrest);

        SumPower += norm;
        const double freq = i * binWidth;
        const double freq2 = freq * freq;
        const double index = static_cast<double>(i);
        sumFreq += freq;
        sumFreqSq += freq2;

        weightedSumFreqs += freq * norm;
        sumIndex += index * norm;
        sumIndex2 += index * index * norm;
        irregularityDenominator += norm * norm;

        const double v = std::max(amin, p);
        logSumPower += std::log(v);
        linSumPower += v;

        const double diff = mag - m_PreviousSpectralPower[i];
        Desc.SpectralFlux += diff * diff;

        sumFreq2 += freq2 * norm;
        sumFreq3 += freq2 * freq * norm;
        sumFreq4 += freq2 * freq2 * norm;

        if (i >= hfStart) {
            highFreqEnergy += norm;
        }
        if (i > 0) {
            const double d = Desc.SpectralMagnitudeNorm[i - 1] - norm;
            irregularityNumerator += d * d;
            harmonicitySum += norm;
            harmonicityPeak = std::max(harmonicityPeak, norm);
        }

        m_PreviousSpectralPower[i] = mag;
    }

    Desc.SpectralFlux = std::sqrt(Desc.SpectralFlux);

    const double SumPowerEps = SumPower + 1e-12;
    const double invSum = 1.0 / SumPowerEps;

    Desc.SpectralCentroid = weightedSumFreqs / SumPowerEps;
    Desc.CentroidVelocity = std::abs(Desc.SpectralCentroid - m_PrevCentroid);
    m_PrevCentroid = Desc.SpectralCentroid;

    const double Mean = 1.0 / static_cast<double>(NHalf);
    double Variance = 0.0;

    for (size_t i = 0; i < NHalf; ++i) {
        const double normSp = (Desc.SpectralMagnitudeNorm[i] + 1e-12) / SumPowerEps;
        Desc.SpectralMagnitudeFrameNorm[i] = normSp;

        const double diff = normSp - Mean;
        Variance += diff * diff;
    }

    const double E1 = Desc.SpectralCentroid;
    const double E2 = sumFreq2 * invSum;
    const double E3 = sumFreq3 * invSum;
    const double E4 = sumFreq4 * invSum;
    const double mu3 = E3 - 3.0 * E1 * E2 + 2.0 * E1 * E1 * E1;
    const double mu4 = E4 - 4.0 * E1 * E3 + 6.0 * E1 * E1 * E2 - 3.0 * E1 * E1 * E1 * E1;

    const double spectralSpreadHzVariance = std::max(0.0, E2 - E1 * E1);
    Desc.SpectralSpreadHz = std::sqrt(spectralSpreadHzVariance);

    const double EIndex = sumIndex * invSum;
    const double EIndex2 = sumIndex2 * invSum;
    const double rawIndexVariance = std::max(0.0, EIndex2 - EIndex * EIndex);
    const double indexRange = static_cast<double>(std::max<size_t>(1, NHalf - 1));
    const double invIndexRange = 1.0 / indexRange;
    Desc.SpectralSpreadVariance = rawIndexVariance * invIndexRange * invIndexRange;

    if (SumPower <= 0.0) {
        Desc.SpectralSpreadHz = 0.0;
        Desc.SpectralSpreadVariance = 0.0;
    }

    const double spreadEps = Desc.SpectralSpreadHz + 1e-12;
    Desc.SpectralSkewness = mu3 / (spreadEps * spreadEps * spreadEps);
    Desc.SpectralKurtosis = mu4 / (spreadEps * spreadEps * spreadEps * spreadEps) - 3.0;

    Desc.StdDev = std::sqrt(Variance / static_cast<double>(NHalf));
    Desc.SpectralIrregularity = irregularityDenominator > 0.0 ? (irregularityNumerator / irregularityDenominator) : 0.0;
    if (maxMagCrest == 0.0f) {
        Desc.SpectralCrest = 0.0;
    } else {
        const float meanMagCrest = sumMagCrest / static_cast<float>(NHalf);
        Desc.SpectralCrest = meanMagCrest > 0.0f ? static_cast<double>(maxMagCrest / meanMagCrest) : 0.0;
    }

    Desc.SpectralFlatness = std::exp(logSumPower / NHalf) / (linSumPower / NHalf);
    Desc.Harmonicity = harmonicityPeak / (harmonicitySum + 1e-12);
    Desc.HighFreqRatio = highFreqEnergy / SumPowerEps;

    const double K = static_cast<double>(NHalf);
    const double denom = (K * sumFreqSq - sumFreq * sumFreq) + 1e-12;
    const double slope = (K * weightedSumFreqs - sumFreq * SumPower) / denom;
    Desc.SpectralSlope = slope / SumPowerEps;

    m_SpectralPrefix.resize(NHalf + 1);
    m_SpectralPrefix[0] = 0.0;

    for (size_t i = 0; i < NHalf; ++i) {
        m_SpectralPrefix[i + 1] = m_SpectralPrefix[i] + Desc.Power[i];
    }

    OnsetExec(Desc);
    MFCCExec(Desc);
    SpectralChromaExec(Desc);
}

// ╭─────────────────────────────────────╮
// │                MFCC                 │
// ╰─────────────────────────────────────╯
void MIR::MFCCInit() {
    const int n_f = m_FFTSize / 2 + 1;
    const double fmin = 0.0;
    const double fmax = static_cast<double>(m_Sr) * 0.5;

    // FFT bin frequencies (match np.fft.rfftfreq)
    std::vector<double> fft_freqs(n_f);
    for (int i = 0; i < n_f; ++i) {
        fft_freqs[i] = static_cast<double>(i) * static_cast<double>(m_Sr) / static_cast<double>(m_FFTSize);
    }

    // Slaney mel scale (htk=False)
    const double f_sp = 200.0 / 3.0;
    const double min_log_hz = 1000.0;
    const double min_log_mel = (min_log_hz - fmin) / f_sp;
    const double logstep = std::log(6.4) / 27.0;

    auto hz_to_mel = [&](double hz) {
        if (hz < min_log_hz)
            return (hz - fmin) / f_sp;
        return min_log_mel + std::log(hz / min_log_hz) / logstep;
    };

    auto mel_to_hz = [&](double mel) {
        if (mel < min_log_mel)
            return mel * f_sp + fmin;
        return min_log_hz * std::exp((mel - min_log_mel) * logstep);
    };

    const double mel_min = hz_to_mel(fmin);
    const double mel_max = hz_to_mel(fmax);

    std::vector<double> mel_pts(m_MFCCMels + 2);
    for (int i = 0; i < m_MFCCMels + 2; ++i) {
        mel_pts[i] = mel_min + (mel_max - mel_min) * static_cast<double>(i) / static_cast<double>(m_MFCCMels + 1);
    }

    std::vector<double> hz_pts(m_MFCCMels + 2);
    for (int i = 0; i < m_MFCCMels + 2; ++i) {
        hz_pts[i] = mel_to_hz(mel_pts[i]);
    }

    // Slaney triangular filters (area normalized), matching librosa.filters.mel
    m_MFCCFilter.assign(m_MFCCMels, std::vector<double>(n_f, 0.0));
    m_MFCCActiveBins.assign(m_MFCCMels, {0, -1});

    std::vector<double> fdiff(m_MFCCMels + 1);
    for (int i = 0; i < m_MFCCMels + 1; ++i) {
        fdiff[i] = hz_pts[i + 1] - hz_pts[i];
    }

    for (int m = 0; m < m_MFCCMels; ++m) {
        const double enorm = 2.0 / (hz_pts[m + 2] - hz_pts[m]);

        int first = -1;
        int last = -1;
        for (int k = 0; k < n_f; ++k) {
            const double ramp = hz_pts[m] - fft_freqs[k];
            const double lower = -ramp / fdiff[m];
            const double upper = (hz_pts[m + 2] - fft_freqs[k]) / fdiff[m + 1];
            const double w = std::max(0.0, std::min(lower, upper));
            const double v = w * enorm;
            m_MFCCFilter[m][k] = v;

            if (v > 0.0) {
                if (first < 0)
                    first = k;
                last = k;
            }
        }

        m_MFCCActiveBins[m] = {first < 0 ? 0 : first, last};
    }

    // Orthonormal DCT-II (librosa: dct_type=2, norm="ortho")
    m_DCTBasis.assign(m_MFCCCount, std::vector<double>(m_MFCCMels));

    const double scale0 = std::sqrt(1.0 / m_MFCCMels);
    const double scale = std::sqrt(2.0 / m_MFCCMels);

    for (int k = 0; k < m_MFCCCount; ++k) {
        for (int n = 0; n < m_MFCCMels; ++n) {
            m_DCTBasis[k][n] = (k == 0 ? scale0 : scale) * std::cos(std::numbers::pi * (n + 0.5) * k / m_MFCCMels);
        }
    }

    m_MFCCEnergy.resize(m_MFCCMels);
}

// ─────────────────────────────────────
void MIR::MFCCExec(Description &Desc) {
    const int n_mels = std::min<int>(m_MFCCMels, static_cast<int>(m_MFCCFilter.size()));
    const int n_mfcc = std::min<int>(m_MFCCCount, static_cast<int>(m_DCTBasis.size()));

    if (n_mels <= 0 || n_mfcc <= 0 || Desc.Power.empty()) {
        Desc.MFCC.assign(std::max(0, n_mfcc), 0.0);
        Desc.LogMelSpectrum.assign(n_mels, 0.0);
        return;
    }

    const int n_f = std::min<int>(static_cast<int>(Desc.Power.size()), static_cast<int>(m_MFCCFilter[0].size()));
    if (n_f <= 0) {
        Desc.MFCC.assign(n_mfcc, 0.0);
        Desc.LogMelSpectrum.assign(n_mels, 0.0);
        return;
    }

    constexpr double kAmin = 1e-10f;
    constexpr double kTopDb = 80.0f;

    if (static_cast<int>(m_MFCCEnergy.size()) != n_mels)
        m_MFCCEnergy.resize(n_mels);

    // Mel projection (power domain) -> fill m_MFCCEnergy and Melogram.
    Desc.LogMelSpectrum.resize(n_mels);
    for (int m = 0; m < n_mels; ++m) {
        const double *filter = m_MFCCFilter[m].data();
        double melEnergy = 0.0;
        for (int k = 0; k < n_f; ++k) {
            melEnergy += filter[k] * Desc.Power[k];
        }

        m_MFCCEnergy[m] = melEnergy;

        // Match librosa.power_to_db(ref=1.0): 10 * log10(max(amin, S)).
        Desc.LogMelSpectrum[m] = 10.0 * std::log10(std::max(kAmin, melEnergy));
    }

    // Apply librosa-style top_db floor.
    float maxLog = -std::numeric_limits<float>::infinity();
    for (int m = 0; m < n_mels; ++m)
        if (Desc.LogMelSpectrum[m] > maxLog)
            maxLog = static_cast<float>(Desc.LogMelSpectrum[m]);

    if (std::isfinite(maxLog)) {
        const float floor = maxLog - kTopDb;
        for (int m = 0; m < n_mels; ++m) {
            if (Desc.LogMelSpectrum[m] < floor)
                Desc.LogMelSpectrum[m] = floor;
        }
    }

    // MFCCs are DCT-II over the log-mel (dB) spectrum.
    Desc.MFCC.assign(n_mfcc, 0.0);
    for (int k = 0; k < n_mfcc; ++k) {
        const double *basis = m_DCTBasis[k].data();
        double coeff = 0.0;
        for (int n = 0; n < n_mels; ++n)
            coeff += basis[n] * Desc.LogMelSpectrum[n];
        Desc.MFCC[k] = coeff;
    }
}

// ╭─────────────────────────────────────╮
// │               Chroma                │
// ╰─────────────────────────────────────╯
double MIR::HzToOcts(double frequency, double tuning, int binsPerOctave) const {
    const double a440 = m_ChromaA440 * std::pow(2.0, tuning / static_cast<double>(binsPerOctave));
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
    const size_t nFft = static_cast<size_t>(m_FFTSize);
    const size_t nHalf = nFft / 2 + 1;

    m_ChromaFilter.assign(m_ChromaSize, std::vector<double>(nHalf, 0.0));
    if (nFft == 0 || nHalf == 0 || m_Sr <= 0.0f) {
        return;
    }

    std::vector<double> frqbins(nFft, 0.0);
    if (nFft > 1) {
        for (size_t k = 1; k < nFft; ++k) {
            const double frequency = static_cast<double>(k) * static_cast<double>(m_Sr) / static_cast<double>(nFft);
            frqbins[k] =
                static_cast<double>(m_ChromaSize) * HzToOcts(frequency, m_ChromaTuning, static_cast<int>(m_ChromaSize));
        }
        frqbins[0] = frqbins[1] - 1.5 * static_cast<double>(m_ChromaSize);
    } else {
        frqbins[0] = -1.5 * static_cast<double>(m_ChromaSize);
    }

    std::vector<double> binwidthbins(nFft, 1.0);
    for (size_t k = 0; k + 1 < nFft; ++k) {
        binwidthbins[k] = std::max(frqbins[k + 1] - frqbins[k], 1.0);
    }

    const double nChroma2 = std::round(static_cast<double>(m_ChromaSize) / 2.0);
    for (size_t k = 0; k < nHalf; ++k) {
        double columnNorm = 0.0;
        for (int chroma = 0; chroma < m_ChromaSize; ++chroma) {
            const double distance = PositiveRemainder(frqbins[k] - static_cast<double>(chroma) + nChroma2 +
                                                          10.0 * static_cast<double>(m_ChromaSize),
                                                      static_cast<double>(m_ChromaSize)) -
                                    nChroma2;
            const double weight = std::exp(-0.5 * std::pow(2.0 * distance / binwidthbins[k], 2.0));
            m_ChromaFilter[chroma][k] = weight;
            columnNorm += weight * weight;
        }

        if (columnNorm > 0.0) {
            const double invNorm = 1.0 / std::sqrt(columnNorm);
            for (int chroma = 0; chroma < m_ChromaSize; ++chroma) {
                m_ChromaFilter[chroma][k] *= invNorm;
            }
        }

        const double octaveWeight =
            std::exp(-0.5 * std::pow((frqbins[k] / static_cast<double>(m_ChromaSize) - m_ChromaCenterOctave) /
                                         m_ChromaOctaveWidth,
                                     2.0));
        for (int chroma = 0; chroma < m_ChromaSize; ++chroma) {
            m_ChromaFilter[chroma][k] *= octaveWeight;
        }
    }

    const int chromaShift = 3 * (m_ChromaSize / 12);
    if (chromaShift > 0 && chromaShift < m_ChromaSize) {
        Matrix rolled(m_ChromaSize, std::vector<double>(nHalf, 0.0));
        for (int chroma = 0; chroma < m_ChromaSize; ++chroma) {
            rolled[chroma] = m_ChromaFilter[(chroma + chromaShift) % m_ChromaSize];
        }
        m_ChromaFilter.swap(rolled);
    }
}

// ─────────────────────────────────────
void MIR::SpectralChromaExec(Description &Desc) {

    // TODO:
    if (Desc.Chroma.size() != static_cast<size_t>(m_ChromaSize))
        Desc.Chroma.resize(m_ChromaSize);

    std::fill(Desc.Chroma.begin(), Desc.Chroma.end(), 0.0);

    if (m_ChromaFilter.empty() || m_ChromaFilter[0].empty() || Desc.Power.empty()) {
        return;
    }

    const size_t nHalf = std::min(Desc.Power.size(), m_ChromaFilter[0].size());
    for (int chroma = 0; chroma < m_ChromaSize; ++chroma) {
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
    m_ZCRFrameLength = static_cast<int>(m_FFTSize);
    m_ZCRHopLength = static_cast<int>(m_HopSize);

    if (m_ZCRFrameLength <= 0) {
        m_ZCRFrameLength = 1;
    }
    if (m_ZCRHopLength <= 0) {
        m_ZCRHopLength = 1;
    }

    // Reuse this buffer across calls to avoid per-block allocations.
    const size_t pad = m_ZCRCenter ? static_cast<size_t>(m_ZCRFrameLength / 2) : 0;
    m_ZCRScratch.resize(static_cast<size_t>(m_ZCRFrameLength) + (2 * pad));
}

// ─────────────────────────────────────
void MIR::ZeroCrossingRateExec(std::vector<double> &In, Description &Desc) {
    if (In.empty()) {
        Desc.ZeroCrossingRate = 0.0;
        return;
    }

    const int frameLength = std::max(1, m_ZCRFrameLength);
    const size_t frameLengthSz = static_cast<size_t>(frameLength);
    const double *yData = nullptr;
    size_t ySize = 0;

    if (m_ZCRCenter) {
        const size_t pad = frameLengthSz / 2;
        const size_t inSize = In.size();
        ySize = inSize + (2 * pad);
        if (m_ZCRScratch.size() < ySize) {
            m_ZCRScratch.resize(ySize);
        }

        double *dst = m_ZCRScratch.data();
        const double edgeLeft = In.front();
        const double edgeRight = In.back();
        std::fill_n(dst, pad, edgeLeft);
        std::copy(In.begin(), In.end(), dst + pad);
        std::fill_n(dst + pad + inSize, pad, edgeRight);
        yData = dst;
    } else {
        yData = In.data();
        ySize = In.size();
    }

    if (ySize < frameLengthSz) {
        Desc.ZeroCrossingRate = 0.0;
        return;
    }

    size_t crossings = 0;
    if (m_ZCRPad) {
        crossings += 1;
    }

    const double threshold = m_ZCRThreshold;
    if (m_ZCRZeroPos) {
        double prev = yData[0];
        if (std::abs(prev) <= threshold) {
            prev = 0.0;
        }

        for (int i = 1; i < frameLength; ++i) {
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
        for (int i = 1; i < frameLength; ++i) {
            double curr = yData[static_cast<size_t>(i)];
            if (std::abs(curr) <= threshold)
                curr = 0.0;

            const int currSign = (curr > 0.0) - (curr < 0.0);
            crossings += static_cast<size_t>(prevSign != currSign);
            prevSign = currSign;
        }
    }

    Desc.ZeroCrossingRate = static_cast<double>(crossings) / static_cast<double>(frameLengthSz);
}

// ╭─────────────────────────────────────╮
// │        Spectral Descriptions        │
// ╰─────────────────────────────────────╯
void MIR::SpectralFluxExec(Description &Desc) {

    // Must be magnitude spectrum, not power
    const auto &S = Desc.SpectralMagnitudeNorm;
    const size_t N = S.size();

    if (m_PreviousSpectralPower.size() != N) {
        m_PreviousSpectralPower.assign(N, 0.0);
    }

    double flux = 0.0;

    // Essentia uses all bins
    for (size_t k = 0; k < N; ++k) {
        double diff = S[k] - m_PreviousSpectralPower[k];
        flux += diff * diff;
    }

    flux = std::sqrt(flux);

    // store spectrum for next frame
    m_PreviousSpectralPower = S;

    Desc.SpectralFlux = flux;
}

// ─────────────────────────────────────
void MIR::SpectralFlatnessExec(Description &Desc) {
    const auto &S = Desc.Power; // usar potência
    const size_t N = S.size();

    constexpr double amin = 1e-10; // mesmo default do librosa

    double logSum = 0.0;
    double linSum = 0.0;

    for (size_t k = 0; k < N; ++k) {
        double v = std::max(amin, S[k]); // já está |X|^2
        logSum += std::log(v);
        linSum += v;
    }

    const double geoMean = std::exp(logSum / double(N));
    const double arithMean = linSum / double(N);

    Desc.SpectralFlatness = geoMean / arithMean;
}

// ─────────────────────────────────────
void MIR::SpectralHarmonicityExec(Description &Desc) {
    const auto &S = Desc.SpectralMagnitudeNorm;
    const size_t N = S.size();

    double sum = 0.0;
    double peak = 0.0;

    // skip DC
    for (size_t k = 1; k < N; ++k) {
        double v = S[k];
        sum += v;
        if (v > peak)
            peak = v;
    }

    constexpr double eps = 1e-12;
    Desc.Harmonicity = peak / (sum + eps);
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
void MIR::GetDescription(std::vector<double> &In, Description &Desc) {
    // 1. Temporal Domain
    GetSignalPower(In, Desc);
    YINExec(In, Desc);
    ZeroCrossingRateExec(In, Desc);

    // 2. Frequency Domain (windowing + FFT)
    const double *x = In.data();
    const double *w = m_FullWindowingFunc.data();
    for (size_t i = 0; i < m_FFTSize; ++i) {
        m_FullFFTIn[i] = x[i] * w[i];
    }

    GetSpectralDescriptions(Desc);

    // Perc
    ExtendedTechExec(Desc);

    // ONNX
    ONNXExec(Desc);
}

} // namespace OpenScofo
