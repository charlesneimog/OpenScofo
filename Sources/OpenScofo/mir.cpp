#include "OpenScofo.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
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
    if (m_FFTPlan != nullptr) {
        fftw_destroy_plan(m_FFTPlan);
        m_FFTPlan = nullptr;
    }

    if (m_FFTIn != nullptr) {
        fftw_free(m_FFTIn);
        m_FFTIn = nullptr;
    }
    if (m_FFTOut != nullptr) {
        fftw_free(m_FFTOut);
        m_FFTOut = nullptr;
    }

    if (m_OnsetInit) {
        delete[] m_ODSData;
        delete m_ODS;
    }

    if (m_ONNXModelLoaded) {
        if (m_OnnxCTX) {
            onnx_context_free(m_OnnxCTX);
            m_OnnxCTX = nullptr;
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
    if (sameAudioConfig && m_FFTPlan != nullptr && m_FFTIn != nullptr && m_FFTOut != nullptr &&
        m_WindowingFunc.size() == static_cast<size_t>(FftSize)) {
        return;
    }

    m_HopSize = HopSize;
    m_FFTSize = FftSize;
    m_Sr = Sr;
    m_BlockSize = HopSize;
    m_Accum = 0;
    m_PrevCentroid = 0.0;

    if (m_FFTPlan != nullptr) {
        fftw_destroy_plan(m_FFTPlan);
        m_FFTPlan = nullptr;
    }
    if (m_FFTIn != nullptr) {
        fftw_free(m_FFTIn);
        m_FFTIn = nullptr;
    }
    if (m_FFTOut != nullptr) {
        fftw_free(m_FFTOut);
        m_FFTOut = nullptr;
    }

    m_PreviousSpectralPower.assign(static_cast<size_t>(round(m_FFTSize / 2)), 0.0);

    FFTWInit();
    OnsetInit();
    MFCCInit();
    SpectralChromaInit();
    ZeroCrossingRateInit();
    YINInit();

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
    m_FFTIn = fftw_alloc_real(m_FFTSize);
    if (!m_FFTIn) {
        spdlog::critical("fftw_alloc_real failed");
        return;
    }

    m_FFTOut = fftw_alloc_complex(WindowHalf + 1);
    if (!m_FFTOut) {
        fftw_free(m_FFTIn);
        spdlog::critical("fftw_alloc_complex failed");
        return;
    }

    m_FFTPlan = fftw_plan_dft_r2c_1d((int)m_FFTSize, m_FFTIn, m_FFTOut, FFTW_PATIENT);

    // Match librosa/scipy get_window('hann', N, fftbins=True): periodic Hann.
    m_WindowingFunc.resize(m_FFTSize);
    for (size_t i = 0; i < m_FFTSize; i++) {
        m_WindowingFunc[i] = 0.5 * (1.0 - cos(2.0 * std::numbers::pi * i / m_FFTSize));
    }
}

// ╭─────────────────────────────────────╮
// │          Machine Learning           │
// ╰─────────────────────────────────────╯
void MIR::LoadONNXModel(fs::path path) {
    auto u8 = path.u8string();
    std::string path_utf8(u8.begin(), u8.end());

    m_OnnxCTX = onnx_context_alloc_from_file(path_utf8.c_str(), nullptr, 0);

    if (m_OnnxCTX == nullptr) {
        spdlog::error("Failed to load ONNX model: {}.", path.string());
        return;
    }

    if (m_OnnxCTX && m_OnnxCTX->g != nullptr) {
        struct onnx_graph_t *g = m_OnnxCTX->g;
        for (int i = 0; i < g->nlen; i++) {
            struct onnx_node_t *n = &g->nodes[i];
            if (n->opset > CURRENT_ONNX_OPSET) {
                spdlog::error("Unsupported opset => {} {}.", n->proto->op_type, n->opset);
                return;
            }
        }
    }

    // Get labels
    bool TreeEnsembleClassifierFound = false;
    struct onnx_graph_t *g = m_OnnxCTX->g;
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
                    m_ONNXResults[label] = 0.0f;
                }
            }
        }
    }

    if (!TreeEnsembleClassifierFound) {
        spdlog::error("TreeEnsembleClassifier not found in model, please use py.train to train models for OpenScofo");
        return;
    }

    m_ONNXModelLoaded = true;
}

// ╭─────────────────────────────────────╮
// │           Onset Detector            │
// ╰─────────────────────────────────────╯
void MIR::OnsetInit() {
    size_t nbytes = onsetsds_memneeded(ODS_ODF_MKL, m_FFTSize, m_MedSpan);
    delete[] m_ODSData;
    delete m_ODS;
    m_ODSData = new float[nbytes / sizeof(float)];
    m_ODS = new OnsetsDS();
    onsetsds_init(m_ODS, m_ODSData, ODS_FFT_FFTW3_R2C, ODS_ODF_MKL, m_OnsetFFTSize, m_MedSpan, m_Sr);
    if (!m_ODS) {
        spdlog::critical("Not possible to initialize the onset detector");
        return;
    }
    m_OnsetInit = true;
}

// ─────────────────────────────────────
void MIR::OnsetExec(Description &Desc) {
    if (!m_OnsetInit) {
        return;
    }
    Desc.Onset = onsetsds_process(m_ODS, reinterpret_cast<float *>(m_FFTOut));
}

// ╭─────────────────────────────────────╮
// │          Audio Descriptors          │
// ╰─────────────────────────────────────╯
void MIR::GetSignalPower(std::vector<double> &In, Description &Desc) {
    // TODO: Recalculate this based on the sample rate
    const std::array<double, 3> b1 = {1.53512485958697, -2.69169618940638, 1.19839281085285};
    const std::array<double, 3> a1 = {1.0, -1.69065929318241, 0.73248077421585};
    const std::array<double, 3> b2 = {1.0, -2.0, 1.0};
    const std::array<double, 3> a2 = {1.0, -1.99004745483398, 0.99007225036621};

    // Filter states (local variables)
    double x1_1 = 0.0, x2_1 = 0.0;
    double y1_1 = 0.0, y2_1 = 0.0;
    double x1_2 = 0.0, x2_2 = 0.0;
    double y1_2 = 0.0, y2_2 = 0.0;

    double z = 0.0;
    double z_loudness = 0.0;
    for (double sample : In) {
        double s1 = b1[0] * sample + b1[1] * x1_1 + b1[2] * x2_1 - a1[1] * y1_1 - a1[2] * y2_1;
        x2_1 = x1_1;
        x1_1 = sample;
        y2_1 = y1_1;
        y1_1 = s1;
        double s2 = b2[0] * s1 + b2[1] * x1_2 + b2[2] * x2_2 - a2[1] * y1_2 - a2[2] * y2_2;
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

void MIR::YINExec(std::vector<double> &In, Description &Desc) {
    const size_t frame = In.size();
    if (frame < 2 || m_YINDifference.empty() || m_YINCMNDF.empty()) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    const double sampleRate = std::max(1.0, static_cast<double>(m_Sr));
    const size_t minTau = std::max<size_t>(1, static_cast<size_t>(sampleRate / m_YINMaxFrequency));
    const size_t maxTauByPitch = std::max(minTau + 1, static_cast<size_t>(std::ceil(sampleRate / m_YINMinFrequency)));
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

    const double confidence = std::clamp(1.0 - bestValue, 0.0, 1.0);
    if (refinedTau <= 0.0 || confidence <= 0.0) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    const double pitch = sampleRate / refinedTau;
    if (pitch < m_YINMinFrequency || pitch > m_YINMaxFrequency) {
        Desc.Pitch = 0.0;
        Desc.PitchConfidence = 0.0;
        return;
    }

    Desc.Pitch = pitch;
    Desc.PitchConfidence = confidence;
}

// ─────────────────────────────────────
void MIR::GetSpectralDescriptions(Description &Desc) {
    const size_t NHalf = m_FFTSize / 2 + 1;
    const double binWidth = static_cast<double>(m_Sr) / static_cast<double>(m_FFTSize);
    const size_t hfStart = NHalf / 4;

    fftw_execute(m_FFTPlan);

    Desc.MaxAmp = 0;
    double SumPower = 0.0;
    const double invN = 1.0 / static_cast<double>(m_FFTSize);

    // FUSE 1: Calculate power-domain arrays and descriptors in one pass.
    double logSumPower = 0.0;
    double linSumPower = 0.0;
    double flux = 0.0;
    double irregularityNumerator = 0.0;
    double irregularityDenominator = 0.0;
    double harmonicityPeak = 0.0;
    double harmonicitySum = 0.0;
    double highFreqEnergy = 0.0;
    double weightedSumFreqs = 0.0;
    constexpr double amin = 1e-10;
    double previousSpectralBin = 0.0;
    bool hasPreviousSpectralBin = false;

    if (m_PreviousSpectralPower.size() != NHalf) {
        m_PreviousSpectralPower.assign(NHalf, 0.0);
    }

    for (size_t i = 0; i < NHalf; ++i) {
        const double re = m_FFTOut[i][0];
        const double im = m_FFTOut[i][1];
        const double p = re * re + im * im;
        const double sp = std::sqrt(p) * invN;

        Desc.Power[i] = p;
        Desc.SpectralPower[i] = sp;

        if (sp > Desc.MaxAmp)
            Desc.MaxAmp = sp;

        SumPower += sp;
        weightedSumFreqs += (static_cast<double>(i) * binWidth) * sp;
        irregularityDenominator += sp * sp;
        if (i >= hfStart) {
            highFreqEnergy += sp;
        }

        if (hasPreviousSpectralBin) {
            const double binDelta = previousSpectralBin - sp;
            irregularityNumerator += binDelta * binDelta;
        }
        previousSpectralBin = sp;
        hasPreviousSpectralBin = true;

        // Flatness on power spectrum.
        const double v = std::max(amin, p);
        logSumPower += std::log(v);
        linSumPower += v;

        // Flux and harmonicity skip DC to match previous behavior.
        if (i > 0) {
            const double diff = sp - m_PreviousSpectralPower[i];
            if (diff > 0.0) {
                flux += diff;
            }

            harmonicitySum += sp;
            if (sp > harmonicityPeak) {
                harmonicityPeak = sp;
            }
        }
        m_PreviousSpectralPower[i] = sp;
    }

    const bool hasSpectralEnergy = SumPower > 1e-12;
    const double centroid = hasSpectralEnergy ? (weightedSumFreqs / SumPower) : 0.0;
    Desc.SpectralCentroid = centroid;
    Desc.CentroidVelocity = std::abs(centroid - m_PrevCentroid);
    m_PrevCentroid = centroid;

    // FUSE 2: Normalize and calculate mean/variance simultaneously.
    const double SumPowerEps = SumPower + 1e-12;
    const double Mean = 1.0 / NHalf;
    double Variance = 0.0;
    double weightedSpreadVariance = 0.0;

    for (size_t i = 0; i < NHalf; ++i) {
        double normSp = (Desc.SpectralPower[i] + 1e-12) / SumPowerEps;
        Desc.NormSpectralPower[i] = normSp;
        double Diff = normSp - Mean;
        Variance += Diff * Diff;

        const double freqDiff = (static_cast<double>(i) * binWidth) - centroid;
        weightedSpreadVariance += (freqDiff * freqDiff) * Desc.SpectralPower[i];
    }

    Desc.StdDev = std::sqrt(Variance / NHalf);
    Desc.SpectralSpread = hasSpectralEnergy ? std::sqrt(weightedSpreadVariance / SumPower) : 0.0;

    // Finalize descriptors derived from the fused pass.
    Desc.SpectralFlux = flux;
    Desc.SpectralIrregularity = irregularityDenominator > 0.0 ? (irregularityNumerator / irregularityDenominator) : 0.0;
    Desc.SpectralCrest = Desc.MaxAmp / ((SumPower / static_cast<double>(NHalf)) + 1e-12);
    Desc.SpectralFlatness =
        std::exp(logSumPower / static_cast<double>(NHalf)) / (linSumPower / static_cast<double>(NHalf));
    Desc.Harmonicity = harmonicityPeak / (harmonicitySum + 1e-12);
    Desc.HighFreqRatio = highFreqEnergy / (SumPower + 1e-12);

    const size_t prefixSize = NHalf + 1;
    if (m_SpectralPrefix.size() != prefixSize) {
        m_SpectralPrefix.resize(prefixSize);
    }

    m_SpectralPrefix[0] = 0.0;
    for (size_t i = 0; i < NHalf; ++i) {
        m_SpectralPrefix[i + 1] = m_SpectralPrefix[i] + Desc.SpectralPower[i];
    }

    // Remaining descriptors that require their own transforms.
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
    m_DCTBasis.assign(m_MFCC, std::vector<double>(m_MFCCMels));

    const double scale0 = std::sqrt(1.0 / m_MFCCMels);
    const double scale = std::sqrt(2.0 / m_MFCCMels);

    for (int k = 0; k < m_MFCC; ++k) {
        for (int n = 0; n < m_MFCCMels; ++n) {
            m_DCTBasis[k][n] = (k == 0 ? scale0 : scale) * std::cos(std::numbers::pi * (n + 0.5) * k / m_MFCCMels);
        }
    }

    m_MFCCEnergy.resize(m_MFCCMels);
}

// ─────────────────────────────────────
void MIR::MFCCExec(Description &Desc) {
    const int n_mels = std::min<int>(m_MFCCMels, static_cast<int>(m_MFCCFilter.size()));
    const int n_mfcc = std::min<int>(m_MFCC, static_cast<int>(m_DCTBasis.size()));

    if (n_mels <= 0 || n_mfcc <= 0 || Desc.Power.empty()) {
        Desc.MFCC.assign(std::max(0, n_mfcc), 0.0);
        return;
    }

    const int n_f = std::min<int>(static_cast<int>(Desc.Power.size()), static_cast<int>(m_MFCCFilter[0].size()));

    if (n_f <= 0) {
        Desc.MFCC.assign(n_mfcc, 0.0);
        return;
    }

    constexpr float kAmin = 1e-10f;
    constexpr float kTopDb = 80.0f;

    if (static_cast<int>(m_MFCCEnergy.size()) != n_mels)
        m_MFCCEnergy.resize(n_mels);

    // Mel projection (power domain)
    for (int m = 0; m < n_mels; ++m) {
        const double *filter = m_MFCCFilter[m].data();
        float mel = 0.0f;
        for (int k = 0; k < n_f; ++k) {
            mel += static_cast<float>(filter[k]) * static_cast<float>(Desc.Power[k]);
        }

        m_MFCCEnergy[m] = static_cast<double>(mel);
    }

    // power_to_db(ref=1.0, top_db=80)
    float maxLog = -std::numeric_limits<float>::infinity();

    for (int m = 0; m < n_mels; ++m) {
        const float v = 10.0f * std::log10(std::max(kAmin, static_cast<float>(m_MFCCEnergy[m])));

        m_MFCCEnergy[m] = static_cast<double>(v);
        if (v > maxLog)
            maxLog = v;
    }

    if (std::isfinite(maxLog)) {
        const float floor = maxLog - kTopDb;
        for (int m = 0; m < n_mels; ++m) {
            if (static_cast<float>(m_MFCCEnergy[m]) < floor)
                m_MFCCEnergy[m] = static_cast<double>(floor);
        }
    }

    Desc.MFCC.assign(n_mfcc, 0.0);

    // DCT-II
    for (int k = 0; k < n_mfcc; ++k) {
        const double *basis = m_DCTBasis[k].data();
        float coeff = 0.0f;

        for (int n = 0; n < n_mels; ++n) {
            coeff += static_cast<float>(basis[n]) * static_cast<float>(m_MFCCEnergy[n]);
        }

        Desc.MFCC[k] = static_cast<double>(coeff);
    }
}

// ╭─────────────────────────────────────╮
// │               Chroma                │
// ╰─────────────────────────────────────╯
double MIR::HzToOcts(double frequency, double tuning, int binsPerOctave) const {
    const double a440 = m_ChromaA440 * std::pow(2.0, tuning / static_cast<double>(binsPerOctave));
    return std::log2(frequency / (a440 / 16.0));
}

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
        for (size_t chroma = 0; chroma < m_ChromaSize; ++chroma) {
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
            for (size_t chroma = 0; chroma < m_ChromaSize; ++chroma) {
                m_ChromaFilter[chroma][k] *= invNorm;
            }
        }

        const double octaveWeight = std::exp(
            -0.5 * std::pow((frqbins[k] / static_cast<double>(m_ChromaSize) - m_ChromaCenterOctave) / m_ChromaOctaveWidth,
                            2.0));
        for (size_t chroma = 0; chroma < m_ChromaSize; ++chroma) {
            m_ChromaFilter[chroma][k] *= octaveWeight;
        }
    }

    const size_t chromaShift = 3 * (m_ChromaSize / 12);
    if (chromaShift > 0 && chromaShift < m_ChromaSize) {
        Matrix rolled(m_ChromaSize, std::vector<double>(nHalf, 0.0));
        for (size_t chroma = 0; chroma < m_ChromaSize; ++chroma) {
            rolled[chroma] = m_ChromaFilter[(chroma + chromaShift) % m_ChromaSize];
        }
        m_ChromaFilter.swap(rolled);
    }
}

// ─────────────────────────────────────
void MIR::SpectralChromaExec(Description &Desc) {
    if (Desc.Chroma.size() != m_ChromaSize)
        Desc.Chroma.resize(m_ChromaSize);

    std::fill(Desc.Chroma.begin(), Desc.Chroma.end(), 0.0);

    if (m_ChromaFilter.empty() || m_ChromaFilter[0].empty() || Desc.Power.empty()) {
        return;
    }

    const size_t nHalf = std::min(Desc.Power.size(), m_ChromaFilter[0].size());

    for (size_t chroma = 0; chroma < m_ChromaSize; ++chroma) {
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
    const auto &S = Desc.SpectralPower;
    const size_t N = S.size();

    if (m_PreviousSpectralPower.size() != N) {
        m_PreviousSpectralPower.assign(N, 0.0);
    }

    double flux = 0.0;
    for (size_t k = 1; k < N; ++k) {
        double diff = S[k] - m_PreviousSpectralPower[k];
        if (diff > 0.0) {
            flux += diff;
        }
        m_PreviousSpectralPower[k] = S[k];
    }
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
    const auto &S = Desc.SpectralPower;
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
    for (size_t i = 0; i < Desc.NormSpectralPower.size(); i++) {
        Desc.ReverbSpectralPower[i] = 0; //(Desc.ReverbSpectralPower[i] * decay) * (Desc.NormSpectralPower[i] * decay);
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
    const double *w = m_WindowingFunc.data();
    for (size_t i = 0; i < m_FFTSize; ++i) {
        m_FFTIn[i] = x[i] * w[i];
    }
    GetSpectralDescriptions(Desc);

    // Perc vs Pitch
    // Take the strongest percussive indicator: energy burst OR noise burst
    double inst_perc = Desc.SpectralFlux / (Desc.Harmonicity + 1e-6);

    inst_perc *= (1.0 - Desc.SilenceProb * Desc.Harmonicity);
    inst_perc = std::min(1.0, std::max(1e-300, inst_perc));

    Desc.PercussiveProb = inst_perc;
    m_PrevPercussiveProb = Desc.PercussiveProb;
}

} // namespace OpenScofo
