#include "mir.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <math.h>
#include <vector>
#include <numeric>
#include <complex>

namespace OScofo {

// ╭─────────────────────────────────────╮
// │Constructor and Destructor Functions │
// ╰─────────────────────────────────────╯
MIR::MIR(float Sr, float FftSize, float HopSize) {
    LOGE() << "OScofoMIR::OScofoMIR";

    m_HopSize = HopSize;
    m_FFTSize = FftSize;
    m_Sr = Sr;

    float WindowHalf = FftSize / 2;
    m_FFTIn = (double *)fftw_alloc_real((size_t)FftSize);
    if (!m_FFTIn) {
        SetError("OScofoMIR::OScofoMIR fftw_alloc_real failed");
        return;
    }

    m_FFTOut = (fftw_complex *)fftw_alloc_complex((size_t)WindowHalf + 1);
    if (!m_FFTOut) {
        fftw_free(m_FFTIn); // Free previously allocated memory
        SetError("OScofoMIR::OScofoMIR fftw_alloc_complex failed");
        return;
    }

#ifdef __EMSCRIPTEN__
    m_FFTPlan = fftw_plan_dft_r2c_1d((int)m_FftSize, m_FFTIn, m_FFTOut, FFTW_ESTIMATE);
#else
    m_FFTPlan = fftw_plan_dft_r2c_1d((int)m_FFTSize, m_FFTIn, m_FFTOut, FFTW_MEASURE);
#endif

    // hanning
    m_WindowingFunc.resize((size_t)m_FFTSize);
    for (size_t i = 0; i < (size_t)m_FFTSize; i++) {
        m_WindowingFunc[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (m_FFTSize - 1)));
    }

    // Init Onset
    OnsetInit();
    MFCCInit();
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
    }
}

// ╭─────────────────────────────────────╮
// │               Errors                │
// ╰─────────────────────────────────────╯
bool MIR::HasErrors() {
    return m_HasErrors;
}

// ─────────────────────────────────────
std::vector<std::string> MIR::GetErrorMessage() {
    return m_Errors;
}

// ─────────────────────────────────────
void MIR::SetError(const std::string &message) {
    m_HasErrors = true;
    m_Errors.push_back(message);
}

// ─────────────────────────────────────
void MIR::ClearError() {
    m_HasErrors = false;
    m_Errors.clear();
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

// ╭─────────────────────────────────────╮
// │          Set|Get Functions          │
// ╰─────────────────────────────────────╯
void MIR::SetdBTreshold(double dB) {
    m_dBTreshold = dB;
}

// ╭─────────────────────────────────────╮
// │          Machine Learning           │
// ╰─────────────────────────────────────╯
void MIR::LoadONNXModel(fs::path path) {
    m_OnnxCTX = onnx_context_alloc_from_file(path.c_str(), NULL, 0);
    if (m_OnnxCTX == nullptr) {
        SetError("Failed to load ONNX model: " + path.string());
        return;
    }

    if (m_OnnxCTX && m_OnnxCTX->g != nullptr) {
        struct onnx_graph_t *g = m_OnnxCTX->g;
        for (int i = 0; i < g->nlen; i++) {
            struct onnx_node_t *n = &g->nodes[i];
            if (n->opset > CURRENT_ONNX_OPSET) {
                SetError("Unsupported opset => " + std::string(n->proto->op_type) + " " + std::to_string(n->opset));
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
        SetError("TreeEnsembleClassifier not found in model, please use py.train to train models for OScofo");
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
    m_OnsetInit = true;
}

// ─────────────────────────────────────
void MIR::OnsetExec(Description &Desc) {
    if (!m_OnsetInit) {
        return;
    }

    const size_t numBins = m_OnsetFFTSize / 2 + 1;
    std::vector<float> fft_data(2 * numBins);

    for (size_t i = 0; i < numBins; ++i) {
        fft_data[2 * i] = m_FFTOut[i][0];
        fft_data[2 * i + 1] = m_FFTOut[i][1];
    }
    Desc.Onset = onsetsds_process(m_ODS, fft_data.data());
}

// ╭─────────────────────────────────────╮
// │           Time Coherence            │
// ╰─────────────────────────────────────╯
void MIR::BuildTimeCoherenceTemplate(States &ScoreStates) {
    const double dt = static_cast<double>(m_HopSize) / m_Sr;

    double totalTime = 0.0;
    for (const auto &s : ScoreStates) {
        totalTime += s.Duration;
    }

    const size_t N = static_cast<size_t>(std::ceil(totalTime / dt));
    std::vector<double> probability(N, 0.0);

    // Gaussian width (seconds)
    const double sigma = 0.005;
    const double inv2Sigma2 = 1.0 / (2.0 * sigma * sigma);

    // 2. event onsets
    double onset = 0.0;
    for (const auto &state : ScoreStates) {
        if (state.Type == REST) {
            continue;
        }
        const size_t center = static_cast<size_t>(onset / dt);

        // limit influence window to ±3σ
        const int radius = static_cast<int>(std::ceil(3.0 * sigma / dt));

        for (int i = -radius; i <= radius; ++i) {
            const int idx = static_cast<int>(center) + i;
            if (idx < 0 || idx >= static_cast<int>(N))
                continue;

            const double t = (idx * dt) - onset;
            const double g = std::exp(-(t * t) * inv2Sigma2);

            probability[idx] += g;
        }

        onset += state.Duration;
    }

    for (auto &p : probability)
        p = std::min(1.0, p);
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

    // Silence detection
    Desc.Silence = Desc.dB < m_dBTreshold;

    // Loudness (based on sum of squares)
    double meanSquare = z_loudness / In.size();
    if (meanSquare <= 0.0) {
        Desc.Loudness = -100.0; // substitui -inf
    } else {
        Desc.Loudness = -0.691 + 10.0 * std::log10(meanSquare);
    }

    // Compute silence probability
    const double L0 = -70.0;
    const double alpha = 0.25;
    Desc.SilenceProb = 1.0 / (1.0 + std::exp(alpha * (Desc.Loudness - L0)));
}

// ─────────────────────────────────────
void MIR::GetFFTDescriptions(std::vector<double> &In, Description &Desc) {
    // real audio analisys
    size_t N = In.size();
    size_t NHalf = N / 2 + 1; // include Nyquist to match rfft/librosa

    if (NHalf != Desc.Power.size()) {
        Desc.Power.resize(NHalf);
        Desc.SpectralPower.resize(NHalf);
        Desc.NormSpectralPower.resize(NHalf);
    }

    std::copy(In.begin(), In.end(), m_FFTIn);
    fftw_execute(m_FFTPlan);

    // FFT Mag
    Desc.MaxAmp = 0;
    const double invN = 1.0 / static_cast<double>(N);
    for (size_t i = 0; i < NHalf; ++i) {
        const double re = m_FFTOut[i][0];
        const double im = m_FFTOut[i][1];
        const double p = re * re + im * im;
        const double sp = std::sqrt(p) * invN;

        Desc.Power[i] = p;
        Desc.SpectralPower[i] = sp;
        Desc.MaxAmp = std::max(Desc.MaxAmp, sp);
    }

    OnsetExec(Desc);
    MFCCExec(Desc);

    // Normalize Spectral Power
    double SumPower = std::accumulate(Desc.SpectralPower.begin(), Desc.SpectralPower.end(), 0.0);

    for (size_t i = 0; i < NHalf; i++) {
        Desc.NormSpectralPower[i] = (Desc.SpectralPower[i] + 1e-12) / (SumPower + 1e-12);
    }

    const double Mean = 1.0 / NHalf; // Sum is 1.0, so mean = 1/N
    double Variance = 0.0;
    for (size_t i = 0; i < NHalf; i++) {
        double Diff = Desc.NormSpectralPower[i] - Mean;
        Variance += Diff * Diff;
    }
    Variance /= NHalf;
    Desc.StdDev = std::sqrt(Variance);
}

// ╭─────────────────────────────────────╮
// │                MFCC                 │
// ╰─────────────────────────────────────╯
void MIR::MFCCInit() {
    const int n_f = m_FFTSize / 2 + 1; // match librosa rfft length
    const double fmin = 0.0;
    const double fmax = m_Sr * 0.5;

    // FFT bin frequencies
    std::vector<double> fft_freqs(n_f);
    for (int i = 0; i < n_f; ++i) {
        fft_freqs[i] = i * m_Sr / static_cast<double>(m_FFTSize);
    }

    // Slaney-style mel scale used by librosa (htk=False)
    const double f_sp = 200.0 / 3.0;
    const double min_log_hz = 1000.0;
    const double min_log_mel = (min_log_hz - fmin) / f_sp;
    const double logstep = std::log(6.4) / 27.0;

    auto hz_to_mel = [&](double hz) {
        if (hz < min_log_hz) {
            return (hz - fmin) / f_sp;
        }
        return min_log_mel + std::log(hz / min_log_hz) / logstep;
    };

    auto mel_to_hz = [&](double mel) {
        if (mel < min_log_mel) {
            return mel * f_sp + fmin;
        }
        return min_log_hz * std::exp((mel - min_log_mel) * logstep);
    };

    const double mel_min = hz_to_mel(fmin);
    const double mel_max = hz_to_mel(fmax);

    std::vector<double> mel_pts(m_MFCCMels + 2);
    for (int i = 0; i < m_MFCCMels + 2; i++) {
        mel_pts[i] = mel_min + (mel_max - mel_min) * i / (m_MFCCMels + 1);
    }

    std::vector<double> hz_pts(m_MFCCMels + 2);
    for (int i = 0; i < m_MFCCMels + 2; i++) {
        hz_pts[i] = mel_to_hz(mel_pts[i]);
    }

    // Triangular Slaney mel filters normalized to match librosa
    m_MFCCFilter.assign(m_MFCCMels, std::vector<double>(n_f, 0.0));
    for (int m = 0; m < m_MFCCMels; m++) {
        const double f_left = hz_pts[m];
        const double f_center = hz_pts[m + 1];
        const double f_right = hz_pts[m + 2];
        const double norm = 2.0 / (f_right - f_left);

        for (int k = 0; k < n_f; k++) {
            const double f = fft_freqs[k];
            double w = 0.0;
            if (f >= f_left && f <= f_center) {
                w = (f - f_left) / (f_center - f_left);
            } else if (f > f_center && f <= f_right) {
                w = (f_right - f) / (f_right - f_center);
            }
            m_MFCCFilter[m][k] = std::max(0.0, w) * norm;
        }
    }

    // Orthonormal DCT-II basis (librosa dct_type=2, norm="ortho")
    m_DCTBasis.assign(m_MFCC, std::vector<double>(m_MFCCMels));
    const double scale0 = std::sqrt(1.0 / m_MFCCMels);
    const double scale = std::sqrt(2.0 / m_MFCCMels);
    for (int k = 0; k < m_MFCC; k++) {
        for (int n = 0; n < m_MFCCMels; n++) {
            m_DCTBasis[k][n] = (k == 0 ? scale0 : scale) * std::cos(M_PI * (n + 0.5) * k / m_MFCCMels);
        }
    }

    m_MFCCEnergy.resize(m_MFCCMels);
}

// ─────────────────────────────────────
void MIR::MFCCExec(Description &Desc) {
    const int n_f = m_FFTSize / 2 + 1;
    const int n_mels = m_MFCCMels;
    const int n_mfcc = m_MFCC;

    // Mel energies from power spectrum (librosa: np.maximum(..., 1e-30))
    for (int m = 0; m < n_mels; m++) {
        double sum = 0.0;
        for (int k = 0; k < n_f; k++) {
            sum += m_MFCCFilter[m][k] * Desc.Power[k];
        }
        m_MFCCEnergy[m] = std::max(sum, 1e-30);
    }

    // Power to dB and dynamic range compression (80 dB)
    double max_log = -1e30;
    for (int m = 0; m < n_mels; m++) {
        const double v = 10.0 * std::log10(m_MFCCEnergy[m]);
        m_MFCCEnergy[m] = v;
        if (v > max_log) {
            max_log = v;
        }
    }

    const double floor = max_log - 80.0;
    for (int m = 0; m < n_mels; m++) {
        if (m_MFCCEnergy[m] < floor) {
            m_MFCCEnergy[m] = floor;
        }
    }

    if (static_cast<int>(Desc.MFCC.size()) != n_mfcc) {
        Desc.MFCC.resize(n_mfcc);
    }

    // DCT-II with ortho normalization (librosa dct_type=2)
    for (int k = 0; k < n_mfcc; k++) {
        double sum = 0.0;
        for (int n = 0; n < n_mels; n++) {
            sum += m_DCTBasis[k][n] * m_MFCCEnergy[n];
        }
        Desc.MFCC[k] = sum;
    }
}

// ╭─────────────────────────────────────╮
// │            Main Function            │
// ╰─────────────────────────────────────╯
void MIR::GetDescription(std::vector<double> &In, Description &Desc) {
    // apply windowing function
    for (size_t i = 0; i < (size_t)m_FFTSize; i++) {
        In[i] *= m_WindowingFunc[i];
    }

    GetSignalPower(In, Desc);
    GetFFTDescriptions(In, Desc);
}
} // namespace OScofo
