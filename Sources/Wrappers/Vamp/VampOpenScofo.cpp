#include "VampOpenScofo.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>

namespace {

enum class VampOutputKind {
    Scalar,
    Magnitude,
    Power,
    SpectralMagnitudeNorm,
    SpectralMagnitudeFrameNorm,
    ReverbSpectralPower,
    MFCC,
    LogMel,
    Chroma,
};

struct VampOutputSpec {
    const char *identifier;
    const char *name;
    const char *description;
    VampOutputKind kind;
    OpenScofo::Descriptors descriptor;
};

const std::vector<VampOutputSpec> &GetVampOutputSpecs() {
    static const std::vector<VampOutputSpec> specs = {
        // Spectral arrays
        {"magnitude", "Magnitude", "Spectral array: magnitude spectrum", VampOutputKind::Magnitude,
         OpenScofo::Descriptors::MAGNITUDE},
        {"power", "Power", "Spectral array: power spectrum", VampOutputKind::Power,
         OpenScofo::Descriptors::POWERARRAY},
        {"spectral_magnitude_norm", "Spectral Magnitude Norm", "Spectral array: bin-normalized magnitude",
         VampOutputKind::SpectralMagnitudeNorm, OpenScofo::Descriptors::INVALID},
        {"spectral_magnitude_frame_norm", "Spectral Magnitude Frame Norm",
         "Spectral array: frame-normalized magnitude", VampOutputKind::SpectralMagnitudeFrameNorm,
         OpenScofo::Descriptors::INVALID},
        {"reverb_spectral_power", "Reverb Spectral Power", "Spectral array: reverberant spectral power",
         VampOutputKind::ReverbSpectralPower, OpenScofo::Descriptors::INVALID},

        // Compact spectral arrays
        {"mfcc", "MFCC", "Spectral array: mel-frequency cepstral coefficients", VampOutputKind::MFCC,
         OpenScofo::Descriptors::MFCC},
        {"logmel", "Log Mel Spectrum", "Spectral array: log mel spectrum", VampOutputKind::LogMel,
         OpenScofo::Descriptors::MELOGRAM},
        {"chroma", "Chroma", "Spectral array: chroma energy", VampOutputKind::Chroma,
         OpenScofo::Descriptors::CHROMA},

        // Event/probability
        {"onset", "Onset", "Event/probability: onset detection function", VampOutputKind::Scalar,
         OpenScofo::Descriptors::ODSONSET},
        {"silence", "Silence Probability", "Event/probability: silence probability", VampOutputKind::Scalar,
         OpenScofo::Descriptors::SILENCEPROB},
        {"ext", "Extended Technique Probability", "Event/probability: extended technique probability",
         VampOutputKind::Scalar, OpenScofo::Descriptors::EXTENDEDTECHNIQUE},

        // Amplitude
        {"maxamp", "Max Amplitude", "Amplitude: maximum normalized amplitude", VampOutputKind::Scalar,
         OpenScofo::Descriptors::MAXAMP},
        {"loudness", "Loudness", "Amplitude: loudness", VampOutputKind::Scalar, OpenScofo::Descriptors::LOUDNESS},
        {"db", "dB", "Amplitude: decibels", VampOutputKind::Scalar, OpenScofo::Descriptors::DB},
        {"rms", "RMS", "Amplitude: root mean square", VampOutputKind::Scalar, OpenScofo::Descriptors::RMS},
        {"stddev", "Amplitude Standard Deviation", "Amplitude: standard deviation", VampOutputKind::Scalar,
         OpenScofo::Descriptors::STDDEV},

        // Spectral scalars
        {"flux", "Spectral Flux", "Spectral scalar: flux", VampOutputKind::Scalar, OpenScofo::Descriptors::FLUX},
        {"irregularity", "Spectral Irregularity", "Spectral scalar: irregularity", VampOutputKind::Scalar,
         OpenScofo::Descriptors::IRREGULARITY},
        {"crest", "Spectral Crest", "Spectral scalar: crest", VampOutputKind::Scalar,
         OpenScofo::Descriptors::CREST},
        {"centroid", "Spectral Centroid", "Spectral scalar: centroid", VampOutputKind::Scalar,
         OpenScofo::Descriptors::CENTROID},
        {"centroid_velocity", "Centroid Velocity", "Spectral scalar: centroid velocity", VampOutputKind::Scalar,
         OpenScofo::Descriptors::CENTROIDVEL},
        {"spreadhz", "Spectral Spread Hz", "Spectral scalar: spread in Hz", VampOutputKind::Scalar,
         OpenScofo::Descriptors::SPREADHZ},
        {"spread_variance", "Spectral Spread Variance", "Spectral scalar: normalized spread variance",
         VampOutputKind::Scalar, OpenScofo::Descriptors::SPREADVARIANCE},
        {"flatness", "Spectral Flatness", "Spectral scalar: flatness", VampOutputKind::Scalar,
         OpenScofo::Descriptors::FLATNESS},
        {"entropy", "Spectral Entropy", "Spectral scalar: entropy", VampOutputKind::Scalar,
         OpenScofo::Descriptors::ENTROPY},
        {"rolloff", "Spectral Rolloff", "Spectral scalar: rolloff frequency", VampOutputKind::Scalar,
         OpenScofo::Descriptors::ROLLOFF},
        {"hfr", "High Frequency Ratio", "Spectral scalar: high-frequency ratio", VampOutputKind::Scalar,
         OpenScofo::Descriptors::HFR},
        {"harmonicity", "Harmonicity", "Spectral scalar: harmonicity", VampOutputKind::Scalar,
         OpenScofo::Descriptors::HARMONICITY},
        {"zcr", "Zero Crossing Rate", "Spectral scalar: zero crossing rate", VampOutputKind::Scalar,
         OpenScofo::Descriptors::ZCR},
        {"skewness", "Spectral Skewness", "Spectral scalar: skewness", VampOutputKind::Scalar,
         OpenScofo::Descriptors::SKEWNESS},
        {"slope", "Spectral Slope", "Spectral scalar: slope", VampOutputKind::Scalar,
         OpenScofo::Descriptors::SLOPE},
        {"kurtosis", "Spectral Kurtosis", "Spectral scalar: kurtosis", VampOutputKind::Scalar,
         OpenScofo::Descriptors::KURTOSIS},

        // Pitch
        {"yin", "Pitch", "Pitch: YIN frequency estimate", VampOutputKind::Scalar, OpenScofo::Descriptors::YIN},
        {"yin_confidence", "Pitch Confidence", "Pitch: YIN confidence", VampOutputKind::Scalar,
         OpenScofo::Descriptors::YINCONFIDENCE},
    };

    return specs;
}

bool IsVectorOutput(VampOutputKind kind) {
    return kind != VampOutputKind::Scalar;
}

size_t GetSpectralBinCount(size_t blockSize) {
    return blockSize / 2 + 1;
}

size_t GetBinCount(const VampOutputSpec &spec, const OpenScofo::Configuration &config, size_t blockSize) {
    switch (spec.kind) {
    case VampOutputKind::Magnitude:
    case VampOutputKind::Power:
    case VampOutputKind::SpectralMagnitudeNorm:
    case VampOutputKind::SpectralMagnitudeFrameNorm:
    case VampOutputKind::ReverbSpectralPower:
        return GetSpectralBinCount(blockSize > 0 ? blockSize : static_cast<size_t>(config.FFTSize));
    case VampOutputKind::MFCC:
        return static_cast<size_t>(config.MFCCCount);
    case VampOutputKind::LogMel:
        return static_cast<size_t>(config.MFCCMels);
    case VampOutputKind::Chroma:
        return static_cast<size_t>(config.ChromaSize);
    case VampOutputKind::Scalar:
        return 1;
    }

    return 1;
}

const std::vector<double> &GetVectorValue(const OpenScofo::Description &desc, VampOutputKind kind) {
    switch (kind) {
    case VampOutputKind::Magnitude:
        return desc.Magnitude;
    case VampOutputKind::Power:
        return desc.Power;
    case VampOutputKind::SpectralMagnitudeNorm:
        return desc.SpectralMagnitudeNorm;
    case VampOutputKind::SpectralMagnitudeFrameNorm:
        return desc.SpectralMagnitudeFrameNorm;
    case VampOutputKind::ReverbSpectralPower:
        return desc.ReverbSpectralPower;
    case VampOutputKind::MFCC:
        return desc.MFCC;
    case VampOutputKind::LogMel:
        return desc.LogMelSpectrum;
    case VampOutputKind::Chroma:
        return desc.Chroma;
    case VampOutputKind::Scalar:
        break;
    }

    static const std::vector<double> empty;
    return empty;
}

} // namespace

VampOpenScofo::VampOpenScofo(float inputSampleRate)
    : Plugin(inputSampleRate), m_blockSize(0), m_stepSize(0), m_selectedScoreIndex(0), m_OpenScofo(nullptr) {

    //
}

// ─────────────────────────────────────
VampOpenScofo::~VampOpenScofo() {
    delete m_OpenScofo;
}

// ─────────────────────────────────────
std::string VampOpenScofo::getIdentifier() const {
    return "openscorefollower";
}

// ─────────────────────────────────────
std::string VampOpenScofo::getName() const {
    return "OpenSource Score Follower";
}

// ─────────────────────────────────────
std::string VampOpenScofo::getDescription() const {
    return "A Machine Listening System for Contemporary Music (just descriptors work on vamp plugins)";
}

// ─────────────────────────────────────
std::string VampOpenScofo::getType() const {
    return "Score Following";
}

// ─────────────────────────────────────
std::string VampOpenScofo::getMaker() const {
    return "Charles K. Neimog";
}

// ─────────────────────────────────────
int VampOpenScofo::getPluginVersion() const {
    return 1;
}

// ─────────────────────────────────────
std::string VampOpenScofo::getCopyright() const {
    return "GPL3";
}

// ─────────────────────────────────────
Vamp::Plugin::InputDomain VampOpenScofo::getInputDomain() const {
    return TimeDomain;
}

// ─────────────────────────────────────
size_t VampOpenScofo::getPreferredBlockSize() const {
    return 2048;
}

// ─────────────────────────────────────
size_t VampOpenScofo::getPreferredStepSize() const {
    return 512;
}

// ─────────────────────────────────────
size_t VampOpenScofo::getMinChannelCount() const {
    return 1;
}

// ─────────────────────────────────────
size_t VampOpenScofo::getMaxChannelCount() const {
    return 1;
}

// ─────────────────────────────────────
Vamp::Plugin::ParameterList VampOpenScofo::getParameterDescriptors() const {
    ParameterList list;

    ParameterDescriptor scorePath;
    scorePath.identifier = "score_path";
    scorePath.name = "Score Path";
    scorePath.description = "Select score index from Programs, or set OPENSCOFO_VAMP_SCORE_PATH to a full path.";
    scorePath.unit = "index";
    scorePath.minValue = 0;
    scorePath.maxValue = m_scorePaths.empty() ? 0 : float(m_scorePaths.size() - 1);
    scorePath.defaultValue = 0;
    scorePath.isQuantized = true;
    scorePath.quantizeStep = 1;
    list.push_back(scorePath);

    return list;
}

// ─────────────────────────────────────
float VampOpenScofo::getParameter(std::string identifier) const {
    if (identifier == "score_path") {
        return float(m_selectedScoreIndex);
    }

    return 0.f;
}

// ─────────────────────────────────────
void VampOpenScofo::setParameter(std::string identifier, float value) {
    if (identifier != "score_path")
        return;

    if (m_scorePaths.empty()) {
        m_selectedScoreIndex = 0;
        return;
    }

    int idx = int(std::lround(value));
    idx = std::max(0, std::min(idx, int(m_scorePaths.size() - 1)));
    m_selectedScoreIndex = idx;

    if (m_OpenScofo) {
        LoadScoreAtIndex(m_selectedScoreIndex);
    }
}

// ─────────────────────────────────────
Vamp::Plugin::ProgramList VampOpenScofo::getPrograms() const {
    ProgramList list;
    for (const std::string &path : m_scorePaths) {
        list.push_back(path);
    }
    return list;
}

// ─────────────────────────────────────
std::string VampOpenScofo::getCurrentProgram() const {
    return m_currentProgram;
}

// ─────────────────────────────────────
void VampOpenScofo::selectProgram(std::string name) {
    if (name.empty())
        return;

    auto it = std::find(m_scorePaths.begin(), m_scorePaths.end(), name);
    if (it == m_scorePaths.end()) {
        m_scorePaths.push_back(name);
        m_selectedScoreIndex = int(m_scorePaths.size() - 1);
    } else {
        m_selectedScoreIndex = int(std::distance(m_scorePaths.begin(), it));
    }

    if (m_OpenScofo && m_OpenScofo->LoadScore(name)) {
        m_currentProgram = name;
    }
}

// ─────────────────────────────────────
Vamp::Plugin::FeatureSet VampOpenScofo::getRemainingFeatures() {
    return {};
}

// ─────────────────────────────────────
bool VampOpenScofo::initialise(size_t channels, size_t stepSize, size_t blockSize) {
    if (channels != 1) {
        return false;
    }

    m_blockSize = blockSize;
    m_stepSize = stepSize;

    delete m_OpenScofo;
    m_OpenScofo = new OpenScofo::OpenScofo(m_inputSampleRate, m_blockSize, m_blockSize);

    if (!m_scorePaths.empty()) {
        LoadScoreAtIndex(m_selectedScoreIndex);
    }

    return true;
}

// ─────────────────────────────────────
void VampOpenScofo::reset() {
}

// ─────────────────────────────────────
void VampOpenScofo::RefreshScorePathsFromEnv() {
}

// ─────────────────────────────────────
bool VampOpenScofo::LoadScoreAtIndex(int index) {
    if (!m_OpenScofo)
        return false;

    if (index < 0 || index >= int(m_scorePaths.size()))
        return false;

    const std::string &scorePath = m_scorePaths[size_t(index)];
    if (!m_OpenScofo->LoadScore(scorePath))
        return false;

    m_currentProgram = scorePath;
    return true;
}

// ─────────────────────────────────────
Vamp::Plugin::FeatureSet VampOpenScofo::process(const float *const *inputBuffers, Vamp::RealTime) {
    FeatureSet fs;

    bool ok = m_OpenScofo->ProcessBlock(inputBuffers[0], m_blockSize);
    if (!ok) {
    }
    OpenScofo::Description Desc = m_OpenScofo->GetDescription();

    auto pushScalar = [&](int idx, float value) {
        Feature f;
        f.values.push_back(value);
        fs[idx].push_back(f);
    };
    auto pushVector = [&](int idx, const std::vector<double> &vec) {
        Feature f;
        f.values.reserve(vec.size());
        for (double v : vec) {
            f.values.push_back(static_cast<float>(v));
        }

        fs[idx].push_back(f);
    };

    int i = 0;
    for (const auto &spec : GetVampOutputSpecs()) {
        if (IsVectorOutput(spec.kind)) {
            pushVector(i++, GetVectorValue(Desc, spec.kind));
        } else {
            pushScalar(i++, static_cast<float>(m_OpenScofo->GetDescriptionFloat(Desc, spec.descriptor)));
        }
    }

    return fs;
}

// ─────────────────────────────────────
Vamp::Plugin::OutputList VampOpenScofo::getOutputDescriptors() const {
    OutputList list;

    OpenScofo::Configuration config;
    if (m_OpenScofo) {
        config = m_OpenScofo->GetConfiguration();
    }

    auto makeDescriptor = [&](const VampOutputSpec &spec) {
        OutputDescriptor d;
        d.identifier = spec.identifier;
        d.name = spec.name;
        d.description = spec.description;
        d.hasFixedBinCount = true;
        d.binCount = GetBinCount(spec, config, m_blockSize > 0 ? m_blockSize : getPreferredBlockSize());
        d.sampleType = OutputDescriptor::OneSamplePerStep;
        return d;
    };

    for (const auto &spec : GetVampOutputSpecs()) {
        list.push_back(makeDescriptor(spec));
    }

    return list;
}

static Vamp::PluginAdapter<VampOpenScofo> adapter;

// Windows hosts often fail to discover descriptors unless the entry point is
// explicitly exported from the DLL.
#if defined(_WIN32)
#define VAMP_PLUGIN_EXPORT extern "C" __declspec(dllexport)
#else
#define VAMP_PLUGIN_EXPORT extern "C"
#endif

VAMP_PLUGIN_EXPORT const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int version, unsigned int index) {
    if (version < 1)
        return nullptr;

    if (index > 0)
        return nullptr;

    return adapter.getDescriptor();
}
