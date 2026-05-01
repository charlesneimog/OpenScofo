#include "VampOpenScofo.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>

namespace {

std::string trim(const std::string &text) {
    size_t start = 0;
    while (start < text.size() && std::isspace(static_cast<unsigned char>(text[start]))) {
        start++;
    }

    size_t end = text.size();
    while (end > start && std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        end--;
    }

    return text.substr(start, end - start);
}

std::vector<std::string> splitSemicolonList(const std::string &text) {
    std::vector<std::string> out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ';')) {
        std::string cleaned = trim(item);
        if (!cleaned.empty()) {
            out.push_back(cleaned);
        }
    }
    return out;
}

} // namespace

// ─────────────────────────────────────
VampOpenScofo::VampOpenScofo(float inputSampleRate)
    : Plugin(inputSampleRate), m_blockSize(0), m_stepSize(0), m_selectedScoreIndex(0), m_OScofo(nullptr) {
    RefreshScorePathsFromEnv();
}

// ─────────────────────────────────────
VampOpenScofo::~VampOpenScofo() {
    delete m_OScofo;
}

// ─────────────────────────────────────
std::string VampOpenScofo::getIdentifier() const {
    return "openscorefollower";
}
std::string VampOpenScofo::getName() const {
    return "Open Score Follower";
}
std::string VampOpenScofo::getDescription() const {
    return "OpenScofo Score Follower";
}
std::string VampOpenScofo::getType() const {
    return "Score Following";
}
std::string VampOpenScofo::getMaker() const {
    return "Charles K. Neimog";
}
int VampOpenScofo::getPluginVersion() const {
    return 1;
}

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

float VampOpenScofo::getParameter(std::string identifier) const {
    if (identifier == "score_path") {
        return float(m_selectedScoreIndex);
    }

    return 0.f;
}

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

    if (m_OScofo) {
        LoadScoreAtIndex(m_selectedScoreIndex);
    }
}

Vamp::Plugin::ProgramList VampOpenScofo::getPrograms() const {
    ProgramList list;
    for (const std::string &path : m_scorePaths) {
        list.push_back(path);
    }
    return list;
}

std::string VampOpenScofo::getCurrentProgram() const {
    return m_currentProgram;
}

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

    if (m_OScofo && m_OScofo->LoadScore(name)) {
        m_currentProgram = name;
    }
}

// ─────────────────────────────────────
Vamp::Plugin::FeatureSet VampOpenScofo::getRemainingFeatures() {
    return {};
}

// ─────────────────────────────────────
bool VampOpenScofo::initialise(size_t channels, size_t stepSize, size_t blockSize) {
    if (channels != 1)
        return false;

    m_blockSize = blockSize;
    m_stepSize = stepSize;

    delete m_OScofo;
    m_OScofo = new OpenScofo::OpenScofo(m_inputSampleRate, m_blockSize, m_blockSize);

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
    m_scorePaths.clear();

    const char *singlePath = std::getenv("OPENSCOFO_VAMP_SCORE_PATH");
    if (singlePath && *singlePath) {
        std::string path = trim(std::string(singlePath));
        if (!path.empty()) {
            m_scorePaths.push_back(path);
            m_selectedScoreIndex = 0;
            return;
        }
    }

    const char *env = std::getenv("OPENSCOFO_VAMP_SCORE_PATHS");
    if (!env || !*env) {
        m_selectedScoreIndex = 0;
        return;
    }

    m_scorePaths = splitSemicolonList(env);
    if (m_scorePaths.empty() || m_selectedScoreIndex >= int(m_scorePaths.size())) {
        m_selectedScoreIndex = 0;
    }
}

// ─────────────────────────────────────
bool VampOpenScofo::LoadScoreAtIndex(int index) {
    if (!m_OScofo)
        return false;

    if (index < 0 || index >= int(m_scorePaths.size()))
        return false;

    const std::string &scorePath = m_scorePaths[size_t(index)];
    if (!m_OScofo->LoadScore(scorePath))
        return false;

    m_currentProgram = scorePath;
    return true;
}

// ─────────────────────────────────────
Vamp::Plugin::FeatureSet VampOpenScofo::process(const float *const *inputBuffers, Vamp::RealTime) {
    FeatureSet fs;

    bool ok = m_OScofo->ProcessBlock(inputBuffers[0], m_blockSize);
    if (!ok) {
    }
    OpenScofo::Description Desc = m_OScofo->GetDescription();

    // lambda
    auto pushScalar = [&](int idx, double value) {
        Feature f;
        f.values.push_back((float)value);
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

    // vectors
    pushVector(i++, Desc.MFCC);
    pushVector(i++, Desc.LogMelSpectrum);
    pushVector(i++, Desc.Chroma);

    // scalars
    pushScalar(i++, Desc.Onset);
    pushScalar(i++, Desc.SilenceProb);
    pushScalar(i++, Desc.ExtendedTechProb);

    pushScalar(i++, Desc.MaxAmp);
    pushScalar(i++, Desc.Loudness);
    pushScalar(i++, Desc.SpectralFlux);
    pushScalar(i++, Desc.SpectralIrregularity);
    pushScalar(i++, Desc.SpectralCrest);
    pushScalar(i++, Desc.SpectralCentroid);
    pushScalar(i++, Desc.SpectralSpreadHz);
    pushScalar(i++, Desc.SpectralSpreadVariance);
    pushScalar(i++, Desc.SpectralFlatness);

    pushScalar(i++, Desc.CentroidVelocity);

    pushScalar(i++, Desc.HighFreqRatio);
    pushScalar(i++, Desc.Harmonicity);
    pushScalar(i++, Desc.ZeroCrossingRate);
    pushScalar(i++, Desc.StdDev);
    pushScalar(i++, Desc.Pitch);
    pushScalar(i++, Desc.PitchConfidence);

    pushScalar(i++, Desc.dB);
    pushScalar(i++, Desc.RMS);

    return fs;
}

// ─────────────────────────────────────
Vamp::Plugin::OutputList VampOpenScofo::getOutputDescriptors() const {
    OutputList list;

    auto makeScalar = [](const std::string &id, const std::string &name) {
        OutputDescriptor d;
        d.identifier = id;
        d.name = name;
        d.hasFixedBinCount = true;
        d.binCount = 1;
        d.sampleType = OutputDescriptor::OneSamplePerStep;
        return d;
    };

    auto makeVector = [](const std::string &id, const std::string &name, size_t bins) {
        OutputDescriptor d;
        d.identifier = id;
        d.name = name;
        d.hasFixedBinCount = true;
        d.binCount = bins;
        d.sampleType = OutputDescriptor::OneSamplePerStep;
        return d;
    };

    // ---- vectors (set fixed sizes!) ----
    list.push_back(makeVector("mfcc", "MFCC", 13));               // adjust if needed
    list.push_back(makeVector("melogram", "LogMelSpectrum", 40)); // adjust if needed
    list.push_back(makeVector("chroma", "Chroma", 12));           // standard

    // ---- scalars ----
    list.push_back(makeScalar("onset", "Onset"));
    list.push_back(makeScalar("silence", "Silence"));
    list.push_back(makeScalar("ext", "Extended Technique Probability"));

    list.push_back(makeScalar("max_amp", "Max Amplitude"));
    list.push_back(makeScalar("loudness", "Loudness"));
    list.push_back(makeScalar("flux", "Spectral Flux"));
    list.push_back(makeScalar("irregularity", "Spectral Irregularity"));
    list.push_back(makeScalar("crest", "Spectral Crest"));
    list.push_back(makeScalar("centroid", "Spectral Centroid"));
    list.push_back(makeScalar("spreadhz", "Spectral Spread Hz (librosa)"));
    list.push_back(makeScalar("spread_variance", "Spectral Spread Variance (Essentia)"));
    list.push_back(makeScalar("flatness", "Spectral Flatness"));

    list.push_back(makeScalar("centroid_velocity", "Centroid Velocity"));

    list.push_back(makeScalar("hfr", "High Frequency Ratio"));
    list.push_back(makeScalar("harmonicity", "Harmonicity"));
    list.push_back(makeScalar("zcr", "Zero Crossing Rate"));
    list.push_back(makeScalar("stddev", "Ampltiude Standard Deviation"));
    list.push_back(makeScalar("pitch", "Pitch"));
    list.push_back(makeScalar("pitch_confidence", "Pitch Confidence"));

    list.push_back(makeScalar("db", "dB"));
    list.push_back(makeScalar("rms", "RMS"));

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
