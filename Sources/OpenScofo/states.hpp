#pragma once

#include <string>
#include <vector>
#include <variant>
#include <span>
#include <unordered_map>

namespace OpenScofo {

enum AudioDescType {
    PITCH,
    SILENCE,
    LABEL,
};

// ─────────────────────────────────────
enum Descriptors {
    MFCC,
    LOUDNESS,
    RMS,
    POWER,
    CHROMA,
    ZCR,
    HFR,
    CENTROID,
    SPREAD,
    FLATNESS,
    FLUX,
    IRREGULARITY,
    HARMONICITY,
    YIN,

    //
    SILENCEPROB,
    PERCUSSIVEPROB,
    ONSET,
};

// ─────────────────────────────────────
enum EventType {
    // Tradicional
    REST,
    NOTE,
    CHORD,
    TRILL,
    MULTI,

    //  AI model
    EXTENDED,

    // Cort Lippe event (event is just count)
    EVENT,
};

enum HMMType { SEMIMARKOV, MARKOV };

// ─────────────────────────────────────
class Action {
  public:
    bool isLua;
    std::string Lua;
    std::string Receiver;
    std::vector<std::variant<float, int, std::string>> Args;
    bool AbsoluteTime;
    double Time;
};

using ActionVec = std::vector<Action>;

// ─────────────────────────────────────
class AudioState {
  public:
    AudioDescType Type;
    double Freq;
    double Midi;
    std::string Label;
    unsigned Index;
};

// ─────────────────────────────────────
class MarkovState {
  public:
    int Index;
    int ScorePos;
    int SubStateIndex = 0;
    EventType Type;
    HMMType HSMMType;
    int MarkovIndex = -1;

    // States Actions
    ActionVec Actions;
    std::vector<AudioState> AudioStates;

    // Forward Algorithm — Right-Censored Forward (Cuvillier 2014)
    // F_j(t) = sum_u obs_prod(u) * F_j^i(t-u) * D_j(u)
    // F_j^o(t) = sum_u obs_prod(u) * F_j^i(t-u) * d_j(u)   [d_j = D_j(u)-D_j(u+1)]
    // F_j^i(t) = sum_{i!=j} p_ij * F_i^o(t)                [= F_{j-1}^o(t) for linear chain]
    // obs_prod(u) = prod_{v=0}^{u-1} b_j(x_{t-v}) / N(t-v)  [normalized, stable]
    double InitProb;
    std::vector<double> Forward;  // F_j(t): current-state prob, circular buffer
    std::vector<double> ExitProb; // F_j^o(t): exit prob, circular buffer
    std::vector<double> BestObs;  // b_j(x_t): observation, circular buffer

    // Time
    int UpperBound;
    double BPMExpected;
    double BPMObserved = 0;
    double OnsetExpected = 0.0;
    double OnsetObserved = 0;
    double PhaseExpected;
    double PhaseObserved = 0;
    double IOIPhiN;
    double IOIHatPhiN;
    double Duration = 0.0;
    double PhaseCoupling;
    double SyncStrength;
    double TimeProb;

    // Error Handling
    int Line;
};

using States = std::vector<MarkovState>;

// ─────────────────────────────────────
class Description {
  public:
    bool Onset;

    // Probability
    double SilenceProb;
    double ExtendedTechProb;

    // Amplitude
    double dB;
    double RMS;
    double MaxAmp;
    double Loudness;

    // Spectral
    double Harmonicity;
    double SpectralFlatness;
    double SpectralFlux;
    double SpectralIrregularity = 0.0;
    double SpectralCrest = 0.0;
    double SpectralCentroid = 0.0;
    double CentroidVelocity = 0.0;
    double SpectralSpread = 0.0;
    double HighFreqRatio = 0.0;
    double Peakiness = 0.0;
    double ZeroCrossingRate;
    double StdDev;
    double Pitch = 0.0;
    double PitchConfidence = 0.0;
    std::vector<double> Power;
    std::vector<double> SpectralPower;
    std::vector<double> NormSpectralPower;
    std::vector<double> ReverbSpectralPower;
    std::vector<double> MFCC;
    std::vector<double> Chroma; // size 12

    // ONNX
    std::unordered_map<std::string, float> ONNX;
};

// ─────────────────────────────────────
struct MIRConfig {
    double FFTSize;
    double HopSize;
    double Chroma;
    double SilencedB;
};

// ─────────────────────────────────────
struct MDPConfig {
    double PitchSigma;
    double HopSize;
    double Chroma;
    double SilencedB;
};

} // namespace OpenScofo
