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
    INVALID,

    ONSET,

    // Amplitude
    LOUDNESS,
    DB,
    MAXAMP,
    RMS,
    STDDEV,
    MAGNITUDE,
    POWERARRAY,
    SILENCEPROB,

    // Spectral Arrays
    MFCC,
    CHROMA,
    MELOGRAM,

    // Spectral
    ZCR,
    HFR,
    CENTROID,
    SPREADHZ,
    SPREADVARIANCE,
    CREST,
    FLATNESS,
    ENTROPY,
    ROLLOFF,
    CENTROIDVEL,
    FLUX,
    SKEWNESS,
    SLOPE,
    KURTOSIS,
    IRREGULARITY,
    HARMONICITY,

    // Pitch
    YIN,
    YINCONFIDENCE,

    // Percussive
    EXTENDEDPROB,

    // AI
    ONNX,
};

// ─────────────────────────────────────
enum Mode { SCOREFOLLOWER, DESCRIPTORS };

// ─────────────────────────────────────
enum EventType {
    // Tradicional
    REST,
    NOTE,
    CHORD,
    TRILL,
    MULTI,
    PTECH,
    UTECH,

    // Cort Lippe event (TODO:)
    EVENT,
};

enum HMMType { SEMIMARKOV, MARKOV };

// ─────────────────────────────────────
class ScoreAction {
  public:
    bool isLua;
    std::string Lua;
    std::string Receiver;
    std::vector<std::variant<float, int, std::string>> Args;
    bool AbsoluteTime;
    double Time;
};

using EventActions = std::vector<ScoreAction>;

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
    int MarkovIndex = -1;
    std::vector<AudioState> AudioStates;

    // States Actions
    HMMType HSMMType;
    EventType Type;
    EventActions Actions;

    // Inference
    double InitProb;
    std::vector<double> Forward;
    std::vector<double> ExitProb;
    std::vector<double> BestObs;

    // Time
    int UpperBound = 0;
    double BPMExpected = 0;
    double BPMObserved = 0;
    double OnsetExpected = 0.0;
    double OnsetObserved = 0;
    double PhaseExpected;
    double PhaseObserved = 0;
    double IOIPhiN = 0;
    double IOIHatPhiN = 0;
    double Duration = 0.0;

    // Configuration for each event
    double PhaseCoupling = 0;
    double SyncStrength = 0;

    // Error Handling
    int Line;
};

using States = std::vector<MarkovState>;

// ─────────────────────────────────────
class Description {
  public:
    double Onset;

    // Probability
    double SilenceProb;
    double ExtendedTechProb;
    int WindowLastOnset = 0;

    // Amplitude
    double dB;
    double RMS;
    double MaxAmp;
    double Loudness;

    // Spectral
    double Harmonicity;
    double SpectralFlatness;
    double SpectralEntropy = 0.0;
    double SpectralRolloff = 0.0;
    double SpectralFlux;
    double SpectralIrregularity = 0.0;
    double SpectralCrest = 0.0;
    double SpectralCentroid = 0.0;
    double CentroidVelocity = 0.0;
    double SpectralSpreadHz = 0.0;
    double SpectralSpreadVariance = 0.0;
    double SpectralSkewness = 0.0;
    double SpectralSlope = 0.0;
    double SpectralKurtosis = 0.0;
    double HighFreqRatio = 0.0;
    double ZeroCrossingRate;
    double StdDev;
    double Pitch = 0.0;
    double PitchConfidence = 0.0;

    std::vector<double> Power;
    std::vector<double> Magnitude;
    std::vector<double> SpectralMagnitudeNorm;
    std::vector<double> SpectralMagnitudeFrameNorm;

    std::vector<double> LogMelSpectrum;
    std::vector<double> ReverbSpectralPower;

    std::vector<double> MFCC;
    std::vector<double> Chroma;

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
