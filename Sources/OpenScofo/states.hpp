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

    // AI
    ONNX,

    // Percussive
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

using EventActions = std::vector<Action>;

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
    std::vector<double> Magnitude;
    std::vector<double> SpectralMagnitudeNorm;
    std::vector<double> SpectralMagnitudeFrameNorm;

    std::vector<double> Melogram;
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
