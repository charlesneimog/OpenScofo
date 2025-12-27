#pragma once

#include <string>
#include <vector>
#include <variant>

namespace OScofo {

// ─────────────────────────────────────
enum EventType {
    REST,  // Markov state
    NOTE,  // MacroState, Semimarkov with Markov inside
    CHORD, // MacroState,
    TRILL, // MacroState, SemiMarkov with Markov inside
    MULTI, // MacroState,
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
    EventType Type;
    HMMType Markov;
    double Freq;
    std::vector<double> Obs;
    std::vector<double> Forward;
    unsigned Index;
};

// ─────────────────────────────────────
class MacroState {
  public:
    int Index;
    int ScorePos;
    int SubStateIndex = 0;
    EventType Type;
    HMMType Markov;
    int MarkovIndex = -1;

    // Configs
    double Entropy;

    // States Actions
    ActionVec Actions;
    std::vector<AudioState> SubStates;

    // Forward Algorithm
    double InitProb;
    std::vector<double> Obs;
    std::vector<double> Forward;

    // Time
    double OnsetTime;
    // double Duration;

    double Sigma;
    double PdfStartTime;
    std::vector<double> Pdf;
    std::vector<double> Context;

    // Audio Obs
    std::vector<double> Freqs;

    // Time
    int UpperBound;
    double BPMExpected;
    double BPMObserved = 0;
    double OnsetExpected;
    double OnsetObserved = 0;
    double PhaseExpected;
    double PhaseObserved = 0;
    double IOIPhiN;
    double IOIHatPhiN;
    double Duration;

    double PhaseCoupling;
    double SyncStrength;

    // Error Handling
    int Line;

    std::string __repr__() const {
        std::string oss;
        oss = "<<State(Index=";
        oss += std::to_string(Index);
        oss += ", ScorePosition=";
        oss += std::to_string(ScorePos);
        oss += ", BPMExpected=";
        oss += std::to_string(BPMExpected);
        oss += ")>>";
        return oss;
    }
};

using States = std::vector<MacroState>;

// ─────────────────────────────────────
class Description {
  public:
    bool Silence;
    bool Onset;
    double SilenceProb;

    double dB;
    double RMS;
    double MaxAmp;
    double Loudness;

    double SpectralFlatness;
    double SpectralFlux;
    double StdDev;

    std::vector<double> Power;
    std::vector<double> SpectralPower;
    std::vector<double> NormSpectralPower;
    std::vector<double> MFCC;
};

} // namespace OScofo
