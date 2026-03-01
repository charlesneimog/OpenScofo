#pragma once
#include <unordered_map>
#include <vector>

#include "states.hpp"

namespace OpenScofo {

#ifndef TWO_PI
#define TWO_PI (2 * M_PI)
#endif

using PitchTemplateArray = std::vector<double>;

// ╭─────────────────────────────────────╮
// │     Markov Description Process      │
// ╰─────────────────────────────────────╯
class MDP {
  public:
    MDP(double Sr, double WindowSize, double HopSize);

    // Init Functions
    void SetScoreStates(States States);

    void UpdateAudioTemplate();
    void UpdatePhaseValues();

    // Set
    void SetPitchTemplateSigma(double f);
    void SetHarmonics(int i);
    void SetBPM(double Bpm);
    void SetMinEntropy(double EntropyValue);
    void SetAmplitudeDecay(double decay);
    void SetTunning(double Tunning);
    void SetCurrentEvent(int Event);
    void SetdBTreshold(double dB);

    // Get
    double GetLiveBPM();
    void ResetLiveBpm();

    //
    double GetBlockDuration();

    // Get Functions
    int GetTunning();
    ActionVec GetEventActions(int Index);
    std::vector<MarkovState> GetStates();
    MarkovState GetState(int Index);
    double GetKappa();
    void AddState(MarkovState state);

    int GetStatesSize();
    int GetEvent(Description &Desc);
    double GetPitchSimilarity(double Freq);
    States GetStatesForProcessing();

    // Python For Research
    PitchTemplateArray GetPitchTemplate(double Freq);

    // Set Variables
    void ClearStates();

  private:
    void GetDecodeWindow();

    // Time
    double UpdatePsiN(int StateIndex);
    double InverseA2(double r);
    double ModPhases(double value);
    double CouplingFunction(double Phi, double PhiMu, double Kappa);
    double GetOccupancyDistribution(MarkovState &State, int u);
    double GetSurvivorDistribution(MarkovState &State, int u);
    void InitTimeDecoding();
    void BuildDistributionCache(double ExpectedFrames);

    // Markov and Probabilities
    double GetTransProbability(int i, int j);
    void GetInitialDistribution();
    int GetMaxUForJ(MarkovState &StateJ);

    // Markov
    double GetBestEvent();
    int GetMaxJIndex(int StateIndex);

    void Markov(MarkovState &StateJ, int j, int T, int bufferIndex);
    void SemiMarkov(MarkovState &StateJ, int j, int T, int bufferIndex);

    int Inference(int CurrentState);

    // Pitch Template
    void BuildPitchTemplate(double Freq);

    // Get Audio Obs
    void GetAudioObservations(int T);

  private:
    // Test things
    int m_WinEnd;
    int m_WinStart;
    int m_EventWindowSize;
    std::vector<double> m_Normalization;

    // Config
    double m_MinEntropy = 0;

    // Audio
    double m_Sr;
    double m_FFTSize;
    double m_HopSize;
    double m_dBTreshold = -55;
    int m_BufferSize = 1000;
    bool m_IsSilence = false;

    // Pitch Template
    double m_Harmonics = 8;
    double m_PitchTemplateSigma = 0.387;
    double m_PitchTemplateAmplitudeDecay = 0.87270;

    // Events
    double m_Tunning = 440;
    int m_CurrentStateIndex = 0;

    // Time
    double m_SyncStrength = 0.5;
    double m_PhaseCoupling = 0.5;
    double m_SyncStr = 0;
    double m_TimeInPrevEvent = 0;
    std::unordered_map<double, double> m_KappaCache;
    // Cache for the distributions
    std::unordered_map<int, std::vector<double>> m_OccupancyCache;
    std::unordered_map<int, std::vector<double>> m_SurvivorCache;

    double m_LastTn = 0;
    double m_TauWithSound = 0;
    double m_BlockDur = 0;
    double m_CurrentStateOnset = 0;
    int m_MaxScoreState = 0;

    int m_Tau = 0;
    double m_LastPsiN = 0;
    double m_PsiN = 0;
    double m_PsiN1 = 0;
    double m_BPM = 0;
    double m_Kappa = 1;
    double m_MaxAheadSeconds;
    double m_BeatsAhead = 1;
    double m_NormAlpha = 1;
    double m_SecondsAhead = 2;

    // Time
    std::unordered_map<double, PitchTemplateArray> m_PitchTemplates;

    // Pitch
    std::vector<MarkovState> m_States;
    double m_PitchScalingFactor = 0.5;

    std::unordered_map<double, PitchTemplateArray> m_PitchCQTTemplates;

    // Audio Observations
    Description m_Desc;

    // Markov
    bool m_EventDetected = false;

    // Errors
    bool m_HasErrors = false;
};
} // namespace OpenScofo
