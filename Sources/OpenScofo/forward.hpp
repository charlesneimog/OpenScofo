/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#pragma once
#include <array>
#include <cstddef>
#include <deque>
#include <unordered_map>
#include <vector>

#include <states.hpp>

namespace OpenScofo {

// ─────────────────────────────────────
using PitchTemplateArray = std::vector<double>;

// ─────────────────────────────────────
struct PitchTemplateEntry {
    size_t bin;
    double p;
    double p_log_p;
};

// ╭─────────────────────────────────────╮
// │     Markov Description Process      │
// ╰─────────────────────────────────────╯
class OnlineForward {
  public:
    OnlineForward();
    void UpdateConfiguration(Configuration &Config);
    void SetScoreStates(States States);
    void UpdateAudioTemplate();
    void ResetDecoding();

    // Set
    void SetPitchTemplateSigma(double f);
    void SetHarmonics(int i);
    void SetAmplitudeDecay(double decay);
    void SetTunning(double Tunning);
    void SetDescription(const Description &Desc);
    void SetCurrentEvent(int Event);

    // Get
    double GetCurrentBPM();
    void ResetLiveBpm();

    //
    double GetBlockDuration();

    // Get Functions
    int GetCurrentBufferIndex();
    int GetCurrentStateIndex();
    int GetTunning();
    EventActions GetCurrentEventActions();
    EventActions GetAudioStateChangeActions();
    std::vector<MarkovState> &GetStates();
    MarkovState GetState(int Index);
    void AddState(MarkovState state);

    int GetStatesSize();
    int GetEvent(Description &Desc);
    double GetPitchProbability(double Freq);
    States GetStatesForProcessing();

    // Python For Research
    PitchTemplateArray GetPitchTemplate(double Freq);

    // Set Variables
    void ClearStates();

  private:
    void GetDecodeWindow();

    // Time
    double UpdatePsiN(int StateIndex);
    double A2(double kappa);
    double InverseA2(double r);
    void InitializeA2Table();
    static double CalculateA2(double kappa);
    double ModPhases(double value);
    double CouplingFunction(double Phi, double PhiMu, double Kappa);
    double GetOccupancyDistribution(MarkovState &State, int u);
    double GetSurvivorDistribution(MarkovState &State, int u);
    void InitTimeDecoding();
    void BuildDistributionCache(double ExpectedFrames);
    void ResetCaches();

    // Markov and Probabilities
    double GetTransProbability(int i, int j);
    void GetInitialDistribution();
    int GetMaxUForJ(MarkovState &StateJ);

    // Markov
    double GetBestEvent();
    int GetMaxJIndex(int StateIndex);

    void Markov(MarkovState &StateJ, int j, int bufferIndex);
    void SemiMarkov(MarkovState &StateJ, int j, int bufferIndex);

    int GetAlphaT();

    // Pitch Template
    void BuildPitchTemplate(double Freq);

    // Get Audio Obs
    void GetAudioObservations();
    void NotifyAudioStateChange(int StateIndex);

  private:
    // Test things
    int m_WinEnd;
    int m_WinStart;
    int m_EventWindowSize;
    std::vector<double> m_Normalization;

    int m_BestAudioStateStateIndex;
    int m_BestAudioStateIndex;

    // Config
    double m_MinEntropy = 0;
    AudioDescType m_CurrentAudioState;
    std::string m_AudioStateChangeReceiver;

    // Audio-state change notifications
    int m_LastNotifiedStateIndex = -1;
    int m_LastNotifiedAudioStateIndex = -1;
    std::deque<ScoreAction> m_PendingAudioStateActions;

    // Audio
    double m_Sr;
    double m_FFTSize;
    double m_HopSize;
    double m_dBTreshold = -55;
    int m_BufferSize = 2000;
    bool m_IsSilence = false;

    // Pitch Template
    double m_Harmonics = 8;
    double m_PitchTemplateSigma = 0.5;
    double m_PitchTemplateAmplitudeDecay = 0.5;

    // Events
    double m_Tunning = 440;
    int m_CurrentStateIndex = 0;

    // Time
    double m_TimeInPrevEvent = 0;
    std::unordered_map<int, double> m_KappaCache;
    static constexpr int A2TablePrecision = 100;
    static constexpr int A2TableMaxKappa = 10;
    static constexpr int A2TableSize = A2TablePrecision * A2TableMaxKappa + 1;
    std::array<double, A2TableSize> m_A2Table = {};
    bool m_A2TableInitialized = false;
    // Cache for the distributions
    std::unordered_map<int, std::vector<double>> m_OccupancyPMFCache;
    std::unordered_map<int, std::vector<double>> m_SurvivorCache;

    double m_LastTn = 0;
    double m_BlockDur = 0;
    double m_CurrentStateOnset = 0;
    int m_MaxScoreState = 0;

    int m_Tau = 0;
    double m_LastPsiN = 0;
    double m_PsiN = 0;
    double m_PsiN1 = 0;
    double m_BPM = 0;

    double m_Kappa = 5;
    double m_SyncStrength = 0.5;
    double m_PhaseCoupling = 0.5;
    double m_SyncStr = 0;

    double m_MaxAheadSeconds;
    double m_BeatsAhead = 1;
    double m_NormAlpha = 1;
    double m_SecondsAhead = 2;

    // Time
    std::unordered_map<double, PitchTemplateArray> m_PitchTemplates;
    std::unordered_map<double, std::vector<PitchTemplateEntry>> m_PitchTemplatesPrecomputed;
    std::unordered_map<double, double> m_PitchProbabilityCache;
    int m_PitchProbabilityCacheTau = -1;
    int m_ReverbEnergyCacheTau = -1;
    bool m_ReverbSpectralPowerHasEnergy = false;

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
