#include "mdp.hpp"
#include "log.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <boost/math/special_functions/bessel.hpp>
#include <numeric>

namespace OScofo {

/*
    // ──────────────────────────────── REFERENCES ───────────────────────────────────────

    * GONG, R.; CUVILLIER, P.; OBIN, N.; CONT, A. Real-Time Audio-to-Score Alignment of Sin-
        ging Voice Based on Melody and Lyric Information. In: Interspeech, 2015, Dresde, Germany.
        Anais… [S.l.: s.n.], 2015.

    * CONT, A. A Coupled Duration-Focused Architecture for Real-Time Music-to-Score Align-
        ment. IEEE Transactions on Pattern Analysis and Machine Intelligence, [S.l.], v.32, n.6,
        p.974–987, 2010.

    * CONT, A. Improvement of Observation Modeling for Score Following. 2004.

    * GUÉDON, Y. Hidden Hybrid Markov/Semi-Markov Chains. Computational Statistics & Data
        Analysis, [S.l.], v.49, n.3, p.663–688, 2005.

    * LARGE, E. W.; JONES, M. R. The Dynamics of Attending: How People Track Time-Varying
        Events. Psychological Review, [S.l.], v.106, n.1, p.119–159, 1999.

    * LARGE, E. W.; PALMER, C. Perceiving Temporal Regularity in Music. Cognitive Science,
        [S.l.], v.26, n.1, p.1–37, 2002.
*/

// ╭─────────────────────────────────────╮
// │Constructor and Destructor Functions │
// ╰─────────────────────────────────────╯
MDP::MDP(double Sr, double FFTSize, double HopSize) {
    m_HopSize = HopSize;
    m_FFTSize = FFTSize;
    m_Sr = Sr;

    m_SyncStrength = 0.5;
    m_PhaseCoupling = 0.5;
    m_BlockDur = (1 / m_Sr) * HopSize;
    m_TimeInPrevEvent = 0;

    SetTunning(440);
}

// ╭─────────────────────────────────────╮
// │               Errors                │
// ╰─────────────────────────────────────╯
bool MDP::HasErrors() {
    return m_HasErrors;
}

// ─────────────────────────────────────
std::vector<std::string> MDP::GetErrorMessage() {
    return m_Errors;
}

// ─────────────────────────────────────
void MDP::SetError(const std::string &message) {
    printf("Error: %s.\n", message.c_str());
    m_HasErrors = true;
    m_Errors.push_back(message);
}

// ─────────────────────────────────────
void MDP::ClearError() {
    m_HasErrors = false;
    m_Errors.clear();
}

// ─────────────────────────────────────
ActionVec MDP::GetEventActions(int Index) {
    if (Index < 0 || Index >= (int)m_States.size()) {
        return ActionVec();
    }

    MacroState State = m_States[(size_t)Index];
    return State.Actions;
}

// ─────────────────────────────────────
void MDP::SetScoreStates(States ScoreStates) {
    if (ScoreStates.size() == 0) {
        SetError("ScoreStates is empty, add some events to the score");
        return;
    }

    m_States.clear();
    m_States = ScoreStates;

    for (MacroState &State : m_States) {
        // State.Obs.resize(m_BufferSize + 1, 0);
        State.Forward.resize(m_BufferSize + 1, 0);
        for (AudioState &MicroState : State.AudioStates) {
            MicroState.Obs.resize(m_BufferSize + 1, 0);
            MicroState.Forward.resize(m_BufferSize + 1, 0);
        }
    }

    m_CurrentStateIndex = -1;
    m_Kappa = 1;
    m_BPM = m_States[0].BPMExpected;
    m_PsiN = 60.0f / m_States[0].BPMExpected;
    m_PsiN1 = 60.0f / m_States[0].BPMExpected;
    m_LastPsiN = 60.0f / m_States[0].BPMExpected;
    m_BeatsAhead = m_States[0].BPMExpected / 60 * m_SecondsAhead;
    m_CurrentStateIndex = -1;
    m_SyncStr = 0;

    UpdateAudioTemplate();
    UpdatePhaseValues();
}

// ─────────────────────────────────────
// void MDP::GetPitchTemplate(double Freq) {
//     // Following Gong (2015), eq 5 and 6
//     const double sigmaLog = m_PitchTemplateSigma / 12.0;
//     const double beta = m_PitchTemplateAmplitudeDecay;
//
//     double rootBinFreq = round(Freq / (m_Sr / m_FFTSize));
//     if (m_PitchTemplates.find(rootBinFreq) != m_PitchTemplates.end()) {
//         return;
//     }
//
//     m_PitchTemplates[rootBinFreq].resize(m_FFTSize / 2, 0.0);
//     for (int k = 1; k <= m_Harmonics; k++) {
//         double harmonicFreqHz = Freq * k;
//         double sigmaHz = harmonicFreqHz * (std::pow(2.0, sigmaLog) - 1.0);
//         double envelope = std::exp(-beta * (k - 1));
//         for (size_t i = 0; i < m_FFTSize / 2; i++) {
//             double binFreq = i * (m_Sr / m_FFTSize);
//             double exponent = -0.5 * std::pow((binFreq - harmonicFreqHz) / sigmaHz, 2);
//             double gaussian = (1.0 / (sigmaHz * std::sqrt(2 * M_PI))) * std::exp(exponent);
//             m_PitchTemplates[rootBinFreq][i] += envelope * gaussian;
//         }
//     }
//
//     // Normalize template to sum to 1 (probability distribution)
//     double sum = std::accumulate(m_PitchTemplates[rootBinFreq].begin(), m_PitchTemplates[rootBinFreq].end(), 0.0);
//     if (sum > 0) {
//         for (auto &val : m_PitchTemplates[rootBinFreq]) {
//             val = (val + 1e-24) / (sum + 1e-24); // Avoid zero probabilities
//         }
//     }
// }

// ─────────────────────────────────────
void MDP::BuildPitchTemplate(double Freq) {
    // 1. Setup Constants
    // ------------------
    const double m_MinHarmonicDecay = 0.2;
    const double m_MaxHarmonicDecay = 1.8;
    const double binWidth = m_Sr / m_FFTSize;

    if (Freq > m_Sr / 2.0 || Freq <= 0)
        return;

    // 2. Determine Beta (Harmonic Decay)
    // ------------------
    double m_MinF0 = 32.70;
    double m_MaxF0 = 4186.0;
    double f0Norm = std::log2(Freq / m_MinF0) / std::log2(m_MaxF0 / m_MinF0);
    f0Norm = std::clamp(f0Norm, 0.0, 1.0);

    double pitchDecayCurve = m_PitchTemplateAmplitudeDecay;
    double shaped = std::pow(1.0 - f0Norm, pitchDecayCurve);
    double beta = m_MaxHarmonicDecay - shaped * (m_MaxHarmonicDecay - m_MinHarmonicDecay);

    // 3. Template Initialization
    // ------------------
    double rootBinFreq = round(Freq / binWidth);

    if (m_PitchTemplates.find(rootBinFreq) != m_PitchTemplates.end()) {
        return;
    }

    std::vector<double> templateBins(m_FFTSize / 2, 1e-24);

    // 4. Build Harmonics
    // ------------------
    const double sigmaLog = m_PitchTemplateSigma / 12.0;
    const double sigmaConst = std::pow(2.0, sigmaLog) - 1.0;

    // Inharmonicity factor (B): typical for piano/strings is 0.0001
    // If your library is purely electronic, set to 0.0
    const double B = 0.0001;

    for (int k = 1; k <= m_Harmonics; ++k) {
        // Apply inharmonicity stretch: f_k = k * f0 * sqrt(1 + B * k^2)
        double stretch = std::sqrt(1.0 + B * k * k);
        double harmonicFreqHz = Freq * k * stretch;

        if (harmonicFreqHz >= m_Sr / 2.0)
            break;

        // A. Calculate Sigma
        double sigmaHz = harmonicFreqHz * sigmaConst;

        // Sigma Floor to prevent "picket fence" logic errors
        if (sigmaHz < binWidth * 0.75) {
            sigmaHz = binWidth * 0.75;
        }

        // B. Calculate Amplitude Envelope with Harmonic Boost
        double envelope = std::exp(-beta * (k - 1));

        // INTERVENTION: Boost higher harmonics.
        // Real-world instruments often have a "formant" or resonance
        // that makes harmonics 2-6 stronger than the fundamental.
        if (k > 1) {
            envelope *= 1.25;
        }

        if (envelope < 1e-5)
            break;

        // C. Optimized Loop Bounds
        double rangeHz = 4.0 * sigmaHz;
        int minBin = static_cast<int>(std::floor((harmonicFreqHz - rangeHz) / binWidth));
        int maxBin = static_cast<int>(std::ceil((harmonicFreqHz + rangeHz) / binWidth));

        minBin = std::max(0, minBin);
        maxBin = std::min(static_cast<int>(m_FFTSize / 2) - 1, maxBin);

        double twoSigmaSq = 2.0 * sigmaHz * sigmaHz;
        double normalizationFactor = 1.0 / (sigmaHz * std::sqrt(2.0 * M_PI));

        for (int i = minBin; i <= maxBin; ++i) {
            double binFreq = i * binWidth;
            double diff = binFreq - harmonicFreqHz;
            double exponent = -(diff * diff) / twoSigmaSq;
            double gaussian = normalizationFactor * std::exp(exponent);

            templateBins[i] += envelope * gaussian;
        }
    }

    // 5. Normalization
    // ------------------
    double sum = std::accumulate(templateBins.begin(), templateBins.end(), 0.0);

    if (sum > 1e-12) {
        double invSum = 1.0 / sum;
        for (auto &val : templateBins) {
            val *= invSum;
        }
    }

    m_PitchTemplates[rootBinFreq] = std::move(templateBins);
}

// ─────────────────────────────────────
void MDP::UpdateAudioTemplate() {
    int StateSize = (int)m_States.size();
    m_PitchTemplates.clear();

    for (int h = 0; h < StateSize; h++) {
        if (m_States[h].Type == NOTE || m_States[h].Type == TRILL) {
            for (AudioState &SubState : m_States[h].AudioStates) {
                if (SubState.Type == PITCH) {
                    BuildPitchTemplate(SubState.Freq);
                }
            }
        }
    }
}

// ─────────────────────────────────────
PitchTemplateArray MDP::GetPitchTemplate(double Freq) {
    BuildPitchTemplate(Freq);
    double rootBinFreq = round(Freq / (m_Sr / m_FFTSize));
    return m_PitchTemplates[rootBinFreq];
}

// ─────────────────────────────────────
void MDP::UpdatePhaseValues() {
}

// ╭─────────────────────────────────────╮
// │          Set|Get Functions          │
// ╰─────────────────────────────────────╯
void MDP::ClearStates() {
    m_States.clear();
}
// ─────────────────────────────────────
double MDP::GetLiveBPM() {
    return m_BPM;
}

// ─────────────────────────────────────
double MDP::GetKappa() {
    return m_Kappa;
}

// ─────────────────────────────────────
void MDP::SetBPM(double BPM) {
    m_BPM = BPM;
}

// ─────────────────────────────────────
void MDP::SetdBTreshold(double dB) {
    m_dBTreshold = dB;
}

// ─────────────────────────────────────
void MDP::SetTunning(double Tunning) {
    m_Tunning = Tunning;
}

// ─────────────────────────────────────
void MDP::SetHarmonics(int Harmonics) {
    m_Harmonics = Harmonics;
}

// ─────────────────────────────────────
void MDP::SetMinEntropy(double EntropyValue) {
    m_MinEntropy = EntropyValue;
}

// ─────────────────────────────────────
void MDP::SetAmplitudeDecay(double decay) {
    m_PitchTemplateAmplitudeDecay = decay;
}

// ─────────────────────────────────────
int MDP::GetTunning() {
    return m_Tunning;
}

// ─────────────────────────────────────
void MDP::SetCurrentEvent(int Event) {
    m_CurrentStateIndex = Event;
    m_Tau = 0;
}

// ─────────────────────────────────────
int MDP::GetStatesSize() {
    return m_States.size();
}

// ─────────────────────────────────────
void MDP::AddState(MacroState State) {
    m_States.push_back(State);
}

// ─────────────────────────────────────
MacroState MDP::GetState(int Index) {
    return m_States[Index];
}

// ─────────────────────────────────────
void MDP::SetPitchTemplateSigma(double f) {
    m_PitchTemplateSigma = f;
}

// ╭─────────────────────────────────────╮
// │            Time Decoding            │
// ╰─────────────────────────────────────╯
double MDP::InverseA2(double SyncStrength) {
    // SyncStrength must be between 0 and 1
    if (SyncStrength < 0) {
        return 0;
    }

    // Following Large and Jones (1999, p. 157).
    if (SyncStrength > 0.95) {
        return 10.0f;
    }

    double Low = 0.0;
    double Tol = 1e-16;
    double High = std::max(SyncStrength, 10.0);
    double Mid;

    // In my tests I never reached more than 100 iterations.
    int i;
    for (i = 0; i < 1000; i++) {
        Mid = (Low + High) / 2.0;
        // TODO: Some way to avoid to use boost?
        double I1 = boost::math::cyl_bessel_i(1, Mid);
        double I0 = boost::math::cyl_bessel_i(0, Mid);
        double A2Mid = I1 / I0;
        if (std::fabs(A2Mid - SyncStrength) < Tol) {
            return Mid;
        } else if (A2Mid < SyncStrength) {
            Low = Mid;
        } else {
            High = Mid;
        }
    }
    LOGE() << "InverseA2 not converged after " << i << " iterations.";
    return Mid;
}

// ─────────────────────────────────────
double MDP::CouplingFunction(double Phi, double PhiMu, double Kappa) {
    // Equation 2b from Large and Palmer (2002)
    double ExpKappa = exp(Kappa);
    double PhiNDiff = Phi - PhiMu;
    double CosTerm = cos(TWO_PI * PhiNDiff);
    double SinTerm = sin(TWO_PI * PhiNDiff);
    double PhiN = (1 / (TWO_PI * ExpKappa)) * exp(Kappa * CosTerm) * SinTerm;
    return PhiN;
}

// ─────────────────────────────────────
double MDP::ModPhases(double Phase) {
    // Following Cont (2010) conventions
    Phase = std::fmod(Phase + M_PI, TWO_PI);
    if (Phase < 0) {
        Phase += TWO_PI;
    }
    return Phase - M_PI;
}

// ─────────────────────────────────────
States MDP::GetStatesForProcessing() {
    if (m_CurrentStateIndex < 0)
        m_CurrentStateIndex = 0;

    double EventOnset = m_States[m_CurrentStateIndex].Duration - (m_TimeInPrevEvent + m_BlockDur);
    size_t begin = m_CurrentStateIndex + 1;
    size_t end = begin;

    for (; end < m_States.size(); ++end) {
        if (EventOnset > m_SecondsAhead)
            break;

        EventOnset += m_States[end].Duration;
    }

    return States(m_States.begin() + begin, m_States.begin() + end);
}

// ─────────────────────────────────────
int MDP::GetMaxJIndex(int StateIndex) {
    if (StateIndex == -1) {
        return 1;
    }
    double EventOnset = m_States[StateIndex].Duration - (m_TimeInPrevEvent + m_BlockDur);
    int MaxJ = StateIndex + 1;
    for (size_t i = StateIndex + 1; i < m_States.size(); i++) {
        if (EventOnset > m_SecondsAhead) {
            MaxJ = i;
            break;
        } else {
            EventOnset += m_States[i].Duration;
        }
    }

    return MaxJ;
}

// ─────────────────────────────────────
double MDP::UpdatePsiN(int StateIndex) {
    if (StateIndex == m_CurrentStateIndex) {
        m_TimeInPrevEvent += m_BlockDur;
        m_Tau += 1;
        return m_PsiN;
    } else {
        if (StateIndex == 0) {
            double PsiK = 60 / m_States[0].BPMExpected;
            m_LastPsiN = PsiK;
            m_PsiN = PsiK;
            m_PsiN1 = PsiK;
            m_States[0].OnsetObserved = 0;
            m_BPM = m_States[0].BPMExpected;
            m_CurrentStateOnset = 0;
            m_LastTn = 0;
            m_TimeInPrevEvent = 0;
            m_Tau = 0;
            return m_PsiN;
        } else {
            m_TimeInPrevEvent += m_BlockDur;
            m_LastTn = m_CurrentStateOnset;
            m_CurrentStateOnset += m_TimeInPrevEvent;
        }
    }

    // Cont (2010), Large and Palmer (1999) and Large and Jones (2002)
    MacroState &LastState = m_States[StateIndex - 1];
    MacroState &CurrentState = m_States[StateIndex];
    MacroState &NextState = m_States[StateIndex + 1];

    double IOISeconds = m_CurrentStateOnset - m_LastTn;
    double LastPhiN = LastState.IOIPhiN;
    double LastHatPhiN = LastState.IOIHatPhiN;
    double HatPhiN = CurrentState.IOIHatPhiN;
    double PhiNExpected = LastPhiN + ((m_CurrentStateOnset - m_LastTn) / m_PsiN);
    CurrentState.IOIHatPhiN = PhiNExpected;
    CurrentState.OnsetObserved = m_CurrentStateOnset;

    // Update Variance (Cont, 2010) - Coupling Strength (Large 1999)
    double PhaseDiff = (IOISeconds / m_PsiN) - HatPhiN;
    double SyncStrength = m_SyncStr - m_SyncStrength * (m_SyncStr - cos(TWO_PI * PhaseDiff));
    double Kappa = InverseA2(SyncStrength);
    m_SyncStr = SyncStrength;
    m_Kappa = Kappa;

    // Update and Correct PhiN
    double FValueUpdate = CouplingFunction(LastPhiN, LastHatPhiN, Kappa);
    double PhiN = LastPhiN + (IOISeconds / m_LastPsiN) + (m_PhaseCoupling * FValueUpdate);
    PhiN = ModPhases(PhiN);
    CurrentState.PhaseObserved = PhiN;

    // Prediction for next PsiN+1
    double FValuePrediction = CouplingFunction(PhiN, HatPhiN, Kappa);
    double PsiN1 = m_PsiN * (1 + m_SyncStrength * FValuePrediction);

    // Prediction for Next HatPhiN
    double Tn1 = m_CurrentStateOnset + CurrentState.Duration * PsiN1;
    double PhiN1 = ModPhases((Tn1 - m_CurrentStateOnset) / PsiN1);
    NextState.IOIHatPhiN = PhiN1;

    // Update all next expected onsets
    NextState.OnsetExpected = Tn1;
    double LastOnsetExpected = Tn1;

    // the m_CurrentEvent + 1 already updated, now
    // we update the future events to get the Sojourn Time
    for (int i = m_CurrentStateIndex + 2; i < m_CurrentStateIndex + 20; i++) {
        if ((size_t)i >= m_States.size()) {
            break;
        }
        MacroState &FutureState = m_States[i];
        MacroState &PreviousFutureState = m_States[(i - 1)];
        double Duration = PreviousFutureState.Duration;
        double FutureOnset = LastOnsetExpected + Duration * PsiN1;

        FutureState.OnsetExpected = FutureOnset;
        LastOnsetExpected = FutureOnset;
    }

    // Update Values for next calls
    m_BPM = 60.0f / m_PsiN;
    m_LastPsiN = m_PsiN;

    if (StateIndex != m_CurrentStateIndex) {
        m_TimeInPrevEvent = 0;
        m_Tau = 0;
    }
    return PsiN1;
}

// ╭─────────────────────────────────────╮
// │     Markov Description Process      │
// ╰─────────────────────────────────────╯
void MDP::GetAudioObservations(int FirstStateIndex, int LastStateIndex, int T) {
    std::unordered_map<double, double> PitchObs;

    for (int j = FirstStateIndex; j <= LastStateIndex; j++) {
        if (j < 0) {
            continue;
        }

        MacroState &StateJ = m_States[j];
        int BufferIndex = (T % m_BufferSize);
        if (StateJ.Type == NOTE) {
            for (AudioState &AudioState : StateJ.AudioStates) {
                switch (AudioState.Type) {
                case PITCH: {
                    if (PitchObs.find(AudioState.Freq) != PitchObs.end()) {
                        AudioState.Obs[BufferIndex] = PitchObs[AudioState.Freq];
                        continue;
                    }
                    double KL = GetPitchSimilarity(AudioState.Freq) * (1 - m_Desc.SilenceProb);
                    PitchObs[AudioState.Freq] = KL;
                    AudioState.Obs[BufferIndex] = KL;
                    break;
                }
                case SILENCE: {
                    AudioState.Obs[BufferIndex] = m_Desc.SilenceProb;
                    break;
                }
                default:
                    SetError("Not implemented yet!");
                }
            }
        } else if (StateJ.Type == REST) {
            for (AudioState &AudioState : StateJ.AudioStates) {
                AudioState.Obs[BufferIndex] = m_Desc.SilenceProb;
            }
        } else if (StateJ.Type == TRILL) {
            for (AudioState &AudioState : StateJ.AudioStates) {
                switch (AudioState.Type) {
                case PITCH: {
                    if (PitchObs.find(AudioState.Freq) != PitchObs.end()) {
                        AudioState.Obs[BufferIndex] = PitchObs[AudioState.Freq];
                        continue;
                    }
                    double KL = GetPitchSimilarity(AudioState.Freq) * (1 - m_Desc.SilenceProb);
                    AudioState.Obs[BufferIndex] = KL;
                    break;
                }
                case SILENCE: {
                    AudioState.Obs[BufferIndex] = m_Desc.SilenceProb;
                    break;
                }
                default:
                    SetError("Not implemented yet!");
                }
            }
        }
    }
}

// ─────────────────────────────────────
double MDP::GetPitchSimilarity(double Freq) {
    double KLDiv = 0.0;
    double RootBinFreq = round(Freq / (m_Sr / m_FFTSize));
    PitchTemplateArray PitchTemplate;

    if (m_PitchTemplates.find(RootBinFreq) != m_PitchTemplates.end()) {
        PitchTemplate = m_PitchTemplates[RootBinFreq];
    } else {
        BuildPitchTemplate(Freq);
        PitchTemplate = m_PitchTemplates[RootBinFreq];
    }

    for (size_t i = 0; i < m_FFTSize / 2; i++) {
        double P = PitchTemplate[i];
        double Q = m_Desc.NormSpectralPower[i];
        if (P > 0 && Q > 0) {
            KLDiv += P * log(P / Q);
        } else if (P == 0 && Q >= 0) {
            KLDiv += Q;
        }
    }

    // KLDiv /= (m_Desc.StdDev + 1e-9);
    double noise_robustness = 1.0 / (1.0 + m_Desc.StdDev);
    KLDiv *= noise_robustness;

    KLDiv = exp(-m_PitchScalingFactor * KLDiv);
    return KLDiv;
}

// ─────────────────────────────────────
std::vector<double> MDP::GetInitialDistribution() {
    int Size = m_MaxScoreState - m_CurrentStateIndex;
    std::vector<double> InitialProb(Size);

    double Dur = 0;
    double Sum = 0;

    for (int i = 0; i < Size; i++) {
        if (m_CurrentStateIndex + i < 0) {
            continue;
        }
        double DurProb = exp(-1 * (Dur / m_BeatsAhead));
        InitialProb[i] = DurProb;
        Dur += m_States[m_CurrentStateIndex + i].Duration; // Accumulate duration
        Sum += DurProb;
    }

    // Normalize
    for (int i = 0; i < Size; i++) {
        InitialProb[i] /= Sum;
    }

    return InitialProb;
}

// ─────────────────────────────────────
double MDP::GetTransProbability(int i, int j) {
    // simplest markov
    if (i + 1 == j) {
        return 1;
    } else {
        return 0;
    }
}

// ─────────────────────────────────────
double MDP::GetOccupancyDistribution(MacroState &State, int u) {
    // u: number of time steps spent in current state
    double time_in_seconds = u * m_BlockDur;
    double expected_duration = m_PsiN1 * State.Duration;
    return std::exp(-time_in_seconds / expected_duration);
}
// ─────────────────────────────────────
int MDP::GetMaxUForJ(MacroState &StateJ) {
    double MaxU = StateJ.Duration / m_BlockDur;
    int MaxUInt = round(MaxU);
    return MaxUInt;
}

// ─────────────────────────────────────
double CalculateEntropy(const std::vector<double> &probs) {
    double entropy = 0.0;
    for (double prob : probs) {
        if (prob > 0) { // Avoid log(0) which is undefined
            entropy -= prob * log(prob);
        }
    }
    return entropy;
}

// ─────────────────────────────────────
double MDP::SemiMarkov(MacroState &StateJ, int CurrentState, int j, int T, int bufferIndex) {
    if (T == 0) {
        double Obs = 0;
        for (AudioState &AudioState : StateJ.AudioStates) {
            Obs = std::max(Obs, AudioState.Obs[bufferIndex]);
        }
        return Obs * GetOccupancyDistribution(StateJ, T + 1) * StateJ.InitProb;
    } else {
        auto bestIt = std::max_element(StateJ.AudioStates.begin(), StateJ.AudioStates.end(),
                                       [bufferIndex](const AudioState &a, const AudioState &b) { return a.Obs[bufferIndex] < b.Obs[bufferIndex]; });
        AudioState BestAudioDescription = *bestIt;
        double Obs = bestIt->Obs[bufferIndex];
        double MaxAlpha = 0;
        for (int u = 1; u <= std::min(T, GetMaxUForJ(StateJ)); u++) {
            double ProbPrevObs = 1.0;
            for (int v = 1; v < u; v++) {
                int PrevIndex = (bufferIndex - v + m_BufferSize) % m_BufferSize;
                ProbPrevObs *= BestAudioDescription.Obs[PrevIndex];
            }
            double Sur = GetOccupancyDistribution(StateJ, u);
            double MaxTrans = 0;
            for (int i = CurrentState; i <= j; i++) {
                if (i < 0) {
                    continue;
                }
                MacroState &StateI = m_States[i];
                int PrevIndex = (T - u) % m_BufferSize;
                if (i != j) {
                    MaxTrans = std::max(MaxTrans, GetTransProbability(i, j) * StateI.Forward[PrevIndex]);
                } else {
                    MaxTrans = std::max(MaxTrans, StateJ.Forward[PrevIndex]);
                }
            }

            double MaxResult = ProbPrevObs * Sur * MaxTrans;
            MaxAlpha = std::max(MaxAlpha, MaxResult);
        }
        return Obs * MaxAlpha;
    }
}

// ─────────────────────────────────────
int MDP::Inference(int CurrentState, int MaxState, int T) {
    double MaxValue = 0;
    int BestState = CurrentState;
    int bufferIndex = T % m_BufferSize;

    for (int j = CurrentState; j <= MaxState; j++) {
        if ((j < 0) || (j >= (int)m_States.size() - 1)) {
            continue;
        }
        MacroState &StateJ = m_States[j];
        StateJ.Forward[bufferIndex] = SemiMarkov(StateJ, CurrentState, j, T, bufferIndex);
    }

    // Sum To Normalize
    double SumForward = 0;
    for (int j = CurrentState; j <= MaxState; j++) {
        if ((j < 0) || (j >= (int)m_States.size() - 1))
            continue;
        SumForward += m_States[j].Forward[bufferIndex];
    }

    // Normalization
    std::vector<double> Probs;
    for (int j = CurrentState; j <= MaxState; j++) {
        if ((j < 0) || ((size_t)j >= m_States.size()))
            continue;

        MacroState &StateJ = m_States[j];
        if (T != 0) {
            double Forward = StateJ.Forward[bufferIndex];
            StateJ.Forward[bufferIndex] = Forward / SumForward;
        }
        if (StateJ.Forward[bufferIndex] > MaxValue) {
            MaxValue = StateJ.Forward[bufferIndex];
            BestState = j;
        }
        Probs.push_back(StateJ.Forward[bufferIndex]);
    }

    return BestState;
}

// ─────────────────────────────────────
int MDP::GetEvent(Description &Desc) {
    m_Desc = Desc;
    m_MaxScoreState = GetMaxJIndex(m_CurrentStateIndex);
    GetAudioObservations(m_CurrentStateIndex - 1, m_MaxScoreState, m_Tau);

    if (Desc.Silence || (size_t)m_CurrentStateIndex == m_States.size()) {
        if (m_CurrentStateIndex == -1) {
            return 0;
        }
        return m_States[m_CurrentStateIndex].ScorePos;
    }

    // just one first event block of the current event
    if (m_Tau == 0) {
        std::vector<double> InitialProb = GetInitialDistribution();
        for (int j = m_CurrentStateIndex; j < m_MaxScoreState; j++) {
            if (j < 0) {
                continue;
            }
            MacroState &StateJ = m_States[j];
            StateJ.InitProb = InitialProb[j - m_CurrentStateIndex];
        }
    }

    int StateIndex = Inference(m_CurrentStateIndex, m_MaxScoreState, m_Tau);
    if (StateIndex == -1) {
        return 0;
    }
    m_PsiN = UpdatePsiN(StateIndex); // Time Model Update

    // Return Score Position
    if (m_CurrentStateIndex == StateIndex) {
        m_CurrentStateIndex = StateIndex;
        return m_States[StateIndex].ScorePos;
    } else {
        m_CurrentStateIndex = StateIndex;
        m_MinEntropy = m_States[StateIndex].Entropy;
        m_SyncStrength = m_States[StateIndex].SyncStrength;
        m_PhaseCoupling = m_States[StateIndex].PhaseCoupling;
        return m_States[StateIndex].ScorePos;
    }
}
} // namespace OScofo
