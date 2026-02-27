#include "OpenScofo.hpp"

namespace OpenScofo {

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
    m_EventWindowSize = 20;
    SetTunning(440);

    for (int i = 0; i <= 10000; ++i) { // 0.0000, 0.0001, 0.0002 ...
        double key = i / 10000.0;
        InverseA2(key);
    }
}

// ─────────────────────────────────────
ActionVec MDP::GetEventActions(int Index) {
    if (Index < 0 || Index >= (int)m_States.size()) {
        return ActionVec();
    }

    MarkovState State = m_States[(size_t)Index];
    return State.Actions;
}

// ─────────────────────────────────────
void MDP::SetScoreStates(States ScoreStates) {
    if (ScoreStates.size() == 0) {
        return;
    }

    m_States.clear();
    m_States = ScoreStates;

    m_Normalization.resize(m_BufferSize + 1, 0.0);

    for (MarkovState &State : m_States) {
        State.Forward.resize(m_BufferSize + 1, 0.0);
        State.Forward.resize(m_BufferSize + 1, 0.0); // ADD THIS
        State.BestObs.resize(m_BufferSize + 1, 0.0);
        for (AudioState &AS : State.AudioStates) {
            AS.Obs.resize(m_BufferSize + 1, 0.0);
            AS.Forward.resize(m_BufferSize + 1, 0.0);
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
void MDP::BuildPitchTemplate(double Freq) {
    const double m_MinHarmonicDecay = 0.2;
    const double m_MaxHarmonicDecay = 1.8;
    const double binWidth = m_Sr / m_FFTSize;

    if (Freq > m_Sr / 2.0 || Freq <= 0)
        return;

    double m_MinF0 = 32.70;
    double m_MaxF0 = 4186.0;
    double f0Norm = std::log2(Freq / m_MinF0) / std::log2(m_MaxF0 / m_MinF0);
    f0Norm = std::clamp(f0Norm, 0.0, 1.0);

    double pitchDecayCurve = m_PitchTemplateAmplitudeDecay;
    double shaped = std::pow(1.0 - f0Norm, pitchDecayCurve);
    double beta = m_MaxHarmonicDecay - shaped * (m_MaxHarmonicDecay - m_MinHarmonicDecay);

    double rootBinFreq = round(Freq / binWidth);

    if (m_PitchTemplates.find(rootBinFreq) != m_PitchTemplates.end()) {
        return;
    }

    std::vector<double> templateBins(m_FFTSize / 2, 1e-24);

    const double sigmaLog = m_PitchTemplateSigma / 12.0;
    const double sigmaConst = std::pow(2.0, sigmaLog) - 1.0;

    const double B = 0.0001;
    for (int k = 1; k <= m_Harmonics; ++k) {
        double stretch = std::sqrt(1.0 + B * k * k);
        double harmonicFreqHz = Freq * k * stretch;
        if (harmonicFreqHz >= m_Sr / 2.0)
            break;

        double sigmaHz = harmonicFreqHz * sigmaConst;
        if (sigmaHz < binWidth * 0.75) {
            sigmaHz = binWidth * 0.75;
        }
        double envelope = std::exp(-beta * (k - 1));
        if (k > 1) {
            envelope *= 1.25;
        }
        if (envelope < 1e-5)
            break;

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
double MDP::GetBlockDuration() {
    return m_BlockDur;
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
    spdlog::debug("Current event is {}", Event);
    if (m_States.size() == 0) {
        spdlog::error("There is not Events on Score or the Score was no loaded");
        return;
    }

    if (Event == 0) {
        spdlog::info("Initializing Time Decoding Algorithm");
        InitTimeDecoding();
    }
    m_CurrentStateIndex = Event;
    m_T = 0;
    m_Tau = 0;
    m_TauWithSound = 0;
    std::cout << "\n" << std::endl;
}

// ─────────────────────────────────────
int MDP::GetStatesSize() {
    return m_States.size();
}

// ─────────────────────────────────────
void MDP::AddState(MarkovState State) {
    m_States.push_back(State);
}

// ─────────────────────────────────────
MarkovState MDP::GetState(int Index) {
    return m_States[Index];
}

// ─────────────────────────────────────
void MDP::SetPitchTemplateSigma(double f) {
    m_PitchTemplateSigma = f;
}

// ╭─────────────────────────────────────╮
// │            Time Decoding            │
// ╰─────────────────────────────────────╯
void MDP::InitTimeDecoding(void) {
    double PsiK = 60 / m_States[0].BPMExpected;
    m_LastPsiN = PsiK;
    m_PsiN = PsiK;
    m_PsiN1 = PsiK;
    m_States[0].OnsetObserved = 0;
    m_BPM = m_States[0].BPMExpected;
    m_CurrentStateOnset = 0;
    m_LastTn = 0;
    m_TimeInPrevEvent = 0;
    m_T = 0;
}

// ─────────────────────────────────────
double MDP::InverseA2(double SyncStrength) {
    if (SyncStrength < 0) {
        return 0;
    }

    if (SyncStrength >= 0.95) {
        return 10.0f;
    }

    double key = std::round(SyncStrength * 10000.0) / 10000.0;
    auto it = m_KappaArray.find(key);
    if (it != m_KappaArray.end()) {
        return it->second;
    }

    double Low = 0.0;
    double Tol = 1e-10;
    double High = std::max(SyncStrength, 11.0);
    double Mid;

    int i;
    double A2Mid;
    for (i = 0; i < 8192; i++) {
        Mid = (Low + High) / 2.0;
        double I1 = std::cyl_bessel_i(1, Mid);
        double I0 = std::cyl_bessel_i(0, Mid);
        A2Mid = I1 / I0;
        if (std::fabs(A2Mid - SyncStrength) < Tol) {
            m_KappaArray[key] = Mid;
            return Mid;
        } else if (A2Mid < SyncStrength) {
            Low = Mid;
        } else {
            High = Mid;
        }
    }

    m_KappaArray[key] = Mid;
    spdlog::warn("Kappa for Sync Strength {:.4f} not converged after {} iterations. Returning {:.4f} for A2 = {:.4f}",
                 SyncStrength, i, Mid, A2Mid);
    return Mid;
}

// ─────────────────────────────────────
double MDP::CouplingFunction(double Phi, double PhiMu, double Kappa) {
    double ExpKappa = exp(Kappa);
    double PhiNDiff = Phi - PhiMu;
    double CosTerm = cos(TWO_PI * PhiNDiff);
    double SinTerm = sin(TWO_PI * PhiNDiff);
    double PhiN = (1 / (TWO_PI * ExpKappa)) * exp(Kappa * CosTerm) * SinTerm;
    return PhiN;
}

// ─────────────────────────────────────
double MDP::ModPhases(double Phase) {
    Phase = std::fmod(Phase + M_PI, TWO_PI);
    if (Phase < 0) {
        Phase += TWO_PI;
    }
    return Phase - M_PI;
}

// ─────────────────────────────────────
States MDP::GetStatesForProcessing() {
    double EventOnset = m_States[m_CurrentStateIndex].Duration - (m_TimeInPrevEvent + m_BlockDur);
    size_t begin = m_CurrentStateIndex;
    size_t end = begin;

    for (; end < m_States.size(); ++end) {
        if (EventOnset > m_SecondsAhead)
            break;

        EventOnset += m_States[end].Duration;
    }

    return States(m_States.begin() + begin, m_States.begin() + end);
}

// ─────────────────────────────────────
void MDP::GetDecodeWindow() {
    int half = m_EventWindowSize / 2;
    m_WinStart = std::max(0, (int)m_CurrentStateIndex - half);
    m_WinEnd = std::min((int)m_States.size() - 1, (int)m_CurrentStateIndex + half);
}

// ─────────────────────────────────────
double MDP::UpdatePsiN(int StateIndex) {
    m_Tau += 1;
    m_T += 1;

    if (StateIndex == m_CurrentStateIndex) {
        m_TimeInPrevEvent += m_BlockDur;
        return m_PsiN;
    } else {
        m_TimeInPrevEvent += m_BlockDur;
        m_LastTn = m_CurrentStateOnset;
        m_CurrentStateOnset += m_TimeInPrevEvent;
    }

    MarkovState &LastState = m_States[StateIndex - 1];
    MarkovState &CurrentState = m_States[StateIndex];
    MarkovState &NextState = m_States[StateIndex + 1];

    double IOISeconds = m_CurrentStateOnset - m_LastTn;
    double LastPhiN = LastState.IOIPhiN;
    double LastHatPhiN = LastState.IOIHatPhiN;
    double HatPhiN = CurrentState.IOIHatPhiN;
    double PhiNExpected = LastPhiN + ((m_CurrentStateOnset - m_LastTn) / m_PsiN);
    CurrentState.IOIHatPhiN = PhiNExpected;
    CurrentState.OnsetObserved = m_CurrentStateOnset;

    double PhaseDiff = (IOISeconds / m_PsiN) - HatPhiN;
    double SyncStrength = m_SyncStr - m_SyncStrength * (m_SyncStr - cos(TWO_PI * PhaseDiff));
    double Kappa = InverseA2(SyncStrength);
    m_SyncStr = SyncStrength;
    m_Kappa = Kappa;

    double FValueUpdate = CouplingFunction(LastPhiN, LastHatPhiN, Kappa);
    double PhiN = LastPhiN + (IOISeconds / m_LastPsiN) + (m_PhaseCoupling * FValueUpdate);
    PhiN = ModPhases(PhiN);
    CurrentState.PhaseObserved = PhiN;

    double FValuePrediction = CouplingFunction(PhiN, HatPhiN, Kappa);
    double PsiN1 = m_PsiN * (1 + m_SyncStrength * FValuePrediction);

    double Tn1 = m_CurrentStateOnset + CurrentState.Duration * PsiN1;
    double PhiN1 = ModPhases((Tn1 - m_CurrentStateOnset) / PsiN1);
    NextState.IOIHatPhiN = PhiN1;
    NextState.OnsetExpected = Tn1;

    double LastOnsetExpected = Tn1;
    for (int i = m_CurrentStateIndex + 2; i < m_CurrentStateIndex + 20; i++) {
        if (i >= static_cast<int>(m_States.size()))
            break;
        MarkovState &FutureState = m_States[i];
        MarkovState &PreviousFutureState = m_States[(i - 1)];
        double Duration = PreviousFutureState.Duration;
        double FutureOnset = LastOnsetExpected + Duration * PsiN1;
        FutureState.OnsetExpected = FutureOnset;
        LastOnsetExpected = FutureOnset;
    }

    m_BPM = 60.0f / m_PsiN;
    m_LastPsiN = m_PsiN;

    if (StateIndex != m_CurrentStateIndex) {
        m_TimeInPrevEvent = 0;
        m_T = 0;
    }
    return PsiN1;
}

// ╭─────────────────────────────────────╮
// │     Markov / Semi-Markov Core       │
// ╰─────────────────────────────────────╯
void MDP::GetAudioObservations(int T) {
    std::unordered_map<double, double> PitchObs;
    int bufferIndex = T % m_BufferSize;
    double ObsNoSound = 0;
    double ObsSilence = m_Desc.SilenceProb;

    for (int j = m_WinStart; j <= m_WinEnd; j++) {
        if (j < 0 || j >= (int)m_States.size())
            continue;

        MarkovState &StateJ = m_States[j];
        double BestObs = 1e-300;

        // Time
        {
            int currentBlock = m_TauWithSound;
            int expectedBlock = static_cast<int>(std::round(StateJ.OnsetExpected / m_BlockDur));
            double distance = static_cast<double>(currentBlock - expectedBlock);
            double sigma = 200.0;
            double prob = std::exp(-0.5 * (distance * distance) / (sigma * sigma));
            double eps = 0.05; // piso mínimo
            double beta = 0.46;
            double factor = eps + (1.0 - eps) * (1.0 - std::exp(-beta * m_Kappa)) / (1.0 - std::exp(-beta * 10.0));
            StateJ.TimeProb = prob * factor;
        }

        if (StateJ.Type == NOTE || StateJ.Type == TRILL) {
            for (AudioState &AS : StateJ.AudioStates) {
                switch (AS.Type) {
                case PITCH: {
                    auto it = PitchObs.find(AS.Freq);
                    if (it != PitchObs.end()) {
                        AS.Obs[bufferIndex] = it->second;
                        BestObs = std::max(BestObs, it->second);
                    } else {
                        double kl = GetPitchSimilarity(AS.Freq) * (1.0 - m_Desc.SilenceProb);
                        PitchObs[AS.Freq] = kl;
                        AS.Obs[bufferIndex] = kl;
                        BestObs = std::max(BestObs, kl);
                    }
                    break;
                }
                default:
                    spdlog::error("Not implemented yet!");
                }
            }
            ObsNoSound = std::max(ObsNoSound, BestObs);

        } else if (StateJ.Type == REST) {
            for (AudioState &AS : StateJ.AudioStates) {
                AS.Obs[bufferIndex] = m_Desc.SilenceProb;
                BestObs = std::max(BestObs, m_Desc.SilenceProb);
            }
        }

        StateJ.BestObs[bufferIndex] = BestObs * StateJ.TimeProb;
    }

    if (ObsNoSound > ObsSilence) {
        spdlog::debug("SOUND   | Sound {} | Silence {}", ObsNoSound, ObsSilence);
        m_TauWithSound++;
    } else {
        spdlog::debug("SILENCE | Sound {} | Silence {}", ObsNoSound, ObsSilence);
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
        double P = PitchTemplate[i] + m_Desc.ReverbSpectralPower[i];
        double Q = m_Desc.NormSpectralPower[i];
        if (P > 0 && Q > 0) {
            KLDiv += P * log(P / Q);
        } else if (P == 0 && Q >= 0) {
            KLDiv += Q;
        }
    }

    double noise_robustness = 1.0 / (1.0 + m_Desc.StdDev);
    KLDiv *= noise_robustness;
    KLDiv = exp(-m_PitchScalingFactor * KLDiv);
    return KLDiv;
}

// ─────────────────────────────────────
void MDP::GetInitialDistribution() {
    int Size = m_WinEnd - m_CurrentStateIndex + 1;
    std::vector<double> InitialProb(Size);

    double Dur = 0;
    double Sum = 0;

    for (int i = 0; i < Size; i++) {
        double DurProb = exp(-1.0 * (Dur / m_BeatsAhead));
        InitialProb[i] = DurProb;
        Dur += m_States[m_CurrentStateIndex + i].Duration;
        Sum += DurProb;
    }

    if (Sum > 1e-12) {
        for (int i = 0; i < Size; i++) {
            InitialProb[i] /= Sum;
        }
    }

    for (int j = m_WinStart; j <= m_WinEnd; j++) {
        if (j < 0 || j >= (int)m_States.size())
            continue;
        int idx = j - m_CurrentStateIndex;
        MarkovState &StateJ = m_States[j];
        StateJ.InitProb = (idx < (int)InitialProb.size()) ? InitialProb[idx] : 0.0;
        if (j < m_CurrentStateIndex) {
            m_States[j].InitProb = 0.0;
        }
    }

    return;
}

// ─────────────────────────────────────
double MDP::GetTransProbability(int i, int j) {
    return (i + 1 == j) ? 1.0 : 0.0;
}

// ─────────────────────────────────────
double MDP::GetOccupancyDistribution(MarkovState &State, int u) {
    double expected_frames = (m_PsiN1 * State.Duration) / m_BlockDur;
    if (expected_frames < 1.0)
        expected_frames = 1.0;
    return std::exp(-static_cast<double>(u) / expected_frames);
}

// ─────────────────────────────────────
int MDP::GetMaxUForJ(MarkovState &StateJ) {
    double expected_frames = (m_PsiN1 * StateJ.Duration) / m_BlockDur;
    int cap = static_cast<int>(std::ceil(3.0 * expected_frames));
    return std::max(cap, 1);
}

// ─────────────────────────────────────
double MDP::Markov(MarkovState &StateJ, int j, int T, int bufferIndex) {
    double bj = StateJ.BestObs[bufferIndex];
    if (T == 0) {
        return bj * StateJ.InitProb;
    }
    int prevBuf = (bufferIndex - 1 + m_BufferSize) % m_BufferSize;
    double bestPrev = 1e-300;

    if (j >= m_CurrentStateIndex) {
        bestPrev = std::max(bestPrev, StateJ.Forward[prevBuf]);
    }

    if (j - 1 >= m_CurrentStateIndex && j - 1 >= 0) {
        double trans = GetTransProbability(j - 1, j);
        bestPrev = std::max(bestPrev, trans * m_States[j - 1].Forward[prevBuf]);
    }

    return bj * bestPrev;
}

// ─────────────────────────────────────
double MDP::SemiMarkov(MarkovState &StateJ, int j, int T, int bufferIndex) {
    double bj = StateJ.BestObs[bufferIndex];

    if (m_Tau == 0) {
        return bj * GetOccupancyDistribution(StateJ, 1) * StateJ.InitProb;
    }

    int maxU = std::min(T, GetMaxUForJ(StateJ));
    double bestU = 0.0;

    for (int u = 1; u <= maxU; u++) {
        double prod = 1.0;
        for (int v = 1; v < u; v++) {
            int PrevIndex = (bufferIndex - v + m_BufferSize) % m_BufferSize;
            prod *= StateJ.BestObs[PrevIndex];
        }

        double dju = GetOccupancyDistribution(StateJ, u);

        // as we use linear semi-Markov chains: ∀i, j, pij = δi,i+1.
        int iPrev = j - 1;
        if (iPrev >= 0) {
            int entryBuf = (T - u + m_BufferSize) % m_BufferSize;
            double raw = m_States[iPrev].Forward[entryBuf];
            double bestEntry = GetTransProbability(iPrev, j) * raw;
            double candidate = prod * dju * bestEntry;
            bestU = std::max(bestU, candidate);
        }
    }

    return bj * (bestU + 1e-300);
}

// ─────────────────────────────────────
int MDP::Inference(int T) {
    int bufferIndex = T % m_BufferSize;
    spdlog::debug("WinStart {} | WinFinish {} | BufferSize {} | Tau {} | Kappa {}", m_WinStart, m_WinEnd, bufferIndex,
                  m_Tau, m_Kappa);

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        if (j < 0 || j >= (int)m_States.size())
            continue;

        MarkovState &StateJ = m_States[j];
        if (StateJ.HSMMType == SEMIMARKOV)
            StateJ.Forward[bufferIndex] = SemiMarkov(StateJ, j, T, bufferIndex);
        else
            StateJ.Forward[bufferIndex] = Markov(StateJ, j, T, bufferIndex);
    }

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        if (j < 0 || j >= (int)m_States.size())
            continue;
        m_States[j].ForwardLast = m_States[j].Forward[bufferIndex];
    }

    double maxVal = 1e-300;
    int bestStateIndex = m_CurrentStateIndex;
    for (int j = m_CurrentStateIndex; j <= m_WinEnd; ++j) {
        if (j < 0 || j >= (int)m_States.size())
            continue;

        MarkovState StateJ = m_States[j];
        double fwd = StateJ.Forward[bufferIndex];
        if (fwd > maxVal) {
            maxVal = fwd;
            bestStateIndex = j;
        }
        spdlog::debug("State ({}) | ForwardLast {:.5f}, Obs = {:.5f}, Forward {:.5f}, Time Prob {:.5f}", StateJ.Index,
                      StateJ.ForwardLast, StateJ.BestObs[bufferIndex], StateJ.Forward[bufferIndex], StateJ.TimeProb);
    }

    MarkovState BestState = m_States[bestStateIndex];
    spdlog::debug("State ({}) | ForwardLast {:.5f}, Obs = {:.5f}, Forward {:.5f}", BestState.Index,
                  BestState.ForwardLast, BestState.BestObs[bufferIndex], BestState.Forward[bufferIndex]);

    return bestStateIndex;
}

// ─────────────────────────────────────
int MDP::GetEvent(Description &Desc) {

    spdlog::debug("Starting inference");
    m_Desc = Desc;
    if (m_CurrentStateIndex > (int)m_States.size()) {
        spdlog::debug("Score Finished");
        return m_States.back().ScorePos;
    }

    GetDecodeWindow();
    GetAudioObservations(m_Tau);

    // ── On the very first frame of a new event: set initial distribution ─────
    if (m_T == 0) {
        // TODO:: Save this inside the loop
        GetInitialDistribution();
    }

    // ── Run forward inference ────────────────────────────────────────────────
    int BestState = Inference(m_Tau);
    m_PsiN = UpdatePsiN(BestState); // update tempo model

    // ── Advance the score position if a new event was detected ───────────────
    if (BestState != m_CurrentStateIndex) {
        spdlog::info("New Event Index {:04d}, Score Position {:04d}", BestState, m_States[BestState].ScorePos);
        m_CurrentStateIndex = BestState;
        m_MinEntropy = m_States[BestState].Entropy;
        m_SyncStrength = m_States[BestState].SyncStrength;
        m_PhaseCoupling = m_States[BestState].PhaseCoupling;
    }

    printf("\n");
    return m_States[BestState].ScorePos;
}

} // namespace OpenScofo
