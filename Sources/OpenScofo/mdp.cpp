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

    std::cout << std::numeric_limits<double>::min() << "\n";
    std::cout << std::numeric_limits<double>::denorm_min() << "\n";
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

    // m_Normalization[t % m_BufferSize] = N(t): the per-timestep sum used in the
    // normalized obs product obs_prod *= b_j(x_{t-v}) / N(t-v) to prevent underflow.
    m_Normalization.resize(m_BufferSize + 1, 1.0); // init to 1 so first-step division is safe

    for (MarkovState &State : m_States) {
        State.Forward.resize(m_BufferSize + 1, 0.0);       // F_j(t)
        State.ExitProb.resize(m_BufferSize + 1, 0.0);      // F_j^o(t)
        State.SumIn_History.resize(m_BufferSize + 1, 0.0); // F_j^i(t)
        State.BestObs.resize(m_BufferSize + 1, 1e-300);    // b_j(x_t), floor avoids /0
        for (AudioState &AS : State.AudioStates) {
            AS.Obs.resize(m_BufferSize + 1, 0.0);
            AS.Forward.resize(m_BufferSize + 1, 0.0);
        }
    }

    m_CurrentStateIndex = 0; // ← was -1: caused GetInitialDistribution to index
                             //   m_States[-1] (UB) and give State 0 an InitProb of 0,
                             //   making Forward[0] = 0 even with perfect observation.
                             //   State 0 is the BEGIN/REST state; we ARE there at start.
    m_Kappa = 1;
    m_BPM = m_States[0].BPMExpected;
    m_PsiN = 60.0f / m_States[0].BPMExpected;
    m_PsiN1 = 60.0f / m_States[0].BPMExpected;
    m_LastPsiN = 60.0f / m_States[0].BPMExpected;
    m_BeatsAhead = m_States[0].BPMExpected / 60 * m_SecondsAhead;
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
void MDP::BuildDistributionCache(double expected_frames) {
    if (expected_frames < 1.0)
        expected_frames = 1.0;

    // Round to 1 decimal place to create a cache key (e.g., 12.4 frames -> 124)
    int key = static_cast<int>(std::round(expected_frames * 10.0));

    // If already computed, skip
    if (m_OccupancyCache.find(key) != m_OccupancyCache.end())
        return;

    int max_u = static_cast<int>(std::ceil(5.0 * expected_frames));
    std::vector<double> occ(max_u + 1, 0.0);
    std::vector<double> surv(max_u + 1, 1.0);

    double p = 0.5;
    double r = expected_frames;

    // Base case for u = 0
    occ[0] = std::pow(p, r);
    surv[0] = 1.0;
    double current_surv = 1.0 - occ[0];

    // Compute recursively (O(1) math per step)
    for (int u = 1; u <= max_u; u++) {
        occ[u] = occ[u - 1] * ((u + r - 1.0) / (double)u) * (1.0 - p);
        surv[u] = (current_surv > 0.0) ? current_surv : 0.0;
        current_surv -= occ[u];
    }

    m_OccupancyCache[key] = std::move(occ);
    m_SurvivorCache[key] = std::move(surv);
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
    auto it = m_KappaCache.find(key);
    if (it != m_KappaCache.end()) {
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
            m_KappaCache[key] = Mid;
            return Mid;
        } else if (A2Mid < SyncStrength) {
            Low = Mid;
        } else {
            High = Mid;
        }
    }

    m_KappaCache[key] = Mid;
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

    // \hat{\phi}_n foi previsto no step anterior e já está guardado aqui
    double HatPhiN = CurrentState.IOIHatPhiN;
    CurrentState.OnsetObserved = m_CurrentStateOnset;

    // Equação 1: Atualização de Kappa
    double PhaseDiff = (IOISeconds / m_PsiN) - HatPhiN;
    double SyncStrength = m_SyncStr - m_SyncStrength * (m_SyncStr - cos(TWO_PI * PhaseDiff));
    double Kappa = InverseA2(SyncStrength);
    m_SyncStr = SyncStrength;
    m_Kappa = Kappa;

    // Equação 2: Atualização de Phi_n
    double FValueUpdate = CouplingFunction(LastPhiN, LastHatPhiN, Kappa);
    double PhiN = LastPhiN + (IOISeconds / m_LastPsiN) + (m_PhaseCoupling * FValueUpdate);
    PhiN = ModPhases(PhiN);

    CurrentState.PhaseObserved = PhiN;
    CurrentState.IOIPhiN = PhiN; // CORREÇÃO: Essencial para o "LastPhiN" do próximo loop

    // Equação 3: Previsão do próximo BPM (\hat{\psi}_{n+1})
    double FValuePrediction = CouplingFunction(PhiN, HatPhiN, Kappa);
    double PsiN1 = m_PsiN * (1 + m_SyncStrength * FValuePrediction);

    // Equação 4: Previsão do próximo \hat{\phi}_{n+1}
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
        // Keep this for now
        StateJ.TimeProb = 1;

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

    int key = static_cast<int>(std::round(expected_frames * 10.0));
    BuildDistributionCache(expected_frames); // Builds only if missing

    if (u < static_cast<int>(m_OccupancyCache[key].size())) {
        return m_OccupancyCache[key][u];
    }
    return 0.0; // u is beyond the tail
}

// ─────────────────────────────────────
double MDP::GetSurvivorDistribution(MarkovState &State, int u) {
    double expected_frames = (m_PsiN1 * State.Duration) / m_BlockDur;
    if (expected_frames < 1.0)
        expected_frames = 1.0;

    int key = static_cast<int>(std::round(expected_frames * 10.0));
    BuildDistributionCache(expected_frames); // Builds only if missing

    if (u < static_cast<int>(m_SurvivorCache[key].size())) {
        return m_SurvivorCache[key][u];
    }
    return 0.0; // u is beyond the tail
}

// ─────────────────────────────────────
int MDP::GetMaxUForJ(MarkovState &StateJ) {
    double expected_frames = (m_PsiN1 * StateJ.Duration) / m_BlockDur;
    if (expected_frames < 1.0)
        expected_frames = 1.0;

    // For Negative Binomial, the tail is longer.
    // 5x expected_frames ensures we capture the "ageing" properties
    // mentioned in your text without truncating the probability peak.
    int cap = static_cast<int>(std::ceil(5.0 * expected_frames));
    return std::max(cap, 1);
}

// ─────────────────────────────────────
void MDP::Markov(MarkovState &StateJ, int j, int T, int bufferIndex) {
    double bj = StateJ.BestObs[bufferIndex];

    double fj;
    if (T == 0) {
        fj = bj * StateJ.InitProb;
    } else {
        int prevBuf = (bufferIndex - 1 + m_BufferSize) % m_BufferSize;
        double sumPrev = 0.0;

        // Stay in j (self-loop)
        if (j >= m_CurrentStateIndex) {
            sumPrev += StateJ.Forward[prevBuf];
        }
        // Arrive from j-1
        if (j - 1 >= m_CurrentStateIndex && j - 1 >= 0) {
            double trans = GetTransProbability(j - 1, j);
            sumPrev += trans * m_States[j - 1].ExitProb[prevBuf];
        }
        fj = bj * sumPrev;
    }

    // For Markov states Forward = ExitProb (they can exit at every step)
    StateJ.Forward[bufferIndex] = fj;
    StateJ.ExitProb[bufferIndex] = fj;
}

// ─────────────────────────────────────
void MDP::SemiMarkov(MarkovState &StateJ, int j, int T, int bufferIndex) {
    double bj = StateJ.BestObs[bufferIndex]; // Factored out b_j(o_t)

    double f_tilde_j = 0.0;  // \tilde{f}_j(t)
    double f_tilde_jo = 0.0; // \tilde{f}_j^o(t)
    double ObsProd = 1.0;    // b_j(o_{t-u+1}^{t-1})

    int maxU = GetMaxUForJ(StateJ);

    for (int u = 1; u <= T + 1; u++) {
        double D_bar_ju = GetSurvivorDistribution(StateJ, u); // \overline{D}_j(u)
        double d_ju = GetOccupancyDistribution(StateJ, u);    // d_j(u)

        if (u == T + 1) {
            // INIT TERM: \overline{D}_j(t+1) * b_j(o_0^t) * \pi(j)
            f_tilde_j += D_bar_ju * ObsProd * StateJ.InitProb;
            f_tilde_jo += d_ju * ObsProd * StateJ.InitProb;
            break;
        }

        if (u <= maxU) {
            int entryBuf = ((T - u) % m_BufferSize + m_BufferSize) % m_BufferSize;

            double transition_sum = 0.0;
            for (int i = 0; i < (int)m_States.size(); i++) {
                if (i != j) {
                    double p_ij = GetTransProbability(i, j);
                    if (p_ij > 0.0) {
                        transition_sum += p_ij * m_States[i].ExitProb[entryBuf];
                    }
                }
            }

            f_tilde_j += D_bar_ju * ObsProd * transition_sum;
            f_tilde_jo += d_ju * ObsProd * transition_sum;
        }

        // Advance scaled observation product to prevent floating-point underflow
        int prevBuf = ((T - u) % m_BufferSize + m_BufferSize) % m_BufferSize;
        double prevObs = StateJ.BestObs[prevBuf];
        double prevNorm = std::max(m_Normalization[prevBuf], 1e-300);
        ObsProd *= prevObs / prevNorm;

        if (ObsProd < 1e-15)
            break;
    }

    // Apply the factored-out current observation b_j(o_t)
    StateJ.Forward[bufferIndex] = bj * (f_tilde_j + 1e-300);
    StateJ.ExitProb[bufferIndex] = bj * (f_tilde_jo + 1e-300);
}

// ─────────────────────────────────────
int MDP::Inference(int T) {
    int bIndex = T % m_BufferSize;
    spdlog::debug("WinStart {} | WinFinish {} | BufferSize {} | Tau {} | Kappa {}", m_WinStart, m_WinEnd, bIndex, m_Tau,
                  m_Kappa);

    // ── STEP 1: Compute F_j(t) and F_j^o(t) for all states ─────────
    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        if (j < 0 || j >= (int)m_States.size())
            continue;
        MarkovState &StateJ = m_States[j];
        if (StateJ.HSMMType == SEMIMARKOV)
            SemiMarkov(StateJ, j, T, bIndex); // sets Forward and ExitProb
        else
            Markov(StateJ, j, T, bIndex); // sets Forward and ExitProb
    }

    double N = 0.0;
    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        if (j < 0 || j >= (int)m_States.size())
            continue;
        N += m_States[j].Forward[bIndex];
    }
    if (N < 1e-300)
        N = 1e-300;
    m_Normalization[bIndex] = N;

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        if (j < 0 || j >= (int)m_States.size())
            continue;
        m_States[j].Forward[bIndex] /= N;
        m_States[j].Forward[bIndex] += 1e-300; // floor (reference C code line)
        m_States[j].ExitProb[bIndex] /= N;
        m_States[j].ExitProb[bIndex] += 1e-300;
    }

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        if (j < 0 || j >= (int)m_States.size())
            continue;
        double entry = 0.0;
        if (j > 0 && (j - 1) >= 0 && (j - 1) < (int)m_States.size()) {
            entry = GetTransProbability(j - 1, j) * m_States[j - 1].ExitProb[bIndex];
        }
        m_States[j].SumIn_History[bIndex] = entry;
    }

    double maxVal = 1e-300;
    int bestStateIndex = m_CurrentStateIndex;
    for (int j = m_CurrentStateIndex; j <= m_WinEnd; ++j) {
        if (j < 0 || j >= (int)m_States.size())
            continue;
        MarkovState &StateJ = m_States[j];
        double fwd = StateJ.Forward[bIndex];
        if (fwd > maxVal) {
            maxVal = fwd;
            bestStateIndex = j;
        }
        spdlog::debug("State ({}) | Obs = {:.5f}, Forward {:.5f}, Exit Prob {:.5f}", StateJ.Index,
                      StateJ.BestObs[bIndex], StateJ.Forward[bIndex], StateJ.ExitProb[bIndex]);
    }

    MarkovState &BestState = m_States[bestStateIndex];
    spdlog::debug("Best: State ({}) | Forward {:.5f}", BestState.Index, BestState.Forward[bIndex]);

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

    if (m_T == 0) {
        GetInitialDistribution();
    }

    // Run forward inference
    int BestState = Inference(m_Tau);
    m_PsiN = UpdatePsiN(BestState);

    // Advance the score position if a new event was detected
    if (BestState != m_CurrentStateIndex) {
        spdlog::info("New Event Index {:04d}, Score Position {:04d}", BestState, m_States[BestState].ScorePos);
        m_CurrentStateIndex = BestState;
    }

    printf("\n");
    return m_States[BestState].ScorePos;
}

} // namespace OpenScofo
