#include "OpenScofo.hpp"

#if defined(__APPLE__)
#include <boost/math/special_functions/bessel.hpp>
#define CYL_BESSEL_I(v, x) boost::math::cyl_bessel_i(v, x)
#else
#include <cmath>
#define CYL_BESSEL_I(v, x) std::cyl_bessel_i(v, x)
#endif

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
    m_SyncStrength = 0.5;
    m_PhaseCoupling = 0.5;
    m_TimeInPrevEvent = 0;
    m_EventWindowSize = 20;
    SetTunning(440);

    constexpr int KappaPrecision = 1000;
    m_KappaCache.reserve(static_cast<size_t>(KappaPrecision + 1));
    for (int i = 0; i <= KappaPrecision; i++) {
        double key = i / 1000.0;
        InverseA2(key);
    }

    UpdateAudioParameters(Sr, FFTSize, HopSize);
}

// ─────────────────────────────────────
void MDP::UpdateAudioParameters(double Sr, double FFTSize, double HopSize) {
    m_HopSize = HopSize;
    m_FFTSize = FFTSize;
    m_Sr = Sr;
    m_BlockDur = (1.0 / m_Sr) * m_HopSize;

    m_PitchTemplates.clear();
    m_PitchCQTTemplates.clear();
    m_OccupancyCache.clear();
    m_SurvivorCache.clear();

    m_Normalization.assign(static_cast<size_t>(m_BufferSize + 1), 1.0);
    for (MarkovState &State : m_States) {
        State.Forward.assign(static_cast<size_t>(m_BufferSize + 1), 0.0);
        State.ExitProb.assign(static_cast<size_t>(m_BufferSize + 1), 0.0);
        State.BestObs.assign(static_cast<size_t>(m_BufferSize + 1), 1e-300);
    }

    if (!m_States.empty()) {
        UpdateAudioTemplate();
        UpdatePhaseValues();
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
    spdlog::debug("There is {} score states", m_States.size());

    m_Normalization.resize(m_BufferSize + 1, 1.0); // init to 1 so first-step division is safe
    for (MarkovState &State : m_States) {
        State.Forward.resize(m_BufferSize + 1, 0.0);    // F_j(t)
        State.ExitProb.resize(m_BufferSize + 1, 0.0);   // F_j^o(t)
        State.BestObs.resize(m_BufferSize + 1, 1e-300); // b_j(x_t), floor avoids /0
    }

    m_CurrentStateIndex = 0;
    m_Kappa = 10;
    m_BPM = m_States[0].BPMExpected;
    m_PsiN = 60.0f / m_States[0].BPMExpected;
    m_PsiN1 = 60.0f / m_States[0].BPMExpected;
    m_LastPsiN = 60.0f / m_States[0].BPMExpected;
    m_BeatsAhead = m_States[0].BPMExpected / 60 * m_SecondsAhead;
    m_SyncStr = 0;
    if (std::isfinite(m_States[0].SyncStrength)) {
        m_SyncStrength = std::clamp(m_States[0].SyncStrength, 0.0, 1.0);
    }
    if (std::isfinite(m_States[0].PhaseCoupling)) {
        m_PhaseCoupling = std::clamp(m_States[0].PhaseCoupling, 0.0, 2.0);
    }

    UpdateAudioTemplate();
    UpdatePhaseValues();
}

// ─────────────────────────────────────
// GONG 2015 (adapted)
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
        double normalizationFactor = 1.0 / (sigmaHz * std::sqrt(2.0 * std::numbers::pi));

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
// GONG 2015 (adapted)
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
// GONG 2015 (adapted)
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
    m_Tau = 0;
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
// CONT 2010
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
}

// ─────────────────────────────────────
// https://stat.ethz.ch/R-manual//R-patched/library/stats/html/NegBinomial.html
// CUVILLIER 2016 (Check chapter 3 and 4)
void MDP::BuildDistributionCache(double ExpectedFrames) {
    if (ExpectedFrames < 1.0)
        ExpectedFrames = 1.0;

    int key = static_cast<int>(std::round(ExpectedFrames * 10.0));

    if (m_OccupancyCache.find(key) != m_OccupancyCache.end())
        return;

    int maxU = static_cast<int>(std::ceil(5.0 * ExpectedFrames));
    std::vector<double> occ(maxU + 1, 0.0);
    std::vector<double> surv(maxU + 1, 1.0);

    if (occ.empty() || surv.empty())
        return;

    double p = 0.5;
    double r = ExpectedFrames;
    double currentSurv = 1.0;

    // Log-space constants
    double log_p = std::log(p);
    double log_1_p = std::log(1.0 - p);

    // Base case ln(P(0)) = r * ln(p)
    double current_log_occ = r * log_p;

    for (int u = 0; u <= maxU; u++) {
        // Safely convert log to linear space. Exp(<-708) safely yields 0.0
        occ[u] = std::exp(current_log_occ);

        surv[u] = (currentSurv > 0.0) ? currentSurv : 0.0;
        currentSurv -= occ[u];

        // Advance to next step in log space:
        // ln(P(u+1)) = ln(P(u)) + ln(u + r) - ln(u + 1) + ln(1 - p)
        current_log_occ += std::log(u + r) - std::log(u + 1.0) + log_1_p;
    }

    m_OccupancyCache[key] = std::move(occ);
    m_SurvivorCache[key] = std::move(surv);
}

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
double MDP::A2(double kappa) {
    if (kappa <= 0.0) {
        return 0.0;
    }

    if (kappa > 10.0) {
        return 1.0 - (1.0 / (2.0 * kappa)) - (1.0 / (8.0 * kappa * kappa));
    }

    double I1 = CYL_BESSEL_I(1, kappa);
    double I0 = CYL_BESSEL_I(0, kappa);
    if (!std::isfinite(I1) || !std::isfinite(I0) || I0 <= 0.0) {
        return 1.0 - (1.0 / (2.0 * kappa));
    }
    return I1 / I0;
}

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
double MDP::InverseA2(double SyncStrength) {
    // SyncStrength must be between 0 and 1
    if (SyncStrength < 0) {
        return 0;
    }

    // Following Large and Jones (1999, p. 157).
    if (SyncStrength > 0.95) {
        return 10.0f;
    }

    // Constructor pre-fills this cache at 1e-3 resolution.
    const int cacheKey = static_cast<int>(std::round(std::clamp(SyncStrength, 0.0, 1.0) * 1000.0));
    auto cacheIt = m_KappaCache.find(cacheKey);
    if (cacheIt != m_KappaCache.end()) {
        return cacheIt->second;
    }

    double Low = 0.0;
    double Tol = 1e-16;
    double High = std::max(SyncStrength, 10.0);
    double Mid;

    // In my tests I never reached more than 100 iterations.
    int i;
    for (i = 0; i < 1000; ++i) {
        Mid = (Low + High) / 2.0;
        double I1 = CYL_BESSEL_I(1, Mid);
        double I0 = CYL_BESSEL_I(0, Mid);
        double A2Mid = I1 / I0;
        if (std::fabs(A2Mid - SyncStrength) < Tol) {
            m_KappaCache.emplace(cacheKey, Mid);
            return Mid;
        } else if (A2Mid < SyncStrength) {
            Low = Mid;
        } else {
            High = Mid;
        }
    }
    // LOGE() << "InverseA2 not converged after " << i << " iterations.";
    m_KappaCache.emplace(cacheKey, Mid);
    return Mid;
}

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
double MDP::CouplingFunction(double phi, double phi_hat, double kappa) {
    static constexpr double invTwoPi = 1.0 / (2.0 * std::numbers::pi);
    double diff = 2.0 * std::numbers::pi * (phi - phi_hat);
    double cosDiff = std::cos(diff);
    return invTwoPi * std::exp(kappa * (cosDiff - 1.0)) * std::sin(diff);
}

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
double MDP::ModPhases(double Phase) {
    Phase = std::fmod(Phase + 0.5, 1.0);
    if (Phase < 0.0)
        Phase += 1.0;
    return Phase - 0.5;
}

// ─────────────────────────────────────
// CONT 2010 (Last § of section 4)
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
// CONT 2010 (Last § of section 4)
void MDP::GetDecodeWindow() {
    int half = m_EventWindowSize / 2;
    m_WinStart = std::max(0, static_cast<int>(m_CurrentStateIndex) - half);
    m_WinEnd = std::min(static_cast<int>(m_States.size()) - 1, static_cast<int>(m_CurrentStateIndex) + half);

    if (m_WinStart < 0 || m_WinEnd >= static_cast<int>(m_States.size())) {
        spdlog::critical("MDP::GetDecodeWindow invariant violated: "
                         "window out of bounds "
                         "(winStart={}, winEnd={}, statesSize={}, currentIndex={}, eventWindowSize={})",
                         m_WinStart, m_WinEnd, m_States.size(), m_CurrentStateIndex, m_EventWindowSize);
    }
}

// ─────────────────────────────────────
// CONT 2010 (Section 5, algorithm 1)
double MDP::UpdatePsiN(int StateIndex) {
    m_TimeInPrevEvent += m_BlockDur;
    m_Tau += 1;

    if (StateIndex == m_CurrentStateIndex) {
        m_PsiN1 = m_PsiN;
        return m_PsiN;
    }

    if (StateIndex <= 0 || StateIndex < 2) {
        m_PsiN1 = m_PsiN;
        return m_PsiN;
    }

    m_CurrentStateOnset += m_TimeInPrevEvent;

    // Cont (2010), Large and Palmer (1999) and Large and Jones (2002)
    MarkovState &LastState = m_States[StateIndex - 1];
    MarkovState &CurrentState = m_States[StateIndex];
    MarkovState *NextState = nullptr;
    if ((size_t)(StateIndex + 1) < m_States.size()) {
        NextState = &m_States[StateIndex + 1];
    }

    double IOISeconds = m_CurrentStateOnset - m_LastTn;
    double LastPhiN = LastState.IOIPhiN;
    double LastHatPhiN = LastState.IOIHatPhiN;
    double HatPhiN = CurrentState.IOIHatPhiN;
    CurrentState.OnsetObserved = m_CurrentStateOnset;

    if (std::isfinite(CurrentState.SyncStrength)) {
        m_SyncStrength = std::clamp(CurrentState.SyncStrength, 0.0, 1.0);
    }
    if (std::isfinite(CurrentState.PhaseCoupling)) {
        m_PhaseCoupling = std::clamp(CurrentState.PhaseCoupling, 0.0, 2.0);
    }

    // Correction (1): r_n, kappa_n
    double PhaseDiff = (IOISeconds / m_PsiN) - HatPhiN;
    double SyncStrength = m_SyncStr - m_SyncStrength * (m_SyncStr - cos(2.0 * std::numbers::pi * PhaseDiff));
    double Kappa = InverseA2(SyncStrength);
    m_SyncStr = SyncStrength;
    m_Kappa = Kappa;

    // Correction (2): phi_n
    double FValueUpdate = CouplingFunction(LastPhiN, LastHatPhiN, Kappa);
    double PhiN = LastPhiN + (IOISeconds / m_LastPsiN) + (m_PhaseCoupling * FValueUpdate);
    PhiN = ModPhases(PhiN);
    CurrentState.PhaseObserved = PhiN;
    CurrentState.IOIPhiN = PhiN;

    // Prediction: psi_{n+1}
    double FValuePrediction = CouplingFunction(PhiN, HatPhiN, Kappa);
    double PsiN1 = m_PsiN * (1 + m_SyncStrength * FValuePrediction);
    if (PsiN1 < 1e-6) {
        PsiN1 = 1e-6;
    }

    if (NextState != nullptr) {
        double Tn1 = m_CurrentStateOnset + CurrentState.Duration * PsiN1;
        double PhiN1 = ModPhases((Tn1 - m_CurrentStateOnset) / PsiN1);
        NextState->IOIHatPhiN = PhiN1;
        NextState->OnsetExpected = Tn1;

        double LastOnsetExpected = Tn1;
        for (int i = m_CurrentStateIndex + 2; i < m_CurrentStateIndex + 20; i++) {
            if ((size_t)i >= m_States.size()) {
                break;
            }
            MarkovState &FutureState = m_States[i];
            MarkovState &PreviousFutureState = m_States[(i - 1)];
            double Duration = PreviousFutureState.Duration;
            double FutureOnset = LastOnsetExpected + Duration * PsiN1;

            FutureState.OnsetExpected = FutureOnset;
            LastOnsetExpected = FutureOnset;
        }
    }

    m_BPM = 60.0f / PsiN1;
    m_LastPsiN = m_PsiN;
    m_PsiN1 = PsiN1;

    m_TimeInPrevEvent = 0;
    m_LastTn = m_CurrentStateOnset;

    return PsiN1;
}

// ╭─────────────────────────────────────╮
// │     Markov / Semi-Markov Core       │
// ╰─────────────────────────────────────╯
// Section (CONT 2010) section 3.1 and also CUVILLIER (2016) section 2.2.2
void MDP::GetAudioObservations(int T) {
    std::unordered_map<int, double> pitchObs;
    pitchObs.reserve(static_cast<size_t>(std::max(1, m_WinEnd - m_WinStart + 1)));

    const double binWidth = m_Sr / m_FFTSize;

    int bufferIndex = T % m_BufferSize;
    double ObsNoSound = 0;
    double ObsSilence = m_Desc.SilenceProb;
    double nonSilenceWeight = 1.0 - m_Desc.SilenceProb;

    for (int j = m_WinStart; j <= m_WinEnd; j++) {
        if (j < 0 || j >= (int)m_States.size())
            continue;

        MarkovState &StateJ = m_States[j];
        double BestObs = 1e-300;

        switch (StateJ.Type) {
        case NOTE:
        case TRILL: {
            for (AudioState &AS : StateJ.AudioStates) {
                if (AS.Type != PITCH) {
                    spdlog::error("Memory error on creation of Audio States, please report");
                }
                const int rootBin = static_cast<int>(std::round(AS.Freq / binWidth));
                auto it = pitchObs.find(rootBin);
                if (it != pitchObs.end()) {
                    BestObs = std::max(BestObs, it->second);
                } else {
                    double kl = GetPitchSimilarity(AS.Freq) * nonSilenceWeight;
                    pitchObs.emplace(rootBin, kl);
                    BestObs = std::max(BestObs, kl);
                }
            }
            ObsNoSound = std::max(ObsNoSound, BestObs);
            break;
        }
        case CHORD: {
            double ChordKLObs = 1e-300;
            for (AudioState &AS : StateJ.AudioStates) {
                if (AS.Type != PITCH) {
                    spdlog::error("Memory error on creation of Audio States, please report");
                }

                const int rootBin = static_cast<int>(std::round(AS.Freq / binWidth));
                auto it = pitchObs.find(rootBin);
                if (it != pitchObs.end()) {
                    ChordKLObs += it->second;
                } else {
                    double kl = GetPitchSimilarity(AS.Freq) * nonSilenceWeight;
                    pitchObs.emplace(rootBin, kl);
                    ChordKLObs += kl;
                }
            }
            BestObs = ChordKLObs / StateJ.AudioStates.size();
            ObsNoSound = std::max(ObsNoSound, BestObs);
            break;
        }
        case REST: {
            BestObs = std::max(BestObs, m_Desc.SilenceProb);
            break;
        }
        default:
            spdlog::error("Observation not implemented yet");
        }
        StateJ.BestObs[bufferIndex] = BestObs;
    }

    if (ObsNoSound > ObsSilence) {
        spdlog::debug("SOUND   | Sound {:.4f} | Silence {:.4f}", ObsNoSound, ObsSilence);
        m_IsSilence = false;
    } else {
        spdlog::debug("SILENCE | Sound {:.4f} | Silence {:.4f}", ObsNoSound, ObsSilence);
        m_IsSilence = true;
    }
}

// ─────────────────────────────────────
// CONT (2010) section 3.1;
// CUVILLIER (2016) section 2.2.2;
// GONG (2015)
double MDP::GetPitchSimilarity(double Freq) {
    double KLDiv = 0.0;
    double RootBinFreq = round(Freq / (m_Sr / m_FFTSize));
    auto it = m_PitchTemplates.find(RootBinFreq);
    if (it == m_PitchTemplates.end()) {
        spdlog::critical("Pitch template not found for frequency {:.2f} Hz (root bin {:.2f} Hz). This should not "
                         "happen, please report.",
                         Freq, RootBinFreq);
        return 0.0;
    }

    const PitchTemplateArray &PitchTemplate = it->second;
    const auto &reverbSpectralPower = m_Desc.ReverbSpectralPower;
    const auto &normSpectralPower = m_Desc.NormSpectralPower;
    size_t halfFft = static_cast<size_t>(m_FFTSize / 2);

    for (size_t i = 0; i < halfFft; i++) {
        double P = PitchTemplate[i] + reverbSpectralPower[i];
        double Q = normSpectralPower[i];
        if (Q <= 0.0) {
            continue;
        }
        if (P > 0.0) {
            KLDiv += P * std::log(P / Q);
        } else {
            KLDiv += Q;
        }
    }

    double noise_robustness = 1.0 / (1.0 + m_Desc.StdDev);
    KLDiv *= noise_robustness;
    KLDiv = std::exp(-m_PitchScalingFactor * KLDiv);
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
// CUVILLIER and CONT (2014) section 2.1.
double MDP::GetTransProbability(int i, int j) {
    return (i + 1 == j) ? 1.0 : 0.0;
}

// ─────────────────────────────────────
// CUVILLIER (2015)
// Needs review
double MDP::GetOccupancyDistribution(MarkovState &State, int u) {
    double ExpectedFrames = (m_PsiN1 * State.Duration) / m_BlockDur;
    if (ExpectedFrames < 1.0)
        ExpectedFrames = 1.0;

    int key = static_cast<int>(std::round(ExpectedFrames * 10.0));
    BuildDistributionCache(ExpectedFrames); // Builds only if missing

    if (u < static_cast<int>(m_OccupancyCache[key].size())) {
        return m_OccupancyCache[key][u];
    }

    return 0.0;
}

// ─────────────────────────────────────
// CUVILLIER (2015)
// Needs review
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
// CUVILLIER (2015)
// Needs review
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
// GUÉDON (2005) + CUVILLIER (2016)
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
// GUÉDON (2005) + CUVILLIER (2016)
void MDP::SemiMarkov(MarkovState &StateJ, int j, int T, int bufferIndex) {
    double bj = StateJ.BestObs[bufferIndex];

    double f_tilde_j = 0.0;
    double f_tilde_jo = 0.0;
    double ObsProd = 1.0;

    double expectedFrames = (m_PsiN1 * StateJ.Duration) / m_BlockDur;
    if (expectedFrames < 1.0) {
        expectedFrames = 1.0;
    }
    const int distKey = static_cast<int>(std::round(expectedFrames * 10.0));
    BuildDistributionCache(expectedFrames);

    auto occIt = m_OccupancyCache.find(distKey);
    auto survIt = m_SurvivorCache.find(distKey);
    if (occIt == m_OccupancyCache.end() || survIt == m_SurvivorCache.end()) {
        StateJ.Forward[bufferIndex] = bj * 1e-300;
        StateJ.ExitProb[bufferIndex] = bj * 1e-300;
        return;
    }

    const std::vector<double> &occ = occIt->second;
    const std::vector<double> &surv = survIt->second;
    const int maxU = static_cast<int>(occ.size()) - 1;

    for (int u = 1; u <= T + 1; u++) {
        double D_bar_ju = (u < static_cast<int>(surv.size())) ? surv[static_cast<size_t>(u)] : 0.0;
        double d_ju = (u < static_cast<int>(occ.size())) ? occ[static_cast<size_t>(u)] : 0.0;

        if (u == T + 1) {
            f_tilde_j += D_bar_ju * ObsProd * StateJ.InitProb;
            f_tilde_jo += d_ju * ObsProd * StateJ.InitProb;
            break;
        }

        if (u <= maxU) {
            int entryBuf = ((T - u) % m_BufferSize + m_BufferSize) % m_BufferSize;
            double transition_sum = 0.0;

            if (j > 0) {
                double p_ij = GetTransProbability(j - 1, j);
                transition_sum += p_ij * m_States[j - 1].ExitProb[entryBuf];
            }

            f_tilde_j += D_bar_ju * ObsProd * transition_sum;
            f_tilde_jo += d_ju * ObsProd * transition_sum;
        }

        int prevBuf = ((T - u) % m_BufferSize + m_BufferSize) % m_BufferSize;
        double prevObs = StateJ.BestObs[prevBuf];
        double prevNorm = std::max(m_Normalization[prevBuf], 1e-300);
        ObsProd *= prevObs / prevNorm;

        if (ObsProd < 1e-15)
            break;
    }

    StateJ.Forward[bufferIndex] = bj * (f_tilde_j + 1e-300);
    StateJ.ExitProb[bufferIndex] = bj * (f_tilde_jo + 1e-300);
}

// ─────────────────────────────────────
// GUÉDON (2005) + CUVILLIER (2016)
int MDP::Inference(int T) {
    int bIndex = T % m_BufferSize;
    spdlog::debug("WinStart {:04d} | WinFinish {:04d} | BufferSize {:04d} | Tau {:06d} | Kappa {:.4f}", m_WinStart,
                  m_WinEnd, bIndex, m_Tau, m_Kappa);

    // Compute \tilde{f}_j(t) and \tilde{f}_j^o(t)
    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        MarkovState &StateJ = m_States[j];
        if (StateJ.HSMMType == SEMIMARKOV)
            SemiMarkov(StateJ, j, T, bIndex);
        else
            Markov(StateJ, j, T, bIndex);
    }

    // Calculate the Normalization Denominator
    double N = 0.0;
    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        N += m_States[j].Forward[bIndex];
    }

    if (N < 1e-300)
        N = 1e-300;

    m_Normalization[bIndex] = N;

    // Apply Normalization
    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        m_States[j].Forward[bIndex] /= N;
        m_States[j].Forward[bIndex] += 1e-300;
        m_States[j].ExitProb[bIndex] /= N;
        m_States[j].ExitProb[bIndex] += 1e-300;
    }

    // Find the Argmax (Best State)
    double maxVal = 1e-300;
    int bestStateIndex = m_CurrentStateIndex;

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        MarkovState &StateJ = m_States[j];
        double fwd = StateJ.Forward[bIndex];
        if (fwd > maxVal && j >= m_CurrentStateIndex) {
            maxVal = fwd;
            bestStateIndex = j;
        }

        spdlog::debug("State ({}) | Obs = {:.5f}, Forward {:.5f}, Exit Prob {:.5f}", StateJ.Index,
                      StateJ.BestObs[bIndex], StateJ.Forward[bIndex], StateJ.ExitProb[bIndex]);
    }

    if (m_IsSilence && bestStateIndex != m_CurrentStateIndex && m_States[bestStateIndex].Type != REST) {
        return m_CurrentStateIndex;
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

    if (m_Tau == 0) {
        GetInitialDistribution();
    }

    // Run forward inference
    int BestState = Inference(m_Tau);
    m_PsiN = UpdatePsiN(BestState);

    // Advance the score position if a new event was detected
    if (BestState != m_CurrentStateIndex) {
        spdlog::debug("New Event Index {:04d}, Score Position {:04d}", BestState, m_States[BestState].ScorePos);
        m_CurrentStateIndex = BestState;
    }

    return m_States[BestState].ScorePos;
}

} // namespace OpenScofo
