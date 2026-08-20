/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/
#include "OpenScofo.hpp"

#include <algorithm>

#if defined(__APPLE__) || defined(__EMSCRIPTEN__)
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

    * CUVILLIER, P.. On Temporal Coherency of Probabilistic Models for Audio-to-Score Alignment. Theses,
        Université Pierre et Marie Curie - Paris VI, 2016. (2016PA066532).

    * GUÉDON, Y. Hidden Hybrid Markov/Semi-Markov Chains. Computational Statistics & Data
        Analysis, [S.l.], v.49, n.3, p.663–688, 2005.

    * LARGE, E. W.; JONES, M. R. The Dynamics of Attending: How People Track Time-Varying
        Events. Psychological Review, [S.l.], v.106, n.1, p.119–159, 1999.

    * LARGE, E. W.; PALMER, C. Perceiving Temporal Regularity in Music. Cognitive Science,
        [S.l.], v.26, n.1, p.1–37, 2002.

    * BATSCHELET, E. Circular statistics in biology. London: Academic Press, 1981.


*/

// ╭─────────────────────────────────────╮
// │Constructor and Destructor Functions │
// ╰─────────────────────────────────────╯
OnlineForward::OnlineForward() {
    m_SyncStrength = 0.5;
    m_PhaseCoupling = 0.5;
    m_TimeInPrevEvent = 0;
    m_EventWindowSize = 20;
    SetTunning(440);

    InitializeA2Table();

    constexpr int KappaPrecision = 1000;
    m_KappaCache.reserve(static_cast<size_t>(KappaPrecision + 1));
    for (int i = 0; i <= KappaPrecision; i++) {
        double key = i / 1000.0;
        InverseA2(key);
    }
}

// ─────────────────────────────────────
void OnlineForward::UpdateConfiguration(Configuration &Config) {
    m_Sr = Config.SR;
    m_FFTSize = Config.FFTSize;
    m_HopSize = Config.HOPSize;
    m_BlockDur = (1.0 / m_Sr) * m_HopSize;
    m_Tunning = Config.TuningA4;
    m_PitchTemplateSigma = Config.PitchTemplateSigma;
    m_Harmonics = Config.PitchTemplateHarmonics;
    m_SyncStrength = Config.SyncStrength;
    m_PhaseCoupling = Config.PhaseCoupling;
    m_dBTreshold = Config.dBTreshold;
    m_AudioStateChangeReceiver = Config.AudioStateChangeReceiver;

    m_PitchTemplates.clear();
    m_PitchTemplatesPrecomputed.clear();
    m_PitchCQTTemplates.clear();
    m_OccupancyPMFCache.clear();
    m_SurvivorCache.clear();

    m_Normalization.assign(static_cast<size_t>(m_BufferSize + 1), 1.0);
    for (MarkovState &State : m_States) {
        State.Forward.assign(static_cast<size_t>(m_BufferSize + 1), std::numeric_limits<double>::min());
        State.ExitProb.assign(static_cast<size_t>(m_BufferSize + 1), std::numeric_limits<double>::min());
        State.BestObs.assign(static_cast<size_t>(m_BufferSize + 1), std::numeric_limits<double>::min());
    }

    if (!m_States.empty()) {
        UpdateAudioTemplate();
    }
}

// ─────────────────────────────────────
int OnlineForward::GetCurrentStateIndex() {
    return m_CurrentStateIndex;
}

// ─────────────────────────────────────
EventActions OnlineForward::GetCurrentEventActions() {
    MarkovState State = m_States[m_CurrentStateIndex];
    return State.Actions;
}

// ─────────────────────────────────────
EventActions OnlineForward::GetAudioStateChangeActions() {
    EventActions Actions;
    Actions.reserve(m_PendingAudioStateActions.size());

    while (!m_PendingAudioStateActions.empty()) {
        Actions.push_back(std::move(m_PendingAudioStateActions.front()));
        m_PendingAudioStateActions.pop_front();
    }

    return Actions;
}

// ─────────────────────────────────────
void OnlineForward::NotifyAudioStateChange(int StateIndex) {
    if (m_AudioStateChangeReceiver.empty() || StateIndex < 0 || StateIndex >= static_cast<int>(m_States.size())) {
        return;
    }

    const MarkovState &State = m_States[StateIndex];
    const int AudioStateIndex = State.BestAudioStateIndex;

    if (State.Type != CHORD && (AudioStateIndex < 0 || AudioStateIndex >= static_cast<int>(State.AudioStates.size()))) {
        return;
    }

    if (StateIndex == m_LastNotifiedStateIndex && AudioStateIndex == m_LastNotifiedAudioStateIndex) {
        return;
    }

    m_LastNotifiedStateIndex = StateIndex;
    m_LastNotifiedAudioStateIndex = AudioStateIndex;

    ScoreAction Action{};
    Action.isAudioStateChange = true;
    Action.Receiver = m_AudioStateChangeReceiver;
    Action.AbsoluteTime = true;
    Action.Args.emplace_back(State.ScorePos);

    if (State.Type == CHORD) {
        Action.Args.emplace_back(std::string("chord"));
        for (const AudioState &AudioState : State.AudioStates) {
            if (AudioState.Type == PITCH) {
                Action.Args.emplace_back(static_cast<float>(AudioState.Freq));
            }
        }
    } else {
        const AudioState &AudioState = State.AudioStates[static_cast<size_t>(AudioStateIndex)];
        switch (AudioState.Type) {
        case PITCH:
            Action.Args.emplace_back(static_cast<float>(AudioState.Freq));
            break;
        case SILENCE:
            Action.Args.emplace_back(std::string("silence"));
            break;
        case LABEL:
            Action.Args.emplace_back(AudioState.Label);
            break;
        case ONSET:
            Action.Args.emplace_back(std::string("onset"));
            break;
        }
    }

    m_PendingAudioStateActions.push_back(std::move(Action));
}

// ─────────────────────────────────────
void OnlineForward::ResetCaches() {
    // Clear all caches
    m_PitchTemplates.clear();
    m_PitchTemplatesPrecomputed.clear();
    m_PitchProbabilityCache.clear();
    m_OccupancyPMFCache.clear();
    m_SurvivorCache.clear();
    m_KappaCache.clear();

    // Rehash/compact the maps to free memory
    m_PitchTemplates.rehash(0);
    m_PitchTemplatesPrecomputed.rehash(0);
    m_PitchProbabilityCache.rehash(0);
    m_OccupancyPMFCache.rehash(0);
    m_SurvivorCache.rehash(0);
    m_KappaCache.rehash(0);

    spdlog::debug("All caches cleared and rehashed");
}

// ─────────────────────────────────────
void OnlineForward::SetScoreStates(States ScoreStates) {
    if (ScoreStates.size() == 0) {
        return;
    }

    m_States.clear();
    m_States.shrink_to_fit();
    // m_States = ScoreStates;
    m_States = std::move(ScoreStates);
    spdlog::debug("There is {} score states", m_States.size());
    ResetCaches(); // Clear all caches

    spdlog::debug("BufferSize is {}", m_BufferSize);

    m_Normalization.assign(m_BufferSize, 1.0);
    for (MarkovState &State : m_States) {
        State.Forward.assign(m_BufferSize, std::numeric_limits<double>::min());
        State.ExitProb.assign(m_BufferSize, std::numeric_limits<double>::min());
        State.BestObs.assign(m_BufferSize, std::numeric_limits<double>::min());
    }

    m_CurrentStateIndex = 0;
    m_LastNotifiedStateIndex = -1;
    m_LastNotifiedAudioStateIndex = -1;
    m_PendingAudioStateActions.clear();
    m_Kappa = 10;
    m_BPM = m_States[0].BPMExpected;
    m_PsiN = 60.0f / m_States[0].BPMExpected;
    m_PsiN1 = 60.0f / m_States[0].BPMExpected;
    m_LastPsiN = 60.0f / m_States[0].BPMExpected;
    m_BeatsAhead = m_States[0].BPMExpected / 60 * m_SecondsAhead;
    m_SyncStr = 0;
    m_Tau = 0;

    m_SyncStrength = m_States[0].SyncStrength;
    m_PhaseCoupling = m_States[0].PhaseCoupling;

    UpdateAudioTemplate();
}

// ─────────────────────────────────────
// GONG 2015 (adapted)
void OnlineForward::BuildPitchTemplate(double Freq) {
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

    double rootBinFreq = std::round(Freq / binWidth);

    if (m_PitchTemplates.find(rootBinFreq) != m_PitchTemplates.end()) {
        return;
    }

    // Keep a tiny positive floor across the template. Besides avoiding
    // log(0), this is part of the original KL-divergence model: bins outside
    // the harmonic peaks must not be treated as unmatched observed energy.
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

    std::vector<PitchTemplateEntry> precomputedTemplate;
    precomputedTemplate.reserve(templateBins.size());
    for (size_t i = 0; i < templateBins.size(); ++i) {
        double P = templateBins[i];
        if (P > 0.0) {
            precomputedTemplate.push_back({i, P, P * std::log(P)});
        }
    }

    m_PitchTemplatesPrecomputed[rootBinFreq] = std::move(precomputedTemplate);
    m_PitchTemplates[rootBinFreq] = std::move(templateBins);
}

// ─────────────────────────────────────
// GONG 2015 (adapted)
void OnlineForward::UpdateAudioTemplate() {
    int StateSize = (int)m_States.size();
    m_PitchTemplates.clear();
    m_PitchTemplatesPrecomputed.clear();

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
PitchTemplateArray OnlineForward::GetPitchTemplate(double Freq) {
    BuildPitchTemplate(Freq);
    double rootBinFreq = std::round(Freq / (m_Sr / m_FFTSize));
    return m_PitchTemplates[rootBinFreq];
}

// ╭─────────────────────────────────────╮
// │          Set|Get Functions          │
// ╰─────────────────────────────────────╯
void OnlineForward::ClearStates() {
    m_States.clear();
}
// ─────────────────────────────────────
double OnlineForward::GetCurrentBPM() {
    return m_BPM;
}

// ─────────────────────────────────────
double OnlineForward::GetBlockDuration() {
    return m_BlockDur;
}

// ─────────────────────────────────────
void OnlineForward::SetTunning(double Tunning) {
    m_Tunning = Tunning;
}

// ─────────────────────────────────────
void OnlineForward::SetDescription(const Description &Desc) {
    m_Desc = Desc;
    m_PitchProbabilityCache.clear();
    m_PitchProbabilityCacheTau = -1;
    m_ReverbEnergyCacheTau = -1;
}

// ─────────────────────────────────────
void OnlineForward::SetHarmonics(int Harmonics) {
    m_Harmonics = Harmonics;
}

// ─────────────────────────────────────
void OnlineForward::SetAmplitudeDecay(double decay) {
    m_PitchTemplateAmplitudeDecay = decay;
}

// ─────────────────────────────────────
int OnlineForward::GetCurrentBufferIndex() {
    return m_Tau % m_BufferSize;
}

// ─────────────────────────────────────
int OnlineForward::GetTunning() {
    return m_Tunning;
}

// ─────────────────────────────────────
void OnlineForward::SetCurrentEvent(int Event) {
    spdlog::debug("Current event is {}", Event);
    if (m_States.size() == 0) {
        spdlog::error("There is not Events on Score or the Score was no loaded");
        return;
    }

    if (Event == 0) {
        spdlog::info("Initializing Time Decoding Algorithm");
        ResetDecoding(); // Use full reset instead of InitTimeDecoding
    }
    m_CurrentStateIndex = Event;

    // TODO: CHECK THIS
    m_Tau = 0; // Already set in ResetDecoding, but keep for safety
}

// ─────────────────────────────────────
int OnlineForward::GetStatesSize() {
    return m_States.size();
}

// ─────────────────────────────────────
void OnlineForward::AddState(MarkovState State) {
    m_States.push_back(State);
}

// ─────────────────────────────────────
MarkovState OnlineForward::GetState(int Index) {
    return m_States[Index];
}

// ─────────────────────────────────────
std::vector<MarkovState> &OnlineForward::GetStates() {
    return m_States;
}

// ─────────────────────────────────────
void OnlineForward::SetPitchTemplateSigma(double f) {
    m_PitchTemplateSigma = f;
}

// ╭─────────────────────────────────────╮
// │            Time Decoding            │
// ╰─────────────────────────────────────╯
void OnlineForward::InitTimeDecoding(void) {
    double PsiK = 60 / m_States[0].BPMExpected;
    m_LastPsiN = PsiK;
    m_PsiN = PsiK;
    m_PsiN1 = PsiK;
    m_States[0].OnsetObserved = 0;
    m_BPM = m_States[0].BPMExpected;
    m_CurrentStateOnset = 0;
    m_LastTn = 0;
    m_TimeInPrevEvent = 0;
    m_SyncStr = 0;
    m_Kappa = 10;
    m_Tau = 0;
}

// ─────────────────────────────────────
void OnlineForward::ResetDecoding() {
    // Reset timing variables
    m_Tau = 0;
    m_TimeInPrevEvent = 0;
    m_CurrentStateOnset = 0;
    m_LastTn = 0;

    // Reset synchronization variables
    m_SyncStr = 0;
    m_Kappa = 10;

    // Reset BPM and period predictions
    double PsiK = 60 / m_States[0].BPMExpected;
    m_LastPsiN = PsiK;
    m_PsiN = PsiK;
    m_PsiN1 = PsiK;
    m_BPM = m_States[0].BPMExpected;

    // Reset all state probabilities
    for (MarkovState &State : m_States) {
        std::fill(State.Forward.begin(), State.Forward.end(), 0.0);
        std::fill(State.ExitProb.begin(), State.ExitProb.end(), 0.0);
        std::fill(State.BestObs.begin(), State.BestObs.end(), std::numeric_limits<double>::min());
        State.BestAudioStateIndex = -1;
        State.OnsetObserved = 0;
        State.PhaseObserved = 0;
        State.IOIPhiN = 0;
        State.InitProb = 0.0;
    }

    // Reset normalization buffer
    std::fill(m_Normalization.begin(), m_Normalization.end(), 1.0);

    // Reset first state
    m_States[0].OnsetObserved = 0;

    m_LastNotifiedStateIndex = -1;
    m_LastNotifiedAudioStateIndex = -1;
    m_PendingAudioStateActions.clear();

    spdlog::debug("OnlineForward decoding state fully reset");
}

// ─────────────────────────────────────
// CUVILLIER 2016 + robust right-tail extension.
//
// Main occupancy law:
//     D_l ~ Poisson(l)
//
// A small exponential right tail is added after the expected duration.
// This keeps normal rhythmic passages close to Cuvillier's Poisson model,
// while allowing occasional unexpectedly long events (breaths, phrase
// endings, expressive lengthening).
void OnlineForward::BuildDistributionCache(double ExpectedFrames) {
    ExpectedFrames = std::max(ExpectedFrames, 1.0);
    const int key = static_cast<int>(ExpectedFrames * 10.0 + 0.5);
    if (m_OccupancyPMFCache.contains(key) && m_SurvivorCache.contains(key)) {
        return;
    }
    const double lambda = static_cast<double>(key) * 0.1;

    // TODO: This is a score parameter
    constexpr double timeTolerance = 0.03;

    const double tailScale = std::max(2.0, lambda);
    const double sigma = std::sqrt(lambda);
    const double poissonLimit = lambda + 8.0 * sigma;
    const double tailLimit = lambda + 8.0 * tailScale;

    const int requestedMaxU = static_cast<int>(std::ceil(std::max({5.0 * lambda, poissonLimit, tailLimit})));
    const int maxU = std::clamp(requestedMaxU, 1, m_BufferSize - 1);
    std::vector<double> pmf(static_cast<size_t>(maxU + 1), 0.0);
    std::vector<double> survivor(static_cast<size_t>(maxU + 2), 0.0);

    // Poisson PMF.
    double poissonProbability = std::exp(-lambda);
    double totalProbability = 0.0;

    for (int u = 1; u <= maxU; ++u) {
        poissonProbability *= lambda / static_cast<double>(u);
        double probability = poissonProbability;

        // Heavy RIGHT tail.
        if (static_cast<double>(u) > lambda && timeTolerance > 0.0) {
            const double excess = static_cast<double>(u) - lambda;
            const double tail = std::exp(-excess / tailScale);
            probability = (1.0 - timeTolerance) * probability + timeTolerance * tail;
        }
        pmf[static_cast<size_t>(u)] = probability;
        totalProbability += probability;
    }

    // Normalize.
    if (totalProbability > std::numeric_limits<double>::min()) {
        const double invTotal = 1.0 / totalProbability;
        for (int u = 1; u <= maxU; ++u) {
            pmf[static_cast<size_t>(u)] *= invTotal;
        }
    } else {
        pmf[static_cast<size_t>(maxU)] = 1.0;
    }

    // Survivor:
    // D(u) = P(U >= u)
    double cumulative = 0.0;
    for (int u = maxU; u >= 1; --u) {
        cumulative += pmf[static_cast<size_t>(u)];
        survivor[static_cast<size_t>(u)] = cumulative;
    }
    survivor[0] = 1.0;
    m_OccupancyPMFCache.emplace(key, std::move(pmf));
    m_SurvivorCache.emplace(key, std::move(survivor));
}

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
double OnlineForward::CalculateA2(double kappa) {
    if (kappa <= 0.0) {
        return 0.0;
    }

    if (kappa > 10.0) {
        return 1.0 - (1.0 / (2.0 * kappa)) - (1.0 / (8.0 * kappa * kappa));
    }

    double I1 = CYL_BESSEL_I(1, kappa);
    double I0 = CYL_BESSEL_I(0, kappa);
    return I1 / I0;
}

// ─────────────────────────────────────
// CONT (2010) and BATSCHELET (1981)
void OnlineForward::InitializeA2Table() {
    if (m_A2TableInitialized) {
        return;
    }

    for (int i = 0; i < A2TableSize; ++i) {
        double kappa = static_cast<double>(i) / A2TablePrecision;
        m_A2Table[i] = CalculateA2(kappa);
    }
    m_A2TableInitialized = true;
}

// ─────────────────────────────────────
// CONT 2010
double OnlineForward::A2(double kappa) {
    InitializeA2Table();

    if (kappa <= 0.0) {
        return m_A2Table.front();
    }

    if (kappa >= A2TableMaxKappa) {
        return m_A2Table.back();
    }

    int index = static_cast<int>(std::round(kappa * A2TablePrecision));
    index = std::clamp(index, 0, A2TableSize - 1);
    return m_A2Table[static_cast<size_t>(index)];
}

// ─────────────────────────────────────
// CONT 2010
double OnlineForward::InverseA2(double SyncStrength) {
    InitializeA2Table();

    SyncStrength = std::clamp(SyncStrength, 0.0, 1.0);
    constexpr int KappaPrecision = 1000;
    int key = static_cast<int>(SyncStrength * KappaPrecision);
    auto it = m_KappaCache.find(key);
    if (it != m_KappaCache.end()) {
        return it->second;
    }

    // Following Large and Jones (1999, p. 157).
    if (SyncStrength > 0.95) {
        m_KappaCache[key] = 10.0;
        return 10.0;
    }

    auto lower = std::lower_bound(m_A2Table.begin(), m_A2Table.end(), SyncStrength);
    if (lower == m_A2Table.begin()) {
        m_KappaCache[key] = 0.0;
        return 0.0;
    }

    if (lower == m_A2Table.end()) {
        m_KappaCache[key] = 10.0;
        return 10.0;
    }

    size_t upperIndex = static_cast<size_t>(std::distance(m_A2Table.begin(), lower));
    size_t lowerIndex = upperIndex - 1;
    double lowerA2 = m_A2Table[lowerIndex];
    double upperA2 = m_A2Table[upperIndex];
    size_t nearestIndex = (SyncStrength - lowerA2 <= upperA2 - SyncStrength) ? lowerIndex : upperIndex;
    double Kappa = static_cast<double>(nearestIndex) / A2TablePrecision;
    m_KappaCache[key] = Kappa;
    return Kappa;
}

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
double OnlineForward::CouplingFunction(double phi, double phi_hat, double kappa) {
    static constexpr double invTwoPi = 1.0 / (2.0 * std::numbers::pi);
    double diff = 2.0 * std::numbers::pi * (phi - phi_hat);
    double cosDiff = std::cos(diff);
    return invTwoPi * std::exp(kappa * (cosDiff - 1.0)) * std::sin(diff);
}

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
double OnlineForward::ModPhases(double Phase) {
    Phase = std::fmod(Phase + 0.5, 1.0);
    if (Phase < 0.0) {
        Phase += 1.0;
    }
    return Phase - 0.5;
}

// ─────────────────────────────────────
// CONT 2010 (Last § of section 4)
States OnlineForward::GetStatesForProcessing() {
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
void OnlineForward::GetDecodeWindow() {
    int half = m_EventWindowSize / 2;
    m_WinStart = std::max(0, static_cast<int>(m_CurrentStateIndex) - half);
    m_WinEnd = std::min(static_cast<int>(m_States.size()) - 1, static_cast<int>(m_CurrentStateIndex) + half);

    if (m_WinStart < 0 || m_WinEnd >= static_cast<int>(m_States.size())) {
        spdlog::critical("Inference::GetDecodeWindow invariant violated: "
                         "window out of bounds "
                         "(winStart={}, winEnd={}, statesSize={}, currentIndex={}, eventWindowSize={})",
                         m_WinStart, m_WinEnd, m_States.size(), m_CurrentStateIndex, m_EventWindowSize);
    }
}

// ─────────────────────────────────────
// CONT 2010 (Section 5, algorithm 1)
double OnlineForward::UpdatePsiN(int StateIndex) {
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
    if (StateIndex + 1 > (int)m_States.size()) {
        return m_PsiN;
    }

    // Cont (2010), Large and Palmer (1999) and Large and Jones (2002)
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

    // Update Variance (Cont, 2010) - Coupling Strength (Large 1999)
    double PhaseDiff = (IOISeconds / m_PsiN) - HatPhiN;
    double SyncStrength = m_SyncStr - m_SyncStrength * (m_SyncStr - cos(std::numbers::pi * 2 * PhaseDiff));
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
        MarkovState &FutureState = m_States[i];
        MarkovState &PreviousFutureState = m_States[(i - 1)];
        double Duration = PreviousFutureState.Duration;
        double FutureOnset = LastOnsetExpected + Duration * PsiN1;

        FutureState.OnsetExpected = FutureOnset;
        LastOnsetExpected = FutureOnset;
    }

    // Update Values for next calls
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
void OnlineForward::GetAudioObservations() {
    double soundProb = std::max(0.0, 1.0 - m_Desc.SilenceProb);
    double techWeight = m_Desc.ExtendedTechProb;
    double pitchWeight = 1.0 - m_Desc.ExtendedTechProb;

    EventType CurrentEventType = m_States[m_CurrentStateIndex].Type;
    double maxSoundEvidence = 0.0;
    bool allowSilence = (CurrentEventType != FIRSTEVENT) && (CurrentEventType != REST);

    // Global best AudioState across all semi-Markov states
    double globalBestAudioStateEvidence = 0.0;
    int globalBestStateIndex = -1;
    int globalBestAudioStateIndex = -1;

    for (int j = m_WinStart; j <= m_WinEnd; j++) {
        MarkovState &state = m_States[j];

        double bestPitch = 0.0;
        double bestTech = 0.0;
        double bestOnset = 0.0;
        double bestSilence = 0.0;

        double sumPitch = 0.0;
        int pitchCount = 0;

        for (size_t audioStateIndex = 0; audioStateIndex < state.AudioStates.size(); ++audioStateIndex) {

            const AudioState &as = state.AudioStates[audioStateIndex];
            double audioStateEvidence = 0.0;

            switch (as.Type) {
            case PITCH: {
                double p = GetPitchProbability(as.Freq);
                audioStateEvidence = p;

                bestPitch = std::max(bestPitch, p);
                sumPitch += p;
                pitchCount++;
                break;
            }

            case LABEL: {
                double p = m_Desc.ONNX[as.Label];
                audioStateEvidence = p;

                bestTech = std::max(bestTech, p);
                break;
            }

            case ONSET: {
                audioStateEvidence = m_Desc.Onset;
                bestOnset = std::max(bestOnset, m_Desc.Onset);
                break;
            }

            case SILENCE: {
                if (allowSilence) {
                    audioStateEvidence = m_Desc.SilenceProb;
                    bestSilence = std::max(bestSilence, m_Desc.SilenceProb);
                }
                break;
            }
            }

            // Compare against every AudioState of every MarkovState
            if (audioStateEvidence > globalBestAudioStateEvidence) {
                globalBestAudioStateEvidence = audioStateEvidence;
                globalBestStateIndex = j;
                globalBestAudioStateIndex = static_cast<int>(audioStateIndex);
            }
        }

        double stateLikelihood = 0.0;

        switch (state.Type) {
        case NOTE:
        case TRILL:
            stateLikelihood = std::max(bestPitch * pitchWeight * soundProb, bestSilence);
            break;

        case PTECH: {
            double techObs = bestTech * techWeight * soundProb;
            double pitchObs = bestPitch * pitchWeight * soundProb;

            stateLikelihood = std::max({techObs, pitchObs, bestSilence});
            break;
        }

        case UTECH:
            stateLikelihood = std::max(bestTech * techWeight * soundProb, bestSilence);
            break;

        case CHORD:
            if (pitchCount > 0)
                stateLikelihood = (sumPitch / pitchCount) * pitchWeight * soundProb;
            break;

        case FIRSTEVENT:
        case REST:
            stateLikelihood = m_Desc.SilenceProb;
            break;

        default:
            spdlog::error("Event type of line {} of score file is not implemented, please remove it", state.Line);
            break;
        }

        state.BestObs[m_CircularBufferIndex] = std::max(stateLikelihood, std::numeric_limits<double>::min());

        if (state.Type != REST)
            maxSoundEvidence = std::max(maxSoundEvidence, stateLikelihood);
    }

    m_IsSilence = (m_Desc.SilenceProb > maxSoundEvidence);

    // Save global winner
    m_BestAudioStateStateIndex = globalBestStateIndex;
    m_BestAudioStateIndex = globalBestAudioStateIndex;
}

// ─────────────────────────────────────
// CONT (2010) section 3.1;
// CUVILLIER (2016) section 2.2.2;
// GONG (2015)
double OnlineForward::GetPitchProbability(double Freq) {
    constexpr double minProb = std::numeric_limits<double>::min();

    if (Freq <= 0.0 || m_FFTSize <= 0.0 || m_Sr <= 0.0) {
        return minProb;
    }

    const double binWidth = m_Sr / m_FFTSize;
    const double rootBinFreq = std::round(Freq / binWidth);

    auto it = m_PitchTemplates.find(rootBinFreq);
    if (it == m_PitchTemplates.end()) {
        BuildPitchTemplate(Freq);
        it = m_PitchTemplates.find(rootBinFreq);

        if (it == m_PitchTemplates.end()) {
            return minProb;
        }
    }

    const PitchTemplateArray &pitchTemplate = it->second;
    const auto &reverb = m_Desc.ReverbSpectralPower;
    const auto &spectrum = m_Desc.SpectralMagnitudeFrameNorm;

    const size_t bins = std::min({static_cast<size_t>(m_FFTSize / 2), pitchTemplate.size(), spectrum.size()});

    if (bins == 0) {
        return minProb;
    }

    const double *templateData = pitchTemplate.data();
    const double *spectrumData = spectrum.data();
    const double *reverbData = reverb.data();

    const size_t reverbBins = std::min(bins, reverb.size());

    double klDiv = 0.0;

    // Bins where reverb data exists.
    size_t i = 0;
    for (; i < reverbBins; ++i) {
        const double P = templateData[i] + reverbData[i];
        const double Q = spectrumData[i];

        if (P > 0.0 && Q > 0.0) {
            klDiv += P * std::log(P / Q);
        } else if (P == 0.0 && Q >= 0.0) {
            klDiv += Q;
        }
    }

    // Remaining bins: reverb = 0.
    for (; i < bins; ++i) {
        const double P = templateData[i];
        const double Q = spectrumData[i];

        if (P > 0.0 && Q > 0.0) {
            klDiv += P * std::log(P / Q);
        } else if (P == 0.0 && Q >= 0.0) {
            klDiv += Q;
        }
    }

    klDiv /= 1.0 + m_Desc.StdDev;
    return std::exp(-m_PitchScalingFactor * klDiv);
}

// ─────────────────────────────────────
void OnlineForward::GetInitialDistribution() {
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

    std::vector<double> NonNullInitialProb(Size, 0.0);
    for (int destination = 0; destination < Size; ++destination) {
        const MarkovState &DestinationState = m_States[m_CurrentStateIndex + destination];
        const double destinationExpectedFrames = std::max(1.0, (m_PsiN1 * DestinationState.Duration) / m_BlockDur);
        const double destinationLambda = std::round(destinationExpectedFrames * 10.0) / 10.0;
        const double destinationNonNull =
            DestinationState.HSMMType == SEMIMARKOV ? -std::expm1(-destinationLambda) : 1.0;
        double probability = 0.0;
        for (int source = 0; source <= destination; ++source) {
            double skipped = 1.0;
            for (int k = source; k < destination; ++k) {
                const MarkovState &SkippedState = m_States[m_CurrentStateIndex + k];
                if (SkippedState.HSMMType == SEMIMARKOV) {
                    const double expectedFrames = std::max(1.0, (m_PsiN1 * SkippedState.Duration) / m_BlockDur);
                    const double lambda = std::round(expectedFrames * 10.0) / 10.0;
                    skipped *= std::exp(-lambda);
                } else {
                    skipped = 0.0;
                    break;
                }
            }
            probability += InitialProb[source] * skipped;
        }
        NonNullInitialProb[destination] = destinationNonNull * probability;
    }
    const double nonNullSum = std::accumulate(NonNullInitialProb.begin(), NonNullInitialProb.end(), 0.0);
    if (nonNullSum > std::numeric_limits<double>::min())
        for (double &probability : NonNullInitialProb)
            probability /= nonNullSum;

    for (int j = m_WinStart; j <= m_WinEnd; j++) {
        if (j < 0 || j >= (int)m_States.size())
            continue;
        int idx = j - m_CurrentStateIndex;
        MarkovState &StateJ = m_States[j];
        StateJ.InitProb =
            (idx >= 0 && idx < static_cast<int>(NonNullInitialProb.size())) ? NonNullInitialProb[idx] : 0.0;
        if (j < m_CurrentStateIndex) {
            m_States[j].InitProb = 0.0;
        }
    }
    return;
}

// ─────────────────────────────────────
// CUVILLIER and CONT (2014) section 2.1.
double OnlineForward::GetTransProbability(int i, int j) {
    return (i + 1 == j) ? 1.0 : 0.0;
}

// ─────────────────────────────────────
// CUVILLIER (2015)
// TODO: Needs review
double OnlineForward::GetOccupancyDistribution(MarkovState &State, int u) {
    double ExpectedFrames = (m_PsiN1 * State.Duration) / m_BlockDur;
    if (ExpectedFrames < 1.0) {
        ExpectedFrames = 1.0;
    }

    const int key = static_cast<int>(ExpectedFrames * 10.0 + 0.5);
    BuildDistributionCache(ExpectedFrames);

    const auto &cache = m_OccupancyPMFCache[key];
    if (u >= 0 && u < static_cast<int>(cache.size())) {
        return cache[u];
    }

    return 0.0;
}
// ─────────────────────────────────────
// CUVILLIER (2015)
// TODO: Needs review
double OnlineForward::GetSurvivorDistribution(MarkovState &State, int u) {
    double ExpectedFrames = (m_PsiN1 * State.Duration) / m_BlockDur;
    if (ExpectedFrames < 1.0) {
        ExpectedFrames = 1.0;
    }

    const int key = static_cast<int>(ExpectedFrames * 10.0 + 0.5);
    BuildDistributionCache(ExpectedFrames);

    const auto &cache = m_SurvivorCache[key];
    // survivor[u] defined for u = 0..maxU+1, returns 0 beyond
    if (u >= 0 && u < static_cast<int>(cache.size())) {
        return cache[u];
    }
    return 0.0;
}

// ─────────────────────────────────────
// CUVILLIER (2015)
// TODO: Needs review
int OnlineForward::GetMaxUForJ(MarkovState &StateJ) {
    double Expected_Frames = (m_PsiN1 * StateJ.Duration) / m_BlockDur;
    if (Expected_Frames < 1.0) {
        Expected_Frames = 1.0;
    }
    int Cap = static_cast<int>(std::ceil(5.0 * Expected_Frames));
    return std::max(Cap, 1);
}

// ─────────────────────────────────────
// GUÉDON (2005) + CUVILLIER (2016)
void OnlineForward::Markov(MarkovState &StateJ, int j) {
    double Bj = StateJ.BestObs[m_CircularBufferIndex];
    double Fj;

    if (m_Tau == 0) {
        Fj = Bj * StateJ.InitProb;
    } else {
        int prevBuf = (m_CircularBufferIndex - 1 + m_BufferSize) % m_BufferSize;
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
        Fj = Bj * sumPrev;
    }

    // For Markov states Forward = ExitProb (they can exit at every step)
    StateJ.Forward[m_CircularBufferIndex] = Fj;
    StateJ.ExitProb[m_CircularBufferIndex] = Fj;
}

// ─────────────────────────────────────
// GUÉDON (2005) + CUVILLIER (2016)
void OnlineForward::SemiMarkov(MarkovState &StateJ, int j) {
    double Bj = StateJ.BestObs[m_CircularBufferIndex];

    double FTildeJ = 0.0;
    double FTildeJo = 0.0;
    double ObsProd = 1.0;

    double ExpectedFrames = std::max(1.0, (m_PsiN1 * StateJ.Duration) / m_BlockDur);
    const int key = static_cast<int>(ExpectedFrames * 10.0 + 0.5);
    BuildDistributionCache(ExpectedFrames);
    const auto &surv_cache = m_SurvivorCache[key];
    const auto &occ_cache = m_OccupancyPMFCache[key];
    const int maxU = static_cast<int>(occ_cache.size()) - 1;
    const int observedHistory = std::min(m_Tau, maxU);

    for (int u = 1; u <= observedHistory; ++u) {
        const double Dju = surv_cache[u];
        const double dju = occ_cache[u];
        const int EntryBuf = ((m_Tau - u) % m_BufferSize + m_BufferSize) % m_BufferSize;
        double TransSum = 0.0;
        for (int i = m_WinStart; i < j; ++i)
            TransSum += GetTransProbability(i, j) * m_States[i].ExitProb[EntryBuf];
        FTildeJ += Dju * ObsProd * TransSum;
        FTildeJo += dju * ObsProd * TransSum;

        const double prevObs = StateJ.BestObs[EntryBuf];
        const double prevNorm = m_Normalization[EntryBuf];
        if (prevNorm > std::numeric_limits<double>::min())
            ObsProd *= prevObs / prevNorm;
        else
            ObsProd = 0.0;
        if (ObsProd == 0.0)
            break;
    }

    const int initialDuration = m_Tau + 1;
    if (initialDuration <= maxU) {
        FTildeJ += surv_cache[initialDuration] * ObsProd * StateJ.InitProb;
        FTildeJo += occ_cache[initialDuration] * ObsProd * StateJ.InitProb;
    }

    StateJ.Forward[m_CircularBufferIndex] = Bj * (FTildeJ + std::numeric_limits<double>::min());
    StateJ.ExitProb[m_CircularBufferIndex] = Bj * (FTildeJo + std::numeric_limits<double>::min());
}

// ─────────────────────────────────────
// GUÉDON (2005) + CUVILLIER (2016)
int OnlineForward::GetAlphaT() {
    spdlog::debug("WinStart {:04d} | WinFinish {:04d} | BufferSize {:04d} | Tau {:06d} | Kappa {:.4f}", m_WinStart,
                  m_WinEnd, m_CircularBufferIndex, m_Tau, m_Kappa);

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        MarkovState &StateJ = m_States[j];
        switch (StateJ.HSMMType) {
        case SEMIMARKOV:
            SemiMarkov(StateJ, j);
            break;
        case MARKOV:
            Markov(StateJ, j);
            break;
        }
    }

    // Calculate the Normalization Denominator
    double N = 0.0;
    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        N += m_States[j].Forward[m_CircularBufferIndex];
    }

    if (N < std::numeric_limits<double>::min()) {
        N = std::numeric_limits<double>::min();
    }

    m_Normalization[m_CircularBufferIndex] = N;

    // Apply Normalization
    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        m_States[j].Forward[m_CircularBufferIndex] /= N;
        m_States[j].Forward[m_CircularBufferIndex] += std::numeric_limits<double>::min();
        m_States[j].ExitProb[m_CircularBufferIndex] /= N;
        m_States[j].ExitProb[m_CircularBufferIndex] += std::numeric_limits<double>::min();
    }

    // Find the Argmax (Best State)
    double maxVal = std::numeric_limits<double>::min();
    int BestStateIndex = m_CurrentStateIndex;

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        MarkovState &StateJ = m_States[j];
        double fwd = StateJ.Forward[m_CircularBufferIndex];
        if (fwd > maxVal && j >= m_CurrentStateIndex) {
            maxVal = fwd;
            BestStateIndex = j;
        }

        spdlog::debug("State ({}) | Obs = {:.5f}, Forward {:.5f}, Exit Prob {:.5f}", StateJ.Index,
                      StateJ.BestObs[m_CircularBufferIndex], StateJ.Forward[m_CircularBufferIndex],
                      StateJ.ExitProb[m_CircularBufferIndex]);
    }

    if (m_IsSilence && BestStateIndex != m_CurrentStateIndex && m_States[BestStateIndex].Type != REST) {
        return m_CurrentStateIndex;
    }

    MarkovState &BestState = m_States[BestStateIndex];
    spdlog::debug("Best: State ({}) | Forward {:.5f}", BestState.Index, BestState.Forward[m_CircularBufferIndex]);

    return BestStateIndex;
}

// ─────────────────────────────────────
int OnlineForward::GetEvent(Description &Desc) {
    spdlog::debug("Starting inference");
    m_Desc = Desc;

    m_CircularBufferIndex = m_Tau % m_BufferSize;

    if (m_CurrentStateIndex > (int)m_States.size()) {
        spdlog::debug("Score Finished");
        return m_States.back().ScorePos;
    }

    GetDecodeWindow();
    GetAudioObservations();

    if (m_Tau == 0) {
        GetInitialDistribution();
    }

    // Run forward inference
    int BestState = GetAlphaT();
    m_PsiN = UpdatePsiN(BestState);

    // Advance the score position if a new event was detected
    if (BestState != m_CurrentStateIndex) {
        spdlog::debug("New Event Index {:04d}, Score Position {:04d}", BestState, m_States[BestState].ScorePos);
        m_CurrentStateIndex = BestState;
    }

    NotifyAudioStateChange(BestState);
    return m_States[BestState].ScorePos;
}

} // namespace OpenScofo
