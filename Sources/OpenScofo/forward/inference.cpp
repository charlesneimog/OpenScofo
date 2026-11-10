/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file inference.cpp
 * @brief Online hybrid Markov/semi-Markov score inference and decoding.
 */
#include "OpenScofo.hpp"
#include <algorithm>

namespace OpenScofo {

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
/**
 * @brief Cache occupancy and survivor laws for an expected duration.
 *
 * @param ExpectedFrames Expected event duration in analysis frames.
 *
 * @note Quantizes durations to tenths of a frame and adds an exponential right tail to the Poisson occupancy law.
 * @warning The cached distribution is truncated to the available circular history length.
 */
inline void OnlineForward::BuildDistributionCache(double ExpectedFrames) {
    ExpectedFrames = std::max(ExpectedFrames, 1.0);
    const int key = static_cast<int>(ExpectedFrames * 10.0 + 0.5);
    if (m_OccupancyPMFCache.contains(key) && m_SurvivorCache.contains(key)) {
        return;
    }
    const double lambda = static_cast<double>(key) * 0.1;

    // TODO: Should this be a score parameter?
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
// CONT 2010 (Last § of section 4)
/**
 * @brief Copy upcoming states within the look-ahead horizon.
 *
 * @return Copy of the selected state range.
 *
 * @note Starts at the current state and accumulates durations until the configured horizon is exceeded.
 * @warning Requires a valid current index in a nonempty score.
 */
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
/**
 * @brief Determine the inclusive internal state range for inference.
 *
 * @note Uses all contiguous states of the current section when restricted, otherwise a centered event window.
 * @warning Requires a valid current index in a nonempty score.
 */
void OnlineForward::GetDecodeWindow() {
    if (m_SectionRestrict) {
        const std::string &Section = m_States[static_cast<size_t>(m_CurrentStateIndex)].Section;
        m_WinStart = m_CurrentStateIndex;
        while (m_WinStart > 0 && m_States[static_cast<size_t>(m_WinStart - 1)].Section == Section) {
            m_WinStart--;
        }

        m_WinEnd = m_CurrentStateIndex;
        while (m_WinEnd + 1 < static_cast<int>(m_States.size()) &&
               m_States[static_cast<size_t>(m_WinEnd + 1)].Section == Section) {

            m_WinEnd++;
        }
        return;
    }

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
/**
 * @brief Initialize the score-position prior in the active window.
 *
 * @note This is a OpenScofo implementation
 */
void OnlineForward::GetInitialDistribution() {
    double duration = 0.0;
    double sum = 0.0;

    // Weight states by their distance from the current position.
    for (int j = m_CurrentStateIndex; j <= m_WinEnd; ++j) {
        ScoreState &state = m_States[j];
        state.InitProb = state.IsInterEventSilence ? 0.0 : std::exp(-duration / m_BeatsAhead);
        duration += state.Duration;
        sum += state.InitProb;
    }

    // Normalize the active prior.
    if (sum > 0.0) {
        for (int j = m_CurrentStateIndex; j <= m_WinEnd; ++j) {
            m_States[j].InitProb /= sum;
        }
    }

    // Clear states behind the current position.
    for (int j = m_WinStart; j < m_CurrentStateIndex; ++j) {
        m_States[j].InitProb = 0.0;
    }
}

// ─────────────────────────────────────
/**
 * @brief Use a linear kind of left-to-right score topology.
 *
 * @note Initialize the score-position prior in the active window. This is just difference because I believe that should
 * be 50% have silence (on non-legato)
 */
double OnlineForward::GetSemiMarkovTransitionProbability(int i, int j) {
    if (i < 0 || j <= i || j >= static_cast<int>(m_States.size())) {
        return 0.0;
    }

    // Inter-event silence is optional: the decoder may either enter the
    // silence state or skip directly to the next musical event.
    //
    // Ideally this prior could be controlled by score information, e.g.
    // articulation/legato. Since OpenScofo currently has no such parameter,
    // assign equal prior probability (0.5) to both paths. The observation
    // likelihoods will subsequently determine which path is favored.
    // This is something that is implemented different from Antescofo
    if (i + 1 < static_cast<int>(m_States.size()) && m_States[i + 1].IsInterEventSilence) {
        return (j == i + 1 || j == i + 2) ? 0.5 : 0.0;
    }
    return j == i + 1 ? 1.0 : 0.0;
}

// ─────────────────────────────────────
// CUVILLIER (2015)
// TODO: Needs review
/**
 * @brief Read the probability of an event lasting u analysis frames.
 *
 * @param State Score state to read or update.
 * @param u Event duration or age in analysis frames.
 *
 * @return Occupancy probability, or zero outside the cached duration range.
 *
 * @note Builds the quantized duration cache on demand from the predicted beat period.
 */
double OnlineForward::GetOccupancyDistribution(ScoreState &State, int u) {
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
/**
 * @brief Read the probability of an event surviving at least u frames.
 *
 * @param State Score state to read or update.
 * @param u Event duration or age in analysis frames.
 *
 * @return Survivor probability, or zero outside the cached range.
 *
 * @note Uses the survivor law corresponding to the cached occupancy distribution.
 */
double OnlineForward::GetSurvivorDistribution(ScoreState &State, int u) {
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
/**
 * @brief Estimate a duration-history cap for a score state.
 *
 * @param StateJ Score state whose forward probabilities or duration are evaluated.
 *
 * @return Estimated maximum duration in analysis frames.
 *
 * @note Returns five times the expected duration, rounded up, with a minimum of one frame.
 */
int OnlineForward::GetMaxUForJ(ScoreState &StateJ) {
    double Expected_Frames = (m_PsiN1 * StateJ.Duration) / m_BlockDur;
    if (Expected_Frames < 1.0) {
        Expected_Frames = 1.0;
    }
    int Cap = static_cast<int>(std::ceil(5.0 * Expected_Frames));
    return std::max(Cap, 1);
}

// ─────────────────────────────────────
// GUÉDON (2005) + CUVILLIER (2016)
/**
 * @brief Update score-level Markov forward and exit probabilities.
 *
 * @param StateJ Score state whose forward probabilities or duration are evaluated.
 * @param j Destination state index.
 *
 * @note Uses current evidence and predecessor exits; inserted silence states have equal stay and leave priors.
 * @warning j must identify StateJ and all circular probability buffers must be initialized.
 */
void OnlineForward::Markov(ScoreState &StateJ, int j) {
    double Bj = StateJ.BestObs[m_CircularBufferIndex];
    double Fj;

    if (m_Tau == 0) {
        Fj = Bj * StateJ.InitProb;
    } else {
        int prevBuf = (m_CircularBufferIndex - 1 + m_BufferSize) % m_BufferSize;
        double sumPrev = 0.0;

        // Stay in j (self-loop)
        if (j >= m_CurrentStateIndex) {
            sumPrev += (StateJ.IsInterEventSilence ? 0.5 : 1.0) * StateJ.Forward[prevBuf];
        }

        // Arrive from j-1
        if (j - 1 >= m_CurrentStateIndex && j - 1 >= 0) {
            double trans = GetSemiMarkovTransitionProbability(j - 1, j);
            sumPrev += trans * m_States[j - 1].ExitProb[prevBuf];
        }
        Fj = Bj * sumPrev;
    }

    // Gap duration is geometric: stay or leave with equal probability. Exit
    // mass is consumed once by the next state's transition, not counted twice.
    StateJ.Forward[m_CircularBufferIndex] = Fj;
    StateJ.ExitProb[m_CircularBufferIndex] = (StateJ.IsInterEventSilence ? 0.5 : 1.0) * Fj;
}

// ─────────────────────────────────────
// GUÉDON (2005) + CUVILLIER (2016)
/**
 * @brief Update duration-conditioned forward and exit probabilities.
 *
 * @param StateJ Score state whose forward probabilities or duration are evaluated.
 * @param j Destination state index.
 *
 * @note Uses prepared segment likelihoods for both predecessor entries and the initial segment.
 * @warning j must identify StateJ and all circular probability buffers must be initialized.
 */
void OnlineForward::SemiMarkov(ScoreState &StateJ, int j) {
    const double ExpectedFrames = std::max(1.0, (m_PsiN1 * StateJ.Duration) / m_BlockDur);

    const int key = static_cast<int>(ExpectedFrames * 10.0 + 0.5);
    BuildDistributionCache(ExpectedFrames);

    const auto &surv_cache = m_SurvivorCache[key];
    const auto &occ_cache = m_OccupancyPMFCache[key];

    const int maxU = static_cast<int>(occ_cache.size()) - 1;
    const int observedHistory = std::min(m_Tau, maxU);

    PrepareSegmentLikelihoods(StateJ, std::min(m_Tau + 1, maxU));

    // Precompute the only two possible predecessor transitions.
    const int Prev1 = j - 1;
    const int Prev2 = j - 2;

    const double P1 = Prev1 >= m_WinStart ? GetSemiMarkovTransitionProbability(Prev1, j) : 0.0;
    const double P2 = Prev2 >= m_WinStart ? GetSemiMarkovTransitionProbability(Prev2, j) : 0.0;

    double FTildeJ = 0.0;
    double FTildeJo = 0.0;
    const bool IsLeftRight = StateJ.MicroTopologyType == LEFT_RIGHT;

    // This is a optimization, once that we
    for (int u = 1; u <= observedHistory; ++u) {
        const int EntryBuf = ((m_Tau - u) % m_BufferSize + m_BufferSize) % m_BufferSize;
        double Incoming = 0.0;

        if (P2 != 0.0) {
            Incoming += P2 * m_States[Prev2].ExitProb[EntryBuf];
        }

        if (P1 != 0.0) {
            Incoming += P1 * m_States[Prev1].ExitProb[EntryBuf];
        }

        const double segment = GetSegmentLikelihood(u, Incoming);
        FTildeJ += surv_cache[u] * segment;
        FTildeJo += occ_cache[u] * segment;

        if (IsLeftRight) {
            m_SegmentLogForwardWeights[u] = LogMultiply(LogProbability(surv_cache[u]), LogProbability(Incoming));
        }
    }

    // A state already active at the start of decoding
    // uses the same observation model.
    const int InitDur = m_Tau + 1;

    if (InitDur <= maxU) {
        const double segment = GetSegmentLikelihood(InitDur, StateJ.InitProb);
        FTildeJ += surv_cache[InitDur] * segment;
        FTildeJo += occ_cache[InitDur] * segment;

        if (IsLeftRight) {
            m_SegmentLogForwardWeights[InitDur] =
                LogMultiply(LogProbability(surv_cache[InitDur]), LogProbability(StateJ.InitProb));
        }
    }

    StateJ.Forward[m_CircularBufferIndex] = FTildeJ + m_SegmentProbabilityFloor;
    StateJ.ExitProb[m_CircularBufferIndex] = FTildeJo + m_SegmentProbabilityFloor;
}

// ─────────────────────────────────────
// GUÉDON (2005) + CUVILLIER (2016)
/**
 * @brief Advance and normalize the decode-window forward probabilities.
 *
 * @return Index of the best eligible internal state.
 *
 * @note Shares normalization with ordered microstate histories and chooses an eligible state without moving
 * backward.
 * @warning Requires current-frame evidence and initialized decoder history.
 */
int OnlineForward::GetAlphaT() {
    spdlog::debug("WinStart {:04d} | WinFinish {:04d} | BufferSize {:04d} | Tau {:06d} | Kappa {:.4f}", m_WinStart,
                  m_WinEnd, m_CircularBufferIndex, m_Tau, m_Kappa);

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        ScoreState &StateJ = m_States[j];
        switch (StateJ.HSMMType) {
        case SEMIMARKOV:
            SemiMarkov(StateJ, j);
            UpdateMicroStatePosterior(StateJ);
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
    const double logN = std::log(N);
    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        m_States[j].Forward[m_CircularBufferIndex] /= N;
        m_States[j].Forward[m_CircularBufferIndex] += std::numeric_limits<double>::min();
        m_States[j].ExitProb[m_CircularBufferIndex] /= N;
        m_States[j].ExitProb[m_CircularBufferIndex] += std::numeric_limits<double>::min();

        if (m_States[j].MicroTopologyType == LEFT_RIGHT) {
            // Shared score-level scaling: each path retains the normalization
            // of its segment history for the next HSMM update.
            for (MarkovMicroState &Micro : m_States[j].MicroStates) {
                for (double &LogAlpha : Micro.LogForwardByAge) {
                    if (LogAlpha != LogZero) {
                        LogAlpha -= logN;
                    }
                }
            }
        }
    }

    // Find the Argmax (Best State)
    double maxVal = std::numeric_limits<double>::min();
    int BestStateIndex = m_CurrentStateIndex;

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        ScoreState &StateJ = m_States[j];
        double fwd = StateJ.Forward[m_CircularBufferIndex];
        bool canAdvance = j >= m_CurrentStateIndex;
        if (m_IsSilence && j > m_CurrentStateIndex) {
            // During silence, cross REST states only if every intervening
            // state is also REST. Sounded skips remain available once audio
            // returns; no internal state is removed from the forward model.
            canAdvance = StateJ.Type == REST;
            for (int k = m_CurrentStateIndex + 1; canAdvance && k < j; ++k) {
                if (m_States[k].Type != REST) {
                    canAdvance = false;
                }
            }
        }
        if (fwd > maxVal && canAdvance) {
            maxVal = fwd;
            BestStateIndex = j;
        }

        spdlog::debug("State ({}) | Obs = {:.5f}, Forward {:.5f}, Exit Prob {:.5f}", StateJ.Index,
                      StateJ.BestObs[m_CircularBufferIndex], StateJ.Forward[m_CircularBufferIndex],
                      StateJ.ExitProb[m_CircularBufferIndex]);
    }

    ScoreState &BestState = m_States[BestStateIndex];
    spdlog::debug("Best: State ({}) | Forward {:.5f}", BestState.Index, BestState.Forward[m_CircularBufferIndex]);

    return BestStateIndex;
}

// ─────────────────────────────────────
/**
 * @brief Decode one audio frame and update the current event and tempo.
 *
 * @param Desc Audio description for the current analysis frame.
 *
 * @return Public score position of the winning internal state.
 *
 * @note Updates observation evidence, performs forward inference, and queues audio-state-change notifications.
 * @warning Requires a nonempty initialized score and a valid current state index; serialize with state changes.
 */
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
