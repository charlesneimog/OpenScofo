/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file microstates.cpp
 * @brief Ordered microstate inference and segment likelihoods.
 */
#include "OpenScofo.hpp"
#include <algorithm>

namespace OpenScofo {

// ─────────────────────────────────────
/**
 * @brief Clear microstate emissions, winners, and age-conditioned history.
 *
 * @param State Score state to read or update.
 *
 * @note Retains the score-defined observations and microstate topology.
 */
void OnlineForward::ResetMicroStateRuntime(ScoreState &State) {
    State.MicroForwardLastFrame = -1;
    State.BestMicroStateIndex = -1;
    State.BestMicroObservationIndex = -1;
    for (MarkovMicroState &Micro : State.MicroStates) {
        Micro.CurrentEmission = 0.0;
        Micro.BestObservationIndex = -1;
        Micro.LogForwardByAge.clear();
    }
}

// ─────────────────────────────────────
/**
 * @brief Advance the internal monotonic path sums for every parent entry time.
 *
 * @note Cont (2010), section 5.2.2, specifies the MULTI hierarchy and ordering.
 * OpenScofo's conditional path model assigns unit topological weights to staying
 * and advancing one microstate; these are not transition probabilities. The
 * forward sum-product structure is that of Guédon (2005), equation (7).
 * Histories retain outer frame scaling and no parent duration or tempo enters
 * this update. Descending microstates and ages preserve previous-frame inputs.
 */
void OnlineForward::UpdateMicroStateForward(ScoreState &State) {
    if (State.MicroTopologyType != LEFT_RIGHT) {
        return;
    }

    // Retain the decoder's available history independently of the current
    // occupancy support, so later duration/tempo changes can reuse older ages.
    const int MaxAge = std::min(m_Tau + 1, m_BufferSize - 1);
    for (MarkovMicroState &Micro : State.MicroStates) {
        if (State.MicroForwardLastFrame != m_Tau - 1) {
            Micro.LogForwardByAge.clear();
        }
        Micro.LogForwardByAge.resize(static_cast<size_t>(MaxAge + 1), LogZero);
    }

    for (size_t k = State.MicroStates.size(); k-- > 0;) {
        MarkovMicroState &Micro = State.MicroStates[k];
        const double logEmission = LogProbability(Micro.CurrentEmission);

        for (int u = MaxAge; u >= 1; --u) {
            double incoming = k == 0 ? 0.0 : LogZero;
            if (u > 1) {
                incoming = Micro.LogForwardByAge[u - 1]; // Stay, including at the final microstate.
                if (k > 0) {
                    incoming = LogAdd(incoming, State.MicroStates[k - 1].LogForwardByAge[u - 1]);
                }
            }
            Micro.LogForwardByAge[u] = LogMultiply(logEmission, incoming);
        }
    }
    State.MicroForwardLastFrame = m_Tau;
}

// ─────────────────────────────────────
/**
 * @brief Prepare normalized observation likelihoods for each observed segment age.
 *
 * @note Includes the current emission exactly once. Past ordinary observations use
 * stored frame normalization; ordered histories already carry that same scaling.
 * Only the observation model dispatches on microstate topology.
 */
void OnlineForward::PrepareSegmentLikelihoods(ScoreState &State, int MaxAge) {
    m_SegmentLogLikelihoods.assign(static_cast<size_t>(MaxAge + 1), LogZero);
    m_SegmentLogForwardWeights.assign(static_cast<size_t>(MaxAge + 1), LogZero);
    m_SegmentProbabilityFloor = 0.0;

    if (State.MicroTopologyType == LEFT_RIGHT) {
        PrepareInternalMarkovSegmentLikelihoods(State, MaxAge);
    } else {
        PrepareOrdinarySegmentLikelihoods(State, MaxAge);
    }
}

// ─────────────────────────────────────
/**
 * @brief Prepare the stationary observation product with existing outer scaling.
 */
void OnlineForward::PrepareOrdinarySegmentLikelihoods(const ScoreState &State, int MaxAge) {
    const double Bj = State.BestObs[m_CircularBufferIndex];
    const double logBj = LogProbability(Bj);
    // Preserve the ordinary recurrence's Bj * min() numerical floor.
    m_SegmentProbabilityFloor = Bj * std::numeric_limits<double>::min();
    double ObsProd = 1.0;
    for (int u = 1; u <= MaxAge; ++u) {
        m_SegmentLogLikelihoods[u] = LogMultiply(logBj, LogProbability(ObsProd));
        if (u == MaxAge) {
            break;
        }
        const int PreviousBuf = ((m_Tau - u) % m_BufferSize + m_BufferSize) % m_BufferSize;
        const double previousNormalization = m_Normalization[PreviousBuf];
        if (previousNormalization > std::numeric_limits<double>::min()) {
            ObsProd *= State.BestObs[PreviousBuf] / previousNormalization;
        } else {
            ObsProd = 0.0;
        }
        if (ObsProd == 0.0) {
            break;
        }
    }
}

// ─────────────────────────────────────
/**
 * @brief Cache log C(u, K) for the conditional monotonic path model.
 *
 * @note OpenScofo implementation choice: average uniformly over all admissible
 * paths for the candidate segment length. C(u, K) = sum_r binomial(u-1, r),
 * 0 <= r <= min(K-1, u-1). This is not specified by the three reference papers.
 * Counts never use emissions, score normalization, parent duration, or tempo.
 */
const std::vector<double> &OnlineForward::GetInternalLogPathCounts(size_t MicroStateCount, int MaxAge) {
    auto &counts = m_InternalLogPathCountCache[MicroStateCount];
    if (counts.empty()) {
        counts.push_back(LogZero); // Age zero is unused.
    }
    for (int u = static_cast<int>(counts.size()); u <= MaxAge; ++u) {
        double logCount = LogZero;
        double logBinomial = 0.0;
        for (size_t r = 0; r < MicroStateCount && r < static_cast<size_t>(u); ++r) {
            if (r > 0) {
                logBinomial += std::log(static_cast<double>(u - r)) - std::log(static_cast<double>(r));
            }
            logCount = LogAdd(logCount, logBinomial);
        }
        counts.push_back(logCount);
    }
    return counts;
}

// ─────────────────────────────────────
/**
 * @brief Prepare MULTI segment likelihoods by averaging over monotonic paths.
 *
 * @note OpenScofo's endpoint policy permits ending in any current microstate
 * for both parent occupancy and exit. Divide the summed forward weights by
 * C(u, K), not by their observation-dependent sum: likelihood magnitude must
 * remain available to the outer HSMM. Unit emissions therefore give unit
 * conditional likelihood at every age, regardless of the number of paths.
 */
void OnlineForward::PrepareInternalMarkovSegmentLikelihoods(ScoreState &State, int MaxAge) {
    UpdateMicroStateForward(State);
    if (State.MicroStates.empty()) {
        return;
    }
    const auto &counts = GetInternalLogPathCounts(State.MicroStates.size(), MaxAge);
    for (int u = 1; u <= MaxAge; ++u) {
        double logLikelihood = LogZero;
        for (const MarkovMicroState &Micro : State.MicroStates) {
            logLikelihood = LogAdd(logLikelihood, Micro.LogForwardByAge[u]);
        }
        if (logLikelihood != LogZero) {
            m_SegmentLogLikelihoods[u] = logLikelihood - counts[u];
        }
    }
}

// ─────────────────────────────────────
/**
 * @brief Weight a prepared segment likelihood by its entry probability.
 *
 * @note Combines entry mass before exponentiation so scaled microstate histories
 * need not be representable as standalone linear probabilities.
 */
double OnlineForward::GetSegmentLikelihood(int Age, double IncomingProbability) const {
    return std::exp(LogMultiply(m_SegmentLogLikelihoods[Age], LogProbability(IncomingProbability)));
}

// ─────────────────────────────────────
/**
 * @brief Select the ordered microstate winner from the outer HSMM occupancy mass.
 *
 * @note Marginalizes each endpoint over all parent entry hypotheses, using
 * the same path-count divisor, survivor and incoming mass as SemiMarkov.
 * Comparing log masses affects only reporting, never inference.
 */
void OnlineForward::UpdateMicroStatePosterior(ScoreState &State) {
    if (State.MicroTopologyType != LEFT_RIGHT) {
        return;
    }
    State.BestMicroStateIndex = -1;
    State.BestMicroObservationIndex = -1;
    const int MaxAge = static_cast<int>(m_SegmentLogForwardWeights.size()) - 1;
    const auto &counts = GetInternalLogPathCounts(State.MicroStates.size(), MaxAge);
    double best = LogZero;
    for (size_t k = 0; k < State.MicroStates.size(); ++k) {
        double mass = LogZero;
        for (int u = 1; u <= MaxAge; ++u) {
            double contribution = LogMultiply(State.MicroStates[k].LogForwardByAge[u], m_SegmentLogForwardWeights[u]);
            if (contribution != LogZero) {
                mass = LogAdd(mass, contribution - counts[u]);
            }
        }
        if (mass > best) {
            best = mass;
            State.BestMicroStateIndex = static_cast<int>(k);
            State.BestMicroObservationIndex = State.MicroStates[k].BestObservationIndex;
        }
    }
}

} // namespace OpenScofo
