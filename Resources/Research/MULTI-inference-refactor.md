# MULTI inference refactor

The implemented model is the conditional monotonic-path average explicitly requested by the user. All internal endpoints contribute to both parent occupancy and parent exit. There are no internal transition-probability parameters and no microstate duration allocation.

## 1. Previous models

The committed HEAD allocated the parent's expected frame count across microstates using `DurationWeight`, converted each allocation into a geometric advance probability, and made the final microstate absorbing. It already stored age-conditioned histories, but embedded the MULTI computation inside a separate branch of `SemiMarkov()`.

The working tree at the start of this task had removed the duration allocation and replaced the histories with one `LogForward` value per microstate. It advanced exactly one microstate per frame and selected microstate `u-1` as the likelihood for age `u`. It gave zero segment likelihood when `u > K`.

## 2. Problems corrected

The committed model coupled the internal observation trajectory to the parent's duration prior. The initial working tree instead identified segment age with microstate index. Both models conflict with the requested separation between parent occupancy and conditional internal observations. Summing paths without dividing by their count would introduce a third problem: multiplicity alone would increase the observation term.

## 3. New model and sources

For each candidate segment, entry is at microstate zero. A path can stay or advance exactly one microstate, with unit topological weights. The last microstate can stay. Backward moves, skips, and advances beyond the last state are absent. All current endpoints are allowed.

The conditional likelihood is the mean observation product over those paths:

`L(u,t) = sum_k alpha_k(u,t) / C(u,K)`

`C(u,K) = sum_{r=0}^{min(K-1,u-1)} binomial(u-1,r)`

The source responsibilities are deliberately distinct:

- Cont (2010), section 5.2.2 and figure 6, p. 980, specifies the parent semi-Markov event and ordered internal Markov observation states. Section 5.2.1 specifies the separate TRILL framewise maximum. Section 4 presents max-product/Viterbi recursions, rather than the sum-product algorithm used here.
- Guédon (2005), section 4.1, equations (5)-(10), pp. 673-675, supplies the sum-product recursion and interaction of Markov and semi-Markov quantities, with shared observation normalization. His basic hybrid chain has one state space; it does not specify this nested MULTI conditional path model.
- Philippe Cuvillier (2016), appendix A.3.3, pp. 191-192, states the right-censored occupancy and exit recursions, including initial terms and normalization. Survivor weights apply to occupancy and PMF weights apply to exit.
- Unit topological weights, a uniform distribution over paths conditional on segment length, division by `C(u,K)`, and unrestricted internal endpoints are the user's explicit OpenScofo model choices. None is attributed to the papers as a numerical prescription.

## 4. Functions removed

This implementation removes `GetMarkovTransitionProbability()`, which misleadingly described the deterministic internal moves as probabilities. `UpdateMicroStateForward()` and `UpdateMicroStatePosterior()` are rewritten. Observation preparation is separated into ordinary and internal helpers.

`PrepareMicroStateDurations()` was already absent in the initial working tree and remains removed from the final implementation. Score-level `Markov()`, `SemiMarkov()`, and `GetSemiMarkovTransitionProbability()` retain their function bodies from that working tree.

## 5. Fields removed

The scalar `MarkovMicroState::LogForward` is removed in this pass. The duration machinery already absent at the start remains absent: `DurationWeight`, `m_ActiveMarkovScoreStateIndex`, `m_MicroExpectedFrames`, `m_MicroLogSelfProb`, `m_MicroLogAdvanceProb`, `m_MicroLogEmission`, and `m_MicroPosterior`.

## 6. Retained and restored fields

- `MicroStates` and `MicroTopologyType` retain the hierarchy and distinguish ordered MULTI from unordered alternatives.
- Restored `Micro.LogForwardByAge[u]` stores the log sum of observation weights of monotonic paths ending in that microstate, for the segment of age `u` ending at the current frame. This numerator is not divided by `C(u,K)` in the history itself. It carries shared score frame scaling. Age zero is unused; `LogZero` means zero weight.
- `MicroForwardLastFrame` detects missing observation frames; a gap invalidates the old internal trajectories.
- `CurrentEmission` and `BestObservationIndex` keep their existing observation meanings.
- `m_SegmentLogLikelihoods[u]` stores the log conditional segment likelihood, including past frame scaling but excluding the current frame's normalization until the outer update finishes.
- `m_SegmentLogForwardWeights[u]` stores log survivor times incoming/initial mass, used only for reporting.
- `m_SegmentProbabilityFloor` preserves the existing ordinary-state numerical floor.
- New `m_InternalLogPathCountCache[K][u]` stores log `C(u,K)`. It depends only on topology size and segment length, and is cleared with configuration/caches.
- `Forward`, `ExitProb`, `InitProb`, `BestMicroStateIndex`, and `BestMicroObservationIndex` retain their outer inference/reporting roles.

## 7. Final outer recurrence

`PrepareSegmentLikelihoods()` dispatches to the observation helper. `SemiMarkov()` consumes the same prepared likelihood in both survivor and PMF sums. Incoming mass continues to sum predecessor exits through `GetSemiMarkovTransitionProbability()`. The separate initial term at age `t+1` remains. No internal path inference runs inside the outer recurrence.

The complete function is reproduced below from the final source.

## 8. Final internal implementation

`UpdateMicroStateForward()` evaluates the stay-plus-advance sum in log space. Both microstate and age loops run backwards to preserve the previous frame's inputs during the in-place update. `PrepareInternalMarkovSegmentLikelihoods()` log-sums all current endpoints and subtracts the cached log path count. The complete functions are reproduced below.

## 9. Evaluation of a segment of length u

At frame `t`, the age-`u` vector corresponds exactly to observations `x[t-u+1] ... x[t]`. The preceding frame's age-`u-1` vector corresponds to the same entry time. Each update multiplies the current emission by the sum of staying and advancing path weights. Summing all endpoints and dividing by `C(u,K)` gives the requested observation term. No relationship between `u` and a particular microstate is imposed.

## 10. Normalization

Before current-frame outer normalization, a cached numerator contains:

`sum_paths product_{s=t-u+1}^{t} emission(path_s, x_s) / product_{s=t-u+1}^{t-1} N_s`.

After `GetAlphaT()` computes the shared `N_t`, every finite age/microstate log weight subtracts `log(N_t)`. Thus the stored history is ready for the next frame. No age vector is normalized to unit observation mass. Path-count averaging and score normalization are separate operations. Entry mass is combined with the log likelihood before exponentiation to avoid overflowing a standalone scaled likelihood.

## 11. Reporting

For endpoint `k`, `UpdateMicroStatePosterior()` log-sums:

`survivor(u) * incoming(u) * alpha_k(u,t) / C(u,K)`.

The initial hypothesis uses `InitProb` as its incoming weight. `BestMicroStateIndex` is the endpoint with the largest summed occupancy mass, and `BestMicroObservationIndex` is that microstate's existing current emission winner. The common current frame normalization does not change the argmax. Reporting does not modify inference.

## 12. Why age histories remain necessary

Simultaneous parent entry hypotheses see different observation segments and require different endpoint distributions. One scalar per microstate cannot represent them. The implementation stores at most the existing decoder history bound, `m_BufferSize-1` ages, independently of the current duration-cache support. Keeping that history permits later changes in parent duration or tempo to reuse older observed segments. Work and state storage per active MULTI are O(KH), where H is the retained history bound; segment preparation and reporting are O(KU) for the current outer support U.

## 13. OpenScofo-specific choices and preserved policies

The explicit choices introduced or retained for this refactor are: conditional uniform path averaging; unit stay/advance weights; start at microstate zero; all-endpoint occupancy and exit; log-domain arithmetic; caching path counts by K; in-place descending updates; retaining histories to the existing buffer bound; invalidating histories after missing frames; and choosing the reporting winner by summed outer occupancy mass with the existing earliest-index tie behavior.

Existing OpenScofo policies remain: score/parser pitch interpolation; observation definitions and gates; unordered maxima; ordinary observation products and numerical floors; duration distributions, their cache construction and truncation; score initialization; optional-silence branches; decoder window; score-level normalization floors; and actions. These retained policies are not newly attributed to the three references. There are no unresolved transition or endpoint parameters after the user's explicit conditional-model specification. This is an OpenScofo extension, not a claim to reproduce an unspecified numerical MULTI model from Cont.

## 14. Build results

Native Debug build in the existing `build` directory, with Clang, succeeds using `cmake --build build --parallel 2`. It builds the core, enabled Python/PureData/SuperCollider/CSound/Vamp wrappers, tests, and configured performance targets. This does not certify a separate Windows, macOS, Max, or WebAssembly build.

## 15. Test results

All eight registered `openscofo` CTest suites pass through the normal build's `check` target. The MULTI suite contains 33 tests, including six additional tests in this pass and replacements of the deterministic-model expectations.

The checks include ordinary NOTE/REST/CHORD/PTECH/UTECH/TRILL segment products; C4-D4-E4 versus reversed ordering; explicit path enumeration for every entry time; unit likelihood under unit emissions at all tested ages; path counts exceeding linear double range; forbidden skips/backward moves; three microstates over 18 frames; exit before the final microstate; duration-independent internal histories and likelihoods; changed outer occupancy under changed parent duration; reuse of old ages when outer duration support grows; shared normalization against an unnormalized reference; circular histories; reset/missing frames; optional-silence branches; initial contributions; observation winners/actions; and unchanged TRILL emissions.

The previously passing test requiring a one-microstate MULTI to lose all likelihood after one frame is replaced by a long-segment equivalence test against an ordinary state. The synthetic ordered pitch likelihoods are `0.2450745` versus `0.0000745`. A noise-free 18-frame trajectory has one matching path among 154 admissible paths, so its conditional likelihood is `1/154`.

## Final SemiMarkov()

```cpp
void OnlineForward::SemiMarkov(ScoreState &StateJ, int j) {
    const double ExpectedFrames = std::max(1.0, (m_PsiN1 * StateJ.Duration) / m_BlockDur);
    const int key = static_cast<int>(ExpectedFrames * 10.0 + 0.5);
    BuildDistributionCache(ExpectedFrames);
    const auto &surv_cache = m_SurvivorCache[key];
    const auto &occ_cache = m_OccupancyPMFCache[key];
    const int maxU = static_cast<int>(occ_cache.size()) - 1;
    const int observedHistory = std::min(m_Tau, maxU);
    PrepareSegmentLikelihoods(StateJ, std::min(m_Tau + 1, maxU));

    double FTildeJ = 0.0;
    double FTildeJo = 0.0;
    for (int u = 1; u <= observedHistory; ++u) {
        const int EntryBuf = ((m_Tau - u) % m_BufferSize + m_BufferSize) % m_BufferSize;
        double incoming = 0.0;

        for (int i = m_WinStart; i < j; ++i) {
            incoming += GetSemiMarkovTransitionProbability(i, j) * m_States[i].ExitProb[EntryBuf];
        }

        const double segment = GetSegmentLikelihood(u, incoming);
        FTildeJ += surv_cache[u] * segment;
        FTildeJo += occ_cache[u] * segment;
        m_SegmentLogForwardWeights[u] = LogMultiply(LogProbability(surv_cache[u]), LogProbability(incoming));
    }

    // A state already active at the start of decoding uses the same observation model.
    const int initialDuration = m_Tau + 1;
    if (initialDuration <= maxU) {
        const double segment = GetSegmentLikelihood(initialDuration, StateJ.InitProb);
        FTildeJ += surv_cache[initialDuration] * segment;
        FTildeJo += occ_cache[initialDuration] * segment;
        m_SegmentLogForwardWeights[initialDuration] =
            LogMultiply(LogProbability(surv_cache[initialDuration]), LogProbability(StateJ.InitProb));
    }

    StateJ.Forward[m_CircularBufferIndex] = FTildeJ + m_SegmentProbabilityFloor;
    StateJ.ExitProb[m_CircularBufferIndex] = FTildeJo + m_SegmentProbabilityFloor;
}
```

## Observation-model dispatch

```cpp
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
```

## Ordinary segment likelihood

```cpp
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
```

## Internal Forward update

```cpp
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
            double incoming = k == 0 ? 0.0 : LogZero; // Entry is always at microstate zero.
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
```

## Internal path-count cache

```cpp
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
```

## Internal conditional segment likelihood

```cpp
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
```

## Microstate reporting

```cpp
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
```
