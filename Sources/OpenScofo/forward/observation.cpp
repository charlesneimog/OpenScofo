/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file observation.cpp
 * @brief Audio observation likelihoods and emission evaluation.
 */
#include "OpenScofo.hpp"
#include <algorithm>

namespace OpenScofo {

// ─────────────────────────────────────
/**
 * @brief Queue an action when the winning audio observation changes.
 *
 * @param StateIndex Zero-based index in the internal score state vector.
 *
 * @note Uses the configured receiver and suppresses repeated notifications for the same state and observation.
 */
void OnlineForward::NotifyAudioStateChange(int StateIndex) {
    if (m_AudioStateChangeReceiver.empty() || StateIndex < 0 || StateIndex >= static_cast<int>(m_States.size())) {
        return;
    }

    const ScoreState &State = m_States[StateIndex];
    const int MicroStateIndex = State.BestMicroStateIndex;
    const std::vector<Observation> *Observations = &State.Observations;
    int AudioStateIndex = State.BestAudioStateIndex;
    if (State.MicroTopologyType != NO_MICROSTATES) {
        if (MicroStateIndex < 0 || MicroStateIndex >= static_cast<int>(State.MicroStates.size())) {
            return;
        }
        Observations = &State.MicroStates[static_cast<size_t>(MicroStateIndex)].Observations;
        AudioStateIndex = State.BestMicroObservationIndex;
    }

    if (State.Type != CHORD && (AudioStateIndex < 0 || AudioStateIndex >= static_cast<int>(Observations->size()))) {
        return;
    }

    if (StateIndex == m_LastNotifiedStateIndex && AudioStateIndex == m_LastNotifiedAudioStateIndex &&
        MicroStateIndex == m_LastNotifiedMicroStateIndex) {
        return;
    }

    m_LastNotifiedStateIndex = StateIndex;
    m_LastNotifiedAudioStateIndex = AudioStateIndex;
    m_LastNotifiedMicroStateIndex = MicroStateIndex;

    ScoreAction Action{};
    Action.isAudioStateChange = true;
    Action.Receiver = m_AudioStateChangeReceiver;
    Action.AbsoluteTime = true;
    Action.Args.emplace_back(State.ScorePos);

    if (State.Type == CHORD) {
        Action.Args.emplace_back(std::string("chord"));
        for (const Observation &AudioState : State.Observations) {
            if (AudioState.Type == PITCH) {
                Action.Args.emplace_back(static_cast<float>(AudioState.Freq));
            }
        }
    } else {
        const Observation &AudioState = (*Observations)[static_cast<size_t>(AudioStateIndex)];
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

// ╭─────────────────────────────────────╮
// │     Markov / Semi-Markov Core       │
// ╰─────────────────────────────────────╯
/**
 * @brief Evaluate a single audio observation against the current frame.
 *
 * @param Obs Pitch, technique-label, onset, or silence observation to evaluate.
 * @param AllowSilence Whether silence observations may contribute evidence.
 *
 * @return Raw observation evidence.
 *
 * @note Missing ONNX labels and disallowed silence return zero.
 */
double OnlineForward::GetObservationEvidence(const Observation &Obs, bool AllowSilence) {
    switch (Obs.Type) {
    case PITCH:
        return GetPitchProbability(Obs.Freq);
    case LABEL: {
        const auto it = m_Desc.ONNX.find(Obs.Label);
        return it == m_Desc.ONNX.end() ? 0.0 : it->second;
    }
    case ONSET:
        return m_Desc.Onset;
    case SILENCE:
        return AllowSilence ? m_Desc.SilenceProb : 0.0;
    }
    return 0.0;
}

// ─────────────────────────────────────
/**
 * @brief Select the strongest observation evidence for a microstate.
 *
 * @param MicroState Microstate whose observations and winning observation index are evaluated.
 * @param AllowSilence Whether silence observations may contribute evidence.
 * @param ApplyGates Whether to apply sound and technique gates to pitch and label evidence.
 *
 * @return Maximum microstate observation evidence.
 *
 * @note Records the winning observation index and optionally applies sound and technique gates.
 */
double OnlineForward::GetMicroStateEmission(MarkovMicroState &MicroState, bool AllowSilence, bool ApplyGates) {
    const double soundProb = std::max(0.0, 1.0 - m_Desc.SilenceProb);
    double best = 0.0;
    MicroState.BestObservationIndex = -1;
    for (size_t i = 0; i < MicroState.Observations.size(); ++i) {
        const Observation &Obs = MicroState.Observations[i];
        double evidence = GetObservationEvidence(Obs, AllowSilence);
        if (ApplyGates) {
            if (Obs.Type == PITCH) {
                evidence *= soundProb;
            } else if (Obs.Type == LABEL) {
                evidence *= m_Desc.ExtendedTechProb * soundProb;
            }
        }
        if (evidence > best) {
            best = evidence;
            MicroState.BestObservationIndex = static_cast<int>(i);
        }
    }
    return best;
}

// ─────────────────────────────────────
/**
 * @brief Update observation likelihoods within the decode window.
 *
 * @note Combines event and microstate evidence, records winners, and limits REST advancement during silence.
 * @warning Requires initialized decode-window indices and probability history buffers.
 */
void OnlineForward::GetAudioObservations() {
    const double soundProb = std::max(0.0, 1.0 - m_Desc.SilenceProb);
    const double techWeight = m_Desc.ExtendedTechProb;

    EventType CurrentEventType = m_States[m_CurrentStateIndex].Type;

    double maxSoundEvidence = 0.0;

    const bool allowSilence = (CurrentEventType != FIRSTEVENT) && (CurrentEventType != REST);

    // Strongest raw observation across score states.
    double globalBestAudioStateEvidence = 0.0;
    int globalBestStateIndex = -1;
    int globalBestAudioStateIndex = -1;

    for (int j = m_WinStart; j <= m_WinEnd; ++j) {
        ScoreState &state = m_States[j];

        double bestPitch = 0.0;
        double bestTech = 0.0;
        double bestSilence = 0.0;

        double sumPitch = 0.0;
        int pitchCount = 0;
        double bestStateObservationEvidence = 0.0;
        state.BestAudioStateIndex = -1;

        for (size_t audioStateIndex = 0; audioStateIndex < state.Observations.size(); ++audioStateIndex) {

            const Observation &as = state.Observations[audioStateIndex];

            const double audioStateEvidence = GetObservationEvidence(as, state.IsInterEventSilence || allowSilence);
            if (as.Type == PITCH) {
                bestPitch = std::max(bestPitch, audioStateEvidence);
                sumPitch += audioStateEvidence;
                ++pitchCount;
            } else if (as.Type == LABEL) {
                bestTech = std::max(bestTech, audioStateEvidence);
            } else if (as.Type == SILENCE) {
                bestSilence = std::max(bestSilence, audioStateEvidence);
            }

            if (audioStateEvidence > bestStateObservationEvidence) {
                bestStateObservationEvidence = audioStateEvidence;
                state.BestAudioStateIndex = static_cast<int>(audioStateIndex);
            }

            // Keep track of the globally strongest raw AudioState
            // observation. This is intentionally done before the
            // event-level gates below.
            if (audioStateEvidence > globalBestAudioStateEvidence) {
                globalBestAudioStateEvidence = audioStateEvidence;
                globalBestStateIndex = j;
                globalBestAudioStateIndex = static_cast<int>(audioStateIndex);
            }
        }

        double bestMicroLikelihood = 0.0;
        if (state.MicroTopologyType != NO_MICROSTATES) {
            state.BestMicroStateIndex = -1;
            state.BestMicroObservationIndex = -1;
            for (size_t k = 0; k < state.MicroStates.size(); ++k) {
                MarkovMicroState &micro = state.MicroStates[k];
                // TRILL keeps its raw pitch maximum and winner before sound gating.
                micro.CurrentEmission = GetMicroStateEmission(micro, allowSilence, state.Type != TRILL);
                if (micro.CurrentEmission > bestMicroLikelihood) {
                    bestMicroLikelihood = micro.CurrentEmission;
                    if (state.MicroTopologyType == UNORDERED) {
                        state.BestMicroStateIndex = static_cast<int>(k);
                        state.BestMicroObservationIndex = micro.BestObservationIndex;
                    }
                }
                if (micro.CurrentEmission > globalBestAudioStateEvidence) {
                    globalBestAudioStateEvidence = micro.CurrentEmission;
                    globalBestStateIndex = j;
                    globalBestAudioStateIndex = micro.BestObservationIndex;
                }
            }
        }

        double stateLikelihood = 0.0;
        if (state.MicroTopologyType != NO_MICROSTATES) {
            // For LEFT_RIGHT this scalar is diagnostic/silence evidence only.
            // SemiMarkov uses the ordered path for each segment age.
            stateLikelihood = bestMicroLikelihood;
            if (state.Type == TRILL) {
                stateLikelihood = std::max(bestMicroLikelihood * soundProb, bestSilence);
            }
        } else {
            switch (state.Type) {
            case NOTE:
                // Pitch is already an observation probability.
                // Do not gate it using the complement of the
                // extended-technique detector.
                stateLikelihood = bestPitch * soundProb;
                break;

            case PTECH: {
                // Pitched extended technique:
                //
                // The state may be supported independently by:
                //   1. the trained technique classifier, gated by
                //      ExtendedTechProb;
                //   2. its expected pitch.
                //
                // Therefore technique and pitch are not complements.
                const double techObs = bestTech * techWeight * soundProb;

                const double pitchObs = bestPitch * soundProb;

                stateLikelihood = std::max(techObs, pitchObs);

                break;
            }

            case UTECH: {
                // Unpitched extended technique:
                // pitch evidence has no role here.
                const double techObs = bestTech * techWeight * soundProb;

                stateLikelihood = techObs;

                break;
            }

            case CHORD:
                if (pitchCount > 0) {
                    stateLikelihood = (sumPitch / static_cast<double>(pitchCount)) * soundProb;
                }
                break;

            case FIRSTEVENT:
            case REST:
                stateLikelihood = m_Desc.SilenceProb;
                break;

            default:
                spdlog::error("Event type of line {} of score file is not "
                              "implemented, please remove it",
                              state.Line);
                break;
            }
        }

        state.BestObs[m_CircularBufferIndex] = std::max(stateLikelihood, std::numeric_limits<double>::min());

        if (state.Type != REST) {
            maxSoundEvidence = std::max(maxSoundEvidence, stateLikelihood);
        }
    }

    m_IsSilence = (m_Desc.SilenceProb > maxSoundEvidence);

    if (m_IsSilence) {
        // A silent frame can support the next REST, but cannot establish
        // that the performer skipped sounding events to reach a later REST.
        // Keep those hypotheses in the forward window, with no stronger
        // observation than the best sounded event has in this frame.
        bool skippedSoundedState = false;
        const double soundEvidence = std::max(maxSoundEvidence, std::numeric_limits<double>::min());
        for (int j = m_CurrentStateIndex + 1; j <= m_WinEnd; ++j) {
            ScoreState &candidate = m_States[j];
            if (candidate.Type == REST) {
                if (skippedSoundedState) {
                    double &observation = candidate.BestObs[m_CircularBufferIndex];
                    observation = std::min(observation, soundEvidence);
                }
            } else {
                skippedSoundedState = true;
            }
        }
    }

    // Save global winner.
    m_BestAudioStateStateIndex = globalBestStateIndex;
    m_BestAudioStateIndex = globalBestAudioStateIndex;
}

} // namespace OpenScofo
