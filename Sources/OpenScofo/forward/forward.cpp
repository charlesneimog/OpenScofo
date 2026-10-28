/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file forward.cpp
 * @brief Online score-inference configuration, initialization, and state management.
 *
 * @note OnlineForward combines spectral pitch evidence, duration distributions, and optional ordered microstate
 * chains across its implementation translation units.
 * @warning The decoder owns mutable history and caches; serialize state changes and inference.
 */
#include "OpenScofo.hpp"
#include <algorithm>

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
/**
 * @brief Initialize decoder defaults and synchronization lookup caches.
 *
 * @note Sets A4 to 440 Hz and precomputes the A2 table and inverse concentration cache.
 */
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
/**
 * @brief Apply processing parameters and rebuild decoder resources.
 *
 * @param Config Audio, pitch-template, synchronization, and section configuration to apply.
 *
 * @note Clears template and duration caches, resizes probability histories, and rebuilds templates for loaded
 * states.
 * @warning Requires a positive sample rate; reallocates histories and must be serialized with inference.
 */
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
    m_SectionRestrict = Config.SectionRestrict;
    m_AudioStateChangeReceiver = Config.AudioStateChangeReceiver;

    m_PitchTemplates.clear();
    m_PitchTemplatesPrecomputed.clear();
    m_PitchCQTTemplates.clear();
    m_OccupancyPMFCache.clear();
    m_SurvivorCache.clear();
    m_InternalLogPathCountCache.clear();

    m_Normalization.assign(static_cast<size_t>(m_BufferSize + 1), 1.0);
    for (ScoreState &State : m_States) {
        State.Forward.assign(static_cast<size_t>(m_BufferSize + 1), std::numeric_limits<double>::min());
        State.ExitProb.assign(static_cast<size_t>(m_BufferSize + 1), std::numeric_limits<double>::min());
        State.BestObs.assign(static_cast<size_t>(m_BufferSize + 1), std::numeric_limits<double>::min());
        ResetMicroStateRuntime(State);
    }

    if (!m_States.empty()) {
        UpdatePitchTemplate();
    }
}

// ─────────────────────────────────────
/**
 * @brief Read the active internal score state index.
 *
 * @return Zero-based active state index.
 *
 * @note Includes boundary and inter-event silence states in the indexing.
 */
int OnlineForward::GetCurrentStateIndex() {
    return m_CurrentStateIndex;
}

// ─────────────────────────────────────
/**
 * @brief Copy the actions for the active score state.
 *
 * @return Copy of the active state action list.
 *
 * @note Reads the state directly without removing any of its actions.
 * @warning Requires a nonempty score and a valid current state index.
 */
EventActions OnlineForward::GetCurrentEventActions() {
    ScoreState State = m_States[m_CurrentStateIndex];
    return State.Actions;
}

// ─────────────────────────────────────
/**
 * @brief Drain the pending audio-state-change action queue.
 *
 * @return Queued actions, in notification order.
 *
 * @note Preserves queue order while moving actions into the result.
 * @warning Retrieval consumes the queued actions.
 */
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
/**
 * @brief Clear decoder caches and microstate runtime history.
 *
 * @note Rehashes several maps to release storage and clears microstate paths.
 * @warning Can allocate or release memory; serialize with inference.
 */
void OnlineForward::ResetCaches() {
    for (ScoreState &State : m_States) {
        ResetMicroStateRuntime(State);
    }
    // Clear all caches
    m_PitchTemplates.clear();
    m_PitchTemplatesPrecomputed.clear();
    m_PitchProbabilityCache.clear();
    m_OccupancyPMFCache.clear();
    m_SurvivorCache.clear();
    m_InternalLogPathCountCache.clear();
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
/**
 * @brief Install score states and initialize decoder histories.
 *
 * @param ScoreStates Score states to install in the decoder.
 *
 * @note Ignores an empty input; otherwise initializes tempo from the starting state and precomputes pitch
 * templates.
 * @warning Replaces the current states and clears pending notifications; serialize with inference.
 */
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
    for (ScoreState &State : m_States) {
        // States entering the moving decode window after the first frame
        // must not retain an uninitialized initial probability.
        State.InitProb = 0.0;
        State.Forward.assign(m_BufferSize, std::numeric_limits<double>::min());
        State.ExitProb.assign(m_BufferSize, std::numeric_limits<double>::min());
        State.BestObs.assign(m_BufferSize, std::numeric_limits<double>::min());
        State.BestAudioStateIndex = -1;
        State.BestMicroStateIndex = -1;
        State.BestMicroObservationIndex = -1;
    }

    m_CurrentStateIndex = 0;
    if (m_SectionRestrict) {
        auto FirstSectionState = std::find_if(m_States.begin(), m_States.end(),
                                              [](const ScoreState &State) { return !State.Section.empty(); });
        if (FirstSectionState != m_States.end()) {
            m_CurrentStateIndex = static_cast<int>(std::distance(m_States.begin(), FirstSectionState));
        }
    }
    const ScoreState &InitialState = m_States[static_cast<size_t>(m_CurrentStateIndex)];
    m_LastNotifiedStateIndex = -1;
    m_LastNotifiedAudioStateIndex = -1;
    m_LastNotifiedMicroStateIndex = -1;
    m_PendingAudioStateActions.clear();
    m_Kappa = 10;
    m_BPM = InitialState.BPMExpected;
    m_PsiN = 60.0f / InitialState.BPMExpected;
    m_PsiN1 = 60.0f / InitialState.BPMExpected;
    m_LastPsiN = 60.0f / InitialState.BPMExpected;
    m_BeatsAhead = InitialState.BPMExpected / 60 * m_SecondsAhead;
    m_SyncStr = 0;
    m_Tau = 0;

    m_SyncStrength = InitialState.SyncStrength;
    m_PhaseCoupling = InitialState.PhaseCoupling;

    UpdatePitchTemplate();
}

// ╭─────────────────────────────────────╮
// │          Set|Get Functions          │
// ╰─────────────────────────────────────╯
/**
 * @brief Remove all installed score states.
 *
 * @warning Load new states before calling methods that index the score vector.
 */
void OnlineForward::ClearStates() {
    m_States.clear();
}

// ─────────────────────────────────────
/**
 * @brief Read the current tempo estimate.
 *
 * @return Tempo in beats per minute.
 *
 * @note Tempo is initialized from the score and updated on decoded event transitions.
 */
double OnlineForward::GetCurrentBPM() {
    return m_BPM;
}

// ─────────────────────────────────────
/**
 * @brief Read the configured analysis hop duration.
 *
 * @return Analysis hop duration in seconds.
 *
 * @note Computed by UpdateConfiguration() from hop size and sampling rate.
 */
double OnlineForward::GetBlockDuration() {
    return m_BlockDur;
}

// ─────────────────────────────────────
/**
 * @brief Set the reference A4 tuning frequency.
 *
 * @param Tunning Reference A4 tuning frequency in Hz.
 *
 * @note Assigns the tuning value without rebuilding existing templates.
 */
void OnlineForward::SetTunning(double Tunning) {
    m_Tunning = Tunning;
}

// ─────────────────────────────────────
/**
 * @brief Store an audio description and invalidate frame evidence caches.
 *
 * @param Desc Audio description for the current analysis frame.
 *
 * @note Copies the description and resets pitch-probability and reverberation cache markers.
 */
void OnlineForward::SetDescription(const Description &Desc) {
    m_Desc = Desc;
    m_PitchProbabilityCache.clear();
    m_PitchProbabilityCacheTau = -1;
    m_ReverbEnergyCacheTau = -1;
}

// ─────────────────────────────────────
/**
 * @brief Set the harmonic count for future pitch templates.
 *
 * @param Harmonics Number of harmonics used when building pitch templates.
 *
 * @note Existing cached templates are retained until rebuilt.
 */
void OnlineForward::SetHarmonics(int Harmonics) {
    m_Harmonics = Harmonics;
}

// ─────────────────────────────────────
/**
 * @brief Set the pitch-dependent harmonic decay exponent.
 *
 * @param decay Exponent shaping harmonic decay across the pitch range.
 *
 * @note The new exponent affects templates built after this call; existing templates remain cached.
 */
void OnlineForward::SetAmplitudeDecay(double decay) {
    m_PitchTemplateAmplitudeDecay = decay;
}

// ─────────────────────────────────────
/**
 * @brief Read the circular slot for the current inference frame.
 *
 * @return Circular history buffer index.
 *
 * @note Computes the elapsed frame count modulo the history buffer size.
 */
int OnlineForward::GetCurrentBufferIndex() {
    return m_Tau % m_BufferSize;
}

// ─────────────────────────────────────
/**
 * @brief Read the reference A4 tuning as an integer.
 *
 * @return Integer A4 tuning frequency in Hz.
 *
 * @note The stored tuning is floating point, but this accessor truncates it to int.
 */
int OnlineForward::GetTunning() {
    return m_Tunning;
}

// ─────────────────────────────────────
/**
 * @brief Select an internal state and restart decoding there.
 *
 * @param Event Zero-based internal state index to select.
 *
 * @note Logs and ignores invalid indices or an empty score; a valid selection calls ResetDecoding().
 */
void OnlineForward::SetCurrentEvent(int Event) {
    spdlog::debug("Current event is {}", Event);
    if (m_States.size() == 0) {
        spdlog::error("There is not Events on Score or the Score was no loaded");
        return;
    }

    if (Event < 0 || Event >= static_cast<int>(m_States.size())) {
        spdlog::error("Event index {} is outside the score", Event);
        return;
    }
    m_CurrentStateIndex = Event;
    ResetDecoding();
}

// ─────────────────────────────────────
/**
 * @brief Select the first state matching a section and restart decoding.
 *
 * @param Section Section name without surrounding quotes.
 *
 * @return True if the section was selected; false otherwise.
 *
 * @note Logs an error and leaves the current index unchanged when no state matches.
 */
bool OnlineForward::SetCurrentSection(const std::string &Section) {
    auto State = std::find_if(m_States.begin(), m_States.end(),
                              [&](const ScoreState &Candidate) { return Candidate.Section == Section; });
    if (State == m_States.end()) {
        spdlog::error("SECTION '{}' was not found in the loaded score", Section);
        return false;
    }

    m_CurrentStateIndex = static_cast<int>(std::distance(m_States.begin(), State));
    ResetDecoding();
    return true;
}

// ─────────────────────────────────────
/**
 * @brief Read the number of installed internal states.
 *
 * @return Number of internal score states.
 *
 * @note Counts musical events as well as boundaries and inserted silence states.
 */
int OnlineForward::GetStatesSize() {
    return m_States.size();
}

// ─────────────────────────────────────
/**
 * @brief Append a score state to the decoder.
 *
 * @param State Score state to read or update.
 *
 * @note Copies the provided state without initializing its probability buffers or rebuilding templates.
 * @warning The caller must prepare runtime buffers before using the appended state in inference.
 */
void OnlineForward::AddState(ScoreState State) {
    m_States.push_back(State);
}

// ─────────────────────────────────────
/**
 * @brief Copy an internal score state by index.
 *
 * @param Index Zero-based internal score state index.
 *
 * @return Copy of the selected state.
 *
 * @note Uses unchecked vector indexing.
 * @warning Index must be within the installed score.
 */
ScoreState OnlineForward::GetState(int Index) {
    return m_States[Index];
}

// ─────────────────────────────────────
/**
 * @brief Access the mutable internal score vector.
 *
 * @return Reference to the installed states.
 *
 * @note Returns the decoder-owned container without copying it.
 * @warning Mutations can invalidate references and decoder invariants; serialize with inference.
 */
std::vector<ScoreState> &OnlineForward::GetStates() {
    return m_States;
}

// ─────────────────────────────────────
/**
 * @brief Set the width of future spectral pitch-template peaks.
 *
 * @param f Pitch-template peak width in semitones.
 *
 * @note Existing cached templates are not rebuilt by this setter.
 */
void OnlineForward::SetPitchTemplateSigma(double f) {
    m_PitchTemplateSigma = f;
}

// ╭─────────────────────────────────────╮
// │            Time Decoding            │
// ╰─────────────────────────────────────╯
/**
 * @brief Initialize timing from the first score state.
 *
 * @note Resets phase, synchronization, onset, tempo, and elapsed-frame counters.
 * @warning Requires a nonempty score with a positive expected BPM in its first state.
 */
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
/**
 * @brief Reset timing and probabilities at the selected score state.
 *
 * @note Clears microstate histories and pending notifications; returns early when no states are installed.
 */
void OnlineForward::ResetDecoding() {
    if (m_States.empty()) {
        return;
    }

    const int StartStateIndex = std::clamp(m_CurrentStateIndex, 0, static_cast<int>(m_States.size()) - 1);
    ScoreState &StartState = m_States[static_cast<size_t>(StartStateIndex)];

    // Reset timing variables
    m_Tau = 0;
    m_TimeInPrevEvent = 0;
    m_CurrentStateOnset = 0;
    m_LastTn = 0;

    // Reset synchronization variables
    m_SyncStr = 0;
    m_Kappa = 10;

    // Reset BPM and period predictions
    double PsiK = 60 / StartState.BPMExpected;
    m_LastPsiN = PsiK;
    m_PsiN = PsiK;
    m_PsiN1 = PsiK;
    m_BPM = StartState.BPMExpected;
    m_BeatsAhead = StartState.BPMExpected / 60.0 * m_SecondsAhead;
    m_SyncStrength = StartState.SyncStrength;
    m_PhaseCoupling = StartState.PhaseCoupling;

    // Reset all state probabilities
    for (ScoreState &State : m_States) {
        std::fill(State.Forward.begin(), State.Forward.end(), 0.0);
        std::fill(State.ExitProb.begin(), State.ExitProb.end(), 0.0);
        std::fill(State.BestObs.begin(), State.BestObs.end(), std::numeric_limits<double>::min());
        State.BestAudioStateIndex = -1;
        State.BestMicroStateIndex = -1;
        State.BestMicroObservationIndex = -1;
        ResetMicroStateRuntime(State);
        State.OnsetObserved = 0;
        State.PhaseObserved = 0;
        State.IOIPhiN = 0;
        State.InitProb = 0.0;
    }

    // Reset normalization buffer
    std::fill(m_Normalization.begin(), m_Normalization.end(), 1.0);

    // Reset first state
    StartState.OnsetObserved = 0;

    m_LastNotifiedStateIndex = -1;
    m_LastNotifiedAudioStateIndex = -1;
    m_LastNotifiedMicroStateIndex = -1;
    m_PendingAudioStateActions.clear();

    spdlog::debug("OnlineForward decoding state fully reset");
}

} // namespace OpenScofo
