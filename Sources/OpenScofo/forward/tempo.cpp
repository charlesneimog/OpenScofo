/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file tempo.cpp
 * @brief Tempo tracking, phase synchronization, and coupling functions.
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

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
/**
 * @brief Evaluate the circular concentration function I1(kappa)/I0(kappa).
 *
 * @param kappa Circular-distribution concentration parameter.
 *
 * @return Circular concentration function value.
 *
 * @note Uses an asymptotic approximation above kappa = 10 and returns zero for nonpositive kappa.
 */
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
/**
 * @brief Populate the concentration lookup table once.
 *
 * @note Subsequent calls return immediately after initialization.
 */
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
/**
 * @brief Read the concentration function from its lookup table.
 *
 * @param kappa Circular-distribution concentration parameter.
 *
 * @return Approximate A2 value.
 *
 * @note Rounds to the nearest table sample and clamps values outside the table range.
 */
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
/**
 * @brief Approximate concentration from synchronization strength.
 *
 * @param SyncStrength Synchronization strength, clamped to the interval [0, 1].
 *
 * @return Approximate concentration parameter.
 *
 * @note Caches quantized strengths and returns 10 for strengths above 0.95.
 */
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
/**
 * @brief Evaluate the phase-coupling correction.
 *
 * @param phi Observed phase in cycles.
 * @param phi_hat Expected phase in cycles.
 * @param kappa Circular-distribution concentration parameter.
 *
 * @return Signed phase-coupling value.
 *
 * @note Uses a periodic phase difference and an exponential concentration term.
 */
double OnlineForward::CouplingFunction(double phi, double phi_hat, double kappa) {
    static constexpr double invTwoPi = 1.0 / (2.0 * std::numbers::pi);
    double diff = 2.0 * std::numbers::pi * (phi - phi_hat);
    double cosDiff = std::cos(diff);
    return invTwoPi * std::exp(kappa * (cosDiff - 1.0)) * std::sin(diff);
}

// ─────────────────────────────────────
// CONT 2010 (Section 7.1)
/**
 * @brief Wrap a phase to the centered unit interval.
 *
 * @param Phase Phase in cycles.
 *
 * @return Wrapped phase in [-0.5, 0.5).
 *
 * @note Preserves phase equivalence modulo one cycle.
 */
double OnlineForward::ModPhases(double Phase) {
    Phase = std::fmod(Phase + 0.5, 1.0);
    if (Phase < 0.0) {
        Phase += 1.0;
    }
    return Phase - 0.5;
}

// ─────────────────────────────────────
// CONT 2010 (Section 5, algorithm 1)
/**
 * @brief Advance frame time and update beat-period predictions on a state change.
 *
 * @param StateIndex Zero-based index in the internal score state vector.
 *
 * @return Predicted beat period in seconds.
 *
 * @note Inter-event silence does not drive tempo updates; musical transitions update future onsets and phase
 * estimates.
 * @warning StateIndex must identify an installed state.
 */
double OnlineForward::UpdatePsiN(int StateIndex) {
    m_TimeInPrevEvent += m_BlockDur;
    m_Tau += 1;

    if (StateIndex == m_CurrentStateIndex || m_States[StateIndex].IsInterEventSilence) {
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

    // Cont (2010), Large and Palmer (1999) and Large and Jones (2002).
    // Only scored events drive phase and tempo updates.
    int PreviousIndex = StateIndex - 1;
    while (PreviousIndex > 0 && m_States[PreviousIndex].IsInterEventSilence) {
        --PreviousIndex;
    }
    int NextIndex = StateIndex + 1;
    while (NextIndex < static_cast<int>(m_States.size()) && m_States[NextIndex].IsInterEventSilence) {
        ++NextIndex;
    }
    ScoreState &LastState = m_States[PreviousIndex];
    ScoreState &CurrentState = m_States[StateIndex];

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
    if (NextIndex < static_cast<int>(m_States.size())) {
        m_States[NextIndex].IOIHatPhiN = PhiN1;
    }

    // Update all next expected onsets, including zero-beat gap states.
    if (StateIndex + 1 < static_cast<int>(m_States.size())) {
        m_States[StateIndex + 1].OnsetExpected = Tn1;
    }
    double LastOnsetExpected = Tn1;

    // the m_CurrentEvent + 1 already updated, now
    // we update the future events to get the Sojourn Time
    for (int i = StateIndex + 2; i < StateIndex + 20; i++) {
        if ((size_t)i >= m_States.size()) {
            break;
        }
        ScoreState &FutureState = m_States[i];
        ScoreState &PreviousFutureState = m_States[(i - 1)];
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

} // namespace OpenScofo
