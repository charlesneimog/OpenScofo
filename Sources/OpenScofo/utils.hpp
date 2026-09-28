#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace OpenScofo {

// ─────────────────────────────────────
inline constexpr double LogZero = -std::numeric_limits<double>::max();

// ─────────────────────────────────────
inline double LogProbability(double Probability) {
    return Probability > 0.0 ? std::log(Probability) : LogZero;
}

// ─────────────────────────────────────
inline double LogMultiply(double A, double B) {
    return A == LogZero || B == LogZero ? LogZero : A + B;
}

// ─────────────────────────────────────
inline double LogAdd(double A, double B) {
    if (A == LogZero) {
        return B;
    }

    if (B == LogZero) {
        return A;
    }

    const double High = std::max(A, B);
    const double Low = std::min(A, B);

    return High + std::log1p(std::exp(Low - High));
}

} // namespace OpenScofo
