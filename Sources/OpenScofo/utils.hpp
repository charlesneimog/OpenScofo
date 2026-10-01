#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <charconv>
#include <string>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │         String Manipulation         │
// ╰─────────────────────────────────────╯
inline bool ParseDouble(const std::string &text, double &value) {
    const char *begin = text.data();
    const char *end = begin + text.size();
    auto [ptr, ec] = std::from_chars(begin, end, value);
    return ec == std::errc{} && ptr == end;
}

// ─────────────────────────────────────
inline bool ParseFloat(std::string_view text, float &value) {
    const char *begin = text.data();
    const char *end = begin + text.size();
    auto [ptr, ec] = std::from_chars(begin, end, value);
    return ec == std::errc{} && ptr == end;
}

// ╭─────────────────────────────────────╮
// │                Math                 │
// ╰─────────────────────────────────────╯
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
