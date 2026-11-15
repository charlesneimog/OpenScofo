#pragma once

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

namespace OpenScofo {

// ╭─────────────────────────────────────╮
// │         String Manipulation         │
// ╰─────────────────────────────────────╯
inline bool ParseDouble(const std::string &text, double &value) {
    if (text.empty() || std::isspace(static_cast<unsigned char>(text.front()))) {
        return false;
    }

    char *end;
    errno = 0;
    value = std::strtod(text.c_str(), &end);
    return errno != ERANGE && end == text.c_str() + text.size();
}

// ─────────────────────────────────────
inline bool ParseFloat(const std::string &text, float &value) {
    if (text.empty() || std::isspace(static_cast<unsigned char>(text.front()))) {
        return false;
    }

    char *end;
    errno = 0;
    value = std::strtof(text.c_str(), &end);
    return errno != ERANGE && end == text.c_str() + text.size();
}

// ╭─────────────────────────────────────╮
// │                Math                 │
// ╰─────────────────────────────────────╯
inline constexpr double LogZero = -std::numeric_limits<double>::infinity();

// ─────────────────────────────────────
inline double LogProbability(double Probability) {
    return Probability > 0.0 ? std::log(Probability) : LogZero;
}

// ─────────────────────────────────────
inline double LogMultiply(double A, double B) {
    return A + B;
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
