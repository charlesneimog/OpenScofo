#include <gtest/gtest.h>
#include <OpenScofo.hpp>

#include <stdexcept>

double ExplicitOccupancy(double r, double p, int u) {
    if (p <= 0.0 || p >= 1.0) {
        throw std::invalid_argument("Erro: 'p' deve estar entre 0 e 1.");
    }
    if (r <= 0.0) {
        throw std::domain_error("Erro: 'r' deve ser maior que 0.");
    }

    double logCoeff = std::lgamma(u + r) - (std::lgamma(u + 1.0) + std::lgamma(r));
    double logProb = logCoeff + r * std::log(p) + static_cast<double>(u) * std::log(1.0 - p);
    return std::exp(logProb);
}

TEST(MathTest, OccupancyThrowsOnInvalidP) {
    // Verifica se a função realmente lança um invalid_argument
    EXPECT_THROW(ExplicitOccupancy(12.4, 1.5, 2), std::invalid_argument);
}
