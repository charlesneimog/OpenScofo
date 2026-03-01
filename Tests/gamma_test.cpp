#include <iostream>
#include <cmath>
#include <iomanip>

// Method 1: The original way using lgamma and exp
double ExplicitOccupancy(double r, double p, int u) {
    double logCoeff = std::lgamma(u + r) - (std::lgamma(u + 1.0) + std::lgamma(r));
    double logProb = logCoeff + r * std::log(p) + static_cast<double>(u) * std::log(1.0 - p);
    return std::exp(logProb);
}

int main() {
    double r = 12.4; // Example expected frames
    double p = 0.5;
    int maxU = 10; // Test up to u = 10

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "u | Explicit (lgamma) | Recursive (Math) | Difference\n";
    std::cout << "---------------------------------------------------------\n";

    // Method 2: The optimized recursive way
    double current_recursive = std::pow(p, r); // Base case for u = 0

    for (int u = 0; u <= maxU; u++) {
        // Calculate explicitly for comparison
        double current_explicit = ExplicitOccupancy(r, p, u);

        // Print results
        double diff = std::abs(current_explicit - current_recursive);
        std::cout << std::left << std::setw(2) << u << "| " << std::setw(18) << current_explicit << "| "
                  << std::setw(17) << current_recursive << "| " << diff << "\n";

        // Calculate the next recursive value for the next loop iteration
        current_recursive = current_recursive * ((u + 1.0 + r - 1.0) / (u + 1.0)) * (1.0 - p);
    }

    return 0;
}
