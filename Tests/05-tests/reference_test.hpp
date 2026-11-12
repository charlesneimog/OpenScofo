#pragma once

#include <mir/mir.hpp>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <iomanip>

namespace ReferenceTest {

using Json = nlohmann::json;

struct Report {
    int Checks = 0;
    int Failed = 0;
    int MissingOptional = 0;
    double MaximumDifference = 0.0;
    double MaximumPercentError = 0.0;
};

// ─────────────────────────────────────
inline Json LoadResults(const std::filesystem::path &Path) {
    std::ifstream file(Path);
    if (!file) {
        throw std::runtime_error(std::string("Cannot open reference JSON: ") + Path.string());
    }

    Json data = Json::parse(file);
    if (data.at("schema_version").get<int>() != 1) {
        throw std::runtime_error("Unsupported descriptor_results.json schema version");
    }
    if (!data.at("cases").is_array() || data.at("cases").empty()) {
        throw std::runtime_error("Reference JSON has no test cases");
    }
    return data;
}

// ─────────────────────────────────────
inline bool IsClose(double Reference, double Result, double AbsoluteTolerance, double RelativeTolerancePercent) {
    if (!std::isfinite(Reference) || !std::isfinite(Result)) {
        return false;
    }

    const double Difference = std::abs(Reference - Result);
    const double RelativeTolerance = RelativeTolerancePercent / 100.0;

    return Difference <= AbsoluteTolerance + RelativeTolerance * std::abs(Reference);
}

// ─────────────────────────────────────
inline OpenScofo::Configuration ReadConfiguration(const Json &Data) {
    const Json &stored = Data.at("configuration");
    OpenScofo::Configuration config;
    config.SR = stored.at("sample_rate").get<float>();
    config.FFTSize = stored.at("fft_size").get<float>();
    config.HOPSize = stored.at("openscofo_hop").get<float>();

    if (!std::isfinite(config.SR) || config.SR <= 0 || !std::isfinite(config.FFTSize) || config.FFTSize < 512 ||
        !std::isfinite(config.HOPSize) || config.HOPSize <= 0 || std::floor(config.FFTSize) != config.FFTSize) {
        throw std::runtime_error("Invalid sampling configuration in reference JSON");
    }
    // Current reference generator uses these explicit parameters.
    config.MFCCMels = 40;
    config.MFCCCount = 13;
    config.ChromaSize = 12;
    config.ZCRCenter = true;
    config.ZCRPad = false;
    config.ZCRZeroPos = true;
    config.ZCRThreshold = 1e-10;
    config.SpectralRolloffCutoff = 0.85;
    return config;
}

// ─────────────────────────────────────
inline void ResizeDescription(const OpenScofo::Configuration &Config, OpenScofo::Description &Desc) {
    const size_t nBins = static_cast<size_t>(Config.FFTSize) / 2 + 1;
    Desc.Power.resize(nBins);
    Desc.Magnitude.resize(nBins);
    Desc.SpectralMagnitudeNorm.resize(nBins);
    Desc.SpectralMagnitudeFrameNorm.resize(nBins);
    Desc.SpectralPowerFrameNorm.resize(nBins);
    Desc.ReverbSpectralPower.resize(nBins);
    Desc.LogMelSpectrum.resize(static_cast<size_t>(Config.MFCCMels));
    Desc.MFCC.resize(static_cast<size_t>(Config.MFCCCount));
    Desc.Chroma.resize(static_cast<size_t>(Config.ChromaSize));
}

// ─────────────────────────────────────
inline std::vector<double> ReadWindow(const Json &Case, const char *Key, size_t FFTSize) {
    const Json &array = Case.at("input").at(Key);
    if (!array.is_array() || array.size() != FFTSize) {
        throw std::runtime_error(std::string("Input window '") + Key + "' has incorrect length");
    }

    std::vector<double> values;
    values.reserve(array.size());
    for (const Json &value : array) {
        if (!value.is_number()) {
            throw std::runtime_error(std::string("Input window '") + Key + "' contains null/non-numeric data");
        }
        const double number = value.get<double>();
        if (!std::isfinite(number)) {
            throw std::runtime_error(std::string("Input window '") + Key + "' contains NaN/Infinity");
        }
        values.push_back(number);
    }
    return values;
}

// ─────────────────────────────────────
inline void Fail(Report &Summary, long long Start, const char *Name, const std::string &Reason) {
    ++Summary.Checks;
    ++Summary.Failed;
    std::cout << "FAIL | start=" << Start << " | " << Name << " | " << Reason << '\n';
}

// ─────────────────────────────────────

// ─────────────────────────────────────
inline void Compare(Report &Summary, long long Start, const char *Name, const char *Library, const Json &Expected,
                    const Json &Actual, double TolerancePercent, double AbsoluteTolerance = 1e-5) {
    if (Expected.is_array() != Actual.is_array()) {
        Fail(Summary, Start, Name, "scalar/vector shape mismatch");
        return;
    }

    const bool vectorValue = Expected.is_array();
    const size_t count = vectorValue ? Expected.size() : 1;
    const size_t actualCount = vectorValue ? Actual.size() : 1;

    if (count == 0 || count != actualCount) {
        Fail(Summary, Start, Name,
             "length mismatch (expected=" + std::to_string(count) + ", actual=" + std::to_string(actualCount) + ")");
        return;
    }

    double maximumPercent = 0.0;
    double maximumAbsolute = 0.0;
    double differenceAtMaximumPercent = 0.0;
    double referenceAtMaximum = 0.0;
    double actualAtMaximum = 0.0;
    size_t maximumIndex = 0;

    bool passed = true;
    bool hasFailure = false;

    double failurePercent = 0.0;
    double failureDifference = 0.0;
    double failureReference = 0.0;
    double failureActual = 0.0;
    size_t failureIndex = 0;

    for (size_t i = 0; i < count; ++i) {
        const Json &referenceValue = vectorValue ? Expected[i] : Expected;
        const Json &actualValue = vectorValue ? Actual[i] : Actual;

        if (!referenceValue.is_number() || !actualValue.is_number()) {
            Fail(Summary, Start, Name, "null or non-numeric value at index " + std::to_string(i));
            return;
        }

        const double reference = referenceValue.get<double>();
        const double actual = actualValue.get<double>();

        if (!std::isfinite(reference) || !std::isfinite(actual)) {
            Fail(Summary, Start, Name, "non-finite value at index " + std::to_string(i));
            return;
        }

        const double difference = std::abs(reference - actual);

        const double percentage = reference == 0.0 ? (difference == 0.0 ? 0.0 : std::numeric_limits<double>::infinity())
                                                   : 100.0 * (difference / std::abs(reference));

        maximumAbsolute = std::max(maximumAbsolute, difference);

        if (i == 0 || percentage > maximumPercent) {
            maximumPercent = percentage;
            differenceAtMaximumPercent = difference;
            maximumIndex = i;
            referenceAtMaximum = reference;
            actualAtMaximum = actual;
        }

        if (!IsClose(reference, actual, AbsoluteTolerance, TolerancePercent)) {
            passed = false;

            if (!hasFailure) {
                hasFailure = true;
                failurePercent = percentage;
                failureDifference = difference;
                failureReference = reference;
                failureActual = actual;
                failureIndex = i;
            }
        }
    }

    ++Summary.Checks;
    Summary.MaximumDifference = std::max(Summary.MaximumDifference, maximumAbsolute);
    Summary.MaximumPercentError = std::max(Summary.MaximumPercentError, maximumPercent);

    if (!passed) {
        ++Summary.Failed;

        // Display the failing coefficient rather than an unrelated maximum.
        maximumPercent = failurePercent;
        differenceAtMaximumPercent = failureDifference;
        referenceAtMaximum = failureReference;
        actualAtMaximum = failureActual;
        maximumIndex = failureIndex;

        std::cout << "\033[31m";
    }

    std::cout << std::left << std::setw(5) << (passed ? "PASS" : "FAIL") << " " << Library
              << " | start=" << std::setw(8) << Start << " | " << std::setw(10) << Name << " | ref=" << std::setw(13)
              << std::setprecision(8) << referenceAtMaximum << " | mir=" << std::setw(13) << actualAtMaximum
              << " | diff=" << std::setw(13) << differenceAtMaximumPercent << " | error=";

    if (std::isinf(maximumPercent)) {
        std::cout << std::setw(20) << "undefined (ref=0)";
    } else {
        std::ostringstream error;
        error << std::setprecision(10) << maximumPercent << "%";
        std::cout << std::setw(20) << error.str();
    }

    if (vectorValue) {
        std::cout << " | index=" << std::setw(8) << maximumIndex << " | size=" << std::setw(8) << count;
    } else {
        std::cout << " | " << std::setw(14) << ""
                  << " | " << std::setw(13) << "";
    }

    if (!passed) {
        std::cout << "\033[0m";
    }

    std::cout << '\n';
}

// ─────────────────────────────────────
inline void CompareDescriptor(Report &Summary, long long Start, const Json &Descriptors, const char *Library,
                              const char *Name, const Json &Actual, double TolerancePercent, bool Optional = false,
                              double AbsoluteTolerance = 1e-5) {
    if (!Descriptors.contains(Name)) {
        if (Optional) {
            ++Summary.MissingOptional;
        } else {
            Fail(Summary, Start, Name, "reference descriptor missing from JSON");
        }
        return;
    }

    const Json &record = Descriptors.at(Name);

    if (record.at("library") != Library) {
        Fail(Summary, Start, Name, "reference library mismatch");
        return;
    }

    Compare(Summary, Start, Name, Library, record.at("reference"), Actual, TolerancePercent, AbsoluteTolerance);
}

// ─────────────────────────────────────
inline void PrintSummary(const Report &Summary, const char *Library) {
    std::cout << "\n"
              << Library << ": " << (Summary.Checks - Summary.Failed) << " passed, " << Summary.Failed << " failed, "
              << Summary.Checks << " checked"
              << ", max percentage error=" << std::setprecision(10);
    if (std::isinf(Summary.MaximumPercentError)) {
        std::cout << "undefined (ref=0)";
    } else {
        std::cout << Summary.MaximumPercentError << "%";
    }
    std::cout << ", max absolute difference=" << Summary.MaximumDifference << '\n';
    if (Summary.MissingOptional > 0) {
        std::cout << "Optional records not present: " << Summary.MissingOptional
                  << " (rolloff/entropy are disabled in the Python generator)\n";
    }
}

} // namespace ReferenceTest
