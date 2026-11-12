#include "reference_test.hpp"

using namespace OpenScofo;
using namespace ReferenceTest;

const fs::path Data = BENCHMARK_DATA_DIR;

// ─────────────────────────────────────
int main() {
    const fs::path fileName = Data / "05-tests" / "descriptor_results.json";

    try {
        const Json data = LoadResults(fileName);
        Configuration config = ReadConfiguration(data);
        config.RequestedDescriptors = {MAGNITUDE, CREST, FLUX, SPREADVARIANCE, SKEWNESS, KURTOSIS, ROLLOFF, ENTROPY};

        // Percentage points: 0.003 means 0.003%, NOT 0.3%.
        // The JSON's configuration.threshold is an absolute tolerance and is
        // deliberately not used for percentage-based comparisons.
        const double tolerancePercent = 0.013;

        if (data.at("configuration").at("essentia_spectrum_source") != "openscofo.magnitude") {
            throw std::runtime_error("Unexpected Essentia spectrum source in JSON");
        }

        Report report;
        const size_t fftSize = static_cast<size_t>(config.FFTSize);
        for (const Json &testCase : data.at("cases")) {
            const long long start = testCase.at("start_sample").get<long long>();
            const std::vector<double> previous = ReadWindow(testCase, "previous", fftSize);
            const std::vector<double> current = ReadWindow(testCase, "current", fftSize);
            const Json &inputs = testCase.at("input");
            const Json &descriptors = testCase.at("descriptors");

            // Python's Essentia references were computed from OpenScofo's raw
            // Magnitude. Verify that spectrum before comparing its descriptors.
            MIR mir;
            mir.UpdateConfiguration(config);
            Description desc{};
            ResizeDescription(config, desc);
            mir.GetDescription(previous, desc);
            Compare(report, start, "magnitude_previous", "Essentia", inputs.at("previous_magnitude"),
                    Json(desc.Magnitude), tolerancePercent);

            mir.GetDescription(current, desc);
            Compare(report, start, "magnitude_current", "Essentia", inputs.at("current_magnitude"),
                    Json(desc.Magnitude), tolerancePercent);

            CompareDescriptor(report, start, descriptors, "essentia", "crest", desc.SpectralCrest, tolerancePercent);
            CompareDescriptor(report, start, descriptors, "essentia", "flux", desc.SpectralFlux, tolerancePercent);
            CompareDescriptor(report, start, descriptors, "essentia", "spread_v", desc.SpectralSpreadVariance,
                              tolerancePercent);
            CompareDescriptor(report, start, descriptors, "essentia", "skewness", desc.SpectralSkewness,
                              tolerancePercent);
            CompareDescriptor(report, start, descriptors, "essentia", "kurtosis", desc.SpectralKurtosis,
                              tolerancePercent);

            // Only tested when uncommented in the Python reference generator.
            CompareDescriptor(report, start, descriptors, "essentia", "rolloff", desc.SpectralRolloff, tolerancePercent,
                              true);
            CompareDescriptor(report, start, descriptors, "essentia", "entropy", desc.SpectralEntropy, tolerancePercent,
                              true);
        }

        PrintSummary(report, "Essentia + spectra");
        return report.Checks > 0 && report.Failed == 0 ? 0 : 1;
    } catch (const std::exception &error) {
        std::cerr << "test_essentia: " << error.what() << '\n';
        return 1;
    }
}
