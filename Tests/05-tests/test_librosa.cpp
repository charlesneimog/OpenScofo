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
        config.RequestedDescriptors = {RMS, ZCR, CENTROID, SPREADHZ, FLATNESS, MFCC, CHROMA, LOGMEL};

        // Percentage points: 0.003 means 0.003%, NOT 0.3%.
        // The JSON's configuration.threshold is an absolute tolerance and is
        // deliberately not used for percentage-based comparisons.
        //
        const double tolerancePercent = 0.013;

        Report report;
        const size_t fftSize = static_cast<size_t>(config.FFTSize);
        for (const Json &testCase : data.at("cases")) {
            const long long start = testCase.at("start_sample").get<long long>();
            const std::vector<double> previous = ReadWindow(testCase, "previous", fftSize);
            const std::vector<double> current = ReadWindow(testCase, "current", fftSize);
            const Json &descriptors = testCase.at("descriptors");

            // Fresh MIR per case: frames at random positions do not share history.
            MIR mir;
            mir.UpdateConfiguration(config);
            Description desc{};
            ResizeDescription(config, desc);
            mir.GetDescription(previous, desc);
            mir.GetDescription(current, desc);

            CompareDescriptor(report, start, descriptors, "librosa", "logmel", desc.LogMelSpectrum, tolerancePercent);
            CompareDescriptor(report, start, descriptors, "librosa", "mfcc", desc.MFCC, tolerancePercent);
            CompareDescriptor(report, start, descriptors, "librosa", "flatness", desc.SpectralFlatness,
                              tolerancePercent);
            CompareDescriptor(report, start, descriptors, "librosa", "rms", desc.RMS, tolerancePercent);
            CompareDescriptor(report, start, descriptors, "librosa", "zcr", desc.ZeroCrossingRate, tolerancePercent);
            CompareDescriptor(report, start, descriptors, "librosa", "spreadhz", desc.SpectralSpreadHz,
                              tolerancePercent);
            CompareDescriptor(report, start, descriptors, "librosa", "centroid", desc.SpectralCentroid,
                              tolerancePercent);
            CompareDescriptor(report, start, descriptors, "librosa", "chroma", desc.Chroma, tolerancePercent);
        }

        PrintSummary(report, "Librosa");
        return report.Checks > 0 && report.Failed == 0 ? 0 : 1;
    } catch (const std::exception &error) {
        std::cerr << "test_librosa: " << error.what() << '\n';
        return 1;
    }
}
