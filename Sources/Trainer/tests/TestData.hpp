#pragma once

#include <AudioFile.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

class TrainingTempDirectory {
  public:
    TrainingTempDirectory() {
        const auto stamp = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        Path = std::filesystem::temp_directory_path() / ("openscofo-training-test-" + std::to_string(stamp));
        std::filesystem::create_directories(Path);
    }

    ~TrainingTempDirectory() {
        std::error_code error;
        std::filesystem::remove_all(Path, error);
    }

    std::filesystem::path Path;
};

inline std::vector<float> TestTone(int sampleRate, std::size_t size, double frequency, int seed) {
    std::mt19937 random(static_cast<std::mt19937::result_type>(seed));
    std::normal_distribution<float> noise(0.0F, 0.02F);
    std::vector<float> samples(size);
    for (std::size_t i = 0; i < size; ++i) {
        const double phase =
            2.0 * std::numbers::pi * frequency * static_cast<double>(i) / static_cast<double>(sampleRate);
        samples[i] = static_cast<float>(0.45 * std::sin(phase) + 0.15 * std::sin(2.07 * phase)) + noise(random);
    }
    return samples;
}

inline void WriteTestAudio(const std::filesystem::path &path, const std::vector<float> &samples, int sampleRate,
                           bool aiff = false) {
    AudioFile<float> file;
    file.setSampleRate(static_cast<std::uint32_t>(sampleRate));
    file.setBitDepth(16);
    file.setNumChannels(1);
    file.setNumSamplesPerChannel(static_cast<int>(samples.size()));
    file.samples.front() = samples;
    if (!file.save(path.string(), aiff ? AudioFileFormat::Aiff : AudioFileFormat::Wave)) {
        throw std::runtime_error("Could not write test audio file: " + path.string());
    }
}

inline void CreateTrainingDataset(const std::filesystem::path &root, int sampleRate) {
    const std::filesystem::path low = root / "low";
    const std::filesystem::path high = root / "high";
    std::filesystem::create_directories(low);
    std::filesystem::create_directories(high);
    for (int i = 0; i < 3; ++i) {
        WriteTestAudio(low / ("low-" + std::to_string(i) + (i == 1 ? ".aiff" : ".wav")),
                       TestTone(sampleRate, 4096, 220.0 + i * 3.0, 100 + i), sampleRate, i == 1);
        WriteTestAudio(high / ("high-" + std::to_string(i) + (i == 2 ? ".aif" : ".wav")),
                       TestTone(sampleRate, 2048, 1500.0 + i * 7.0, 200 + i), sampleRate, i == 2);
    }
}
