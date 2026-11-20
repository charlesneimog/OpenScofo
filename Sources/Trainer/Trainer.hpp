/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace OpenScofo::Training {

struct Configuration {
    int SampleRate = 48000;
    std::vector<int> TrainingSampleRates = {44100, 48000};
    int FFTSize = 2048;
    int HopSize = 512;

    double MinDB = -60.0;
    double MaxPitchConfidence = 0.9;
    double TestFraction = 0.1;
    double AugmentationFactor = 1.0;
    int MaxAugmentationLoops = 128;

    int Iterations = 1000;
    double LearningRate = 0.05;
    int Depth = 6;
    int EarlyStoppingRounds = 70;
    int RandomSeed = 80;
    // -1 lets CatBoost select the available CPU threads automatically.
    int ThreadCount = -1;

    std::filesystem::path BasePath;
    std::filesystem::path DatasetFolder;
    std::vector<std::filesystem::path> ImpulseResponseFolders;
    std::filesystem::path CachePath;
    std::vector<std::string> Descriptors;
};

struct SampleOrigin {
    std::uint64_t Recording = 0;
    bool Augmented = false;
};

struct Partition {
    std::vector<float> Features;
    std::vector<std::int32_t> Labels;
    std::vector<SampleOrigin> Origins;
    std::size_t FeatureCount = 0;

    std::size_t SampleCount() const;
    bool Empty() const;
    void Clear();
};

struct Dataset {
    Partition Train;
    Partition Validation;
    Partition Test;
    std::vector<std::string> Labels;
    std::vector<std::string> Descriptors;
    std::map<std::string, std::size_t> OriginalTrainCounts;
    std::map<std::string, std::size_t> AugmentedTrainCounts;
    std::map<std::string, std::size_t> ValidationCounts;
    std::map<std::string, std::size_t> TestCounts;
    bool LoadedFromCache = false;

    void Clear();
};

struct Evaluation {
    double Accuracy = 0.0;
    double BalancedAccuracy = 0.0;
    std::map<std::string, double> ClassAccuracy;
    std::size_t TestedSamples = 0;
};

using ProgressCallback = std::function<void(const std::string &status, double progress)>;
using LogCallback = std::function<void(const std::string &message)>;

class Trainer {
  public:
    Trainer();
    explicit Trainer(const Configuration &configuration);
    ~Trainer();

    Trainer(const Trainer &) = delete;
    Trainer &operator=(const Trainer &) = delete;
    Trainer(Trainer &&) noexcept;
    Trainer &operator=(Trainer &&) noexcept;

    void SetConfiguration(const Configuration &configuration);
    const Configuration &GetConfiguration() const;

    void SetProgressCallback(ProgressCallback callback);
    void SetLogCallback(LogCallback callback);

    void Analyze();
    void Train();
    void ExportModel(const std::filesystem::path &path);
    void ClearData();
    void Cancel();

    bool HasDataset() const;
    bool HasModel() const;
    const Dataset &GetDataset() const;
    const Evaluation &GetEvaluation() const;

  private:
    class Impl;
    std::unique_ptr<Impl> m_Impl;
};

} // namespace OpenScofo::Training
