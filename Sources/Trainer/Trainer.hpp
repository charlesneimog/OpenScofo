/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#pragma once

#include <OpenScofo.hpp>

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iosfwd>
#include <map>
#include <random>
#include <string>
#include <vector>

struct OpenScofoCatBoostHandle;

namespace OpenScofo {

class Trainer {
  public:
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

    struct AudioBuffer {
        std::vector<float> Samples;
        int SampleRate = 0;
    };

    using ProgressCallback = std::function<void(const std::string &status, double progress)>;
    using LogCallback = std::function<void(const std::string &message)>;

    Trainer();
    explicit Trainer(const Configuration &configuration);
    ~Trainer();

    Trainer(const Trainer &) = delete;
    Trainer &operator=(const Trainer &) = delete;
    Trainer(Trainer &&other) noexcept;
    Trainer &operator=(Trainer &&other) noexcept;

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

    AudioBuffer LoadAudio(const std::filesystem::path &path) const;
    std::vector<float> LoadMono(const std::filesystem::path &path, int sampleRate) const;
    std::vector<float> Resample(const std::vector<float> &input, std::size_t outputSize) const;
    std::vector<float> TimeStretch(const std::vector<float> &input, double rate, int fftSize, int hopSize) const;
    std::vector<float> PitchShift(const std::vector<float> &input, double semitones, int fftSize, int hopSize) const;
    std::vector<float> AddGaussianNoise(const std::vector<float> &input, double deviation, std::mt19937 &random) const;
    std::vector<float> Convolve(const std::vector<float> &input, const std::vector<float> &impulseResponse) const;

  private:
    struct Recording;
    struct DiscoveredData;
    class RealFft;

    std::string Lowercase(std::string text) const;
    bool IsAudioFile(const std::filesystem::path &path) const;
    std::uint64_t HashText(const std::string &text) const;
    std::filesystem::path ResolvePath(const std::filesystem::path &base, const std::filesystem::path &path) const;
    bool IsVectorDescriptor(Descriptors descriptor) const;
    void AppendSample(Partition &partition, const std::vector<float> &features, std::size_t label,
                      std::uint64_t recording, bool augmented) const;

    void WriteBytes(std::ostream &stream, const void *data, std::size_t size) const;
    void ReadBytes(std::istream &stream, void *data, std::size_t size) const;
    template <typename T> void WriteValue(std::ostream &stream, const T &value) const;
    template <typename T> T ReadValue(std::istream &stream) const;
    void WriteString(std::ostream &stream, const std::string &value) const;
    std::string ReadString(std::istream &stream) const;
    void WriteStrings(std::ostream &stream, const std::vector<std::string> &values) const;
    std::vector<std::string> ReadStrings(std::istream &stream) const;
    void WritePartition(std::ostream &stream, const Partition &partition) const;
    Partition ReadPartition(std::istream &stream) const;
    void WriteCounts(std::ostream &stream, const std::map<std::string, std::size_t> &counts) const;
    std::map<std::string, std::size_t> ReadCounts(std::istream &stream) const;

    void Report(const std::string &status, double progress) const;
    void Print(const std::string &message) const;
    void CheckCancelled() const;
    std::filesystem::path BasePath() const;
    std::filesystem::path CacheFile() const;
    std::vector<int> SampleRates() const;
    std::vector<Descriptors> GetDescriptors(std::vector<std::string> &canonical) const;
    void Validate() const;
    DiscoveredData Discover() const;
    std::string Signature(const DiscoveredData &discovered, const std::vector<std::string> &descriptors) const;
    bool LoadCache(const std::string &signature);
    void SaveCache(const std::string &signature) const;
    std::vector<float> Flatten(OpenScofo &extractor, Description &description,
                               const std::vector<Descriptors> &descriptors) const;
    std::size_t ProcessSignal(const std::vector<float> &signal, int sampleRate, const Recording &recording,
                              Partition &partition, const std::vector<Descriptors> &descriptors, bool pitchFilter,
                              bool augmented, std::size_t maximumSamples = 0) const;
    std::vector<float> Augment(const std::vector<float> &input, int sampleRate,
                               const std::vector<std::filesystem::path> &impulseResponses, std::mt19937 &random) const;
    void ExtractRecordings(const std::vector<Recording> &recordings, Partition &partition,
                           const std::vector<Descriptors> &descriptors, bool pitchFilter,
                           std::map<std::string, std::size_t> &counts, std::size_t &completed, std::size_t total) const;
    void AssertNoLeakage() const;

    std::string TrainingParameters(bool haveValidation) const;
    std::string JsonString(const std::string &value) const;
    std::string JsonArray(const std::vector<std::string> &values) const;
    void TrainModel();
    Evaluation EvaluateModel() const;
    void ExportTrainedModel(const std::filesystem::path &path) const;
    void ClearModel() noexcept;

    Configuration m_Config;
    Dataset m_Data;
    Evaluation m_Results;
    ProgressCallback m_Progress;
    LogCallback m_Log;
    std::atomic<bool> m_Cancelled{false};
    OpenScofoCatBoostHandle *m_Model = nullptr;
    std::uint32_t m_CacheVersion = 1;
    std::array<char, 8> m_CacheMagic = {'O', 'S', 'C', 'T', 'R', 'N', '1', '\0'};
};

} // namespace OpenScofo
