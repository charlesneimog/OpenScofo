/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#pragma once

#include <OpenScofo.hpp>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

struct OpenScofoCatBoostHandle;
struct PFFFT_Setup;

namespace OpenScofo {

struct CatConfiguration {
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

struct Recording {
    std::filesystem::path Path;
    std::size_t Label = 0;
    std::uint64_t Id = 0;
};

struct DiscoveredData {
    std::vector<std::string> Labels;
    std::vector<Recording> Train;
    std::vector<Recording> Validation;
    std::vector<Recording> Test;
    std::vector<std::filesystem::path> ImpulseResponses;
};

// Owns a PFFFT real-FFT setup and its work buffer. Implemented in Audio.cpp.
class RealFft {
  public:
    RealFft(int size);

    std::unique_ptr<float, void (*)(void *)> AllocateBuffer() const;
    PFFFT_Setup *GetSetup();
    void Forward(const float *input, float *output);
    void Backward(const float *input, float *output);

  private:
    int m_Size;
    std::unique_ptr<PFFFT_Setup, void (*)(PFFFT_Setup *)> m_Setup;
    std::unique_ptr<float, void (*)(void *)> m_Work;
};

class Trainer {
  public:
    Trainer(const CatConfiguration &configuration = CatConfiguration());
    ~Trainer();

    Trainer(Trainer &&other) noexcept;
    Trainer &operator=(Trainer &&other) noexcept;

    void SetConfiguration(const CatConfiguration &configuration);
    const CatConfiguration &GetConfiguration() const;

    void SetProgressCallback(std::function<void(const std::string &status, double progress)> callback);

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
    std::string Lowercase(std::string text) const;
    bool IsAudioFile(const std::filesystem::path &path) const;
    std::uint64_t HashText(const std::string &text) const;
    std::filesystem::path ResolvePath(const std::filesystem::path &base, const std::filesystem::path &path) const;
    bool IsVectorDescriptor(Descriptors descriptor) const;
    void AppendSample(Partition &partition, const std::vector<float> &features, std::size_t label,
                      std::uint64_t recording, bool augmented) const;

    void Report(const std::string &status, double progress) const;
    void CheckCancelled() const;
    std::filesystem::path BasePath() const;
    std::vector<int> SampleRates() const;
    std::vector<Descriptors> GetDescriptors(std::vector<std::string> &canonical) const;
    bool Validate() const;
    DiscoveredData Discover() const;
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
    static int TrainingProgress(void *data, std::size_t iteration, std::size_t iterations);
    void TrainModel();
    Evaluation EvaluateModel() const;
    void ExportTrainedModel(const std::filesystem::path &path) const;
    void ClearModel() noexcept;

    CatConfiguration m_Config;

    Dataset m_Data;
    Evaluation m_Results;
    std::function<void(const std::string &status, double progress)> m_Progress;
    std::atomic<bool> m_Cancelled{false};
    OpenScofoCatBoostHandle *m_Model = nullptr;
};

} // namespace OpenScofo
