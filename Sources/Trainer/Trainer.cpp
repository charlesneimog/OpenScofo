/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "Trainer.hpp"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <random>
#include <set>
#include <utility>

namespace OpenScofo {

namespace fs = std::filesystem;

// ─────────────────────────────────────
std::size_t Partition::SampleCount() const {
    return Labels.size();
}

// ─────────────────────────────────────
bool Partition::Empty() const {
    return Labels.empty();
}

// ─────────────────────────────────────
void Partition::Clear() {
    Features.clear();
    Labels.clear();
    Origins.clear();
    FeatureCount = 0;
}

// ─────────────────────────────────────
void Dataset::Clear() {
    Train.Clear();
    Validation.Clear();
    Test.Clear();
    Labels.clear();
    Descriptors.clear();
    OriginalTrainCounts.clear();
    AugmentedTrainCounts.clear();
    ValidationCounts.clear();
    TestCounts.clear();
}

// ─────────────────────────────────────
std::string Trainer::Lowercase(std::string text) const {
    for (char &character : text) {
        character = static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    }
    return text;
}

// ─────────────────────────────────────
bool Trainer::IsAudioFile(const fs::path &path) const {
    const std::string extension = Lowercase(path.extension().string());
    return extension == ".wav" || extension == ".aif" || extension == ".aiff";
}

// ─────────────────────────────────────
std::uint64_t Trainer::HashText(const std::string &text) const {
    std::uint64_t value = 1469598103934665603ULL;
    for (unsigned char character : text) {
        value ^= character;
        value *= 1099511628211ULL;
    }
    return value;
}

// ─────────────────────────────────────
fs::path Trainer::ResolvePath(const fs::path &base, const fs::path &path) const {
    if (path.empty()) {
        return {};
    }
    if (path.is_absolute()) {
        return path.lexically_normal();
    }
    return (base / path).lexically_normal();
}

// ─────────────────────────────────────
bool Trainer::IsVectorDescriptor(Descriptors descriptor) const {
    return descriptor == ::OpenScofo::MFCC || descriptor == ::OpenScofo::CHROMA || descriptor == ::OpenScofo::LOGMEL ||
           descriptor == ::OpenScofo::MAGNITUDE || descriptor == ::OpenScofo::POWERARRAY;
}

// ─────────────────────────────────────
bool Trainer::AppendSample(Partition &partition, const std::vector<float> &features, std::size_t label,
                           std::uint64_t recording, bool augmented) const {
    if (partition.FeatureCount == 0) {
        partition.FeatureCount = features.size();
    }
    if (features.size() != partition.FeatureCount) {
        m_Logger->error("OpenScofo descriptor width changed while extracting the dataset");
        return false;
    }
    partition.Features.insert(partition.Features.end(), features.begin(), features.end());
    partition.Labels.push_back(static_cast<std::int32_t>(label));
    partition.Origins.push_back({recording, augmented});
    return true;
}

// ─────────────────────────────────────
void Trainer::Report(const std::string &status, double progress) const {
    if (m_Progress) {
        m_Progress(status, std::clamp(progress, 0.0, 1.0));
    }
}

// ─────────────────────────────────────
bool Trainer::CheckCancelled() const {
    return !m_Cancelled.load();
}

// ─────────────────────────────────────
fs::path Trainer::BasePath() const {
    fs::path base = m_Config.BasePath.empty() ? fs::current_path() : m_Config.BasePath;
    if (!base.is_absolute()) {
        base = fs::absolute(base);
    }
    return base.lexically_normal();
}

// ─────────────────────────────────────
std::vector<int> Trainer::SampleRates() const {
    std::vector<int> rates = m_Config.TrainingSampleRates;
    if (rates.empty()) {
        rates.push_back(m_Config.SampleRate);
    }
    std::vector<int> unique;
    for (int rate : rates) {
        if (std::find(unique.begin(), unique.end(), rate) == unique.end()) {
            unique.push_back(rate);
        }
    }
    return unique;
}

// ─────────────────────────────────────
std::vector<Descriptors> Trainer::GetDescriptors(std::vector<std::string> &canonical) const {
    canonical.clear();
    if (m_Config.Descriptors.empty()) {
        m_Logger->error("At least one OpenScofo descriptor is required");
        return {};
    }
    ::OpenScofo::OpenScofo descriptorConverter(static_cast<float>(m_Config.SampleRate),
                                               static_cast<float>(m_Config.FFTSize),
                                               static_cast<float>(m_Config.HopSize));
    std::vector<::OpenScofo::Descriptors> descriptors;
    for (const std::string &name : m_Config.Descriptors) {
        const std::string normalized = Lowercase(name);
        const ::OpenScofo::Descriptors descriptor = descriptorConverter.GetDescriptorsEnum(normalized.c_str());
        if (descriptor == ::OpenScofo::INVALID || descriptor == ::OpenScofo::ONNX) {
            m_Logger->error("Unknown or unsupported training descriptor: {}", name);
            canonical.clear();
            return {};
        }
        descriptors.push_back(descriptor);
        canonical.emplace_back(descriptorConverter.GetDescriptionId(descriptor));
    }
    return descriptors;
}

// ─────────────────────────────────────
bool Trainer::Validate() const {
    if (m_Config.SampleRate <= 0 || m_Config.FFTSize < 512 || m_Config.HopSize <= 0 ||
        m_Config.HopSize > m_Config.FFTSize) {
        m_Logger->error("Invalid sample rate, FFT size, or hop size");
        return false;
    }
    for (int rate : m_Config.TrainingSampleRates) {
        if (rate <= 0) {
            m_Logger->error("Training sample rates must be positive");
            return false;
        }
    }
    if (m_Config.TestFraction < 0.0 || m_Config.TestFraction >= 1.0) {
        m_Logger->error("Test fraction must be in the range [0, 1)");
        return false;
    }
    if (m_Config.MaxPitchConfidence < 0.0 || m_Config.MaxPitchConfidence > 1.0) {
        m_Logger->error("Maximum pitch confidence must be between 0 and 1");
        return false;
    }
    if (m_Config.AugmentationFactor < 1.0) {
        m_Logger->error("Augmentation factor must be at least 1");
        return false;
    }
    if (m_Config.MaxAugmentationLoops < 0) {
        m_Logger->error("Maximum augmentation loops cannot be negative");
        return false;
    }
    if (m_Config.Iterations <= 0 || m_Config.LearningRate <= 0.0 || m_Config.Depth <= 0 ||
        m_Config.EarlyStoppingRounds < 0 || m_Config.ThreadCount == 0 || m_Config.ThreadCount < -1) {
        m_Logger->error("Invalid CatBoost configuration");
        return false;
    }
    return true;
}

// ─────────────────────────────────────
bool Trainer::Discover(DiscoveredData &discovered) const {
    discovered = {};
    const fs::path datasetFolder = ResolvePath(BasePath(), m_Config.DatasetFolder);
    if (datasetFolder.empty() || !fs::is_directory(datasetFolder)) {
        m_Logger->error("Training folder is not a directory: {}", datasetFolder.string());
        return false;
    }

    std::vector<fs::path> classFolders;
    for (const fs::directory_entry &entry : fs::directory_iterator(datasetFolder)) {
        if (entry.is_directory()) {
            classFolders.push_back(entry.path());
        }
    }
    std::sort(classFolders.begin(), classFolders.end());
    if (classFolders.size() < 2) {
        m_Logger->error("The training folder must contain at least two class subdirectories");
        return false;
    }

    std::mt19937 random(static_cast<std::mt19937::result_type>(m_Config.RandomSeed));
    for (const fs::path &classFolder : classFolders) {
        const std::size_t label = discovered.Labels.size();
        discovered.Labels.push_back(classFolder.filename().string());

        std::vector<fs::path> files;
        for (const fs::directory_entry &entry : fs::directory_iterator(classFolder)) {
            if (entry.is_regular_file() && IsAudioFile(entry.path())) {
                files.push_back(fs::absolute(entry.path()).lexically_normal());
            }
        }
        if (files.empty()) {
            m_Logger->error("Class folder contains no WAV, AIFF, or AIF files: {}", classFolder.string());
            return false;
        }
        std::sort(files.begin(), files.end());

        std::vector<Recording> recordings;
        for (const fs::path &file : files) {
            recordings.push_back({file, label, HashText(file.generic_string())});
        }
        std::shuffle(recordings.begin(), recordings.end(), random);

        std::size_t testCount = 0;
        std::size_t validationCount = 0;
        if (m_Config.TestFraction > 0.0 && recordings.size() >= 2) {
            testCount = std::max<std::size_t>(1, static_cast<std::size_t>(recordings.size() * m_Config.TestFraction));
            testCount = std::min(testCount, recordings.size() - 1);
        }
        const std::size_t afterTest = recordings.size() - testCount;
        if (m_Config.TestFraction > 0.0 && afterTest >= 2) {
            validationCount =
                std::max<std::size_t>(1, static_cast<std::size_t>(recordings.size() * m_Config.TestFraction));
            validationCount = std::min(validationCount, afterTest - 1);
        }

        discovered.Test.insert(discovered.Test.end(), recordings.begin(),
                               recordings.begin() + static_cast<std::ptrdiff_t>(testCount));
        discovered.Validation.insert(discovered.Validation.end(),
                                     recordings.begin() + static_cast<std::ptrdiff_t>(testCount),
                                     recordings.begin() + static_cast<std::ptrdiff_t>(testCount + validationCount));
        discovered.Train.insert(discovered.Train.end(),
                                recordings.begin() + static_cast<std::ptrdiff_t>(testCount + validationCount),
                                recordings.end());
    }

    for (const fs::path &folderValue : m_Config.ImpulseResponseFolders) {
        const fs::path folder = ResolvePath(BasePath(), folderValue);
        if (!fs::is_directory(folder)) {
            m_Logger->error("Impulse-response folder is not a directory: {}", folder.string());
            return false;
        }
        for (const fs::directory_entry &entry : fs::directory_iterator(folder)) {
            if (entry.is_regular_file() && IsAudioFile(entry.path())) {
                discovered.ImpulseResponses.push_back(fs::absolute(entry.path()).lexically_normal());
            }
        }
    }
    std::sort(discovered.ImpulseResponses.begin(), discovered.ImpulseResponses.end());
    return true;
}

// ─────────────────────────────────────
std::vector<float> Trainer::Flatten(::OpenScofo::OpenScofo &extractor, ::OpenScofo::Description &description,
                                    const std::vector<::OpenScofo::Descriptors> &descriptors) const {
    std::vector<float> values;
    for (::OpenScofo::Descriptors descriptor : descriptors) {
        if (IsVectorDescriptor(descriptor)) {
            const std::vector<double> &array = extractor.GetDescriptionArray(description, descriptor);
            for (double value : array) {
                values.push_back(static_cast<float>(value));
            }
        } else {
            values.push_back(static_cast<float>(extractor.GetDescriptionFloat(description, descriptor)));
        }
    }
    if (values.empty()) {
        m_Logger->error("The selected descriptors produced an empty feature vector");
    }
    return values;
}

// ─────────────────────────────────────
std::optional<std::size_t> Trainer::ProcessSignal(const std::vector<float> &signal, int sampleRate,
                                                  const Recording &recording, Partition &partition,
                                                  const std::vector<::OpenScofo::Descriptors> &descriptors,
                                                  bool pitchFilter, bool augmented, std::size_t maximumSamples) const {
    if (signal.size() < static_cast<std::size_t>(m_Config.FFTSize)) {
        return 0;
    }
    std::vector<::OpenScofo::Descriptors> requested = descriptors;
    if (std::find(requested.begin(), requested.end(), ::OpenScofo::YINCONFIDENCE) == requested.end()) {
        requested.push_back(::OpenScofo::YINCONFIDENCE);
    }
    ::OpenScofo::OpenScofo extractor(static_cast<float>(sampleRate), static_cast<float>(m_Config.FFTSize),
                                     static_cast<float>(m_Config.HopSize));
    extractor.SetRequestedDescriptors(requested);

    std::size_t appended = 0;
    std::vector<float> window(static_cast<std::size_t>(m_Config.FFTSize));
    for (std::size_t start = 0; start + static_cast<std::size_t>(m_Config.FFTSize) <= signal.size();
         start += static_cast<std::size_t>(m_Config.HopSize)) {
        if (!CheckCancelled()) {
            return std::nullopt;
        }
        std::copy_n(signal.begin() + static_cast<std::ptrdiff_t>(start), m_Config.FFTSize, window.begin());
        extractor.ProcessBlock(window.data(), window.size());
        Description description = extractor.GetDescription();
        if (description.dB < m_Config.MinDB) {
            continue;
        }
        if (pitchFilter && description.PitchConfidence >= m_Config.MaxPitchConfidence) {
            continue;
        }
        const std::vector<float> features = Flatten(extractor, description, descriptors);
        if (features.empty() || !AppendSample(partition, features, recording.Label, recording.Id, augmented)) {
            return std::nullopt;
        }
        appended++;
        if (maximumSamples != 0 && appended >= maximumSamples) {
            break;
        }
    }
    return appended;
}

// ─────────────────────────────────────
std::vector<float> Trainer::Augment(const std::vector<float> &input, int sampleRate,
                                    const std::vector<fs::path> &impulseResponses, std::mt19937 &random) const {
    const int choiceCount = impulseResponses.empty() ? 3 : 4;
    std::uniform_int_distribution<int> choice(0, choiceCount - 1);
    const int selected = choice(random);
    if (selected == 0) {
        std::uniform_real_distribution<double> rate(0.7, 1.2);
        return TimeStretch(input, rate(random), m_Config.FFTSize, std::max(1, m_Config.FFTSize / 4));
    }
    if (selected == 1) {
        std::uniform_real_distribution<double> steps(-2.0, 2.0);
        return PitchShift(input, steps(random), m_Config.FFTSize, std::max(1, m_Config.FFTSize / 4));
    }
    if (selected == 2) {
        std::uniform_real_distribution<double> deviation(0.005, 0.01);
        return AddGaussianNoise(input, deviation(random), random);
    }

    std::uniform_int_distribution<std::size_t> impulseIndex(0, impulseResponses.size() - 1);
    std::vector<float> impulse = LoadMono(impulseResponses[impulseIndex(random)], sampleRate);
    return Convolve(input, impulse);
}

// ─────────────────────────────────────
bool Trainer::ExtractRecordings(const std::vector<Recording> &recordings, Partition &partition,
                                const std::vector<::OpenScofo::Descriptors> &descriptors, bool pitchFilter,
                                std::map<std::string, std::size_t> &counts, std::size_t &completed,
                                std::size_t total) const {
    const std::vector<int> rates = SampleRates();
    for (const Recording &recording : recordings) {
        std::size_t recordingSamples = 0;
        for (int sampleRate : rates) {
            if (!CheckCancelled()) {
                return false;
            }
            const std::vector<float> signal = LoadMono(recording.Path, sampleRate);
            if (signal.empty()) {
                return false;
            }
            const std::optional<std::size_t> processed =
                ProcessSignal(signal, sampleRate, recording, partition, descriptors, pitchFilter, false);
            if (!processed) {
                return false;
            }
            recordingSamples += *processed;
            ++completed;
            Report("analyzing", 0.7 * static_cast<double>(completed) / static_cast<double>(total));
        }
        counts[m_Data.Labels[recording.Label]] += recordingSamples;
    }
    return true;
}

// ─────────────────────────────────────
bool Trainer::AssertNoLeakage() const {
    std::set<std::uint64_t> train;
    std::set<std::uint64_t> validation;
    std::set<std::uint64_t> test;
    for (const SampleOrigin &origin : m_Data.Train.Origins) {
        train.insert(origin.Recording);
    }
    for (const SampleOrigin &origin : m_Data.Validation.Origins) {
        if (origin.Augmented) {
            m_Logger->error("Augmented data entered the validation partition");
            return false;
        }
        validation.insert(origin.Recording);
    }
    for (const SampleOrigin &origin : m_Data.Test.Origins) {
        if (origin.Augmented) {
            m_Logger->error("Augmented data entered the test partition");
            return false;
        }
        test.insert(origin.Recording);
    }
    for (std::uint64_t id : train) {
        if (validation.contains(id) || test.contains(id)) {
            m_Logger->error("A recording occurs in both training and evaluation partitions");
            return false;
        }
    }
    for (std::uint64_t id : validation) {
        if (test.contains(id)) {
            m_Logger->error("A recording occurs in both validation and test partitions");
            return false;
        }
    }
    return true;
}

// ─────────────────────────────────────
bool Trainer::Analyze() {
    m_Cancelled.store(false);
    if (!Validate()) {
        return false;
    }
    Report("analyzing", 0.0);

    std::vector<std::string> canonicalDescriptors;
    const std::vector<::OpenScofo::Descriptors> descriptors = GetDescriptors(canonicalDescriptors);
    if (descriptors.empty()) {
        return false;
    }
    DiscoveredData discovered;
    if (!Discover(discovered)) {
        return false;
    }

    m_Data.Clear();
    m_Data.Labels = discovered.Labels;
    m_Data.Descriptors = canonicalDescriptors;
    ClearModel();
    m_Results = {};

    const std::vector<int> rates = SampleRates();
    const std::size_t total = std::max<std::size_t>(
        1, (discovered.Train.size() + discovered.Validation.size() + discovered.Test.size()) * rates.size());
    std::size_t completed = 0;
    if (!ExtractRecordings(discovered.Test, m_Data.Test, descriptors, false, m_Data.TestCounts, completed, total) ||
        !ExtractRecordings(discovered.Validation, m_Data.Validation, descriptors, true, m_Data.ValidationCounts,
                           completed, total) ||
        !ExtractRecordings(discovered.Train, m_Data.Train, descriptors, true, m_Data.OriginalTrainCounts, completed,
                           total)) {
        return false;
    }

    std::size_t largestClass = 0;
    for (const std::string &label : m_Data.Labels) {
        largestClass = std::max(largestClass, m_Data.OriginalTrainCounts[label]);
    }
    const std::size_t target =
        static_cast<std::size_t>(std::ceil(static_cast<double>(largestClass) * m_Config.AugmentationFactor));
    std::mt19937 random(static_cast<std::mt19937::result_type>(m_Config.RandomSeed));
    for (std::size_t label = 0; label < m_Data.Labels.size(); ++label) {
        const std::string &name = m_Data.Labels[label];
        const std::size_t original = m_Data.OriginalTrainCounts[name];
        std::size_t &augmentedCount = m_Data.AugmentedTrainCounts[name];

        std::vector<Recording> sources;
        for (const Recording &recording : discovered.Train) {
            if (recording.Label == label) {
                sources.push_back(recording);
            }
        }
        const std::size_t maximumAttempts = static_cast<std::size_t>(m_Config.MaxAugmentationLoops) * sources.size();
        std::size_t attempts = 0;
        while (original + augmentedCount < target && attempts < maximumAttempts) {
            if (!CheckCancelled()) {
                return false;
            }
            const Recording &recording = sources[attempts % sources.size()];
            const int sampleRate = rates[attempts % rates.size()];
            const std::vector<float> source = LoadMono(recording.Path, sampleRate);
            if (source.empty()) {
                return false;
            }
            const std::vector<float> augmented = Augment(source, sampleRate, discovered.ImpulseResponses, random);
            if (augmented.empty()) {
                return false;
            }
            const std::size_t remaining = target - original - augmentedCount;
            const std::optional<std::size_t> processed =
                ProcessSignal(augmented, sampleRate, recording, m_Data.Train, descriptors, true, true, remaining);
            if (!processed) {
                return false;
            }
            augmentedCount += *processed;
            ++attempts;
            const double classProgress = static_cast<double>(attempts) / static_cast<double>(maximumAttempts);
            Report("analyzing", 0.7 + 0.25 * (static_cast<double>(label) + classProgress) /
                                          static_cast<double>(m_Data.Labels.size()));
        }
        if (original + augmentedCount < target) {
            m_Logger->warn("Class {} did not reach its balanced augmentation target", name);
        }
    }

    for (const std::string &label : m_Data.Labels) {
        if (m_Data.OriginalTrainCounts[label] == 0) {
            m_Logger->error("No training windows survived filtering for class: {}", label);
            return false;
        }
    }
    if (!AssertNoLeakage() || !CheckCancelled()) {
        return false;
    }
    m_Logger->info("Dataset ready: {} classes, {} training, {} validation, {} test windows", m_Data.Labels.size(),
                   m_Data.Train.SampleCount(), m_Data.Validation.SampleCount(), m_Data.Test.SampleCount());
    Report("analyzing", 1.0);
    return true;
}

// ─────────────────────────────────────
Trainer::Trainer(const CatConfiguration &configuration) {
    m_Log = std::make_shared<OpenScofoLog<std::mutex>>();
    m_Log->SetCallback(nullptr, nullptr, &m_HasErrors);
    m_Log->set_pattern("%v");
    std::vector<spdlog::sink_ptr> sinks{m_Log};
#ifndef NDEBUG
    sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
#endif
    m_Logger = std::make_shared<spdlog::logger>("OpenScofoTrainer", sinks.begin(), sinks.end());
#if defined(NDEBUG)
    m_Logger->set_level(spdlog::level::info);
#else
    m_Logger->set_level(spdlog::level::debug);
#endif
    SetConfiguration(configuration);
}

// ─────────────────────────────────────
Trainer::~Trainer() {
    ClearModel();
}

// ─────────────────────────────────────
void Trainer::SetConfiguration(const CatConfiguration &configuration) {
    m_Config = configuration;
    if (m_Config.BasePath.empty()) {
        m_Config.BasePath = fs::current_path();
    }
    m_Data.Clear();
    m_Results = {};
    ClearModel();
}

// ─────────────────────────────────────
const CatConfiguration &Trainer::GetConfiguration() const {
    return m_Config;
}

// ─────────────────────────────────────
void Trainer::SetProgressCallback(std::function<void(const std::string &status, double progress)> callback) {
    m_Progress = std::move(callback);
}

// ─────────────────────────────────────
void Trainer::SetLogCallback(std::function<void(const spdlog::details::log_msg &, void *data)> callback, void *data) {
    m_Log->SetCallback(std::move(callback), data, &m_HasErrors);
}

// ─────────────────────────────────────
void Trainer::SetLogLevel(spdlog::level::level_enum level) {
    m_Logger->set_level(level);
}

// ─────────────────────────────────────
bool Trainer::Train() {
    m_Cancelled.store(false);
    if (!HasDataset() && !Analyze()) {
        return false;
    }
    Report("training", 0.0);
    if (!TrainModel() || !EvaluateModel(m_Results)) {
        return false;
    }
    if (m_Results.TestedSamples == 0) {
        m_Logger->warn("No independent test samples are available; test evaluation was skipped");
    } else {
        m_Logger->info("Test accuracy {:.4f}, balanced accuracy {:.4f} ({} samples)", m_Results.Accuracy,
                       m_Results.BalancedAccuracy, m_Results.TestedSamples);
    }
    return true;
}

// ─────────────────────────────────────
bool Trainer::ExportModel(const fs::path &path) {
    if (!CheckCancelled()) {
        return false;
    }
    Report("exporting", 0.0);
    if (!ExportTrainedModel(ResolvePath(BasePath(), path))) {
        return false;
    }
    Report("exporting", 1.0);
    return true;
}

// ─────────────────────────────────────
void Trainer::ClearData() {
    m_Cancelled.store(false);
    m_Data.Clear();
    m_Results = {};
    ClearModel();
}

// ─────────────────────────────────────
void Trainer::Cancel() {
    m_Cancelled.store(true);
}

// ─────────────────────────────────────
bool Trainer::HasDataset() const {
    return !m_Data.Train.Empty();
}

// ─────────────────────────────────────
bool Trainer::HasModel() const {
    return m_Model != nullptr;
}

// ─────────────────────────────────────
bool Trainer::IsCancelled() const {
    return m_Cancelled.load();
}

// ─────────────────────────────────────
const Dataset &Trainer::GetDataset() const {
    return m_Data;
}

// ─────────────────────────────────────
const Evaluation &Trainer::GetEvaluation() const {
    return m_Results;
}

} // namespace OpenScofo
