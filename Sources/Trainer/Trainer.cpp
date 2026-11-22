/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "Trainer.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace OpenScofo {

namespace fs = std::filesystem;

// ─────────────────────────────────────
struct Trainer::DiscoveredData {
    std::vector<std::string> Labels;
    std::vector<Recording> Train;
    std::vector<Recording> Validation;
    std::vector<Recording> Test;
    std::vector<fs::path> AudioFiles;
    std::vector<fs::path> ImpulseResponses;
};

// ─────────────────────────────────────
std::size_t Trainer::Partition::SampleCount() const {
    return Labels.size();
}

// ─────────────────────────────────────
bool Trainer::Partition::Empty() const {
    return Labels.empty();
}

// ─────────────────────────────────────
void Trainer::Partition::Clear() {
    Features.clear();
    Labels.clear();
    Origins.clear();
    FeatureCount = 0;
}

// ─────────────────────────────────────
void Trainer::Dataset::Clear() {
    Train.Clear();
    Validation.Clear();
    Test.Clear();
    Labels.clear();
    Descriptors.clear();
    OriginalTrainCounts.clear();
    AugmentedTrainCounts.clear();
    ValidationCounts.clear();
    TestCounts.clear();
    LoadedFromCache = false;
}

// ─────────────────────────────────────
struct Trainer::Recording {
    fs::path Path;
    std::size_t Label = 0;
    std::uint64_t Id = 0;
};

// ─────────────────────────────────────
std::string Trainer::Lowercase(std::string text) const {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
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
void Trainer::AppendSample(Partition &partition, const std::vector<float> &features, std::size_t label,
                           std::uint64_t recording, bool augmented) const {
    if (partition.FeatureCount == 0) {
        partition.FeatureCount = features.size();
    }
    if (features.size() != partition.FeatureCount) {
        throw std::runtime_error("OpenScofo descriptor width changed while extracting the dataset");
    }
    partition.Features.insert(partition.Features.end(), features.begin(), features.end());
    partition.Labels.push_back(static_cast<std::int32_t>(label));
    partition.Origins.push_back({recording, augmented});
}

// ─────────────────────────────────────
void Trainer::WriteBytes(std::ostream &stream, const void *data, std::size_t size) const {
    if (size == 0) {
        return;
    }
    stream.write(static_cast<const char *>(data), static_cast<std::streamsize>(size));
    if (!stream) {
        throw std::runtime_error("Could not write the OpenScofo dataset cache");
    }
}

// ─────────────────────────────────────
void Trainer::ReadBytes(std::istream &stream, void *data, std::size_t size) const {
    if (size == 0) {
        return;
    }
    stream.read(static_cast<char *>(data), static_cast<std::streamsize>(size));
    if (!stream) {
        throw std::runtime_error("OpenScofo dataset cache is truncated");
    }
}

// ─────────────────────────────────────
template <typename T> void Trainer::WriteValue(std::ostream &stream, const T &value) const {
    WriteBytes(stream, &value, sizeof(value));
}

// ─────────────────────────────────────
template <typename T> T Trainer::ReadValue(std::istream &stream) const {
    T value{};
    ReadBytes(stream, &value, sizeof(value));
    return value;
}

// ─────────────────────────────────────
void Trainer::WriteString(std::ostream &stream, const std::string &value) const {
    const std::uint64_t size = value.size();
    WriteValue(stream, size);
    WriteBytes(stream, value.data(), value.size());
}

// ─────────────────────────────────────
std::string Trainer::ReadString(std::istream &stream) const {
    const std::uint64_t size = ReadValue<std::uint64_t>(stream);
    if (size > 16ULL * 1024ULL * 1024ULL) {
        throw std::runtime_error("OpenScofo dataset cache contains an invalid string size");
    }
    std::string value(static_cast<std::size_t>(size), '\0');
    ReadBytes(stream, value.data(), value.size());
    return value;
}

// ─────────────────────────────────────
void Trainer::WriteStrings(std::ostream &stream, const std::vector<std::string> &values) const {
    WriteValue(stream, static_cast<std::uint64_t>(values.size()));
    for (const std::string &value : values) {
        WriteString(stream, value);
    }
}

// ─────────────────────────────────────
std::vector<std::string> Trainer::ReadStrings(std::istream &stream) const {
    const std::uint64_t size = ReadValue<std::uint64_t>(stream);
    if (size > 100000) {
        throw std::runtime_error("OpenScofo dataset cache contains too many strings");
    }
    std::vector<std::string> values;
    values.reserve(static_cast<std::size_t>(size));
    for (std::uint64_t i = 0; i < size; ++i) {
        values.push_back(ReadString(stream));
    }
    return values;
}

// ─────────────────────────────────────
void Trainer::WritePartition(std::ostream &stream, const Partition &partition) const {
    WriteValue(stream, static_cast<std::uint64_t>(partition.FeatureCount));
    WriteValue(stream, static_cast<std::uint64_t>(partition.SampleCount()));
    WriteBytes(stream, partition.Features.data(), partition.Features.size() * sizeof(float));
    WriteBytes(stream, partition.Labels.data(), partition.Labels.size() * sizeof(std::int32_t));
    for (const SampleOrigin &origin : partition.Origins) {
        WriteValue(stream, origin.Recording);
        WriteValue(stream, static_cast<std::uint8_t>(origin.Augmented ? 1 : 0));
    }
}

// ─────────────────────────────────────
Trainer::Partition Trainer::ReadPartition(std::istream &stream) const {
    Partition partition;
    const std::uint64_t featureCount = ReadValue<std::uint64_t>(stream);
    const std::uint64_t sampleCount = ReadValue<std::uint64_t>(stream);
    if (featureCount > 1000000 || sampleCount > 100000000 ||
        (featureCount != 0 && sampleCount > std::numeric_limits<std::uint64_t>::max() / featureCount)) {
        throw std::runtime_error("OpenScofo dataset cache contains invalid dimensions");
    }

    partition.FeatureCount = static_cast<std::size_t>(featureCount);
    partition.Features.resize(static_cast<std::size_t>(featureCount * sampleCount));
    partition.Labels.resize(static_cast<std::size_t>(sampleCount));
    partition.Origins.resize(static_cast<std::size_t>(sampleCount));

    ReadBytes(stream, partition.Features.data(), partition.Features.size() * sizeof(float));
    ReadBytes(stream, partition.Labels.data(), partition.Labels.size() * sizeof(std::int32_t));

    for (SampleOrigin &origin : partition.Origins) {
        origin.Recording = ReadValue<std::uint64_t>(stream);
        origin.Augmented = ReadValue<std::uint8_t>(stream) != 0;
    }

    return partition;
}

// ─────────────────────────────────────
void Trainer::WriteCounts(std::ostream &stream, const std::map<std::string, std::size_t> &counts) const {
    WriteValue(stream, static_cast<std::uint64_t>(counts.size()));
    for (const auto &[label, count] : counts) {
        WriteString(stream, label);
        WriteValue(stream, static_cast<std::uint64_t>(count));
    }
}

// ─────────────────────────────────────
std::map<std::string, std::size_t> Trainer::ReadCounts(std::istream &stream) const {
    const std::uint64_t size = ReadValue<std::uint64_t>(stream);
    if (size > 100000) {
        throw std::runtime_error("OpenScofo dataset cache contains too many class counts");
    }
    std::map<std::string, std::size_t> counts;
    for (std::uint64_t i = 0; i < size; ++i) {
        const std::string label = ReadString(stream);
        counts[label] = static_cast<std::size_t>(ReadValue<std::uint64_t>(stream));
    }
    return counts;
}

// ─────────────────────────────────────
void Trainer::Report(const std::string &status, double progress) const {
    if (m_Progress) {
        m_Progress(status, std::clamp(progress, 0.0, 1.0));
    }
}

// ─────────────────────────────────────
void Trainer::Print(const std::string &message) const {
    if (m_Log) {
        m_Log(message);
    }
}

// ─────────────────────────────────────
void Trainer::CheckCancelled() const {
    if (m_Cancelled.load()) {
        throw std::runtime_error("Training cancelled");
    }
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
fs::path Trainer::CacheFile() const {
    if (m_Config.CachePath.empty()) {
        return BasePath() / "dataset.openscofo-cache";
    }
    return ResolvePath(BasePath(), m_Config.CachePath);
}

// ─────────────────────────────────────
std::vector<int> Trainer::SampleRates() const {
    std::vector<int> rates = m_Config.TrainingSampleRates;
    if (rates.empty()) {
        rates.push_back(m_Config.SampleRate);
    }
    std::vector<int> unique;
    for (int rate : rates) {
        if (rate <= 0) {
            throw std::invalid_argument("Training sample rates must be positive");
        }
        if (std::find(unique.begin(), unique.end(), rate) == unique.end()) {
            unique.push_back(rate);
        }
    }
    return unique;
}

// ─────────────────────────────────────
std::vector<Descriptors> Trainer::GetDescriptors(std::vector<std::string> &canonical) const {
    if (m_Config.Descriptors.empty()) {
        throw std::invalid_argument("At least one OpenScofo descriptor is required");
    }
    ::OpenScofo::OpenScofo descriptorConverter(static_cast<float>(m_Config.SampleRate),
                                               static_cast<float>(m_Config.FFTSize),
                                               static_cast<float>(m_Config.HopSize));
    std::vector<::OpenScofo::Descriptors> descriptors;
    canonical.clear();
    for (const std::string &name : m_Config.Descriptors) {
        const std::string normalized = Lowercase(name);
        const ::OpenScofo::Descriptors descriptor = descriptorConverter.GetDescriptorsEnum(normalized.c_str());
        if (descriptor == ::OpenScofo::INVALID || descriptor == ::OpenScofo::ONNX) {
            throw std::invalid_argument("Unknown or unsupported training descriptor: " + name);
        }
        descriptors.push_back(descriptor);
        canonical.emplace_back(descriptorConverter.GetDescriptionId(descriptor));
    }
    return descriptors;
}

// ─────────────────────────────────────
void Trainer::Validate() const {
    if (m_Config.SampleRate <= 0 || m_Config.FFTSize < 512 || m_Config.HopSize <= 0 ||
        m_Config.HopSize > m_Config.FFTSize) {
        throw std::invalid_argument("Invalid sample rate, FFT size, or hop size");
    }
    if (m_Config.TestFraction < 0.0 || m_Config.TestFraction >= 1.0) {
        throw std::invalid_argument("Test fraction must be in the range [0, 1)");
    }
    if (m_Config.MaxPitchConfidence < 0.0 || m_Config.MaxPitchConfidence > 1.0) {
        throw std::invalid_argument("Maximum pitch confidence must be between 0 and 1");
    }
    if (m_Config.AugmentationFactor < 1.0) {
        throw std::invalid_argument("Augmentation factor must be at least 1");
    }
    if (m_Config.MaxAugmentationLoops < 0) {
        throw std::invalid_argument("Maximum augmentation loops cannot be negative");
    }
    if (m_Config.Iterations <= 0 || m_Config.LearningRate <= 0.0 || m_Config.Depth <= 0 ||
        m_Config.EarlyStoppingRounds < 0 || m_Config.ThreadCount == 0 || m_Config.ThreadCount < -1) {
        throw std::invalid_argument("Invalid CatBoost configuration");
    }
    (void)SampleRates();
}

// ─────────────────────────────────────
Trainer::DiscoveredData Trainer::Discover() const {
    const fs::path datasetFolder = ResolvePath(BasePath(), m_Config.DatasetFolder);
    if (datasetFolder.empty() || !fs::is_directory(datasetFolder)) {
        throw std::runtime_error("Training folder is not a directory: " + datasetFolder.string());
    }

    std::vector<fs::path> classFolders;
    for (const fs::directory_entry &entry : fs::directory_iterator(datasetFolder)) {
        if (entry.is_directory()) {
            classFolders.push_back(entry.path());
        }
    }
    std::sort(classFolders.begin(), classFolders.end());
    if (classFolders.size() < 2) {
        throw std::runtime_error("The training folder must contain at least two class subdirectories");
    }

    DiscoveredData discovered;
    std::mt19937 random(static_cast<std::mt19937::result_type>(m_Config.RandomSeed));
    for (const fs::path &classFolder : classFolders) {
        const std::size_t label = discovered.Labels.size();
        discovered.Labels.push_back(classFolder.filename().string());
        std::vector<Recording> recordings;
        for (const fs::directory_entry &entry : fs::directory_iterator(classFolder)) {
            if (!entry.is_regular_file() || !IsAudioFile(entry.path())) {
                continue;
            }
            const fs::path absolute = fs::absolute(entry.path()).lexically_normal();
            recordings.push_back({absolute, label, HashText(absolute.generic_string())});
            discovered.AudioFiles.push_back(absolute);
        }
        std::sort(recordings.begin(), recordings.end(),
                  [](const Recording &left, const Recording &right) { return left.Path < right.Path; });
        if (recordings.empty()) {
            throw std::runtime_error("Class folder contains no WAV, AIFF, or AIF files: " + classFolder.string());
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
    std::sort(discovered.AudioFiles.begin(), discovered.AudioFiles.end());

    for (const fs::path &folderValue : m_Config.ImpulseResponseFolders) {
        const fs::path folder = ResolvePath(BasePath(), folderValue);
        if (!fs::is_directory(folder)) {
            throw std::runtime_error("Impulse-response folder is not a directory: " + folder.string());
        }
        for (const fs::directory_entry &entry : fs::directory_iterator(folder)) {
            if (entry.is_regular_file() && IsAudioFile(entry.path())) {
                discovered.ImpulseResponses.push_back(fs::absolute(entry.path()).lexically_normal());
            }
        }
    }
    std::sort(discovered.ImpulseResponses.begin(), discovered.ImpulseResponses.end());
    return discovered;
}

// ─────────────────────────────────────
std::string Trainer::Signature(const DiscoveredData &discovered, const std::vector<std::string> &descriptors) const {
    std::ostringstream stream;
    stream.precision(17);
    stream << "cache=" << m_CacheVersion << '\n';
    stream << "openscofo_version=" << OPENSCOFO_VERSION << '\n';
    stream << "sample_rate=" << m_Config.SampleRate << '\n';
    stream << "fft_size=" << m_Config.FFTSize << '\n';
    stream << "hop_size=" << m_Config.HopSize << '\n';
    stream << "min_db=" << m_Config.MinDB << '\n';
    stream << "max_pitch_confidence=" << m_Config.MaxPitchConfidence << '\n';
    stream << "test_fraction=" << m_Config.TestFraction << '\n';
    stream << "augmentation_factor=" << m_Config.AugmentationFactor << '\n';
    stream << "max_augmentation_loops=" << m_Config.MaxAugmentationLoops << '\n';
    stream << "seed=" << m_Config.RandomSeed << '\n';
    stream << "dataset=" << ResolvePath(BasePath(), m_Config.DatasetFolder).generic_string() << '\n';
    for (int rate : SampleRates()) {
        stream << "training_rate=" << rate << '\n';
    }
    for (const std::string &descriptor : descriptors) {
        stream << "descriptor=" << descriptor << '\n';
    }
    for (const fs::path &file : discovered.AudioFiles) {
        stream << "audio=" << file.generic_string() << ':' << fs::file_size(file) << ':'
               << fs::last_write_time(file).time_since_epoch().count() << '\n';
    }
    for (const fs::path &file : discovered.ImpulseResponses) {
        stream << "ir=" << file.generic_string() << ':' << fs::file_size(file) << ':'
               << fs::last_write_time(file).time_since_epoch().count() << '\n';
    }
    return stream.str();
}

// ─────────────────────────────────────
bool Trainer::LoadCache(const std::string &signature) {
    const fs::path path = CacheFile();
    std::ifstream stream(path, std::ios::binary);
    if (!stream) {
        return false;
    }
    try {
        std::array<char, 8> magic{};
        ReadBytes(stream, magic.data(), magic.size());
        if (magic != m_CacheMagic) {
            return false;
        }
        if (ReadValue<std::uint32_t>(stream) != m_CacheVersion || ReadString(stream) != signature) {
            return false;
        }

        Dataset loaded;
        loaded.Labels = ReadStrings(stream);
        loaded.Descriptors = ReadStrings(stream);
        loaded.Train = ReadPartition(stream);
        loaded.Validation = ReadPartition(stream);
        loaded.Test = ReadPartition(stream);
        loaded.OriginalTrainCounts = ReadCounts(stream);
        loaded.AugmentedTrainCounts = ReadCounts(stream);
        loaded.ValidationCounts = ReadCounts(stream);
        loaded.TestCounts = ReadCounts(stream);
        loaded.LoadedFromCache = true;
        m_Data = std::move(loaded);
        Print("Loaded native dataset cache: " + path.string());
        return true;
    } catch (const std::exception &error) {
        Print(std::string("Ignoring invalid dataset cache: ") + error.what());
        return false;
    }
}

// ─────────────────────────────────────
void Trainer::SaveCache(const std::string &signature) const {
    const fs::path path = CacheFile();
    if (path.has_parent_path()) {
        fs::create_directories(path.parent_path());
    }
    fs::path temporary = path;
    temporary += ".tmp";
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) {
        throw std::runtime_error("Could not create dataset cache: " + temporary.string());
    }
    WriteBytes(stream, m_CacheMagic.data(), m_CacheMagic.size());
    WriteValue(stream, m_CacheVersion);
    WriteString(stream, signature);
    WriteStrings(stream, m_Data.Labels);
    WriteStrings(stream, m_Data.Descriptors);
    WritePartition(stream, m_Data.Train);
    WritePartition(stream, m_Data.Validation);
    WritePartition(stream, m_Data.Test);
    WriteCounts(stream, m_Data.OriginalTrainCounts);
    WriteCounts(stream, m_Data.AugmentedTrainCounts);
    WriteCounts(stream, m_Data.ValidationCounts);
    WriteCounts(stream, m_Data.TestCounts);
    stream.close();
    if (!stream) {
        throw std::runtime_error("Could not finish dataset cache: " + temporary.string());
    }
    std::error_code error;
    fs::remove(path, error);
    error.clear();
    fs::rename(temporary, path, error);
    if (error) {
        throw std::runtime_error("Could not install dataset cache: " + error.message());
    }
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
        throw std::runtime_error("The selected descriptors produced an empty feature vector");
    }
    return values;
}

// ─────────────────────────────────────
std::size_t Trainer::ProcessSignal(const std::vector<float> &signal, int sampleRate, const Recording &recording,
                                   Partition &partition, const std::vector<::OpenScofo::Descriptors> &descriptors,
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
        CheckCancelled();
        std::copy_n(signal.begin() + static_cast<std::ptrdiff_t>(start), m_Config.FFTSize, window.begin());
        extractor.ProcessBlock(window.data(), window.size());
        Description description = extractor.GetDescription();
        if (description.dB < m_Config.MinDB) {
            continue;
        }
        if (pitchFilter && description.PitchConfidence >= m_Config.MaxPitchConfidence) {
            continue;
        }
        AppendSample(partition, Flatten(extractor, description, descriptors), recording.Label, recording.Id, augmented);
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
void Trainer::ExtractRecordings(const std::vector<Recording> &recordings, Partition &partition,
                                const std::vector<::OpenScofo::Descriptors> &descriptors, bool pitchFilter,
                                std::map<std::string, std::size_t> &counts, std::size_t &completed,
                                std::size_t total) const {
    const std::vector<int> rates = SampleRates();
    for (const Recording &recording : recordings) {
        std::size_t recordingSamples = 0;
        for (int sampleRate : rates) {
            CheckCancelled();
            const std::vector<float> signal = LoadMono(recording.Path, sampleRate);
            recordingSamples +=
                ProcessSignal(signal, sampleRate, recording, partition, descriptors, pitchFilter, false);
            ++completed;
            Report("analyzing", 0.7 * static_cast<double>(completed) / static_cast<double>(total));
        }
        counts[m_Data.Labels[recording.Label]] += recordingSamples;
    }
}

// ─────────────────────────────────────
void Trainer::AssertNoLeakage() const {
    std::set<std::uint64_t> train;
    std::set<std::uint64_t> validation;
    std::set<std::uint64_t> test;
    for (const SampleOrigin &origin : m_Data.Train.Origins) {
        train.insert(origin.Recording);
    }
    for (const SampleOrigin &origin : m_Data.Validation.Origins) {
        if (origin.Augmented) {
            throw std::runtime_error("Augmented data entered the validation partition");
        }
        validation.insert(origin.Recording);
    }
    for (const SampleOrigin &origin : m_Data.Test.Origins) {
        if (origin.Augmented) {
            throw std::runtime_error("Augmented data entered the test partition");
        }
        test.insert(origin.Recording);
    }
    for (std::uint64_t id : train) {
        if (validation.contains(id) || test.contains(id)) {
            throw std::runtime_error("A recording occurs in both training and evaluation partitions");
        }
    }
    for (std::uint64_t id : validation) {
        if (test.contains(id)) {
            throw std::runtime_error("A recording occurs in both validation and test partitions");
        }
    }
}

// ─────────────────────────────────────
void Trainer::Analyze() {
    m_Cancelled.store(false);
    Validate();
    Report("analyzing", 0.0);
    std::vector<std::string> canonicalDescriptors;
    const std::vector<::OpenScofo::Descriptors> descriptors = GetDescriptors(canonicalDescriptors);
    const DiscoveredData discovered = Discover();
    const std::string signature = Signature(discovered, canonicalDescriptors);
    if (LoadCache(signature)) {
        ClearModel();
        m_Results = {};
        Report("analyzing", 1.0);
        return;
    }

    m_Data.Clear();
    m_Data.Labels = discovered.Labels;
    m_Data.Descriptors = canonicalDescriptors;
    ClearModel();
    m_Results = {};

    const std::size_t total = std::max<std::size_t>(
        1, (discovered.Train.size() + discovered.Validation.size() + discovered.Test.size()) * SampleRates().size());
    std::size_t completed = 0;
    ExtractRecordings(discovered.Test, m_Data.Test, descriptors, false, m_Data.TestCounts, completed, total);
    ExtractRecordings(discovered.Validation, m_Data.Validation, descriptors, true, m_Data.ValidationCounts, completed,
                      total);
    ExtractRecordings(discovered.Train, m_Data.Train, descriptors, true, m_Data.OriginalTrainCounts, completed, total);

    std::size_t largestClass = 0;
    for (const std::string &label : m_Data.Labels) {
        largestClass = std::max(largestClass, m_Data.OriginalTrainCounts[label]);
    }
    const std::size_t target =
        static_cast<std::size_t>(std::ceil(static_cast<double>(largestClass) * m_Config.AugmentationFactor));
    std::mt19937 random(static_cast<std::mt19937::result_type>(m_Config.RandomSeed));
    const std::vector<int> rates = SampleRates();
    for (std::size_t label = 0; label < m_Data.Labels.size(); ++label) {
        std::vector<Recording> sources;
        for (const Recording &recording : discovered.Train) {
            if (recording.Label == label) {
                sources.push_back(recording);
            }
        }
        const std::size_t original = m_Data.OriginalTrainCounts[m_Data.Labels[label]];
        const std::size_t maximumAttempts = static_cast<std::size_t>(m_Config.MaxAugmentationLoops) * sources.size();
        std::size_t attempts = 0;
        while (original + m_Data.AugmentedTrainCounts[m_Data.Labels[label]] < target && attempts < maximumAttempts &&
               !sources.empty()) {
            CheckCancelled();
            const Recording &recording = sources[attempts % sources.size()];
            const int sampleRate = rates[attempts % rates.size()];
            const std::vector<float> source = LoadMono(recording.Path, sampleRate);
            const std::vector<float> augmented = Augment(source, sampleRate, discovered.ImpulseResponses, random);
            const std::size_t remaining = target - original - m_Data.AugmentedTrainCounts[m_Data.Labels[label]];
            m_Data.AugmentedTrainCounts[m_Data.Labels[label]] +=
                ProcessSignal(augmented, sampleRate, recording, m_Data.Train, descriptors, true, true, remaining);
            ++attempts;
            const double classProgress =
                maximumAttempts == 0 ? 1.0 : static_cast<double>(attempts) / static_cast<double>(maximumAttempts);
            Report("analyzing", 0.7 + 0.25 * (static_cast<double>(label) + classProgress) /
                                          static_cast<double>(m_Data.Labels.size()));
        }
        if (original + m_Data.AugmentedTrainCounts[m_Data.Labels[label]] < target) {
            Print("Class " + m_Data.Labels[label] + " did not reach its balanced augmentation target");
        }
    }

    for (const std::string &label : m_Data.Labels) {
        if (m_Data.OriginalTrainCounts[label] == 0) {
            throw std::runtime_error("No training windows survived filtering for class: " + label);
        }
    }
    AssertNoLeakage();
    CheckCancelled();
    SaveCache(signature);
    Report("analyzing", 1.0);
}

// ─────────────────────────────────────
Trainer::Trainer() : Trainer(Configuration{}) {
}

// ─────────────────────────────────────
Trainer::Trainer(const Configuration &configuration) : m_Config(configuration) {
    if (m_Config.BasePath.empty()) {
        m_Config.BasePath = fs::current_path();
    }
}

// ─────────────────────────────────────
Trainer::~Trainer() {
    ClearModel();
}

// ─────────────────────────────────────
Trainer::Trainer(Trainer &&other) noexcept
    : m_Config(std::move(other.m_Config)), m_Data(std::move(other.m_Data)), m_Results(std::move(other.m_Results)),
      m_Progress(std::move(other.m_Progress)), m_Log(std::move(other.m_Log)), m_Cancelled(other.m_Cancelled.load()),
      m_Model(std::exchange(other.m_Model, nullptr)) {
}

// ─────────────────────────────────────
Trainer &Trainer::operator=(Trainer &&other) noexcept {
    if (this != &other) {
        ClearModel();
        m_Config = std::move(other.m_Config);
        m_Data = std::move(other.m_Data);
        m_Results = std::move(other.m_Results);
        m_Progress = std::move(other.m_Progress);
        m_Log = std::move(other.m_Log);
        m_Cancelled.store(other.m_Cancelled.load());
        m_Model = std::exchange(other.m_Model, nullptr);
    }
    return *this;
}

// ─────────────────────────────────────
void Trainer::SetConfiguration(const Configuration &configuration) {
    m_Config = configuration;
    if (m_Config.BasePath.empty()) {
        m_Config.BasePath = fs::current_path();
    }
    m_Data.Clear();
    m_Results = {};
    ClearModel();
}

// ─────────────────────────────────────
const Trainer::Configuration &Trainer::GetConfiguration() const {
    return m_Config;
}

// ─────────────────────────────────────
void Trainer::SetProgressCallback(ProgressCallback callback) {
    m_Progress = std::move(callback);
}

// ─────────────────────────────────────
void Trainer::SetLogCallback(LogCallback callback) {
    m_Log = std::move(callback);
}

// ─────────────────────────────────────
void Trainer::Train() {
    m_Cancelled.store(false);
    if (!HasDataset()) {
        Analyze();
    }
    Report("training", 0.0);
    TrainModel();
    m_Results = EvaluateModel();
    if (m_Results.TestedSamples == 0) {
        Print("No independent test samples are available; test evaluation was skipped");
    }
}

// ─────────────────────────────────────
void Trainer::ExportModel(const fs::path &path) {
    CheckCancelled();
    Report("exporting", 0.0);
    ExportTrainedModel(ResolvePath(BasePath(), path));
    Report("exporting", 1.0);
}

// ─────────────────────────────────────
void Trainer::ClearData() {
    m_Cancelled.store(false);
    const fs::path cache = CacheFile();
    std::error_code error;
    const bool removed = fs::remove(cache, error);
    if (error) {
        throw std::runtime_error("Could not remove dataset cache: " + error.message());
    }
    if (removed) {
        Print("Deleted native dataset cache: " + cache.string());
    }
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
const Trainer::Dataset &Trainer::GetDataset() const {
    return m_Data;
}

// ─────────────────────────────────────
const Trainer::Evaluation &Trainer::GetEvaluation() const {
    return m_Results;
}

} // namespace OpenScofo
