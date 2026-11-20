/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "Trainer.hpp"

#include "Audio.hpp"
#include "CatBoostBackend.hpp"

#include <OpenScofo.hpp>

#include <algorithm>
#include <atomic>
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

namespace OpenScofo::Training {

namespace fs = std::filesystem;

static constexpr std::uint32_t CacheVersion = 1;
static constexpr char CacheMagic[8] = {'O', 'S', 'C', 'T', 'R', 'N', '1', '\0'};

std::size_t Partition::SampleCount() const {
    return Labels.size();
}

bool Partition::Empty() const {
    return Labels.empty();
}

void Partition::Clear() {
    Features.clear();
    Labels.clear();
    Origins.clear();
    FeatureCount = 0;
}

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
    LoadedFromCache = false;
}

enum class PartitionKind { Train, Validation, Test };

struct Recording {
    fs::path Path;
    std::size_t Label = 0;
    std::uint64_t Id = 0;
};

struct DiscoveredData {
    std::vector<std::string> Labels;
    std::vector<Recording> Train;
    std::vector<Recording> Validation;
    std::vector<Recording> Test;
    std::vector<fs::path> AudioFiles;
    std::vector<fs::path> ImpulseResponses;
};

static std::string Lowercase(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
    return text;
}

static bool IsAudioFile(const fs::path &path) {
    const std::string extension = Lowercase(path.extension().string());
    return extension == ".wav" || extension == ".aif" || extension == ".aiff";
}

static std::uint64_t HashText(const std::string &text) {
    std::uint64_t value = 1469598103934665603ULL;
    for (unsigned char character : text) {
        value ^= character;
        value *= 1099511628211ULL;
    }
    return value;
}

static fs::path ResolvePath(const fs::path &base, const fs::path &path) {
    if (path.empty()) {
        return {};
    }
    if (path.is_absolute()) {
        return path.lexically_normal();
    }
    return (base / path).lexically_normal();
}

static ::OpenScofo::Descriptors DescriptorFromName(const std::string &source) {
    const std::string name = Lowercase(source);
    if (name == "onset")
        return ::OpenScofo::ODSONSET;
    if (name == "loudness")
        return ::OpenScofo::LOUDNESS;
    if (name == "db")
        return ::OpenScofo::DB;
    if (name == "maxamp" || name == "max_amp")
        return ::OpenScofo::MAXAMP;
    if (name == "rms")
        return ::OpenScofo::RMS;
    if (name == "stddev")
        return ::OpenScofo::STDDEV;
    if (name == "magnitude")
        return ::OpenScofo::MAGNITUDE;
    if (name == "power" || name == "powerarray")
        return ::OpenScofo::POWERARRAY;
    if (name == "silence")
        return ::OpenScofo::SILENCEPROB;
    if (name == "mfcc")
        return ::OpenScofo::MFCC;
    if (name == "chroma")
        return ::OpenScofo::CHROMA;
    if (name == "logmel")
        return ::OpenScofo::LOGMEL;
    if (name == "zcr")
        return ::OpenScofo::ZCR;
    if (name == "hfr")
        return ::OpenScofo::HFR;
    if (name == "centroid")
        return ::OpenScofo::CENTROID;
    if (name == "spread")
        return ::OpenScofo::SPREADHZ;
    if (name == "spread_variance")
        return ::OpenScofo::SPREADVARIANCE;
    if (name == "crest")
        return ::OpenScofo::CREST;
    if (name == "flatness")
        return ::OpenScofo::FLATNESS;
    if (name == "entropy")
        return ::OpenScofo::ENTROPY;
    if (name == "rolloff")
        return ::OpenScofo::ROLLOFF;
    if (name == "centroid_velocity")
        return ::OpenScofo::CENTROIDVEL;
    if (name == "flux")
        return ::OpenScofo::FLUX;
    if (name == "skewness")
        return ::OpenScofo::SKEWNESS;
    if (name == "slope")
        return ::OpenScofo::SLOPE;
    if (name == "kurtosis")
        return ::OpenScofo::KURTOSIS;
    if (name == "irregularity")
        return ::OpenScofo::IRREGULARITY;
    if (name == "harmonicity")
        return ::OpenScofo::HARMONICITY;
    if (name == "yin")
        return ::OpenScofo::YIN;
    if (name == "yin_confidence" || name == "pitch_confidence")
        return ::OpenScofo::YINCONFIDENCE;
    if (name == "ext")
        return ::OpenScofo::EXTENDEDTECHNIQUE;
    return ::OpenScofo::INVALID;
}

static const char *CanonicalDescriptorName(::OpenScofo::Descriptors descriptor) {
    switch (descriptor) {
    case ::OpenScofo::ODSONSET:
        return "onset";
    case ::OpenScofo::LOUDNESS:
        return "loudness";
    case ::OpenScofo::DB:
        return "db";
    case ::OpenScofo::MAXAMP:
        return "maxamp";
    case ::OpenScofo::RMS:
        return "rms";
    case ::OpenScofo::STDDEV:
        return "stddev";
    case ::OpenScofo::MAGNITUDE:
        return "magnitude";
    case ::OpenScofo::POWERARRAY:
        return "power";
    case ::OpenScofo::SILENCEPROB:
        return "silence";
    case ::OpenScofo::MFCC:
        return "mfcc";
    case ::OpenScofo::CHROMA:
        return "chroma";
    case ::OpenScofo::LOGMEL:
        return "logmel";
    case ::OpenScofo::ZCR:
        return "zcr";
    case ::OpenScofo::HFR:
        return "hfr";
    case ::OpenScofo::CENTROID:
        return "centroid";
    case ::OpenScofo::SPREADHZ:
        return "spread";
    case ::OpenScofo::SPREADVARIANCE:
        return "spread_variance";
    case ::OpenScofo::CREST:
        return "crest";
    case ::OpenScofo::FLATNESS:
        return "flatness";
    case ::OpenScofo::ENTROPY:
        return "entropy";
    case ::OpenScofo::ROLLOFF:
        return "rolloff";
    case ::OpenScofo::CENTROIDVEL:
        return "centroid_velocity";
    case ::OpenScofo::FLUX:
        return "flux";
    case ::OpenScofo::SKEWNESS:
        return "skewness";
    case ::OpenScofo::SLOPE:
        return "slope";
    case ::OpenScofo::KURTOSIS:
        return "kurtosis";
    case ::OpenScofo::IRREGULARITY:
        return "irregularity";
    case ::OpenScofo::HARMONICITY:
        return "harmonicity";
    case ::OpenScofo::YIN:
        return "yin";
    case ::OpenScofo::YINCONFIDENCE:
        return "yin_confidence";
    case ::OpenScofo::EXTENDEDTECHNIQUE:
        return "ext";
    default:
        return "unknown";
    }
}

static bool IsVectorDescriptor(::OpenScofo::Descriptors descriptor) {
    return descriptor == ::OpenScofo::MFCC || descriptor == ::OpenScofo::CHROMA || descriptor == ::OpenScofo::LOGMEL ||
           descriptor == ::OpenScofo::MAGNITUDE || descriptor == ::OpenScofo::POWERARRAY;
}

static void AppendSample(Partition &partition, const std::vector<float> &features, std::size_t label,
                         std::uint64_t recording, bool augmented) {
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

static void WriteBytes(std::ostream &stream, const void *data, std::size_t size) {
    if (size == 0) {
        return;
    }
    stream.write(static_cast<const char *>(data), static_cast<std::streamsize>(size));
    if (!stream) {
        throw std::runtime_error("Could not write the OpenScofo dataset cache");
    }
}

static void ReadBytes(std::istream &stream, void *data, std::size_t size) {
    if (size == 0) {
        return;
    }
    stream.read(static_cast<char *>(data), static_cast<std::streamsize>(size));
    if (!stream) {
        throw std::runtime_error("OpenScofo dataset cache is truncated");
    }
}

template <typename T> static void WriteValue(std::ostream &stream, const T &value) {
    WriteBytes(stream, &value, sizeof(value));
}

template <typename T> static T ReadValue(std::istream &stream) {
    T value{};
    ReadBytes(stream, &value, sizeof(value));
    return value;
}

static void WriteString(std::ostream &stream, const std::string &value) {
    const std::uint64_t size = value.size();
    WriteValue(stream, size);
    WriteBytes(stream, value.data(), value.size());
}

static std::string ReadString(std::istream &stream) {
    const std::uint64_t size = ReadValue<std::uint64_t>(stream);
    if (size > 16ULL * 1024ULL * 1024ULL) {
        throw std::runtime_error("OpenScofo dataset cache contains an invalid string size");
    }
    std::string value(static_cast<std::size_t>(size), '\0');
    ReadBytes(stream, value.data(), value.size());
    return value;
}

static void WriteStrings(std::ostream &stream, const std::vector<std::string> &values) {
    WriteValue(stream, static_cast<std::uint64_t>(values.size()));
    for (const std::string &value : values) {
        WriteString(stream, value);
    }
}

static std::vector<std::string> ReadStrings(std::istream &stream) {
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

static void WritePartition(std::ostream &stream, const Partition &partition) {
    WriteValue(stream, static_cast<std::uint64_t>(partition.FeatureCount));
    WriteValue(stream, static_cast<std::uint64_t>(partition.SampleCount()));
    WriteBytes(stream, partition.Features.data(), partition.Features.size() * sizeof(float));
    WriteBytes(stream, partition.Labels.data(), partition.Labels.size() * sizeof(std::int32_t));
    for (const SampleOrigin &origin : partition.Origins) {
        WriteValue(stream, origin.Recording);
        WriteValue(stream, static_cast<std::uint8_t>(origin.Augmented ? 1 : 0));
    }
}

static Partition ReadPartition(std::istream &stream) {
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

static void WriteCounts(std::ostream &stream, const std::map<std::string, std::size_t> &counts) {
    WriteValue(stream, static_cast<std::uint64_t>(counts.size()));
    for (const auto &[label, count] : counts) {
        WriteString(stream, label);
        WriteValue(stream, static_cast<std::uint64_t>(count));
    }
}

static std::map<std::string, std::size_t> ReadCounts(std::istream &stream) {
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

class Trainer::Impl {
  public:
    Configuration Config;
    Dataset Data;
    Evaluation Results;
    ProgressCallback Progress;
    LogCallback Log;
    CatBoostBackend Backend;
    std::atomic<bool> Cancelled{false};

    explicit Impl(Configuration configuration) : Config(std::move(configuration)) {
        if (Config.BasePath.empty()) {
            Config.BasePath = fs::current_path();
        }
    }

    void Report(const std::string &status, double progress) const {
        if (Progress) {
            Progress(status, std::clamp(progress, 0.0, 1.0));
        }
    }

    void Print(const std::string &message) const {
        if (Log) {
            Log(message);
        }
    }

    void CheckCancelled() const {
        if (Cancelled.load()) {
            throw std::runtime_error("Training cancelled");
        }
    }

    fs::path BasePath() const {
        fs::path base = Config.BasePath.empty() ? fs::current_path() : Config.BasePath;
        if (!base.is_absolute()) {
            base = fs::absolute(base);
        }
        return base.lexically_normal();
    }

    fs::path CacheFile() const {
        if (Config.CachePath.empty()) {
            return BasePath() / "dataset.openscofo-cache";
        }
        return ResolvePath(BasePath(), Config.CachePath);
    }

    std::vector<int> SampleRates() const {
        std::vector<int> rates = Config.TrainingSampleRates;
        if (rates.empty()) {
            rates.push_back(Config.SampleRate);
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

    std::vector<::OpenScofo::Descriptors> Descriptors(std::vector<std::string> &canonical) const {
        if (Config.Descriptors.empty()) {
            throw std::invalid_argument("At least one OpenScofo descriptor is required");
        }
        std::vector<::OpenScofo::Descriptors> descriptors;
        canonical.clear();
        for (const std::string &name : Config.Descriptors) {
            const ::OpenScofo::Descriptors descriptor = DescriptorFromName(name);
            if (descriptor == ::OpenScofo::INVALID || descriptor == ::OpenScofo::ONNX) {
                throw std::invalid_argument("Unknown or unsupported training descriptor: " + name);
            }
            descriptors.push_back(descriptor);
            canonical.emplace_back(CanonicalDescriptorName(descriptor));
        }
        return descriptors;
    }

    void Validate() const {
        if (Config.SampleRate <= 0 || Config.FFTSize < 512 || Config.HopSize <= 0 || Config.HopSize > Config.FFTSize) {
            throw std::invalid_argument("Invalid sample rate, FFT size, or hop size");
        }
        if (Config.TestFraction < 0.0 || Config.TestFraction >= 1.0) {
            throw std::invalid_argument("Test fraction must be in the range [0, 1)");
        }
        if (Config.MaxPitchConfidence < 0.0 || Config.MaxPitchConfidence > 1.0) {
            throw std::invalid_argument("Maximum pitch confidence must be between 0 and 1");
        }
        if (Config.AugmentationFactor < 1.0) {
            throw std::invalid_argument("Augmentation factor must be at least 1");
        }
        if (Config.MaxAugmentationLoops < 0) {
            throw std::invalid_argument("Maximum augmentation loops cannot be negative");
        }
        if (Config.Iterations <= 0 || Config.LearningRate <= 0.0 || Config.Depth <= 0 ||
            Config.EarlyStoppingRounds < 0 || Config.ThreadCount == 0 || Config.ThreadCount < -1) {
            throw std::invalid_argument("Invalid CatBoost configuration");
        }
        (void)SampleRates();
        std::vector<std::string> canonical;
        (void)Descriptors(canonical);
    }

    DiscoveredData Discover() const {
        const fs::path datasetFolder = ResolvePath(BasePath(), Config.DatasetFolder);
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
        std::mt19937 random(static_cast<std::mt19937::result_type>(Config.RandomSeed));
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
            if (Config.TestFraction > 0.0 && recordings.size() >= 2) {
                testCount = std::max<std::size_t>(1, static_cast<std::size_t>(recordings.size() * Config.TestFraction));
                testCount = std::min(testCount, recordings.size() - 1);
            }
            const std::size_t afterTest = recordings.size() - testCount;
            if (Config.TestFraction > 0.0 && afterTest >= 2) {
                validationCount =
                    std::max<std::size_t>(1, static_cast<std::size_t>(recordings.size() * Config.TestFraction));
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

        for (const fs::path &folderValue : Config.ImpulseResponseFolders) {
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

    std::string Signature(const DiscoveredData &discovered, const std::vector<std::string> &descriptors) const {
        std::ostringstream stream;
        stream.precision(17);
        stream << "cache=" << CacheVersion << '\n';
        stream << "openscofo_version=" << OPENSCOFO_VERSION << '\n';
        stream << "sample_rate=" << Config.SampleRate << '\n';
        stream << "fft_size=" << Config.FFTSize << '\n';
        stream << "hop_size=" << Config.HopSize << '\n';
        stream << "min_db=" << Config.MinDB << '\n';
        stream << "max_pitch_confidence=" << Config.MaxPitchConfidence << '\n';
        stream << "test_fraction=" << Config.TestFraction << '\n';
        stream << "augmentation_factor=" << Config.AugmentationFactor << '\n';
        stream << "max_augmentation_loops=" << Config.MaxAugmentationLoops << '\n';
        stream << "seed=" << Config.RandomSeed << '\n';
        stream << "dataset=" << ResolvePath(BasePath(), Config.DatasetFolder).generic_string() << '\n';
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

    bool LoadCache(const std::string &signature) {
        const fs::path path = CacheFile();
        std::ifstream stream(path, std::ios::binary);
        if (!stream) {
            return false;
        }
        try {
            char magic[sizeof(CacheMagic)]{};
            ReadBytes(stream, magic, sizeof(magic));
            if (!std::equal(std::begin(magic), std::end(magic), std::begin(CacheMagic))) {
                return false;
            }
            if (ReadValue<std::uint32_t>(stream) != CacheVersion || ReadString(stream) != signature) {
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
            Data = std::move(loaded);
            Print("Loaded native dataset cache: " + path.string());
            return true;
        } catch (const std::exception &error) {
            Print(std::string("Ignoring invalid dataset cache: ") + error.what());
            return false;
        }
    }

    void SaveCache(const std::string &signature) const {
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
        WriteBytes(stream, CacheMagic, sizeof(CacheMagic));
        WriteValue(stream, CacheVersion);
        WriteString(stream, signature);
        WriteStrings(stream, Data.Labels);
        WriteStrings(stream, Data.Descriptors);
        WritePartition(stream, Data.Train);
        WritePartition(stream, Data.Validation);
        WritePartition(stream, Data.Test);
        WriteCounts(stream, Data.OriginalTrainCounts);
        WriteCounts(stream, Data.AugmentedTrainCounts);
        WriteCounts(stream, Data.ValidationCounts);
        WriteCounts(stream, Data.TestCounts);
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

    std::vector<float> Flatten(::OpenScofo::OpenScofo &extractor, ::OpenScofo::Description &description,
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

    std::size_t ProcessSignal(const std::vector<float> &signal, int sampleRate, const Recording &recording,
                              Partition &partition, const std::vector<::OpenScofo::Descriptors> &descriptors,
                              bool pitchFilter, bool augmented, std::size_t maximumSamples = 0) const {
        if (signal.size() < static_cast<std::size_t>(Config.FFTSize)) {
            return 0;
        }
        std::vector<::OpenScofo::Descriptors> requested = descriptors;
        if (std::find(requested.begin(), requested.end(), ::OpenScofo::YINCONFIDENCE) == requested.end()) {
            requested.push_back(::OpenScofo::YINCONFIDENCE);
        }
        ::OpenScofo::OpenScofo extractor(static_cast<float>(sampleRate), static_cast<float>(Config.FFTSize),
                                         static_cast<float>(Config.HopSize));
        extractor.SetRequestedDescriptors(requested);

        std::size_t appended = 0;
        std::vector<float> window(static_cast<std::size_t>(Config.FFTSize));
        for (std::size_t start = 0; start + static_cast<std::size_t>(Config.FFTSize) <= signal.size();
             start += static_cast<std::size_t>(Config.HopSize)) {
            CheckCancelled();
            std::copy_n(signal.begin() + static_cast<std::ptrdiff_t>(start), Config.FFTSize, window.begin());
            extractor.ProcessBlock(window.data(), window.size());
            ::OpenScofo::Description description = extractor.GetDescription();
            if (description.dB < Config.MinDB) {
                continue;
            }
            if (pitchFilter && description.PitchConfidence >= Config.MaxPitchConfidence) {
                continue;
            }
            AppendSample(partition, Flatten(extractor, description, descriptors), recording.Label, recording.Id,
                         augmented);
            ++appended;
            if (maximumSamples != 0 && appended >= maximumSamples) {
                break;
            }
        }
        return appended;
    }

    std::vector<float> Augment(const std::vector<float> &input, int sampleRate,
                               const std::vector<fs::path> &impulseResponses, std::mt19937 &random) const {
        const int choiceCount = impulseResponses.empty() ? 3 : 4;
        std::uniform_int_distribution<int> choice(0, choiceCount - 1);
        const int selected = choice(random);
        if (selected == 0) {
            std::uniform_real_distribution<double> rate(0.7, 1.2);
            return Audio::TimeStretch(input, rate(random), Config.FFTSize, std::max(1, Config.FFTSize / 4));
        }
        if (selected == 1) {
            std::uniform_real_distribution<double> steps(-2.0, 2.0);
            return Audio::PitchShift(input, steps(random), Config.FFTSize, std::max(1, Config.FFTSize / 4));
        }
        if (selected == 2) {
            std::uniform_real_distribution<double> deviation(0.005, 0.01);
            return Audio::AddGaussianNoise(input, deviation(random), random);
        }

        std::uniform_int_distribution<std::size_t> impulseIndex(0, impulseResponses.size() - 1);
        std::vector<float> impulse = Audio::LoadMono(impulseResponses[impulseIndex(random)], sampleRate);
        return Audio::Convolve(input, impulse);
    }

    void ExtractRecordings(const std::vector<Recording> &recordings, Partition &partition,
                           const std::vector<::OpenScofo::Descriptors> &descriptors, bool pitchFilter,
                           std::map<std::string, std::size_t> &counts, std::size_t &completed,
                           std::size_t total) const {
        const std::vector<int> rates = SampleRates();
        for (const Recording &recording : recordings) {
            std::size_t recordingSamples = 0;
            for (int sampleRate : rates) {
                CheckCancelled();
                const std::vector<float> signal = Audio::LoadMono(recording.Path, sampleRate);
                recordingSamples +=
                    ProcessSignal(signal, sampleRate, recording, partition, descriptors, pitchFilter, false);
                ++completed;
                Report("analyzing", 0.7 * static_cast<double>(completed) / static_cast<double>(total));
            }
            counts[Data.Labels[recording.Label]] += recordingSamples;
        }
    }

    void AssertNoLeakage() const {
        std::set<std::uint64_t> train;
        std::set<std::uint64_t> validation;
        std::set<std::uint64_t> test;
        for (const SampleOrigin &origin : Data.Train.Origins) {
            train.insert(origin.Recording);
        }
        for (const SampleOrigin &origin : Data.Validation.Origins) {
            if (origin.Augmented) {
                throw std::runtime_error("Augmented data entered the validation partition");
            }
            validation.insert(origin.Recording);
        }
        for (const SampleOrigin &origin : Data.Test.Origins) {
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

    void Analyze() {
        Cancelled.store(false);
        Validate();
        Report("analyzing", 0.0);
        std::vector<std::string> canonicalDescriptors;
        const std::vector<::OpenScofo::Descriptors> descriptors = Descriptors(canonicalDescriptors);
        const DiscoveredData discovered = Discover();
        const std::string signature = Signature(discovered, canonicalDescriptors);
        if (LoadCache(signature)) {
            Backend.Clear();
            Results = {};
            Report("analyzing", 1.0);
            return;
        }

        Data.Clear();
        Data.Labels = discovered.Labels;
        Data.Descriptors = canonicalDescriptors;
        Backend.Clear();
        Results = {};

        const std::size_t total =
            std::max<std::size_t>(1, (discovered.Train.size() + discovered.Validation.size() + discovered.Test.size()) *
                                         SampleRates().size());
        std::size_t completed = 0;
        ExtractRecordings(discovered.Test, Data.Test, descriptors, false, Data.TestCounts, completed, total);
        ExtractRecordings(discovered.Validation, Data.Validation, descriptors, true, Data.ValidationCounts, completed,
                          total);
        ExtractRecordings(discovered.Train, Data.Train, descriptors, true, Data.OriginalTrainCounts, completed, total);

        std::size_t largestClass = 0;
        for (const std::string &label : Data.Labels) {
            largestClass = std::max(largestClass, Data.OriginalTrainCounts[label]);
        }
        const std::size_t target =
            static_cast<std::size_t>(std::ceil(static_cast<double>(largestClass) * Config.AugmentationFactor));
        std::mt19937 random(static_cast<std::mt19937::result_type>(Config.RandomSeed));
        const std::vector<int> rates = SampleRates();
        for (std::size_t label = 0; label < Data.Labels.size(); ++label) {
            std::vector<Recording> sources;
            for (const Recording &recording : discovered.Train) {
                if (recording.Label == label) {
                    sources.push_back(recording);
                }
            }
            const std::size_t original = Data.OriginalTrainCounts[Data.Labels[label]];
            const std::size_t maximumAttempts = static_cast<std::size_t>(Config.MaxAugmentationLoops) * sources.size();
            std::size_t attempts = 0;
            while (original + Data.AugmentedTrainCounts[Data.Labels[label]] < target && attempts < maximumAttempts &&
                   !sources.empty()) {
                CheckCancelled();
                const Recording &recording = sources[attempts % sources.size()];
                const int sampleRate = rates[attempts % rates.size()];
                const std::vector<float> source = Audio::LoadMono(recording.Path, sampleRate);
                const std::vector<float> augmented = Augment(source, sampleRate, discovered.ImpulseResponses, random);
                const std::size_t remaining = target - original - Data.AugmentedTrainCounts[Data.Labels[label]];
                Data.AugmentedTrainCounts[Data.Labels[label]] +=
                    ProcessSignal(augmented, sampleRate, recording, Data.Train, descriptors, true, true, remaining);
                ++attempts;
                const double classProgress =
                    maximumAttempts == 0 ? 1.0 : static_cast<double>(attempts) / static_cast<double>(maximumAttempts);
                Report("analyzing", 0.7 + 0.25 * (static_cast<double>(label) + classProgress) /
                                              static_cast<double>(Data.Labels.size()));
            }
            if (original + Data.AugmentedTrainCounts[Data.Labels[label]] < target) {
                Print("Class " + Data.Labels[label] + " did not reach its balanced augmentation target");
            }
        }

        for (const std::string &label : Data.Labels) {
            if (Data.OriginalTrainCounts[label] == 0) {
                throw std::runtime_error("No training windows survived filtering for class: " + label);
            }
        }
        AssertNoLeakage();
        CheckCancelled();
        SaveCache(signature);
        Report("analyzing", 1.0);
    }
};

Trainer::Trainer() : m_Impl(std::make_unique<Impl>(Configuration{})) {
}

Trainer::Trainer(const Configuration &configuration) : m_Impl(std::make_unique<Impl>(configuration)) {
}

Trainer::~Trainer() = default;
Trainer::Trainer(Trainer &&) noexcept = default;
Trainer &Trainer::operator=(Trainer &&) noexcept = default;

void Trainer::SetConfiguration(const Configuration &configuration) {
    m_Impl->Config = configuration;
    if (m_Impl->Config.BasePath.empty()) {
        m_Impl->Config.BasePath = fs::current_path();
    }
    m_Impl->Data.Clear();
    m_Impl->Results = {};
    m_Impl->Backend.Clear();
}

const Configuration &Trainer::GetConfiguration() const {
    return m_Impl->Config;
}

void Trainer::SetProgressCallback(ProgressCallback callback) {
    m_Impl->Progress = std::move(callback);
}

void Trainer::SetLogCallback(LogCallback callback) {
    m_Impl->Log = std::move(callback);
}

void Trainer::Analyze() {
    m_Impl->Analyze();
}

void Trainer::Train() {
    m_Impl->Cancelled.store(false);
    if (!HasDataset()) {
        m_Impl->Analyze();
    }
    m_Impl->Report("training", 0.0);
    m_Impl->Backend.Train(m_Impl->Data, m_Impl->Config, m_Impl->Cancelled, m_Impl->Progress);
    m_Impl->Results = m_Impl->Backend.Evaluate(m_Impl->Data);
    if (m_Impl->Results.TestedSamples == 0) {
        m_Impl->Print("No independent test samples are available; test evaluation was skipped");
    }
}

void Trainer::ExportModel(const fs::path &path) {
    m_Impl->CheckCancelled();
    m_Impl->Report("exporting", 0.0);
    m_Impl->Backend.Export(ResolvePath(m_Impl->BasePath(), path), m_Impl->Data, m_Impl->Config);
    m_Impl->Report("exporting", 1.0);
}

void Trainer::ClearData() {
    m_Impl->Cancelled.store(false);
    const fs::path cache = m_Impl->CacheFile();
    std::error_code error;
    const bool removed = fs::remove(cache, error);
    if (error) {
        throw std::runtime_error("Could not remove dataset cache: " + error.message());
    }
    if (removed) {
        m_Impl->Print("Deleted native dataset cache: " + cache.string());
    }
    m_Impl->Data.Clear();
    m_Impl->Results = {};
    m_Impl->Backend.Clear();
}

void Trainer::Cancel() {
    m_Impl->Cancelled.store(true);
}

bool Trainer::HasDataset() const {
    return !m_Impl->Data.Train.Empty();
}

bool Trainer::HasModel() const {
    return m_Impl->Backend.HasModel();
}

const Dataset &Trainer::GetDataset() const {
    return m_Impl->Data;
}

const Evaluation &Trainer::GetEvaluation() const {
    return m_Impl->Results;
}

} // namespace OpenScofo::Training
