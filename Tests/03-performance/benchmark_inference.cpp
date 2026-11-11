/*
    Copyright (c) 2024-2026 Charles K. Neimog
    SPDX-License-Identifier: GPL-3.0-or-later
*/

#include <OpenScofo.hpp>
#include <AudioFile.h>
#include <nlohmann/json.hpp>

#define MINIMP3_IMPLEMENTATION
#include <minimp3_ex.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;
using Json = nlohmann::json;

constexpr size_t BlockSize = 64;
constexpr int FFTSize = 2048;
constexpr int HopSize = 512;
constexpr double AlignmentToleranceSeconds = 0.250; // Cont (2010), Section 8.2
constexpr double SpeedRegressionTolerance = 0.10;   // 10% tolerance for host scheduling noise
constexpr double OffsetEpsilonMs = 0.001;
const fs::path Data = TEST_DATA_DIR;
const fs::path RecordPath = Data / "03-performance/inference_benchmark_best.json";

struct Piece {
    const char *Name;
    fs::path JsonPath;
    fs::path ScorePath;
    fs::path AudioDirectory;
    int SampleRate = 0;
    std::vector<float> Samples;
    std::map<int, double> Expected;
};

struct Result {
    double Seconds = 0.0;
    int FinalPosition = 0;
    size_t DeadlineMisses = 0;
};

struct Accuracy {
    size_t ReferenceEvents = 0;
    size_t Matched = 0;
    size_t Missed = 0;
    size_t FalsePositives = 0;
    double PrecisionPercent = 0.0;
    double MeanAbsoluteOffsetMs = 0.0;
    double MeanSignedOffsetMs = 0.0;
    double StddevOffsetMs = 0.0;
};

struct Benchmark {
    Accuracy WorstAccuracy;
    double MeanRunMs = 0.0;
    double MedianRunMs = 0.0;
    double MinRunMs = 0.0;
    double MaxRunMs = 0.0;
    double P99BlockUs = 0.0;
    size_t DeadlineMisses = 0;
    size_t BlockCalls = 0;
    int FinalPosition = 0;
};

// ─────────────────────────────────────
static std::vector<Piece> GetPieces() {
    return {{"Miniatura1",
             Data / "04-miniaturas/Extras/miniatura1.json",
             Data / "04-miniaturas/Extras/miniatura1.scofo",
             Data / "04-miniaturas/Audios",
             0,
             {},
             {}},
            {"Miniatura2",
             Data / "04-miniaturas/Extras/miniatura2.json",
             Data / "04-miniaturas/Extras/miniatura2.scofo",
             Data / "04-miniaturas/Audios",
             0,
             {},
             {}},
            {"Miniatura3",
             Data / "04-miniaturas/Extras/miniatura3.json",
             Data / "04-miniaturas/Extras/miniatura3.scofo",
             Data / "04-miniaturas/Audios",
             0,
             {},
             {}},
            {"Canticos",
             Data / "01-benchmark/real/canticos.json",
             Data / "01-benchmark/real/canticos.scofo",
             Data / "01-benchmark/real",
             0,
             {},
             {}}};
}

// ─────────────────────────────────────
static void LoadWav(Piece &P, const fs::path &Path) {
    AudioFile<float> Audio;
    if (!Audio.load(Path.string()) || Audio.getNumChannels() == 0 || Audio.getNumSamplesPerChannel() == 0) {
        throw std::runtime_error("Could not load WAV: " + Path.string());
    }
    P.SampleRate = Audio.getSampleRate();
    P.Samples = Audio.samples[0];
}

// ─────────────────────────────────────
static void LoadMp3(Piece &P, const fs::path &Path) {
    mp3dec_t Decoder;
    mp3dec_file_info_t Info{};
    mp3dec_init(&Decoder);
    const int Error = mp3dec_load(&Decoder, Path.string().c_str(), &Info, nullptr, nullptr);
    std::unique_ptr<mp3d_sample_t, decltype(&std::free)> Buffer(Info.buffer, &std::free);
    if (Error != 0 || !Buffer || Info.channels <= 0 || Info.hz <= 0 || Info.samples == 0) {
        throw std::runtime_error("Could not load MP3: " + Path.string());
    }
    P.SampleRate = Info.hz;
    P.Samples.resize(Info.samples / Info.channels);
    for (size_t i = 0; i < P.Samples.size(); ++i) {
        P.Samples[i] = static_cast<float>(Buffer.get()[i * Info.channels]) / 32768.0f;
    }
}

// ─────────────────────────────────────
static void LoadPiece(Piece &P) {
    std::ifstream File(P.JsonPath);
    if (!File) {
        throw std::runtime_error("Could not open annotations: " + P.JsonPath.string());
    }
    Json Annotations;
    File >> Annotations;

    const fs::path AudioPath = P.AudioDirectory / Annotations.at("audio").get<std::string>();
    if (AudioPath.extension() == ".wav") {
        LoadWav(P, AudioPath);
    } else if (AudioPath.extension() == ".mp3") {
        LoadMp3(P, AudioPath);
    } else {
        throw std::runtime_error("Unsupported audio format: " + AudioPath.string());
    }
    if (P.SampleRate <= 0 || P.Samples.empty()) {
        throw std::runtime_error("Invalid audio: " + AudioPath.string());
    }

    for (const Json &Event : Annotations.at("events")) {
        const int Position = Event.at("event").get<int>();
        const double Timestamp = Event.at("timestamp_seconds").get<double>();
        if (Position <= 0 || !std::isfinite(Timestamp) || Timestamp < 0.0) {
            throw std::runtime_error("Invalid event annotation: " + P.JsonPath.string());
        }
        if (!P.Expected.emplace(Position, Timestamp).second) {
            throw std::runtime_error("Duplicate event annotation: " + P.JsonPath.string());
        }
    }
    if (P.Expected.empty()) {
        throw std::runtime_error("No annotated events: " + P.JsonPath.string());
    }
}

// ─────────────────────────────────────
static void InitializeInference(OpenScofo::OpenScofo &Follower, const Piece &P) {
    if (!Follower.LoadScore(P.ScorePath.string())) {
        throw std::runtime_error("Could not load score: " + P.ScorePath.string());
    }
    Follower.SetCurrentEvent(0);
}

// ─────────────────────────────────────
static void RecordDetection(OpenScofo::OpenScofo &Follower, const Piece &P, size_t Consumed,
                            std::map<int, double> *Detected) {
    if (Detected == nullptr) {
        return;
    }
    const int Position = Follower.GetCurrentScorePosition();
    if (Position > 0) {
        // As in inference_tests, use the first observed time for each event tag.
        Detected->emplace(Position, static_cast<double>(Consumed) / P.SampleRate);
    }
}

// ─────────────────────────────────────
static void ProcessAudio(OpenScofo::OpenScofo &Follower, const Piece &P, std::map<int, double> *Detected) {
    std::array<float, BlockSize> Block{};
    size_t Consumed = 0;
    while (Consumed < P.Samples.size()) {
        const size_t Frames = std::min(BlockSize, P.Samples.size() - Consumed);
        std::copy_n(P.Samples.data() + Consumed, Frames, Block.data());
        if (!Follower.ProcessBlock(Block.data(), Frames)) {
            throw std::runtime_error("ProcessBlock failed: " + std::string(P.Name));
        }
        Consumed += Frames;
        RecordDetection(Follower, P, Consumed, Detected);
    }
}

// ─────────────────────────────────────
static size_t ProcessAudioLatency(OpenScofo::OpenScofo &Follower, const Piece &P, std::vector<double> &Latencies,
                                  std::map<int, double> *Detected) {
    std::array<float, BlockSize> Block{};
    size_t Consumed = 0;
    size_t Misses = 0;
    const double Deadline = 1e6 * BlockSize / P.SampleRate;
    while (Consumed < P.Samples.size()) {
        const size_t Frames = std::min(BlockSize, P.Samples.size() - Consumed);
        std::copy_n(P.Samples.data() + Consumed, Frames, Block.data());
        const Clock::time_point Start = Clock::now();
        const bool Ok = Follower.ProcessBlock(Block.data(), Frames);
        const Clock::time_point End = Clock::now();
        if (!Ok) {
            throw std::runtime_error("ProcessBlock failed: " + std::string(P.Name));
        }
        const double Microseconds = std::chrono::duration<double, std::micro>(End - Start).count();
        Latencies.push_back(Microseconds);
        if (Microseconds > Deadline) {
            ++Misses;
        }
        Consumed += Frames;
        RecordDetection(Follower, P, Consumed, Detected);
    }
    return Misses;
}

// ─────────────────────────────────────
static Result RunInference(const Piece &P, bool MeasureLatency, std::vector<double> &Latencies,
                           std::map<int, double> *Detected) {
    OpenScofo::OpenScofo Follower(P.SampleRate, FFTSize, HopSize);
    InitializeInference(Follower, P);

    Result R;
    const Clock::time_point Start = Clock::now();
    if (MeasureLatency) {
        R.DeadlineMisses = ProcessAudioLatency(Follower, P, Latencies, Detected);
    } else {
        ProcessAudio(Follower, P, Detected);
    }
    const Clock::time_point End = Clock::now();
    R.Seconds = std::chrono::duration<double>(End - Start).count();
    R.FinalPosition = Follower.GetCurrentScorePosition();
    return R;
}

// ─────────────────────────────────────
static double Percentile(const std::vector<double> &Sorted, double P) {
    const size_t Index = static_cast<size_t>(std::ceil(P * Sorted.size())) - 1;
    return Sorted[Index];
}

// ─────────────────────────────────────
static Accuracy EvaluateAccuracy(const Piece &P, const std::map<int, double> &Detected) {
    Accuracy A;
    A.ReferenceEvents = P.Expected.size();
    double AbsoluteSum = 0.0;
    double SignedSum = 0.0;
    double SquaredSum = 0.0;

    for (const auto &Ref : P.Expected) {
        const auto Found = Detected.find(Ref.first);
        if (Found == Detected.end() || std::abs(Found->second - Ref.second) > AlignmentToleranceSeconds) {
            ++A.Missed;
            continue;
        }
        const double OffsetMs = (Found->second - Ref.second) * 1000.0;
        ++A.Matched;
        AbsoluteSum += std::abs(OffsetMs);
        SignedSum += OffsetMs;
        SquaredSum += OffsetMs * OffsetMs;
    }

    // In Cont (2010), false positives are misaligned detections and may overlap with missed events.
    for (const auto &Detection : Detected) {
        const auto Ref = P.Expected.find(Detection.first);
        if (Ref == P.Expected.end() || std::abs(Detection.second - Ref->second) > AlignmentToleranceSeconds) {
            ++A.FalsePositives;
        }
    }

    A.PrecisionPercent = 100.0 * A.Matched / A.ReferenceEvents;
    if (A.Matched != 0) {
        A.MeanAbsoluteOffsetMs = AbsoluteSum / A.Matched;
        A.MeanSignedOffsetMs = SignedSum / A.Matched;
        const double Variance = SquaredSum / A.Matched - A.MeanSignedOffsetMs * A.MeanSignedOffsetMs;
        A.StddevOffsetMs = std::sqrt(std::max(0.0, Variance));
    }
    return A;
}

// ─────────────────────────────────────
static Benchmark BenchmarkPiece(const Piece &P, int Repetitions, bool MeasureLatency) {
    const double AudioSeconds = static_cast<double>(P.Samples.size()) / P.SampleRate;
    const size_t Blocks = (P.Samples.size() - 1) / BlockSize + 1;
    std::vector<double> Latencies;
    std::vector<double> Times;
    Times.reserve(Repetitions);
    if (MeasureLatency) {
        if (Blocks > Latencies.max_size() / static_cast<size_t>(Repetitions)) {
            throw std::runtime_error("Too many latency samples");
        }
        Latencies.reserve(Blocks * static_cast<size_t>(Repetitions));
    }

    // Warm up library code and caches. Not included in the measurements.
    RunInference(P, false, Latencies, nullptr);

    Benchmark B;
    bool HasAccuracy = false;
    bool HasPosition = false;
    for (int i = 0; i < Repetitions; ++i) {
        std::map<int, double> Detected;
        const Result R = RunInference(P, MeasureLatency, Latencies, &Detected);
        const Accuracy A = EvaluateAccuracy(P, Detected);
        if (!HasAccuracy) {
            B.WorstAccuracy = A;
            HasAccuracy = true;
        } else {
            // Track the worst value of EACH metric across repetitions, not just the worst run.
            // Thus one abnormal run cannot be hidden by another abnormal run of a different kind.
            if (A.Matched < B.WorstAccuracy.Matched) {
                B.WorstAccuracy.Matched = A.Matched;
                B.WorstAccuracy.Missed = A.Missed;
                B.WorstAccuracy.PrecisionPercent = A.PrecisionPercent;
            }
            B.WorstAccuracy.FalsePositives = std::max(B.WorstAccuracy.FalsePositives, A.FalsePositives);
            if (A.MeanAbsoluteOffsetMs > B.WorstAccuracy.MeanAbsoluteOffsetMs) {
                B.WorstAccuracy.MeanAbsoluteOffsetMs = A.MeanAbsoluteOffsetMs;
                B.WorstAccuracy.MeanSignedOffsetMs = A.MeanSignedOffsetMs;
                B.WorstAccuracy.StddevOffsetMs = A.StddevOffsetMs;
            }
        }
        if (HasPosition && R.FinalPosition != B.FinalPosition) {
            throw std::runtime_error(std::string(P.Name) + ": inconsistent final position across runs");
        }
        B.FinalPosition = R.FinalPosition;
        HasPosition = true;
        Times.push_back(R.Seconds * 1000.0);
        B.DeadlineMisses += R.DeadlineMisses;
    }

    std::sort(Times.begin(), Times.end());
    B.MeanRunMs = std::accumulate(Times.begin(), Times.end(), 0.0) / Times.size();
    B.MedianRunMs = Percentile(Times, 0.50);
    B.MinRunMs = Times.front();
    B.MaxRunMs = Times.back();
    B.BlockCalls = Latencies.size();
    if (MeasureLatency) {
        std::sort(Latencies.begin(), Latencies.end());
        B.P99BlockUs = Percentile(Latencies, 0.99);
    }

    std::cout << std::fixed << std::setprecision(3) << "\n"
              << P.Name << " (" << Repetitions << " runs)\n"
              << "  Audio:           " << AudioSeconds << " s @ " << P.SampleRate << " Hz\n"
              << "  Mean run:        " << B.MeanRunMs << " ms\n"
              << "  Median run:      " << B.MedianRunMs << " ms\n"
              << "  Min / max:       " << B.MinRunMs << " / " << B.MaxRunMs << " ms\n"
              << "  Real-time factor: " << (B.MeanRunMs / (1000.0 * AudioSeconds)) << "\n"
              << "  Final position:  " << B.FinalPosition << "\n"
              << "  Precision:       " << B.WorstAccuracy.PrecisionPercent << "% (" << B.WorstAccuracy.Matched << "/"
              << B.WorstAccuracy.ReferenceEvents << ")\n"
              << "  Missed / FP:     " << B.WorstAccuracy.Missed << " / " << B.WorstAccuracy.FalsePositives << "\n"
              << "  Mean abs offset: " << B.WorstAccuracy.MeanAbsoluteOffsetMs << " ms\n"
              << "  Mean offset:     " << B.WorstAccuracy.MeanSignedOffsetMs << " ms\n"
              << "  Offset stddev:   " << B.WorstAccuracy.StddevOffsetMs << " ms\n";

    if (MeasureLatency) {
        const double Average = std::accumulate(Latencies.begin(), Latencies.end(), 0.0) / Latencies.size();
        std::cout << "  Calls:           " << Latencies.size() << "\n"
                  << "  Mean block:      " << Average << " us\n"
                  << "  P50 / P95 / P99: " << Percentile(Latencies, 0.50) << " / " << Percentile(Latencies, 0.95)
                  << " / " << B.P99BlockUs << " us\n"
                  << "  Max block:       " << Latencies.back() << " us\n"
                  << "  Block budget:    " << 1e6 * BlockSize / P.SampleRate << " us\n"
                  << "  Budget misses:   " << B.DeadlineMisses << " / " << Latencies.size() << "\n";
    }
    return B;
}

// ─────────────────────────────────────
static Json AccuracyToJson(const Accuracy &A) {
    return {{"reference_events", A.ReferenceEvents},
            {"matched", A.Matched},
            {"missed", A.Missed},
            {"false_positives", A.FalsePositives},
            {"precision_percent", A.PrecisionPercent},
            {"mean_absolute_offset_ms", A.MeanAbsoluteOffsetMs},
            {"mean_signed_offset_ms", A.MeanSignedOffsetMs},
            {"stddev_offset_ms", A.StddevOffsetMs}};
}

// ─────────────────────────────────────
static Json ReadRecord(const fs::path &Path) {
    Json Record;
    if (!fs::exists(Path)) {
        Record = {{"schema_version", 1},
                  {"method", "Cont 2010, section 8.2; first detection per event"},
                  {"aggregation", "worst value per metric across all repetitions"},
                  {"tolerance_ms", AlignmentToleranceSeconds * 1000.0},
                  {"configuration", {{"block_size", BlockSize}, {"fft_size", FFTSize}, {"hop_size", HopSize}}},
                  {"pieces", Json::object()}};
        std::cout << "No previous record; a new one will be created: " << Path << "\n";
        return Record;
    }
    std::ifstream File(Path);
    if (!File) {
        throw std::runtime_error("Could not open previous record: " + Path.string());
    }
    File >> Record; // Invalid JSON fails instead of silently overwriting the record.
    if (!Record.is_object() || Record.at("schema_version") != 1 ||
        Record.at("tolerance_ms").get<double>() != AlignmentToleranceSeconds * 1000.0 ||
        Record.at("configuration").at("block_size").get<size_t>() != BlockSize ||
        Record.at("configuration").at("fft_size").get<int>() != FFTSize ||
        Record.at("configuration").at("hop_size").get<int>() != HopSize || !Record.at("pieces").is_object()) {
        throw std::runtime_error("Incompatible benchmark record: " + Path.string());
    }
    std::cout << "Loaded previous best record: " << Path << "\n";
    return Record;
}

// ─────────────────────────────────────
static bool CompareAndUpdateRecord(Json &Record, const Piece &P, const Benchmark &B, bool MeasureLatency,
                                   int Repetitions) {
    Json &Pieces = Record.at("pieces");
    Json &Entry = Pieces[P.Name];
    const Json CurrentAccuracy = AccuracyToJson(B.WorstAccuracy);
    const std::string Mode = MeasureLatency ? "latency" : "profile";
    bool Regression = false;

    if (Entry.is_null()) {
        Entry = {{"sample_rate", P.SampleRate}, {"accuracy", CurrentAccuracy}, {"performance", Json::object()}};
        std::cout << "  [NEW] " << P.Name << " accuracy baseline\n";
    } else {
        if (Entry.at("sample_rate").get<int>() != P.SampleRate ||
            Entry.at("accuracy").at("reference_events").get<size_t>() != B.WorstAccuracy.ReferenceEvents) {
            throw std::runtime_error(std::string(P.Name) +
                                     ": annotations or sample rate changed; baseline not comparable");
        }
        const Json &Old = Entry.at("accuracy");
        const size_t OldMatched = Old.at("matched").get<size_t>();
        const size_t OldFP = Old.at("false_positives").get<size_t>();
        const double OldOffset = Old.at("mean_absolute_offset_ms").get<double>();

        if (B.WorstAccuracy.Matched < OldMatched) {
            std::cerr << "  [REGRESSION] " << P.Name << " precision: " << Old.at("precision_percent").get<double>()
                      << "% -> " << B.WorstAccuracy.PrecisionPercent << "%\n";
            Regression = true;
        } else if (B.WorstAccuracy.Matched == OldMatched) {
            if (B.WorstAccuracy.FalsePositives > OldFP) {
                std::cerr << "  [REGRESSION] " << P.Name << " false positives: " << OldFP << " -> "
                          << B.WorstAccuracy.FalsePositives << "\n";
                Regression = true;
            }
            if (B.WorstAccuracy.Matched != 0 && B.WorstAccuracy.MeanAbsoluteOffsetMs > OldOffset + OffsetEpsilonMs) {
                std::cerr << "  [REGRESSION] " << P.Name << " mean absolute offset: " << OldOffset << " -> "
                          << B.WorstAccuracy.MeanAbsoluteOffsetMs << " ms\n";
                Regression = true;
            }
        }

        // Accuracy takes priority; with equal hits, keep the lowest FP count and temporal error.
        if (B.WorstAccuracy.Matched > OldMatched ||
            (B.WorstAccuracy.Matched == OldMatched && B.WorstAccuracy.FalsePositives <= OldFP &&
             (B.WorstAccuracy.Matched == 0 || B.WorstAccuracy.MeanAbsoluteOffsetMs + OffsetEpsilonMs < OldOffset))) {
            Entry["accuracy"] = CurrentAccuracy;
            std::cout << "  [BEST] " << P.Name << " accuracy/offset record improved\n";
        }
    }

    Json &Performance = Entry.at("performance");
    const bool HadMode = Performance.contains(Mode);
    if (!HadMode) {
        Performance[Mode] = {{"best_median_run_ms", B.MedianRunMs},
                             {"mean_run_ms_at_best", B.MeanRunMs},
                             {"min_run_ms_at_best", B.MinRunMs},
                             {"max_run_ms_at_best", B.MaxRunMs},
                             {"repetitions_at_best", Repetitions}};
        if (MeasureLatency) {
            Performance[Mode]["best_p99_block_us"] = B.P99BlockUs;
            Performance[Mode]["least_budget_misses"] = B.DeadlineMisses;
        }
        std::cout << "  [NEW] " << P.Name << " " << Mode << " performance baseline\n";
    } else {
        Json &Previous = Performance.at(Mode);
        const double BestMedian = Previous.at("best_median_run_ms").get<double>();
        if (B.MedianRunMs > BestMedian * (1.0 + SpeedRegressionTolerance)) {
            std::cerr << "  [REGRESSION] " << P.Name << " " << Mode << " median: " << BestMedian << " -> "
                      << B.MedianRunMs << " ms (more than 10% slower)\n";
            Regression = true;
        }
        if (B.MedianRunMs < BestMedian) {
            Previous["best_median_run_ms"] = B.MedianRunMs;
            Previous["mean_run_ms_at_best"] = B.MeanRunMs;
            Previous["min_run_ms_at_best"] = B.MinRunMs;
            Previous["max_run_ms_at_best"] = B.MaxRunMs;
            Previous["repetitions_at_best"] = Repetitions;
            std::cout << "  [BEST] " << P.Name << " " << Mode << " median improved\n";
        }
        if (MeasureLatency) {
            const double BestP99 = Previous.at("best_p99_block_us").get<double>();
            const size_t LeastMisses = Previous.at("least_budget_misses").get<size_t>();
            if (B.P99BlockUs > BestP99 * (1.0 + SpeedRegressionTolerance)) {
                std::cerr << "  [REGRESSION] " << P.Name << " P99 block: " << BestP99 << " -> " << B.P99BlockUs
                          << " us\n";
                Regression = true;
            }
            if (B.DeadlineMisses > LeastMisses) {
                std::cerr << "  [REGRESSION] " << P.Name << " deadline misses: " << LeastMisses << " -> "
                          << B.DeadlineMisses << "\n";
                Regression = true;
            }
            if (B.P99BlockUs < BestP99) {
                Previous["best_p99_block_us"] = B.P99BlockUs;
            }
            if (B.DeadlineMisses < LeastMisses) {
                Previous["least_budget_misses"] = B.DeadlineMisses;
            }
        }
    }
    return Regression;
}

// ─────────────────────────────────────
static void WriteRecord(const fs::path &Path, const Json &Record) {
    fs::create_directories(Path.parent_path());
    const fs::path Temp = Path.string() + ".tmp";
    {
        std::ofstream File(Temp, std::ios::trunc);
        if (!File) {
            throw std::runtime_error("Could not write temporary record: " + Temp.string());
        }
        File << std::setw(2) << Record << '\n';
        File.flush();
        if (!File) {
            throw std::runtime_error("Could not flush temporary record: " + Temp.string());
        }
    }
    fs::rename(Temp, Path); // Atomic replace on POSIX; previous record survives failed writes.
    std::cout << "Saved best record: " << Path << "\n";
}

// ─────────────────────────────────────
int main(int argc, char **argv) {
    try {
        int Repetitions = 1;
        if (argc > 1) {
            size_t Parsed = 0;
            Repetitions = std::stoi(argv[1], &Parsed);
            if (Parsed != std::string(argv[1]).size()) {
                throw std::runtime_error("Invalid repetition count");
            }
        }
        const std::string Selection = argc > 2 ? argv[2] : "all";
        const std::string Mode = argc > 3 ? argv[3] : "profile";
        if (argc > 4 || Repetitions <= 0 || (Mode != "profile" && Mode != "latency")) {
            throw std::runtime_error("Usage: inference_benchmark [repetitions] "
                                     "[all|Miniatura1|Miniatura2|Miniatura3|Canticos] [profile|latency]");
        }

        // Always load and validate the best-known record before the first benchmark run.
        const Json Original = ReadRecord(RecordPath);
        Json Candidate = Original;
        const bool MeasureLatency = Mode == "latency";
        std::vector<Piece> Pieces = GetPieces();
        bool Found = false;
        for (Piece &P : Pieces) {
            if (Selection == "all" || Selection == P.Name) {
                LoadPiece(P);
                Found = true;
            }
        }
        if (!Found) {
            throw std::runtime_error("Unknown piece: " + Selection);
        }

        bool Regression = false;
        for (const Piece &P : Pieces) {
            if (Selection == "all" || Selection == P.Name) {
                const Benchmark B = BenchmarkPiece(P, Repetitions, MeasureLatency);
                if (CompareAndUpdateRecord(Candidate, P, B, MeasureLatency, Repetitions)) {
                    Regression = true;
                }
            }
        }

        if (Regression) {
            std::cerr << "\nFAILED: regression detected. Best record was NOT modified.\n";
            return 1;
        }
        if (Candidate != Original) {
            WriteRecord(RecordPath, Candidate);
        } else {
            std::cout << "\nPASS: no regression; previous best record unchanged.\n";
        }
        return 0;
    } catch (const std::exception &Error) {
        std::cerr << Error.what() << "\n";
        return 1;
    }
}
