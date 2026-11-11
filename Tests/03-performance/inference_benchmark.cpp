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
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

constexpr size_t BlockSize = 64;
constexpr int FFTSize = 2048;
constexpr int HopSize = 512;
const fs::path Data = TEST_DATA_DIR;

struct Piece {
    const char *Name;
    fs::path JsonPath;
    fs::path ScorePath;
    fs::path AudioDirectory;
    int SampleRate = 0;
    std::vector<float> Samples;
};

struct Result {
    double Seconds = 0.0;
    int FinalPosition = 0;
    size_t DeadlineMisses = 0;
};

// The same score/annotation pairs as inference_tests.
static std::vector<Piece> GetPieces() {
    return {{"Miniatura1",
             Data / "04-miniaturas/Extras/miniatura1.json",
             Data / "04-miniaturas/Extras/miniatura1.scofo",
             Data / "04-miniaturas/Audios",
             0,
             {}},
            {"Miniatura2",
             Data / "04-miniaturas/Extras/miniatura2.json",
             Data / "04-miniaturas/Extras/miniatura2.scofo",
             Data / "04-miniaturas/Audios",
             0,
             {}},
            {"Miniatura3",
             Data / "04-miniaturas/Extras/miniatura3.json",
             Data / "04-miniaturas/Extras/miniatura3.scofo",
             Data / "04-miniaturas/Audios",
             0,
             {}},
            {"Canticos",
             Data / "01-benchmark/real/canticos.json",
             Data / "01-benchmark/real/canticos.scofo",
             Data / "01-benchmark/real",
             0,
             {}}};
}

static void LoadWav(Piece &P, const fs::path &Path) {
    AudioFile<float> Audio;
    if (!Audio.load(Path.string()) || Audio.getNumChannels() == 0 || Audio.getNumSamplesPerChannel() == 0) {
        throw std::runtime_error("Could not load WAV: " + Path.string());
    }
    P.SampleRate = Audio.getSampleRate();
    P.Samples = Audio.samples[0];
}

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

static void LoadPiece(Piece &P) {
    std::ifstream File(P.JsonPath);
    if (!File) {
        throw std::runtime_error("Could not open JSON: " + P.JsonPath.string());
    }
    nlohmann::json Annotations;
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
}

// Same initial configuration as test_benchmark.cpp. LoadScore may override it.
static void InitializeInference(OpenScofo::OpenScofo &Follower, const Piece &P) {
    if (!Follower.LoadScore(P.ScorePath.string())) {
        throw std::runtime_error("Could not load score: " + P.ScorePath.string());
    }
    Follower.SetCurrentEvent(0);
}

// No per-block timing or allocation in the profiling loop.
static void ProcessAudio(OpenScofo::OpenScofo &Follower, const Piece &P) {
    std::array<float, BlockSize> Block{};
    size_t Consumed = 0;
    while (Consumed < P.Samples.size()) {
        const size_t Frames = std::min(BlockSize, P.Samples.size() - Consumed);
        std::copy_n(P.Samples.data() + Consumed, Frames, Block.data());
        if (!Follower.ProcessBlock(Block.data(), Frames)) {
            throw std::runtime_error("ProcessBlock failed: " + std::string(P.Name));
        }
        Consumed += Frames;
    }
}

// Measure every 64-sample ProcessBlock call, including expensive hop boundaries.
// The caller preallocates Latencies; no logging or sorting occurs here.
static size_t ProcessAudioLatency(OpenScofo::OpenScofo &Follower, const Piece &P, std::vector<double> &Latencies) {
    std::array<float, BlockSize> Block{};
    size_t Consumed = 0;
    size_t Misses = 0;
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
        // A real Pd callback has a full block period, even for the last short file block.
        const double Deadline = 1e6 * BlockSize / P.SampleRate;
        if (Microseconds > Deadline) {
            ++Misses;
        }
        Consumed += Frames;
    }
    return Misses;
}

static Result RunInference(const Piece &P, bool MeasureLatency, std::vector<double> &Latencies) {
    // A completely new inference object for every run.
    OpenScofo::OpenScofo Follower(P.SampleRate, FFTSize, HopSize);
    InitializeInference(Follower, P);

    Result R;
    const Clock::time_point Start = Clock::now();
    if (MeasureLatency) {
        R.DeadlineMisses = ProcessAudioLatency(Follower, P, Latencies);
    } else {
        ProcessAudio(Follower, P);
    }
    const Clock::time_point End = Clock::now();
    R.Seconds = std::chrono::duration<double>(End - Start).count();
    R.FinalPosition = Follower.GetCurrentScorePosition();
    return R;
}

static double Percentile(const std::vector<double> &Sorted, double P) {
    const size_t Index = static_cast<size_t>(std::ceil(P * Sorted.size())) - 1;
    return Sorted[Index];
}

static void BenchmarkPiece(const Piece &P, int Repetitions, bool MeasureLatency) {
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

    // Warm up library code and memory caches. This run is not timed in the summary.
    std::vector<double> Unused;
    RunInference(P, false, Unused);

    size_t Misses = 0;
    int FinalPosition = 0;
    for (int i = 0; i < Repetitions; ++i) {
        const Result R = RunInference(P, MeasureLatency, Latencies);
        Times.push_back(R.Seconds);
        Misses += R.DeadlineMisses;
        FinalPosition = R.FinalPosition;
    }

    std::sort(Times.begin(), Times.end());
    const double Mean = std::accumulate(Times.begin(), Times.end(), 0.0) / Times.size();
    const double RTF = Mean / AudioSeconds;
    std::cout << std::fixed << std::setprecision(3) << "\n"
              << P.Name << " (" << Repetitions << " runs)\n"
              << "  Audio:          " << AudioSeconds << " s @ " << P.SampleRate << " Hz\n"
              << "  Mean run:       " << Mean * 1000.0 << " ms\n"
              << "  Min / max:      " << Times.front() * 1000.0 << " / " << Times.back() * 1000.0 << " ms\n"
              << "  Real-time factor: " << RTF << " (" << 1.0 / RTF << "x real-time speed)\n"
              << "  Final position: " << FinalPosition << "\n";

    if (MeasureLatency) {
        std::sort(Latencies.begin(), Latencies.end());
        const double Avg = std::accumulate(Latencies.begin(), Latencies.end(), 0.0) / Latencies.size();
        std::cout << "  Calls:          " << Latencies.size() << "\n"
                  << "  Mean block:     " << Avg << " us\n"
                  << "  P50 block:      " << Percentile(Latencies, 0.50) << " us\n"
                  << "  P95 block:      " << Percentile(Latencies, 0.95) << " us\n"
                  << "  P99 block:      " << Percentile(Latencies, 0.99) << " us\n"
                  << "  Max block:      " << Latencies.back() << " us\n"
                  << "  Block budget:   " << 1e6 * BlockSize / P.SampleRate << " us\n"
                  << "  Budget misses:  " << Misses << " / " << Latencies.size() << "\n";
    }
}

int main(int argc, char **argv) {
    try {
        int Repetitions = 20;
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
        const bool MeasureLatency = Mode == "latency";
        std::vector<Piece> Pieces = GetPieces();
        bool Found = false;
        // Load every selected recording before running any inference.
        for (Piece &P : Pieces) {
            if (Selection == "all" || Selection == P.Name) {
                LoadPiece(P);
                Found = true;
            }
        }
        if (!Found) {
            throw std::runtime_error("Unknown piece: " + Selection);
        }
        for (const Piece &P : Pieces) {
            if (Selection == "all" || Selection == P.Name) {
                BenchmarkPiece(P, Repetitions, MeasureLatency);
            }
        }
        return 0;
    } catch (const std::exception &Error) {
        std::cerr << Error.what() << "\n";
        return 1;
    }
}
