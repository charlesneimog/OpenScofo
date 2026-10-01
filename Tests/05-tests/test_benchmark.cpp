
#include <OpenScofo.hpp>
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

#include <AudioFile.h>
#define MINIMP3_IMPLEMENTATION
#include <minimp3_ex.h>

namespace fs = std::filesystem;
using Json = nlohmann::json;

constexpr double Tolerance = 0.250;
constexpr size_t BlockSize = 64;

const fs::path Data = BENCHMARK_DATA_DIR;

// ─────────────────────────────────────
double TestPiece(const fs::path &JsonPath, const fs::path &ScorePath, const fs::path &AudioDirectory) {
    // Read annotations
    std::ifstream File(JsonPath);
    if (!File.is_open()) {
        std::cerr << "Could not open JSON: " << JsonPath << "\n";
        return 0.0;
    }

    Json Annotations;
    File >> Annotations;
    std::map<int, double> Expected;

    for (const auto &Event : Annotations["events"]) {
        Expected[Event["event"].get<int>()] = Event["timestamp_seconds"].get<double>();
    }

    if (Expected.empty()) {
        return 0.0;
    }

    const fs::path AudioPath = AudioDirectory / Annotations["audio"].get<std::string>();
    auto ProcessAudio = [&](int SampleRate, size_t NumSamples, auto GetSample) {
        OpenScofo::OpenScofo Follower(SampleRate, 2048, 256);
        if (!Follower.LoadScore(ScorePath.string())) {
            std::cerr << "Could not load score: " << ScorePath << "\n";
            return 0.0;
        }
        Follower.SetCurrentEvent(0);
        std::map<int, double> Detected;
        std::vector<float> Block(BlockSize);
        size_t Consumed = 0;
        while (Consumed < NumSamples) {
            const size_t Frames = std::min(BlockSize, NumSamples - Consumed);
            for (size_t i = 0; i < Frames; ++i) {
                Block[i] = GetSample(Consumed + i);
            }

            if (!Follower.ProcessBlock(Block.data(), Frames)) {
                std::cerr << "Audio processing failed\n";
                return 0.0;
            }

            Consumed += Frames;
            const int Position = Follower.GetCurrentScorePosition();
            if (Position > 0) {
                const double Time = static_cast<double>(Consumed) / SampleRate;
                Detected.try_emplace(Position, Time);
            }
        }
        size_t Matched = 0;
        for (const auto &[Position, ExpectedTime] : Expected) {
            const auto Found = Detected.find(Position);
            if (Found == Detected.end()) {
                continue;
            }

            const double Error = std::abs(Found->second - ExpectedTime);
            if (Error <= Tolerance) {
                ++Matched;
            }
        }

        return 100.0 * Matched / Expected.size();
    };

    // MP3
    if (AudioPath.extension() == ".mp3") {
        mp3dec_t Decoder;
        mp3dec_file_info_t Info{};
        mp3dec_init(&Decoder);
        const int Error = mp3dec_load(&Decoder, AudioPath.string().c_str(), &Info, nullptr, nullptr);
        if (Error != 0 || Info.buffer == nullptr) {
            std::cerr << "Could not load MP3: " << AudioPath << "\n";
            return 0.0;
        }

        const int Channels = Info.channels;
        const int SampleRate = Info.hz;
        const size_t NumSamples = Info.samples / Channels;

        const double Percentage = ProcessAudio(
            SampleRate, NumSamples, [&](size_t i) { return static_cast<float>(Info.buffer[i * Channels]) / 32768.0f; });
        free(Info.buffer);
        return Percentage;
    }

    // WAV
    if (AudioPath.extension() == ".wav") {
        AudioFile<float> Audio;
        if (!Audio.load(AudioPath.string())) {
            std::cerr << "Could not load WAV: " << AudioPath << "\n";
            return 0.0;
        }

        const int SampleRate = Audio.getSampleRate();
        const size_t NumSamples = static_cast<size_t>(Audio.getNumSamplesPerChannel());
        if (Audio.getNumChannels() == 0 || NumSamples == 0) {
            return 0.0;
        }
        return ProcessAudio(SampleRate, NumSamples, [&](size_t i) { return Audio.samples[0][i]; });
    }

    std::cerr << "Unsupported audio format: " << AudioPath << "\n";
    return 0.0;
}

// ─────────────────────────────────────
TEST(Inference, Miniatura1) {
    const double Percentage = TestPiece(Data / "04-miniaturas/Extras/miniatura1.json",
                                        Data / "04-miniaturas/Extras/miniatura1.scofo", Data / "04-miniaturas/Audios");

    std::cout << "Miniatura I: " << Percentage << "%\n";

    EXPECT_GE(Percentage, 14.81);
}

// ─────────────────────────────────────
TEST(Inference, Miniatura2) {
    const double Percentage = TestPiece(Data / "04-miniaturas/Extras/miniatura2.json",
                                        Data / "04-miniaturas/Extras/miniatura2.scofo", Data / "04-miniaturas/Audios");

    std::cout << "Miniatura II: " << Percentage << "%\n";

    EXPECT_GE(Percentage, 100.0);
}

// ─────────────────────────────────────
TEST(Inference, Miniatura3) {
    const double Percentage = TestPiece(Data / "04-miniaturas/Extras/miniatura3.json",
                                        Data / "04-miniaturas/Extras/miniatura3.scofo", Data / "04-miniaturas/Audios");
    std::cout << "Miniatura III: " << Percentage << "%\n";
    EXPECT_GE(Percentage, 91.93);
}

// ─────────────────────────────────────
TEST(Inference, Canticos) {
    const double Percentage = TestPiece(Data / "01-benchmark/real/canticos.json",
                                        Data / "01-benchmark/real/canticos.scofo", Data / "01-benchmark/real");
    std::cout << "Canticos: " << Percentage << "%\n";
    EXPECT_GE(Percentage, 97.56);
}
