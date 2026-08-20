#include "AudioFile.h"

#include <OpenScofo.hpp>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// ─────────────────────────────────────
// Load WAV using AudioFile
std::vector<double> load_audio_as_wave(const fs::path &path, int &sr, int &ch) {
    AudioFile<double> audioFile;

    if (!audioFile.load(path.string())) {
        std::cerr << "Audio load failed: " << path << "\n";
        std::abort();
    }

    sr = audioFile.getSampleRate();
    ch = audioFile.getNumChannels();

    const int numSamplesPerChannel = audioFile.getNumSamplesPerChannel();

    std::vector<double> wave;
    wave.reserve(static_cast<size_t>(numSamplesPerChannel) * ch);

    // Interleave channels
    for (int i = 0; i < numSamplesPerChannel; ++i) {
        for (int c = 0; c < ch; ++c) {
            wave.push_back(audioFile.samples[c][i]);
        }
    }

    return wave;
}

// ─────────────────────────────────────
void run_scofo(OpenScofo::OpenScofo &scofo, const std::vector<double> &samples) {
    constexpr int HOP = 64;

    std::vector<float> currentSamplesBlock(HOP, 0.0f);

    int currentEvent = -1;

    for (size_t i = 0; i + HOP <= samples.size(); i += HOP) {
        for (int j = 0; j < HOP; ++j) {
            currentSamplesBlock[static_cast<size_t>(j)] = static_cast<float>(samples[i + static_cast<size_t>(j)]);
        }

        const bool ok = scofo.ProcessBlock(currentSamplesBlock.data(), HOP);

        if (!ok) {
            std::cerr << "Score processing failed\n";
            return;
        }

        const int scoreEvent = scofo.GetCurrentScorePosition();

        if (scoreEvent != currentEvent) {
            currentEvent = scoreEvent;
        }
    }
}

// ─────────────────────────────────────
bool process_score_pair(const fs::path &wavPath, const fs::path &scorePath, int repetitions) {
    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "Audio: " << wavPath.filename() << "\n";
    std::cout << "Score: " << scorePath.filename() << "\n";
    std::cout << "========================================\n";

    int sr = 0;
    int ch = 0;

    // Create a fresh OpenScofo instance for each score.
    OpenScofo::OpenScofo scofo(48000, 2048, 512);

    if (!scofo.LoadScore(scorePath.string())) {
        std::cerr << "LoadScore failed: " << scorePath << "\n";
        return false;
    }

    std::vector<double> samples = load_audio_as_wave(wavPath, sr, ch);

    if (sr != 48000) {
        std::cerr << "Warning: samplerate = " << sr << " Hz, expected 48000 Hz\n";
    }

    std::cout << "Samples: " << samples.size() << "\n";
    std::cout << "Channels: " << ch << "\n";

    for (int run = 0; run < repetitions; ++run) {
        std::cout << "Run " << (run + 1) << "/" << repetitions << "\n";

        scofo.SetCurrentEvent(0);
        run_scofo(scofo, samples);
    }

    return true;
}

// ─────────────────────────────────────
int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <directory> [repetitions]\n";

        return 1;
    }

    const fs::path directory = argv[1];

    int repetitions = 20;

    if (argc >= 3) {
        repetitions = std::max(1, std::atoi(argv[2]));
    }

    if (!fs::exists(directory)) {
        std::cerr << "Directory does not exist: " << directory << "\n";
        return 1;
    }

    if (!fs::is_directory(directory)) {
        std::cerr << "Path is not a directory: " << directory << "\n";
        return 1;
    }

    std::vector<fs::path> scoreFiles;

    // Find all .scofo files
    for (const auto &entry : fs::directory_iterator(directory)) {
        if (!entry.is_regular_file()) {
            continue;
        }

        const fs::path &path = entry.path();

        if (path.extension() == ".scofo") {
            scoreFiles.push_back(path);
        }
    }

    // Natural/simple alphabetical ordering:
    // score1.scofo, score10.scofo, score2.scofo...
    std::sort(scoreFiles.begin(), scoreFiles.end());

    if (scoreFiles.empty()) {
        std::cerr << "No .scofo files found in " << directory << "\n";
        return 1;
    }

    int processed = 0;
    int failed = 0;

    for (const fs::path &scorePath : scoreFiles) {
        fs::path wavPath = scorePath;

        // score1.scofo -> score1.wav
        wavPath.replace_extension(".wav");

        if (!fs::exists(wavPath)) {
            std::cerr << "Missing WAV for " << scorePath.filename() << ": expected " << wavPath.filename() << "\n";

            ++failed;
            continue;
        }

        if (process_score_pair(wavPath, scorePath, repetitions)) {
            ++processed;
        } else {
            ++failed;
        }
    }

    std::cout << "\n";
    std::cout << "========================================\n";
    std::cout << "Finished\n";
    std::cout << "Processed: " << processed << "\n";
    std::cout << "Failed:    " << failed << "\n";
    std::cout << "========================================\n";

    return failed == 0 ? 0 : 1;
}
