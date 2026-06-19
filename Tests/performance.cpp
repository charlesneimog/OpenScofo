#include "AudioFile.h"

#include <iostream>
#include <vector>
#include <cstring> // memmove, memcpy
#include <stdexcept>

#include <OpenScofo.hpp>
#include <spdlog/spdlog.h>

// ─────────────────────────────────────
// Load WAV (or AIFF) using AudioFile
std::vector<double> load_audio_as_wave(const char *path, int &sr, int &ch) {
    AudioFile<double> audioFile;

    if (!audioFile.load(path)) {
        throw std::runtime_error("audio load failed");
    }

    sr = audioFile.getSampleRate();
    ch = audioFile.getNumChannels();
    int numSamplesPerChannel = audioFile.getNumSamplesPerChannel();

    std::vector<double> wave;
    wave.reserve(numSamplesPerChannel * ch);

    // Interleave channels to match your previous MP3 layout
    for (int i = 0; i < numSamplesPerChannel; i++) {
        for (int c = 0; c < ch; c++) {
            // Scale to match previous int16-based pipeline
            wave.push_back(audioFile.samples[c][i]);
        }
    }

    return wave;
}

// ─────────────────────────────────────
void run_scofo(OpenScofo::OpenScofo &scofo, const std::vector<double> &samples) {
    const int WINDOW = 2048;
    const int SCORE_HOP = 512;
    const int HOP = 64;

    std::vector<double> window(WINDOW, 0.0);
    int currentEvent = -1;
    std::vector<float> CurrSamplesBlock(HOP, 0);
    int pos = 0;

    for (size_t i = 0; i + HOP <= samples.size(); i += HOP) {
        pos = 0;
        for (size_t j = i; j < i + HOP; j++) {
            CurrSamplesBlock[pos] = samples[j];
            pos++;
        }
        bool ok = scofo.ProcessBlock(CurrSamplesBlock.data(), HOP);
        if (!ok) {
            printf("score failed\n");
            exit(-1);
        }
        int score_event = scofo.GetCurrentScorePosition();
        if (score_event != currentEvent) {
            printf("currentEvent %d\n", score_event);
            currentEvent = score_event;
        }
    }
}

// ─────────────────────────────────────
int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <audio.wav> <score.txt>\n";
        return 1;
    }

    const char *audio_path = argv[1];
    const char *score_path = argv[2];

    int sr = 0;
    int ch = 0;

    // Initialize OpenScofo
    OpenScofo::OpenScofo scofo(48000, 2048, 512);

    if (!scofo.LoadScore(score_path)) {
        std::cerr << "ParseScore failed\n";
        return -1;
    }

    // Load audio (WAV)
    std::vector<double> samples = load_audio_as_wave(audio_path, sr, ch);

    if (sr != 48000) {
        std::cerr << "warning: samplerate = " << sr << " (expected 48000)\n";
    }

    std::cout << "Samples in audio: " << samples.size() << "\n";
    std::cout << "Channels: " << ch << "\n";

    // Run processing
    run_scofo(scofo, samples);

    return 0;
}
