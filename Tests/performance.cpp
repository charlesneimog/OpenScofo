#define MINIMP3_IMPLEMENTATION
#include <minimp3_ex.h>

#include <iostream>
#include <vector>
#include <OpenScofo.hpp>

// ─────────────────────────────────────
std::vector<double> load_mp3_as_wave(const char *path, int &sr, int &ch) {
    mp3dec_t dec;
    mp3dec_file_info_t info;
    mp3dec_init(&dec);

    if (mp3dec_load(&dec, path, &info, NULL, NULL) != 0)
        throw std::runtime_error("mp3 decode failed");

    sr = info.hz;
    ch = info.channels;

    std::vector<double> wave(info.samples);
    for (int i = 0; i < (int)info.samples; i++) {
        wave[i] = static_cast<double>(info.buffer[i]);
    }

    free(info.buffer);
    return wave;
}

// ─────────────────────────────────────
void run_scofo(OpenScofo::OpenScofo &scofo, const std::vector<double> &samples) {
    const int WINDOW = 2048;
    std::vector<double> window(WINDOW, 0.0); // contiguous sliding window
    int currentEvent = -1;

    size_t blockIndex = 0;
    const int SCORE_HOP = 512; // or whatever OpenScofo expects

    for (size_t i = 0; i < samples.size(); i += 64) {
        // shift window
        memmove(window.data(), window.data() + 64, (WINDOW - 64) * sizeof(double));
        memcpy(window.data() + (WINDOW - 64), samples.data() + i, 64 * sizeof(double));
        blockIndex += 64;

        if (blockIndex >= SCORE_HOP) {
            blockIndex = 0;
            bool ok = scofo.ProcessBlock(window);
            int event = scofo.GetEventIndex();
            if (event != currentEvent) {
                spdlog::info("Current Event is {}", event);
                currentEvent = event;
            }
            if (event > 100)
                break;
            if (!ok) {
                std::cerr << "ProcessBlock failed\n";
                break;
            }
        }
    }
}

// ─────────────────────────────────────
int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <audio.mp3> <score.txt>\n";
        return 1;
    }

    const char *audio_path = argv[1];
    const char *score_path = argv[2];
    int sr = 0;
    int ch = 0;

    // Initialize OpenScofo
    OpenScofo::OpenScofo scofo(48000, 2048, 512);
    bool ok = scofo.ParseScore(score_path);
    if (!ok) {
        std::cerr << "ParseScore failed\n";
        exit(-1);
    }

    // Load audio
    std::vector<double> samples = load_mp3_as_wave(audio_path, sr, ch);
    if (sr != 48000) {
        std::cerr << "warning: samplerate = " << sr << "\n";
    }
    std::cout << "Samples in audio " << samples.size() << "\n";

    // Measure performance
    run_scofo(scofo, samples); // run processing

    return 0;
}
