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

    constexpr double inv = 1.0 / 32768.0;
    std::vector<double> wave(info.samples);
    for (int i = 0; i < info.samples; i++) {
        wave[i] = static_cast<double>(info.buffer[i]) * inv;
    }

    double maxAbs = 0.0;
    for (double s : wave)
        maxAbs = std::max(maxAbs, std::abs(s));

    free(info.buffer);
    return wave;
}

// ─────────────────────────────────────
int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <audio.mp3> <score.txt>\n";
        return 1;
    }

    const int WINDOW = 4096;
    const int HOP = 1024;
    const char *audio_path = argv[1];
    const char *score_path = argv[2];
    int sr = 0;
    int ch = 0;

    OpenScofo::OpenScofo scofo(48000, 4096, 512);
    bool ok = scofo.ParseScore(score_path);
    if (!ok) {
        exit(-1);
    }

    std::vector<double> samples = load_mp3_as_wave(audio_path, sr, ch);
    if (sr != 48000) {
        std::cerr << "warning: samplerate = " << sr << "\n";
    }

    std::cout << "Samples in audio " << samples.size() << "\n";

    std::vector<double> window(WINDOW, 0.0f);
    int currentEvent = -1;
    for (size_t pos = 0; pos + HOP <= samples.size(); pos += HOP) {
        memmove(window.data(), window.data() + HOP, (WINDOW - HOP) * sizeof(float));
        memcpy(window.data() + (WINDOW - HOP), samples.data() + pos, HOP * sizeof(float));
        bool ok = scofo.ProcessBlock(window);
        int event = scofo.GetEventIndex();
        if (event != currentEvent) {
            std::cerr << "Current Event is " << event << "\n";
            currentEvent = event;
        }
        if (event > 100) {
            break;
        }

        if (!ok) {
            std::cerr << "ProcessBlock failed\n";
            break;
        }
    }

    return 0;
}
