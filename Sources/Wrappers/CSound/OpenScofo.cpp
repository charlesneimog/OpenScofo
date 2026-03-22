#define opaddr opadr
#include <plugin.h>
#undef opaddr

#ifdef _CR
#undef _CR
#endif

#include <OpenScofo.hpp>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

struct OpenScofoState {
    std::unique_ptr<OpenScofo::OpenScofo> scofo;
    std::vector<double> buffer;
    int fftSize = 2048;
    int hopSize = 512;
    int blockIndex = 0;
    bool scoreMode = true;
    OpenScofo::Descriptors descriptor = OpenScofo::Descriptors::RMS;
};

static std::string oscofo_to_lower_trim(const std::string &text)
{
    const auto first = std::find_if_not(text.begin(), text.end(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    });

    const auto last = std::find_if_not(text.rbegin(), text.rend(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    }).base();

    if (first >= last) {
        return {};
    }

    std::string result(first, last);
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return result;
}

static std::string oscofo_trim(const std::string &text)
{
    const auto first = std::find_if_not(text.begin(), text.end(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    });

    const auto last = std::find_if_not(text.rbegin(), text.rend(), [](unsigned char ch) {
        return std::isspace(ch) != 0;
    }).base();

    if (first >= last) {
        return {};
    }

    return std::string(first, last);
}

static std::string oscofo_string_arg(STRINGDAT *arg)
{
    if (arg == nullptr || arg->data == nullptr || arg->size <= 0) {
        return {};
    }

    const char *ptr = reinterpret_cast<const char *>(arg->data);
    const int nbytes = arg->size;
    const int len = static_cast<int>(strnlen(ptr, static_cast<size_t>(nbytes)));
    return std::string(ptr, ptr + len);
}

static bool oscofo_is_scalar_descriptor(OpenScofo::Descriptors descriptor)
{
    switch (descriptor) {
    case OpenScofo::Descriptors::MFCC:
    case OpenScofo::Descriptors::CHROMA:
    case OpenScofo::Descriptors::MELOGRAM:
    case OpenScofo::Descriptors::MAGNITUDE:
    case OpenScofo::Descriptors::POWER:
    case OpenScofo::Descriptors::ONNX:
    case OpenScofo::Descriptors::INVALID:
        return false;
    default:
        return true;
    }
}

struct OpenScofoOpcode : csnd::Plugin<4, 3> {
    OpenScofoState *state = nullptr;

    int init()
    {
        state = new OpenScofoState();
        csound->plugin_deinit(this);

        state->scofo = std::make_unique<OpenScofo::OpenScofo>(static_cast<float>(sr()), 2048.0f, 512.0f);
        state->fftSize = static_cast<int>(state->scofo->GetFFTSize());
        state->hopSize = static_cast<int>(state->scofo->GetHopSize());
        state->buffer.assign(static_cast<size_t>(state->fftSize), 0.0);

        const std::string mode = oscofo_to_lower_trim(oscofo_string_arg(&inargs.str_data(1)));
        const std::string arg = oscofo_trim(oscofo_string_arg(&inargs.str_data(2)));

        if (mode == "desc" || mode == "descriptor" || mode == "descriptors") {
            state->scoreMode = false;
            if (!arg.empty()) {
                const std::string descriptorId = oscofo_to_lower_trim(arg);
                state->descriptor = state->scofo->GetDescriptorsEnum(descriptorId.c_str());
            }
            if (!oscofo_is_scalar_descriptor(state->descriptor)) {
                state->descriptor = OpenScofo::Descriptors::RMS;
                csound->warning("[oscofo] descriptor mode expects a scalar descriptor; falling back to 'rms'.");
            }
        } else {
            state->scoreMode = true;
            if (!arg.empty()) {
                if (!state->scofo->ParseScore(arg)) {
                    csound->warning("[oscofo] failed to parse score: " + arg);
                } else {
                    state->scofo->SetCurrentEvent(0);
                }
            }
        }

        outargs[0] = FL(0.0);
        outargs[1] = FL(0.0);
        outargs[2] = FL(0.0);
        outargs[3] = FL(0.0);

        if (state->scoreMode && state->scofo->ScoreIsLoaded()) {
            outargs[1] = static_cast<MYFLT>(state->scofo->GetLiveBPM());
        }

        return OK;
    }

    int aperf()
    {
        outargs[2] = FL(0.0);

        if (state == nullptr || state->scofo == nullptr || state->buffer.empty() || nsmps == 0) {
            return OK;
        }

        csnd::AudioSig in(this, inargs(0));
        const int ksmps = static_cast<int>(nsmps);
        const int fftSize = state->fftSize;

        if (ksmps >= fftSize) {
            const int inputStart = static_cast<int>(offset + nsmps) - fftSize;
            for (int i = 0; i < fftSize; ++i) {
                state->buffer[static_cast<size_t>(i)] = static_cast<double>(in[inputStart + i]);
            }
        } else {
            std::move(state->buffer.begin() + ksmps, state->buffer.end(), state->buffer.begin());
            const int writeStart = fftSize - ksmps;
            for (int i = 0; i < ksmps; ++i) {
                state->buffer[static_cast<size_t>(writeStart + i)] = static_cast<double>(in[static_cast<int>(offset) + i]);
            }
        }

        state->blockIndex += ksmps;
        while (state->blockIndex >= state->hopSize) {
            state->blockIndex -= state->hopSize;

            if (state->scoreMode) {
                if (!state->scofo->ScoreIsLoaded()) {
                    continue;
                }

                if (state->scofo->ProcessBlock(state->buffer)) {
                    outargs[0] = static_cast<MYFLT>(state->scofo->GetEventIndex());
                    outargs[1] = static_cast<MYFLT>(state->scofo->GetLiveBPM());
                    outargs[2] = FL(1.0);
                }
            } else {
                OpenScofo::Description desc = state->scofo->GetAudioDescription(state->buffer);
                outargs[3] = static_cast<MYFLT>(state->scofo->GetDescriptionFloat(desc, state->descriptor));
                outargs[2] = FL(1.0);
            }
        }

        return OK;
    }

    int deinit()
    {
        delete state;
        state = nullptr;
        return OK;
    }
};

#include <modload.h>
void csnd::on_load(csnd::Csound *csound)
{
    csnd::plugin<OpenScofoOpcode>(csound, "oscofo", "kkkk", "aSS", csnd::thread::a);
}
