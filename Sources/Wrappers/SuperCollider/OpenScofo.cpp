#include <OpenScofo.hpp>
#include <SC_PlugIn.hpp>
#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

#ifdef SUPERNOVA
#error "Supernova is not supported by OpenScofo"
#endif

// Callback de erro adaptado para SC
static void OpenScofo_error_callback(const spdlog::details::log_msg &log, void *data) {
    // ignorar logs abaixo de warn
    if (log.level < spdlog::level::warn)
        return;

    std::string text(log.payload.data(), log.payload.size());

    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
        printf("\n[OpenScofo][ERR] %s\n", text.c_str());
        break;
    case spdlog::level::warn:
        printf("\n[OpenScofo][WARN] %s\n", text.c_str());
        break;
    case spdlog::level::info:
    case spdlog::level::debug:
    case spdlog::level::trace:
        printf("\n[OpenScofo][INFO] %s\n", text.c_str());
        break;
    default:
        break;
    }
}

static InterfaceTable *ft;
struct ScOpenScofo : public SCUnit {
  public:
    ScOpenScofo() {
        const float sampleRate = in0(1) > 0.0f ? in0(1) : 48000.0f;
        const float fftSize = in0(2) > 0.0f ? in0(2) : 2048.0f;
        const float hopSize = in0(3) > 0.0f ? in0(3) : 512.0f;

        m_FFTSize = std::max(1, static_cast<int>(fftSize));
        m_HopSize = std::max(1, static_cast<int>(hopSize));
        m_InBuffer.assign(m_FFTSize, 0.0);
        m_OScofo = new OpenScofo::OpenScofo(sampleRate, fftSize, hopSize);
        m_OScofo->SetErrorCallback(OpenScofo_error_callback, nullptr);

        if (isAudioRateIn(0)) {
            set_calc_function<ScOpenScofo, &ScOpenScofo::next_a>();
            next_a(1);
        }
    }

    ~ScOpenScofo() {
        delete m_OScofo;
    }

    void ParseScore(const char *path) {
        if (m_OScofo) {
            bool ok = m_OScofo->ParseScore(path);
            if (!ok) {
                return;
            }
            m_LastEventIndex = -1;
        }
    }

    void SetEventNotifications(bool enabled) {
        m_EmitEventNotifications = enabled;
    }

    bool LoadOnnxModel(const char *modelPath, const char *descriptorIdsCsv) {
        if (!m_OScofo || !modelPath || !descriptorIdsCsv) {
            return false;
        }

        std::vector<OpenScofo::Descriptors> descriptors;
        descriptors.reserve(16);

        std::string csv(descriptorIdsCsv);
        std::string item;
        for (size_t i = 0; i <= csv.size(); ++i) {
            const bool split = (i == csv.size()) || (csv[i] == ',');
            if (!split) {
                item.push_back(csv[i]);
                continue;
            }

            const auto begin =
                std::find_if_not(item.begin(), item.end(), [](unsigned char c) { return std::isspace(c) != 0; });
            const auto end = std::find_if_not(item.rbegin(), item.rend(), [](unsigned char c) {
                                 return std::isspace(c) != 0;
                             }).base();

            if (begin < end) {
                std::string descriptorId(begin, end);
                const auto descriptor = m_OScofo->GetDescriptorsEnum(descriptorId.c_str());
                if (descriptor != OpenScofo::Descriptors::INVALID) {
                    descriptors.push_back(descriptor);
                }
            }

            item.clear();
        }

        if (descriptors.empty()) {
            return false;
        }

        m_OScofo->LoadONNXModel(modelPath, descriptors);
        return true;
    }

    int GetEventIndex() {
        if (m_OScofo) {
            return m_OScofo->GetCurrentScorePosition();
        }
        return -1;
    }

  private:
    OpenScofo::OpenScofo *m_OScofo = nullptr;
    int m_FFTSize = 2048;
    int m_HopSize = 512;
    bool m_FollowScore = true;
    bool m_EmitEventNotifications = false;
    std::vector<double> m_InBuffer;
    OpenScofo::Description m_LastDesc;
    int m_LastEventIndex = -1;

    void EmitCurrentEventIfChanged() {
        if (!m_EmitEventNotifications || !m_OScofo) {
            return;
        }

        if (!m_FollowScore || !m_OScofo->ScoreIsLoaded()) {
            return;
        }

        const int currentEvent = m_OScofo->GetCurrentScorePosition();
        if (currentEvent == m_LastEventIndex) {
            return;
        }

        m_LastEventIndex = currentEvent;
        float currentEventAsFloat = static_cast<float>(currentEvent);
        SendNodeReply(&mParent->mNode, 0, "/OpenScofo/currentEvent", 1, &currentEventAsFloat);
    }

    void next_a(int inNumSamples) {
        (void)inNumSamples;
        int n = mWorld->mFullRate.mBufLength;
        const float *inBuf = in(0);
        bool ok = m_OScofo->ProcessBlock(inBuf, n);
        if (!ok) {
            // set fail
        }

        EmitCurrentEventIfChanged();

        out(0)[0] = 0.0f;
    }
};

// ─────────────────────────────────────
void cmdGetCurrentEvent(ScOpenScofo *unit, sc_msg_iter *args) {
    (void)args;
    float currentEvent = static_cast<float>(unit->GetEventIndex());
    SendNodeReply(&unit->mParent->mNode, 0, "/OpenScofo/currentEvent", 1, &currentEvent);
}

// ─────────────────────────────────────
void cmdParseScore(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *path = args->gets(); // Get the string argument
    if (path) {
        unit->ParseScore(path);
    }
}

// ─────────────────────────────────────
void cmdSetEventNotifications(ScOpenScofo *unit, sc_msg_iter *args) {
    const int enabled = args->geti();
    unit->SetEventNotifications(enabled != 0);
}

// ─────────────────────────────────────
void cmdLoadOnnxModel(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *modelPath = args->gets();
    const char *descriptorIdsCsv = args->gets();
    unit->LoadOnnxModel(modelPath, descriptorIdsCsv);
}

// ─────────────────────────────────────
PluginLoad(OpenScofoUGens) {
    ft = inTable;

    registerUnit<ScOpenScofo>(ft, "OpenScofo");
    DefineUnitCmd("OpenScofo", "parseScore", (UnitCmdFunc)&cmdParseScore);
    DefineUnitCmd("OpenScofo", "getCurrentEvent", (UnitCmdFunc)&cmdGetCurrentEvent);
    DefineUnitCmd("OpenScofo", "setEventNotifications", (UnitCmdFunc)&cmdSetEventNotifications);
    DefineUnitCmd("OpenScofo", "loadOnnxModel", (UnitCmdFunc)&cmdLoadOnnxModel);

    printf("\nOpenScofo version %s (%s), by Charles K. Neimog\n\n", OPENSCOFO_VERSION, OPENSCOFO_BUILD_TIME);
}
