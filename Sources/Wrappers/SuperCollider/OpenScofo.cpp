#include <OpenScofo.hpp>
#include <SC_PlugIn.hpp>
#include <algorithm>
#include <cctype>
#include <exception>
#include <string>
#include <vector>

#ifdef SUPERNOVA
#error "Supernova is not supported by OpenScofo"
#endif

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

        if (isAudioRateIn(0)) {
            set_calc_function<ScOpenScofo, &ScOpenScofo::next_a>();
            next_a(1);
        } else {
            set_calc_function<ScOpenScofo, &ScOpenScofo::next_k>();
            next_k(1);
        }
    }

    ~ScOpenScofo() {
        delete m_OScofo;
    }

    void ParseScore(const char *path) {
        if (m_OScofo) {
            m_OScofo->ParseScore(path);
            m_LastEventIndex = -1;
        }
    }

    void SetFollowScore(bool follow) {
        m_FollowScore = follow;
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

            const auto begin = std::find_if_not(item.begin(), item.end(), [](unsigned char c) {
                return std::isspace(c) != 0;
            });
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
            return m_OScofo->GetEventIndex();
        }
        return -1;
    }

    bool GetDescriptorValues(const char *descriptorId, std::vector<float> &values) {
        if (!m_OScofo || !descriptorId) {
            return false;
        }

        if (!m_LastDescValid) {
            m_LastDesc = m_OScofo->GetAudioDescription(m_InBuffer);
            m_LastDescValid = true;
        }

        OpenScofo::Description desc = (m_FollowScore && m_OScofo->ScoreIsLoaded()) ? m_OScofo->GetDescription()
                                                                                      : m_LastDesc;

        const auto descriptor = m_OScofo->GetDescriptorsEnum(descriptorId);
        if (descriptor == OpenScofo::Descriptors::INVALID || descriptor == OpenScofo::Descriptors::ONNX) {
            return false;
        }

        switch (descriptor) {
        case OpenScofo::Descriptors::MFCC:
        case OpenScofo::Descriptors::CHROMA:
        case OpenScofo::Descriptors::MELOGRAM:
        case OpenScofo::Descriptors::MAGNITUDE:
        case OpenScofo::Descriptors::POWER: {
            auto &arrayValues = m_OScofo->GetDescriptionArray(desc, descriptor);
            values.assign(arrayValues.begin(), arrayValues.end());
            return true;
        }
        default: {
            const double value = m_OScofo->GetDescriptionFloat(desc, descriptor);
            values.assign(1, static_cast<float>(value));
            return true;
        }
        }
    }

  private:
    OpenScofo::OpenScofo *m_OScofo = nullptr;
    int m_BlockIndex = 0;
    int m_FFTSize = 2048;
    int m_HopSize = 512;
    bool m_FollowScore = true;
    bool m_EmitEventNotifications = false;
    bool m_LastDescValid = false;
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

        const int currentEvent = m_OScofo->GetEventIndex();
        if (currentEvent == m_LastEventIndex) {
            return;
        }

        m_LastEventIndex = currentEvent;
        float currentEventAsFloat = static_cast<float>(currentEvent);
        SendNodeReply(&mParent->mNode, 0, "/oscofo/currentEvent", 1, &currentEventAsFloat);
    }

    void next_a(int inNumSamples) {
        (void)inNumSamples;
        int audioSamples = mWorld->mFullRate.mBufLength;
        const float *inBuf = in(0);
        const int size = (int)m_InBuffer.size();

        // 2. Shift your internal FFT buffer left to make room
        if (audioSamples < size) {
            std::copy(m_InBuffer.begin() + audioSamples, m_InBuffer.end(), m_InBuffer.begin());
            // 3. Insert the new audio block at the end
            for (int i = 0; i < audioSamples; ++i) {
                m_InBuffer[size - audioSamples + i] = inBuf[i];
            }
        }

        m_BlockIndex += audioSamples;

        // 4. If we hit the 512 hop size, process it
        if (m_BlockIndex >= m_HopSize) {
            m_BlockIndex -= m_HopSize; // Reset for next hop

            if (m_OScofo) {
                if (m_FollowScore && m_OScofo->ScoreIsLoaded()) {
                    m_OScofo->ProcessBlock(m_InBuffer);
                    m_LastDesc = m_OScofo->GetDescription();
                } else {
                    m_LastDesc = m_OScofo->GetAudioDescription(m_InBuffer);
                }
                m_LastDescValid = true;
                EmitCurrentEventIfChanged();
            }
        }

        out(0)[0] = 0.0f;
    }

    void next_k(int inNumSamples) {
        (void)inNumSamples;
        out(0)[0] = 0.0f;
    }
};

// ─────────────────────────────────────
void cmdGetCurrentEvent(ScOpenScofo *unit, sc_msg_iter *args) {
    (void)args;
    float currentEvent = static_cast<float>(unit->GetEventIndex());
    SendNodeReply(&unit->mParent->mNode, 0, "/oscofo/currentEvent", 1, &currentEvent);
}

// ─────────────────────────────────────
void cmdParseScore(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *path = args->gets(); // Get the string argument
    if (path) {
        unit->ParseScore(path);
    }
}

// ─────────────────────────────────────
void cmdSetFollowScore(ScOpenScofo *unit, sc_msg_iter *args) {
    const int follow = args->geti();
    unit->SetFollowScore(follow != 0);
}

// ─────────────────────────────────────
void cmdSetEventNotifications(ScOpenScofo *unit, sc_msg_iter *args) {
    const int enabled = args->geti();
    unit->SetEventNotifications(enabled != 0);
}

// ─────────────────────────────────────
void cmdGetDescriptor(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *descriptorId = args->gets();
    if (!descriptorId) {
        return;
    }

    std::vector<float> values;
    if (!unit->GetDescriptorValues(descriptorId, values) || values.empty()) {
        return;
    }

    std::string address = std::string("/oscofo/descriptor/") + descriptorId;
    SendNodeReply(&unit->mParent->mNode, 0, address.c_str(), static_cast<int>(values.size()), values.data());
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
    DefineUnitCmd("OpenScofo", "setFollowScore", (UnitCmdFunc)&cmdSetFollowScore);
    DefineUnitCmd("OpenScofo", "setEventNotifications", (UnitCmdFunc)&cmdSetEventNotifications);
    DefineUnitCmd("OpenScofo", "getDescriptor", (UnitCmdFunc)&cmdGetDescriptor);
    DefineUnitCmd("OpenScofo", "loadOnnxModel", (UnitCmdFunc)&cmdLoadOnnxModel);

    Print("\nOpenScofo version %d.%d.%d (build on %s), by Charles K. Neimog\n\n", OSCOFO_VERSION_MAJOR,
          OSCOFO_VERSION_MINOR, OSCOFO_VERSION_PATCH, OSCOFO_BUILD_TIME);
}
