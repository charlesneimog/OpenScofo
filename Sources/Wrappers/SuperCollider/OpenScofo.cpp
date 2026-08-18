/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include <OpenScofo.hpp>
#include <SC_PlugIn.hpp>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cmath>
#include <atomic>
#include <exception>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_set>
#include <variant>
#include <vector>

#ifdef SUPERNOVA
#error "Supernova is not supported by OpenScofo"
#endif

static InterfaceTable *ft;

// ─────────────────────────────────────
static void OpenScofo_error_callback(const spdlog::details::log_msg &log, void *data) {
    (void)data;
    // ignorar logs abaixo de warn
    if (log.level < spdlog::level::warn)
        return;

    std::string text(log.payload.data(), log.payload.size());

    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
        Print("[OpenScofo][ERR] %s\n", text.c_str());
        break;
    case spdlog::level::warn:
        Print("[OpenScofo][WARN] %s\n", text.c_str());
        break;
    case spdlog::level::info:
    case spdlog::level::debug:
    case spdlog::level::trace:
        Print("[OpenScofo][INFO] %s\n", text.c_str());
        break;
    default:
        break;
    }
}

static constexpr float kEncodedActionMagic = -900001.0f;
static constexpr float kActionNumberTag = 0.0f;
static constexpr float kActionStringTag = 1.0f;
static constexpr const char *kDefaultOSCNamespace = "openscofo";

static std::string NormalizeOSCNamespace(const char *namespaceName) {
    std::string result = namespaceName ? namespaceName : kDefaultOSCNamespace;

    while (!result.empty() && result.front() == '/') {
        result.erase(result.begin());
    }
    while (!result.empty() && result.back() == '/') {
        result.pop_back();
    }

    return result.empty() ? kDefaultOSCNamespace : result;
}

// ─────────────────────────────────────
struct ScOpenScofo : public SCUnit {
  public:
    ScOpenScofo() {
        const float sampleRate = in0(1) > 0.0f ? in0(1) : 48000.0f;
        OpenScofo::Configuration config;

        m_OpenScofo = new OpenScofo::OpenScofo(sampleRate, config.FFTSize, config.HOPSize);
        m_OpenScofo->SetErrorCallback(OpenScofo_error_callback, nullptr);

        if (isAudioRateIn(0)) {
            set_calc_function<ScOpenScofo, &ScOpenScofo::next_a>();
            next_a(1);
        }
    }

    // ─────────────────────────────────────
    ~ScOpenScofo() {
        if (m_LoadThread.joinable()) {
            m_LoadThread.join();
        }
        std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
        delete m_OpenScofo;
    }

    // ─────────────────────────────────────
    void LoadScore(const char *path) {
        if (!path || path[0] == '\0') {
            return;
        }

        if (m_LoadInProgress.exchange(true)) {
            Print("[OpenScofo][WARN] SuperCollider score load ignored because another load is still running.\n");
            return;
        }

        if (m_LoadThread.joinable()) {
            m_LoadThread.join();
        }

        std::string scorePath(path);
        m_LoadThread = std::thread([this, scorePath]() {
            try {
                std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
                if (m_OpenScofo) {
                    bool ok = m_OpenScofo->LoadScore(scorePath);
                    if (ok) {
                        m_LastEventIndex = -1;
                        m_LastActionStateIndex = -1;
                        m_ScheduledActions.clear();
                    }
                }
            } catch (const std::exception &e) {
                Print("[OpenScofo][ERR] SuperCollider score load failed: %s\n", e.what());
            } catch (...) {
                Print("[OpenScofo][ERR] SuperCollider score load failed with an unknown exception.\n");
            }
            m_LoadInProgress = false;
        });
    }

    // ─────────────────────────────────────
    void SetEventNotifications(bool enabled) {
        m_EmitEventNotifications = enabled;
    }

    // ─────────────────────────────────────
    void SetFollowScore(bool enabled) {
        m_FollowScore = enabled;
    }

    // ─────────────────────────────────────
    void SetNamespace(const char *namespaceName) {
        std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
        m_OSCNamespace = NormalizeOSCNamespace(namespaceName);
        m_WarnedMissingReceivers.clear();
    }

    // ─────────────────────────────────────
    void RegisterActionReceiver(const char *receiver) {
        std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
        if (receiver && receiver[0] != '\0') {
            m_ActionReceivers.insert(receiver);
            m_WarnedMissingReceivers.erase(receiver);
        }
    }

    // ─────────────────────────────────────
    void UnregisterActionReceiver(const char *receiver) {
        std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
        if (receiver && receiver[0] != '\0') {
            m_ActionReceivers.erase(receiver);
        }
    }

    // ─────────────────────────────────────
    bool LoadOnnxModel(const char *modelPath, const char *descriptorIdsCsv) {
        std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
        if (!m_OpenScofo || !modelPath || !descriptorIdsCsv) {
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
                const auto descriptor = m_OpenScofo->GetDescriptorsEnum(descriptorId.c_str());
                if (descriptor != OpenScofo::Descriptors::INVALID) {
                    descriptors.push_back(descriptor);
                }
            }

            item.clear();
        }

        if (descriptors.empty()) {
            return false;
        }

        m_OpenScofo->LoadONNXModel(modelPath, descriptors);
        return true;
    }

    // ─────────────────────────────────────
    void SendDescriptor(const char *descriptorId) {
        std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
        if (!m_OpenScofo || !descriptorId) {
            return;
        }

        const auto descriptor = m_OpenScofo->GetDescriptorsEnum(descriptorId);
        if (descriptor == OpenScofo::Descriptors::INVALID) {
            return;
        }

        m_OpenScofo->RequestDescriptor(descriptor);
        OpenScofo::Description desc = m_OpenScofo->GetDescription();
        const std::string replyAddress = DescriptorAddress(m_OpenScofo->GetDescriptionId(descriptor));

        try {
            switch (descriptor) {
            case OpenScofo::Descriptors::MFCC:
            case OpenScofo::Descriptors::CHROMA:
            case OpenScofo::Descriptors::LOGMEL:
            case OpenScofo::Descriptors::MAGNITUDE:
            case OpenScofo::Descriptors::POWERARRAY: {
                std::vector<double> values = m_OpenScofo->GetDescriptionArray(desc, descriptor);
                std::vector<float> replyValues;
                replyValues.reserve(values.size());
                for (double value : values) {
                    replyValues.push_back(static_cast<float>(value));
                }
                if (!replyValues.empty()) {
                    SendNodeReply(&mParent->mNode, 0, replyAddress.c_str(), static_cast<int>(replyValues.size()),
                                  replyValues.data());
                }
                break;
            }
            case OpenScofo::Descriptors::ONNX:
                for (const auto &onnxDesc : desc.ONNX) {
                    std::string onnxAddress = replyAddress + "/";
                    onnxAddress += onnxDesc.first;
                    float value = static_cast<float>(onnxDesc.second);
                    SendNodeReply(&mParent->mNode, 0, onnxAddress.c_str(), 1, &value);
                }
                break;
            default: {
                float value = static_cast<float>(m_OpenScofo->GetDescriptionFloat(desc, descriptor));
                SendNodeReply(&mParent->mNode, 0, replyAddress.c_str(), 1, &value);
                break;
            }
            }
        } catch (const std::exception &e) {
            Print("[OpenScofo][ERR] Failed to send descriptor '%s': %s\n", descriptorId, e.what());
        }
    }

    // ─────────────────────────────────────
    int GetEventIndex() {
        std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
        if (m_OpenScofo) {
            return m_OpenScofo->GetCurrentScorePosition();
        }
        return -1;
    }

    // ─────────────────────────────────────
    void SendCurrentEvent() {
        std::lock_guard<std::mutex> lock(m_OpenScofoMutex);
        if (!m_OpenScofo) {
            return;
        }

        const float currentEvent = static_cast<float>(m_OpenScofo->GetCurrentScorePosition());
        const std::string replyAddress = EventAddress();
        SendNodeReply(&mParent->mNode, 0, replyAddress.c_str(), 1, &currentEvent);
    }

  private:
    struct ScheduledAction {
        int64_t samplesUntilDispatch = 0;
        OpenScofo::ScoreAction action;
    };

    OpenScofo::OpenScofo *m_OpenScofo = nullptr;
    std::mutex m_OpenScofoMutex;
    std::thread m_LoadThread;
    std::atomic_bool m_LoadInProgress = false;
    bool m_FollowScore = true;
    bool m_EmitEventNotifications = false;
    std::unordered_set<std::string> m_ActionReceivers;
    std::unordered_set<std::string> m_WarnedMissingReceivers;
    std::string m_OSCNamespace = kDefaultOSCNamespace;
    std::vector<ScheduledAction> m_ScheduledActions;
    OpenScofo::Description m_LastDesc;
    int m_LastEventIndex = -1;
    int m_LastActionStateIndex = -1;

    bool EncodeActionArgs(const OpenScofo::ScoreAction &action, std::vector<float> &values) {
        values.clear();

        const bool hasStringArg = std::any_of(action.Args.begin(), action.Args.end(),
                                              [](const auto &arg) { return std::holds_alternative<std::string>(arg); });

        if (!hasStringArg) {
            values.reserve(action.Args.size());
            for (const auto &arg : action.Args) {
                if (std::holds_alternative<float>(arg)) {
                    values.push_back(std::get<float>(arg));
                } else if (std::holds_alternative<int>(arg)) {
                    values.push_back(static_cast<float>(std::get<int>(arg)));
                }
            }
            return true;
        }

        values.push_back(kEncodedActionMagic);
        values.push_back(static_cast<float>(action.Args.size()));

        for (const auto &arg : action.Args) {
            if (std::holds_alternative<float>(arg)) {
                values.push_back(kActionNumberTag);
                values.push_back(std::get<float>(arg));
            } else if (std::holds_alternative<int>(arg)) {
                values.push_back(kActionNumberTag);
                values.push_back(static_cast<float>(std::get<int>(arg)));
            } else if (std::holds_alternative<std::string>(arg)) {
                const auto &text = std::get<std::string>(arg);
                values.push_back(kActionStringTag);
                values.push_back(static_cast<float>(text.size()));
                for (const unsigned char ch : text) {
                    values.push_back(static_cast<float>(ch));
                }
            }
        }

        return true;
    }

    std::string ActionAddress(const std::string &receiver) const {
        std::string address = "/";
        address += m_OSCNamespace;
        address += "/";
        address += receiver;
        return address;
    }

    std::string EventAddress() const {
        return ActionAddress("currentEvent");
    }

    std::string DescriptorAddress(const char *descriptorId) const {
        std::string address = ActionAddress("descriptor");
        address += "/";
        address += descriptorId;
        return address;
    }

    bool HasActionReceiver(const std::string &receiver) {
        if (m_ActionReceivers.find(receiver) != m_ActionReceivers.end()) {
            return true;
        }

        if (m_WarnedMissingReceivers.insert(receiver).second) {
            const std::string address = ActionAddress(receiver);
            Print("[OpenScofo][WARN] SuperCollider sendto receiver '%s' has no registered listener. "
                  "Register it with listen(...) or registerActionReceiver(...).\n",
                  address.c_str());
        }

        return false;
    }

    // ─────────────────────────────────────
    void DispatchAction(const OpenScofo::ScoreAction &action) {
        if (action.isLua) {
            Print("[OpenScofo][WARN] SuperCollider wrapper does not execute Lua score actions yet.\n");
            return;
        }

        std::vector<float> values;
        if (!EncodeActionArgs(action, values)) {
            return;
        }

        if (action.Receiver.empty() || !HasActionReceiver(action.Receiver)) {
            return;
        }

        const std::string replyAddress = ActionAddress(action.Receiver);
        SendNodeReply(&mParent->mNode, 0, replyAddress.c_str(), static_cast<int>(values.size()),
                      values.empty() ? nullptr : values.data());
    }

    // ─────────────────────────────────────
    int64_t ActionDelaySamples(const OpenScofo::ScoreAction &action) {
        double timeMs = action.Time;
        if (!action.AbsoluteTime && m_OpenScofo) {
            timeMs = 60.0 / m_OpenScofo->GetCurrentBPM() * action.Time * 1000.0;
        }

        if (timeMs <= 0.0 || !m_OpenScofo) {
            return 0;
        }

        return static_cast<int64_t>(std::ceil(timeMs * m_OpenScofo->GetSr() / 1000.0));
    }

    // ─────────────────────────────────────
    void QueueOrDispatchAction(const OpenScofo::ScoreAction &action) {
        const int64_t delaySamples = ActionDelaySamples(action);
        if (delaySamples <= 0) {
            DispatchAction(action);
            return;
        }

        m_ScheduledActions.push_back({delaySamples, action});
    }

    // ─────────────────────────────────────
    void ProcessEventActionsIfChanged() {
        if (!m_OpenScofo || !m_FollowScore || !m_OpenScofo->ScoreIsLoaded()) {
            return;
        }

        OpenScofo::EventActions audioStateActions = m_OpenScofo->GetAudioStateChangeActions();
        for (const auto &action : audioStateActions) {
            DispatchAction(action);
        }

        const int currentStateIndex = m_OpenScofo->GetCurrentStateIndex();
        if (currentStateIndex == m_LastActionStateIndex) {
            return;
        }

        m_LastActionStateIndex = currentStateIndex;
        OpenScofo::EventActions actions = m_OpenScofo->GetCurrentEventActions();
        for (const auto &action : actions) {
            QueueOrDispatchAction(action);
        }
    }

    // ─────────────────────────────────────
    void ProcessScheduledActions(int samplesElapsed) {
        auto it = m_ScheduledActions.begin();
        while (it != m_ScheduledActions.end()) {
            it->samplesUntilDispatch -= samplesElapsed;
            if (it->samplesUntilDispatch <= 0) {
                DispatchAction(it->action);
                it = m_ScheduledActions.erase(it);
            } else {
                ++it;
            }
        }
    }

    // ─────────────────────────────────────
    void EmitCurrentEventIfChanged() {
        if (!m_EmitEventNotifications || !m_OpenScofo) {
            return;
        }

        if (!m_FollowScore || !m_OpenScofo->ScoreIsLoaded()) {
            return;
        }

        const int currentEvent = m_OpenScofo->GetCurrentScorePosition();
        if (currentEvent == m_LastEventIndex) {
            return;
        }

        m_LastEventIndex = currentEvent;
        float currentEventAsFloat = static_cast<float>(currentEvent);
        const std::string replyAddress = EventAddress();
        SendNodeReply(&mParent->mNode, 0, replyAddress.c_str(), 1, &currentEventAsFloat);
    }

    void next_a(int inNumSamples) {
        (void)inNumSamples;
        int n = mWorld->mFullRate.mBufLength;
        if (m_LoadInProgress) {
            std::fill_n(out(0), n, 0.0f);
            return;
        }

        std::unique_lock<std::mutex> lock(m_OpenScofoMutex, std::try_to_lock);
        if (!lock.owns_lock() || !m_OpenScofo) {
            std::fill_n(out(0), n, 0.0f);
            return;
        }

        const float *inBuf = in(0);
        bool ok = m_OpenScofo->ProcessBlock(inBuf, n);
        if (!ok) {
            // set fail
        }

        ProcessEventActionsIfChanged();
        ProcessScheduledActions(n);
        EmitCurrentEventIfChanged();

        std::fill_n(out(0), n, 0.0f);
    }
};

// ─────────────────────────────────────
void cmdGetCurrentEvent(ScOpenScofo *unit, sc_msg_iter *args) {
    (void)args;
    unit->SendCurrentEvent();
}

// ─────────────────────────────────────
void cmdLoadScore(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *path = args->gets(); // Get the string argument
    if (path) {
        unit->LoadScore(path);
    }
}

// ─────────────────────────────────────
void cmdSetEventNotifications(ScOpenScofo *unit, sc_msg_iter *args) {
    const int enabled = args->geti();
    unit->SetEventNotifications(enabled != 0);
}

// ─────────────────────────────────────
void cmdSetFollowScore(ScOpenScofo *unit, sc_msg_iter *args) {
    const int enabled = args->geti();
    unit->SetFollowScore(enabled != 0);
}

// ─────────────────────────────────────
void cmdSetNamespace(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *namespaceName = args->gets();
    unit->SetNamespace(namespaceName);
}

// ─────────────────────────────────────
void cmdRegisterActionReceiver(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *receiver = args->gets();
    unit->RegisterActionReceiver(receiver);
}

// ─────────────────────────────────────
void cmdUnregisterActionReceiver(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *receiver = args->gets();
    unit->UnregisterActionReceiver(receiver);
}

// ─────────────────────────────────────
void cmdGetDescriptor(ScOpenScofo *unit, sc_msg_iter *args) {
    const char *descriptorId = args->gets();
    unit->SendDescriptor(descriptorId);
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

    registerUnit<ScOpenScofo>(ft, "OpenScofoUGen");
    DefineUnitCmd("OpenScofoUGen", "loadScore", (UnitCmdFunc)&cmdLoadScore);
    DefineUnitCmd("OpenScofoUGen", "getCurrentEvent", (UnitCmdFunc)&cmdGetCurrentEvent);
    DefineUnitCmd("OpenScofoUGen", "setFollowScore", (UnitCmdFunc)&cmdSetFollowScore);
    DefineUnitCmd("OpenScofoUGen", "setEventNotifications", (UnitCmdFunc)&cmdSetEventNotifications);
    DefineUnitCmd("OpenScofoUGen", "setNamespace", (UnitCmdFunc)&cmdSetNamespace);
    DefineUnitCmd("OpenScofoUGen", "registerActionReceiver", (UnitCmdFunc)&cmdRegisterActionReceiver);
    DefineUnitCmd("OpenScofoUGen", "unregisterActionReceiver", (UnitCmdFunc)&cmdUnregisterActionReceiver);
    DefineUnitCmd("OpenScofoUGen", "getDescriptor", (UnitCmdFunc)&cmdGetDescriptor);
    DefineUnitCmd("OpenScofoUGen", "loadOnnxModel", (UnitCmdFunc)&cmdLoadOnnxModel);

    Print("\nOpenScofo version %s (%s), by Charles K. Neimog\n\n", OPENSCOFO_VERSION, OPENSCOFO_BUILD_TIME);
}
