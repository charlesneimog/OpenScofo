/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include <OpenScofo.hpp>
#include <cmath>
#include <cstdint>
#include <string>
#include <variant>
#include <vector>

#undef TWOPI

#include <csound.h>

#define opaddr opadr
#include <plugin.h>
#undef opaddr

#include <modload.h>

namespace csnd {

// ─────────────────────────────────────
// Logging callback
// ─────────────────────────────────────
static void oscofo_error_callback(const spdlog::details::log_msg &log, void *data) {
    Csound *csound = static_cast<Csound *>(data);

    if (log.level < spdlog::level::warn)
        return;

    std::string text(log.payload.data(), log.payload.size());

    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
        csound->init_error("[OpenScofo] " + text);
        break;
    case spdlog::level::warn:
        csound->warning("[OpenScofo] " + text);
        break;
    case spdlog::level::info:
    case spdlog::level::debug:
    case spdlog::level::trace:
        csound->message("[OpenScofo] " + text);
        break;
    default:
        break;
    }
}

// ╭─────────────────────────────────────╮
// │               Opcode                │
// ╰─────────────────────────────────────╯
struct CSoundOpenScofo : Plugin<3, 4> {
    struct ScheduledAction {
        int64_t samplesUntilDispatch = 0;
        OpenScofo::ScoreAction action;
    };

    OpenScofo::OpenScofo *oscofo = nullptr;
    float m_SR;
    float m_FFT;
    float m_HOP;
    uint32_t m_BlockIndex = 0;
    int m_LastEvent = 0;
    int m_LastActionStateIndex = -1;
    std::vector<ScheduledAction> m_ScheduledActions;

    // ─────────────────────────────────────
    bool SetControlChannel(const std::string &name, MYFLT value) {
        CSOUND *cs = csound->get_csound();
        MYFLT *channel = nullptr;
        const int result = cs->GetChannelPtr(cs, &channel, name.c_str(), CSOUND_CONTROL_CHANNEL | CSOUND_INPUT_CHANNEL);
        if (result != CSOUND_SUCCESS || channel == nullptr) {
            csound->warning("[OpenScofo] Failed to set Csound control channel '" + name + "'.");
            return false;
        }

        *channel = value;
        return true;
    }

    // ─────────────────────────────────────
    bool ConvertActionArgs(const OpenScofo::ScoreAction &action, std::vector<double> &values) {
        values.clear();
        values.reserve(action.Args.size());

        for (const auto &arg : action.Args) {
            if (std::holds_alternative<float>(arg)) {
                values.push_back(std::get<float>(arg));
            } else if (std::holds_alternative<int>(arg)) {
                values.push_back(static_cast<double>(std::get<int>(arg)));
            } else if (std::holds_alternative<std::string>(arg)) {
                auto msg = "[OpenScofo] Csound sendto receiver '" + action.Receiver +
                           "' was not sent: string argument '" + std::get<std::string>(arg) +
                           "' is not supported. Use numeric arguments only.";
                csound->warning(msg);
                return false;
            }
        }

        return true;
    }

    // ─────────────────────────────────────
    void DispatchAction(const OpenScofo::ScoreAction &action) {
        if (action.isLua) {
            csound->warning("[OpenScofo] Csound wrapper does not execute Lua score actions yet.");
            return;
        }

        if (action.Receiver.empty()) {
            return;
        }

        std::vector<double> values;
        if (!ConvertActionArgs(action, values)) {
            return;
        }

        const std::string base = "OpenScofo/" + action.Receiver;
        SetControlChannel(base + "/count", static_cast<MYFLT>(values.size()));
        for (size_t i = 0; i < values.size(); ++i) {
            SetControlChannel(base + "/" + std::to_string(i), static_cast<MYFLT>(values[i]));
        }
        SetControlChannel(base + "/trigger", static_cast<MYFLT>(1.0));
    }

    // ─────────────────────────────────────
    int64_t ActionDelaySamples(const OpenScofo::ScoreAction &action) {
        double timeMs = action.Time;
        if (!action.AbsoluteTime && oscofo) {
            timeMs = 60.0 / oscofo->GetCurrentBPM() * action.Time * 1000.0;
        }

        if (timeMs <= 0.0 || !oscofo) {
            return 0;
        }

        return static_cast<int64_t>(std::ceil(timeMs * oscofo->GetSr() / 1000.0));
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
        if (!oscofo || !oscofo->ScoreIsLoaded()) {
            return;
        }

        const int currentStateIndex = oscofo->GetCurrentStateIndex();
        if (currentStateIndex == m_LastActionStateIndex) {
            return;
        }

        m_LastActionStateIndex = currentStateIndex;
        OpenScofo::EventActions actions = oscofo->GetCurrentEventActions();
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
    int init() {
        m_SR = static_cast<MYFLT>(csound->sr());
        m_FFT = static_cast<MYFLT>(inargs[2]);
        m_HOP = static_cast<MYFLT>(inargs[3]);

        if (m_SR <= 0 || m_FFT <= 0 || m_HOP <= 0) {
            csound->init_error("FFT or HOP not defined, defined fft = 2048, hop == 512");
            m_SR = 48000;
            m_FFT = 2048;
            m_HOP = 512;
        }

        oscofo = new OpenScofo::OpenScofo(m_SR, m_FFT, m_HOP);
        oscofo->SetErrorCallback(oscofo_error_callback, static_cast<void *>(csound));

        STRINGDAT &scorePath = inargs.str_data(1);
        const char *score = reinterpret_cast<const char *>(scorePath.data);
        if (!oscofo->LoadScore(score)) {
            csound->init_error("Failed to parse score");
            return NOTOK;
        }

        OpenScofo::States states = oscofo->GetStates();
        auto msg = "\nScore has " + std::to_string(states.size()) + " states\n";
        csound->message(msg.c_str());

        outargs[0] = static_cast<MYFLT>(0.0);
        outargs[1] = static_cast<MYFLT>(oscofo->GetCurrentBPM());
        outargs[2] = static_cast<MYFLT>(0.0);
        m_LastEvent = -1;
        m_LastActionStateIndex = -1;

        return OK;
    }

    // ─────────────────────────────────────
    int kperf() {
        csnd::AudioSig in(this, inargs(0));
        MYFLT *audio = in.begin();
        int n = in.GetNsmps();

        outargs[0] = static_cast<MYFLT>(0.0);
        outargs[1] = static_cast<MYFLT>(0.0);
        outargs[2] = static_cast<MYFLT>(0.0);

        if (!oscofo->ProcessBlock(audio, n)) {
            // how to get OPDS?
            csound->warning("[OpenScofo] ProcessBlock failed");
            return NOTOK;
        }

        ProcessEventActionsIfChanged();
        ProcessScheduledActions(n);

        int event = oscofo->GetCurrentScorePosition();
        outargs[0] = static_cast<MYFLT>(event);
        outargs[1] = static_cast<MYFLT>(oscofo->GetCurrentBPM());

        if (event != m_LastEvent) {
            outargs[2] = static_cast<MYFLT>(1.0);
        } else {
            outargs[2] = static_cast<MYFLT>(0.0);
        }
        m_LastEvent = event;
        return OK;
    }

    // ─────────────────────────────────────
    int deinit() {
        if (oscofo) {
            delete oscofo;
            oscofo = nullptr;
        }
        m_ScheduledActions.clear();
        return OK;
    }
};

} // namespace csnd

// Registration
void csnd::on_load(Csound *csound) {
    auto msg = "\n[OpenScofo] version " + std::string(OPENSCOFO_VERSION) + " (" + std::string(OPENSCOFO_BUILD_TIME) +
               "), by Charles K. Neimog\n\n";

    csound->message(msg.c_str());
    csnd::plugin<csnd::CSoundOpenScofo>(csound, "OpenScofoScore", "kkk", "aSii", csnd::thread::ik);
}
