#include <OpenScofo.hpp>
#undef TWOPI

// TODO: Check how to define this without extra define
#define opaddr opadr
#include <plugin.h>
#undef opaddr

#include <modload.h>

namespace csnd {

// ─────────────────────────────────────
static void oscofo_error_callback(const spdlog::details::log_msg &log, void *data) {
    Csound *csound = static_cast<Csound *>(data);
    spdlog::level::level_enum pdlevel = spdlog::level::warn;
    if (log.level < pdlevel) {
        return;
    }

    std::string text(log.payload.data(), log.payload.size());
    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
        csound->warning(std::format("[OpenScofo] {}", text));
        break;
    case spdlog::level::info:
    case spdlog::level::warn:
        csound->warning(std::format("[OpenScofo] {}", text));
        break;
    case spdlog::level::debug:
    case spdlog::level::trace:
        csound->message(std::format("[OpenScofo] {}", text));
        break;
    default:
        break;
    }
}

// ─────────────────────────────────────
struct CSoundOpenScofo : Plugin<3, 2> {
    OpenScofo::OpenScofo *oscofo;
    uint32_t m_FFT;
    uint32_t m_HOP;
    uint32_t m_BlockIndex = 0;
    uint32_t m_ScoreEventIndex;
    std::vector<double> m_InBuffer;
    MYFLT m_CurrentEvent = FL(0.0);
    MYFLT m_CurrentBPM = FL(0.0);
    bool m_Following = false; // Add a flag to prevent the -nan silence issue

    int init() {
        double sr = csound->sr();
        m_FFT = 2048;
        m_HOP = 512;
        m_BlockIndex = 0;
        m_ScoreEventIndex = 0;
        if (m_FFT < m_HOP) {
            csound->warning(std::format("ksmps (defined as {}) must be less then {}", m_HOP, m_FFT));
            return NOTOK;
        }

        m_ScoreEventIndex = -1;
        oscofo = new OpenScofo::OpenScofo(sr, m_FFT, m_HOP);
        oscofo->SetErrorCallback(oscofo_error_callback, static_cast<void *>(csound));

        STRINGDAT ScorePath = inargs.str_data(1);
        const char *score = reinterpret_cast<const char *>(ScorePath.data);
        bool ok = oscofo->ParseScore(score);
        if (!ok) {
            csound->warning("Failed to parse the score");
            return NOTOK;
        }

        OpenScofo::States states = oscofo->GetStates();
        csound->message(std::format("There is {} states", states.size()));

        outargs[0] = FL(0.0);
        outargs[1] = FL(0.0);
        outargs[2] = FL(0.0);

        m_InBuffer.resize(m_FFT);
        return OK;
    }

    // ─────────────────────────────────────
    int aperf() {
        csnd::AudioSig in(this, &inargs[0]);
        int n = in.GetNsmps();

        // 1. Only reset the trigger output, not the state
        outargs[2] = FL(0.0);

        // [ Keep your existing buffer copy logic here ]

        if (m_BlockIndex >= m_HOP) {
            // Optional: Only process if you've received a trigger/start command
            // if (m_Following) {
            bool ok = oscofo->ProcessBlock(m_InBuffer);
            // ... error handling ...

            int event = oscofo->GetEventIndex();
            m_CurrentEvent = static_cast<MYFLT>(event);
            m_CurrentBPM = static_cast<MYFLT>(oscofo->GetLiveBPM());

            if (event >= 0 && m_ScoreEventIndex != static_cast<uint32_t>(event)) {
                m_ScoreEventIndex = event;
                outargs[2] = FL(1.0); // Trigger high for exactly one k-cycle
            }
            // }

            // 2. Subtract instead of forcing to 0 to preserve remainder samples
            m_BlockIndex -= m_HOP;
        }

        // 3. Constantly output the held state
        outargs[0] = m_CurrentEvent;
        outargs[1] = m_CurrentBPM;

        return OK;
    }

    // ─────────────────────────────────────
    int deinit() {
        if (oscofo) {
            delete oscofo;
            oscofo = nullptr;
        }
        return OK;
    }
};
} // namespace csnd

void csnd::on_load(Csound *csound) {
    csound->message(
        std::format("\n[OpenScofo] version {} ({}), by Charles K. Neimog\n\n", OPENSCOFO_VERSION, OSCOFO_BUILD_TIME));
    csnd::plugin<CSoundOpenScofo>(csound, "OpenScofoScore", "kkk", "aS", csnd::thread::ia);
}
