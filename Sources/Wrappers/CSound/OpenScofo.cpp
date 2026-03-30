#include <OpenScofo.hpp>
#undef TWOPI

#define opaddr opadr
#include <plugin.h>
#undef opaddr

#include <modload.h>

namespace csnd {

// ─────────────────────────────────────
// Logging bridge
// ─────────────────────────────────────
static void oscofo_error_callback(const spdlog::details::log_msg &log, void *data) {
    Csound *csound = static_cast<Csound *>(data);

    if (log.level < spdlog::level::warn)
        return;

    std::string text(log.payload.data(), log.payload.size());

    switch (log.level) {
    case spdlog::level::critical:
    case spdlog::level::err:
    case spdlog::level::warn:
        csound->warning(std::format("[OpenScofo] {}", text));
        break;
    case spdlog::level::info:
    case spdlog::level::debug:
    case spdlog::level::trace:
        csound->message(std::format("[OpenScofo] {}", text));
        break;
    default:
        break;
    }
}

// ╭─────────────────────────────────────╮
// │               Opcode                │
// ╰─────────────────────────────────────╯
struct CSoundOpenScofo : Plugin<3, 2> {
    OpenScofo::OpenScofo *oscofo = nullptr;
    float m_SR;
    float m_FFT;
    float m_HOP;
    uint32_t m_BlockIndex = 0;
    int m_LastEvent = 0;

    // ─────────────────────────────────────
    int init() {
        m_SR = sr();
        m_FFT = 2048;
        m_HOP = 512;
        oscofo = new OpenScofo::OpenScofo(48000, 2048, 512);
        oscofo->SetErrorCallback(oscofo_error_callback, static_cast<void *>(csound));

        STRINGDAT &scorePath = inargs.str_data(1);
        const char *score = reinterpret_cast<const char *>(scorePath.data);
        if (!oscofo->ParseScore(score)) {
            csound->warning("Failed to parse score");
            return NOTOK;
        }

        outargs[0] = FL(0.0);
        outargs[1] = oscofo->GetLiveBPM();
        outargs[2] = FL(0.0);
        m_LastEvent = -1;

        return OK;
    }

    // ─────────────────────────────────────
    int kperf() {
        csnd::AudioSig in(this, inargs(0));
        MYFLT *audio = in.begin();
        int n = in.GetNsmps();

        outargs[0] = FL(0.0);
        outargs[1] = FL(0.0);
        outargs[2] = FL(0.0);

        if (!oscofo->ProcessBlock(audio, n)) {
            csound->warning("[OpenScofo] ProcessBlock failed");
            return NOTOK;
        }

        int event = oscofo->GetEventIndex();
        outargs[0] = static_cast<MYFLT>(event);
        outargs[1] = static_cast<MYFLT>(oscofo->GetLiveBPM());

        if (event != m_LastEvent) {
            outargs[2] = FL(1.0);
            csound->message(std::format("Event is {}", event));
        } else {
            outargs[2] = FL(0.0);
        }

        m_LastEvent = event;

        return OK;
    }

    // ─────────────────────────────────────
    // FIX 1: now actually called because plugin_deinit was registered above
    int deinit() {
        if (oscofo) {
            delete oscofo;
            oscofo = nullptr;
        }
        return OK;
    }
};

} // namespace csnd

// Registration
void csnd::on_load(Csound *csound) {
    csound->message(
        std::format("\n[OpenScofo] version {} ({}), by Charles K. Neimog\n\n", OPENSCOFO_VERSION, OSCOFO_BUILD_TIME));
    csnd::plugin<csnd::CSoundOpenScofo>(csound, "OpenScofoScore", "kkk", "aS", csnd::thread::ik);
}
