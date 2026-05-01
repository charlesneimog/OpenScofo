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
struct CSoundOpenScofo : Plugin<3, 4> {
    OpenScofo::OpenScofo *oscofo = nullptr;
    float m_SR;
    float m_FFT;
    float m_HOP;
    uint32_t m_BlockIndex = 0;
    int m_LastEvent = 0;

    // ─────────────────────────────────────
    int init() {
        m_SR = static_cast<MYFLT>(csound->sr());
        m_FFT = static_cast<MYFLT>(inargs[2]);
        m_HOP = static_cast<MYFLT>(inargs[3]);

        if (m_SR <= 0 || m_FFT <= 0 || m_HOP <= 0) {
            csound->warning("FFT or HOP not defined, defined fft = 2048, hop == 512");
            m_SR = 48000;
            m_FFT = 2048;
            m_HOP = 512;
        }

        oscofo = new OpenScofo::OpenScofo(m_SR, m_FFT, m_HOP);
        oscofo->SetErrorCallback(oscofo_error_callback, static_cast<void *>(csound));

        STRINGDAT &scorePath = inargs.str_data(1);
        const char *score = reinterpret_cast<const char *>(scorePath.data);
        if (!oscofo->LoadScore(score)) {
            csound->warning("Failed to parse score");
            return NOTOK;
        }

        OpenScofo::States states = oscofo->GetStates();
        csound->message(std::format("\nScore has {} states\n", states.size()));

        outargs[0] = FL(0.0);
        outargs[1] = FL(oscofo->GetLiveBPM());
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

        int event = oscofo->GetCurrentScorePosition();
        outargs[0] = static_cast<MYFLT>(event);
        outargs[1] = static_cast<MYFLT>(oscofo->GetLiveBPM());

        if (event != m_LastEvent) {
            outargs[2] = FL(1.0);
        } else {
            outargs[2] = FL(0.0);
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
        return OK;
    }
};

} // namespace csnd

// Registration
void csnd::on_load(Csound *csound) {
    csound->message(std::format("\n[OpenScofo] version {} ({}), by Charles K. Neimog\n\n", OPENSCOFO_VERSION,
                                OPENSCOFO_BUILD_TIME));
    csnd::plugin<csnd::CSoundOpenScofo>(csound, "OpenScofoScore", "kkk", "aSii", csnd::thread::ik);
}
