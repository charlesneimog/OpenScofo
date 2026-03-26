#include <OpenScofo.hpp>

// TODO: Check how to define this without extra define
#define opaddr opadr
#include <plugin.h>
#undef opaddr

#include <modload.h>

namespace csnd {

struct CSoundOpenScofo : Plugin<4, 3> {
    OpenScofo::OpenScofo *oscofo;
    uint32_t m_FFT;
    uint32_t m_HOP;
    int m_BlockIndex = 0;
    int m_ScoreEventIndex;
    std::vector<double> m_InBuffer;

    int init() {
        double sr = csound->sr();
        m_FFT = 2048;
        m_HOP = 512;
        if (m_FFT < m_HOP) {
            csound->warning(std::format("ksmps (defined as {}) must be less then {}", m_HOP, m_FFT));
            return NOTOK;
        }

        m_ScoreEventIndex = -1;
        oscofo = new OpenScofo::OpenScofo(sr, m_FFT, m_HOP);

        STRINGDAT ScorePath = inargs.str_data(1);
        const char *score = reinterpret_cast<const char *>(ScorePath.data);
        bool ok = oscofo->ParseScore(score);
        if (!ok) {
            csound->warning("Failed to parse the score");
            return NOTOK;
        }

        outargs[0] = FL(0.0);
        outargs[1] = FL(0.0);
        outargs[2] = FL(0.0);

        m_InBuffer.resize(m_FFT);
        return OK;
    }

    int aperf() {
        csnd::AudioSig in(this, inargs(0));
        outargs[0] = FL(0.0);
        outargs[1] = FL(0.0);
        outargs[2] = FL(0.0);

        if (nsmps != m_HOP) {
            csound->warning(std::format("ksmps must be {} for now, received {}", m_HOP, nsmps));
            return NOTOK;
        }

        m_BlockIndex += m_HOP;
        std::memmove(m_InBuffer.data(), m_InBuffer.data() + m_HOP, (m_InBuffer.size() - m_HOP) * sizeof(double));
        for (uint32_t i = 0; i < m_HOP; ++i) {
            m_InBuffer[m_InBuffer.size() - m_HOP + i] = in[i];
        }

        bool ok = oscofo->ProcessBlock(m_InBuffer);
        if (!ok) {
            csound->message("Failed to process block");
            return NOTOK;
        }
        int event = oscofo->GetEventIndex();
        outargs[0] = static_cast<MYFLT>(event);
        outargs[1] = static_cast<MYFLT>(oscofo->GetLiveBPM());
        if (event >= 0 && m_ScoreEventIndex != event) {
            m_ScoreEventIndex = event;
            outargs[2] = FL(1.0);
        } else {
            outargs[2] = FL(0.0);
        }
        return OK;
    }

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
    csnd::plugin<CSoundOpenScofo>(csound, "OpenScofoScore", "kkk", "aS", csnd::thread::ia);
}
