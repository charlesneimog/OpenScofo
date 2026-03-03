#include <OpenScofo.hpp>

// ─────────────────────────────────────
template <typename Mutex> class OpenScofoRaise : public spdlog::sinks::base_sink<Mutex> {
  protected:
    void sink_it_(const spdlog::details::log_msg &msg) override {
        std::string text(msg.payload.data(), msg.payload.size());
        if ((msg.level == spdlog::level::err || msg.level == spdlog::level::critical)) {
            throw std::runtime_error(text);
        }
    }

    void flush_() override {
        fflush(stdout); // garante que tudo seja impresso
    }
};
