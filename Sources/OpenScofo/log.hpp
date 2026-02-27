#pragma once

#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/base_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

// ─────────────────────────────────────
template <typename T> inline std::string VectorToString(const std::vector<T> &v) {
    std::ostringstream oss;
    oss << "[";

    for (size_t i = 0; i < v.size(); ++i) {
        oss << std::fixed << std::setprecision(6) << v[i];
        if (i + 1 < v.size())
            oss << ", ";
    }

    oss << "]";
    return oss.str();
}

// ─────────────────────────────────────
template <typename Mutex> class OpenScofoLog : public spdlog::sinks::base_sink<Mutex> {
  private:
    std::function<void(const spdlog::details::log_msg &, void *)> m_Callback;
    void *m_Data = nullptr;
    bool *m_Errors = nullptr;

  public:
    void SetCallback(std::function<void(const spdlog::details::log_msg &, void *data)> cb, void *data, bool *Errors) {
        m_Callback = std::move(cb);
        m_Data = data;
        m_Errors = Errors;
    }

  protected:
    void sink_it_(const spdlog::details::log_msg &msg) override {
        if (!m_Callback)
            return;

        // Set error flag only if pointer is valid
        if (m_Errors && (msg.level == spdlog::level::err || msg.level == spdlog::level::critical)) {
            *m_Errors = true;
        }

        m_Callback(msg, m_Data);
    }

    void flush_() override {
    }
};
