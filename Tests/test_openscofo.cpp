#include <gtest/gtest.h>
#include <OpenScofo.hpp>
#include "exception.hpp"

// Test fixture for OpenScofo
class OpenScofoTest : public ::testing::Test {
  protected:
    std::shared_ptr<OpenScofoRaise<std::mutex>> sink;
    OpenScofo::OpenScofo *oscofo;

    void SetUp() override {
        sink = std::make_shared<OpenScofoRaise<std::mutex>>();
        auto logger = std::make_shared<spdlog::logger>("openscofo", sink);
        spdlog::set_default_logger(logger);
        spdlog::set_level(spdlog::level::trace); // or desired level
        spdlog::flush_on(spdlog::level::info);

        // Initialize with typical Sample Rate, FFT Size, and Hop Size
        oscofo = new OpenScofo::OpenScofo(44100.0f, 1024.0f, 256.0f);
    }

    void TearDown() override {
        delete oscofo;
    }
};

// 1. Test basic initialization
TEST_F(OpenScofoTest, Initialization) {
    ASSERT_NE(oscofo, nullptr);
}

// 2. Test default parameters using the Description struct
TEST_F(OpenScofoTest, DefaultDescription) {
    OpenScofo::Description desc = oscofo->GetDescription();
    EXPECT_FALSE(desc.Onset);
}

// 3. Test Audio Processing with an empty buffer
