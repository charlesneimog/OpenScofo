#include <gtest/gtest.h>
#include <OpenScofo.hpp>
#include "exception.hpp"

std::filesystem::path assets = TEST_DATA_DIR;

// Test fixture for OpenScofo
class Score : public ::testing::Test {
  protected:
    OpenScofo::Score *score;
    std::shared_ptr<OpenScofoRaise<std::mutex>> sink;

    void SetUp() override {
        sink = std::make_shared<OpenScofoRaise<std::mutex>>();
        auto logger = std::make_shared<spdlog::logger>("openscofo", sink);

        logger->set_error_handler([](const std::string &msg) { throw std::runtime_error(msg); });

        spdlog::set_default_logger(logger);
        spdlog::set_level(spdlog::level::trace);
        spdlog::flush_on(spdlog::level::info);

        score = new OpenScofo::Score(1024.0f, 256.0f);
    }

    void TearDown() override {
        delete score;
    }
};

// 1. Test basic initialization
TEST_F(Score, Initialization) {
    ASSERT_NE(score, nullptr);
}

// WRong trill
TEST_F(Score, WrongTRILL) {
    EXPECT_THROW(score->Parse(assets / "wrongtrill.txt"), std::runtime_error);
}

// 1. Score Parse
TEST_F(Score, ScoreParse) {
    OpenScofo::States states = score->Parse(assets / "bwv-1013.txt");
}

TEST_F(Score, WrongFFTSize) {
    EXPECT_THROW(score->Parse(assets / "wrongfftsize.txt"), std::runtime_error);
}

TEST_F(Score, WrongHopSize) {
    EXPECT_THROW(score->Parse(assets / "wronghopsize.txt"), std::runtime_error);
}
