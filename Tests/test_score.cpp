#include <gtest/gtest.h>

#include <OpenScofo.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace {

const std::filesystem::path Assets = TEST_DATA_DIR;

TEST(ScoreSections, ParsesNamesConfigurationAndPerSectionTiming) {
    OpenScofo::Score Score;
    auto [Config, States] = Score.Parse(Assets / "sections.scofo");

    ASSERT_TRUE(Config.SectionRestrict);
    ASSERT_EQ(States.size(), 9U);

    const std::vector<std::string> ExpectedSections = {"1", "1", "1", "A", "A", "12B2", "12B2", "D", "D"};
    for (size_t Index = 0; Index < States.size(); ++Index) {
        EXPECT_EQ(States[Index].Section, ExpectedSections[Index]);
    }

    EXPECT_DOUBLE_EQ(States[1].BPMExpected, 60.0);
    EXPECT_DOUBLE_EQ(States[4].BPMExpected, 120.0);
    EXPECT_DOUBLE_EQ(States[6].BPMExpected, 90.0);
    EXPECT_DOUBLE_EQ(States[1].OnsetExpected, 0.0);
    EXPECT_DOUBLE_EQ(States[4].OnsetExpected, 0.0);
    EXPECT_DOUBLE_EQ(States[6].OnsetExpected, 0.0);
    EXPECT_EQ(States[7].Type, OpenScofo::FIRSTEVENT);
    EXPECT_DOUBLE_EQ(States[7].BPMExpected, 90.0);
    EXPECT_DOUBLE_EQ(States[8].OnsetExpected, 0.0);
}

TEST(ScoreSections, SelectsSectionAndRestoresItsTemporalConfiguration) {
    OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    ASSERT_TRUE(Scofo.LoadScore(Assets / "sections.scofo"));

    ASSERT_TRUE(Scofo.SetCurrentSection("1"));
    EXPECT_EQ(Scofo.GetStates()[static_cast<size_t>(Scofo.GetCurrentStateIndex())].Type, OpenScofo::FIRSTEVENT);
    EXPECT_EQ(Scofo.GetStates()[static_cast<size_t>(Scofo.GetCurrentStateIndex())].ScorePos, 0);

    ASSERT_TRUE(Scofo.SetCurrentSection("A"));
    EXPECT_EQ(Scofo.GetStates()[static_cast<size_t>(Scofo.GetCurrentStateIndex())].Section, "A");
    EXPECT_EQ(Scofo.GetStates()[static_cast<size_t>(Scofo.GetCurrentStateIndex())].Type, OpenScofo::FIRSTEVENT);
    EXPECT_DOUBLE_EQ(Scofo.GetCurrentBPM(), 120.0);

    ASSERT_TRUE(Scofo.SetCurrentSection("12B2"));
    EXPECT_EQ(Scofo.GetStates()[static_cast<size_t>(Scofo.GetCurrentStateIndex())].Section, "12B2");
    EXPECT_EQ(Scofo.GetStates()[static_cast<size_t>(Scofo.GetCurrentStateIndex())].Type, OpenScofo::FIRSTEVENT);
    EXPECT_DOUBLE_EQ(Scofo.GetCurrentBPM(), 90.0);

    ASSERT_TRUE(Scofo.SetCurrentSection("D"));
    EXPECT_EQ(Scofo.GetStates()[static_cast<size_t>(Scofo.GetCurrentStateIndex())].Type, OpenScofo::FIRSTEVENT);
    EXPECT_DOUBLE_EQ(Scofo.GetCurrentBPM(), 90.0);

    EXPECT_FALSE(Scofo.SetCurrentSection("missing"));
}

TEST(ScoreSections, RestrictsForwardInferenceToSelectedSection) {
    OpenScofo::Configuration Config;
    Config.SectionRestrict = true;

    OpenScofo::States States;
    for (int Index = 0; Index < 4; ++Index) {
        OpenScofo::MarkovState State{};
        State.Index = Index;
        State.ScorePos = Index + 1;
        State.Section = Index < 2 ? "A" : "B";
        State.Type = OpenScofo::REST;
        State.HSMMType = OpenScofo::MARKOV;
        State.Duration = 1.0;
        State.BPMExpected = Index < 2 ? 60.0 : 120.0;
        State.SyncStrength = 0.5;
        State.PhaseCoupling = 0.5;

        OpenScofo::AudioState Silence{};
        Silence.Type = OpenScofo::SILENCE;
        State.AudioStates.push_back(Silence);
        States.push_back(std::move(State));
    }

    OpenScofo::OnlineForward Forward;
    Forward.UpdateConfiguration(Config);
    Forward.SetScoreStates(std::move(States));
    ASSERT_TRUE(Forward.SetCurrentSection("B"));

    OpenScofo::Description Description{};
    Description.SilenceProb = 1.0;
    Forward.GetEvent(Description);

    const auto &ProcessedStates = Forward.GetStates();
    EXPECT_DOUBLE_EQ(ProcessedStates[0].Forward[0], 0.0);
    EXPECT_DOUBLE_EQ(ProcessedStates[1].Forward[0], 0.0);
    EXPECT_GT(ProcessedStates[2].Forward[0], 0.0);
    EXPECT_GT(ProcessedStates[3].Forward[0], 0.0);
}

} // namespace
