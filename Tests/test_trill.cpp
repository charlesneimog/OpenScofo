#include <gtest/gtest.h>
#include <OpenScofo.hpp>

#include <algorithm>
#include <filesystem>
#include <limits>
#include <numeric>

namespace {

const std::filesystem::path TrillScore = std::filesystem::path(TEST_DATA_DIR) / "trill.scofo";

TEST(TrillMicroStates, ParsesOnlyTrillPitchesAsMicroStates) {
    OpenScofo::Score Parser;
    auto [Config, States] = Parser.Parse(TrillScore);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ASSERT_EQ(States.size(), 6U);
    EXPECT_EQ(States[0].Type, OpenScofo::FIRSTEVENT);
    EXPECT_EQ(States[1].Type, OpenScofo::NOTE);
    EXPECT_EQ(States[2].Type, OpenScofo::CHORD);
    EXPECT_EQ(States[3].Type, OpenScofo::TRILL);
    EXPECT_EQ(States[4].Type, OpenScofo::REST);
    for (size_t Index : {0U, 1U, 2U, 4U, 5U}) {
        EXPECT_FALSE(States[Index].Observations.empty());
        EXPECT_TRUE(States[Index].MicroStates.empty());
        EXPECT_EQ(States[Index].MicroTopologyType, OpenScofo::NO_MICROSTATES);
    }
    ASSERT_EQ(States[1].Observations.size(), 1U);
    EXPECT_DOUBLE_EQ(States[1].Observations[0].Midi, 60.0);
    ASSERT_EQ(States[2].Observations.size(), 3U);
    EXPECT_DOUBLE_EQ(States[2].Observations[0].Midi, 60.0);
    EXPECT_DOUBLE_EQ(States[2].Observations[1].Midi, 64.0);
    EXPECT_DOUBLE_EQ(States[2].Observations[2].Midi, 67.0);

    const auto &Trill = States[3];
    EXPECT_EQ(Trill.HSMMType, OpenScofo::SEMIMARKOV);
    EXPECT_DOUBLE_EQ(Trill.Duration, 4.0);
    EXPECT_EQ(Trill.MicroTopologyType, OpenScofo::UNORDERED);
    EXPECT_TRUE(Trill.Observations.empty());
    ASSERT_EQ(Trill.MicroStates.size(), 2U);
    for (size_t Index = 0; Index < Trill.MicroStates.size(); ++Index) {
        ASSERT_EQ(Trill.MicroStates[Index].Observations.size(), 1U);
        EXPECT_EQ(Trill.MicroStates[Index].Observations[0].Type, OpenScofo::PITCH);
        EXPECT_DOUBLE_EQ(Trill.MicroStates[Index].Observations[0].Midi, 60.0 + 2.0 * Index);
    }
}

OpenScofo::Description PitchDescription(OpenScofo::OnlineForward &Forward, double Frequency) {
    OpenScofo::Description Desc{};
    Desc.SpectralMagnitudeFrameNorm = Forward.GetPitchTemplate(Frequency);
    // A small noise floor keeps both pitch likelihoods finite and distinct.
    for (double &Bin : Desc.SpectralMagnitudeFrameNorm) {
        Bin += 1e-8;
    }
    const double Sum =
        std::accumulate(Desc.SpectralMagnitudeFrameNorm.begin(), Desc.SpectralMagnitudeFrameNorm.end(), 0.0);
    for (double &Bin : Desc.SpectralMagnitudeFrameNorm) {
        Bin /= Sum;
    }
    return Desc;
}

TEST(TrillMicroStates, RejectsEmptyTrill) {
    OpenScofo::Score Parser;
    const auto [Config, States] = Parser.Parse(TrillScore.parent_path() / "trill-invalid.scofo");
    for (const auto &State : States) {
        EXPECT_EQ(State.Type, OpenScofo::FIRSTEVENT);
        EXPECT_TRUE(State.MicroStates.empty());
    }
}

TEST(TrillMicroStates, PreservesMaxEmissionAndOrdinaryEventFormulas) {
    OpenScofo::Score Parser;
    auto [Config, States] = Parser.Parse(TrillScore);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ASSERT_EQ(States.size(), 6U);
    const double C4 = States[3].MicroStates[0].Observations[0].Freq;
    const double D4 = States[3].MicroStates[1].Observations[0].Freq;
    for (size_t Index : {1U, 2U, 3U, 4U, 5U}) {
        OpenScofo::OnlineForward Forward;
        Forward.UpdateConfiguration(Config);
        Forward.SetScoreStates({States[Index]});
        for (double Frequency : {C4, D4}) {
            for (double Silence : {0.0, 0.25, 1.0}) {
                for (double Technique : {0.0, 0.95}) {
                    auto Desc = PitchDescription(Forward, Frequency);
                    Desc.SilenceProb = Silence;
                    Desc.ExtendedTechProb = Technique;
                    Forward.SetDescription(Desc);
                    const double P0 = Forward.GetPitchProbability(C4);
                    const double P1 = Forward.GetPitchProbability(D4);
                    double Expected = 0.0;
                    if (Index == 3) {
                        Expected = std::max(P0, P1) * (1.0 - Silence);
                    } else if (Index == 2) {
                        for (const auto &Obs : States[Index].Observations) {
                            Expected += Forward.GetPitchProbability(Obs.Freq);
                        }
                        Expected = Expected / 3.0 * (1.0 - Silence);
                    } else if (Index == 4) {
                        Expected = Silence;
                    } else {
                        Expected = P0 * (1.0 - Silence);
                    }
                    const int Buffer = Forward.GetCurrentBufferIndex();
                    Forward.GetEvent(Desc);
                    EXPECT_DOUBLE_EQ(Forward.GetStates()[0].BestObs[Buffer],
                                     std::max(Expected, std::numeric_limits<double>::min()));
                    if (Index == 3) {
                        EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, P1 > P0 ? 1 : 0);
                        EXPECT_EQ(Forward.GetStates()[0].BestMicroObservationIndex, 0);
                    }
                }
            }
        }
    }
}

TEST(TrillMicroStates, ReportsPitchChangesAndResetsWinner) {
    OpenScofo::Score Parser;
    auto [Config, States] = Parser.Parse(TrillScore);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ASSERT_EQ(States.size(), 6U);
    OpenScofo::OnlineForward Forward;
    Forward.UpdateConfiguration(Config);
    Forward.SetScoreStates({States[3]});
    for (size_t Index : {0U, 1U, 0U}) {
        const double Frequency = States[3].MicroStates[Index].Observations[0].Freq;
        auto Desc = PitchDescription(Forward, Frequency);
        Forward.GetEvent(Desc);
        EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, static_cast<int>(Index));
        auto Actions = Forward.GetAudioStateChangeActions();
        ASSERT_EQ(Actions.size(), 1U);
        EXPECT_TRUE(Actions[0].isAudioStateChange);
        EXPECT_EQ(Actions[0].Receiver, "trill_pitch");
        ASSERT_EQ(Actions[0].Args.size(), 2U);
        EXPECT_EQ(std::get<int>(Actions[0].Args[0]), States[3].ScorePos);
        EXPECT_FLOAT_EQ(std::get<float>(Actions[0].Args[1]), static_cast<float>(Frequency));
        Forward.GetEvent(Desc);
        EXPECT_TRUE(Forward.GetAudioStateChangeActions().empty());
    }
    Forward.ResetDecoding();
    EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, -1);
    EXPECT_EQ(Forward.GetStates()[0].BestMicroObservationIndex, -1);
    auto Desc = PitchDescription(Forward, States[3].MicroStates[0].Observations[0].Freq);
    Forward.GetEvent(Desc);
    EXPECT_EQ(Forward.GetAudioStateChangeActions().size(), 1U);
}

TEST(ObservationNotifications, TracksEachNoteWinnerAndReportsPitch) {
    OpenScofo::Score Parser;
    auto [Config, States] = Parser.Parse(TrillScore);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ASSERT_EQ(States.size(), 6U);
    auto OtherNote = States[1];
    OtherNote.Observations[0] = States[3].MicroStates[1].Observations[0];
    OtherNote.Index = 1;
    OtherNote.ScorePos = 2;
    OpenScofo::OnlineForward Forward;
    Forward.UpdateConfiguration(Config);
    Forward.SetScoreStates({States[1], OtherNote});
    const double Frequency = States[1].Observations[0].Freq;
    auto Desc = PitchDescription(Forward, Frequency);
    Forward.GetEvent(Desc);
    // The weaker candidate must have its own winner, too.
    EXPECT_EQ(Forward.GetStates()[0].BestAudioStateIndex, 0);
    EXPECT_EQ(Forward.GetStates()[1].BestAudioStateIndex, 0);
    auto Actions = Forward.GetAudioStateChangeActions();
    ASSERT_EQ(Actions.size(), 1U);
    EXPECT_TRUE(Actions[0].isAudioStateChange);
    EXPECT_EQ(Actions[0].Receiver, Config.AudioStateChangeReceiver);
    ASSERT_EQ(Actions[0].Args.size(), 2U);
    EXPECT_EQ(std::get<int>(Actions[0].Args[0]), States[1].ScorePos);
    EXPECT_FLOAT_EQ(std::get<float>(Actions[0].Args[1]), static_cast<float>(Frequency));
    Forward.GetEvent(Desc);
    EXPECT_TRUE(Forward.GetAudioStateChangeActions().empty());
    Forward.ResetDecoding();
    EXPECT_EQ(Forward.GetStates()[0].BestAudioStateIndex, -1);
    EXPECT_EQ(Forward.GetStates()[1].BestAudioStateIndex, -1);
    Forward.GetEvent(Desc);
    EXPECT_EQ(Forward.GetAudioStateChangeActions().size(), 1U);
}

TEST(ObservationNotifications, ReportsDirectWinnerChangesAndClearsZeroEvidence) {
    OpenScofo::Score Parser;
    auto [Config, States] = Parser.Parse(TrillScore);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ASSERT_EQ(States.size(), 6U);
    auto State = States[1];
    State.Type = OpenScofo::UTECH;
    State.Observations = {{OpenScofo::LABEL, 0, 0, "first"},
                          {OpenScofo::LABEL, 0, 0, "second"},
                          {OpenScofo::ONSET},
                          {OpenScofo::SILENCE}};
    State.BestAudioStateIndex = 3;
    OpenScofo::OnlineForward Forward;
    Forward.UpdateConfiguration(Config);
    Forward.SetScoreStates({State});
    EXPECT_EQ(Forward.GetStates()[0].BestAudioStateIndex, -1);
    const std::vector<std::string> Expected = {"first", "second", "onset", "silence"};
    for (int Index = 0; Index < 4; ++Index) {
        OpenScofo::Description Desc{};
        Desc.ExtendedTechProb = 1.0;
        Desc.ONNX["first"] = Index == 0 ? 0.8f : 0.1f;
        Desc.ONNX["second"] = Index == 1 ? 0.9f : 0.1f;
        Desc.Onset = Index == 2 ? 1.0 : 0.0;
        Desc.SilenceProb = Index == 3 ? 0.9 : 0.0;
        Forward.GetEvent(Desc);
        EXPECT_EQ(Forward.GetStates()[0].BestAudioStateIndex, Index);
        auto Actions = Forward.GetAudioStateChangeActions();
        ASSERT_EQ(Actions.size(), 1U);
        ASSERT_EQ(Actions[0].Args.size(), 2U);
        EXPECT_EQ(std::get<std::string>(Actions[0].Args[1]), Expected[Index]);
        Forward.GetEvent(Desc);
        EXPECT_TRUE(Forward.GetAudioStateChangeActions().empty());
    }
    OpenScofo::Description NoEvidence{};
    Forward.GetEvent(NoEvidence);
    EXPECT_EQ(Forward.GetStates()[0].BestAudioStateIndex, -1);
    EXPECT_TRUE(Forward.GetAudioStateChangeActions().empty());
}

#if defined(OPENSCOFO_LUA)
TEST(TrillMicroStates, ExportsNestedObservationsToLua) {
    OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    ASSERT_TRUE(Scofo.LoadScore(TrillScore));
    EXPECT_TRUE(Scofo.LuaExecute(R"(
        local all = require('OpenScofo').get_states()
        local states = {}
        for _, state in ipairs(all) do
            if not state.inter_event_silence then table.insert(states, state) end
        end
        assert(#states[2].audiostates == 1 and #states[2].microstates == 0)
        assert(#states[3].audiostates == 3 and #states[3].microstates == 0)
        local trill = states[4]
        assert(#trill.audiostates == 0 and #trill.microstates == 2)
        assert(trill.micro_topology == 1)
        assert(trill.microstates[1].observations[1].midi == 60)
        assert(trill.microstates[2].observations[1].midi == 62)
    )")) << Scofo.LuaGetError();
}
#endif

} // namespace
