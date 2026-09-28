#include <gtest/gtest.h>
#include <OpenScofo.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numeric>

namespace OpenScofo {

// Keep deterministic inference tests independent of FFT/ONNX and public APIs.
struct OnlineForwardTestAccess {
    static void Setup(OnlineForward &Forward, States ScoreStates) {
        Configuration Config;
        Config.AudioStateChangeReceiver = "microstate";
        Forward.UpdateConfiguration(Config);
        Forward.SetScoreStates(std::move(ScoreStates));
        Forward.m_BlockDur = 1.0;
        Forward.m_PsiN1 = 1.0;
    }

    static void Frame(OnlineForward &Forward, int Time) {
        Forward.m_Tau = Time;
        Forward.m_CircularBufferIndex = Time % Forward.m_BufferSize;
        Forward.m_WinStart = 0;
        Forward.m_WinEnd = static_cast<int>(Forward.m_States.size()) - 1;
    }

    static double Transition(OnlineForward &Forward, int StateIndex, int I, int J) {
        Forward.m_ActiveMarkovScoreStateIndex = StateIndex;
        Forward.PrepareMicroStateDurations(Forward.m_States[StateIndex]);
        const double Result = Forward.GetMarkovTransitionProbability(I, J);
        Forward.m_ActiveMarkovScoreStateIndex = -1;
        return Result;
    }

    static void Emissions(OnlineForward &Forward, int StateIndex, const std::vector<double> &Values) {
        auto &State = Forward.m_States[StateIndex];
        for (size_t K = 0; K < Values.size(); ++K) {
            State.MicroStates[K].CurrentEmission = Values[K];
            State.MicroStates[K].BestObservationIndex = 0;
        }
    }

    static void Advance(OnlineForward &Forward, int Time, int MaxAge, const std::vector<double> &Values) {
        Frame(Forward, Time);
        Emissions(Forward, 0, Values);
        Forward.UpdateMicroStateForward(Forward.m_States[0], 0, MaxAge);
        EXPECT_EQ(Forward.m_ActiveMarkovScoreStateIndex, -1);
    }

    static void SemiMarkov(OnlineForward &Forward, int Index) {
        Forward.SemiMarkov(Forward.m_States[Index], Index);
    }

    static double Normalize(OnlineForward &Forward) {
        Forward.GetAlphaT();
        return Forward.m_Normalization[Forward.m_CircularBufferIndex];
    }

    static void Observe(OnlineForward &Forward, const Description &Desc) {
        Forward.SetDescription(Desc);
        Forward.GetAudioObservations();
    }

    static double Emission(OnlineForward &Forward, MarkovMicroState &Micro, const Description &Desc) {
        Forward.SetDescription(Desc);
        return Forward.GetMicroStateEmission(Micro, true, true);
    }

    static double Occupancy(OnlineForward &Forward, int Index, int Age, bool Survivor) {
        return Survivor ? Forward.GetSurvivorDistribution(Forward.m_States[Index], Age)
                        : Forward.GetOccupancyDistribution(Forward.m_States[Index], Age);
    }

    static void Notify(OnlineForward &Forward, int Index) {
        Forward.NotifyAudioStateChange(Index);
    }

    static double ScoreTransition(OnlineForward &Forward, int I, int J) {
        return Forward.GetSemiMarkovTransitionProbability(I, J);
    }

    static void Tempo(OnlineForward &Forward, double Period) {
        Forward.m_PsiN1 = Period;
    }
};

} // namespace OpenScofo

namespace {

using namespace OpenScofo;
using Access = OnlineForwardTestAccess;
const std::filesystem::path Asset = std::filesystem::path(TEST_DATA_DIR) / "microstates.scofo";
const double LogZero = -std::numeric_limits<double>::max();

ScoreState Chain(EventType Type, size_t Count, double Duration) {
    ScoreState State{};
    State.Type = Type;
    State.HSMMType = SEMIMARKOV;
    State.MicroTopologyType = LEFT_RIGHT;
    State.Duration = Duration;
    State.BPMExpected = 60;
    State.ScorePos = 1;
    State.MicroStates.resize(Count);
    for (auto &Micro : State.MicroStates) {
        Micro.Observations.push_back({PITCH, 440.0});
    }
    return State;
}

double Alpha(const MarkovMicroState &Micro, int Age) {
    return Micro.LogForwardByAge[Age] == LogZero ? 0.0 : std::exp(Micro.LogForwardByAge[Age]);
}

// Independent exhaustive enumeration of all allowed paths (small tests only).
void EnumeratePaths(const std::vector<std::vector<double>> &Emissions, const std::vector<double> &Advance, int Time,
                    int End, size_t K, double Mass, std::vector<double> &Result) {
    Mass *= Emissions[Time][K];
    if (Time == End) {
        Result[K] += Mass;
        return;
    }
    EnumeratePaths(Emissions, Advance, Time + 1, End, K, Mass * (1.0 - Advance[K]), Result);
    if (K + 1 < Advance.size()) {
        EnumeratePaths(Emissions, Advance, Time + 1, End, K + 1, Mass * Advance[K], Result);
    }
}

std::vector<double> Segment(const std::vector<std::vector<double>> &Emissions, const std::vector<double> &Advance,
                            int Start, int End) {
    std::vector<double> Result(Advance.size(), 0.0);
    EnumeratePaths(Emissions, Advance, Start, End, 0, 1.0, Result);
    return Result;
}

TEST(InterEventSilence, ParserSeparatesUnscoredGapsFromSoundAndScoredRests) {
    Score Parser;
    const auto [Config, States] = Parser.Parse(Asset.parent_path() / "markov-silence.scofo");
    ASSERT_EQ(States.size(), 14U); // Eight events, two section starts, four gaps.
    int Gaps = 0, Events = 0, ScorePosition = 0;
    double Beats = 0.0;
    for (size_t I = 0; I < States.size(); ++I) {
        const auto &State = States[I];
        EXPECT_EQ(State.Index, I);
        if (State.IsInterEventSilence) {
            ++Gaps;
            ASSERT_GT(I, 0U);
            ASSERT_LT(I + 1, States.size());
            EXPECT_EQ(State.HSMMType, MARKOV);
            EXPECT_EQ(State.Type, REST);
            EXPECT_DOUBLE_EQ(State.Duration, 0.0);
            EXPECT_EQ(State.ScorePos, States[I - 1].ScorePos);
            EXPECT_EQ(State.Section, States[I + 1].Section);
            EXPECT_DOUBLE_EQ(State.OnsetExpected, States[I + 1].OnsetExpected);
            EXPECT_TRUE(State.Actions.empty());
            ASSERT_EQ(State.Observations.size(), 1U);
            EXPECT_EQ(State.Observations[0].Type, SILENCE);
        } else if (State.Type != FIRSTEVENT) {
            ++Events;
            if (State.Type != REST) ++ScorePosition;
            EXPECT_EQ(State.ScorePos, ScorePosition);
            Beats += State.Duration;
            if (State.Type == REST) {
                EXPECT_EQ(State.HSMMType, SEMIMARKOV);
                EXPECT_DOUBLE_EQ(State.Duration, 2.0);
            } else {
                for (const auto &Obs : State.Observations) EXPECT_NE(Obs.Type, SILENCE);
                for (const auto &Micro : State.MicroStates)
                    for (const auto &Obs : Micro.Observations) EXPECT_NE(Obs.Type, SILENCE);
            }
        }
    }
    EXPECT_EQ(Gaps, 4);
    EXPECT_EQ(Events, 8);
    EXPECT_DOUBLE_EQ(Beats, 10.0);
}

States GapChain(bool Nested = false) {
    ScoreState Sound = Chain(UTECH, 0, 2.0);
    Sound.MicroTopologyType = NO_MICROSTATES;
    Sound.Observations = {{LABEL, 0, 0, "first"}};
    ScoreState Gap{};
    Gap.Type = REST;
    Gap.HSMMType = MARKOV;
    Gap.IsInterEventSilence = true;
    Gap.ScorePos = 1;
    Gap.BPMExpected = 60;
    Gap.Observations = {{SILENCE}};
    ScoreState Next = Sound;
    Next.ScorePos = 2;
    Next.Observations = {{LABEL, 0, 0, "second"}};
    if (Nested) {
        Next.MicroTopologyType = LEFT_RIGHT;
        Next.MicroStates.resize(1);
        Next.MicroStates[0].Observations = Next.Observations;
        Next.Observations.clear();
    }
    States Result = {Sound, Gap, Next};
    for (size_t I = 0; I < Result.size(); ++I) Result[I].Index = I;
    return Result;
}

TEST(InterEventSilence, NormalizedBranchesAndGeometricOccupancy) {
    OnlineForward Forward;
    Access::Setup(Forward, GapChain());
    EXPECT_DOUBLE_EQ(Access::ScoreTransition(Forward, 0, 1), 0.5);
    EXPECT_DOUBLE_EQ(Access::ScoreTransition(Forward, 0, 2), 0.5);
    EXPECT_DOUBLE_EQ(Access::ScoreTransition(Forward, 1, 2), 1.0);
    EXPECT_DOUBLE_EQ(Access::ScoreTransition(Forward, 2, 0), 0.0);
    Forward.GetStates()[1].InitProb = 1.0;
    double Evidence = 1.0;
    for (int T = 0; T < 12; ++T) {
        Access::Frame(Forward, T);
        for (auto &State : Forward.GetStates()) State.BestObs[T] = 0.0;
        Forward.GetStates()[1].BestObs[T] = 1.0;
        Evidence *= Access::Normalize(Forward);
        EXPECT_NEAR(Evidence, std::pow(0.5, T), 1e-14);
        EXPECT_NEAR(Forward.GetStates()[1].Forward[T], 1.0, 1e-14);
        EXPECT_NEAR(Forward.GetStates()[1].ExitProb[T], 0.5, 1e-14);
    }
}

TEST(InterEventSilence, BothSemiMarkovRecursionsAcceptTheDirectBypass) {
    for (bool Nested : {false, true}) {
        OnlineForward Forward;
        Access::Setup(Forward, GapChain(Nested));
        Access::Frame(Forward, 1);
        Forward.GetStates()[0].ExitProb[0] = 0.8;
        Forward.GetStates()[1].ExitProb[0] = 0.0;
        Forward.GetStates()[2].BestObs[1] = 0.7;
        if (Nested) Access::Emissions(Forward, 2, {0.7});
        Access::SemiMarkov(Forward, 2);
        EXPECT_NEAR(Forward.GetStates()[2].Forward[1],
                    Access::Occupancy(Forward, 2, 1, true) * 0.8 * 0.5 * 0.7, 1e-14);
    }
}

TEST(InterEventSilence, DecodesWithAndWithoutPauseAndDoesNotAdvanceOnSilence) {
    for (int Pause : {0, 12}) {
        OnlineForward Forward;
        Access::Setup(Forward, GapChain());
        Description Desc{};
        Desc.ExtendedTechProb = 1.0;
        Desc.ONNX = {{"first", 1.0f}};
        EXPECT_EQ(Forward.GetEvent(Desc), 1);
        EXPECT_EQ(Forward.GetEvent(Desc), 1);
        Desc.ONNX.clear();
        Desc.SilenceProb = 1.0;
        const double BPM = Forward.GetCurrentBPM();
        for (int T = 0; T < Pause; ++T) {
            EXPECT_EQ(Forward.GetEvent(Desc), 1);
            EXPECT_EQ(Forward.GetCurrentStateIndex(), 1);
            EXPECT_DOUBLE_EQ(Forward.GetCurrentBPM(), BPM);
            EXPECT_TRUE(Forward.GetCurrentEventActions().empty());
        }
        Desc.SilenceProb = 0.0;
        Desc.ONNX = {{"second", 1.0f}};
        EXPECT_EQ(Forward.GetEvent(Desc), 2);
        EXPECT_EQ(Forward.GetCurrentStateIndex(), 2);
        EXPECT_TRUE(std::isfinite(Forward.GetCurrentBPM()));
    }
}

TEST(MicroStateParsing, PreservesPitchOrderAndGroupsTechniqueLabels) {
    Score Parser;
    auto [Config, States] = Parser.Parse(Asset);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ASSERT_EQ(States.size(), 6U);
    const auto &Multi = States[1];
    EXPECT_EQ(Multi.Type, MULTI);
    EXPECT_EQ(Multi.HSMMType, SEMIMARKOV);
    EXPECT_EQ(Multi.MicroTopologyType, LEFT_RIGHT);
    EXPECT_TRUE(Multi.Observations.empty());
    EXPECT_DOUBLE_EQ(Multi.Duration, 2.0);
    ASSERT_EQ(Multi.MicroStates.size(), 9U);
    for (size_t K = 0; K < 9; ++K) {
        ASSERT_EQ(Multi.MicroStates[K].Observations.size(), 1U);
        EXPECT_EQ(Multi.MicroStates[K].Observations[0].Type, PITCH);
        EXPECT_DOUBLE_EQ(Multi.MicroStates[K].Observations[0].Midi, 60.0 + 0.5 * K);
        EXPECT_DOUBLE_EQ(Multi.MicroStates[K].DurationWeight, 1.0);
    }
    for (size_t Index : {2U, 3U, 4U, 5U}) {
        const auto &State = States[Index];
        const bool Pitched = Index == 2 || Index == 4;
        EXPECT_EQ(State.Type, Pitched ? PTECH : UTECH);
        EXPECT_EQ(State.HSMMType, SEMIMARKOV);
        EXPECT_EQ(State.MicroTopologyType, UNORDERED);
        EXPECT_TRUE(State.Observations.empty());
        ASSERT_EQ(State.MicroStates.size(), Pitched ? 2U : 1U);
        ASSERT_EQ(State.MicroStates[0].Observations.size(), Index < 4 ? 2U : 1U);
        for (const auto &Obs : State.MicroStates[0].Observations) {
            EXPECT_EQ(Obs.Type, LABEL);
            EXPECT_FALSE(Obs.Label.empty());
        }
        if (Pitched) {
            EXPECT_EQ(State.MicroStates[1].Observations[0].Type, PITCH);
        }
    }
    EXPECT_EQ(States[2].MicroStates[0].Observations[1].Label, "key_click");
    EXPECT_EQ(States[3].MicroStates[0].Observations[1].Label, "aeolian");
}

TEST(MultiGlissParsing, ExpandsIntervalsOnlyForMultiAndPreservesTuning) {
    Score Parser;
    auto [Config, States] = Parser.Parse(Asset.parent_path() / "multi-gliss.scofo");
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence || State.Type == FIRSTEVENT; });
    ASSERT_EQ(States.size(), 8U);
    const std::vector<std::vector<double>> Expected = {
        {60, 60.5, 61, 61.5, 62, 62.5, 63, 63.5, 64, 64.5, 65, 65.5, 66, 66.5, 67},
        {67, 66.5, 66, 65.5, 65, 64.5, 64, 63.5, 63, 62.5, 62, 61.5, 61, 60.5, 60},
        {60, 60.5, 61, 61.5, 62, 61.5, 61, 60.5, 60},
        {60}, {60}};
    for (size_t I = 0; I < Expected.size(); ++I) {
        const auto &State = States[I];
        EXPECT_EQ(State.Type, MULTI);
        EXPECT_EQ(State.MicroTopologyType, LEFT_RIGHT);
        EXPECT_TRUE(State.Observations.empty());
        EXPECT_DOUBLE_EQ(State.Duration, 2.0);
        EXPECT_DOUBLE_EQ(State.OnsetExpected, 2.0 * I);
        ASSERT_EQ(State.MicroStates.size(), Expected[I].size());
        for (size_t K = 0; K < Expected[I].size(); ++K) {
            ASSERT_EQ(State.MicroStates[K].Observations.size(), 1U);
            const auto &Pitch = State.MicroStates[K].Observations[0];
            EXPECT_EQ(Pitch.Type, PITCH);
            EXPECT_DOUBLE_EQ(Pitch.Midi, Expected[I][K]);
            EXPECT_NEAR(Pitch.Freq, 440.0 * std::pow(2.0, (Expected[I][K] - 69.0) / 12.0), 1e-10);
        }
    }
    const auto &Chord = States[5];
    EXPECT_EQ(Chord.Type, CHORD);
    EXPECT_EQ(Chord.MicroTopologyType, NO_MICROSTATES);
    EXPECT_TRUE(Chord.MicroStates.empty());
    ASSERT_EQ(Chord.Observations.size(), 2U);
    EXPECT_DOUBLE_EQ(Chord.Observations[0].Midi, 60);
    EXPECT_DOUBLE_EQ(Chord.Observations[1].Midi, 67);
    const auto &Trill = States[6];
    EXPECT_EQ(Trill.Type, TRILL);
    EXPECT_EQ(Trill.MicroTopologyType, UNORDERED);
    ASSERT_EQ(Trill.MicroStates.size(), 2U);
    EXPECT_DOUBLE_EQ(Trill.MicroStates[0].Observations[0].Midi, 60);
    EXPECT_DOUBLE_EQ(Trill.MicroStates[1].Observations[0].Midi, 67);
    ASSERT_EQ(States[7].MicroStates.size(), 5U);
    for (size_t K = 0; K < 5; ++K) {
        const auto &Pitch = States[7].MicroStates[K].Observations[0];
        EXPECT_DOUBLE_EQ(Pitch.Midi, 60.25 + K * 0.5);
        EXPECT_NEAR(Pitch.Freq, 442.0 * std::pow(2.0, (Pitch.Midi - 69.0) / 12.0), 1e-10);
    }
}

TEST(MicroStateTransitions, TopologyDurationWeightsAndActiveParent) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(MULTI, 4, 8), Chain(PTECH, 4, 10)});
    for (int I = 0; I < 3; ++I) {
        EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, I, I), 0.5);
        EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, I, I) + Access::Transition(Forward, 0, I, I + 1), 1.0);
    }
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, 3, 3), 1.0);
    for (const auto &Pair : {std::pair{1, 0}, {0, 2}, {2, 0}, {-1, 0}, {0, 4}, {4, 4}}) {
        EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, Pair.first, Pair.second), 0.0);
    }
    // Parent 1 must be used even though the decoded score state is still 0.
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 1, 0, 1), 0.4);
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 1, 1, 2), 0.4);
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 1, 2, 3), 0.4);
    Forward.GetStates()[1].MicroStates[1].DurationWeight = 3;
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 1, 1, 2), 0.2);
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 1, 2, 3), 0.6);
    Access::Tempo(Forward, 2.0);
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, 0, 1), 0.25);
    Access::Tempo(Forward, 0.01);
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, 0, 1), 1.0);
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 1, 1, 2), 1.0);
    Forward.GetStates()[0].MicroTopologyType = UNORDERED;
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, 0, 1), 0.0);
}

TEST(MicroStateParsing, RejectsEmptyPitchAndLabelGroups) {
    Score Parser;
    const auto [Config, States] = Parser.Parse(Asset.parent_path() / "microstates-invalid.scofo");
    for (const auto &State : States) {
        EXPECT_EQ(State.Type, FIRSTEVENT);
        EXPECT_TRUE(State.MicroStates.empty());
    }
}

TEST(MicroStateForward, MatchesExhaustivePathsForEveryEntryTime) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(MULTI, 3, 6)});
    const std::vector<std::vector<double>> Emissions = {
        {0.8, 0.2, 0.1}, {0.3, 0.9, 0.2}, {0.4, 0.7, 0.8}, {0.1, 0.2, 0.9}, {0.7, 0.4, 0.3}};
    for (int T = 0; T < static_cast<int>(Emissions.size()); ++T) {
        Access::Advance(Forward, T, 5, Emissions[T]);
        for (int U = 1; U <= T + 1; ++U) {
            const auto Expected = Segment(Emissions, {0.5, 0.5, 0.0}, T - U + 1, T);
            for (size_t K = 0; K < 3; ++K) {
                EXPECT_NEAR(Alpha(Forward.GetStates()[0].MicroStates[K], U), Expected[K], 1e-14);
            }
        }
    }
}

TEST(MicroStateForward, OrderedEvidenceWinsForMultiAndPtech) {
    for (EventType Type : {MULTI, PTECH}) {
        double Likelihood[2] = {};
        for (int Reverse = 0; Reverse < 2; ++Reverse) {
            OnlineForward Forward;
            Access::Setup(Forward, {Chain(Type, 4, 4)});
            for (int T = 0; T < 4; ++T) {
                std::vector<double> Emissions(4, 0.01);
                Emissions[Reverse ? 3 - T : T] = 0.99;
                Access::Advance(Forward, T, 4, Emissions);
            }
            for (const auto &Micro : Forward.GetStates()[0].MicroStates) {
                Likelihood[Reverse] += Alpha(Micro, 4);
            }
        }
        EXPECT_GT(Likelihood[0], Likelihood[1] * 1000.0);
    }
}

TEST(MicroStateForward, SemiMarkovUsesSegmentsAndPathPosterior) {
    for (EventType Type : {MULTI, PTECH}) {
        OnlineForward Forward;
        Access::Setup(Forward, {Chain(MULTI, 1, 8), Chain(Type, 4, 8)});
        Forward.GetStates()[1].InitProb = 0.3;
        const std::vector<double> Incoming = {0.2, 0.4, 0.1, 0.5};
        const std::vector<std::vector<double>> Emissions = {
            {0.6, 0.7, 0.8, 0.99}, {0.3, 0.8, 0.2, 0.9}, {0.2, 0.7, 0.9, 0.1}, {0.1, 0.3, 0.4, 0.95}};
        const std::vector<double> Q = {0.5, 0.5, 0.5, 0.0};
        for (int T = 0; T < 4; ++T) {
            Access::Frame(Forward, T);
            Forward.GetStates()[0].ExitProb[T] = Incoming[T];
            Access::Emissions(Forward, 1, Emissions[T]);
            Forward.GetStates()[1].BestObs[T] = 1e-50; // Must not enter the recursion.
            Access::SemiMarkov(Forward, 1);
            double ExpectedForward = 0.0, ExpectedExit = 0.0;
            std::vector<double> Posterior(4, 0.0);
            for (int U = 1; U <= T + 1; ++U) {
                const auto Paths = Segment(Emissions, Q, T - U + 1, T);
                const double Entry = U == T + 1 ? 0.3 : Incoming[T - U];
                for (size_t K = 0; K < 4; ++K) {
                    const double Contribution = Access::Occupancy(Forward, 1, U, true) * Paths[K] * Entry;
                    ExpectedForward += Contribution;
                    ExpectedExit += Access::Occupancy(Forward, 1, U, false) * Paths[K] * Entry;
                    Posterior[K] += Contribution;
                }
            }
            EXPECT_NEAR(Forward.GetStates()[1].Forward[T], ExpectedForward, 1e-13);
            EXPECT_NEAR(Forward.GetStates()[1].ExitProb[T], ExpectedExit, 1e-13);
            EXPECT_EQ(Forward.GetStates()[1].BestMicroStateIndex,
                      std::distance(Posterior.begin(), std::max_element(Posterior.begin(), Posterior.end())));
            if (T == 0) {
                EXPECT_EQ(Forward.GetStates()[1].BestMicroStateIndex, 0); // Instantaneous winner is 3.
            }
        }
    }
}

TEST(MicroStateForward, GlobalScalingMatchesUnnormalizedReference) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(MULTI, 2, 4), Chain(PTECH, 4, 6)});
    Forward.GetStates()[0].InitProb = 0.7;
    Forward.GetStates()[1].InitProb = 0.3;
    const std::vector<std::vector<std::vector<double>>> Emissions = {
        {{0.8, 0.2}, {0.3, 0.9}, {0.4, 0.8}, {0.1, 0.9}, {0.7, 0.3}},
        {{0.7, 0.2, 0.1, 0.05},
         {0.1, 0.8, 0.3, 0.1},
         {0.2, 0.3, 0.8, 0.1},
         {0.1, 0.2, 0.3, 0.9},
         {0.7, 0.3, 0.2, 0.1}}};
    const std::vector<std::vector<double>> Q = {{0.5, 0.0}, {2.0 / 3, 2.0 / 3, 2.0 / 3, 0.0}};
    double Exit[2][5] = {};
    double Evidence[5] = {};
    for (int T = 0; T < 5; ++T) {
        Access::Frame(Forward, T);
        double Expected[2] = {};
        for (int J = 0; J < 2; ++J) {
            Access::Emissions(Forward, J, Emissions[J][T]);
            for (int U = 1; U <= T + 1; ++U) {
                const auto Paths = Segment(Emissions[J], Q[J], T - U + 1, T);
                const double B = std::accumulate(Paths.begin(), Paths.end(), 0.0);
                const double Entry = U == T + 1 ? (J == 0 ? 0.7 : 0.3) : (J == 0 ? 0.0 : Exit[0][T - U]);
                Expected[J] += Access::Occupancy(Forward, J, U, true) * B * Entry;
                Exit[J][T] += Access::Occupancy(Forward, J, U, false) * B * Entry;
            }
        }
        Evidence[T] = Expected[0] + Expected[1];
        const double N = Access::Normalize(Forward);
        EXPECT_NEAR(N, Evidence[T] / (T == 0 ? 1.0 : Evidence[T - 1]), 1e-13);
        for (int J = 0; J < 2; ++J) {
            EXPECT_NEAR(Forward.GetStates()[J].Forward[T], Expected[J] / Evidence[T], 1e-13);
            EXPECT_NEAR(Forward.GetStates()[J].ExitProb[T], Exit[J][T] / Evidence[T], 1e-13);
            for (int U = 1; U <= T + 1; ++U) {
                const auto Paths = Segment(Emissions[J], Q[J], T - U + 1, T);
                const double Scale = (T - U < 0 ? 1.0 : Evidence[T - U]) / Evidence[T];
                for (size_t K = 0; K < Paths.size(); ++K) {
                    EXPECT_NEAR(Alpha(Forward.GetStates()[J].MicroStates[K], U), Paths[K] * Scale, 1e-12);
                }
            }
        }
    }
}

TEST(MicroStateForward, ResizesResetsAndDoesNotReuseMissingFrames) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(MULTI, 2, 4)});
    Access::Advance(Forward, 0, 2, {0.8, 0.2});
    Access::Advance(Forward, 1, 4, {0.3, 0.9});
    EXPECT_NEAR(Alpha(Forward.GetStates()[0].MicroStates[1], 2), 0.8 * 0.5 * 0.9, 1e-14);
    EXPECT_EQ(Forward.GetStates()[0].MicroStates[0].LogForwardByAge.size(), 5U);
    Access::Advance(Forward, 2, 1, {0.4, 0.7});
    EXPECT_EQ(Forward.GetStates()[0].MicroStates[0].LogForwardByAge.size(), 2U);
    Access::Advance(Forward, 4, 4, {0.5, 0.8});
    EXPECT_DOUBLE_EQ(Alpha(Forward.GetStates()[0].MicroStates[0], 2), 0.0);
    EXPECT_DOUBLE_EQ(Alpha(Forward.GetStates()[0].MicroStates[1], 2), 0.0);
    Forward.ResetDecoding();
    for (const auto &Micro : Forward.GetStates()[0].MicroStates) {
        EXPECT_TRUE(Micro.LogForwardByAge.empty());
        EXPECT_DOUBLE_EQ(Micro.CurrentEmission, 0.0);
        EXPECT_EQ(Micro.BestObservationIndex, -1);
    }
    EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, -1);
    Forward.ClearStates();
    EXPECT_TRUE(Forward.GetStates().empty());
}

TEST(MicroStateForward, HandlesEmptySingleAndLongImpossibleHypotheses) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(MULTI, 0, 0)});
    Access::Advance(Forward, 0, 0, {});
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, 0, 0), 0.0);
    Access::Setup(Forward, {Chain(MULTI, 1, 0)});
    EXPECT_DOUBLE_EQ(Access::Transition(Forward, 0, 0, 0), 1.0);
    Access::Advance(Forward, 0, 0, {0.5});
    Access::Advance(Forward, 1, 2, {0.8});
    EXPECT_DOUBLE_EQ(Alpha(Forward.GetStates()[0].MicroStates[0], 2), 0.0);
    Access::Setup(Forward, {Chain(MULTI, 1, 1)});
    // No possible entry, but diagnostic age hypotheses continue to be evaluated.
    for (int T = 0; T < 80; ++T) {
        Access::Frame(Forward, T);
        Access::Emissions(Forward, 0, {0.5});
        Access::Normalize(Forward);
        EXPECT_TRUE(std::isfinite(Forward.GetStates()[0].Forward[T]));
        for (double Value : Forward.GetStates()[0].MicroStates[0].LogForwardByAge) {
            EXPECT_TRUE(std::isfinite(Value));
        }
    }
}

TEST(MicroStateEmissions, UtechMaxGatesLabelsAndNotifiesAlternatives) {
    Score Parser;
    auto [Config, States] = Parser.Parse(Asset);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ASSERT_EQ(States.size(), 6U);
    OnlineForward Forward;
    Access::Setup(Forward, {States[3]});
    Access::Frame(Forward, 0);
    Description Desc{};
    Desc.Onset = 0.1;
    Desc.SilenceProb = 0.2;
    Desc.ExtendedTechProb = 0.5;
    Desc.ONNX = {{"jet_whistle", 0.6f}, {"aeolian", 0.9f}};
    Access::Observe(Forward, Desc);
    auto &State = Forward.GetStates()[0];
    EXPECT_DOUBLE_EQ(State.BestObs[0], static_cast<double>(0.9f) * (0.5 * 0.8));
    EXPECT_EQ(State.BestMicroStateIndex, 0);
    EXPECT_EQ(State.BestMicroObservationIndex, 1);
    Access::Notify(Forward, 0);
    auto Actions = Forward.GetAudioStateChangeActions();
    ASSERT_EQ(Actions.size(), 1U);
    EXPECT_EQ(std::get<std::string>(Actions[0].Args[1]), "aeolian");
    Desc.Onset = 0.8;
    Access::Observe(Forward, Desc);
    EXPECT_DOUBLE_EQ(State.BestObs[0], static_cast<double>(0.9f) * (0.5 * 0.8));
    EXPECT_EQ(State.BestMicroStateIndex, 0);
    Desc.Onset = 0;
    Desc.ONNX.clear();
    Access::Observe(Forward, Desc);
    EXPECT_DOUBLE_EQ(State.BestObs[0], std::numeric_limits<double>::min());
    EXPECT_EQ(State.BestMicroStateIndex, -1);
    for (const auto &Micro : State.MicroStates) {
        EXPECT_TRUE(Micro.LogForwardByAge.empty());
    }
    MarkovMicroState Empty;
    EXPECT_DOUBLE_EQ(Access::Emission(Forward, Empty, Desc), 0.0);
}

TEST(MicroStateEmissions, PtechReportsWinningLabelAmongSoundedAlternatives) {
    Score Parser;
    auto [Config, States] = Parser.Parse(Asset);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ASSERT_EQ(States.size(), 6U);
    OnlineForward Forward;
    Access::Setup(Forward, {States[2]});
    Forward.GetStates()[0].InitProb = 1.0;
    Description Desc{};
    Desc.Onset = 1.0;
    Desc.SilenceProb = 1.0;
    Desc.ExtendedTechProb = 1.0;
    Access::Frame(Forward, 0);
    Access::Observe(Forward, Desc);
    Access::Normalize(Forward);
    EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, -1);
    Desc.Onset = 0;
    Desc.SilenceProb = 0.0;
    Desc.ONNX = {{"slap", 0.3f}, {"key_click", 0.8f}};
    Access::Frame(Forward, 1);
    Access::Observe(Forward, Desc);
    Access::Normalize(Forward);
    EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, 0);
    EXPECT_EQ(Forward.GetStates()[0].BestMicroObservationIndex, 1);
    Access::Notify(Forward, 0);
    const auto Actions = Forward.GetAudioStateChangeActions();
    ASSERT_EQ(Actions.size(), 1U);
    EXPECT_EQ(std::get<std::string>(Actions[0].Args[1]), "key_click");
}

TEST(MicroStateForward, SingleMicroStateMatchesOrdinarySemiMarkov) {
    OnlineForward Nested, Direct;
    auto Note = Chain(NOTE, 0, 8);
    Note.MicroTopologyType = NO_MICROSTATES;
    Note.Observations.push_back({PITCH, 440.0});
    Access::Setup(Nested, {Chain(MULTI, 1, 8)});
    Access::Setup(Direct, {Note});
    Nested.GetStates()[0].InitProb = 1.0;
    Direct.GetStates()[0].InitProb = 1.0;
    for (int T = 0; T < 6; ++T) {
        Access::Frame(Nested, T);
        Access::Frame(Direct, T);
        const double Emission = 0.2 + 0.1 * T;
        Access::Emissions(Nested, 0, {Emission});
        Direct.GetStates()[0].BestObs[T] = Emission;
        EXPECT_NEAR(Access::Normalize(Nested), Access::Normalize(Direct), 1e-13);
        EXPECT_NEAR(Nested.GetStates()[0].ExitProb[T], Direct.GetStates()[0].ExitProb[T], 1e-13);
    }
}

TEST(MicroStateIntegration, DiscoversNestedDescriptorsAndValidatesLabels) {
    ::OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    std::string Errors;
    Scofo.SetErrorCallback([&](const spdlog::details::log_msg &Log, void *) {
        if (Log.level == spdlog::level::err) {
            Errors.append(Log.payload.data(), Log.payload.size());
        }
    });
    // The fixture has labels but intentionally no ONNX model.
    EXPECT_FALSE(Scofo.LoadScore(Asset));
    EXPECT_NE(Errors.find("slap"), std::string::npos);
    const auto Config = Scofo.GetConfiguration();
    for (Descriptors Descriptor : {ONNX, EXTENDEDTECHNIQUE}) {
        EXPECT_NE(std::find(Config.RequestedDescriptors.begin(), Config.RequestedDescriptors.end(), Descriptor),
                  Config.RequestedDescriptors.end());
    }
}

#if defined(OPENSCOFO_LUA)
TEST(MicroStateIntegration, LuaExportsPhasesLabelsAndDurationWeights) {
    Score Parser;
    auto [Config, States] = Parser.Parse(Asset);
    std::erase_if(States, [](const auto &State) { return State.IsInterEventSilence; });
    ::OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    Scofo.GetStates() = States;
    EXPECT_TRUE(Scofo.LuaExecute(R"(
        local s = require('OpenScofo').get_states()
        assert(#s[2].microstates == 9 and #s[2].audiostates == 0)
        assert(s[2].micro_topology == 2)
        assert(#s[3].microstates == 2 and s[3].micro_topology == 1)
        assert(s[3].microstates[1].observations[2].label == 'key_click')
        assert(s[3].microstates[1].duration_weight == 1)
        assert(#s[4].microstates == 1 and s[4].micro_topology == 1)
        assert(s[4].microstates[1].observations[2].label == 'aeolian')
    )")) << Scofo.LuaGetError();
}
#endif

} // namespace
