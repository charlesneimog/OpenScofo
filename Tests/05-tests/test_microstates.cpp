#include <gtest/gtest.h>
#include <OpenScofo.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <limits>
#include <numeric>

namespace OpenScofo {

// Keep inference tests independent of FFT/ONNX and public APIs.
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

    static void Emissions(OnlineForward &Forward, int StateIndex, const std::vector<double> &Values) {
        auto &State = Forward.m_States[StateIndex];
        for (size_t K = 0; K < Values.size(); ++K) {
            State.MicroStates[K].CurrentEmission = Values[K];
            State.MicroStates[K].BestObservationIndex = 0;
        }
    }

    static void Advance(OnlineForward &Forward, int Time, const std::vector<double> &Values) {
        Frame(Forward, Time);
        Emissions(Forward, 0, Values);
        Forward.UpdateMicroStateForward(Forward.m_States[0]);
    }

    static void SemiMarkov(OnlineForward &Forward, int Index) {
        Forward.SemiMarkov(Forward.m_States[Index], Index);
        Forward.UpdateMicroStatePosterior(Forward.m_States[Index]);
    }

    static double SegmentLikelihood(OnlineForward &Forward, int Age, double Incoming = 1.0) {
        return Forward.GetSegmentLikelihood(Age, Incoming);
    }

    static double LogPathCount(OnlineForward &Forward, size_t Count, int Age) {
        return Forward.GetInternalLogPathCounts(Count, Age)[Age];
    }

    static void Normalization(OnlineForward &Forward, int Time, double Value) {
        Forward.m_Normalization[Time % Forward.m_BufferSize] = Value;
    }

    static int BufferSize(const OnlineForward &Forward) {
        return Forward.m_BufferSize;
    }

    static void UseBufferSize(OnlineForward &Forward, int Size) {
        Forward.m_BufferSize = Size;
        // Poison unused history so reading a nonexistent frame changes the result.
        Forward.m_Normalization.assign(Size, 0.25);
        for (auto &State : Forward.m_States) {
            State.BestObs.assign(Size, 0.9);
            State.ExitProb.assign(Size, 1000.0);
            State.Forward.assign(Size, 0.0);
        }
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

    static void ResetCaches(OnlineForward &Forward) {
        Forward.ResetCaches();
    }

    static void Notify(OnlineForward &Forward, int Index) {
        Forward.NotifyAudioStateChange(Index);
    }

    static double ScoreTransition(OnlineForward &Forward, int I, int J) {
        return Forward.GetSemiMarkovTransitionProbability(I, J);
    }

    static void Timing(OnlineForward &Forward, double Period, double BlockDuration) {
        Forward.m_PsiN1 = Period;
        Forward.m_BlockDur = BlockDuration;
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

double PathCount(int Age, size_t Count) {
    double Sum = 0.0, Binomial = 1.0;
    for (size_t R = 0; R < Count && R < static_cast<size_t>(Age); ++R) {
        if (R > 0) {
            Binomial *= static_cast<double>(Age - R) / static_cast<double>(R);
        }
        Sum += Binomial;
    }
    return Sum;
}

// Endpoint contribution to the path-averaged likelihood, including outer scaling.
double Alpha(const ScoreState &State, size_t K, int Age) {
    const auto &History = State.MicroStates[K].LogForwardByAge;
    if (Age >= static_cast<int>(History.size()) || History[Age] == LogZero) {
        return 0.0;
    }
    return std::exp(History[Age] - std::log(PathCount(Age, State.MicroStates.size())));
}

// Independent reference: enumerate complete paths instead of using the decoder's
// dynamic programming or binomial-count cache. Count even zero-emission paths.
std::vector<double> Segment(const std::vector<std::vector<double>> &Emissions, int Start, int End,
                            bool Ordered = true) {
    const size_t Count = Emissions[0].size();
    std::vector<double> Result(Count, 0.0);
    if (!Ordered) {
        Result[0] = 1.0;
        for (int T = Start; T <= End; ++T) {
            Result[0] *= Emissions[T][0];
        }
        return Result;
    }
    size_t Paths = 0;
    std::function<void(int, size_t, double)> Visit = [&](int T, size_t K, double Product) {
        Product *= Emissions[T][K];
        if (T == End) {
            Result[K] += Product;
            ++Paths;
        } else {
            Visit(T + 1, K, Product);
            if (K + 1 < Count) {
                Visit(T + 1, K + 1, Product);
            }
        }
    };
    if (Count > 0) {
        Visit(Start, 0, 1.0);
        for (double &Weight : Result) {
            Weight /= static_cast<double>(Paths);
        }
    }
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

TEST(InterEventSilence, BothSegmentModelsAcceptTheDirectBypass) {
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
    EXPECT_EQ(Multi.Type, GLISS);
    EXPECT_EQ(Multi.HSMMType, SEMIMARKOV);
    EXPECT_EQ(Multi.MicroTopologyType, LEFT_RIGHT);
    EXPECT_TRUE(Multi.Observations.empty());
    EXPECT_DOUBLE_EQ(Multi.Duration, 2.0);
    ASSERT_EQ(Multi.MicroStates.size(), 9U);
    for (size_t K = 0; K < 9; ++K) {
        ASSERT_EQ(Multi.MicroStates[K].Observations.size(), 1U);
        EXPECT_EQ(Multi.MicroStates[K].Observations[0].Type, PITCH);
        EXPECT_DOUBLE_EQ(Multi.MicroStates[K].Observations[0].Midi, 60.0 + 0.5 * K);
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
        EXPECT_EQ(State.Type, GLISS);
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

TEST(MicroStatePaths, CountsAllAdmissibleMonotonicPathsIncludingEarlyEndpoints) {
    OnlineForward Forward;
    for (size_t Count : {1U, 2U, 3U, 6U}) {
        const std::vector<std::vector<double>> Emissions(9, std::vector<double>(Count, 1.0));
        for (int Age = 1; Age <= 9; ++Age) {
            const auto Endpoints = Segment(Emissions, 0, Age - 1);
            EXPECT_NEAR(std::accumulate(Endpoints.begin(), Endpoints.end(), 0.0), 1.0, 1e-14);
            EXPECT_NEAR(std::exp(Access::LogPathCount(Forward, Count, Age)), PathCount(Age, Count), 1e-10);
        }
    }
    EXPECT_NEAR(std::exp(Access::LogPathCount(Forward, 3, 12)), 67.0, 1e-12);
    EXPECT_NEAR(std::exp(Access::LogPathCount(Forward, 2, 12)), 12.0, 1e-12);
    EXPECT_DOUBLE_EQ(Access::LogPathCount(Forward, 1, 100), 0.0);
    // Counts larger than double's range remain usable without exponentiation.
    const double LogCount = Access::LogPathCount(Forward, 1500, 1500);
    EXPECT_GT(LogCount, std::log(std::numeric_limits<double>::max()));
    EXPECT_NEAR(LogCount, 1499.0 * std::log(2.0), 1e-10);
}

TEST(MicroStateParsing, RejectsEmptyPitchAndLabelGroups) {
    Score Parser;
    const auto [Config, States] = Parser.Parse(Asset.parent_path() / "microstates-invalid.scofo");
    for (const auto &State : States) {
        EXPECT_EQ(State.Type, FIRSTEVENT);
        EXPECT_TRUE(State.MicroStates.empty());
    }
}

TEST(MicroStateForward, MatchesEnumeratedPathAveragesForEveryEntryTime) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(GLISS, 3, 6)});
    const std::vector<std::vector<double>> Emissions = {
        {0.8, 0.2, 0.1}, {0.3, 0.9, 0.2}, {0.4, 0.7, 0.8}, {0.1, 0.2, 0.9}, {0.7, 0.4, 0.3}};
    for (int T = 0; T < static_cast<int>(Emissions.size()); ++T) {
        Access::Advance(Forward, T, Emissions[T]);
        for (int U = 1; U <= T + 1; ++U) {
            const auto Expected = Segment(Emissions, T - U + 1, T);
            for (size_t K = 0; K < 3; ++K) {
                EXPECT_NEAR(Alpha(Forward.GetStates()[0], K, U), Expected[K], 1e-14);
            }
        }
    }
}

TEST(MicroStatePaths, UnitEmissionsHaveUnitLikelihoodAtEverySegmentAge) {
    for (size_t Count : {1U, 2U, 3U, 6U}) {
        OnlineForward Forward;
        Access::Setup(Forward, {Chain(GLISS, Count, 100.0)});
        Forward.GetStates()[0].InitProb = 1.0;
        for (int T = 0; T < 160; ++T) {
            Access::Frame(Forward, T);
            Access::Emissions(Forward, 0, std::vector<double>(Count, 1.0));
            Access::SemiMarkov(Forward, 0);
            for (int Age = 1; Age <= T + 1; ++Age) {
                EXPECT_NEAR(Access::SegmentLikelihood(Forward, Age), 1.0, 1e-12);
            }
            EXPECT_NEAR(Forward.GetStates()[0].Forward[T], Access::Occupancy(Forward, 0, T + 1, true), 1e-12);
            EXPECT_NEAR(Forward.GetStates()[0].ExitProb[T], Access::Occupancy(Forward, 0, T + 1, false), 1e-12);
        }
    }
}

TEST(MicroStatePaths, AllowsStayingAndAdvancingButForbidsSkippingAndGoingBackward) {
    const std::vector<std::vector<std::vector<double>>> Examples = {
        {{1, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 1}},
        {{1, 0, 0}, {0, 0, 1}},
        {{1, 0, 0}, {0, 1, 0}, {1, 0, 0}},
        {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}}};
    for (size_t Example = 0; Example < Examples.size(); ++Example) {
        OnlineForward Forward;
        Access::Setup(Forward, {Chain(GLISS, 3, 8.0)});
        Forward.GetStates()[0].InitProb = 1.0;
        const auto &Emissions = Examples[Example];
        for (int T = 0; T < static_cast<int>(Emissions.size()); ++T) {
            Access::Frame(Forward, T);
            Access::Emissions(Forward, 0, Emissions[T]);
            Access::SemiMarkov(Forward, 0);
        }
        const double Likelihood = Access::SegmentLikelihood(Forward, Emissions.size());
        EXPECT_NEAR(Likelihood, Example == 0 ? 1.0 / 16.0 : 0.0, 1e-14);
    }
}

TEST(MicroStatePaths, ThreeMicrostatesExplainAnEighteenFrameTrajectory) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(GLISS, 3, 18.0)});
    Forward.GetStates()[0].InitProb = 1.0;
    for (int T = 0; T < 18; ++T) {
        Access::Frame(Forward, T);
        std::vector<double> Emissions(3, 0.0);
        Emissions[T / 6] = 1.0; // C C C C C C -> D D D D D D -> E E E E E E.
        Access::Emissions(Forward, 0, Emissions);
        Access::SemiMarkov(Forward, 0);
    }
    EXPECT_NEAR(Access::SegmentLikelihood(Forward, 18), 1.0 / 154.0, 1e-14);
    EXPECT_NEAR(Forward.GetStates()[0].Forward[17], Access::Occupancy(Forward, 0, 18, true) / 154.0, 1e-14);
    EXPECT_NEAR(Forward.GetStates()[0].ExitProb[17], Access::Occupancy(Forward, 0, 18, false) / 154.0, 1e-14);
    EXPECT_GT(Forward.GetStates()[0].Forward[17], 0.0);
    EXPECT_GT(Forward.GetStates()[0].ExitProb[17], 0.0);
    EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, 2);
}

TEST(MicroStatePaths, ParentCanExitWithoutReachingTheFinalMicrostate) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(GLISS, 3, 8.0)});
    Forward.GetStates()[0].InitProb = 0.4;
    for (int T = 0; T < 10; ++T) {
        Access::Frame(Forward, T);
        Access::Emissions(Forward, 0, {0.8, 0.0, 0.0});
        Access::SemiMarkov(Forward, 0);
        const int Age = T + 1;
        const double InitialSegment = 0.4 * std::pow(0.8, Age) / PathCount(Age, 3);
        // Both initial terms use every current endpoint, including microstate zero.
        EXPECT_NEAR(Forward.GetStates()[0].Forward[T], Access::Occupancy(Forward, 0, Age, true) * InitialSegment,
                    1e-14);
        EXPECT_NEAR(Forward.GetStates()[0].ExitProb[T], Access::Occupancy(Forward, 0, Age, false) * InitialSegment,
                    1e-14);
        EXPECT_GT(Forward.GetStates()[0].ExitProb[T], 0.0);
        EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, 0);
    }
}

TEST(MicroStatePaths, ParentDurationChangesOuterOccupancyButNotConditionalLikelihood) {
    OnlineForward Short, Long;
    Access::Setup(Short, {Chain(GLISS, 3, 4.0)});
    Access::Setup(Long, {Chain(GLISS, 3, 12.0)});
    Short.GetStates()[0].InitProb = Long.GetStates()[0].InitProb = 1.0;
    const std::vector<std::vector<double>> Emissions = {{0.9, 0.1, 0.2}, {0.8, 0.2, 0.1}, {0.7, 0.5, 0.1},
                                                        {0.3, 0.9, 0.2}, {0.2, 0.9, 0.4}, {0.1, 0.7, 0.8},
                                                        {0.1, 0.2, 0.9}, {0.1, 0.2, 0.8}};
    for (int T = 0; T < static_cast<int>(Emissions.size()); ++T) {
        for (OnlineForward *Forward : {&Short, &Long}) {
            Access::Frame(*Forward, T);
            Access::Emissions(*Forward, 0, Emissions[T]);
            Access::SemiMarkov(*Forward, 0);
        }
        for (int Age = 1; Age <= T + 1; ++Age) {
            EXPECT_NEAR(Access::SegmentLikelihood(Short, Age), Access::SegmentLikelihood(Long, Age), 1e-14);
        }
        for (size_t K = 0; K < 3; ++K) {
            EXPECT_EQ(Short.GetStates()[0].MicroStates[K].LogForwardByAge,
                      Long.GetStates()[0].MicroStates[K].LogForwardByAge);
        }
    }
    EXPECT_GT(std::abs(Access::Occupancy(Short, 0, 8, true) - Access::Occupancy(Long, 0, 8, true)), 1e-3);
    EXPECT_GT(std::abs(Short.GetStates()[0].Forward[7] - Long.GetStates()[0].Forward[7]), 1e-8);
    EXPECT_GT(std::abs(Short.GetStates()[0].ExitProb[7] - Long.GetStates()[0].ExitProb[7]), 1e-8);
}

TEST(MicroStatePaths, IncreasingParentDurationCanReuseOlderObservedSegments) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(GLISS, 3, 0.01)});
    Forward.GetStates()[0].InitProb = 1.0;
    for (int T = 0; T < 12; ++T) {
        Access::Frame(Forward, T);
        Access::Emissions(Forward, 0, {0.8, 0.8, 0.8});
        Access::SemiMarkov(Forward, 0);
    }
    Forward.GetStates()[0].Duration = 100.0;
    Access::Frame(Forward, 12);
    Access::Emissions(Forward, 0, {0.8, 0.8, 0.8});
    Access::SemiMarkov(Forward, 0);
    EXPECT_NEAR(Access::SegmentLikelihood(Forward, 13), std::pow(0.8, 13), 1e-14);
    EXPECT_NEAR(Forward.GetStates()[0].Forward[12], Access::Occupancy(Forward, 0, 13, true) * std::pow(0.8, 13), 1e-14);
}

TEST(MicroStateForward, ObservationPathsDoNotDependOnParentDurationTempoOrBlockSize) {
    for (double Duration : {0.01, 8.0, 1000.0}) {
        for (double Period : {0.01, 1.0, 2.0}) {
            for (double BlockDuration : {0.001, 1.0, 4.0}) {
                OnlineForward Forward;
                Access::Setup(Forward, {Chain(GLISS, 3, Duration)});
                Access::Timing(Forward, Period, BlockDuration);
                Access::Advance(Forward, 0, {0.8, 0.2, 0.1});
                Access::Advance(Forward, 1, {0.3, 0.9, 0.2});
                Access::Advance(Forward, 2, {0.4, 0.7, 0.8});
                EXPECT_NEAR(Alpha(Forward.GetStates()[0], 2, 3), 0.8 * 0.9 * 0.8 / 4.0, 1e-14);
            }
        }
    }
}

TEST(MicroStateForward, OrderedEvidenceWinsForMultiAndPtech) {
    for (EventType Type : {GLISS, PTECH}) {
        double Likelihood[2] = {};
        for (int Reverse = 0; Reverse < 2; ++Reverse) {
            OnlineForward Forward;
            Access::Setup(Forward, {Chain(Type, 4, 4)});
            for (int T = 0; T < 4; ++T) {
                std::vector<double> Emissions(4, 0.01);
                Emissions[Reverse ? 3 - T : T] = 0.99;
                Access::Advance(Forward, T, Emissions);
            }
            Likelihood[Reverse] = Alpha(Forward.GetStates()[0], 3, 4);
        }
        EXPECT_GT(Likelihood[0], Likelihood[1] * 1000.0);
    }
}

TEST(MicroStateForward, SemiMarkovUsesSegmentsAndPathPosterior) {
    for (EventType Type : {GLISS, PTECH}) {
        OnlineForward Forward;
        Access::Setup(Forward, {Chain(GLISS, 1, 8), Chain(Type, 4, 8)});
        Forward.GetStates()[1].InitProb = 0.3;
        const std::vector<double> Incoming = {0.2, 0.4, 0.1, 0.5};
        const std::vector<std::vector<double>> Emissions = {
            {0.6, 0.7, 0.8, 0.99}, {0.3, 0.8, 0.2, 0.9}, {0.2, 0.7, 0.9, 0.1}, {0.1, 0.3, 0.4, 0.95}};
        for (int T = 0; T < 4; ++T) {
            Access::Frame(Forward, T);
            Forward.GetStates()[0].ExitProb[T] = Incoming[T];
            Access::Emissions(Forward, 1, Emissions[T]);
            Forward.GetStates()[1].BestObs[T] = 1e-50; // Must not enter the recursion.
            Access::SemiMarkov(Forward, 1);
            double ExpectedForward = 0.0, ExpectedExit = 0.0;
            std::vector<double> Posterior(4, 0.0);
            for (int U = 1; U <= T + 1; ++U) {
                const auto Paths = Segment(Emissions, T - U + 1, T);
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
    Access::Setup(Forward, {Chain(GLISS, 2, 4), Chain(PTECH, 4, 6)});
    Forward.GetStates()[0].InitProb = 0.7;
    Forward.GetStates()[1].InitProb = 0.3;
    const std::vector<std::vector<std::vector<double>>> Emissions = {
        {{0.8, 0.2}, {0.3, 0.9}, {0.4, 0.8}, {0.1, 0.9}, {0.7, 0.3}},
        {{0.7, 0.2, 0.1, 0.05},
         {0.1, 0.8, 0.3, 0.1},
         {0.2, 0.3, 0.8, 0.1},
         {0.1, 0.2, 0.3, 0.9},
         {0.7, 0.3, 0.2, 0.1}}};
    double Exit[2][5] = {};
    double Evidence[5] = {};
    for (int T = 0; T < 5; ++T) {
        Access::Frame(Forward, T);
        double Expected[2] = {};
        for (int J = 0; J < 2; ++J) {
            Access::Emissions(Forward, J, Emissions[J][T]);
            for (int U = 1; U <= T + 1; ++U) {
                const auto Paths = Segment(Emissions[J], T - U + 1, T);
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
                const auto Paths = Segment(Emissions[J], T - U + 1, T);
                const double Scale = (T - U < 0 ? 1.0 : Evidence[T - U]) / Evidence[T];
                for (size_t K = 0; K < Paths.size(); ++K) {
                    EXPECT_NEAR(Alpha(Forward.GetStates()[J], K, U), Paths[K] * Scale, 1e-12);
                }
            }
        }
    }
}

TEST(MicroStateForward, ResetsAndDoesNotReuseMissingFrames) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(GLISS, 2, 4)});
    Access::Advance(Forward, 0, {0.8, 0.2});
    Access::Advance(Forward, 1, {0.3, 0.9});
    EXPECT_NEAR(Alpha(Forward.GetStates()[0], 1, 2), 0.8 * 0.9 / 2.0, 1e-14);
    Access::Advance(Forward, 2, {0.4, 0.7});
    EXPECT_NEAR(Alpha(Forward.GetStates()[0], 1, 2), 0.3 * 0.7 / 2.0, 1e-14);
    Access::Advance(Forward, 4, {0.5, 0.8});
    EXPECT_DOUBLE_EQ(Alpha(Forward.GetStates()[0], 0, 2), 0.0);
    EXPECT_DOUBLE_EQ(Alpha(Forward.GetStates()[0], 1, 2), 0.0);
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

TEST(MicroStateForward, ConfigurationCacheAndScoreResetsClearPaths) {
    for (int Reset = 0; Reset < 4; ++Reset) {
        OnlineForward Forward;
        Access::Setup(Forward, {Chain(GLISS, 3, 6)});
        Access::Advance(Forward, 0, {0.7, 0.1, 0.1});
        Access::Advance(Forward, 1, {0.1, 0.8, 0.1});
        EXPECT_GT(Alpha(Forward.GetStates()[0], 1, 2), 0.0);
        switch (Reset) {
        case 0:
            Forward.ResetDecoding();
            break;
        case 1:
            Access::ResetCaches(Forward);
            break;
        case 2: {
            Configuration Config;
            Forward.UpdateConfiguration(Config);
            break;
        }
        case 3:
            Forward.SetScoreStates(Forward.GetStates());
            break;
        }
        const auto &State = Forward.GetStates()[0];
        EXPECT_EQ(State.MicroForwardLastFrame, -1);
        for (const auto &Micro : State.MicroStates) {
            EXPECT_TRUE(Micro.LogForwardByAge.empty());
            EXPECT_DOUBLE_EQ(Micro.CurrentEmission, 0.0);
            EXPECT_EQ(Micro.BestObservationIndex, -1);
        }
        Access::Advance(Forward, 2, {0.5, 0.9, 0.9});
        EXPECT_DOUBLE_EQ(Alpha(State, 0, 1), 0.5);
        EXPECT_DOUBLE_EQ(Alpha(State, 1, 2), 0.0);
        EXPECT_DOUBLE_EQ(Alpha(State, 2, 3), 0.0);
    }
}

TEST(MicroStateForward, HandlesEmptySingleAndLongHypothesesWithoutIncomingMass) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(GLISS, 0, 0)});
    Access::Advance(Forward, 0, {});
    EXPECT_TRUE(Forward.GetStates()[0].MicroStates.empty());
    Access::Setup(Forward, {Chain(GLISS, 1, 0)});
    Access::Advance(Forward, 0, {0.5});
    EXPECT_DOUBLE_EQ(Alpha(Forward.GetStates()[0], 0, 1), 0.5);
    Access::Advance(Forward, 1, {0.8});
    EXPECT_DOUBLE_EQ(Alpha(Forward.GetStates()[0], 0, 2), 0.5 * 0.8);
    Access::Setup(Forward, {Chain(GLISS, 1, 1)});
    // No possible entry, but diagnostic age hypotheses continue to be evaluated.
    for (int T = 0; T < 80; ++T) {
        Access::Frame(Forward, T);
        Access::Emissions(Forward, 0, {0.5});
        Access::Normalize(Forward);
        EXPECT_TRUE(std::isfinite(Forward.GetStates()[0].Forward[T]));
        EXPECT_TRUE(std::isfinite(Forward.GetStates()[0].MicroStates[0].LogForwardByAge[1]));
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

TEST(MicroStateForward, SingleMicroStateMatchesOrdinaryLikelihoodForLongParentSegments) {
    OnlineForward Nested, Direct;
    auto Note = Chain(NOTE, 0, 8);
    Note.MicroTopologyType = NO_MICROSTATES;
    Note.Observations.push_back({PITCH, 440.0});
    Access::Setup(Nested, {Chain(GLISS, 1, 8)});
    Access::Setup(Direct, {Note});
    Nested.GetStates()[0].InitProb = 1.0;
    Direct.GetStates()[0].InitProb = 1.0;
    for (int T = 0; T < 6; ++T) {
        Access::Frame(Nested, T);
        Access::Frame(Direct, T);
        Access::Emissions(Nested, 0, {0.5});
        Direct.GetStates()[0].BestObs[T] = 0.5;
        Access::SemiMarkov(Nested, 0);
        Access::SemiMarkov(Direct, 0);
        EXPECT_NEAR(Nested.GetStates()[0].Forward[T], Direct.GetStates()[0].Forward[T], 1e-14);
        EXPECT_NEAR(Nested.GetStates()[0].ExitProb[T], Direct.GetStates()[0].ExitProb[T], 1e-14);
        EXPECT_NEAR(Access::SegmentLikelihood(Nested, T + 1), std::pow(0.5, T + 1), 1e-14);
        EXPECT_GT(Nested.GetStates()[0].Forward[T], 0.0);
        EXPECT_EQ(Nested.GetStates()[0].BestMicroStateIndex, 0);
    }
}

TEST(SemiMarkovSegments, OrdinaryProductsPreserveNormalizationInitializationAndCircularHistory) {
    for (EventType Type : {NOTE, REST, CHORD, PTECH, UTECH, TRILL}) {
        for (MicroTopology Topology : {NO_MICROSTATES, UNORDERED}) {
            OnlineForward Forward;
            auto State = Chain(Type, 0, 2.0);
            State.MicroTopologyType = Topology;
            Access::Setup(Forward, {State, State});
            Access::UseBufferSize(Forward, 5);
            Forward.GetStates()[1].InitProb = 0.3;
            std::vector<double> Observations, Normalizations, Incoming;
            for (int T = 0; T < 12; ++T) {
                Access::Frame(Forward, T);
                const int Buf = T % Access::BufferSize(Forward);
                Observations.push_back(0.2 + 0.03 * T);
                Normalizations.push_back(0.4 + 0.01 * T);
                Incoming.push_back(0.1 + 0.02 * T);
                Forward.GetStates()[1].BestObs[Buf] = Observations[T];
                Forward.GetStates()[0].ExitProb[Buf] = Incoming[T];
                Access::Normalization(Forward, T, Normalizations[T]);
                Access::SemiMarkov(Forward, 1);

                double ExpectedForward = 0.0, ExpectedExit = 0.0;
                for (int U = 1; U <= std::min(T + 1, 4); ++U) {
                    double Product = Observations[T];
                    for (int S = T - U + 1; S < T; ++S) {
                        Product *= Observations[S] / Normalizations[S];
                    }
                    const double Entry = U == T + 1 ? 0.3 : Incoming[T - U];
                    EXPECT_NEAR(Access::SegmentLikelihood(Forward, U), Product, 1e-14);
                    ExpectedForward += Access::Occupancy(Forward, 1, U, true) * Product * Entry;
                    ExpectedExit += Access::Occupancy(Forward, 1, U, false) * Product * Entry;
                }
                EXPECT_NEAR(Forward.GetStates()[1].Forward[Buf], ExpectedForward, 1e-14);
                EXPECT_NEAR(Forward.GetStates()[1].ExitProb[Buf], ExpectedExit, 1e-14);
            }
        }
    }
}

TEST(SemiMarkovSegments, OrdinaryZeroNormalizationAndObservationKeepExistingFloor) {
    OnlineForward Forward;
    auto State = Chain(NOTE, 0, 2.0);
    State.MicroTopologyType = NO_MICROSTATES;
    Access::Setup(Forward, {State});
    Forward.GetStates()[0].InitProb = 1.0;
    for (double Normalization : {0.0, std::numeric_limits<double>::min()}) {
        for (double PreviousObservation : {0.0, 0.5}) {
            Access::Frame(Forward, 1);
            Forward.GetStates()[0].BestObs[0] = PreviousObservation;
            Forward.GetStates()[0].BestObs[1] = 0.5;
            Access::Normalization(Forward, 0, Normalization);
            Access::SemiMarkov(Forward, 0);
            EXPECT_DOUBLE_EQ(Access::SegmentLikelihood(Forward, 2), 0.0);
            EXPECT_DOUBLE_EQ(Forward.GetStates()[0].Forward[1], 0.5 * std::numeric_limits<double>::min());
            EXPECT_DOUBLE_EQ(Forward.GetStates()[0].ExitProb[1], 0.5 * std::numeric_limits<double>::min());
        }
    }
}

TEST(SemiMarkovSegments, MixedModelsMatchUnnormalizedPathsAcrossCircularWraps) {
    OnlineForward Forward;
    auto Ordinary = Chain(NOTE, 0, 2.0);
    Ordinary.MicroTopologyType = NO_MICROSTATES;
    Access::Setup(Forward, {Ordinary, Chain(GLISS, 2, 4.0), Ordinary, Ordinary});
    Access::UseBufferSize(Forward, 5);
    Forward.GetStates()[0].InitProb = 0.4;
    Forward.GetStates()[1].InitProb = 0.3;
    Forward.GetStates()[2].InitProb = 0.2;
    Forward.GetStates()[3].InitProb = 0.1;
    const std::vector<double> Initial = {0.4, 0.3, 0.2, 0.1};
    std::vector<std::vector<std::vector<double>>> Emissions(4);
    std::vector<std::vector<double>> Exit(4);
    std::vector<double> Evidence;
    for (int T = 0; T < 12; ++T) {
        Access::Frame(Forward, T);
        const int Buf = T % Access::BufferSize(Forward);
        Emissions[0].push_back({0.3 + 0.03 * T});
        Emissions[1].push_back({0.7 - 0.04 * T, 0.2 + 0.05 * T});
        Emissions[2].push_back({0.6 - 0.02 * T});
        Emissions[3].push_back({0.4 + 0.02 * T});
        Forward.GetStates()[0].BestObs[Buf] = Emissions[0][T][0];
        Forward.GetStates()[1].BestObs[Buf] = 1e-50;
        Forward.GetStates()[2].BestObs[Buf] = Emissions[2][T][0];
        Forward.GetStates()[3].BestObs[Buf] = Emissions[3][T][0];
        Access::Emissions(Forward, 1, Emissions[1][T]);
        double Expected[4] = {}, ExpectedExit[4] = {};
        for (int J = 0; J < 4; ++J) {
            for (int U = 1; U <= std::min(T + 1, 4); ++U) {
                const auto Paths = Segment(Emissions[J], T - U + 1, T, J == 1);
                const double Product = std::accumulate(Paths.begin(), Paths.end(), 0.0);
                const double Entry = U == T + 1 ? Initial[J] : (J == 0 ? 0.0 : Exit[J - 1][T - U]);
                Expected[J] += Access::Occupancy(Forward, J, U, true) * Product * Entry;
                ExpectedExit[J] += Access::Occupancy(Forward, J, U, false) * Product * Entry;
            }
        }
        for (int J = 0; J < 4; ++J) {
            Exit[J].push_back(ExpectedExit[J]);
        }
        Evidence.push_back(std::accumulate(std::begin(Expected), std::end(Expected), 0.0));
        ASSERT_GT(Evidence[T], 0.0);
        const double N = Access::Normalize(Forward);
        EXPECT_NEAR(N, Evidence[T] / (T == 0 ? 1.0 : Evidence[T - 1]), 1e-13);
        for (int J = 0; J < 4; ++J) {
            EXPECT_NEAR(Forward.GetStates()[J].Forward[Buf], Expected[J] / Evidence[T], 1e-13);
            EXPECT_NEAR(Forward.GetStates()[J].ExitProb[Buf], ExpectedExit[J] / Evidence[T], 1e-13);
        }
    }
}

TEST(SemiMarkovSegments, C4D4E4HasDifferentLikelihoodFromE4D4C4) {
    const std::vector<double> Pitches = {60.0, 62.0, 64.0};
    double Likelihood[2] = {};
    for (int Reverse = 0; Reverse < 2; ++Reverse) {
        OnlineForward Forward;
        auto State = Chain(GLISS, 3, 6.0);
        for (size_t K = 0; K < Pitches.size(); ++K) {
            State.MicroStates[K].Observations[0].Midi = Pitches[K];
        }
        Access::Setup(Forward, {State});
        Forward.GetStates()[0].InitProb = 1.0;
        for (int T = 0; T < 3; ++T) {
            Access::Frame(Forward, T);
            std::vector<double> Emissions(3, 0.01);
            Emissions[Reverse ? 2 - T : T] = 0.99;
            Access::Emissions(Forward, 0, Emissions);
            Forward.GetStates()[0].BestObs[T] = 0.99; // Identical diagnostic maxima in either order.
            Access::SemiMarkov(Forward, 0);
        }
        Likelihood[Reverse] = Access::SegmentLikelihood(Forward, 3);
        EXPECT_NEAR(Forward.GetStates()[0].Forward[2], Access::Occupancy(Forward, 0, 3, true) * Likelihood[Reverse],
                    1e-14);
    }
    EXPECT_NEAR(Likelihood[0], 0.2450745, 1e-14);
    EXPECT_NEAR(Likelihood[1], 0.0000745, 1e-14);
    EXPECT_GT(Likelihood[0], 1000.0 * Likelihood[1]);
}

TEST(SemiMarkovSegments, PitchObservationsPreserveMultiOrderingAndUnorderedTrillLikelihood) {
    for (MicroTopology Topology : {LEFT_RIGHT, UNORDERED}) {
        double Likelihood[2] = {};
        for (int Reverse = 0; Reverse < 2; ++Reverse) {
            OnlineForward Forward;
            auto State = Chain(Topology == LEFT_RIGHT ? GLISS : TRILL, 3, 6.0);
            State.MicroTopologyType = Topology;
            for (size_t K = 0; K < 3; ++K) {
                auto &Pitch = State.MicroStates[K].Observations[0];
                Pitch.Midi = 60.0 + 2.0 * K; // C4, D4, E4.
                Pitch.Freq = 440.0 * std::pow(2.0, (Pitch.Midi - 69.0) / 12.0);
            }
            Access::Setup(Forward, {State});
            Forward.GetStates()[0].InitProb = 1.0;
            double Product = 1.0;
            std::vector<std::vector<double>> Emissions;
            for (int T = 0; T < 3; ++T) {
                Access::Frame(Forward, T);
                Description Desc{};
                Desc.SpectralMagnitudeFrameNorm = Forward.GetPitchTemplate(
                    State.MicroStates[Reverse ? 2 - T : T].Observations[0].Freq);
                for (double &Bin : Desc.SpectralMagnitudeFrameNorm) {
                    Bin += 1e-8;
                }
                const double Sum = std::accumulate(Desc.SpectralMagnitudeFrameNorm.begin(),
                                                    Desc.SpectralMagnitudeFrameNorm.end(), 0.0);
                for (double &Bin : Desc.SpectralMagnitudeFrameNorm) {
                    Bin /= Sum;
                }
                Access::Observe(Forward, Desc);
                const auto &Current = Forward.GetStates()[0];
                Product *= Current.BestObs[T];
                std::vector<double> FrameEmissions;
                for (const auto &Micro : Current.MicroStates) {
                    FrameEmissions.push_back(Micro.CurrentEmission);
                }
                Emissions.push_back(std::move(FrameEmissions));
                Access::SemiMarkov(Forward, 0);
            }
            Likelihood[Reverse] = Access::SegmentLikelihood(Forward, 3);
            if (Topology == LEFT_RIGHT) {
                const auto Endpoints = Segment(Emissions, 0, 2);
                EXPECT_NEAR(Likelihood[Reverse], std::accumulate(Endpoints.begin(), Endpoints.end(), 0.0), 1e-14);
                const int Best = std::distance(Endpoints.begin(), std::max_element(Endpoints.begin(), Endpoints.end()));
                EXPECT_EQ(Forward.GetStates()[0].BestMicroStateIndex, Best);
                EXPECT_EQ(Forward.GetStates()[0].BestMicroObservationIndex, 0);
                Access::Notify(Forward, 0);
                const auto Actions = Forward.GetAudioStateChangeActions();
                ASSERT_EQ(Actions.size(), 1U);
                EXPECT_FLOAT_EQ(std::get<float>(Actions[0].Args[1]),
                                static_cast<float>(State.MicroStates[Best].Observations[0].Freq));
            } else {
                EXPECT_NEAR(Likelihood[Reverse], Product, 1e-14);
            }
        }
        if (Topology == LEFT_RIGHT) {
            EXPECT_GT(Likelihood[0], Likelihood[1]);
        } else {
            EXPECT_NEAR(Likelihood[0], Likelihood[1], 1e-14);
        }
    }
}

TEST(SemiMarkovSegments, OrderedLikelihoodCombinesTinyEntryBeforeExponentiation) {
    OnlineForward Forward;
    Access::Setup(Forward, {Chain(GLISS, 1, 2.0), Chain(GLISS, 2, 2.0)});
    Access::Frame(Forward, 0);
    Access::Emissions(Forward, 1, {0.5, 0.5});
    Access::SemiMarkov(Forward, 1);
    // Model a heavily scaled past age hypothesis without converting it to linear space.
    Forward.GetStates()[1].MicroStates[0].LogForwardByAge[1] = 720.0;
    Access::Frame(Forward, 1);
    Access::Emissions(Forward, 1, {0.5, 0.5});
    Forward.GetStates()[0].ExitProb[0] = 0.0;
    Forward.GetStates()[1].InitProb = std::exp(-720.0);
    Access::SemiMarkov(Forward, 1);
    EXPECT_NEAR(Access::SegmentLikelihood(Forward, 2, std::exp(-720.0)), 0.5, 1e-10);
    EXPECT_NEAR(Forward.GetStates()[1].Forward[1], Access::Occupancy(Forward, 1, 2, true) * 0.5, 1e-10);
    EXPECT_NEAR(Forward.GetStates()[1].ExitProb[1], Access::Occupancy(Forward, 1, 2, false) * 0.5, 1e-10);
    EXPECT_TRUE(std::isfinite(Forward.GetStates()[1].Forward[1]));
    EXPECT_TRUE(std::isfinite(Forward.GetStates()[1].ExitProb[1]));
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
TEST(MicroStateIntegration, LuaExportsPhasesAndLabelsWithoutDurationWeights) {
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
        assert(s[3].microstates[1].duration_weight == nil)
        assert(#s[4].microstates == 1 and s[4].micro_topology == 1)
        assert(s[4].microstates[1].observations[2].label == 'aeolian')
    )")) << Scofo.LuaGetError();
}
#endif

} // namespace
