#include <gtest/gtest.h>
#include <OpenScofo.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <numbers>
#include <tuple>

namespace OpenScofo {

struct MIRConfigurationTestAccess {
    static void CheckOnset(MIR &Mir) {
        ASSERT_NE(Mir.m_ODS, nullptr);
        EXPECT_EQ(Mir.m_ODS->odftype, ODS_ODF_COMPLEX);
        EXPECT_EQ(Mir.m_ODS->medspan, 31U);
        EXPECT_EQ(Mir.m_ODS->fftsize, 4096U);
        EXPECT_FLOAT_EQ(Mir.m_ODS->srate, 44100.0f);
    }
};

// This executable is separate from microstate_tests, which has its own accessor.
struct OnlineForwardTestAccess {
    static void CheckTiming(const OnlineForward &Forward) {
        EXPECT_DOUBLE_EQ(Forward.m_SyncStrength, 0.75);
        EXPECT_DOUBLE_EQ(Forward.m_PhaseCoupling, 1.25);
    }

    static double Occupancy(OnlineForward &Forward, ScoreState State, int Age) {
        return Forward.GetOccupancyDistribution(State, Age);
    }

    static void Notify(OnlineForward &Forward, int Index) {
        Forward.m_States[Index].BestAudioStateIndex = 0;
        Forward.NotifyAudioStateChange(Index);
    }

    static std::pair<double, double> AdvanceTiming(OnlineForward &Forward) {
        // A deterministic late entrance, with a nonzero preceding phase error.
        Forward.m_CurrentStateIndex = 1;
        Forward.m_CurrentStateOnset = 0.0;
        Forward.m_TimeInPrevEvent = 0.8;
        Forward.m_LastTn = 0.0;
        Forward.m_States[1].IOIPhiN = 0.15;
        Forward.m_States[1].IOIHatPhiN = 0.0;
        Forward.UpdatePsiN(2);
        return {Forward.m_States[2].PhaseObserved, Forward.m_PsiN1};
    }
};

} // namespace OpenScofo

namespace {

const std::filesystem::path ScorePath = std::filesystem::path(TEST_DATA_DIR) / "configuration.scofo";

void ExpectOverrides(const OpenScofo::Configuration &Config) {
    const OpenScofo::Configuration Defaults;
    // Check both the parsed result and that the fixture really overrides every default.
#define CHECK_OVERRIDE(Field, Value)                                                                                   \
    EXPECT_EQ(Config.Field, Value) << #Field;                                                                          \
    EXPECT_NE(Config.Field, Defaults.Field) << #Field << " was not overridden"
    CHECK_OVERRIDE(SectionRestrict, true);
    CHECK_OVERRIDE(SR, 44100);
    CHECK_OVERRIDE(FFTSize, 4096);
    CHECK_OVERRIDE(HOPSize, 256);
    CHECK_OVERRIDE(TuningA4, 442);
    CHECK_OVERRIDE(PitchTemplateSigma, 0.75f);
    CHECK_OVERRIDE(PitchTemplateHarmonics, 6);
    CHECK_OVERRIDE(MFCCMels, 32);
    CHECK_OVERRIDE(MFCCCount, 16);
    CHECK_OVERRIDE(OnsetType, ODS_ODF_COMPLEX);
    CHECK_OVERRIDE(MedSpan, 31);
    CHECK_OVERRIDE(dBTreshold, -48);
    CHECK_OVERRIDE(SpectralRolloffCutoff, 0.65);
    CHECK_OVERRIDE(YINThreshold, 0.25);
    CHECK_OVERRIDE(YINMinFrequency, 80);
    CHECK_OVERRIDE(YINMaxFrequency, 1200);
    CHECK_OVERRIDE(ChromaSize, 24);
    CHECK_OVERRIDE(ChromaCenterOctave, 4);
    CHECK_OVERRIDE(ChromaOctaveWidth, 1.5);
    CHECK_OVERRIDE(ZCRCenter, false);
    CHECK_OVERRIDE(ZCRPad, true);
    CHECK_OVERRIDE(ZCRZeroPos, false);
    CHECK_OVERRIDE(ZCRThreshold, 0.002);
    CHECK_OVERRIDE(SyncStrength, 0.75f);
    CHECK_OVERRIDE(PhaseCoupling, 1.25f);
    CHECK_OVERRIDE(AudioStateChangeReceiver, "configuration_state");
    CHECK_OVERRIDE(TimbreONNXModel, ScorePath.parent_path() / "configuration.onnx");
    const std::vector<std::string> Descriptors = {"rms", "centroid"};
    CHECK_OVERRIDE(ONNXDescriptors, Descriptors);
#undef CHECK_OVERRIDE
}

class ScoreConfiguration : public ::testing::Test {
  protected:
    OpenScofo::Configuration Config;
    OpenScofo::States States;

    void SetUp() override {
        ASSERT_TRUE(std::filesystem::is_regular_file(ScorePath));
        OpenScofo::Score Parser;
        std::tie(Config, States) = Parser.Parse(ScorePath);
        ASSERT_FALSE(States.empty());
    }
};

std::vector<double> Tone(const OpenScofo::Configuration &Config, double Frequency, double Amplitude = 0.4) {
    std::vector<double> Audio(static_cast<size_t>(Config.FFTSize));
    for (size_t I = 0; I < Audio.size(); ++I) {
        Audio[I] = Amplitude * std::sin(2.0 * std::numbers::pi * Frequency * I / Config.SR);
    }
    return Audio;
}

OpenScofo::Description Describe(OpenScofo::Configuration Config, const std::vector<double> &Audio) {
    Config.RequestedDescriptors = {OpenScofo::MFCC, OpenScofo::LOGMEL, OpenScofo::CHROMA,
                                   OpenScofo::YIN,  OpenScofo::ZCR,    OpenScofo::ODSONSET};
    OpenScofo::OpenScofo Scofo(Config.SR, Config.FFTSize, Config.HOPSize);
    Scofo.SetConfiguration(Config);
    EXPECT_TRUE(Scofo.ProcessBlock(Audio.data(), Audio.size()));
    return Scofo.GetDescription();
}

TEST_F(ScoreConfiguration, ParsesEveryOverrideAndAppliesPitchAndEventSettings) {
    ExpectOverrides(Config);
    bool SawNote = false;
    for (const auto &State : States) {
        SCOPED_TRACE(State.Index);
        EXPECT_DOUBLE_EQ(State.BPMExpected, 90.0);
        EXPECT_DOUBLE_EQ(State.SyncStrength, 0.75);
        EXPECT_DOUBLE_EQ(State.PhaseCoupling, 1.25);
        EXPECT_DOUBLE_EQ(State.TimeTolerance, 49.0);
        if (State.Type == OpenScofo::NOTE && State.Section == "configured") {
            ASSERT_EQ(State.Observations.size(), 1U);
            EXPECT_DOUBLE_EQ(State.Observations[0].Midi, 81.0);  // A4 + 12 semitones
            EXPECT_DOUBLE_EQ(State.Observations[0].Freq, 884.0); // A4 = 442 Hz
            SawNote = true;
        }
    }
    EXPECT_TRUE(SawNote);
}

TEST_F(ScoreConfiguration, LoadScoreReconfiguresAnalysisAndRunsModel) {
    OpenScofo::OpenScofo Scofo(44100, 2048, 512);
    Scofo.SetRequestedDescriptors(
        {OpenScofo::MFCC, OpenScofo::LOGMEL, OpenScofo::CHROMA, OpenScofo::YIN, OpenScofo::ZCR, OpenScofo::ODSONSET});
    ASSERT_TRUE(Scofo.LoadScore(ScorePath));
    ExpectOverrides(Scofo.GetConfiguration());
    EXPECT_EQ(Scofo.GetSr(), 44100);
    EXPECT_EQ(Scofo.GetFFTSize(), 4096);
    EXPECT_EQ(Scofo.GetHopSize(), 256);
    EXPECT_DOUBLE_EQ(Scofo.GetBlockDuration(), 256.0 / 44100.0);
    EXPECT_DOUBLE_EQ(Scofo.GetCurrentBPM(), 90.0);
    const auto Audio = Tone(Config, 884);
    ASSERT_TRUE(Scofo.ProcessBlock(Audio.data(), Audio.size()));
    const auto Desc = Scofo.GetDescription();
    EXPECT_EQ(Desc.Power.size(), 2049U);
    EXPECT_EQ(Desc.MFCC.size(), 16U);
    EXPECT_EQ(Desc.LogMelSpectrum.size(), 32U);
    EXPECT_EQ(Desc.Chroma.size(), 24U);
    EXPECT_NEAR(Desc.Pitch, 884.0, 2.0);
    ASSERT_TRUE(Desc.ONNX.contains("loud"));
    EXPECT_GT(Desc.ONNX.at("loud"), 0.99f);
    ASSERT_TRUE(Scofo.SetCurrentSection("other"));
    EXPECT_EQ(Scofo.GetStates()[Scofo.GetCurrentStateIndex()].Section, "other");
}

TEST_F(ScoreConfiguration, ModelUsesScoreDescriptorOrderAndBothLabels) {
    OpenScofo::ONNXModel Model;
    ASSERT_TRUE(Model.Load(Config.TimbreONNXModel, {OpenScofo::RMS, OpenScofo::CENTROID}, Config));
    EXPECT_EQ(Model.GetDescriptors(), (std::vector<OpenScofo::Descriptors>{OpenScofo::RMS, OpenScofo::CENTROID}));
    EXPECT_EQ(Model.GetLabels(), (std::vector<std::string>{"quiet", "loud"}));
    OpenScofo::Description Desc{};
    Desc.SpectralCentroid = 884;
    Desc.RMS = 0.01;
    Model.Execute(Desc);
    ASSERT_TRUE(Desc.ONNX.contains("quiet"));
    EXPECT_GT(Desc.ONNX.at("quiet"), 0.99f);
    Desc.RMS = 0.4;
    Model.Execute(Desc);
    ASSERT_TRUE(Desc.ONNX.contains("loud"));
    EXPECT_GT(Desc.ONNX.at("loud"), 0.99f);
}

TEST_F(ScoreConfiguration, ConfiguresActualOnsetDetector) {
    Config.RequestedDescriptors = {OpenScofo::ODSONSET};
    OpenScofo::MIR Mir;
    Mir.UpdateConfiguration(Config);
    OpenScofo::MIRConfigurationTestAccess::CheckOnset(Mir);
}

TEST_F(ScoreConfiguration, HarmonicCountAndSigmaChangePitchTemplates) {
    OpenScofo::OnlineForward Forward;
    Forward.UpdateConfiguration(Config);
    const auto Configured = Forward.GetPitchTemplate(884);
    ASSERT_EQ(Configured.size(), 2048U);
    auto Changed = Config;
    Changed.PitchTemplateHarmonics = 1;
    Forward.UpdateConfiguration(Changed);
    EXPECT_NE(Configured, Forward.GetPitchTemplate(884));
    Changed = Config;
    Changed.PitchTemplateSigma = 0.1;
    Forward.UpdateConfiguration(Changed);
    EXPECT_NE(Configured, Forward.GetPitchTemplate(884));
    EXPECT_EQ(Forward.GetTunning(), 442);
}

TEST_F(ScoreConfiguration, TimingSettingsReachFollowerAndChangePhaseAndTempo) {
    OpenScofo::OnlineForward Forward;
    Forward.UpdateConfiguration(Config);
    Forward.SetScoreStates(States);
    ASSERT_TRUE(Forward.SetCurrentSection("configured"));
    OpenScofo::OnlineForwardTestAccess::CheckTiming(Forward);
    const auto Configured = OpenScofo::OnlineForwardTestAccess::AdvanceTiming(Forward);
    for (bool ChangePhase : {true, false}) {
        auto ChangedStates = States;
        for (auto &State : ChangedStates) {
            if (ChangePhase)
                State.PhaseCoupling = 0;
            else
                State.SyncStrength = 0;
        }
        OpenScofo::OnlineForward Other;
        Other.UpdateConfiguration(Config);
        Other.SetScoreStates(ChangedStates);
        const auto Changed = OpenScofo::OnlineForwardTestAccess::AdvanceTiming(Other);
        if (ChangePhase)
            EXPECT_NE(Configured.first, Changed.first);
        else
            EXPECT_NE(Configured.second, Changed.second);
    }
}

TEST_F(ScoreConfiguration, AudioStateReceiverIsUsedForNotifications) {
    OpenScofo::OnlineForward Forward;
    Forward.UpdateConfiguration(Config);
    Forward.SetScoreStates(States);
    const auto Note =
        std::find_if(States.begin(), States.end(), [](const auto &State) { return State.Type == OpenScofo::NOTE; });
    ASSERT_NE(Note, States.end());
    OpenScofo::OnlineForwardTestAccess::Notify(Forward, static_cast<int>(Note - States.begin()));
    const auto Actions = Forward.GetAudioStateChangeActions();
    ASSERT_EQ(Actions.size(), 1U);
    EXPECT_EQ(Actions.front().Receiver, "configuration_state");
}

TEST_F(ScoreConfiguration, MelAndChromaSettingsChangeComputedFeatures) {
    const auto Audio = Tone(Config, 884);
    const auto Configured = Describe(Config, Audio);
    auto Changed = Config;
    Changed.MFCCMels = 40;
    const auto Mel = Describe(Changed, Audio);
    EXPECT_EQ(Mel.LogMelSpectrum.size(), 40U);
    EXPECT_NE(Configured.MFCC, Mel.MFCC);
    Changed = Config;
    Changed.MFCCCount = 13;
    EXPECT_EQ(Describe(Changed, Audio).MFCC.size(), 13U);
    Changed = Config;
    Changed.ChromaSize = 12;
    EXPECT_EQ(Describe(Changed, Audio).Chroma.size(), 12U);
    Changed = Config;
    Changed.ChromaCenterOctave = 5;
    EXPECT_NE(Configured.Chroma, Describe(Changed, Audio).Chroma);
    Changed = Config;
    Changed.ChromaOctaveWidth = 2;
    EXPECT_NE(Configured.Chroma, Describe(Changed, Audio).Chroma);
}

TEST_F(ScoreConfiguration, RolloffCutoffChangesReportedFrequency) {
    auto Audio = Tone(Config, 300);
    const auto High = Tone(Config, 3000, 0.3);
    for (size_t I = 0; I < Audio.size(); ++I)
        Audio[I] += High[I];
    auto Changed = Config;
    Changed.SpectralRolloffCutoff = 0.95;
    EXPECT_LT(Describe(Config, Audio).SpectralRolloff, Describe(Changed, Audio).SpectralRolloff);
}

TEST_F(ScoreConfiguration, YinFrequencyLimitsConstrainPitchEstimates) {
    const auto Audio = Tone(Config, 440);
    EXPECT_NEAR(Describe(Config, Audio).Pitch, 440, 2);
    auto Changed = Config;
    Changed.YINMinFrequency = 700;
    const auto LowLimit = Describe(Changed, Audio).Pitch;
    EXPECT_TRUE(LowLimit == 0 || LowLimit >= 700);
    Changed = Config;
    Changed.YINMaxFrequency = 300;
    const auto HighLimit = Describe(Changed, Audio).Pitch;
    EXPECT_TRUE(HighLimit == 0 || HighLimit <= 300);
}

TEST_F(ScoreConfiguration, YinThresholdChangesWhichPeriodIsSelected) {
    // A weak fundamental and strong second harmonic have two candidate minima.
    auto Audio = Tone(Config, 220, 0.1);
    const auto Harmonic = Tone(Config, 440, 0.4);
    for (size_t I = 0; I < Audio.size(); ++I) {
        Audio[I] += Harmonic[I];
    }
    auto Changed = Config;
    Changed.YINThreshold = 0.01;
    EXPECT_NEAR(Describe(Config, Audio).Pitch, 440, 10);
    EXPECT_NEAR(Describe(Changed, Audio).Pitch, 220, 2);
}

TEST_F(ScoreConfiguration, HopSizeControlsWhenAudioAnalysisUpdates) {
    OpenScofo::OpenScofo Scofo(Config.SR, Config.FFTSize, Config.HOPSize);
    Scofo.SetConfiguration(Config);
    const auto Audio = Tone(Config, 440);
    // Feed whole hops so the update counter has no accumulated remainder.
    for (size_t Offset = 0; Offset < Audio.size(); Offset += 256) {
        ASSERT_TRUE(Scofo.ProcessBlock(Audio.data() + Offset, 256));
    }
    const double PreviousRMS = Scofo.GetDescription().RMS;
    ASSERT_GT(PreviousRMS, 0);
    const std::vector<double> Silence(256, 0.0);
    ASSERT_TRUE(Scofo.ProcessBlock(Silence.data(), 255));
    EXPECT_DOUBLE_EQ(Scofo.GetDescription().RMS, PreviousRMS);
    ASSERT_TRUE(Scofo.ProcessBlock(Silence.data(), 1));
    EXPECT_LT(Scofo.GetDescription().RMS, PreviousRMS);
}

TEST_F(ScoreConfiguration, ZeroCrossingOptionsAffectCrossingCount) {
    std::vector<double> Audio(4096);
    for (size_t I = 0; I < Audio.size(); ++I)
        Audio[I] = I % 2 == 0 ? 0.001 : -0.01;
    // With zero treated as a separate sign, every pair crosses, plus initial padding.
    EXPECT_DOUBLE_EQ(Describe(Config, Audio).ZeroCrossingRate, 1.0);
    auto Changed = Config;
    Changed.ZCRPad = false;
    EXPECT_DOUBLE_EQ(Describe(Changed, Audio).ZeroCrossingRate, 4095.0 / 4096.0);
    Changed = Config;
    Changed.ZCRCenter = true;
    EXPECT_LT(Describe(Changed, Audio).ZeroCrossingRate, 1.0);
    // Both tiny positive and negative values become zero at a higher threshold.
    Changed = Config;
    Changed.ZCRThreshold = 0.02;
    EXPECT_DOUBLE_EQ(Describe(Changed, Audio).ZeroCrossingRate, 1.0 / 4096.0);
    for (size_t I = 0; I < Audio.size(); ++I)
        Audio[I] = I % 2 == 0 ? 0.001 : 0.01;
    EXPECT_DOUBLE_EQ(Describe(Config, Audio).ZeroCrossingRate, 1.0);
    Changed = Config;
    Changed.ZCRZeroPos = true;
    EXPECT_DOUBLE_EQ(Describe(Changed, Audio).ZeroCrossingRate, 1.0 / 4096.0);
}

// These are characterization tests for documented gaps, NOT evidence that the
// settings work. Replace with behavior tests when the runtime wiring is fixed.
TEST_F(ScoreConfiguration, KnownLimitationDbThresholdDoesNotAffectSilenceProbability) {
    const auto Audio = Tone(Config, 440, 0.0001);
    const auto Original = Describe(Config, Audio);
    auto Changed = Config;
    Changed.dBTreshold = -20;
    const auto Other = Describe(Changed, Audio);
    EXPECT_DOUBLE_EQ(Original.SilenceProb, Other.SilenceProb);
    EXPECT_NEAR(Original.SilenceProb, 1.0 / (1.0 + std::exp(0.25 * (Original.Loudness + 60))), 1e-12);
}

TEST_F(ScoreConfiguration, KnownLimitationTimeToleranceDoesNotAffectDurationDistribution) {
    OpenScofo::OnlineForward Forward;
    Forward.UpdateConfiguration(Config);
    Forward.SetScoreStates(States);
    auto State = States[1];
    State.TimeTolerance = 4;
    auto Other = State;
    Other.TimeTolerance = 64;
    for (int Age : {1, 50, 100, 200}) {
        EXPECT_DOUBLE_EQ(OpenScofo::OnlineForwardTestAccess::Occupancy(Forward, State, Age),
                         OpenScofo::OnlineForwardTestAccess::Occupancy(Forward, Other, Age));
    }
}

} // namespace
