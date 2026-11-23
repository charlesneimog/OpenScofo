#include "../Trainer.hpp"
#include "TestData.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <set>

using Trainer = OpenScofo::Trainer;
using Configuration = OpenScofo::CatConfiguration;
using Dataset = OpenScofo::Dataset;
using SampleOrigin = OpenScofo::SampleOrigin;

double EstimateFrequency(const std::vector<float> &samples, int sampleRate) {
    const std::size_t begin = samples.size() / 4;
    const std::size_t end = samples.size() * 3 / 4;
    std::size_t risingCrossings = 0;
    for (std::size_t i = begin + 1; i < end; ++i) {
        if (samples[i - 1] <= 0.0F && samples[i] > 0.0F) {
            ++risingCrossings;
        }
    }
    return static_cast<double>(risingCrossings) * static_cast<double>(sampleRate) / static_cast<double>(end - begin);
}

TEST(TrainingAudio, LoadsWavAiffAndResamples) {
    Trainer trainer;
    TrainingTempDirectory directory;
    const int sampleRate = 8000;
    const std::vector<float> tone = TestTone(sampleRate, 800, 440.0, 1);
    WriteTestAudio(directory.Path / "tone.wav", tone, sampleRate);
    WriteTestAudio(directory.Path / "tone.aiff", tone, sampleRate, true);
    WriteTestAudio(directory.Path / "tone.aif", tone, sampleRate, true);

    EXPECT_EQ(trainer.LoadMono(directory.Path / "tone.wav", sampleRate).size(), tone.size());
    EXPECT_EQ(trainer.LoadMono(directory.Path / "tone.aiff", sampleRate).size(), tone.size());
    EXPECT_EQ(trainer.LoadMono(directory.Path / "tone.aif", sampleRate).size(), tone.size());
    EXPECT_NEAR(trainer.LoadMono(directory.Path / "tone.wav", sampleRate * 2).size(), tone.size() * 2, 1);
    EXPECT_THROW(trainer.LoadMono(directory.Path / "tone.mp3", sampleRate), std::invalid_argument);
}

TEST(TrainingAudio, AugmentationsPreserveExpectedDimensions) {
    Trainer trainer;
    const std::vector<float> tone = TestTone(16000, 4096, 440.0, 2);
    std::mt19937 random(80);
    const std::vector<float> noise = trainer.AddGaussianNoise(tone, 0.01, random);
    const std::vector<float> stretched = trainer.TimeStretch(tone, 0.8, 512, 128);
    const std::vector<float> reverberated = trainer.Convolve(tone, {1.0F, 0.5F, 0.25F});

    EXPECT_EQ(noise.size(), tone.size());
    EXPECT_NE(noise, tone);
    EXPECT_EQ(stretched.size(), 5120U);
    EXPECT_EQ(reverberated.size(), tone.size());
}

TEST(TrainingAudio, ConvolutionPreservesAUnitImpulseInput) {
    Trainer trainer;
    const std::vector<float> input = {0.25F, -0.5F, 1.0F, 0.125F};
    const std::vector<float> output = trainer.Convolve(input, {1.0F});
    ASSERT_EQ(output.size(), input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        EXPECT_NEAR(output[i], input[i], 1.0e-5F);
    }
}

TEST(TrainingAudio, PitchShiftChangesPitchAndPreservesLength) {
    Trainer trainer;
    const int sampleRate = 16000;
    const std::vector<float> tone = TestTone(sampleRate, 16000, 440.0, 3);
    const std::vector<float> shifted = trainer.PitchShift(tone, 12.0, 1024, 256);
    ASSERT_EQ(shifted.size(), tone.size());
    EXPECT_NEAR(EstimateFrequency(shifted, sampleRate), 880.0, 35.0);
}

TEST(TrainingApi, UsesRecordingLevelPartitionsAndCache) {
    TrainingTempDirectory directory;
    const int sampleRate = 8000;
    CreateTrainingDataset(directory.Path / "dataset", sampleRate);

    Configuration configuration;
    configuration.BasePath = directory.Path;
    configuration.DatasetFolder = "dataset";
    // configuration.CachePath = "native.cache";
    configuration.SampleRate = sampleRate;
    configuration.TrainingSampleRates = {sampleRate};
    configuration.FFTSize = 512;
    configuration.HopSize = 256;
    configuration.MinDB = -100.0;
    configuration.MaxPitchConfidence = 1.0;
    configuration.TestFraction = 0.34;
    configuration.Descriptors = {"rms", "mfcc", "chroma"};

    Trainer trainer(configuration);
    trainer.Analyze();
    const Dataset &dataset = trainer.GetDataset();
    ASSERT_FALSE(dataset.Train.Empty());
    ASSERT_FALSE(dataset.Validation.Empty());
    ASSERT_FALSE(dataset.Test.Empty());
    // EXPECT_FALSE(dataset.LoadedFromCache);
    EXPECT_EQ(dataset.Labels, (std::vector<std::string>{"high", "low"}));
    EXPECT_EQ(dataset.Train.FeatureCount, 26U);

    std::size_t largestOriginalClass = 0;
    for (const auto &[label, count] : dataset.OriginalTrainCounts) {
        (void)label;
        largestOriginalClass = std::max(largestOriginalClass, count);
    }
    for (const std::string &label : dataset.Labels) {
        EXPECT_EQ(dataset.OriginalTrainCounts.at(label) + dataset.AugmentedTrainCounts.at(label), largestOriginalClass);
    }

    std::set<std::uint64_t> train;
    std::set<std::uint64_t> validation;
    std::set<std::uint64_t> test;
    for (const SampleOrigin &origin : dataset.Train.Origins)
        train.insert(origin.Recording);
    for (const SampleOrigin &origin : dataset.Validation.Origins) {
        EXPECT_FALSE(origin.Augmented);
        validation.insert(origin.Recording);
    }
    for (const SampleOrigin &origin : dataset.Test.Origins) {
        EXPECT_FALSE(origin.Augmented);
        test.insert(origin.Recording);
    }
    for (std::uint64_t id : train) {
        EXPECT_FALSE(validation.contains(id));
        EXPECT_FALSE(test.contains(id));
    }
    for (std::uint64_t id : validation)
        EXPECT_FALSE(test.contains(id));

    Trainer cached(configuration);
    cached.Analyze();
    // EXPECT_TRUE(cached.GetDataset().LoadedFromCache);
    EXPECT_EQ(cached.GetDataset().Train.SampleCount(), dataset.Train.SampleCount());

    Configuration changed = configuration;
    changed.MinDB = -99.0;
    Trainer invalidated(changed);
    invalidated.Analyze();
    // EXPECT_FALSE(invalidated.GetDataset().LoadedFromCache);

    invalidated.ClearData();
    EXPECT_FALSE(std::filesystem::exists(directory.Path / "native.cache"));
}
