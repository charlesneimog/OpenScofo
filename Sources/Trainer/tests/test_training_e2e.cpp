#include "../Trainer.hpp"
#include "TestData.hpp"

#include <OpenScofo.hpp>
#include <gtest/gtest.h>
#include <onnx.h>

using Trainer = OpenScofo::Trainer;
using Configuration = OpenScofo::CatConfiguration;

TEST(TrainingEndToEnd, ExtractsTrainsExportsAndRunsInOpenScofo) {
    TrainingTempDirectory directory;
    const int sampleRate = 8000;
    CreateTrainingDataset(directory.Path / "dataset", sampleRate);

    Configuration configuration;
    configuration.BasePath = directory.Path;
    configuration.DatasetFolder = "dataset";
    configuration.SampleRate = sampleRate;
    configuration.TrainingSampleRates = {sampleRate};
    configuration.FFTSize = 512;
    configuration.HopSize = 256;
    configuration.MinDB = -100.0;
    configuration.MaxPitchConfidence = 1.0;
    configuration.TestFraction = 0.34;
    configuration.Descriptors = {"rms", "centroid", "flatness"};
    configuration.Iterations = 30;
    configuration.Depth = 3;
    configuration.EarlyStoppingRounds = 5;
    configuration.ThreadCount = -1;

    Trainer trainer(configuration);
    ASSERT_TRUE(trainer.Analyze());
    ASSERT_TRUE(trainer.Train());
    ASSERT_TRUE(trainer.HasModel());
    EXPECT_GT(trainer.GetEvaluation().TestedSamples, 0U);

    const std::filesystem::path onnxPath = directory.Path / "model.onnx";
    const std::filesystem::path nativePath = directory.Path / "model.cbm";
    ASSERT_TRUE(trainer.ExportModel(onnxPath));
    ASSERT_TRUE(trainer.ExportModel(nativePath));
    ASSERT_TRUE(std::filesystem::exists(onnxPath));
    ASSERT_TRUE(std::filesystem::exists(nativePath));

    onnx_context_t *context = onnx_context_alloc_from_file(onnxPath.string().c_str(), nullptr, 0);
    ASSERT_NE(context, nullptr);
    EXPECT_STREQ(onnx_metadata_get(context, "openscofo.descriptors"), "[\"rms\",\"centroid\",\"flatness\"]");
    EXPECT_STREQ(onnx_metadata_get(context, "openscofo.labels"), "[\"high\",\"low\"]");
    EXPECT_STREQ(onnx_metadata_get(context, "openscofo.sample_rate"), "8000");
    EXPECT_STREQ(onnx_metadata_get(context, "openscofo.fft_size"), "512");
    EXPECT_STREQ(onnx_metadata_get(context, "openscofo.hop_size"), "256");
    EXPECT_NE(onnx_metadata_get(context, "openscofo.version"), nullptr);
    onnx_context_free(context);

    OpenScofo::OpenScofo inference(static_cast<float>(sampleRate), 512.0F, 256.0F);
    inference.LoadONNXModel(onnxPath, {OpenScofo::RMS, OpenScofo::CENTROID, OpenScofo::FLATNESS});
    inference.SetRequestedDescriptors({OpenScofo::RMS, OpenScofo::CENTROID, OpenScofo::FLATNESS, OpenScofo::ONNX});
    const std::vector<float> audio = TestTone(sampleRate, 512, 220.0, 88);
    ASSERT_TRUE(inference.ProcessBlock(audio.data(), 256));
    ASSERT_TRUE(inference.ProcessBlock(audio.data() + 256, 256));
    const OpenScofo::Description description = inference.GetDescription();
    EXPECT_EQ(description.ONNX.size(), 2U);
    EXPECT_TRUE(description.ONNX.contains("high"));
    EXPECT_TRUE(description.ONNX.contains("low"));
}
