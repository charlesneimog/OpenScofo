#include <gtest/gtest.h>
#include <onnx.h>

#include <filesystem>

TEST(OnnxMetadata, ReadsPortableConfigurationFixture) {
    const auto Path = std::filesystem::path(TEST_DATA_DIR) / "configuration.onnx";
    auto *Context = onnx_context_alloc_from_file(Path.string().c_str(), nullptr, 0);
    ASSERT_NE(Context, nullptr) << Path;
    EXPECT_STREQ(onnx_metadata_get(Context, "openscofo.sample_rate"), "44100");
    EXPECT_STREQ(onnx_metadata_get(Context, "openscofo.fft_size"), "4096");
    EXPECT_STREQ(onnx_metadata_get(Context, "openscofo.hop_size"), "256");
    EXPECT_STREQ(onnx_metadata_get(Context, "openscofo.labels"), "[\"quiet\", \"loud\"]");
    // Deliberately absent so configuration tests exercise ONNXDESCRIPTORS.
    EXPECT_EQ(onnx_metadata_get(Context, "openscofo.descriptors"), nullptr);
    onnx_context_free(Context);
}
