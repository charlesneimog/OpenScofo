#include <gtest/gtest.h>
#include <onnx.h>

// clang-format off
TEST(OnnxMetadata, OnnxMetadata) {
    fprintf(stderr, "before context allocation\n");
    struct onnx_context_t *ctx = onnx_context_alloc_from_file("/home/neimog/Documents/Git/OpenScofo/Tests/models/flute.onnx", NULL, 0);
    fprintf(stderr, "after context allocation: %p\n", (void *)ctx);

    if (!ctx) {
        return ;
    }

    const char *keys[] = 
        {
            "openscofo.sample_rate", 
            "openscofo.fft_size", 
            "openscofo.hop_size", 
            "openscofo.descriptors",
            "openscofo.labels"
        };

    for (int i = 0; i < 5; ++i) {
        const char *s = onnx_metadata_get(ctx, keys[i]);
        fprintf(stderr, "%s: %s\n", keys[i], s ? s : "(absent)");
    }

    onnx_context_free(ctx);
    return ;
}
