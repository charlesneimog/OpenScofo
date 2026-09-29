#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <pffft/pffft.h> // Ensure pffft.h and pffft.c are compiled with your project

struct PFFFTContext {
    int N;
    PFFFT_Setup *setup;
    float *in;
    float *out;
    float *work;
};

// Allocation and setup phase (Excluded from execution benchmark)
PFFFTContext fft_init(int N) {
    PFFFTContext ctx;
    ctx.N = N;

    ctx.setup = pffft_new_setup(N, PFFFT_REAL);

    // PFFFT requires strict 16-byte alignment for vector instructions
    ctx.in = (float *)pffft_aligned_malloc(N * sizeof(float));
    ctx.out = (float *)pffft_aligned_malloc(N * sizeof(float));
    ctx.work = (float *)pffft_aligned_malloc(N * sizeof(float)); // Scratch space

    return ctx;
}

// Pure execution phase (The target for benchmarking)
void fft_execute(PFFFTContext &ctx) {
    // Computes forward transform and natively reorders it back to normal frequency bins
    pffft_transform_ordered(ctx.setup, ctx.in, ctx.out, ctx.work, PFFFT_FORWARD);
}

void fft_cleanup(PFFFTContext &ctx) {
    pffft_aligned_free(ctx.in);
    pffft_aligned_free(ctx.out);
    pffft_aligned_free(ctx.work);

    // CORRECTED: The PFFFT API uses destroy, not free, for the setup struct.
    pffft_destroy_setup(ctx.setup);
}

void run_main(int N) {
    const int ITERATIONS = 100000;
    const float SAMPLING_FREQ = 8000.0f;
    const float TARGET_FREQ = 440.0f;

    // 1. Initialize PFFFT Plan
    PFFFTContext ctx = fft_init(N);

    // 2. Generate Input Signal (440 Hz Sinusoid)
    for (int i = 0; i < N; ++i) {
        ctx.in[i] = std::sin(2.0f * M_PI * TARGET_FREQ * i / SAMPLING_FREQ);
    }

    // 3. Benchmark the pure execution loop
    std::cout << "Running PFFFT Benchmark for " << ITERATIONS << " iterations..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < ITERATIONS; ++i) {
        fft_execute(ctx);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << "PFFFT Total Time: " << elapsed.count() << " ms" << std::endl;
    std::cout << "PFFFT Average Time per FFT: " << (elapsed.count() / ITERATIONS) * 1000.0 << " us" << std::endl;

    // Clean up resources
    fft_cleanup(ctx);
}

int main() {
    run_main(512);
    run_main(1024);
    run_main(2048);
    run_main(4096);
    return 0;
}
