#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fftw3.h>

// Global or context structure to hold the persistent plan and buffers
struct FFTWContext {
    int N;
    float *in;
    fftwf_complex *out;
    fftwf_plan plan;
};

// Allocation and planning phase (Excluded from execution benchmark)
FFTWContext fft_init(int N) {
    FFTWContext ctx;
    ctx.N = N;

    // FFTW-specific SIMD-aligned allocation
    ctx.in = fftwf_alloc_real(N);
    ctx.out = fftwf_alloc_complex(N / 2 + 1);

    // FFTW_MEASURE intensely analyzes hardware topology at startup
    ctx.plan = fftwf_plan_dft_r2c_1d(N, ctx.in, ctx.out, FFTW_MEASURE);

    return ctx;
}

// Pure execution phase (The target for benchmarking)
void fft_execute(FFTWContext &ctx) {
    fftwf_execute(ctx.plan);
}

void fft_cleanup(FFTWContext &ctx) {
    fftwf_destroy_plan(ctx.plan);
    fftwf_free(ctx.in);
    fftwf_free(ctx.out);
}

int main() {
    const int N = 2048;
    const int ITERATIONS = 1000000;
    const float SAMPLING_FREQ = 8000.0f;
    const float TARGET_FREQ = 440.0f;

    // 1. Initialize FFTW Plan
    FFTWContext ctx = fft_init(N);

    // 2. Generate Input Signal (440 Hz Sinusoid)
    for (int i = 0; i < N; ++i) {
        ctx.in[i] = std::sin(2.0f * M_PI * TARGET_FREQ * i / SAMPLING_FREQ);
    }

    // 3. Benchmark the pure execution loop
    std::cout << "Running FFTW3 (Float) Benchmark for " << ITERATIONS << " iterations..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < ITERATIONS; ++i) {
        fft_execute(ctx);
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << "FFTW3 Total Time: " << elapsed.count() << " ms" << std::endl;
    std::cout << "FFTW3 Average Time per FFT: " << (elapsed.count() / ITERATIONS) * 1000.0 << " us" << std::endl;

    // Clean up resources
    fft_cleanup(ctx);
    return 0;
}
