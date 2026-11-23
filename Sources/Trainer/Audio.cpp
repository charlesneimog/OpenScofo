/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "Trainer.hpp"

#include <AudioFile.h>
#include <pffft/pffft.h>
#include <signalsmith-stretch/signalsmith-stretch.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <new>
#include <stdexcept>
#include <string>

namespace OpenScofo {

// ─────────────────────────────────────
RealFft::RealFft(int size)
    : m_Size(size), m_Setup(pffft_new_setup(size, PFFFT_REAL), pffft_destroy_setup), m_Work(AllocateBuffer()) {
    if (m_Setup == nullptr) {
        throw std::runtime_error("PFFFT does not support FFT size " + std::to_string(size));
    }
}

// ─────────────────────────────────────
std::unique_ptr<float, void (*)(void *)> RealFft::AllocateBuffer() const {
    float *buffer = static_cast<float *>(pffft_aligned_malloc(static_cast<std::size_t>(m_Size) * sizeof(float)));
    if (buffer == nullptr) {
        throw std::bad_alloc();
    }
    return {buffer, pffft_aligned_free};
}

// ─────────────────────────────────────
PFFFT_Setup *RealFft::GetSetup() {
    return m_Setup.get();
}

// ─────────────────────────────────────
void RealFft::Forward(const float *input, float *output) {
    pffft_transform(m_Setup.get(), input, output, m_Work.get(), PFFFT_FORWARD);
}

// ─────────────────────────────────────
void RealFft::Backward(const float *input, float *output) {
    pffft_transform(m_Setup.get(), input, output, m_Work.get(), PFFFT_BACKWARD);
}

// ─────────────────────────────────────
AudioBuffer Trainer::LoadAudio(const std::filesystem::path &path) const {
    const std::string pathString = path.string();
    if (!IsAudioFile(path)) {
        spdlog::error("Unsupported audio file format '{}': {}", path.extension().string(), pathString);
        return {};
    }

    AudioFile<float> file;
    file.shouldLogErrorsToConsole(false);
    if (!file.load(pathString)) {
        throw std::runtime_error("Could not read WAV/AIFF audio file: " + pathString);
    }

    const int channels = file.getNumChannels();
    const int frames = file.getNumSamplesPerChannel();
    const std::uint32_t sampleRate = file.getSampleRate();
    if (frames <= 0 || channels <= 0 || file.samples.size() != static_cast<std::size_t>(channels) || sampleRate == 0 ||
        sampleRate > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("Audio file has invalid dimensions: " + pathString);
    }

    for (const std::vector<float> &channel : file.samples) {
        if (channel.size() < static_cast<std::size_t>(frames)) {
            throw std::runtime_error("Audio file has inconsistent channel lengths: " + pathString);
        }
    }

    AudioBuffer result;
    result.SampleRate = static_cast<int>(sampleRate);
    result.Samples.resize(static_cast<std::size_t>(frames));
    const float channelScale = 1.0F / channels;

    for (std::size_t frame = 0; frame < result.Samples.size(); ++frame) {
        float sample = 0.0F;
        for (int channel = 0; channel < channels; ++channel) {
            sample += file.samples[channel][frame];
        }
        result.Samples[frame] = sample * channelScale;
    }
    return result;
}

// ─────────────────────────────────────
std::vector<float> Trainer::Resample(const std::vector<float> &input, std::size_t outputSize) const {
    if (input.empty() || outputSize == 0) {
        return {};
    }
    if (input.size() == outputSize) {
        return input;
    }
    if (input.size() == 1) {
        return std::vector<float>(outputSize, input.front());
    }
    if (outputSize == 1) {
        return {input.front()};
    }

    std::vector<float> output(outputSize);
    const double scale = static_cast<double>(input.size() - 1) / static_cast<double>(outputSize - 1);
    for (std::size_t i = 0; i < outputSize; ++i) {
        const double position = static_cast<double>(i) * scale;
        const std::size_t left = static_cast<std::size_t>(position);
        const std::size_t right = std::min(left + 1, input.size() - 1);
        const double fraction = position - static_cast<double>(left);
        output[i] = static_cast<float>(input[left] + (input[right] - input[left]) * fraction);
    }
    return output;
}

// ─────────────────────────────────────
std::vector<float> Trainer::LoadMono(const std::filesystem::path &path, int sampleRate) const {
    if (sampleRate <= 0) {
        spdlog::error("Training sample rates must be positive");
        return {};
    }
    AudioBuffer source = LoadAudio(path);
    if (source.Samples.empty()) {
        return {};
    }
    if (source.SampleRate == sampleRate) {
        return std::move(source.Samples);
    }

    const double ratio = static_cast<double>(sampleRate) / static_cast<double>(source.SampleRate);
    const std::size_t outputSize = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::llround(static_cast<double>(source.Samples.size()) * ratio)));
    return Resample(source.Samples, outputSize);
}

// ─────────────────────────────────────
std::vector<float> Trainer::TimeStretch(const std::vector<float> &input, double rate, int fftSize, int hopSize) const {
    if (input.empty()) {
        return {};
    }
    if (!std::isfinite(rate) || rate <= 0.0) {
        spdlog::error("Time-stretch rate must be positive");
        return {};
    }
    if (fftSize < 32 || hopSize <= 0 || hopSize > fftSize) {
        spdlog::error("Invalid FFT or hop size for time stretching");
        return {};
    }

    const double outputSizeValue = std::round(static_cast<double>(input.size()) / rate);
    if (outputSizeValue > static_cast<double>(std::numeric_limits<int>::max()) ||
        input.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        spdlog::error("Time-stretch input or output exceeds Signalsmith Stretch's supported size");
        return {};
    }

    const std::size_t outputSize = std::max<std::size_t>(1, static_cast<std::size_t>(outputSizeValue));
    std::vector<float> output(outputSize);
    const float *inputBuffers[] = {input.data()};
    float *outputBuffers[] = {output.data()};

    signalsmith::stretch::SignalsmithStretch<float> stretch;
    stretch.configure(1, fftSize, hopSize);
    if (!stretch.exact(inputBuffers, static_cast<int>(input.size()), outputBuffers, static_cast<int>(output.size()))) {
        spdlog::warn("Signalsmith Stretch could not time-stretch the signal; using linear resampling instead");
        return Resample(input, outputSize);
    }
    return output;
}

// ─────────────────────────────────────
std::vector<float> Trainer::PitchShift(const std::vector<float> &input, double semitones, int fftSize,
                                       int hopSize) const {
    if (input.empty()) {
        return {};
    }
    if (!std::isfinite(semitones)) {
        spdlog::error("Pitch-shift semitones must be finite");
        return {};
    }
    if (fftSize < 32 || hopSize <= 0 || hopSize > fftSize) {
        spdlog::error("Invalid FFT or hop size for pitch shifting");
        return {};
    }
    if (input.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        spdlog::error("Pitch-shift input exceeds Signalsmith Stretch's supported size");
        return {};
    }

    std::vector<float> output(input.size());
    const float *inputBuffers[] = {input.data()};
    float *outputBuffers[] = {output.data()};

    signalsmith::stretch::SignalsmithStretch<float> stretch;
    stretch.configure(1, fftSize, hopSize);
    stretch.setTransposeSemitones(static_cast<float>(semitones));
    if (!stretch.exact(inputBuffers, static_cast<int>(input.size()), outputBuffers, static_cast<int>(output.size()))) {
        spdlog::warn("Signalsmith Stretch could not pitch-shift the signal; returning it unchanged");
        return input;
    }
    return output;
}

// ─────────────────────────────────────
std::vector<float> Trainer::AddGaussianNoise(const std::vector<float> &input, double deviation,
                                             std::mt19937 &random) const {
    if (deviation < 0.0 || !std::isfinite(deviation)) {
        spdlog::error("Noise deviation must be finite and non-negative");
        return {};
    }
    std::normal_distribution<float> noise(0.0F, static_cast<float>(deviation));
    std::vector<float> output(input);
    for (float &sample : output) {
        sample += noise(random);
    }
    return output;
}

// ─────────────────────────────────────
std::vector<float> Trainer::Convolve(const std::vector<float> &input, const std::vector<float> &impulseResponse) const {
    if (input.empty() || impulseResponse.empty()) {
        return input;
    }

    float impulsePeak = 0.0F;
    for (float sample : impulseResponse) {
        impulsePeak = std::max(impulsePeak, std::abs(sample));
    }
    if (impulsePeak <= 1.0e-12F) {
        return input;
    }

    const std::size_t convolutionSize = input.size() + impulseResponse.size() - 1;
    std::size_t fftSize = static_cast<std::size_t>(pffft_min_fft_size(PFFFT_REAL));
    while (fftSize < convolutionSize) {
        fftSize <<= 1U;
    }
    if (fftSize > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        spdlog::error("Convolution of {} samples exceeds PFFFT's supported FFT size", convolutionSize);
        return {};
    }

    RealFft transform(static_cast<int>(fftSize));
    auto time = transform.AllocateBuffer();
    auto firstSpectrum = transform.AllocateBuffer();
    auto secondSpectrum = transform.AllocateBuffer();

    std::fill_n(time.get(), fftSize, 0.0F);
    std::copy(input.begin(), input.end(), time.get());
    transform.Forward(time.get(), firstSpectrum.get());

    std::fill_n(time.get(), fftSize, 0.0F);
    for (std::size_t i = 0; i < impulseResponse.size(); ++i) {
        time.get()[i] = impulseResponse[i] / impulsePeak;
    }
    transform.Forward(time.get(), secondSpectrum.get());

    pffft_zconvolve_no_accu(transform.GetSetup(), firstSpectrum.get(), secondSpectrum.get(), firstSpectrum.get(),
                            1.0F / static_cast<float>(fftSize));
    transform.Backward(firstSpectrum.get(), time.get());

    std::vector<float> result(input.size());
    float resultPeak = 0.0F;
    for (std::size_t i = 0; i < input.size(); ++i) {
        result[i] = time.get()[i];
        resultPeak = std::max(resultPeak, std::abs(result[i]));
    }
    const float scale = resultPeak > 1.0e-12F ? 1.0F / resultPeak : 1.0F;
    for (float &sample : result) {
        sample *= scale;
    }
    return result;
}

} // namespace OpenScofo
