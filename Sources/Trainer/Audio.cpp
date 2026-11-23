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
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace OpenScofo {

// ─────────────────────────────────────
RealFft::RealFft(int size, std::shared_ptr<spdlog::logger> logger)
    : m_Size(size), m_Logger(std::move(logger)), m_Setup(pffft_new_setup(size, PFFFT_REAL), pffft_destroy_setup),
      m_Work(AllocateBuffer()) {
    if (m_Setup == nullptr) {
        m_Logger->error("PFFFT does not support FFT size {}", size);
    }
}

// ─────────────────────────────────────
std::unique_ptr<float, void (*)(void *)> RealFft::AllocateBuffer() const {
    float *buffer = static_cast<float *>(pffft_aligned_malloc(static_cast<std::size_t>(m_Size) * sizeof(float)));
    if (buffer == nullptr) {
        m_Logger->critical("PFFFT could not allocate an aligned buffer of {} samples", m_Size);
    }
    return {buffer, pffft_aligned_free};
}

// ─────────────────────────────────────
bool RealFft::IsValid() const {
    return m_Setup != nullptr && m_Work != nullptr;
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
        m_Logger->error("Unsupported audio file format '{}': {}", path.extension().string(), pathString);
        return {};
    }

    AudioFile<float> file;
    file.shouldLogErrorsToConsole(false);
    if (!file.load(pathString)) {
        m_Logger->error("Could not read WAV/AIFF audio file: {}", pathString);
        return {};
    }

    const int channels = file.getNumChannels();
    const int frames = file.getNumSamplesPerChannel();
    const std::uint32_t sampleRate = file.getSampleRate();
    if (frames <= 0 || channels <= 0 || file.samples.size() != static_cast<std::size_t>(channels) || sampleRate == 0 ||
        sampleRate > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
        m_Logger->error("Audio file has invalid dimensions: {}", pathString);
        return {};
    }

    for (const std::vector<float> &channel : file.samples) {
        if (channel.size() < static_cast<std::size_t>(frames)) {
            m_Logger->error("Audio file has inconsistent channel lengths: {}", pathString);
            return {};
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
        m_Logger->error("Training sample rates must be positive");
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
std::optional<std::size_t> Trainer::Cut(const std::vector<float> &samples, int sampleRate,
                                        const std::filesystem::path &label) {
    constexpr int silenceMilliseconds = 100;
    constexpr int preOnsetMilliseconds = 20;
    constexpr int extraMilliseconds = 250;
    constexpr std::size_t fadeSize = 64;

    m_Cancelled.store(false);
    if (label.empty()) {
        m_Logger->error("A non-empty label is required when cutting audio");
        return std::nullopt;
    }
    if (sampleRate <= 0 || m_Config.FFTSize < 512 || m_Config.HopSize <= 0 || m_Config.HopSize > m_Config.FFTSize) {
        m_Logger->error("Invalid sample rate, FFT size, or hop size for cutting audio");
        return std::nullopt;
    }
    if (samples.size() < static_cast<std::size_t>(m_Config.FFTSize)) {
        m_Logger->error("Audio is too short to cut: {} samples", samples.size());
        return std::nullopt;
    }

    const std::filesystem::path directory = ResolvePath(BasePath(), label);
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error) {
        m_Logger->error("Could not create label directory '{}': {}", directory.string(), error.message());
        return std::nullopt;
    }

    ::OpenScofo::OpenScofo extractor(static_cast<float>(sampleRate), static_cast<float>(m_Config.FFTSize),
                                     static_cast<float>(m_Config.HopSize));
    extractor.RequestDescriptor(extractor.GetDescriptorsEnum("db"));
    extractor.RequestDescriptor(extractor.GetDescriptorsEnum("onset"));

    const std::size_t hopSize = static_cast<std::size_t>(m_Config.HopSize);
    const std::size_t silenceRequired =
        static_cast<std::size_t>(silenceMilliseconds) * static_cast<std::size_t>(sampleRate) / 1000U;
    const std::size_t preOnset =
        static_cast<std::size_t>(preOnsetMilliseconds) * static_cast<std::size_t>(sampleRate) / 1000U;
    std::vector<std::pair<std::size_t, std::size_t>> segments;
    bool active = false;
    std::size_t begin = 0;
    std::size_t silenceSamples = 0;

    Report("cutting", 0.0);
    for (std::size_t index = 0; index + hopSize <= samples.size(); index += hopSize) {
        if (!CheckCancelled()) {
            return std::nullopt;
        }
        if (!extractor.ProcessBlock(samples.data() + index, hopSize)) {
            continue;
        }
        const ::OpenScofo::Description description = extractor.GetDescription();
        if (!active) {
            if (description.Onset && description.dB > m_Config.MinDB) {
                begin = index > preOnset ? index - preOnset : 0;
                active = true;
                silenceSamples = 0;
            }
        } else if (description.dB < m_Config.MinDB) {
            silenceSamples += hopSize;
            if (silenceSamples >= silenceRequired) {
                segments.emplace_back(begin, std::min(index + hopSize, samples.size()));
                active = false;
                silenceSamples = 0;
            }
        } else {
            silenceSamples = 0;
        }
        Report("cutting", static_cast<double>(index + hopSize) / static_cast<double>(samples.size()));
    }
    if (active) {
        segments.emplace_back(begin, samples.size());
    }

    const std::size_t extra =
        static_cast<std::size_t>(extraMilliseconds) * static_cast<std::size_t>(sampleRate) / 1000U;
    const std::string prefix = directory.filename().empty() ? "segment" : directory.filename().string();
    std::size_t exported = 0;
    for (std::size_t index = 0; index < segments.size(); ++index) {
        auto [start, end] = segments[index];
        const std::size_t nextStart = index + 1 < segments.size() ? segments[index + 1].first : samples.size();
        end = std::min({end + extra, nextStart, samples.size()});
        if (end <= start) {
            continue;
        }
        if (end - start > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            m_Logger->error("Cut audio segment exceeds the WAV writer's supported size");
            return std::nullopt;
        }

        std::vector<float> chunk(samples.begin() + static_cast<std::ptrdiff_t>(start),
                                 samples.begin() + static_cast<std::ptrdiff_t>(end));
        const std::size_t fade = std::min(fadeSize, chunk.size() / 2);
        for (std::size_t sample = 0; sample < fade; ++sample) {
            const float gain = fade > 1 ? static_cast<float>(sample) / static_cast<float>(fade - 1) : 0.0F;
            chunk[sample] *= gain;
            chunk[chunk.size() - 1 - sample] *= gain;
        }

        std::ostringstream filename;
        filename << prefix << '_' << std::setw(4) << std::setfill('0') << index << ".wav";
        AudioFile<float> output;
        output.setSampleRate(static_cast<std::uint32_t>(sampleRate));
        output.setBitDepth(16);
        output.setNumChannels(1);
        output.setNumSamplesPerChannel(static_cast<int>(chunk.size()));
        output.samples.front() = std::move(chunk);
        const std::filesystem::path outputPath = directory / filename.str();
        if (!output.save(outputPath.string(), AudioFileFormat::Wave)) {
            m_Logger->error("Could not write cut audio file: {}", outputPath.string());
            return std::nullopt;
        }
        ++exported;
    }

    Report("cutting", 1.0);
    m_Logger->info("Exported {} segment(s) to {}", exported, directory.string());
    return exported;
}

// ─────────────────────────────────────
std::vector<float> Trainer::TimeStretch(const std::vector<float> &input, double rate, int fftSize, int hopSize) const {
    if (input.empty()) {
        return {};
    }
    if (!std::isfinite(rate) || rate <= 0.0) {
        m_Logger->error("Time-stretch rate must be positive");
        return {};
    }
    if (fftSize < 32 || hopSize <= 0 || hopSize > fftSize) {
        m_Logger->error("Invalid FFT or hop size for time stretching");
        return {};
    }

    const double outputSizeValue = std::round(static_cast<double>(input.size()) / rate);
    if (outputSizeValue > static_cast<double>(std::numeric_limits<int>::max()) ||
        input.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        m_Logger->error("Time-stretch input or output exceeds Signalsmith Stretch's supported size");
        return {};
    }

    const std::size_t outputSize = std::max<std::size_t>(1, static_cast<std::size_t>(outputSizeValue));
    std::vector<float> output(outputSize);
    const float *inputBuffers[] = {input.data()};
    float *outputBuffers[] = {output.data()};

    signalsmith::stretch::SignalsmithStretch<float> stretch;
    stretch.configure(1, fftSize, hopSize);
    if (!stretch.exact(inputBuffers, static_cast<int>(input.size()), outputBuffers, static_cast<int>(output.size()))) {
        m_Logger->warn("Signalsmith Stretch could not time-stretch the signal; using linear resampling instead");
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
        m_Logger->error("Pitch-shift semitones must be finite");
        return {};
    }
    if (fftSize < 32 || hopSize <= 0 || hopSize > fftSize) {
        m_Logger->error("Invalid FFT or hop size for pitch shifting");
        return {};
    }
    if (input.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        m_Logger->error("Pitch-shift input exceeds Signalsmith Stretch's supported size");
        return {};
    }

    std::vector<float> output(input.size());
    const float *inputBuffers[] = {input.data()};
    float *outputBuffers[] = {output.data()};

    signalsmith::stretch::SignalsmithStretch<float> stretch;
    stretch.configure(1, fftSize, hopSize);
    stretch.setTransposeSemitones(static_cast<float>(semitones));
    if (!stretch.exact(inputBuffers, static_cast<int>(input.size()), outputBuffers, static_cast<int>(output.size()))) {
        m_Logger->warn("Signalsmith Stretch could not pitch-shift the signal; returning it unchanged");
        return input;
    }
    return output;
}

// ─────────────────────────────────────
std::vector<float> Trainer::AddGaussianNoise(const std::vector<float> &input, double deviation,
                                             std::mt19937 &random) const {
    if (deviation < 0.0 || !std::isfinite(deviation)) {
        m_Logger->error("Noise deviation must be finite and non-negative");
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
        m_Logger->error("Convolution of {} samples exceeds PFFFT's supported FFT size", convolutionSize);
        return {};
    }

    RealFft transform(static_cast<int>(fftSize), m_Logger);
    auto time = transform.AllocateBuffer();
    auto firstSpectrum = transform.AllocateBuffer();
    auto secondSpectrum = transform.AllocateBuffer();
    if (!transform.IsValid() || !time || !firstSpectrum || !secondSpectrum) {
        return {};
    }

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
