
#include <mir/mir.hpp>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <vector>

using namespace OpenScofo;

int main() {
    MIR mir;
    Configuration config;
    Description desc;
    config.SR = 48000;
    config.FFTSize = 2048;
    config.RequestedDescriptors = {ODSONSET,
                                   LOUDNESS,
                                   DB,
                                   MAXAMP,
                                   RMS,
                                   STDDEV,
                                   MAGNITUDE,
                                   POWERARRAY,
                                   SILENCEPROB,
                                   MFCC,
                                   CHROMA,
                                   LOGMEL,
                                   ZCR,
                                   HFR,
                                   CENTROID,
                                   SPREADHZ,
                                   SPREADVARIANCE,
                                   CREST,
                                   FLATNESS,
                                   ENTROPY,
                                   ROLLOFF,
                                   CENTROIDVEL,
                                   FLUX,
                                   SKEWNESS,
                                   SLOPE,
                                   KURTOSIS,
                                   IRREGULARITY,
                                   HARMONICITY,
                                   YIN,
                                   YINCONFIDENCE,
                                   EXTENDEDTECHNIQUE,
                                   ONNX};

    mir.UpdateConfiguration(config);

    const size_t nBins = config.FFTSize / 2 + 1;

    desc.Power.resize(nBins);
    desc.Magnitude.resize(nBins);
    desc.SpectralMagnitudeNorm.resize(nBins);
    desc.SpectralMagnitudeFrameNorm.resize(nBins);
    desc.SpectralPowerFrameNorm.resize(nBins);
    desc.ReverbSpectralPower.resize(nBins);
    desc.MFCC.resize(config.MFCCCount);
    desc.LogMelSpectrum.resize(config.MFCCMels);
    desc.Chroma.resize(config.ChromaSize);
    desc.Power.resize(nBins);
    desc.Magnitude.resize(nBins);
    desc.SpectralMagnitudeNorm.resize(nBins);
    desc.SpectralMagnitudeFrameNorm.resize(nBins);
    desc.SpectralPowerFrameNorm.resize(nBins);
    desc.ReverbSpectralPower.resize(nBins);

    desc.LogMelSpectrum.resize(config.MFCCMels);
    desc.MFCC.resize(config.MFCCCount);
    desc.Chroma.resize(config.ChromaSize);

    // Generate 440 Hz sine wave with random noise
    std::vector<double> in(config.FFTSize);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> noise(-0.01, 0.01);

    const double frequency = 440.0;
    const double twoPi = 2.0 * std::numbers::pi;

    for (size_t i = 0; i < in.size(); ++i) {
        const double time = static_cast<double>(i) / config.SR;
        in[i] = std::sin(twoPi * frequency * time) + noise(rng);
    }

    for (size_t i = 0; i < 100000; ++i) {
        mir.GetDescription(in, desc);
        // if (fmod(i, 1000) == 0) {
        //     std::cout << "Done " << i << " \n";
        // }
    }

    std::cout << "Expected pitch: " << frequency << " Hz\n";
    std::cout << "Detected pitch: " << desc.Pitch << " Hz\n";
    std::cout << "Confidence: " << desc.PitchConfidence << '\n';
    return 0;
}
