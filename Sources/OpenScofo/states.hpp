/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#pragma once

#include <string>
#include <vector>
#include <variant>
#include <span>
#include <unordered_map>
#include <filesystem>

#include <onsetsds.h>

namespace OpenScofo {

namespace fs = std::filesystem;

enum AudioDescType {
    PITCH,
    SILENCE,
    LABEL,
    ONSET,
};

// ─────────────────────────────────────
enum Descriptors {
    INVALID = -1,

    ODSONSET,

    // Amplitude
    LOUDNESS,
    DB,
    MAXAMP,
    RMS,
    STDDEV,
    MAGNITUDE,
    POWERARRAY,
    SILENCEPROB,

    // Spectral Arrays
    MFCC,
    CHROMA,
    LOGMEL,

    // Spectral
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

    // Pitch
    YIN,
    YINCONFIDENCE,

    // Percussive
    EXTENDEDTECHNIQUE,

    // AI
    ONNX,
};

// ─────────────────────────────────────
enum Mode { SCOREFOLLOWER, DESCRIPTORS };

// ─────────────────────────────────────
enum EventType {
    FIRSTEVENT,

    // Tradicional
    REST,
    NOTE,
    CHORD,
    TRILL,
    MULTI,
    PTECH,
    UTECH,

    // Cort Lippe event (TODO:)
    EVENT,
};

enum HMMType { SEMIMARKOV, MARKOV };

// ─────────────────────────────────────
struct ScoreAction {
    bool isLua;
    bool isAudioStateChange = false;
    std::string Lua;
    std::string Receiver;
    std::vector<std::variant<float, int, std::string>> Args;
    bool AbsoluteTime;
    double Time;
};

using EventActions = std::vector<ScoreAction>;

// ─────────────────────────────────────
struct AudioState {
    AudioDescType Type;
    double Freq = 0;
    double Midi = 0;
    std::string Label = "";
    unsigned Index = 0;
};

// ─────────────────────────────────────
struct MarkovState {
    int Index;
    int ScorePos;
    int MarkovIndex = -1;
    std::vector<AudioState> AudioStates;
    int BestAudioStateIndex = -1;

    // States Actions
    HMMType HSMMType;
    EventType Type;
    EventActions Actions;

    // Inference
    double InitProb;
    std::vector<double> Forward;
    std::vector<double> ExitProb;
    std::vector<double> BestObs;

    // Time
    int UpperBound = 0;
    double BPMExpected = 0;
    double BPMObserved = 0;
    double OnsetExpected = 0.0;
    double OnsetObserved = 0;
    double PhaseExpected;
    double PhaseObserved = 0;
    double IOIPhiN = 0;
    double IOIHatPhiN = 0;
    double Duration = 0.0;

    // Configuration for each event
    double PhaseCoupling = 0.5;
    double SyncStrength = 0.5;

    // Error Handling
    int Line;
};

using States = std::vector<MarkovState>;

// ─────────────────────────────────────
struct Description {
    double Onset;

    // Probability
    double SilenceProb;
    double ExtendedTechProb;
    int WindowLastOnset = 0;

    // Amplitude
    double dB;
    double RMS;
    double MaxAmp;
    double Loudness;

    // Spectral
    double Harmonicity;
    double SpectralFlatness;
    double SpectralEntropy = 0.0;
    double SpectralRolloff = 0.0;
    double SpectralFlux;
    double SpectralIrregularity = 0.0;
    double SpectralIrregularityJensen = 0.0;
    double SpectralIrregularityKrimphoff = 0.0;
    double SpectralCrest = 0.0;
    double SpectralCentroid = 0.0;
    double CentroidVelocity = 0.0;
    double SpectralSpreadHz = 0.0;
    double SpectralSpreadVariance = 0.0;
    double SpectralSkewness = 0.0;
    double SpectralSlope = 0.0;
    double SpectralKurtosis = 0.0;
    double HighFreqRatio = 0.0;
    double ZeroCrossingRate;
    double StdDev;
    double Pitch = 0.0;
    double PitchConfidence = 0.0;

    std::vector<double> Power;
    std::vector<double> Magnitude;
    std::vector<double> SpectralMagnitudeNorm;
    std::vector<double> SpectralMagnitudeFrameNorm;
    std::vector<double> SpectralPowerFrameNorm;

    std::vector<double> LogMelSpectrum;
    std::vector<double> ReverbSpectralPower;

    std::vector<double> MFCC;
    std::vector<double> Chroma;

    // ONNX
    std::unordered_map<std::string, float> ONNX;
};

// ─────────────────────────────────────
struct Configuration {
    // Audio Parameters
    float SR = 48000;
    float FFTSize = 2048;
    float HOPSize = 512;

    double TuningA4 = 440.0;

    // Pitch Template
    float PitchTemplateSigma = 0.5;
    int PitchTemplateHarmonics = 10;

    // MFCC
    int MFCCMels = 40;
    int MFCCCount = 13;

    // Onset
    onsetsds_odf_types OnsetType = ODS_ODF_MKL;
    int MedSpan = 50;

    // Silence Threshold
    float dBTreshold = -60;

    // Yin
    double SpectralRolloffCutoff = 0.85;
    double YINThreshold = 0.15;
    double YINMinFrequency = 50.0;
    double YINMaxFrequency = 2000.0;

    // Chroma
    int ChromaSize = 12;
    double ChromaCenterOctave = 5.0;
    double ChromaOctaveWidth = 2.0;

    // ZCR
    bool ZCRCenter = true;
    bool ZCRPad = false;
    bool ZCRZeroPos = true;
    double ZCRThreshold = 1e-10;

    // Temporal model
    float SyncStrength = 0.5;
    float PhaseCoupling = 0.5;

    // Audio state listener
    std::string AudioStateChangeReceiver;

    // ONNX
    fs::path TimbreONNXModel;
    std::vector<std::string> ONNXDescriptors;
    std::vector<Descriptors> RequestedDescriptors;
};

// ─────────────────────────────────────
struct SpectralAccumulators {
    double SumPower = 0.0;
    double WSumFreqs = 0.0;
    double irregularityJensenNumerator = 0.0;
    double irregularityKrimphoffSum = 0.0;
    double irregularityDenominator = 0.0;
    double logSumPower = 0.0;
    double linSumPower = 0.0;
    double spectralEnergySum = 0.0;
    double harmonicityPeak = 0.0;
    double harmonicitySum = 0.0;
    double highFreqEnergy = 0.0;
    double sumFreq2 = 0.0, sumFreq3 = 0.0, sumFreq4 = 0.0;
    double sumIndex = 0.0, sumIndex2 = 0.0;
    double sumFreq = 0.0, sumFreqSq = 0.0;
    double maxMag = 0.0;
    float sumMagCrest = 0.0f;
    float maxMagCrest = 0.0f;
};

} // namespace OpenScofo
