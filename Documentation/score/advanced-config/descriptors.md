---
tags:
  - Advanced Configuration
  - Audio Analysis
---

# Descriptors

These settings shape the audio measurements used for onset detection, pitch estimation, and AI model input. Keep them consistent with model training. Setting a descriptor parameter does not itself request that descriptor; the host API, score events, and model inputs determine which optional analyses are needed.

## Audio analysis

| Keyword | Default | Accepted value | Effect |
| --- | --- | --- | --- |
| `SR` | `48000` Hz in parsed configuration | Positive integer | Analysis sample rate. Match the host audio rate and model training rate. |
| `FFTSIZE` | `2048` samples | Positive power of two | Analysis window size. Larger windows cover more audio and provide finer frequency bins. |
| `HOPSIZE` | `512` samples | Positive power of two | Spacing between analysis updates. Smaller hops request more frequent updates. |

`SR` does not resample the audio. Score loading reports a mismatch with the engine's existing rate. Set it explicitly when using a rate other than 48000 Hz.

Integer settings are converted from numeric values before validation. Use whole numbers. The parser validates each option separately; it does not check every relationship between settings.

## Onset and articulation

| Keyword | Default | Accepted value | Effect |
| --- | --- | --- | --- |
| `ONSETFUNCTION` | `mkl` | One of the lowercase names below | Selects the onset detection function. |
| `MEDSPAN` | `50` | Positive integer | Median span used by the onset detector. |

| Value | Implementation selected by the score parser |
| --- | --- |
| `pow` | Power (`ODS_ODF_POWER`) |
| `pd` | Phase deviation (`ODS_ODF_PHASE`) |
| `wpd` | Weighted phase deviation (`ODS_ODF_WPHASE`) |
| `sf`, `hfc` | Both select magnitude sum (`ODS_ODF_MAGSUM`) |
| `cd` | Complex domain (`ODS_ODF_COMPLEX`) |
| `rcd` | Rectified complex domain (`ODS_ODF_RCOMPLEX`) |
| `mkl` | Modified Kullback–Leibler (`ODS_ODF_MKL`) |

In the current parser, `sf` and `hfc` are equivalent selections. Switching between them does not select different algorithms.

Test onset changes with the exact articulation that fails. Compare attack detection before changing `FFTSIZE` or `HOPSIZE`, which also affect other analyses. See [Onset Descriptor](../../descriptors/onset.md) for the output.

## MFCC and logmel

| Keyword | Default | Accepted value | Effect |
| --- | --- | --- | --- |
| `MFCCMELS` | `40` | Positive integer | Number of mel bands used by MFCC and logmel. |
| `MFCCCOUNT` | `13` | Positive integer | Number of MFCC coefficients. |

These settings change model input dimensions. Match the trained feature extraction setup; changing them is not a general remedy for poor technique recognition.

## YIN pitch estimation

| Keyword | Default | Accepted value | Effect |
| --- | --- | --- | --- |
| `YINTHRESHOLD` | `0.15` | `0` to `1`, inclusive | Threshold for the first candidate minimum in YIN's normalized difference function. |
| `YINMINFREQUENCY` | `50` Hz | Positive number | Lower frequency limit for estimated pitch. |
| `YINMAXFREQUENCY` | `2000` Hz | Positive number | Upper frequency limit for estimated pitch. |

Choose a minimum below the maximum and a range supported by the sample rate and window size. The parser checks positivity separately, without enforcing the relationship between the limits. YIN supplies pitch descriptors; harmonic template matching has its own [Pitch Following](pitch-following.md) settings.

## Spectral rolloff and chroma

| Keyword | Default | Accepted value | Effect |
| --- | --- | --- | --- |
| `SPECTRALROLLOFFCUTOFF` | `0.85` | `0` to `1`, inclusive | Fraction of spectral power used to locate the rolloff frequency. |
| `CHROMASIZE` | `12` | Positive integer | Number of chroma bins. |
| `CHROMACENTEROCTAVE` | `5` | Any number | Center octave of the chroma weighting. |
| `CHROMAOCTAVEWIDTH` | `2` | Positive number | Width of the chroma octave weighting. |

## Zero-crossing rate

| Keyword | Default | Accepted value | Effect |
| --- | --- | --- | --- |
| `ZCRCENTER` | `true` | Boolean | Enables edge-value padding for centered ZCR analysis. |
| `ZCRPAD` | `false` | Boolean | Counts an initial crossing when enabled. |
| `ZCRZEROPOS` | `true` | Boolean | Uses sign-bit comparisons; when disabled, zero is treated as a separate sign value. |
| `ZCRTHRESHOLD` | `0.0000000001` | Non-negative number | Samples with absolute magnitude at or below this value are treated as zero. |

For these three boolean options, the parser accepts `true`, `on`, `yes`, or `1`, and `false`, `off`, `no`, or `0`. Word values are case-insensitive. `ZCRPAD` controls the initial crossing count; `ZCRCENTER` controls edge padding.

## Example

```openscofo
BPM 72
SR 48000
FFTSIZE 2048
HOPSIZE 512
ONSETFUNCTION mkl
MEDSPAN 50
MFCCMELS 40
MFCCCOUNT 13
ZCRCENTER ON
ZCRPAD OFF
ZCRZEROPOS ON
ZCRTHRESHOLD 0.0000000001

NOTE C4 1
```

This illustrates configuration syntax using defaults. For descriptor meanings and output values, see [Audio Analysis](../../descriptors/index.md). For input selection and metadata precedence, see [AI Model Recognition](ai-model-recognition.md).
