---
icon: lucide/waves
tags:
  - Pitch Detection
---

# Time-Domain and Pitch

Use this page for waveform and pitch-related descriptor definitions. These descriptors complement spectral features with noisiness, estimated fundamental frequency, and pitch confidence.

## Reference Table

### Variables

| Symbol | Meaning | Implementation |
| --- | --- | --- |
| $x[n]$ | Input audio sample at index $n$. | `In[n]` |
| $n$ | Time-domain sample index. | Loop variable `i` in `ZeroCrossingRateExec` |
| $N$ | Number of samples in the analysis frame. | `m_Config.FFTSize` |
| $\mathbb{I}[\cdot]$ | Indicator function: 1 when the condition is true, otherwise 0. | Increment of local variable `crossings` |
| $sign(x[n])$ | Sign of a sample after thresholding near-zero values. | `std::signbit` when `ZCRZeroPos = true`; explicit `-1, 0, 1` sign otherwise |
| $SR$ | Sample rate in Hz. | `m_Config.SR` |
| $k$ | FFT-bin index of a local spectral peak. | Candidate index in `SigmundExec` |
| $\delta$ | Sub-bin offset of a spectral peak. | Hann peak interpolation in `sigmund.cpp` |
| $A_p$ | Interpolated amplitude of peak $p$. | `SigmundPeak::Amplitude` |
| $Pitch$ | Estimated fundamental frequency in Hz. | `SigmundEstimate::Frequency`, stored in `Desc.Pitch` |
| $Confidence$ | Fraction of salient-peak loudness supporting the winning pitch, in $[0,1]$. | `SigmundEstimate::Quality`, stored in `Desc.PitchConfidence` |

---

## Zero Crossing Rate

**ID**: `zcr` :custom-librosa:[^5]

Zero Crossing Rate counts how often the waveform crosses the zero-amplitude line, indicating noisiness or percussiveness.

The current implementation optionally pads the frame when `ZCRCenter` is enabled, applies the threshold `ZCRThreshold`, and then counts sign changes. With the default `ZCRZeroPos = true`, zero is treated as non-negative through `std::signbit`.

$$ZCR = \frac{1}{N}\sum_{n=1}^{N-1}\mathbb{I}\left[sign(x[n]) \ne sign(x[n-1])\right]$$

---

## Pitch & PitchConfidence

**IDs**: `yin`, `yin_confidence` (legacy public names)

Estimated fundamental frequency and confidence are calculated with Sigmund-style sinusoidal-peak and harmonic voting. The detector reuses the same $N$-point Hann-windowed PFFFT spectrum used by the spectral descriptors; it performs no second FFT and no zero-padding.

For neighboring magnitudes $L$, $C$, and $R$ around a local maximum, the Hann peak offset and frequency are:

$$\delta = \operatorname{clamp}\left(\frac{2(R-L)}{L+2C+R},-\frac{1}{2},\frac{1}{2}\right)$$

$$f_p = (k+\delta)\frac{SR}{N}$$

Up to 20 masked peaks are found, the 12 most salient vote across 16 possible harmonic numbers, and the strongest accepted pitch is refined from its supporting harmonics. Confidence is the supporting loudness divided by total salient-peak loudness. `YINMinFrequency` and `YINMaxFrequency` remain the final output limits for compatibility.

[^5]: `Descriptor` fully compatible with [`librosa`](https://librosa.org/).
