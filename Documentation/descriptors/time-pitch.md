---
tags:
  - Pitch Detection
---

# Time-Domain and Pitch Descriptors

Use this page for waveform and pitch-related descriptor definitions. These descriptors complement spectral features with noisiness, estimated fundamental frequency, and pitch confidence.

## Reference Table

### Variables

| Symbol | Meaning | Implementation |
| --- | --- | --- |
| $x[n]$ | Input audio sample at index $n$. | `In[n]` |
| $n$ | Time-domain sample index. | Loop variable `i` in `ZeroCrossingRateExec`, and input index in `YINExec` |
| $N$ | Number of samples in the analysis frame. | `m_Config.FFTSize`; `In.size()` in `YINExec` |
| $\mathbb{I}[\cdot]$ | Indicator function: 1 when the condition is true, otherwise 0. | Increment of local variable `crossings` |
| $sign(x[n])$ | Sign of a sample after thresholding near-zero values. | `std::signbit` when `ZCRZeroPos = true`; explicit `-1, 0, 1` sign otherwise |
| $\tau$ | YIN lag, measured in samples. | Loop variable `tau` |
| $\tau_{min}$ | Smallest lag searched by YIN. | `minTau = SR / YINMaxFrequency` |
| $\tau_{max}$ | Largest lag searched by YIN. | `maxTau`, limited by frame size, YIN buffer size, and `YINMinFrequency` |
| $d(\tau)$ | YIN difference function for lag $\tau$. | `Diff[tau]`, stored in `m_YINDifference` |
| $d'(\tau)$ | Cumulative mean normalized difference function (CMNDF). | `Cmnfg[tau]`, stored in `m_YINCMNDF` |
| $j$ | Auxiliary summation index in the CMNDF denominator. | Iteration accumulated into local variable `cumulative` |
| $\hat{\tau}$ | Selected pitch lag before interpolation. | `tauEstimate` |
| $\tau_{refined}$ | Lag after parabolic interpolation. | `refinedTau` |
| $SR$ | Sample rate in Hz. | `m_Config.SR` |
| $Pitch$ | Estimated fundamental frequency in Hz. | Local variable `Pitch`, stored in `Desc.Pitch` |
| $Confidence$ | Pitch confidence in the range $[0,1]$. | Local variable `Confidence`, stored in `Desc.PitchConfidence` |

---

## Zero Crossing Rate

**ID**: `zcr` :custom-librosa:[^5]

Zero Crossing Rate counts how often the waveform crosses the zero-amplitude line, indicating noisiness or percussiveness.

The current implementation optionally pads the frame when `ZCRCenter` is enabled, applies the threshold `ZCRThreshold`, and then counts sign changes. With the default `ZCRZeroPos = true`, zero is treated as non-negative through `std::signbit`.

$$ZCR = \frac{1}{N}\sum_{n=1}^{N-1}\mathbb{I}\left[sign(x[n]) \ne sign(x[n-1])\right]$$

---

## Pitch & PitchConfidence

**ID**: `pitch`

Estimated fundamental frequency and confidence are calculated using the YIN algorithm's cumulative mean normalized difference function (CMNDF).

The current implementation searches for the first CMNDF value below `YINThreshold`, refines the lag by parabolic interpolation, and reports confidence as $1-d'(\tau)$ clamped to $[0,1]$.

$$d'(\tau) = \begin{cases} 1, & \tau = 0 \\ \frac{d(\tau)}{\frac{1}{\tau}\sum_{j=1}^{\tau}d(j)}, & \tau > 0 \end{cases}$$

$$Pitch = \frac{SR}{\tau_{refined}}$$

$$Confidence = \max(0, \min(1, 1-d'(\hat{\tau})))$$

[^5]: `Descriptor` fully compatible with [`librosa`](https://librosa.org/).
