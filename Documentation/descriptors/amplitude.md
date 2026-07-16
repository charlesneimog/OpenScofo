---
tags:
  - Amplitude Descriptors
---

# Amplitude Descriptors

Use this page for the technical definitions of loudness, RMS, dB, and related amplitude features.

The equations below use $x[n]$ for the input audio samples in the current frame. In the implementation, this frame is passed to `MIR::GetDescription` as `In`, and amplitude descriptors are computed first by `MIR::GetSignalPower` in `Sources/OpenScofo/mir.cpp`.

## Reference Table

### Variables

| Symbol | Meaning | Implementation |
| --- | --- | --- |
| $x[n]$ | Input audio sample at index $n$. | `sample` while iterating over `In` |
| $n$ | Time-domain sample index. | Loop iteration over `In` |
| $N$ | Number of samples in the current frame. | `In.size()`; normally `m_Config.FFTSize` |
| $\sum x[n]^2$ | Sum of squared input samples. | Local accumulator `z` |
| $RMS$ | Root mean square amplitude of the current frame. | Local variable `rms`, stored in `Desc.RMS` |
| $L_{dB}$ | Log-amplitude level derived from RMS. | `Desc.dB` |
| $y[n]$ | Sample after the ITU-R BS.1770 filter stages. | Local variable `s2` |
| $\sum y[n]^2$ | Sum of squared filtered samples used for loudness. | Local accumulator `z_loudness` |
| $L$ | Loudness value after ITU-R filtering and logarithmic conversion. | `Desc.Loudness` |
| $L_0$ | Loudness midpoint used by the silence logistic curve. | Local constant `L0 = -60.0` |
| $\alpha$ | Slope of the silence logistic curve. | Local constant `alpha = 0.25` |
| $P_{silence}$ | Silence probability for the current frame. | `Desc.SilenceProb` |
| $X[k]$ | Complex FFT value at bin $k$. | Read from `m_FullFFTOut` in `GetSpectralDescriptions` |
| $X_R[k]$ | Real part of $X[k]$. | Local variable `re` |
| $X_I[k]$ | Imaginary part of $X[k]$. | Local variable `im` |
| $M[k]$ | Raw spectral magnitude, $\sqrt{X_R[k]^2 + X_I[k]^2}$. | Local variable `mag`, stored in `Desc.Magnitude[k]` |
| $A[k]$ | FFT-size-normalized magnitude, $M[k] / N$. | Local variable `norm`, stored in `Desc.SpectralMagnitudeNorm[k]` |
| $MaxAmp$ | Maximum normalized spectral magnitude in the frame. | `Desc.MaxAmp` |

---

## Root Mean Square (RMS)

**ID**: `rms` :custom-librosa:[^1] :custom-essentia:

RMS measures frame energy. Higher values indicate louder sound; lower values indicate quiet sound or silence.

`OpenScofo` implements the following equation:

$$RMS = \sqrt{\frac{1}{N} \sum_{n=0}^{N-1} x[n]^2}$$


This is the root mean square of the time-domain signal $x[n]$ of length $N$.

---

## dB

**ID**: `db`

dB is a logarithmic sound-level measure derived from RMS.

The equation implemented is:

$$L_{dB} = 20 \log_{10}(RMS)$$

If $RMS = 0$, `OpenScofo` returns $-100$ instead of $-\infty$.

`librosa` and `essentia` do not expose this exact `dB` descriptor, but converting compatible `RMS` values gives the same scale.

---

## Max Amplitude

**ID**: `maxamp`

Maximum normalized spectral amplitude in the current frame:

$$MaxAmp = \max_k A[k]$$

---

## Loudness

**ID**: `loudness`

Loudness estimates perceived sound intensity after perceptual filtering.

$$L = -0.691 + 10 \log_{10}\left(\frac{1}{N}\sum_{n=0}^{N-1} y[n]^2\right)$$

Here \(y[n]\) is the audio after the filtering stage defined by [`ITU‑R BS.1770`](https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.1770-5-202311-I!!PDF-E.pdf){:target="_blank"}.

The term \(N\) is the number of samples in the analyzed frame.

`OpenScofo` follows the ITU-R BS.1770 method rather than Essentia's simplified loudness descriptor. The implementation is based on [klangfreund/LUFSMeter](https://github.com/klangfreund/LUFSMeter).


---

## Silence Probability

**ID**: `silence`

Probability that the current frame is silence, derived from loudness ($L$) with $\alpha = 0.25$ and $L_0 = -60.0$:

$$P_{silence} = \frac{1}{1 + e^{\alpha (L - L_0)}}$$

---
[^1]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-9}$.
