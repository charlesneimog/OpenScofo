# Amplitude Descriptors

The equations below use $x[n]$ for the input audio samples in the current frame. In the implementation, this frame is passed to `MIR::GetDescription` as `In`, and amplitude descriptors are computed first by `MIR::GetSignalPower` in `Sources/OpenScofo/mir.cpp`.

## Variable Reference

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

A measure of the amplitude (energy) of the current audio frame. It represents how loud the sound is within a short time window. Higher RMS values indicate louder sounds, while lower values indicate quieter sounds or silence.

`OpenScofo` implements the following equation:

$$RMS = \sqrt{\frac{1}{N} \sum_{n=0}^{N-1} x[n]^2}$$


Where root mean square of the amplitude representing the energy of the time-domain signal $x[n]$ of length $N$.

---

## dB

**ID**: `db` 

A logarithmic measure of sound level derived from the signal amplitude. Unlike RMS, which measures the raw energy of the signal, dB expresses this energy on a logarithmic scale that better reflects how humans perceive changes in loudness.

The equation implemented is:

$$L_{dB} = 20 \log_{10}(RMS)$$

If $RMS = 0$, `OpenScofo` returns $-100$ instead of $-\infty$.
        
Note that neither `librosa` or `essentia` implemented `dB`, but as it uses `RMS`, once you convert `RMS` to `dB` it is compatible with both.

---

## Max Amplitude

**ID**: `maxamp` 

Maximum normalized spectral amplitude detected in the current frame of audio. It is computed during the spectral pass from the FFT-size-normalized magnitude spectrum:

$$MaxAmp = \max_k A[k]$$

---

## Loudness

**ID**: `loudness` 

An estimate of perceived sound intensity based on psychoacoustic models of human hearing. Unlike dB, it applies perceptual models and filters derived from psychoacoustic studies to approximate how humans actually perceive loudness.

Equation used is:

$$L = -0.691 + 10 \log_{10}\left(\frac{1}{N}\sum_{n=0}^{N-1} y[n]^2\right)$$

where \(y[n]\) denotes the audio samples after applying the filtering stage defined in the [`ITU‑R BS.1770`](https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.1770-5-202311-I!!PDF-E.pdf){:target="_blank"} recommendation.  

The term \(N\) represents the number of samples in the analyzed frame. This formulation corresponds to the energy-based loudness estimate used in the loudness measurement procedure defined by the standard.

The `loudness` descriptor implemented in `essentia` is based on a simplified perceptual model derived from signal energy with a power-law compression. While computationally inexpensive, it does not incorporate perceptual frequency weighting or the measurement procedure defined in modern broadcast loudness standards.

`OpenScofo` instead implements loudness estimation following the methodology described in [`ITU‑R BS.1770`](https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.1770-5-202311-I!!PDF-E.pdf){:target="_blank"}. This approach applies perceptual filtering prior to the energy calculation and expresses the result in a logarithmic scale, which aligns with the methodology adopted in contemporary loudness measurement practices for audio production and broadcasting.

As reference, `OpenScofo` implements the code implemented in [klangfreund/LUFSMeter](https://github.com/klangfreund/LUFSMeter).


---

## Silence Probability

**ID**: `silence` 

Probability that the current frame corresponds to silence, derived from Loudness ($L$) via a logistic function where $\alpha = 0.25$ and $L_0 = -60.0$:

$$P_{silence} = \frac{1}{1 + e^{\alpha (L - L_0)}}$$

---
[^1]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-9}$.
