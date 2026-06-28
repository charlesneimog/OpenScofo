# Spectral Descriptors

Spectral descriptors describe aspects of tone color: brightness, noisiness, concentration of energy, changes between moments, and the presence of pitch-class information. In musical terms, they help `OpenScofo` distinguish sounds such as sustained tones, noisy attacks, percussive gestures, bright spectra, dark spectra, and timbral changes over time.

The equations below show how these descriptors are implemented in `OpenScofo`. They are included as a technical reference, but you do not need to read every formula to use the descriptors musically. Most of them are computed from the spectrum of a short audio frame after a Hann window and FFT.

For the equations, $X[k] = X_R[k] + iX_I[k]$ is one frequency bin of the FFT, with real part $X_R[k]$ and imaginary part $X_I[k]$. Its raw magnitude is $M[k] = \sqrt{X_R[k]^2 + X_I[k]^2}$, its normalized magnitude is $A[k] = M[k] / N$, its power is $P[k] = X_R[k]^2 + X_I[k]^2$, and $f_k$ is the frequency represented by that bin. $N$ is the FFT size and $K$ is the number of one-sided FFT bins.

## Variable Reference

The table below maps the notation used in this page to the implementation in `Sources/OpenScofo/mir.cpp`.

| Symbol | Meaning | Implementation |
| --- | --- | --- |
| $N$ | FFT size. | `m_Config.FFTSize` |
| $K$ | Number of one-sided FFT bins, including DC and Nyquist. | `NHalf = FFTSize / 2 + 1` |
| $k$ | FFT-bin index. | Loop variable `i` or `bin` |
| $i$ | Generic index; in MFCC equations, the MFCC coefficient index. | Loop variable `i` or `k`, depending on the function |
| $j$ | Auxiliary summation index. | Local loop variable when needed |
| $m$ | Mel-band index. | Loop variable `m` |
| $X[k]$ | Complex FFT value at bin $k$. | Read from `m_FullFFTOut` after `pffft_transform_ordered` |
| $X_R[k]$ | Real part of $X[k]$. | Local variable `re` |
| $X_I[k]$ | Imaginary part of $X[k]$. | Local variable `im` |
| $P[k]$ | Raw power spectrum, $X_R[k]^2 + X_I[k]^2$. | Local variable `p`, stored in `Desc.Power[k]` |
| $M[k]$ | Raw magnitude spectrum, $\sqrt{P[k]}$. | Local variable `mag`, stored in `Desc.Magnitude[k]` |
| $A[k]$ | FFT-size-normalized magnitude, $M[k] / N$. | Local variable `norm`, stored in `Desc.SpectralMagnitudeNorm[k]` |
| $A_{norm}[k]$ | Frame-normalized magnitude, $A[k] / \sum_j A[j]$. | Local variable `normSp`, stored in `Desc.SpectralMagnitudeFrameNorm[k]` |
| $f_k$ | Frequency in Hz represented by bin $k$. | `i * binWidth`, where `binWidth = SR / FFTSize` |
| $\epsilon$ | Small positive value used to avoid division by zero or logarithms of zero. | `1e-10` for log/power floors, `1e-12` for denominator guards |
| $t$ | Analysis-frame index. | Current call to `GetSpectralDescriptions`; previous-frame values are stored in `m_PreviousSpectralPower` and `m_PrevCentroid` |
| $\mu$ | Frequency centroid in Hz. | `Desc.SpectralCentroid` |
| $\sigma$ | Spectral spread in Hz. | `Desc.SpectralSpreadHz` |
| $\mu_k$ | Centroid over FFT-bin indices, used for `spread_variance`. | Local variable `EIndex` |
| $r$ | Spectral rolloff cutoff ratio. | `m_Config.SpectralRolloffCutoff` |
| $b$ | First bin whose cumulative power reaches the rolloff cutoff. | Local variable `rolloffBin` |
| $p_k$ | Power-normalized probability for entropy. | Local variable `prob` |
| $H_m[k]$ | Mel-filter weight for mel band $m$ and FFT bin $k$. | `m_MFCCFilter[m][k]` |
| $E_m$ | Mel-band energy before logarithmic compression. | Local variable `MelEnergy`, stored in `m_MFCCEnergy[m]` |
| $L_m$ | Log-mel value before the 80 dB floor is applied. | Local variable `LogMel` |
| $LogMel[m]$ | Log-mel spectrum after applying the 80 dB floor. | `Desc.LogMelSpectrum[m]` |
| $B$ | Number of mel bands. | `m_Config.MFCCMels` |
| $D_{i,m}$ | DCT-II basis coefficient for MFCC coefficient $i$ and mel band $m$. | `m_DCTBasis[i][m]` |
| $\alpha_i$ | DCT normalization factor. | Computed when `m_DCTBasis` is initialized |
| $W_c[k]$ | Chroma-filter weight for pitch class $c$ and FFT bin $k$. | `m_ChromaFilter[c][k]` |
| $c$ | Chroma pitch-class index. | Loop variable `chroma` |

---

## Spectral Flatness

**ID**: `flatness` :custom-librosa:[^2]

Spectral flatness indicates how noisy versus tonal a sound is. A high flatness means the spectrum is uniform like white noise, while a low flatness shows clear peaks, like a sustained musical note.

The current implementation is based on the equation below, where $K$ is the number of one-sided FFT bins and $\epsilon = 10^{-10}$ prevents logarithms of zero.

$$Flatness = \frac{\exp\left(\frac{1}{K}\sum_{k=0}^{K-1}\ln(\max(\epsilon, P[k]))\right)}{\frac{1}{K}\sum_{k=0}^{K-1}\max(\epsilon, P[k])}$$

---

## Spectral Flux

**ID**: `flux` :custom-essentia:[^6]

Spectral flux measures how quickly the spectrum of a sound changes over time. High flux indicates sudden changes or transients, like drum hits, while low flux corresponds to steady, continuous sounds.

The current implementation compares the current magnitude spectrum with the magnitude spectrum stored from the previous analysis frame.

$$Flux_t = \sqrt{\sum_{k=0}^{K-1}\left(M_t[k] - M_{t-1}[k]\right)^2}$$

---

## Spectral Irregularity

**ID**: `irregularity`

Spectral irregularity quantifies how uneven or jagged a spectrum is between adjacent frequency bins. High irregularity indicates complex, inharmonic, or noisy timbres, while low values suggest smooth, harmonic sounds.

`OpenScofo` computes both Jensen's bin-based algorithm and Krimphoff's running-average strategy. The `irregularity` descriptor uses the Krimphoff value by default.

Jensen's algorithm, cited by Brent, divides squared adjacent-bin magnitude differences by the total spectral power.[^7] With $K = N/2 + 1$ one-sided FFT bins, this is:

$$Irregularity = \frac{\sum_{k=0}^{K-2}\left(M[k] - M[k+1]\right)^2}{\sum_{k=0}^{K-1}M[k]^2}$$

Krimphoff's strategy compares each interior bin against the average of itself and its two neighbors:

$$Irregularity = \log_{10}\left(\sum_{k=1}^{K-2}\left|M[k] - \frac{M[k-1] + M[k] + M[k+1]}{3}\right|\right)$$

---

## Spectral Crest

**ID**: `crest` :custom-essentia:[^6]

Spectral crest measures the ratio of the highest spectral peak to the average spectral amplitude. A high crest indicates a tone dominated by strong harmonics or transients, while a low crest corresponds to more even, noise-like spectra.

The current implementation uses the raw magnitude spectrum.

$$Crest = \frac{\max_k M[k]}{\frac{1}{K}\sum_{k=0}^{K-1}M[k]}$$

---

## Spectral Skewness

**ID**: `skewness` :custom-essentia:[^6]

Spectral skewness measures the asymmetry of the spectral distribution around its centroid. It indicates whether the energy is biased toward low or high frequencies.

The current implementation uses normalized magnitude as the spectral weight. Let $\mu$ be the spectral centroid and $\sigma$ be the spectral spread in Hz.

$$\mu = \frac{\sum_{k=0}^{K-1} f_k A[k]}{\sum_{k=0}^{K-1}A[k]}$$

$$\sigma = \sqrt{\frac{\sum_{k=0}^{K-1} f_k^2 A[k]}{\sum_{k=0}^{K-1}A[k]} - \mu^2}$$

$$Skewness = \frac{\frac{\sum_{k=0}^{K-1}(f_k-\mu)^3A[k]}{\sum_{k=0}^{K-1}A[k]}}{(\sigma + \epsilon)^3}$$

---

## Spectral Kurtosis

**ID**: `kurtosis` :custom-essentia:[^6]

Spectral kurtosis measures the peakedness or tailedness of the spectral distribution around its centroid. It quantifies how concentrated the spectral energy is in a few frequency bins versus being evenly spread.

The current implementation uses the fourth central moment and returns excess kurtosis by subtracting 3.

$$Kurtosis = \frac{\frac{\sum_{k=0}^{K-1}(f_k-\mu)^4A[k]}{\sum_{k=0}^{K-1}A[k]}}{(\sigma + \epsilon)^4} - 3$$

---

## Spectral RollOff

**ID**: `rolloff` :custom-essentia:[^6]

Spectral rolloff indicates the frequency below which a fixed percentage of a sound's spectral energy is contained. Higher values make the sound perceptually brighter or sharper, while lower values make it darker or warmer.

The current implementation uses a cutoff ratio $r$ from `SpectralRolloffCutoff` (default $0.85$). It finds the first bin $b$ whose cumulative power reaches $r$ times the stabilized total power.

$$b = \min\left\{i : \sum_{k=0}^{i}P[k] \ge r\sum_{k=0}^{K-1}\max(\epsilon, P[k])\right\}$$

$$Rolloff = f_b$$

---

## Spectral Entropy

**ID**: `entropy` :custom-essentia:[^6]

Spectral entropy indicates how uniformly a sound's spectral energy is distributed across frequencies. Higher values make the sound perceptually more noisy or disordered, while lower values make it more tonal or structured.

The current implementation normalizes the power spectrum into a probability distribution and computes Shannon entropy with base-2 logarithms.

$$p_k = \frac{P[k]}{\sum_{j=0}^{K-1}P[j]}$$

$$Entropy = -\sum_{k=0}^{K-1}p_k\log_2(p_k)$$

---

## Spectral Centroid

**ID**: `centroid` :custom-librosa:[^3]

Spectral centroid indicates the center of mass of a sound's spectrum. Higher values make the sound perceptually brighter, while lower values make it darker or warmer.

The current implementation uses normalized magnitudes as weights and returns the result in Hz.

$$Centroid = \frac{\sum_{k=0}^{K-1} f_k A[k]}{\sum_{k=0}^{K-1}A[k]}$$

---

## Centroid Velocity

**ID**: `velocity` :custom-essentia:[^6]

Centroid velocity measures how quickly the spectral centroid changes over time, reflecting dynamic shifts in brightness or timbre.

The current implementation stores the previous centroid and returns the absolute difference.

$$Velocity_t = |Centroid_t - Centroid_{t-1}|$$

---

## Spectral Spread

**ID**: `spread` :custom-librosa:[^3]

Spectral spread quantifies how dispersed the energy is around the spectral centroid, indicating whether the sound is focused (narrow) or diffuse (wide) in frequency.

The current implementation computes the standard deviation of frequency around the centroid, using normalized magnitudes as weights.

$$Spread = \sqrt{\frac{\sum_{k=0}^{K-1} f_k^2 A[k]}{\sum_{k=0}^{K-1}A[k]} - Centroid^2}$$

The `spread_variance` descriptor follows Essentia's `CentralMoments` + `DistributionShape` spread value. It treats the normalized magnitude spectrum as a probability distribution over FFT-bin indices and returns the normalized second central moment:

$$SpreadVariance = \frac{\sum_{k=0}^{K-1}(k-\mu_k)^2 A[k]}{(K-1)^2\sum_{k=0}^{K-1}A[k]}$$

---

## High Frequency Ratio

**ID**: `hfr`

High Frequency Ratio measures the proportion of energy in the upper part of the spectrum, reflecting the brightness or presence of high-pitched content in a sound.

The current implementation starts the high-frequency region at $K/4$ and uses normalized magnitudes.

$$HFR = \frac{\sum_{k=\lfloor K/4 \rfloor}^{K-1}A[k]}{\sum_{k=0}^{K-1}A[k]}$$

---

## Standard Deviation

**ID**: `std`

Standard deviation describes how far the frame-normalized spectral magnitude is from a uniform distribution.

The current implementation uses $A_{norm}[k] = A[k] / \sum_j A[j]$ and compares it to the uniform mean $\mu = 1/K$.

$$StdDev = \sqrt{\frac{1}{K}\sum_{k=0}^{K-1}\left(A_{norm}[k]-\frac{1}{K}\right)^2}$$

---

## Harmonicity

**ID**: `harmonicity`

Harmonicity measures how concentrated the normalized magnitude spectrum is around its strongest non-DC bin. High harmonicity indicates a clear dominant component, while low harmonicity indicates a more distributed or noise-like spectrum.

The current implementation is intentionally simple: it ignores bin 0 and returns the ratio between the largest normalized magnitude and the sum of normalized magnitudes over bins $1$ to $K-1$.

$$Harmonicity = \frac{\max_{1 \le k < K} A[k]}{\sum_{k=1}^{K-1}A[k] + \epsilon}$$

---

## Log-Mel Spectrogram

**ID**: `logmel` :custom-librosa:[^5]

The log-mel spectrum represents how the energy of a sound is distributed across perceptual frequency bands, using the mel scale and a logarithmic dB compression to approximate human loudness perception.

The current implementation builds Slaney-normalized mel filters, computes mel-band energies, converts them with `power_to_db(ref=1.0)`, and applies an 80 dB top range.

$$E_m = \sum_{k=0}^{K-1}H_m[k]P[k]$$

$$L_m = 10\log_{10}(\max(10^{-10}, E_m))$$

$$LogMel[m] = \max(L_m, \max_j L_j - 80)$$

---

## MFCC

**ID**: `mfcc` :custom-librosa:[^3]

MFCCs summarize the shape of a sound's spectrum on a perceptual, mel-based scale, giving a compact representation of timbre and tone color as humans hear it.

The current implementation applies an orthonormal DCT-II basis to the log-mel spectrum.

$$MFCC[i] = \sum_{m=0}^{B-1}D_{i,m}LogMel[m]$$

where

$$D_{i,m} = \alpha_i\cos\left(\frac{\pi}{B}(m+0.5)i\right), \quad \alpha_0=\sqrt{\frac{1}{B}}, \quad \alpha_i=\sqrt{\frac{2}{B}}\;\text{for}\;i>0$$

---

## Chroma

**ID**: `chroma` :custom-librosa:[^4]

Chroma features capture the intensity of the twelve pitch classes (C, C sharp, ..., B) in a sound, representing its harmonic and melodic content independently of octave.

The current implementation builds chroma filters from FFT-bin frequencies, normalizes each filter column, applies octave weighting, rolls the chroma basis by three pitch classes, and sums weighted power.

$$Chroma[c] = \sum_{k=0}^{K-1}W_c[k]P[k]$$

where $W_c[k]$ is the chroma filter weight for pitch class $c$ and FFT bin $k$.

[^2]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-9}$.
[^3]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-5}$.
[^4]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-3}$.
[^5]: `Descriptor` full compatible with [`librosa`](https://librosa.org/).
[^6]: `Descriptor` compatible with [`essentia`](https://essentia.upf.edu/) in order of $10^{-4}$.
[^7]: Brent, W. “Physical and Perceptual Aspects of Percussive Timbre.” PhD thesis, UC San Diego, 2010. https://escholarship.org/uc/item/5bx4j1fj.
