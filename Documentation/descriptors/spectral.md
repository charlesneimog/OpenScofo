# Spectral Descriptors

Spectral descriptors are used after a `FFT`. For `FFT` we use a Hann window. 

---


## `Spectral Flatness` 
**ID**: `flatness` | :custom-librosa:[^2]

Spectral flatness indicates how noisy versus tonal a sound is. A high flatness means the spectrum is uniform like white noise, while a low flatness shows clear peaks, like a sustained musical note.


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $$
        Flatness = \frac{\exp\left( \frac{1}{K} \sum_{k=0}^{K-1} \ln(|X[k]|^2) \right)}{\frac{1}{K} \sum_{k=0}^{K-1} |X[k]|^2}
        $$

-   ??? note "Notes"


</div>

---

## `Spectral Flux`
**ID**: `flux` |  :custom-essentia:[^6]

Spectral flux measures how quickly the spectrum of a sound changes over time. High flux indicates sudden changes or transients, like drum hits, while low flux corresponds to steady, continuous sounds.


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Flux = \sum_{k=1}^{K-1} \max(0, |X[k]| - |X_{prev}[k]|)$

-   ??? note "Notes"


</div>

---

## `Spectral Irregularity`
**ID**: `irregularity`


Spectral irregularity quantifies how uneven or jagged a spectrum is between adjacent frequency bins. High irregularity indicates complex, inharmonic, or noisy timbres, while low values suggest smooth, harmonic sounds.


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Irregularity = \frac{\sum_{k=1}^{K-1} (|X[k]| - |X[k-1]|)^2}{\sum_{k=0}^{K-1} |X[k]|^2}$

-   ??? note "Notes"

</div>

---

## `Spectral Crest`
**ID**: `crest` | :custom-essentia:[^6]

Spectral crest measures the ratio of the highest spectral peak to the average spectral amplitude. A high crest indicates a tone dominated by strong harmonics or transients, while a low crest corresponds to more even, noise-like spectra.

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Crest = \frac{\max_{k} |X[k]|}{\frac{1}{K} \sum_{k=0}^{K-1} |X[k]|}$

-   ??? note "Notes"

</div>

## `Spectral Skewness`
**ID**: `skewness` | :custom-essentia:[^6]

Spectral kurtosis measures the _peakedness_ or _tailedness_ of the spectral distribution around its centroid. It quantifies how concentrated the spectral energy is in a few frequency bins versus being evenly spread.

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Irregularity = \frac{\sum_{k=1}^{K-1} (|X[k]| - |X[k-1]|)^2}{\sum_{k=0}^{K-1} |X[k]|^2}$

-   ??? note "Notes"

</div>

## `Spectral Kurtois`
**ID**: `kurtosis` | :custom-essentia:[^6]
    
Spectral skewness measures the asymmetry of the spectral distribution around its centroid. It indicates whether the energy is biased toward low or high frequencies.


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Irregularity = \frac{\sum_{k=1}^{K-1} (|X[k]| - |X[k-1]|)^2}{\sum_{k=0}^{K-1} |X[k]|^2}$

-   ??? note "Notes"

</div>

---

## `Spectral RollOff` 
**ID**: `rolloff` | :custom-essentia:[^6]

Spectral rolloff indicates the frequency below which a fixed percentage of a sound’s spectral energy is contained. Higher values make the sound perceptually brighter or sharper, while lower values make it darker or warmer.


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        TODO

-   ??? note "Notes"

        TODO

</div>

---

## `Spectral Entropy` 
**ID**: `entropy` | :custom-essentia:[^6]

Spectral entropy indicates how uniformly a sound’s spectral energy is distributed across frequencies. Higher values make the sound perceptually more noisy or disordered, while lower values make it more tonal or structured.

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        TODO

-   ??? note "Notes"

        TODO

</div>

## `Spectral Centroid` 
**ID**: `centroid` | :custom-librosa:[^3]

Spectral centroid indicates the “center of mass” of a sound’s spectrum. Higher values make the sound perceptually brighter, while lower values make it darker or warmer.

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $$
        Centroid = \frac{\sum_{k=0}^{K-1} f_k |X[k]|}{\sum_{k=0}^{K-1} |X[k]|}
        $$

-   ??? note "Notes"


</div>

---

## `Centroid Velocity`
**ID**: `velocity` | :custom-essentia:[^6]

Centroid velocity measures how quickly the spectral centroid changes over time, reflecting dynamic shifts in brightness or timbre.

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Velocity = |Centroid_t - Centroid_{t-1}|$

-   ??? note "Notes"


</div>

---

## `Spectral Spread` 
**ID**: `spread` | :custom-librosa:[^3]

Spectral spread quantifies how dispersed the energy is around the spectral centroid, indicating whether the sound is focused (narrow) or diffuse (wide) in frequency.

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Spread = \sqrt{\frac{\sum_{k=0}^{K-1} (f_k - Centroid)^2 |X[k]|}{\sum_{k=0}^{K-1} |X[k]|}}$

-   ??? note "Notes"


</div>

---

## `High Frequency Ratio`
**ID**: `hfr`


High Frequency Ratio measures the proportion of energy in the upper part of the spectrum, reflecting the brightness or presence of high-pitched content in a sound.


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $HFR = \frac{\sum_{k=K/4}^{K-1} |X[k]|}{\sum_{k=0}^{K-1} |X[k]|}$

-   ??? note "Notes"


</div>

---

## `Zero Crossing Rate`
**ID**: `zcr` | :custom-librosa:[^5]

Zero Crossing Rate counts how often the waveform crosses the zero amplitude line, indicating the noisiness or percussiveness of a sound.

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $ZCR = \frac{1}{N} \sum_{n=1}^{N-1} \mathbb{I}\{\text{sgn}(x[n]) \neq \text{sgn}(x[n-1])\}$

-   ??? note "Notes"


</div>

---

## `Standard Deviation`
**ID**: `std`

Standard deviation of the normalized spectral power compared to the mean ($\mu = \frac{1}{K}$):

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $StdDev = \sqrt{\frac{1}{K} \sum_{k=0}^{K-1} \left( |X_{norm}[k]| - \mu \right)^2}$

-   ??? note "Notes"


</div>

---

## `Pitch` & `PitchConfidence`
**ID**: `pitch`


Estimated fundamental frequency and confidence, calculated using the YIN algorithm's cumulative mean normalized difference function (CMNDF):


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $d'_t(\tau) = \begin{cases} 1 & \text{if } \tau = 0 \\ \frac{d_t(\tau)}{\frac{1}{\tau} \sum_{j=1}^{\tau} d_t(j)} & \text{otherwise} \end{cases}$

-   ??? note "Notes"


</div>

---

## `Normalized Magnitude`

Vector of magnitude values normalized by the FFT size $N$:

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $|X[k]| = \frac{Power[k]}{N}$

-   ??? note "Notes"


</div>


## `Harmonicity` 
**ID**: `harmonicity` 

!!! warning "Very simple version implemented for now."

Harmonicity measures how well a sound’s spectrum aligns with a harmonic series (integer multiples of a fundamental frequency). High harmonicity indicates a clear pitched sound with harmonically related partials, while low harmonicity indicates inharmonic or noise-like spectra.


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Harmonicity = \frac{\max_{k>0} |X[k]|}{\sum_{k>0} |X[k]|}$

-   ??? note "Notes"

        Probably will be Based on [Yu](https://doi.org/10.1121/1.2356838) (2006);

</div>

---

## `Log-Mel Spectrogram`
**ID**: `logmel` | :custom-librosa:[^5]

The log-mel spectrum represents how the **energy of a sound is distributed across perceptual frequency bands**, using the mel scale and a logarithmic (dB) compression to approximate human loudness perception.

<div class="grid cards" style="font-weigth:bold" markdown>

-    ??? equation "Equation"

        $$
        E_m = \sum_{k=0}^{K-1} H_m(k)\,P[k]
        $$

        $$
        \text{LogMel}[m] = \max\left(L_{min},\; 10 \log_{10}(E_m)\right)
        $$

-    ??? note "Notes"

        - $P[k]$ is the **power spectrum** ($|X[k]|^2$)  
        - $H_m(k)$ is the **mel filterbank** (triangular filters)  
        - $E_m$ is the **energy in mel band $m$**  
        - Log scaling (dB) approximates **human loudness perception**  
        - $L_{min}$ prevents $-\infty$ (numerical stability)  
        - Often followed by a **top-dB clipping** (e.g., 80 dB range)  
        - This is computed **per frame**; stacking over time forms a **mel spectrogram**

</div>


---

## `MFCC` 
**ID**: `mfcc` | :custom-librosa:[^3]

MFCCs summarize the shape of a sound’s spectrum on a perceptual, mel-based scale, giving a compact representation of timbre and tone color as humans hear it.

<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $MFCC[i] = \sum_{m=0}^{M-1} \cos\left( \frac{\pi i (m + 0.5)}{M} \right) \max(L_{min}, 10 \log_{10}(E_m))$

-   ??? note "Notes"

</div>


---

## `Chroma`
**ID**: `chroma` | :custom-librosa:[^4]

Chroma features capture the intensity of the twelve pitch classes (C, C♯, …, B) in a sound, representing its harmonic and melodic content independently of octave.


<div class="grid cards" style="font-weigth:bold" markdown>

-   ??? equation "Equation"

        $Chroma[c] = \sum_{k=0}^{K-1} W_{c}[k] \cdot Power[k]$

-   ??? note "Notes"


</div>


[^2]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-9}$.
[^3]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-5}$.
[^4]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-3}$.
[^5]: `Descriptor` full compatible with [`librosa`](https://librosa.org/).

[^6]: `Descriptor` compatible with [`essentia`](https://essentia.upf.edu/) in order of $10^{-4}$.
