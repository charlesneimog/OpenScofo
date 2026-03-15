# Amplitude Descriptors


## `Root Mean Square (RMS)` 

:custom-librosa:[^1] :custom-essentia:

A measure of the amplitude (energy) of the current audio frame. It represents how loud the sound is within a short time window. Higher RMS values indicate louder sounds, while lower values indicate quieter sounds or silence.


<div class="grid cards" style="font-weigth:bold" markdown>
-   ??? equation "Equation"

        $RMS = \sqrt{\frac{1}{N} \sum_{n=0}^{N-1} x[n]^2}$

-   ??? note "Notes"

        Root mean square of the amplitude representing the energy of the time-domain signal $x[n]$ of length $N$:
</div>

---

## `dB`

A logarithmic measure of sound level derived from the signal amplitude. Unlike RMS, which measures the raw energy of the signal, dB expresses this energy on a logarithmic scale that better reflects how humans perceive changes in loudness.

<div class="grid cards" style="font-weigth:bold" markdown>
-   ??? equation "Equation"

        $$L_{dB} = 20 \log_{10}(RMS)$$

-   ??? note "Notes"
        
        Neither `librosa` or `essentia` implemented `dB`, but as it uses `RMS`, once you convert `RMS` to `dB` it is compatible with both.

    
</div>

---

## `Max Amplitude `

Maximum normalized spectral amplitude detected in the current frame of audio.

<div class="grid cards" style="font-weigth:bold" markdown>
-   ??? equation "Equation"

        $$MaxAmp = \max_{k} |X[k]|$$
        Where $X[k]$ is an FFT Bin.

-   ??? note "Notes"
</div>


## `Loudness`

An estimate of perceived sound intensity based on psychoacoustic models of human hearing. Unlike dB, it applies perceptual models and filters derived from psychoacoustic studies to approximate how humans actually perceive loudness.

<div class="grid cards" markdown>
-   ??? equation "Equation"

        $L = -0.691 + 10 \log_{10}\left(\frac{1}{N}\sum_{n=0}^{N-1} y[n]^2\right)$

        where \(y[n]\) denotes the audio samples after applying the filtering stage defined in the [`ITU‑R BS.1770`](https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.1770-5-202311-I!!PDF-E.pdf){:target="_blank"} recommendation.  

        The term \(N\) represents the number of samples in the analyzed frame.  
        This formulation corresponds to the energy-based loudness estimate used in the loudness measurement procedure defined by the standard.
    
-   ??? note "Notes"

        The `loudness` descriptor implemented in `essentia` is based on a simplified perceptual model derived from signal energy with a power-law compression. While computationally inexpensive, it does not incorporate perceptual frequency weighting or the measurement procedure defined in modern broadcast loudness standards.

        `OpenScofo` instead implements loudness estimation following the methodology described in [`ITU‑R BS.1770`](https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.1770-5-202311-I!!PDF-E.pdf){:target="_blank"}. This approach applies perceptual filtering prior to the energy calculation and expresses the result in a logarithmic scale, which aligns with the methodology adopted in contemporary loudness measurement practices for audio production and broadcasting.

</div>


---

## `Silence Probability`

Probability that the current frame corresponds to silence, derived from Loudness ($L$) via a logistic function where $\alpha = 0.25$ and $L_0 = -60.0$:

<div class="grid cards" markdown>
-   ??? equation "Equation"

        $P_{silence} = \frac{1}{1 + e^{\alpha (L - L_0)}}$
    
-   ??? note "Notes"


</div>


---
[^1]: `Descriptor` compatible with [`librosa`](https://librosa.org/) in order of $10^{-9}$.


