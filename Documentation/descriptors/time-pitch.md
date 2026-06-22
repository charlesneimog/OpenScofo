# Time-Domain and Pitch Descriptors

These descriptors are not spectral descriptors in the strict sense. They are computed from the waveform itself, or from pitch-period estimation over the audio frame, rather than from measurements such as spectral brightness, spread, or energy distribution across FFT bins.

They are still useful for machine listening because they describe musical qualities that complement the spectral descriptors: noisiness in the waveform, estimated fundamental frequency, and the reliability of that pitch estimate.

---

## `Zero Crossing Rate`
**ID**: `zcr` | :custom-librosa:[^5]

Zero Crossing Rate counts how often the waveform crosses the zero amplitude line, indicating the noisiness or percussiveness of a sound.

The current implementation optionally pads the frame when `ZCRCenter` is enabled, applies the threshold `ZCRThreshold`, and then counts sign changes. With the default `ZCRZeroPos = true`, zero is treated as non-negative through `std::signbit`.

$$ZCR = \frac{1}{N}\sum_{n=1}^{N-1}\mathbb{I}\left[sign(x[n]) \ne sign(x[n-1])\right]$$

---

## `Pitch` & `PitchConfidence`
**ID**: `pitch`

Estimated fundamental frequency and confidence are calculated using the YIN algorithm's cumulative mean normalized difference function (CMNDF).

The current implementation searches for the first CMNDF value below `YINThreshold`, refines the lag by parabolic interpolation, and reports confidence as $1-d'(\tau)$ clamped to $[0,1]$.

$$d'(\tau) = \begin{cases} 1, & \tau = 0 \\ \frac{d(\tau)}{\frac{1}{\tau}\sum_{j=1}^{\tau}d(j)}, & \tau > 0 \end{cases}$$

$$Pitch = \frac{SR}{\tau}$$

$$Confidence = \max(0, \min(1, 1-d'(\tau)))$$

[^5]: `Descriptor` full compatible with [`librosa`](https://librosa.org/).
