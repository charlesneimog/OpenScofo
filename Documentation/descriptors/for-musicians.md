---
icon: lucide/music
tags:
  - Audio Analysis
hide:
  - toc
---

# Descriptors for Musicians

An audio descriptor turns a quality of sound, such as loudness, brightness, or noisiness, into numbers. OpenScofo measures these qualities in short moments of audio, so the values change as you play.

The table uses OpenScofo's core descriptor IDs. Some descriptors range from 0 to 1, others use Hz or dB, and some return a list of numbers. “Higher” and “lower” describe tendencies: compare sounds recorded with the same microphone, gain, sample rate, FFT size, and analysis settings. Interpret timbre descriptors while sound is present; silence and very quiet input can give misleading values.

| Descriptor name | OpenScofo ID | What it tells you musically |
| --- | --- | --- |
| Zero-crossing rate | `zcr` | How often the waveform crosses zero, on a **0–1** scale. Low, smooth tones tend to give lower values; breath, hiss, and higher tones often give higher values. Sample rate, near-zero threshold, and padding settings affect the result, so noise has no single fixed value. |
| RMS level | `rms` | The overall strength of the sound. A crescendo usually raises it; a diminuendo lowers it. **0** means silence. Microphone gain also changes it, and it is not capped at 1. |
| Level in decibels | `db` | Sound level on a logarithmic scale. A value of **−20 dB** is louder than **−60 dB**. These are digital signal levels, not the sound pressure level in the room. |
| Loudness | `loudness` | A perceptually weighted level estimate based on ITU-R BS.1770 weighting. Higher values generally suggest louder sound; lower values suggest quieter sound. Values are often negative and depend on recording gain. |
| Maximum spectral amplitude | `maxamp` | The strength of the strongest frequency component. Higher values mean a stronger peak; lower values mean a weaker peak. This is not the peak of the audio waveform. |
| Silence probability | `silence` | A **0–1** estimate based on weighted sound level. Near **1** means likely silence; near **0** means sound is present. Very soft playing can resemble silence, and microphone gain affects the result. |
| Pitch | `yin` | The estimated pitch of a single note in **Hz**: A4 is about **440 Hz**. Higher values mean higher notes. Chords and noisy techniques may not give a clear single pitch. |
| Pitch confidence | `yin_confidence` | A **0–1** estimate of how clearly the sound repeats like a pitched tone. Near **1** means a strong repeating pattern; near **0** means an uncertain pitch. High confidence does not guarantee the correct octave. |
| Onset detection signal | `onset` | Helps locate attacks such as plucks, tongue attacks, and drum hits. In the current implementation, peaks suggest possible attacks; the output can be negative or exceed 1. It is the changing detection signal, not the detector's separate yes/no result. |
| Spectral centroid | `centroid` | A useful clue to **brightness**, measured in Hz. Higher values often mean a brighter or sharper sound; lower values often mean a darker or rounder sound. It is not the played note's pitch. |
| Centroid movement | `centroid_velocity` | How much the spectral centroid changes between two moments, in **Hz per analysis step**. Near **0** means a stable centroid; higher values suggest a stronger brightness change in either direction. The time between analysis steps affects the value. |
| Spectral rolloff | `rolloff` | The frequency below which a chosen share of the sound's energy lies, **85% by default**. Higher values often suggest more high-frequency content; lower values often suggest a darker sound. Measured in Hz; the cutoff setting changes the result. |
| High-frequency ratio | `hfr` | The share of the spectrum's magnitude above OpenScofo's high-frequency boundary, on a **0–1** scale. Near **0** means little content above it; nearer **1** means most content is there. The boundary depends on the sample rate and FFT size. |
| Spectral spread | `spread` | How widely the sound's frequency content is spread, in Hz. Lower values mean a narrow concentration, like a pure tone; higher values mean a wider range, like a broad noisy sound. |
| Spectral spread variance | `spread_variance` | Another measure of frequency width, scaled to the analyzed frequency range. Lower values mean concentrated content; higher values mean more widely separated content. |
| Spectral flatness | `flatness` | Whether the spectrum has distinct peaks or is more even, on a **0–1** scale. Low values mean distinct peaks dominate, as in tonal or resonant sounds; high values mean a more even, noise-like spectrum. Low flatness does not require a harmonic series. Silence can also give a value near 1. |
| Spectral entropy | `entropy` | How widely energy is shared among frequencies. Lower values mean energy concentrated in a few places; higher values mean energy distributed across many frequencies, often in noisy sounds. It can exceed 1, and its upper limit depends on FFT size. |
| Harmonicity | `harmonicity` | OpenScofo measures how strongly one frequency component dominates the others, on a **0–1** scale. Near **1** means strong dominance; low values mean little dominance. It does not measure how well components form a harmonic series. A rich, clearly pitched sound may still have a low value. |
| Spectral crest | `crest` | How much the strongest frequency peak stands above the average. High values mean a prominent peak; values near **1** mean a nearly even spectrum. It can exceed 1 and depends on FFT size; silence returns 0. |
| Spectral flux | `flux` | How much the spectrum changes from one moment to the next. Low values suggest a steady tone; high values suggest an attack, a sudden dynamic change, or a timbral change. Louder input can increase it. |
| Spectral irregularity | `irregularity` | How uneven neighboring parts of the spectrum are. At the same recording level, higher values suggest more pronounced local peaks and dips; lower values suggest a smoother shape. Values can be negative and change with gain and FFT size; this is not a direct noise meter. |
| Spectral standard deviation | `stddev` | How unevenly frequency content is distributed. During audible sound, values near **0** suggest an even spectrum; higher values suggest some frequencies stand out more strongly. The scale depends on FFT size, and near-silent input can distort this reading. |
| Spectral slope | `slope` | The overall tilt from low to high frequencies. Negative values mean content generally falls toward the highs; positive values mean it rises. Near **0** means little overall tilt. |
| Spectral skewness | `skewness` | Which way the spectrum stretches away from its main concentration. Positive values mean a longer tail toward high frequencies; negative values mean a longer tail toward low frequencies. Useful for comparing timbres with similar brightness. |
| Spectral kurtosis | `kurtosis` | How strongly distant frequency content stands out around the main concentration. Higher values can reveal a concentrated sound with faraway components; lower values suggest a broader, flatter distribution. It can be negative and has no simple “bright/dark” meaning. |
| Magnitude spectrum | `magnitude` | A list showing the strength of each frequency region. Peaks can reveal the fundamental, overtones, or resonances; a broad spread can suggest noise. Values depend on gain and FFT size and can exceed 1. |
| Power spectrum | `power` | Like the magnitude spectrum, but with strengths squared, giving stronger components more weight. Returns a list, not one overall loudness value. Values depend on gain and FFT size and can exceed 1. |
| Log-mel spectrum | `logmel` | A list of levels in dB for frequency bands arranged on a scale related to hearing. Larger values in the upper bands suggest more high-frequency content; larger values in the lower bands suggest more low-frequency content. Values can be negative and change with gain and analysis settings. |
| MFCCs | `mfcc` | A compact set of numbers describing tone color. Useful for learning differences between normal tone, breath, pizzicato, or other techniques. Values can be positive or negative; there is no single “high MFCC means brighter” rule. The pattern matters. |
| Chroma | `chroma` | A profile of pitch classes, usually C through B, combining octaves. A stronger entry means more evidence for that pitch class. Values depend on gain and analysis settings and are not limited to 0–1. It does not directly identify the octave or name a chord. |
| Extended-technique estimate | `ext` | A built-in estimate combining spectral change with pitch uncertainty, mapped onto a bounded **0–1** scale. Higher values favor changing sounds without a clear pitch; lower values favor stable pitched sounds. Gain affects it. This is separate from ONNX classification and does not identify named techniques or give a trained class probability. |
| Model output | `onnx` | The values returned by your trained recognition model. Their meaning and order depend on the model: for a classifier, a higher class score may indicate stronger evidence for that technique. Outputs are not always probabilities. |

!!! tip "Listen while watching the values"
    Try a sustained note, a crescendo, a breath sound, and a short attack in the [online descriptor test](testing/index.html). Compare how the values move before choosing a threshold for your own playing.

For formulas and settings, see [Amplitude](amplitude.md), [Time and Pitch](time-pitch.md), [Onset](onset.md), and [Spectral Descriptors](spectral.md). The [librosa zero-crossing reference](https://librosa.org/doc/0.11.0/generated/librosa.feature.zero_crossing_rate.html) also explains the fraction measured by `zcr`.
