---
tags:
  - Audio Analysis
---

# Descriptors for Musicians

An audio descriptor turns a quality of sound, such as loudness, brightness, or noisiness, into numbers. OpenScofo measures these qualities in short moments of audio, so the values change as you play.

The table uses OpenScofo's core descriptor IDs. Some descriptors range from 0 to 1, others use Hz or dB, and some return a list of numbers. “Higher” and “lower” describe tendencies: compare sounds recorded with the same microphone, gain, and analysis settings.

| Descriptor name | OpenScofo ID | What it tells you musically |
| --- | --- | --- |
| Zero-crossing rate | `zcr` | How often the waveform crosses zero. Low, smooth tones tend toward **0**; breath, hiss, and higher tones usually give higher values. Near **1** means the waveform changes sign almost every sample; ordinary white noise is typically around **0.5**, not 1. |
| RMS level | `rms` | The overall strength of the sound. A crescendo raises it; a diminuendo lowers it. **0** means silence. Microphone gain also changes it. |
| Level in decibels | `db` | Sound level on a logarithmic scale. A value of **−20 dB** is louder than **−60 dB**. These are digital signal levels, not the sound pressure level in the room. |
| Loudness | `loudness` | A level estimate that gives different frequencies different weight to better reflect hearing. Higher, less negative values mean louder sound; lower values mean quieter sound. |
| Maximum spectral amplitude | `maxamp` | The strength of the strongest frequency component. Higher values mean a stronger peak; lower values mean a weaker peak. This is not the peak of the audio waveform. |
| Silence probability | `silence` | How quiet the sound is according to OpenScofo. Near **1** means likely silence; near **0** means sound is present. Very soft playing can resemble silence. |
| Pitch | `yin` | The estimated pitch of a single note in **Hz**: A4 is about **440 Hz**. Higher values mean higher notes. Chords and noisy techniques may not give a clear single pitch. |
| Pitch confidence | `yin_confidence` | How clearly the sound repeats like a pitched tone. Near **1** means a strong repeating pattern; near **0** means an uncertain pitch. High confidence does not guarantee the correct octave. |
| Onset strength | `onset` | Evidence of a new attack, such as a pluck, tongue attack, or drum hit. Peaks suggest new events; lower values suggest less attack activity. The current descriptor is an attack-strength value, not simply a 0/1 switch. |
| Spectral centroid | `centroid` | A useful clue to **brightness**, measured in Hz. Higher values often mean a brighter or sharper sound; lower values often mean a darker or rounder sound. It is not the played note's pitch. |
| Centroid movement | `centroid_velocity` | How much brightness changes between two moments. Near **0** means stable brightness; higher values mean a stronger change, whether brighter or darker. |
| Spectral rolloff | `rolloff` | The frequency below which most of the sound's energy lies. Higher values suggest more energy extending into high frequencies; lower values suggest a darker sound. Measured in Hz. |
| High-frequency ratio | `hfr` | How much of the spectrum's magnitude lies in its upper frequency region. Near **0** means little high-frequency content; nearer **1** means most content is there, as in a very high hiss. The dividing frequency depends on the sample rate. |
| Spectral spread | `spread` | How widely the sound's frequency content is spread, in Hz. Lower values mean a narrow concentration, like a pure tone; higher values mean a wider range, like a broad noisy sound. |
| Spectral spread variance | `spread_variance` | Another measure of frequency width, scaled to the analyzed frequency range. Lower values mean concentrated content; higher values mean more widely separated content. |
| Spectral flatness | `flatness` | Whether the spectrum has distinct peaks or resembles noise. Near **0** means strong tonal peaks (clear harmonic serie in the sound); toward **1** means a more even, noise-like spectrum. Useful for comparing a clear note with an airy or noisy technique. |
| Spectral entropy | `entropy` | How widely energy is shared among frequencies. Lower values mean energy concentrated in a few places; higher values mean energy distributed across many frequencies, often in noisy sounds. It is not a 0–1 scale. |
| Harmonicity | `harmonicity` | In OpenScofo, how much one frequency component dominates. Near **1** means one very strong component; near **0** means many components share the sound. A rich, clearly pitched tone can still have a low value. |
| Spectral crest | `crest` | How much the strongest frequency peak stands above the average. High values mean a prominent peak; values near **1** mean a nearly even spectrum. |
| Spectral flux | `flux` | How much the spectrum changes from one moment to the next. Low values suggest a steady tone; high values suggest an attack, a sudden dynamic change, or a timbral change. Louder input can increase it. |
| Spectral irregularity | `irregularity` | How uneven neighboring parts of the spectrum are. Higher values mean more pronounced local peaks and dips; lower values mean a smoother shape. It also changes with recording level, so it is not a direct noise meter. |
| Spectral standard deviation | `stddev` | How unevenly frequency content is distributed. Near **0** means a very even spectrum; higher values mean some frequencies stand out more strongly. |
| Spectral slope | `slope` | The overall tilt from low to high frequencies. Negative values mean content generally falls toward the highs; positive values mean it rises. Near **0** means little overall tilt. |
| Spectral skewness | `skewness` | Which way the spectrum stretches away from its main concentration. Positive values mean a longer tail toward high frequencies; negative values mean a longer tail toward low frequencies. Useful for comparing timbres with similar brightness. |
| Spectral kurtosis | `kurtosis` | How strongly distant frequency content stands out around the main concentration. Higher values can reveal a concentrated sound with faraway components; lower values suggest a broader, flatter distribution. It can be negative and has no simple “bright/dark” meaning. |
| Magnitude spectrum | `magnitude` | A list showing the strength of each frequency region. Peaks reveal the fundamental, overtones, or resonances; a broad spread can reveal noise. Each entry describes a different frequency. |
| Power spectrum | `power` | Like the magnitude spectrum, but with strengths squared, giving stronger components more weight. Returns a list, not one overall loudness value. |
| Log-mel spectrum | `logmel` | A list of levels in frequency bands arranged on a scale related to hearing. Larger values in the upper bands suggest more high-frequency content; larger values in the lower bands suggest more low-frequency content. |
| MFCCs | `mfcc` | A compact set of numbers describing tone color. Useful for learning differences between normal tone, breath, pizzicato, or other techniques. There is no single “high MFCC means brighter” rule; the pattern matters. |
| Chroma | `chroma` | A profile of pitch classes, usually C through B, combining octaves. A stronger entry means more evidence for that pitch class. Useful for harmonic content; it does not directly identify the octave or name a chord. |
| Extended-technique estimate | `ext` | A built-in estimate combining spectral change with pitch uncertainty. Higher values favor changing sounds without a clear pitch; lower values favor stable pitched sounds. It does not identify a particular technique by name. |
| Model output | `onnx` | The values returned by your trained recognition model. Their meaning and order depend on the model: for a classifier, a higher class score may indicate stronger evidence for that technique. Outputs are not always probabilities. |

!!! tip "Listen while watching the values"
    Try a sustained note, a crescendo, a breath sound, and a short attack in the [online descriptor test](testing/index.html). Compare how the values move before choosing a threshold for your own playing.

For formulas and settings, see [Amplitude](amplitude.md), [Time and Pitch](time-pitch.md), [Onset](onset.md), and [Spectral Descriptors](spectral.md). The [librosa zero-crossing reference](https://librosa.org/doc/0.11.0/generated/librosa.feature.zero_crossing_rate.html) also explains the fraction measured by `zcr`.
