---
icon: lucide/audio-lines
tags:
  - Advanced Configuration
  - Pitch Following
---

# Pitch Following

Use these settings when pitched events are matched too narrowly or too loosely, or when the score needs a different tuning or transposition. OpenScofo compares the audio spectrum with harmonic pitch templates.

## Options

| Keyword | Default | Accepted value | Effect |
| --- | --- | --- | --- |
| `TUNINGA4` | `440` Hz | Positive number | Reference frequency for converting subsequent written pitches to Hz. |
| `TRANSPOSE` | `0` | Number of semitones | Adds a semitone offset to subsequent named pitches. Values outside `-36` to `36` produce a warning but are still stored. |
| `PITCHTEMPLATESIGMA` | `0.5` | `0` to `1`, inclusive | Width of the harmonic peaks, expressed on a semitone scale. Larger values broaden the pitch template. |
| `PITCHTEMPLATEHARMONICS` | `10` | Positive integer | Maximum number of harmonics used to build each pitch template, including the fundamental. |

Tuning and transposition are applied while reading pitches. Set them before the affected events. Pitch-template settings are global, so keep them at the beginning of the score.

The actual template width has a minimum determined by FFT bin spacing: `PITCHTEMPLATESIGMA 0` does not create a zero-width template. Harmonics above Nyquist or with negligible weight are omitted, so the harmonic count is an upper limit.

## Practical start

Confirm `TUNINGA4` and `TRANSPOSE` first. If expressive intonation or wide vibrato causes missed notes, test a slightly larger sigma. A wider template can also make neighboring pitches harder to distinguish. Compare the same passage before and after each change.

```openscofo
BPM 72
TUNINGA4 442
TRANSPOSE -12
PITCHTEMPLATESIGMA 0.8
PITCHTEMPLATEHARMONICS 8

NOTE C5 1
NOTE D5 1
```

Here the named notes are followed one octave below their written pitches, using A4 = 442 Hz.

YIN pitch estimation and spectral rolloff are separate [descriptor settings](descriptors.md#yin-pitch-estimation). They do not set the width of the follower's harmonic pitch templates. For analysis window and update sizes, see [Audio analysis](descriptors.md#audio-analysis).
