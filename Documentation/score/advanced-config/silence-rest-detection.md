---
icon: lucide/volume-x
tags:
  - Advanced Configuration
  - Silence Detection
---

# Silence/rest Detection

Written rests and short gaps between sounds use silence evidence from the audio analysis. Check this behavior when room noise obscures rests or very quiet playing is mistaken for silence.

## Configuration and current limitation

| Keyword | Default | Accepted value | Current behavior |
| --- | --- | --- | --- |
| `DBTHRESHOLD` | `-60` | Any numeric value; no range check | Stored in configuration and copied into the follower, but not used by the current silence-probability calculation. |

`Score::NewConfig` also handles the legacy spelling `DBTRESHOLD`, but the current language grammar only declares `DBTHRESHOLD`. Use the correctly spelled keyword in scores.

Changing `DBTHRESHOLD` currently does not adjust silence detection. The audio analysis calculates silence probability from filtered loudness `L` using:

```text
silence_probability = 1 / (1 + exp(0.25 * (L + 60)))
```

The midpoint is fixed at `-60`: the probability is `0.5` there, rises for quieter input, and falls for louder input. The follower compares this probability with evidence for sounding events; it is not a hard dB gate.

## Written rests and gaps

A `REST` has a duration in beats and a silence observation. The parser can also insert a zero-duration internal silence state after a `NOTE`, `CHORD`, `PTECH`, or `UTECH` when another sounding event follows in the same section. These internal states retain the preceding score position and allow gaps between events without requiring an extra written rest.

```openscofo
BPM 72

NOTE C4 1
REST 2
NOTE D4 1
```

Use explicit `REST` events for notated silence. Test with the actual microphone, room, and input gain: those affect the loudness used by the current calculation. For the available measurements, see [Amplitude Descriptors](../../descriptors/amplitude.md).
