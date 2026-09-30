---
tags:
  - Advanced Configuration
  - Tempo
---

# Tempo and Synchronization

Use these settings when OpenScofo recognizes the right events but its timing is too rigid or reacts too strongly to the performer.

## Options

| Keyword | Default | Accepted value | Effect |
| --- | --- | --- | --- |
| `BPM` | Must be declared before the first event | Use a value of at least `1` | Sets the expected tempo in beats per minute and inserts an internal `FIRSTEVENT` boundary. |
| `PHASECOUPLING` | `0.5` | `0` to `2`, inclusive | Strength of the phase correction in the temporal model. |
| `SYNCSTRENGTH` | `0.5` | `0` to `1`, inclusive | Strength of synchronization and tempo prediction updates. |
| `TIMETOLERANCE` | Internal value `16`, equivalent to input `0.8` | Clamped to `0` to `1` | Stored on subsequent states, but currently unused by the runtime duration distribution. |

`BPM` below `1` logs an error and clamps the parser's current tempo to `1`; the inserted boundary still receives the supplied value. Use a valid positive tempo rather than relying on this recovery behavior. Out-of-range `PHASECOUPLING` and `SYNCSTRENGTH` values log an error and leave the previous configuration value in place.

## Practical start

Set a realistic `BPM` first. Test a higher `SYNCSTRENGTH` when tempo adaptation is too weak, or a lower value when the follower overreacts. `PHASECOUPLING` controls phase correction; change it separately so you can hear which adjustment helps.

Use sections for passages with different timing settings. The parser records temporal settings on events and section starts; selecting a section restores its starting tempo and temporal model. See [Sections](../config.md#sections) for selection and `SECTIONRESTRICT`.

```openscofo
SECTIONRESTRICT ON

SECTION "A"
BPM 96
SYNCSTRENGTH 0.4
PHASECOUPLING 0.5
NOTE C4 1
NOTE D4 1

SECTION "B"
BPM 72
SYNCSTRENGTH 0.6
PHASECOUPLING 0.5
NOTE E4 2
```

## Current time-tolerance behavior

The parser maps `TIMETOLERANCE t` to the internal value `64 - 60 * clamp(t, 0, 1)`. The initial internal value is `16`.

The current follower's `BuildDistributionCache` does not read that state value. It uses a fixed duration-tail mixture weight of `0.03`. Changing `TIMETOLERANCE` therefore does not currently change duration tolerance during following.
