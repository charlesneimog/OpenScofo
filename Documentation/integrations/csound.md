---
icon: custom/csound
tags:
  - Host Integration
  - Csound
---

# Csound

## Overview

Use the Csound opcode for real-time score following. In Csound, `sendto` schedules instruments.

## Installation

Load the compiled plugin in `<CsOptions>`:

```csound
<CsOptions>
--opcode-lib=/path/to/OpenScofo.so
</CsOptions>
```

Use `.dylib` on macOS and `.dll` on Windows.

!!! tip "Install using the installation executable"
    If you use the `OpenScofo` [installer](https://github.com/charlesneimog/OpenScofo/releases){:target = "_blank"} you don't need to use the `CsOptions`.

## Minimal Example

```csound
kEvent, kBPM, kTrig OpenScofoScore aIn, "/path/to/score.scofo", 2048, 512
```

## Reference Table

### Messages / API

| Argument / output | Type | Meaning |
| --- | --- | --- |
| `aIn` | a-rate input | mono signal to analyze |
| `SScorePath` | string | absolute path to a `.scofo` score |
| `iFFTSize` | i-rate | FFT size, usually `2048` |
| `iHopSize` | i-rate | hop size, usually `512` |
| `kEvent` | k-rate output | current score event index |
| `kBPM` | k-rate output | estimated live BPM |
| `kTrig` | k-rate output | `1` when a new event is detected |

## Score Actions

The `sendto` receiver is the Csound instrument name or number. The bracket list becomes p-fields after `p1`.

| Score action | Csound event |
| --- | --- |
| `sendto 2 [0 0.25]` | `i 2 0 0.25` |
| `sendto 2 [0 0.25 440 0.2]` | `i 2 0 0.25 440 0.2` |
| `sendto namedPing [0 0.5 880]` | `i "namedPing" 0 0.5 880` |

`p2` and `p3` are required. Use OpenScofo `delay` for score-relative timing and Csound `p2` for local offsets. For shared action syntax, see [Computer Actions](../score/actions/).

```openscofo
NOTE Bb4 1.5
    sendto 2 [0 0.10 466.16 0.10]
    delay 0.25 tempo sendto namedPing [0 0.09 932.33 0.07]
```

## Descriptors

Descriptor output is not the primary Csound workflow. Use the Vamp, Python, JavaScript, or C++ integrations for descriptor-centered work.

## Complete Example

```csound
<CsoundSynthesizer>
<CsOptions>
--opcode-lib=/path/to/OpenScofo.so
</CsOptions>

<CsInstruments>
sr = 48000
ksmps = 64
nchnls = 1
0dbfs = 1

instr 1
    a1 diskin2 "/path/to/audio.wav", 1, 0, 0
    kEvent, kBPM, kTrig OpenScofoScore a1, "/path/to/score.scofo", 2048, 512
    printf "Event: %03d | BPM: %.2f\n", kTrig, kEvent, kBPM
    out a1
endin

instr 2
    iFreq = p4
    iAmp = p5
    aEnv linseg 0, 0.01, iAmp, p3 - 0.02, iAmp, 0.01, 0
    aTone poscil aEnv, iFreq
    out aTone
endin
</CsInstruments>

<CsScore>
i1 0 30
</CsScore>
</CsoundSynthesizer>
```

!!! warning "`p2` and `p3` are required"
    If an action omits `p2` or `p3`, OpenScofo prints a warning and does not schedule the event.
