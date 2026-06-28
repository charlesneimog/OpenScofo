# CSound

`OpenScofo` provides a Csound opcode for real-time score following. The opcode tracks an incoming audio stream against an OpenScofo score and can trigger Csound instruments from score `sendto` actions.

!!! danger "This is a pre-alpha version!"
    The Csound opcode is still under active development and is **not yet ready for production use**. Expect missing features, bugs, and breaking changes between releases. Please report any issues you encounter.

-----

## Loading `OpenScofo`

Load the compiled plugin in the `<CsOptions>` block:

```csound
<CsOptions>
--opcode-lib=/path/to/OpenScofo.so
</CsOptions>
```

!!! warning "Use `.dylib` for macOS and `.dll` for Windows."
    For Windows use `--opcode-lib=/path/to/OpenScofo.dll`, for macOS use `--opcode-lib=/path/to/OpenScofo.dylib`.

-----

## Opcode Syntax

```csound
kEvent, kBPM, kTrig OpenScofoScore aIn, SScorePath, iFFTSize, iHopSize
```

### Inputs

* **`aIn`** *(a-rate)*: Audio input signal, from a sound file or real-time input.
* **`SScorePath`** *(string)*: Absolute path to an OpenScofo score.
* **`iFFTSize`** *(i-rate)*: FFT size. Recommended value: `2048`.
* **`iHopSize`** *(i-rate)*: Hop size. Recommended value: `512`.

### Outputs

* **`kEvent`** *(k-rate)*: Current score event index. Rests do not count as events.
* **`kBPM`** *(k-rate)*: Estimated live BPM.
* **`kTrig`** *(k-rate)*: Outputs `1` when a new event is detected, otherwise `0`.

-----

## Score Actions

In the Csound wrapper, `sendto` schedules a Csound `i` event. The `sendto` receiver is the Csound instrument name or number.

```text
sendto 2 [0 0.25]
sendto namedPing [0 0.5]
```

These schedule:

```csound
i 2 0 0.25
i "namedPing" 0 0.5
```

Csound `sendto` actions must include at least two arguments: `p2` and `p3`. If either is missing, OpenScofo prints a warning and does not schedule the Csound event.

### Passing p-fields

Use bracket arguments to pass Csound p-fields after the instrument.

```text
sendto 2 [0 0.25 440 0.2]
sendto namedPing [0 0.5 880 0.1]
```

These schedule:

```csound
i 2 0 0.25 440 0.2
i "namedPing" 0 0.5 880 0.1
```

In Csound:

* **`p1`** is the instrument number or name.
* **`p2`** is the start time, relative to now for real-time events.
* **`p3`** is the duration.
* **`p4`, `p5`, ...** are custom parameters for the instrument.

So yes: if you want to control the start time and duration, the first two arguments in the `sendto` list should be `p2` and `p3`.

```text
sendto 2 [p2 p3 p4 p5 ...]
```

For example:

```text
sendto 2 [0 0.25 660 0.15]
```

means: trigger `instr 2` now, play for `0.25` seconds, with `p4=660` and `p5=0.15`.

### Delayed Actions

OpenScofo score delays still work, but they are different from Csound `p2`.

There are two timing layers:

* **OpenScofo `delay`** waits before sending the event to Csound.
* **Csound `p2`** waits after Csound receives the event, before starting the instrument.

So the sounding time is:

```text
OpenScofo event time + OpenScofo delay + Csound p2
```

If you want the Csound instrument to start exactly when the OpenScofo action fires, use `p2 = 0`.

```text
delay 0.5 tempo sendto 2 [0 0.25 660 0.15]
```

This means:

1. OpenScofo waits `0.5` beats at the current followed tempo.
2. OpenScofo sends this Csound event:

```csound
i 2 0 0.25 660 0.15
```

3. Because `p2 = 0`, `instr 2` starts immediately when Csound receives it.
4. Because `p3 = 0.25`, `instr 2` plays for `0.25` seconds.

If you also set a non-zero `p2`, the two delays add together:

```text
delay 0.5 tempo sendto 2 [0.1 0.25 660 0.15]
```

This means:

1. OpenScofo waits `0.5` beats.
2. OpenScofo sends:

```csound
i 2 0.1 0.25 660 0.15
```

3. Csound waits another `0.1` seconds because `p2 = 0.1`.
4. Then `instr 2` starts and plays for `0.25` seconds.

As a rule of thumb: use OpenScofo `delay` for score-relative or tempo-relative timing, and use Csound `p2` only when you want a small local offset inside Csound.

-----

## Example Score Actions

```text
BPM 100

NOTE Bb4 1.5
    sendto 2 [0 0.10 466.16 0.10]
    delay 0.25 tempo sendto namedPing [0 0.09 932.33 0.07]
    delay 0.5 tempo sendto 2 [0 0.08 587.33 0.08]

REST 1
    // Space for the flute: no Csound action here.

NOTE D4 2
    sendto 2 [0 0.10 293.66 0.12]
    delay 0.125 tempo sendto 2 [0 0.08 369.99 0.10]
    delay 0.25 tempo sendto namedPing [0 0.08 739.99 0.06]
```

-----

## Example `.csd`

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
    if iFreq <= 0 then
        iFreq = 440
    endif

    iAmp = p5
    if iAmp <= 0 then
        iAmp = 0.2
    endif

    aEnv linseg 0, 0.01, iAmp, p3 - 0.02, iAmp, 0.01, 0
    aTone poscil aEnv, iFreq
    out aTone
endin

instr namedPing
    iFreq = p4
    if iFreq <= 0 then
        iFreq = 880
    endif

    iAmp = p5
    if iAmp <= 0 then
        iAmp = 0.12
    endif

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
