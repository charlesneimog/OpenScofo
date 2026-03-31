# CSound

<release interface="CSound">Loading Releases</release>

`OpenScofo` features a Csound opcode dedicated to real-time score following.

!!! danger "This is a pre-alpha version!"

## Loading `OpenScofo`

To load the opcode, reference the compiled library in your `<CsOptions>` block:

```csound
<CsOptions>
--opcode-lib=/path/to/OpenScofo.so
</CsOptions>
```

!!! warning "Use `.dylib` for macOS and `.dll` for Windows."

-----

## Syntax

```csound
kEvent, kBPM, kTrig OpenScofoScore aIn, SScorePath, iFFTSize, iHopSize 
```

### Inputs

  * **`aIn`** *(a-rate)*: The input audio signal (from a file or real-time microphone input).
  * **`SScorePath`** *(String)*: The absolute file path to your OpenScofo text score.
  * **`iFFTSize`** *(i-rate)*: The FFT size. Recommended value is `2048`.
  * **`iHopSize`** *(i-rate)*: The Hop size. Recommended value is `512`.

### Outputs

  * **`kEvent`** *(k-rate)*: The current event index in the score. Rests do not count as events.
  * **`kBPM`** *(k-rate)*: The estimated live BPM based on the player's execution of the score.
  * **`kTrig`** *(k-rate)*: A trigger signal that outputs `1` when a new event index is detected, and `0` otherwise.

-----

## Example Usage

This example loads a score file and tracks the incoming audio, printing the current event and BPM whenever a new event is detected.

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
    ; Audio input from file (or real-time MIC)
    a1 diskin2 "/path/to/audio.wav", 1, 0, 0

    ; Run OpenScofoScore
    kEvent, kBPM, kTrig OpenScofoScore a1, "/path/to/score.txt", 2048, 512
    
    ; Print updates only when kTrig is 1
    printf "Event: %03d | BPM: %.2f\n", kTrig, kEvent, kBPM
    
    out a1
endin

</CsInstruments>

<CsScore>
i1 0 120
</CsScore>
</CsoundSynthesizer>
```
