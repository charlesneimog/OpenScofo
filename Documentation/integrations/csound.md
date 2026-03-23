# CSound

`OpenScofo` is a Csound opcode for real-time score following and audio descriptor extraction.

!!! danger "This is a pre-alpha version!"

## Syntax

```csound
kEvent, kBPM, kTrig, kDesc OpenScofo ain, SMode, SArg
```

## Description

The `OpenScofo` opcode operates in two distinct modes depending on the `SMode` argument:

1. **Score Mode:** Tracks an audio signal against a text-based score, outputting the current event index and live BPM estimation.
2. **Descriptor Mode:** Analyzes the audio signal and outputs a specific scalar audio descriptor (e.g., `RMS`, `crest`, and others).

## Loading `OpenScofo`

To load the `OpenScofo` opcode, you use the following code:

``` csound
<CsOptions>
--opcode-lib=/path/to/OScofoCSound.so
</CsOptions>
```

For MacOS extension is `dylib`, for Windows extension is `dll`.

## Initialization (Inputs)

* **`ain`** *(a-rate)*: The input audio signal to be analyzed.
* **`SMode`** *(String)*: The operational mode. Use `"score"` for score following, or `"desc"` (or `"descriptor"`) for audio feature extraction.
* **`SArg`** *(String)*: 
    * If `SMode` is `"score"`, this must be the **file path** to your OpenScofo text score.
    * If `SMode` is `"desc"`, this must be the **name ID of the descriptor** (e.g., `"rms"`, `"pitch"`, etc.), check complet list [Descriptors](../descriptors/index.md). 


!!! warning "Just scalar descriptors for now"

    If a non-scalar descriptor (like MFCC or Chroma) is requested, it safely falls back to `rms`.

## Performance (Outputs)

* **`kEvent`** *(k-rate)*: The current event index in the score (Active in Score Mode).
* **`kBPM`** *(k-rate)*: The estimated live BPM (Active in Score Mode).
* **`kTrig`** *(k-rate)*: A trigger signal that outputs `1` when a processing block finishes/updates, and `0` otherwise. Useful for triggering sub-instruments.
* **`kDesc`** *(k-rate)*: The float value of the requested audio descriptor (Active in Descriptor Mode).

---

## Examples

### 1. Score Following Mode

Loads a score file and tracks the audio, printing the event and BPM whenever it updates.

```csound
instr 1
    a1 diskin2 "/path/to/audio.wav", 1, 0, 0
    
    ; Run OpenScofo in "score" mode
    kEvent, kBPM, kTrig, kDesc OpenScofo a1, "score", "/path/to/score.txt"
    
    ; Print updates only when the trigger is 1
    if (kTrig == 1) then
        printf "Event: %d | BPM: %.2f\n", changed(kEvent), kEvent, kBPM
    endif
    
    out a1
endin
```

### 2. Audio Descriptor Mode
Extracts the RMS value from the incoming audio signal.

```csound
instr 2
    a1 diskin2 "/path/to/audio.wav", 1, 0, 1
    
    ; Run OpenScofo in "desc" mode, requesting Spectral Crest descriptor
    kEvent, kBPM, kTrig, kDesc OpenScofo a1, "desc", "crest"
    
    ; Print the descriptor value when the trigger is 1
    if (kTrig == 1) then
        printf "RMS Value: %f\n", changed(kDesc), kDesc
    endif
    
    out a1
endin
```
