---
icon: material/language-javascript
tags:
  - Host Integration
  - JavaScript
---

# JavaScript

## Overview

Use the JavaScript/Emscripten bindings for browser-based score following, real-time processing, and offline descriptor analysis.

<release interface="JavaScript">Loading Releases</release>

## Installation

Load the generated OpenScofo JavaScript/WASM module in your web project. See the [descriptor test page](../descriptors/testing){:target="_blank"} for a working browser example.

## Minimal Example

```javascript
OpenScofo().then((mod) => {
    const scofo = new mod.OpenScofo(48000, 2048, 512);
    scofo.load_score(scoreText);
});
```

## Reference Table

### API

| API | Purpose |
| --- | --- |
| `new mod.OpenScofo(sampleRate, blockSize, hopSize)` | create an engine |
| `load_score(scoreText)` | load score text |
| `process_block(audioBuffer)` | process one audio frame |
| `set_db_threshold(value)` | set silence threshold |
| `set_tuning(a4)` | set A4 tuning |
| `set_current_event(index)` | force score event |
| `set_current_section(section)` | reset to the first event of a named section |
| `load_onnx_model(path)` | load ONNX model |
| `get_live_bpm()` | return estimated BPM |
| `get_event_index()` | return current event index |
| `get_states()` | return Markov states |
| `get_audio_description()` | return descriptors |

### Types

| Type | Use |
| --- | --- |
| `VectorDouble` | numeric buffers passed to the engine |
| `Description` | descriptor values for a block |
| `MarkovState` | score-following state |

## Score Actions

Use the JavaScript API to inspect score action data after processing blocks. Host behavior is application-defined; see [Platform Integrations](index/#sendto-behavior).

## Descriptors

`Description` exposes values such as `mfcc`, `logmelspectrum`, `chroma`, `onset`, `silence`, `loudness`, `flux`, `centroid`, `flatness`, `zcr`, `pitch`, `pitch_confidence`, `db`, and `rms`.

## Complete Example

```javascript
OpenScofo().then((mod) => {
    const scofo = new mod.OpenScofo(44100, 1024, 512);
    scofo.load_score(scoreText);

    const buffer = new mod.VectorDouble();
    // Fill buffer with audio samples.

    scofo.process_block(buffer);

    const bpm = scofo.get_live_bpm();
    const event = scofo.get_event_index();
    const desc = scofo.get_audio_description();
});
```

## Remarks

Errors from the internal logger are forwarded to JavaScript: `error` and `critical` become exceptions; `info` and `debug` print to the console. ONNX loading requires correct file access in the Emscripten environment.
