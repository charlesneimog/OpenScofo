# JavaScript

<release interface="JavaScript">Loading Releases</release>

Check the JavaScript module using the [OpenScofo descriptors](../descriptors/testing){target="_blank"}.


This module exposes `OpenScofo` to JavaScript via Emscripten bindings. It enables browser-based score following, real-time processing, and offline analysis.

## Overview

The bindings wrap the core C++ API and expose:

- Core `OpenScofo` engine
- Audio feature descriptors
- Markov states and score tracking
- Configuration parameters
- Processing functions

Errors from the internal logger are forwarded to JavaScript:
- `error` / `critical` → thrown as exceptions
- `info` / `debug` → printed to console

---

## Module Usage

After compiling with Emscripten, the module can be used as:

```javascript
const scofo = new Module.OpenScofo(sampleRate, blockSize, hopSize);
```

---

## Types

### VectorDouble

```js
let v = new Module.VectorDouble();
```

Used for arrays of numeric data exchanged with the engine.

## Description

Audio descriptor container (per block).

```js
desc.mfcc
desc.logmelspectrum
desc.chroma
desc.onset
desc.silence
desc.ext
desc.max_amp
desc.loudness
desc.flux
desc.irregularity
desc.crest
desc.centroid
desc.magnitude
desc.spreadhz
desc.spread_variance
desc.flatness
desc.skewness
desc.slope
desc.kurtosis
desc.centroid_velocity
desc.hfr
desc.harmonicity
desc.zcr
desc.stddev
desc.pitch
desc.pitch_confidence
desc.db
desc.rms
desc.magnitude
```

---

## MarkovState

Represents a state in the score-following model.

```js
state.position
state.type
state.markov
state.forward
state.bpm_expected
state.bpm_observed
state.onset_expected
state.onset_observed
state.phase_expected
state.phase_observed
state.ioi_phi_n
state.ioi_hat_phi_n
state.audio_states
state.duration
state.line
```

---

## OpenScofo Class

### Constructor

```js
var oscofo = null;
OpenScofo().then(async (mod) => {
    wasmModule = mod;
    oscofo = new mod.OpenScofo(48000, 2048, 512);
});
```

---

### Score

```js
scofo.parse_score(scoreText)
```

---

### Configuration

```js
scofo.set_db_threshold(value)
scofo.set_tuning(a4)
scofo.set_current_event(index)
scofo.load_onnx_model(path)
scofo.set_amplitude_decay(value)
scofo.set_harmonics(count)
scofo.set_pitch_template_sigma(value)
```

---

### Processing

```js
scofo.process_block(audioBuffer)
```

* `audioBuffer`: `VectorDouble` or compatible array
* Processes one audio frame

---

### Getters

```js
scofo.get_live_bpm()
scofo.get_event_index()
scofo.get_states()
scofo.get_pitch_template()
scofo.get_block_duration()
scofo.get_audio_description()
```

---

## Minimal Example

```js
const scofo = new Module.OpenScofo(44100, 1024, 512);

scofo.parse_score(scoreText);

let buffer = new Module.VectorDouble();
// fill buffer with audio samples

scofo.process_block(buffer);

let bpm = scofo.get_live_bpm();
let state = scofo.get_event_index();
let desc = scofo.get_audio_description();
```

---

## Notes

* Designed for real-time and browser contexts
* Compatible with Web Audio pipelines
* Can be used with AudioWorklets for low-latency processing
* ONNX model loading requires proper file access setup in Emscripten

The API closely mirrors the C++ interface, enabling consistent behavior across environments.
