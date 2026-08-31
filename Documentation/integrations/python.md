---
icon: material/language-python
tags:
  - Host Integration
  - Python
---

# Python

## Overview

Use the Python bindings for development, validation, descriptor extraction, and research workflows.

## Installation

```bash
pip install OpenScofo
```

## Minimal Example

```python
from OpenScofo import OpenScofo

scofo = OpenScofo(48000, 4096, 1024)
scofo.load_score("score.scofo")
```

## Reference Table

### Constructor

| Argument | Meaning |
| --- | --- |
| `sr` | sample rate |
| `fft_size` | FFT/window size |
| `hop` | hop size |

### Core methods

| Method | Purpose |
| --- | --- |
| `load_score(path)` | load a score file |
| `process_block(audio)` | process one NumPy audio block |
| `set_db_threshold(value)` | set silence threshold |
| `set_tuning(value)` | set A4 tuning |
| `set_current_event(event)` | force score position |
| `set_current_section(section)` | reset to the first event of a named section |
| `set_harmonics(value)` | set pitch-template harmonics |
| `set_pitch_template_sigma(value)` | set pitch tolerance |
| `get_live_bpm()` | return estimated BPM |
| `get_event_index()` | return current event index |
| `get_states()` | return score states |
| `get_pitch_template(freq)` | return pitch template |
| `get_block_duration()` | return block duration in seconds |
| `get_audio_description(audio)` | return descriptors for one block |

## Score Actions

Use the binding's score state/action APIs when you need host-side action handling. See the [C++ integration](cpp/) for the underlying action structure.

## Descriptors

```python
desc = scofo.get_audio_description(segment)
```

Common `Description` attributes include `mfcc`, `chroma`, `onset`, `silence_prob`, `loudness`, `spectral_flux`, `spectral_flatness`, `harmonicity`, `db`, `rms`, and `power`.

## Complete Example

```python
import librosa
from OpenScofo import OpenScofo

scofo = OpenScofo(48000, 4096, 1024)
scofo.load_score("score.scofo")

y, _ = librosa.load("audio.wav", sr=48000)
fftsize = 4096
hopsize = 1024

for pos in range(0, len(y) - fftsize, hopsize):
    block = y[pos:pos + fftsize]
    scofo.process_block(block)
```

## Remarks

Use Python for offline validation and training workflows. For real-time performance, prefer Pd, Max, Csound, SuperCollider, or an embedded C++ host.
