# Python

`OpenScofo` provides Python bindings (via `pybind11`) for development, validation, and research workflows.

## Import and create an object

```py
from OpenScofo import OpenScofo

scofo = OpenScofo(48000, 4096, 1024)
```

Constructor arguments:

- `sr` (float): sample rate.
- `fft_size` (float): FFT/window size.
- `hop` (float): hop size.

## `OpenScofo` methods

### Score

#### `parse_score(path)`

- Input: score file path.
- Output: `bool`.

```py
ok = scofo.parse_score("myscore.txt")
```

### Processing

#### `process_block(audio)`

- Input: 1D NumPy array (`float64`) with one audio block.
- Output: `bool` indicating processing success.

```py
import librosa

y, _ = librosa.load(path, sr=48000)
fftsize = 4096
hopsize = 1024
pos = 0
while pos + fftsize <= len(y):
    segment = y[pos: pos + fftsize]
    ok = scofo.process_block(segment)
    pos += hopsize
```

### Configuration

#### `set_db_threshold(value)`

```py
scofo.set_db_threshold(-80)
```

#### `set_tuning(value)`

```py
scofo.set_tuning(442)
```

#### `set_current_event(event)`

```py
scofo.set_current_event(2)
```

#### `set_amplitude_decay(value)`

Sets the amplitude decay for harmonic weighting.

```py
scofo.set_amplitude_decay(0.7)
```

#### `set_harmonics(value)`

Sets the number of harmonics used to build pitch templates.

```py
scofo.set_harmonics(8)
```

#### `set_pitch_template_sigma(value)`

Sets the sigma used for pitch template generation.

```py
scofo.set_pitch_template_sigma(1.2)
```

### Information

#### `get_live_bpm()`

```py
bpm = scofo.get_live_bpm()
```

#### `get_event_index()`

```py
score_index = scofo.get_event_index()
```

#### `get_states()`

```py
states = scofo.get_states()
```

#### `get_pitch_template(freq)`

```py
template = scofo.get_pitch_template(440)
```

#### `get_cqt_template(freq)`

```py
cqt_template = scofo.get_cqt_template(440)
```

#### `get_block_duration()`

```py
seconds = scofo.get_block_duration()
```

#### `get_audio_description(audio)`

```py
y, _ = librosa.load(path, sr=48000)
segment = y[0:4096]
desc = scofo.get_audio_description(segment)
```

## Exposed types

### `Description`

Attributes:

- `mfcc`
- `chroma`
- `onset`
- `silence_prob`
- `spectral_power`
- `norm_spectral_power`
- `pseudo_cqt`
- `loudness`
- `spectral_flux`
- `spectral_flatness`
- `harmonicity`
- `db`
- `rms`
- `power`

### `State`

Attributes:

- `position`
- `type`
- `markov`
- `forward`
- `bpm_expected`
- `bpm_observed`
- `onset_expected`
- `onset_observed`
- `phase_expected`
- `phase_observed`
- `ioi_phi_n`
- `ioi_hat_phi_n`
- `audio_states`
- `duration`
- `line`

### `AudioState`

Attributes:

- `frequency`
- `index`

### `EventType`

Enum values:

- `EventType.REST`
- `EventType.NOTE`
- `EventType.CHORD`
- `EventType.TRILL`
- `EventType.MULTI`

