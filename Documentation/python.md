# Python

`OpenScofo` also provides a Python module, which is primarily used for development and research purposes. In addition, this module can be employed for validation tasks, prototyping, and other auxiliary applications, such as experimental evaluation and rapid testing of algorithmic components.


## OpenScofo Class

### `OpenScofo`

``` py
from OpenScofo import OpenScofo

scofo = OpenScofo(48000, 4096, 1024)

```

#### `parse_score`

``` py
scofo.parse_score("myscore.txt")
```

#### `process_block`

``` py
y, _ = librosa.load(path, sr=48000)
pos = 0
fftsize = 4096
hopsize = 1024
while pos + 4096 < len(y):
    segment = y[pos : pos + fftsize]
    desc = scofo.process_block(segment)
    pos += hopsize
```

#### `set_db_threshold`

``` py
scofo.set_db_threshold(-80)
```

#### `set_tuning`

``` py
scofo.set_tuning(442)
```

#### `set_current_event`

``` py
scofo.set_current_event(2)
```

#### `set_amplitude_decay`

Set the the amplitude decay in each harmonics increate.

``` py
scofo.set_amplitude_decay(0.7)
```

#### `set_harmonics`

Set the the number of harmonics used in the pitch template generation.

``` py
scofo.set_harmonics(8)
```

#### `set_pitch_template_sigma`

Set the pitch template sigma for the pitch template generation.

``` py
scofo.set_pitch_template_sigma(1.2)
```

#### `get_live_bpm`

Get the current `BPM` measure by the Large (1999) algorithm.

``` py
bpm = scofo.get_live_bpm()
```

#### `get_event_index`

``` py
score_index = scofo.get_event_index()
```

#### `get_states`

``` py
score_states = scofo.get_states()
```

#### `get_pitch_template`

``` py
score_states = scofo.get_pitch_template(440)
```

#### `get_cqt_template`

``` py
score_states = scofo.get_cqt_template(440)
```

#### `get_audio_description`

``` py
y, _ = librosa.load(path, sr=48000)
segment = y[0: 0 + 4096]
desc = scofo.get_audio_description(segment)
```

