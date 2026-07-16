---
icon: material/music
tags:
  - Examples
---

# Examples

Copy a pattern into a score and rename receivers to match your host patch.

## Overview

| Area | Pattern | Receiver names |
| --- | --- | --- |
| Interactive electronics | delay answer | `delay`, `echo` |
| Interactive electronics | processing change | `reverb`, `freeze`, `granular` |
| Interactive electronics | crossfade | `dry`, `wet` |
| Performance control | cue sections | `section`, `cue`, `texture` |
| Performance control | extended techniques | `sample`, `noise_texture` |
| Media | video, lights, projection | `video`, `lights`, `projection` |
| AI | technique recognition | `technique`, `model_state` |

## Interactive Electronics

### Delay

```openscofo
BPM 60

NOTE C4 1
    sendto delay [1]

NOTE D4 1
    delay 1 tempo sendto echo [1]
```

### Processing

```openscofo
BPM 60

CHORD (E4 G4 B4) 2
    sendto reverb [0.8]

NOTE A4 1
    sendto freeze [1]

REST 1
    sendto granular [close]
```

### Crossfade

```openscofo
BPM 72

NOTE C4 1
    sendto dry [0.2]
    sendto wet [0.8]

NOTE G4 1
    delay 0.5 tempo sendto wet [0.0]
```

## Performance Control

### Extended Techniques

```openscofo
ONNXMODEL "flute.onnx"
ONNXDESCRIPTORS mfcc zcr centroid spread
BPM 60

PTECH pizz C4 1
    sendto sample [start pizz_echo]

UTECH jet-whistle 2
    sendto noise_texture [open]
```

Requires labels from a trained model. See [AI Models](../ai/).

## Media

```openscofo
BPM 60

NOTE C5 1
    sendto video [scene 2]

NOTE F4 2
    delay 1 tempo sendto lights [fade blue 3]

NOTE A4 1
    sendto projection [/scene 4]
```

Your host patch can translate `projection`, `video`, or `lights` into OSC, MIDI, DMX, or another protocol.

## AI

```openscofo
ONNXMODEL "flute.onnx"
ONNXDESCRIPTORS mfcc logmel zcr centroid flatness hfr
BPM 60

PTECH key-click D4 0.5
    sendto technique [key_click]

UTECH breath 1
    sendto model_state [breath]
```

See also: [Computer Actions](../score/actions/), [Platform Integrations](../integrations/).
