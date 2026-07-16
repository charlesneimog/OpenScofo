---
icon: material/robot
tags:
  - AI Models
---

# AI Models

Use AI models when you need to recognize sounds such as extended techniques, breath, key clicks, or percussive gestures.

## Use a Model

```openscofo
ONNXMODEL "flute.onnx"
ONNXDESCRIPTORS mfcc zcr centroid spread
BPM 60

PTECH pizz C4 1
    sendto sample [start pizz_echo]

UTECH jet-whistle 2
    sendto noise_texture [open]
```

Rules:

- event labels must match labels used during training;
- descriptor order must match the training order;
- model paths are relative to the score file.

## Prepare Training Data

```text
Flute/
  pizz/
    pizz_C4_01.wav
    pizz_D4_01.wav
  jet-whistle/
    jet_01.wav
    jet_02.wav
```

Each subfolder is a class label. Use short, isolated, representative recordings.

## Train and Use

1. Record examples.
2. Organize them by label.
3. Extract descriptors.
4. Train an ONNX model.
5. Load it with `ONNXMODEL` and `ONNXDESCRIPTORS`.

See also: [Audio Analysis](../descriptors/), [Advanced Configuration](../score/advanced-config/), [AI Model Reference](../descriptors/ai/).
