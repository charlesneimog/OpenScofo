---
tags:
  - Advanced Configuration
---

# Advanced Configuration

Use these settings when the default following or audio analysis needs adjustment. Start with [Configuring a Score](config.md) for the basic syntax and sections, then choose the area you need:

| Topic | Settings and behavior |
| --- | --- |
| [Pitch Following](advanced-config/pitch-following.md) | Tuning, transposition, and harmonic pitch templates |
| [Tempo and Synchronization](advanced-config/tempo-and-synchronization.md) | BPM, phase coupling, synchronization strength, and time tolerance |
| [Silence/rest Detection](advanced-config/silence-rest-detection.md) | Written rests, silence probability, and the current threshold limitation |
| [AI Model Recognition](advanced-config/ai-model-recognition.md) | ONNX models, technique labels, descriptor order, and model metadata |
| [Descriptors](advanced-config/descriptors.md) | Analysis sizes, onset detection, MFCC/logmel, YIN, chroma, and zero crossings |

## Where settings apply

Write configuration keywords in uppercase, followed by a value. Put global analysis and model settings before the first event. `BPM` must be declared before any event, even though the score class has an initial tempo of 60.

`TRANSPOSE` and `TUNINGA4` affect pitches parsed after them. Tempo and temporal parameters are associated with score states; use section boundaries for deliberate timing changes. Analysis sizes, pitch-template settings, descriptors, and model selection are global configuration, so placing them between events does not schedule a change during performance.

The reference pages follow `Score::NewConfig` in `Sources/OpenScofo/score.cpp`, parser state in `score.hpp`, and the defaults in `states.hpp`. Runtime notes distinguish accepted configuration from settings actually used by the follower. Numeric ranges describe parser validation; acceptance alone does not guarantee that a combination is useful for your instrument or model.
