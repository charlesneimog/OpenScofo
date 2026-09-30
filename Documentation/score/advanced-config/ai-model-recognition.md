---
tags:
  - Advanced Configuration
  - AI Models
---

# AI Model Recognition

Use an ONNX model for technique labels in `PTECH` and `UTECH` events. Recognition depends on matching the model's labels, input descriptors, and training configuration.

## Options

| Keyword | Default | Value | Effect |
| --- | --- | --- | --- |
| `ONNXMODEL` | No model | Model file path | Selects the ONNX technique-recognition model. Relative paths are resolved from the score file's directory. |
| `TIMBREMODEL` | No model | Model file path | Alias of `ONNXMODEL`. Prefer `ONNXMODEL` in new scores. |
| `ONNXDESCRIPTORS` | Empty list | Space-separated descriptor names | Supplies the model input descriptors in order, unless valid model metadata overrides them. |

Quote paths containing spaces. Model selection and descriptor configuration apply globally when the score is loaded.

The current score grammar accepts these descriptor names:

```text
mfcc logmel loudness rms power chroma zcr hfr centroid spread
flatness flux irregularity kurtosis harmonicity yin
```

Choose only the descriptors the model expects, in the order used for training. Vector descriptors such as MFCC and chroma contribute multiple input values, so their sizes also matter.

## Model metadata

Valid, nonempty `openscofo.descriptors` metadata takes precedence over `ONNXDESCRIPTORS`. When that metadata is absent or cannot be parsed, the loader uses the score/API list. A model needs a valid descriptor list from one of these sources.

The loader also reads `openscofo.sample_rate`, `openscofo.fft_size`, and `openscofo.hop_size`. Mismatches produce warnings; these metadata fields do not automatically reconfigure the analysis. Match `SR`, `FFTSIZE`, `HOPSIZE`, and the [descriptor settings](descriptors.md) to training.

Labels come from `openscofo.labels`, with a fallback to classifier labels. Technique names in the score must match available model labels; unknown labels cause score loading to fail.

## Example

This example assumes the model was trained at 48000 Hz with a 2048-sample window, a 512-sample hop, the listed descriptors, and labels `jet_whistle` and `pizzicato`. Replace these values with your model's actual training setup.

```openscofo
BPM 72
SR 48000
FFTSIZE 2048
HOPSIZE 512
ONNXMODEL "flute-techniques.onnx"
ONNXDESCRIPTORS mfcc logmel centroid flatness hfr flux zcr irregularity kurtosis
MFCCMELS 40
MFCCCOUNT 13

UTECH jet_whistle 1
PTECH pizzicato D4 0.5
```

See [AI Models](../../ai/index.md) and [Training](../../ai/train.md) for creating and exporting models, and [Events](../events.md) for technique-event syntax.
