---
icon: material/robot
tags:
  - AI Models
---

# AI Models

Use AI models to recognize non pitched sounds, this include: 

- breath; 
- key clicks, tongue-ram, jet-whistle;
- percussive sounds;
- phonems;


I would make the workflow more explicit and avoid separating “How to train” from “Train and Use”:

## Prepare Training Data

To train a model for your piece, organize audio recordings of each extended technique into separate folders. Each folder represents one class that the model will learn to recognize.

For example, if a flute piece uses `pizz`, `jet-whistle`, and `tongue-ram`, the training data could be organized as:

```text
Flute/
├── pizz/
│   ├── pizz_C4_recording_01.wav
│   ├── pizz_C4_recording_02.wav
│   ├── pizz_D4_recording_01.wav
│   └── pizz_D4_recording_02.wav
├── jet-whistle/
│   ├── jet_recording_01.wav
│   └── jet_recording_02.wav
└── tongue-ram/
    ├── tongue-ram_C3_recording_01.wav
    ├── tongue-ram_C3_recording_02.wav
    ├── tongue-ram_D3_recording_01.wav
    └── tongue-ram_D3_recording_02.wav
```

The folder names become the **class labels** of the trained model. In this example, the model learns to distinguish between `pizz`, `jet-whistle`, and `tongue-ram`.

!!! tip
    Record several examples of each technique. Variations in pitch, dynamics, articulation, and recording conditions can help the model recognize the technique more reliably.

## Train the Model

Once you have the `Flute` folder with all class-label folders inside, compress it into a `.zip` file. Then open the [OpenScofo Model Trainer](https://openscofo-production.up.railway.app/){target="_blank"} and upload the ZIP file.

The trainer will:

1. After upload, click on `Train Model`;
2. Wait some moments;
3. Download the model;
4. Rename it;
5. Load it in your score;

For example, if you rename to `flute.onnx`, you must use:

``` openscofo
ONNXMODEL flute.onnx
```

The resulting model can then be loaded directly into an OpenScofo score.

## Use the Model in a Score

Load the exported model with `ONNXMODEL` and use the trained class labels in your technique events:

```openscofo
ONNXMODEL flute.onnx

// Pitch is still written in the score, but the model is used
// to recognize the timbral/technical class.
PTECH pizz A4 2
PTECH tongue-ram C3 2
UTECH jet-whistle 2
```

When the model was trained using the OpenScofo training tools, information about the audio descriptors is stored in the ONNX model metadata, so `ONNXDESCRIPTORS` does not need to be specified manually.

See also: [Audio Analysis](../descriptors/), [Advanced Configuration](../score/advanced-config/), and [AI Model Reference](../descriptors/ai/).

