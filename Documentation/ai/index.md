---
icon: material/robot
tags:
  - AI Models
---

# AI Models

Use AI models to recognize **unpitched or timbre-based sounds and playing techniques**, including:

- breath sounds;
- key clicks, tongue rams, and jet whistles;
- percussive sounds;
- phonemes.

OpenScofo represents these sounds using `PTECH` and `UTECH` events. See [Score Events](../score/events/#event-syntax-reference) for more details.

## Prepare Training Data

To train a model for your piece, organize audio recordings of each extended technique into separate folders. Each folder represents one **class** that the model will learn to recognize.

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

!!! tip "More recordings better recognition!"
    Record several examples of each technique. Variations in pitch, dynamics, articulation, and recording conditions can help the model recognize each technique more reliably.

## Train the Model

Once the `Flute` folder contains all the class-label folders, compress it into a `.zip` file. Then open the [OpenScofo Model Trainer](https://openscofo-production.up.railway.app/){target="_blank"}.

Then:

1. Upload the `.zip` file;
2. Click **Train Model**;
3. Wait for the training to finish;
4. Download the generated `.onnx` model;
5. Give the model a descriptive name, such as `flute.onnx`;
6. Place it alongside your score and load it with `ONNXMODEL`.

For example, if you rename the model to `flute.onnx`, use:

```openscofo
ONNXMODEL flute.onnx
```

!!! warning "**Use exactly the same class labels that you used as folder names during training.**"

    For example, if the training folder is named `pizz`, then `pizz` is the class label stored in the model. You must therefore write:

    ```openscofo
    PTECH pizz C4 2
    ```

    Using a different label, such as:

    ```openscofo
    PTECH pizzicato C4 2
    ```

    will not work because `pizzicato` was not one of the class labels used during training.

## Use the Model in a Score

Load the exported model with `ONNXMODEL` and use the trained class labels in your technique events:

```openscofo
ONNXMODEL flute.onnx

// Pitch is still written in the score, while the model is used
// to recognize the playing-technique class.
PTECH pizz A4 2
PTECH tongue-ram C3 2
UTECH jet-whistle 2
```

When a model is trained using the OpenScofo training tools, information about the audio descriptors is stored in the ONNX model metadata. Therefore, `ONNXDESCRIPTORS` does not need to be specified manually.

See also: [Audio Analysis](../descriptors/), [Advanced Configuration](../score/advanced-config/), and [AI Model Reference](../descriptors/ai/).
