# AI Model

`OpenScofo` uses a machine learning model based on a **Random Forest classifier** to detect extended techniques.

The model operates on **spectral descriptors** extracted from audio recordings (see previous section). To enable detection in real time, the model must first be trained on labeled examples of the techniques you want to recognize in real time.

---

## Training Overview

Training consists of three steps:

1. **Collect audio samples**
2. **Organize them into labeled folders**
3. **Extract features and train the model**

---

### 1. Dataset Structure

Audio files must be organized in a directory structure where:

* The **top-level folder** (any good name for it).
* Each **subfolder represents a technique label**
* Each subfolder contains `.wav` or `.aif` audio examples of that technique, small samples, not a record of 30 minutes full of tongue-ram.

Example structure

```python
Flute/
├── jet_whistle/
│   ├── jet-whistle_01.wav
│   ├── jet-whistle_02.wav
│   ├── jet-whistle_03.wav
│   └── ...
├── key_click/
│   ├── Fl-key_click_A#4.wav
│   ├── Fl-key_click_A4.wav
│   ├── Fl-key_click_F#4.wav
│   ├── Fl-key_click_F4.wav
│   ├── Fl-key_click_G#4.wav
│   └── ...
├── pizzicato/
│   ├── Fl-pizzicato_A#4.wav
│   ├── Fl-pizzicato_A4.wav
│   ├── Fl-pizzicato_B3.wav
│   └── ...
└── tongue_ram/
    ├── Fl-tongue_ram_A3.wav
    ├── Fl-tongue_ram_B3.wav
    ├── Fl-tongue_ram_C#3.wav
    └── ...
```

---

### 2. Feature Extraction

After preparing the dataset, select the spectral descriptors used for training.

Commonly used feature set in my pieces:

* MFCC
* Log-mel spectrogram features
* Spectral centroid
* Spectral flatness
* High-frequency ratio (HFR)
* Spectral flux
* Zero-crossing rate (ZCR)
* Irregularity

You may adjust this set depending on the target instrument and recording conditions, but consistency between training and inference is required.

---

### 3. Training Procedure

Once the dataset and feature set are defined:

1. Load all audio files from the dataset structure
2. Extract the selected spectral descriptors
3. Train a **Random Forest classifier**
4. Save the trained model for inference in `OpenScofo`

---

### Practical Note

* The quality of classification depends more on **dataset quality and consistency** than on model complexity.
* Balanced representation across techniques is strongly recommended to avoid bias. For example, 80 samples of `tongue-ram` and just one `jet-whistle` is very bad.


## Training 

For now you can use Pure Data or Python to train these models.

### Pure Data

For Pure Data you can use the `py4pd` object + `py.o.train`. With these objects, you can easily train your model. To install it, you need to follow the steps:

1. Install last version of Python. You can check this link to install: [https://www.python.org/downloads/](https://www.python.org/downloads/).
2. Open Pure Data, Go to Tools, Find External, search for `py4pd` and install it;
3. Them add `declare -lib py4pd`  and create the object `py.o.train`.
4. Open the help patch of `py.o.train`. Read and explore it.

<img src="../assets/py.o.train.png" style="display:block; margin: 0 auto; width:80%;">


### Python

On Python, the `OpenScofo` module provide a class to allow this training. You can use it with the following script.

``` python
import OpenScofo

# sample_rate, fft_size and hop_size must be the same you will use in the score
trainer = OpenScofo.ExtendedTechniqueClassifier(
    sample_rate=48000,
    fft_size=2048,
    hop_size=512,
)

# check the ids using in spectral descriptors
# you must set then in the same order in the score:
#       ONNXDESCRIPTORS mfcc logmel centroid flatness hfr flux zcr irregularity kurtosis
trainer.set_descriptors(
    [
        "mfcc",
        "logmel",
        "centroid",
        "flatness",
        "hfr",
        "flux",
        "zcr",
        "irregularity",
        "kurtosis",
    ]
)

# folder where there is sub-folder with audios inside, for example:
'''
Flute
├── jet_whistle
│   ├── jet-whistle1.wav
│   ├── jet-whistle2.wav
│   └── jet-whistle-3.wav
│   └── ...
├── key_click
│   ├── Fl-key_cl-A#4-f-N-N.wav
│   ├── Fl-key_cl-A4-f-N-N.wav
│   ├── Fl-key_cl-F#4-f-N-N.wav
│   ├── Fl-key_cl-F4-f-N-N.wav
│   ├── Fl-key_cl-G#4-f-N-N.wav
│   └── ...
├── pizzicato
│   ├── Fl-pizz-A#4-f-N-N.wav
│   ├── Fl-pizz-A4-f-N-N.wav
│   ├── Fl-pizz-B3-f-N-N.wav
│   ├── Fl-pizz-B3-f-N-N.wav
│   └── ...
└── tongue_ram
    ├── Fl-tng_ram-A3-mf-N-N.wav
    ├── Fl-tng_ram-B3-mf-N-N.wav
    ├── Fl-tng_ram-C#3-mf-N-N.wav
│   └── ...
'''

trainer.set_train_folder("/home/neimog/Downloads/Flute")

# Impulse Responses are good to prevent overfit.
trainer.set_ir_folders(["/home/neimog/Nextcloud/MusicData/Impulse_Responses/05_Halls/"])

trainer.analyze()
trainer.train()

# name of the model, you import it on OpenScofo using 
# ONNXMODEL flute-v5.onnx

trainer.export_model("flute-v5.onnx")
```

## Score Example

Here one score example using a trained model.

<img src="../assets/score-extended-techniques.png" style="display:block; margin: 0 auto; width:80%;">

This will be the score for OpenScofo:

```
/* Generated by OpenScofo online editor */

BPM 80

// Model Exported
ONNXMODEL flute-v5.onnx

// Descriptors used for train (exact same order)
ONNXDESCRIPTORS mfcc logmel centroid flatness hfr flux zcr irregularity kurtosis

// Measure number 1
UTECH jet_whistle 1
REST 0.5	
NOTE Bb4 1.5 // tied
NOTE A4 0.5
REST 0.5
	

// Measure number 2
UTECH jet_whistle 1
REST 0.5
NOTE Bb4 1.5 // tied
NOTE A4 0.5
REST 0.5

// Measure number 3
PTECH pizzicato D4 0.5
PTECH pizzicato A4 0.5
REST 0.5
PTECH pizzicato D4 0.5
PTECH pizzicato A4 0.5
REST 0.5
PTECH pizzicato D4 0.5
PTECH pizzicato A4 0.5

// Measure number 4
REST 0.5
PTECH pizzicato D4 0.5
PTECH pizzicato A4 0.5
REST 0.5
UTECH jet_whistle 1
REST 1
	
	

// Measure number 5
PTECH pizzicato D4 0.5
PTECH pizzicato A4 0.5
REST 0.5
PTECH pizzicato D4 0.5
PTECH pizzicato Bb4 0.5
REST 0.5
PTECH pizzicato D4 0.5
PTECH pizzicato B4 0.5

// Measure number 6
REST 0.5
PTECH pizzicato D4 0.5
PTECH pizzicato Bb4 0.5
REST 0.5
UTECH jet_whistle 1
UTECH jet_whistle 1

// Measure number 7
UTECH jet_whistle 1
REST 1
REST 2

// Measure number 8
NOTE D4 2
REST 2
``` 
