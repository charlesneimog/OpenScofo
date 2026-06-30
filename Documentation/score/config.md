# Score Configuration

`OpenScofo` is configured directly through the score, which simplifies its use across different environments.

Instead of writing separate methods for each platform — such as Pd, Max, SuperCollider or others — a simple interface can be created. The composer defines all necessary settings within the score itself.

This approach allows development to focus on improving the core functionality of `OpenScofo`, rather than maintaining platform-specific integration layers.

---
## Audio Configuration
---

### `SR`

- `Default: same as the host`

Defines the sample rate expected by the score. This setting is optional, but it is useful when an AI model has been trained for a specific sample rate (for example, `48000` Hz). In that case, you can specify:

```text
SR 48000
```
If OpenScofo is running in a host with a different sample rate (for example, `44100` Hz), a warning will be showed.

About AI models check [How to use AI models?](../descriptors/ai.md).

### `FFTSIZE`

- `Default is 2048`
- `Must be a power of 2.`


Define the FFT Size used in decoding.

```
FFTSIZE 2048
```

---

### `HOPSIZE`

- `Default is 512`
- `Must be a power of 2.`

Define the Hop Size used in decoding.

```
HOPSIZE 256
```

---

## Time Configuration
---

### `BPM`

- `Default is 60`

`BPM` keyword set the `BPM` value for the next event.

```
BPM 50
```

---
### `PHASECOUPLING`

- `Default is 0.5`
- `Range 0-2`

This value must be between 0 and 2 and captures the amount of force exerted on attentional rhythm, determining, among other factors, the speed at which the coupled system relaxes in relation to the attractor (Large and Jones, 1999, p. 128).

---
### `SYNCSTRENGTH`

- `Default is 0.5`
- `Range 0-1`

Set the value of $\eta_s$ (Sync Strength) corresponds to what [Large and Jones](https://psycnet.apa.org/doi/10.1037/0033-295X.106.1.119){:target="_blank"} (1999) refer to as the adaptation rate. *This value determines how much of the previous predictions will be considered for the next BPM prediction and must be between 0 and 1*. According to [Large and Jones](https://psycnet.apa.org/doi/10.1037/0033-295X.106.1.119) (1999, p. 131): 
> if it is set to 1, each estimate of $\kappa$ will be based solely on the current onsets. If $\eta_s$ < 1, the focus adapts more slowly because the previous context is taken into account.

```
SYNCSTRENGTH 0.4
```

---
## Listening Module Configuration
---

### `PITCHTEMPLATESIGMA`

- `Default is 0.5`
- `Range -12 - 12`

Defines the width of the pitch template: wider values increase flexibility, narrower values increase precision. This value scale using `MIDI`, so using `0.5` will give you a flexibility of half-tone more or less. In anothers words, even if the player play a half-tone low or upper the pitch will yet 'match', which means have a good observation probability.

```
PITCHTEMPLATESIGMA 0.8
```

---
### `ONNXMODEL`

Defines a path to a `.onnx` model trained with `py.o.train` for identification of extended techniques. These paths are relative to the score file, so `ONNXMODEL "flute.onnx"` expected a file `flute.onnx` side by side of the score file loaded.

```
ONNXMODEL "flute.onnx"
```

!!! tip "Train models the object `py.o.train`."
    The entire model of `OpenScofo` is designed to be used with models trained by the Pd Object `py.o.train`. More details about AI models check [How to use AI models?](../descriptors/ai.md)


---
### `ONNXDESCRIPTORS`

Defines which descriptors did you use for the train of the `ONNX` model.

```

ONNXDESCRIPTORS mfcc zcr centroid spread
```

!!! warning "Order matters!"
    Use descriptors in the same order as your training. For example, using `zcr mfcc` when the model was trained with `mfcc zcr` will **change the result**. More details about AI models check [How to use AI models?](../descriptors/ai.md).


To detect Flute extended techniques, like tongue-ram, pizz, key-click, jet-whistle, etc... A good list would be `mfcc logmel zcr centroid flatness hfr`.


Check the `OpenScofo` [descriptors](./../descriptors/spectral.md).

---
### `ONSETFUNCTION`

Defines the onset detection function (ODF), each emphasizing different signal characteristics.

* `pow` detects onsets through changes in signal energy and works well for signals with strong amplitude attacks, such as percussion.
* `pd` measures phase deviation between frames, detecting irregularities in phase progression that often occur at note attacks.
* `wpd` is a weighted version of phase deviation where phase changes are weighted by spectral magnitude, emphasizing stronger spectral components.
* `sf` (spectral flux) measures the positive change in the magnitude spectrum between consecutive frames and is widely used for detecting new spectral energy introduced at onsets.
* `cd` (complex domain) combines magnitude and phase information to detect deviations from the expected spectral evolution, making it robust across different musical signals.
* `rcd` is a rectified version of the complex-domain method that counts only increases in deviation, helping reduce false detections.
* `hfc` (high frequency content) emphasizes changes in high-frequency spectral bins, which are often strong during percussive attacks.
* `mkl` (modified Kullback–Leibler divergence) measures changes in the spectral distribution between frames and is effective for detecting structural spectral changes in pitched or harmonic material.

!!! warning "Change this only if you know what you are doing"
    This parameter is a core component of the extended-techniques model. It determines when the listener model should rely on the pitch model and when it should switch to the extended-techniques model. Changing it can significantly affect score-following performance.
