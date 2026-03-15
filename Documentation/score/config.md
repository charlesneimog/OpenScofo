# Configuration

`OpenScofo` is configured directly through the score, which simplifies its use across different environments.

Instead of writing separate methods for each platform — such as Pd, Max, SuperCollider or others — a simple interface can be created. The composer defines all necessary settings within the score itself.

This approach allows development to focus on improving the core functionality of `OpenScofo`, rather than maintaining platform-specific integration layers.

!!! tip "Always try on `OpenScofo` Editor"
    Always try the examples on [OpenScofo Online Score Editor](https://charlesneimog.github.io/OpenScofo/Editor){:target="_blank"}, with color highlight. Writing scores is easier there.

---
## <h2 align="center">Audio Configuration</h2>
---

### `FFTSIZE`

- `Default is 2048`
- `Must be a power of 2.`


Define the FFT Size used in decoding.

```
FFTSIZE 2048
```

---

### `HOPSIZE`

- `Default is 1024`
- `Must be a power of 2.`

Define the Hop Size used in decoding.

```
HOPSIZE 256
```

---

## <h2 align="center">Time Configuration</h2>
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
## <h2 align="center">Listening Module Configuration</h2>
---

### `PITCHTEMPLATESIGMA`

- `Default is 0.5`
- `Range -12 - 12`

Defines the width of the pitch template: wider values increase flexibility, narrower values increase precision. This value scale using `MIDI`, so using `0.5` will give you a flexibility of half-tone more or less. In anothers words, even if the player play a half-tone low or upper the pitch will yet 'match', which means have a good observation probability.

```
PITCHTEMPLATESIGMA 0.8
```

---
### `TIMBREMODEL`

Defines a path to a `.onnx` model trained with `py.train-onnx` for identification of extended techniques. These paths are relative to the score file, so `TIMBREMODEL "flute.onnx"` expected a file `flute.onnx` side by side of the score file loaded.

```
TIMBREMODEL "flute.onnx"
```

!!! tip "Train models the object `py.train-onnx`."
    The entire model of `OpenScofo` is designed to be used with models trained by the Pd Object `py.train-onnx`. Check it on Resources

---
### `ONNXDESCRIPTORS`

Defines which descriptors did you use for the train of the `ONNX` model.

```

ONNXDESCRIPTORS mfcc zcr centroid spread
```

!!! warning "Order matters!"
    Use descriptors in the same order as your training. For example, using `zcr mfcc` when the model was trained with `mfcc zcr` will **change the result**.


Check the `OpenScofo` [descriptors](../listening.md).

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




??? question "More details for Onset Detection"






