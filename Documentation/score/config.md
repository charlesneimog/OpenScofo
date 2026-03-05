# Score Configuration

`OpenScofo` is configured directly through the score, which simplifies its use across different environments.

Instead of writing separate methods for each platform — such as Pd, Max, SuperCollider or others — a simple interface can be created. The composer defines all necessary settings within the score itself.

This approach allows development to focus on improving the core functionality of `OpenScofo`, rather than maintaining platform-specific integration layers.

!!! tip "Always try on `OpenScofo` Editor"
    Always try the examples on [`OpenScofo` Online Score Editor](./../../Editor), with color highlight, write scores its easier.

---
## <h2 align="center">Audio Configuration</h2>
---

### `FFTSIZE`

- `Default is 2048`
- `Must be a power of 2, I recomend 1024, 2048 or 4096`


Define the FFT Size used in decoding.

```
FFTSIZE 2048
```

---

### `HOPSIZE`

- `Default is 1024`
- `Must be a power of 2, I recomend 512 or 1024`

Define the Hop Size used in decoding.

```
HOPSIZE 4096
```

---



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

Set the value of $\eta_s$ (Sync Strength) corresponds to what [Large and Jones](https://psycnet.apa.org/doi/10.1037/0033-295X.106.1.119){:target="_blank"} (1999) refer to as the adaptation rate. This value determines how much of the previous predictions will be considered for the next BPM prediction and must be between 0 and 1. According to [Large and Jones](https://psycnet.apa.org/doi/10.1037/0033-295X.106.1.119) (1999, p. 131): 
> if it is set to 1, each estimate of $\kappa$ will be based solely on the current onsets. If $\eta_s$ < 1, the focus adapts more slowly because the previous context is taken into account.

```
SYNCSTRENGTH 0.4
```

---
## <h2 align="center">Listening Module Configuration</h2>
---

### `PITCHSIGMA`

- `Default is 0.5`
- `Range -12 - 12`

Defines the width of the pitch template: wider values increase flexibility, narrower values increase precision. This value scale using `MIDI`, so using `0.5` will give you a flexibility of half-tone more or less.

### `TIMBREMODEL`

Defines a path to a `.onnx` model trained with `py.train` for identification of extended techniques. These paths are relative to the score file, so `TIMBREMODEL "flute.onnx"` expected a file `flute.onnx` side by side of the score file loaded.

```
TIMBREMODEL "flute.onnx"
```

