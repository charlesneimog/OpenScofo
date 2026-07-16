---
icon: material/glasses
tags:
  - Advanced Configuration
---

# Advanced Configuration

Use these settings only when the default follower or listening behavior needs adjustment.

Most scores should only change `BPM`, `TRANSPOSE`, `TUNINGA4`, `ONNXMODEL`, `ONNXDESCRIPTORS`, and, for silence-sensitive pieces, `DBTHRESHOLD`. The basic score options `BPM`, `TRANSPOSE`, `TUNINGA4`, and `SR` are introduced in [Configuring a Score](config/).

Options such as `FFTSIZE`, `HOPSIZE`, MFCC, chroma, YIN, and ZCR settings are advanced. Change them only while testing a specific listening problem, for example missed low attacks, noisy extended techniques, or unstable pitch sources.

There is no universal best setting: each change trades sensitivity, stability, latency, or model compatibility against something else.

## When to use?

__Use this Advanced Setting when__ OpenScofo needs to match the written tempo, transposition, tuning reference, audio engine sample rate, or different MIR behaviors.

!!! danger "Risk if misconfigured!"
    The follower may expect events too early or too late, listen for the wrong pitch, or disagree with a model trained at another sample rate.

## Pitch following

### When to use?

__Use this Advanced Setting when__ ordinary `NOTE` events are followed too narrowly, too loosely, or incorrectly in a passage with vibrato, unstable intonation, multiphonics, or strong harmonics.

!!! danger "Risk if misconfigured!"
    Too much tolerance can make different notes look similar; too little tolerance can reject expressive intonation. YIN and rolloff changes can improve one register while making another less reliable.

### Options

`PITCHTEMPLATESIGMA`, `PITCHTEMPLATEHARMONICS`, `TUNINGA4`, `TRANSPOSE`, `YINTHRESHOLD`, `YINMINFREQUENCY`, `YINMAXFREQUENCY`, `SPECTRALROLLOFFCUTOFF`.

### Practical start

First confirm that `TUNINGA4` and `TRANSPOSE` are correct. If the performer uses wide vibrato or unstable pitch, try a slightly larger `PITCHTEMPLATESIGMA`. If the follower is attracted to harmonics or neighboring notes, return toward the default.

Touch YIN and spectral options only while testing a concrete pitch-detection failure.

```openscofo
TUNINGA4 442
PITCHTEMPLATESIGMA 0.8
PITCHTEMPLATEHARMONICS 8
```

## Tempo and synchronization

### When to use?

__Use this Advanced Setting when__ OpenScofo follows the right events but its timing feels too rigid, too nervous, or too slow to adapt to the performer.

!!! danger "Risk if misconfigured!"
    High synchronization can make the follower chase expressive rubato or false detections; low synchronization can make it feel heavy and late.

### Options

`BPM`, `PHASECOUPLING`, `SYNCSTRENGTH`.

### Practical start

Set a realistic `BPM` first. Raise `SYNCSTRENGTH` only when the performer changes tempo and the follower lags behind. Lower it when the follower overreacts to small expressive timing changes.

Adjust `PHASECOUPLING` when attacks are recognized but the internal beat phase feels consistently pulled ahead or behind.

```openscofo
BPM 96
SYNCSTRENGTH 0.4
PHASECOUPLING 0.5
```

## Silence/rest detection

### When to use?

__Use this Advanced Setting when__ rests are skipped because stage noise is interpreted as playing, or soft entrances are missed because they are treated as silence.

!!! danger "Risk if misconfigured!"
    A threshold that is too sensitive can follow room noise; a threshold that is not sensitive enough can miss quiet playing, breath sounds, or delicate attacks.

### Options

`DBTHRESHOLD`.

### Practical start

Change `DBTHRESHOLD` only after testing with the real microphone, performer, room, and electronics. If background noise triggers events, try a less sensitive threshold. If soft sounds disappear, try a more sensitive threshold.

```openscofo
DBTHRESHOLD -55

REST 2
NOTE C4 1
```

## Onset/articulation detection

### When to use?

__Use this Advanced Setting when__ attacks, articulations, repeated notes, or percussive gestures are consistently detected too late, too early, or not at all.

!!! danger "Risk if misconfigured!"
    A setting that catches sharp attacks may become jumpy with noisy or sustained sounds. Smaller analysis windows can feel faster but less stable; larger windows can be steadier but later.

### Options

`ONSETFUNCTION`, `MEDSPAN`, `FFTSIZE`, `HOPSIZE`.

### Practical start

Try `ONSETFUNCTION` before changing analysis sizes. `mkl` is the default; alternatives include `pow`, `pd`, `wpd`, `sf`, `cd`, `rcd`, and `hfc`. Use a short test passage with the exact articulation that fails.

Change `MEDSPAN`, `FFTSIZE`, or `HOPSIZE` only when deliberately testing onset sensitivity against latency and stability.

```openscofo
ONSETFUNCTION hfc

NOTE C4 0.25
NOTE C4 0.25
NOTE C4 0.5
```

## Extended technique / AI model recognition

### When to use?

__Use this Advanced Setting when__ the score contains `PTECH` or `UTECH` events and needs extended technique recognition from a trained ONNX model.

!!! danger "Risk if misconfigured!"
    If the wrong model, labels, descriptor order, or sample rate is used, extended technique recognition can become meaningless even when the score syntax is correct.

### Options

`ONNXMODEL`, `TIMBREMODEL`, `ONNXDESCRIPTORS`, `SR`.

### Practical start

Load the model that was trained for the instrument and techniques in the score. Use `ONNXMODEL` with the same descriptor list and order used during training. `TIMBREMODEL` is accepted as an alias, but `ONNXMODEL` is clearer in new scores.

Keep `SR` consistent with the model training sample rate when the model expects one.

```openscofo
ONNXMODEL flute-techniques.onnx
ONNXDESCRIPTORS mfcc logmel centroid flatness hfr flux zcr irregularity kurtosis

UTECH jet_whistle 1
PTECH pizzicato D4 0.5
```

## Timbre and descriptor tuning

### When to use?

__Use this Advanced Setting when__ an ONNX model depends on timbral descriptors and a specific descriptor is known to be causing a recognition problem.

!!! danger "Risk if misconfigured!"
    Changing descriptor settings after training can make a good model behave like the wrong model. Chroma changes can affect pitch-class balance; MFCC changes can alter timbre recognition; ZCR changes can affect noisy or percussive technique recognition.

### Options

`MFCCMELS`, `MFCCCOUNT`, `CHROMASIZE`, `CHROMACENTEROCTAVE`, `CHROMAOCTAVEWIDTH`, `ZCRCENTER`, `ZCRPAD`, `ZCRZEROPOS`, `ZCRTHRESHOLD`, `SPECTRALROLLOFFCUTOFF`, `FFTSIZE`, `HOPSIZE`.

### Practical start

Do not change these for ordinary score preparation. They shape the numbers sent to descriptor-based listening and ONNX models, so they should normally match the training setup.

Adjust them only in a controlled test where you compare recognition before and after the change.

```openscofo
ONNXMODEL flute-techniques.onnx
ONNXDESCRIPTORS mfcc zcr chroma
MFCCMELS 40
MFCCCOUNT 13
ZCRTHRESHOLD 0.0000000001
```

## Compact reference

| Option | Composer-facing meaning | When to change it | Typical risk |
| --- | --- | --- | --- |
| `SR` | Audio sample rate expected by the score or model | Only when the host or trained model requires a specific rate | Model mismatch or incorrect analysis timing |
| `FFTSIZE` | Listening window size | When testing latency versus pitch/timbre stability | Smaller can be unstable; larger can react late |
| `HOPSIZE` | How often OpenScofo updates its listening analysis | When testing responsiveness versus stability | Smaller can be jumpy; larger can miss fast detail |
| `TUNINGA4` | Reference tuning for written pitches | When the performer or ensemble uses another A4 | Notes may be heard as out of tune or wrong |
| `PHASECOUPLING` | How strongly timing locks to the internal beat phase | When events are recognized but feel phase-shifted | Follower can pull ahead or drag behind |
| `SYNCSTRENGTH` | How strongly tempo adapts to the performer | When rubato or tempo drift is followed too slowly or too nervously | Too high chases mistakes; too low feels rigid |
| `BPM` | Starting tempo and beat duration | For the notated or expected tempo | Delays and following timing drift |
| `PITCHTEMPLATESIGMA` | Pitch tolerance around written notes | For wide vibrato, unstable pitch, or overly loose pitch matching | Too wide confuses notes; too narrow rejects expression |
| `PITCHTEMPLATEHARMONICS` | Number of harmonics considered in pitch matching | When harmonic content helps or hurts note following | Harmonics may reinforce the wrong pitch |
| `TRANSPOSE` | Semitone transposition of following pitches | For transposing instruments or octave-shifted setups | The follower listens for the wrong notes |
| `ONNXMODEL` | ONNX model used for extended technique recognition | When using `PTECH` or `UTECH` | Wrong model or labels break recognition |
| `TIMBREMODEL` | Alias for loading the timbre/ONNX model | When maintaining older scores that use this name | Same risks as `ONNXMODEL` |
| `ONNXDESCRIPTORS` | Descriptor list and order sent to the model | When matching the model training setup | Wrong order makes model output unreliable |
| `ONSETFUNCTION` | Method used to detect attacks | When articulation or repeated notes are missed | May become too sensitive or miss other attacks |
| `MFCCMELS` | Mel-band resolution for MFCC/logmel timbre descriptors | Only to match or test model training settings | Timbre model input changes |
| `MFCCCOUNT` | Number of MFCC values | Only to match or test model training settings | Descriptor size or timbre balance changes |
| `MEDSPAN` | Smoothing span for onset detection | When onset detection is too jumpy or too dull | Can smear attacks or create false triggers |
| `DBTHRESHOLD` | Level below which sound is treated as silence | For noisy stages or very soft playing | Noise may trigger events, or soft notes may disappear |
| `SPECTRALROLLOFFCUTOFF` | Brightness/energy cutoff used by spectral descriptors | When testing descriptor behavior for timbre or pitch failures | Timbre balance changes unexpectedly |
| `YINTHRESHOLD` | Confidence threshold for YIN pitch tracking | Only when testing a specific pitch-tracking failure | Pitch may become too eager or too reluctant |
| `YINMINFREQUENCY` | Lowest pitch YIN will consider | For instruments outside the default low range | Low notes may be ignored |
| `YINMAXFREQUENCY` | Highest pitch YIN will consider | For instruments outside the default high range | High notes may be ignored |
| `CHROMASIZE` | Number of pitch-class bins for chroma descriptors | Only to match a model or descriptor experiment | Chroma model input changes |
| `CHROMACENTEROCTAVE` | Octave region emphasized by chroma | When testing register-sensitive chroma behavior | Wrong register gets emphasized |
| `CHROMAOCTAVEWIDTH` | Width of octave emphasis in chroma | When chroma is too narrow or too broad in register | Can overfocus or blur register information |
| `ZCRCENTER` | Whether ZCR analysis is centered in the frame | Only for ZCR/model troubleshooting | Noisy-technique descriptor changes |
| `ZCRPAD` | Whether ZCR analysis pads the frame | Only for ZCR/model troubleshooting | ZCR values may no longer match training |
| `ZCRZEROPOS` | How zero crossings are counted around zero | Only for ZCR/model troubleshooting | Small waveform changes alter ZCR |
| `ZCRTHRESHOLD` | Minimum crossing level for ZCR | When noisy or quiet techniques confuse ZCR | Too low follows noise; too high loses detail |
