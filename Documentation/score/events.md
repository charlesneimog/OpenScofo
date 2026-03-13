# Events

`OpenScofo` for now supports four events, `NOTE`, `TRILL`, `CHORD` and `REST`. Before define then, I need to explain how to represent `PITCHES` and `DURATIONS`.

!!! tip "Check the MusicXML Importer"
    `OpenScofo` has a [MusicXML importer](../../Editor), it is very accurate and tested on MuseScore and Sibelius (but should work for all major programs). 


--- 
## Pitches and Durations

### `PITCH`

On OpenScofo pitches are represented by `name` or `MIDI`.

- `NOTENAME`: Things like `C4`, `Db5`, `C#5`, `Gb5`, etc...
- `MIDI`: Normal `MIDI`. 60 for `C4`, 67 for `G4`, etc...

---
### `TIME`

On `OpenScofo` you represent `<DURATION>` using **number of tempos** relative to the Time Unit. Fractions like `(1/2)`, `(1/1)`, `(1/8)` are not supported.

Because that `OpenScofo` uses time representation relative to the value of `BPM` define previously. For example, if in your score you have the `BPM` set as **:material-music-note-quarter: = 60** this means that

- :material-music-note-half: duration value is `2`;
- :material-music-note-quarter: duration value is `1`;
- :material-music-note-eighth: is `0.5`,
- :material-music-note-sixteenth: is `0.25`
- :material-music-note-eighth-dotted: is `0.75`

if in your score you have the `BPM` set as **:material-music-note-eighth: = 60** this means that

- :material-music-note-half: duration value is `4`;
- :material-music-note-quarter: duration value is `2`;
- :material-music-note-eighth: is `1`,
- :material-music-note-sixteenth: is `0.5`

---
## Events
---

### `NOTE`

`NOTE` events describes tradicional pitches. It must be defined as `NOTE <PITCH> <DURATION>`.

- `NOTE C4 2`
- `NOTE 60 2`
- `NOTE C#5 0.3`
- `NOTE Bb3 0.25`

---
### `TRILL`

`TRILL` events describes trill and tremolo events. It must be defined as `TRILL (<PITCH1> <PITCH2>) <DURATION>`.

- `TRILL (C4 D4) 2`
- `TRILL (60 67) 2`
- `TRILL (C#5 E5) 0.3`
- `TRILL (Bb4 D5) 0.25`

On [Cânticos de Silício I](https://charlesneimog.github.io/Canticos-de-Silicio-I/) I have this:

<figure markdown="span">
  ![Image title](assets/canticos_example-1.png){ width="300" }
  <figcaption>Image caption</figcaption>
</figure>

---
### `CHORD`

`CHORD` events describes chord, stable multiphonics, and others events. It must be defined as `CHORD (<PITCH1> <PITCH2>) <DURATION>`.

- `CHORD (C4 E4 G4) 2`
- `CHORD (60 64 67) 2`
- `CHORD (C#5 E5 Ab5) 0.3`
- `CHORD (Bb4 D5 F4) 0.25`

---
### `REST`

`REST` events describes rests. It must be defined as `REST <DURATION>`.

- `REST 0.2`
- `REST 4`


### `TECH`

`TECH` events describes extended techniques events. It must be defined as `TECH <LABEL> <DURATION>`. These events require that you train a ONNX model. It is a simple process describe in [Training a ONNX model](onnx.md). Basically you put audio from specific techniques inside folders with labels (for example **tongue-ram**) and execute a Pd patch, this will generate a `.onnx` file that you load using `TIMBREMODEL` in the score.

- `TECH tongue-ram C4 1`
- `TECH jet-whiste 1`


---
