# Conventions

This section defines the basic notation conventions used in `OpenScofo`, including how pitches, time values, and comments are represented in a score.

## Pitch 

On `OpenScofo` pitches are represented by `MIDI` or `PitchName`.

### MIDI

When using `MIDI`, each note corresponds to a specific MIDI number, where middle C (`C4`) is **60**.

* MIDI numbers increase by one per semitone:
  * C4 → 60
  * C#4 → 61
  * D4 → 62

* **Microtonal pitches** are supported using decimal values. For example:
  * `60.5` → a quarter-tone `C4`;
  * `61.25` → a eight-tone sharp `C#4`;

### Pitch Name 

Pitch names are organized into three components: **class**, **accidental modifiers**, and **octave**.

#### Class 

We follow the standard seven note classes of Western music: `C`, `D`, `E`, `F`, `G`, `A`, `B`. These are the natural notes, forming the diatonic scale. Each note class can be modified by accidentals to indicate pitch alterations.

#### Accidental Modifiers

| Symbol | Effect |
| ------ | --------------------------------------- |
| `#`    | Increases class by 1 (sharp)                  |
| `b`    | Decreases class by 1 (flat)                   |
| `##`   | Increases class by 2 (double sharp)           |
| `bb`   | Decreases class by 2 (double flat)            |
| `+`    | Increases class by 0.5 (quarter-tone sharp)   |
| `#+`   | Increases class by 1.5 (sharp + quarter-tone) |
| `-`    | Decreases class by 0.5 (quarter-tone flat)    |
| `b-`   | Decreases class by 1.5 (flat + quarter-tone)  |

#### Octave

The octave in a pitch name defines its position on the keyboard. `OpenScofo` follows the _Scientific Pitch Notation_ ([SPN](https://rwu.pressbooks.pub/musictheory/chapter/aspn/){:target="_black"}) convention, where `C4` represents middle C (the central C on a piano).

!!! warning "C4 notation"
    `C4` refers to middle C (the central C on a piano).

---
## Time Definition

On `OpenScofo` you represent `<DURATION>` using **number of tempos** relative to the Time Unit. Fractions like `(1/2)`, `(1/1)`, `(1/8)` are not supported.

`OpenScofo` uses time representation relative to the value of `BPM`. For example, if in your score you have the `BPM` set as **:material-music-note-quarter: = 60** this means that

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

!!! warning "Tied notes are not supported"
    `OpenScofo` does not currently support tied notes. Instead, represent the combined duration as a single value. For example, if **:material-music-note-quarter: = 60**, a dotted quarter note should be written as `1.5`.

---
## Comments

You can also create comments on `OpenScofo` language. For one line use `//` and for multiple lines use `/*` to start and `*/` to finish it.


``` openscofo
// This is a piece of Charles K. Neimog

/*
    In this piece we have a lot of different events and actios that happens 
    with the computer.
*/

BPM 80

NOTE C4 1
REST 2
NOTE D5 1
```
