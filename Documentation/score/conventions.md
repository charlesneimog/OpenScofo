---
icon: material/piano
tags:
  - Notation
---

# Pitch and Notation Conventions

Use this page as a lookup for pitch names, MIDI pitches, durations, and comments.

## Reference Table

### Pitch

| Form | Example | Meaning |
| --- | --- | --- |
| Pitch name | `C4`, `F#4`, `Bb3` | Scientific pitch notation; `C4` is middle C. |
| Quarter-tone | `C+4`, `D-4` | Quarter-tone sharp / flat. |
| Compound accidental | `#+`, `b-`, `##`, `bb` | Sharp+quarter, flat+quarter, double sharp, double flat. |
| MIDI | `60`, `61`, `60.5` | `60` is `C4`; decimals allow microtones. |

### Durations

Durations are beats relative to the current `BPM`.

| Value | Meaning when :material-music-note-quarter: = 100 |
| --- | --- |
| `2` | half note (:material-music-note-half:) |
| `1.5` | Dotted quarter (:material-music-note-quarter-dotted:) |
| `1` | quarter note (:material-music-note-quarter:) |
| `0.75` | dotted eighth (:material-music-note-eighth-dotted:) |
| `0.5` | eighth note (:material-music-note-eighth:) |
| `0.25` | sixteenth note (:material-music-note-quarter:) |

Be carefull with Compound time signature (`6/8`, `9/8`, etc...). If you use **:material-music-note-quarter-dotted: = 80**, the `DURATION` for :material-music-note-eighth: will be **0.333**.

### Comments

```openscofo
// one line

/*
multiple
lines
*/
```

!!! warning "Fractions not supported"
    Fractions such as `(1/2)` are not supported. Write tied notes as one combined duration.
