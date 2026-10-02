---
icon: material/music-note
tags:
  - Musical Events
---

# Musical Events

Musical Events define what OpenScofo listens for. Write one event per line, and indent any associated actions below it using one tab (or four spaces).

```openscofo
NOTE C4 1
    sendto delay [1000]
```

## Event Syntax Reference

| Event | Purpose | Syntax | Example | Remarks |
| --- | --- | --- | --- | --- |
| `NOTE` | Single pitch | `NOTE <PITCH> <DURATION>` | `#!openscofo NOTE C4 1` | Pitch name or MIDI number. |
| `CHORD` | Simultaneous pitches | `CHORD (<PITCH...>) <DURATION>` | `#!openscofo CHORD (C4 E4 G4) 2` | For chords and stable multiphonics. |
| `TRILL` | Alternating pitches | `TRILL (<PITCH...>) <DURATION>` | `#!openscofo TRILL (D4 E4) 4` | For trills and tremolos. |
| `GLISS` | Ordered pitch trajectory | `GLISS (<PITCH...>) <DURATION>` | `#!openscofo GLISS (C4 C#4 D4 D#4 E4) 2` | Glissando in quarter-tone steps between written pitches. |
| `REST` | Silence | `REST <DURATION>` | `#!openscofo REST 1` | Keeps score time moving. |
| `PTECH` | Pitched technique | `PTECH <LABEL> <PITCH> <DURATION>` | `#!openscofo PTECH pizz C4 1` | For extended techniques **with** pitch. Check [AI](../ai/index.md)! |
| `UTECH` | Unpitched technique | `UTECH <LABEL> <DURATION>` | `#!openscofo UTECH jet-whistle 2` | For extended techniques **without** pitch. Check [AI](../ai/index.md)! |

## Example

`TRILL`, `PTECH`, and `UTECH` use the strongest current internal observation, without an ordering constraint.
`UTECH` uses alternative technique labels. `PTECH` uses alternative technique labels or its expected pitch.
`GLISS` expands each interval into an ordered chain of quarter-tone (0.5 MIDI) steps during score parsing.
For example, `GLISS (C4 G4) 2` creates 15 pitch microstates from MIDI 60 through 67, including both endpoints.
Descending intervals use descending steps; shared endpoints and consecutive repeated pitches appear once.
These microstates share the event's overall duration equally
by default, and its final internal state is absorbing until the outer duration model exits the event.

`NOTE`, `CHORD`, `PTECH`, and `UTECH` do not use silence as an internal observation. After these events,
OpenScofo inserts an optional Markov silence state before the next sounded event in the same section.
The follower can go directly to that event for legato playing, or enter silence and wait for sound to return.
No extra state is inserted before an explicit `REST`, at a section boundary, or after the final event.

These unscored silence states have geometric occupancy, following the distinction between scored
semi-Markov events and atemporal Markov events in [Cont's hybrid model](https://doi.org/10.1109/TPAMI.2009.106).
OpenScofo currently uses equal probabilities for taking/bypassing the silence and for staying/leaving it;
these numerical priors are implementation defaults, not values prescribed by Cont.
A gap adds no scored beats, retains the preceding event's score position, has no score actions, and does
not trigger a tempo update. Explicit `REST` events retain their written durations and score positions.

State lists include these internal states, so state indices can differ from score positions.
Python and Lua expose `inter_event_silence` to identify them. Audio-state listeners can report `silence`
with the preceding event's score position while the follower occupies a gap.

```openscofo hl_lines="3 6 9"
BPM 60

NOTE C4 1
    sendto delay [1]

CHORD (E4 G4 B4) 2
    sendto reverb [0.8]

PTECH pizz C4 1
    sendto sample [start pizz_echo]
```

!!! warning "`TIMEDEVENT` and `LUAEVENT`"
    `TIMEDEVENT` and `LUAEVENT` are planned but not implemented. 


## Pitch and Notation Conventions

Use this page as a lookup for pitch names, MIDI pitches, durations, and comments.

### Pitch

| Form | Example | Meaning |
| --- | --- | --- |
| Pitch name | `C4`, `F#4`, `Bb3` | Scientific pitch notation; `C4` is middle C. |
| Quarter-tone | `C+4`, `D-4` | Microtonal/Quarter-tone sharp / flat. |
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
