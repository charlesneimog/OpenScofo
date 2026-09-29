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

## Reference Table

| Event | Purpose | Syntax | Example | Remarks |
| --- | --- | --- | --- | --- |
| `NOTE` | Single pitch | `NOTE <PITCH> <DURATION>` | `NOTE C4 1` | Pitch name or MIDI number. |
| `CHORD` | Simultaneous pitches | `CHORD (<PITCH...>) <DURATION>` | `CHORD (C4 E4 G4) 2` | For chords and stable multiphonics. |
| `TRILL` | Alternating pitches | `TRILL (<PITCH...>) <DURATION>` | `TRILL (D4 E4) 4` | For trills and tremolos. |
| `GLISS` | Ordered pitch trajectory | `GLISS (<PITCH...>) <DURATION>` | `GLISS (C4 C#4 D4 D#4 E4) 2` | Glissando in quarter-tone steps between written pitches. |
| `REST` | Silence | `REST <DURATION>` | `REST 1` | Keeps score time moving. |
| `PTECH` | Pitched technique | `PTECH <LABEL> <PITCH> <DURATION>` | `PTECH pizz C4 1` | For extended techniques **with** pitch. Check [AI](../ai/index.md)! |
| `UTECH` | Unpitched technique | `UTECH <LABEL> <DURATION>` | `UTECH jet-whistle 2` | For extended techniques **without** pitch. Check [AI](../ai/index.md)! |

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
