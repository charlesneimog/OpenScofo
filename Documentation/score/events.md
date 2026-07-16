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
| `REST` | Silence | `REST <DURATION>` | `REST 1` | Keeps score time moving. |
| `PTECH` | Pitched technique | `PTECH <LABEL> <PITCH> <DURATION>` | `PTECH pizz C4 1` | For extended techniques **with** pitch. Check [AI](../ai/index.md)! |
| `UTECH` | Unpitched technique | `UTECH <LABEL> <DURATION>` | `UTECH jet-whistle 2` | For extended techniques **without** pitch. Check [AI](../ai/index.md)! |

## Example

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
