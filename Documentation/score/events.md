# Events

`OpenScofo` currently supports six event types: `NOTE`, `TRILL`, `CHORD`, `REST`, `PTECH`, and `UTECH`. Two additional event types, `LUAEVENT` and `TIMEDEVENT`, are planned but have not yet been implemented.

!!! tip "Check the MusicXML Importer"
    `OpenScofo` has a [MusicXML importer](https://charlesneimog.github.io/OpenScofo/Editor/), it is very accurate and tested on MuseScore (but should work for all major programs). 

---
## `NOTE`

`NOTE` events describes tradicional pitches. It must be defined as `NOTE <PITCH> <DURATION>`.

<score id="note_event" timeSig="4/4" notes="C5/q, F#4/8, G4/8, A5/h"></score>

``` openscofo
NOTE C5 1
NOTE F#4 0.5
NOTE G4 0.5
NOTE A5 2
```

---
## `CHORD`

`CHORD` events describes chords, stable multiphonics, and others events. It must be defined as `CHORD (<PITCH1> <PITCH2>) <DURATION>`. 

<score id="chord_event" timeSig="3/4" notes="(C#5 E5 G#5)/h, (C4 D#4 Ab5)/q"></score>

``` openscofo
CHORD (C#5 E5 G#5) 2
CHORD (C4 Ab5 D4) 1
```

---
## `TRILL`

`TRILL` events describes trill and tremolo events. It must be defined as `TRILL (<PITCH1> <PITCH2>) <DURATION>`.

<score id="tremolo_test" timeSig="4/4" notes="G4/h:trem3, D#4/h"></score>

``` openscofo
TRILL (G4 D#4) 4
```

---
## `REST`

`REST` events describes rests. It must be defined as `REST <DURATION>`.

<score id="rest_event" timeSig="4/4" notes="B4/hr, G4/q, D#5/8, C5/16r, C5/16"></score>

``` openscofo
REST 2
NOTE G4 1
NOTE D#5 0.5
REST 0.25
NOTE C5 0.25
```

---
## `PTECH` 


`PTECH` describe events that are non tradicional events but yet have pitches on it, for example: `tongue-ram`, `key-click`. Events must be defined as:

```text
PTECH <LABEL> <PITCH> <DURATION>
```

- `<LABEL>`: Name of the technique, such as `pizz` or `tongue-ram`. This label must match one of the labels used during model training. See [How to use AI models?](../descriptors/ai.md).
- `<PITCH>`: Pitch of the note, such as `C4` or `D#5`. Required for pitched techniques.
- `<DURATION>`: Duration in beats, written as an integer or floating-point value.

!!! tip "`PTECH` events require an ONNX model"
    `PTECH` require an ONNX model capable of recognizing the corresponding playing techniques. See [How to use AI models?](../descriptors/ai.md).


<score id="ptech_example" timeSig="4/4" notes="C4/h:x, D4/h:tongue-ram"></score>

```openscofo
PTECH pizz C4 2
PTECH tongue-ram D4 2
```

---
## `UTECH` 

`UTECH` describe events that are non tradicional events and do not have pitches on it, for example: `jet-whistle`, `sounds-behind-the-bridge`. Events must be defined as:

<score id="utech_example" timeSig="4/4" notes="C4/h:x, C4/h:tongue-ram"></score>

`UTECH` events must be defined as:

```
UTECH <LABEL> <DURATION>
```

* `<LABEL>`: Name of the technique (e.g., `slap`, `jet-whistle`).
* `<DURATION>`: Duration in beats (integer or float).

UTECH events are ideal for percussive or unpitched extended techniques. 

```openscofo
UTECH slap 2
UTECH jet-whistle 1
```

!!! tip "`UTECH` events require an ONNX model"
    `UTECH` require an ONNX model capable of recognizing the corresponding playing techniques. See [How to use AI models?](../descriptors/ai.md).

---
## `TIMEDEVENT`

!!! danger "Not implemented yet"
    Planned to next research steps.

This event is adapted from patches by Cort Lippe. In a Max patch, a simple score is displayed for the player, with a slider beneath it that fills progressively during the event. This feature is implemented using TIMEEVENT `<TIME>`, where `<TIME>` is specified in seconds.


---
## `LUAEVENT`

!!! danger "Not implemented yet"
    Planned to next research steps.

In this event, you will be able to define a Lua function to be called on each evaluation of the model. The syntax will be something like `<LUAEVENT> <FUNCTION_TO_ENTER> <FUNCTION_TO_EXIT>`. The `FUNCTION_TO_ENTER` you set a Lua function that receive a Description as argument. This can help the model to enter in the event (note that there is a `TIME` model running too), lets say you will enter when, after sometime, the is a `C4` note. And you will exit when there is a `G4` for at least 2 seconds.
