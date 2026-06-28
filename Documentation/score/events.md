# Events

`OpenScofo` for now supports four events, `NOTE`, `TRILL`, `CHORD`, `REST`, `TECH`, `LUAEVENT`.

!!! tip "Check the MusicXML Importer"
    `OpenScofo` has a [MusicXML importer](https://charlesneimog.github.io/OpenScofo/Editor/), it is very accurate and tested on MuseScore and Sibelius (but should work for all major programs). 


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

!!! warning "Experimental"
    Part of ongoing research at the University of São Paulo. These events represent extended techniques **with pitch**, like `tongue-ran`, `key-click`, etc...

`PTECH` events must be defined as:

```
PTECH <LABEL> <PITCH> <DURATION>
```

* `<LABEL>`: Name of the technique (e.g., `pizz`, `tongue-ram`).
* `<PITCH>`: The pitch of the note (e.g., `C4`, `D#5`). Required for pitched techniques.
* `<DURATION>`: Duration in beats (integer or float).

These events require an ONNX model to recognize the technique. Training the model involves placing audio examples in labeled folders and running the `py.train-onnx` object, generating a `.onnx` file to load using `TIMBREMODEL`. Check the [`xlab`](https://github.com/charlesneimog/pd-xlab){:target="_blank"} library for Pure Data.


!!! tip "`<PITCH>` allows differentiation"
    When consecutive techniques occur, the pitch allows the model to distinguish them (e.g., consecutive `tongue-ram` on different notes).

<score id="ptech_example" timeSig="4/4" notes="C4/h:x, D4/h:tongue-ram"></score>

```openscofo
PTECH pizz C4 2
PTECH tongue-ram D4 2
```

---
## `UTECH` 

<score id="utech_example" timeSig="4/4" notes="C4/h:x, C4/h:tongue-ram"></score>

!!! warning "Experimental"
    These events represent extended techniques **without pitch**. They are rendered as glyphs or symbols instead of traditional notes.

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
