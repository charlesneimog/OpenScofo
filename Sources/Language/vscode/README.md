# OpenScofo Language for VS Code

Syntax highlighting for `.scofo` scores written in the OpenScofo score language.

OpenScofo scores describe musical events and the actions that should happen when
those events are followed in performance. This extension makes those files easier
to read and edit in Visual Studio Code.

## Features

- Syntax highlighting for OpenScofo score files (`.scofo`)
- Highlighting for configuration keywords such as `BPM`, `FFTSIZE`, `HOPSIZE`,
  `TRANSPOSE`, `ONNXMODEL`, `ONNXDESCRIPTORS`, and `ONSETFUNCTION`
- Event highlighting for `NOTE`, `REST`, `CHORD`, `TRILL`, `PTECH`, `UTECH`, and
  `LUAEVENT`
- Action highlighting for `sendto`, `delay`, and `luacall`
- Receiver and argument highlighting for OpenScofo actions
- Pitch, duration, time unit, descriptor, path, and string highlighting
- Embedded Lua highlighting inside `LUA { ... }` blocks
- Line comments (`//`) and block comments (`/* ... */`)

## Example

```scofo
// Configuration
BPM 120
FFTSIZE 2048
HOPSIZE 512
ONSETFUNCTION hfc

LUA {
    function announce(note)
        print("event: " .. note)
    end
}

NOTE C4 1
    sendto synth [440]

NOTE E4 1
    delay 100 ms sendto synth [440]

CHORD (C4 E4 G4) 2
    sendto harmony [major]

PTECH [aeolian, sing] C4 2
    sendto electronics [aeolian],
    luacall(announce("aeolian"))
```

## Getting Started

1. Install **OpenScofo Language Parse** from the VS Code Marketplace.
2. Open a file with the `.scofo` extension.
3. VS Code will automatically use the `OpenScofo` language mode.

If the language mode is not selected automatically, open the language selector in
the lower-right corner of VS Code and choose **OpenScofo**.

## About OpenScofo Scores

An OpenScofo score is a plain text file made of three main parts:

- **Configuration**, for audio, tempo, pitch, onset detection, and AI descriptor
  settings.
- **Events**, such as notes, rests, chords, trills, instrumental techniques, and
  Lua-defined events.
- **Actions**, such as sending messages to receivers, scheduling delayed
  commands, or calling Lua functions.

The recommended file extension is `.scofo`.

## Documentation

- OpenScofo documentation: https://charlesneimog.github.io/OpenScofo/
- Score language introduction:
  https://charlesneimog.github.io/OpenScofo/score/intro/
- Source repository: https://github.com/charlesneimog/OpenScofo

## Current Scope

This extension currently provides syntax highlighting and language association
for `.scofo` files. It does not yet provide diagnostics, autocomplete, formatting,
or score validation.

## License

GPL-3.0-only
