---
icon: octicons/book-16
tags:
  - Language Reference
---

# Language Reference

An OpenScofo score is a plain-text `.scofo` file. It combines configuration, musical events, and computer actions.

```openscofo
/* Minimal score example */

BPM 60

NOTE C4 1
    sendto delay [1]

// Extended technique
PTECH tongue-ram D3 1
    delay 1 tempo sendto granular [open]

```

## What's in an OpenScofo Score?

| Element | Purpose | Reference |
| --- | --- | --- |
| Comments | Human-readable notes in the score | [Comments](#comments) |
| Configuration | Tempo, sample rate, and low-level settings | [Configuring a Score](config/) |
| Score Events | What OpenScofo listens for | [Musical Events](events/) |
| Actions | What the computer does when an event is detected | [Computer Actions](actions/) |
| Lua | Complex events that are easier to implement in Lua than in Pd, Max, Csound, etc. | [Lua](lua/) |

## File Extension and Editors

Use `.scofo`.

- VS Code: [OpenScofo Language Parse extension](https://marketplace.visualstudio.com/items?itemName=charlesneimog.openscofo-language-parse){:target="_blank"}
- Neovim: [OpenScofo Neovim configuration](https://github.com/charlesneimog/OpenScofo/tree/main/Sources/Language/nvim){:target="_blank"}
- Browser: [OpenScofo Online Editor](https://charlesneimog.github.io/OpenScofo/Editor/){:target="_blank"}

```openscofo
BPM 120

PHASECOUPLING 0.75

SR 48000
FFTSIZE 2048
HOPSIZE 512

NOTE C4 2
TRILL (G4 C5) 2
    // send 1 to the receiver switch
    sendto action1 [switch 1]
    delay 1.5 tempo sendto action2 [1 2 3]

PTECH tongue-ram C3 0.33
PTECH tongue-ram C#3 0.33
UTECH jet-whistle 0.33

```

See also: [Core Language Concepts](../concepts/core-language-concepts/).
