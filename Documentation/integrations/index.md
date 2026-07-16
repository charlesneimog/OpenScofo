---
icon: material/hub
tags:
  - Platform Integrations
---

# Platform Integrations

OpenScofo runs in real-time hosts, offline analysis tools, browser contexts, and language bindings.

```openscofo
BPM 60

NOTE C4 1
    sendto delay [1]
```

## Reference Table

| Environment | Use |
| --- | --- |
| [Pure Data](puredata/) | live electronics patches |
| [Max/MSP](max/) | live electronics and media patches |
| [Csound](csound/) | instrument scheduling |
| [SuperCollider](supercollider/) | synths, patterns, OSC-style cues |
| [Vamp](vamp/) | offline descriptor analysis |
| [Python](python/) | training and scripting |
| [JavaScript](javascript/) | browser contexts |
| [C++](cpp/) | embedding and development |

## `sendto` Behavior

The score language is shared, but `sendto` is host-dependent.

| Host | Behavior | Example |
| --- | --- | --- |
| Pure Data | sends to `[r receiver]` | `sendto delay [1]` -> `[r delay]` |
| Max | sends to `[receive receiver]` | `sendto delay [1]` -> `[receive delay]` |
| Csound | schedules an instrument event | `sendto 2 [0 0.25 440]` -> `i 2 0 0.25 440` |
| SuperCollider | sends to `/<namespace>/receiver` | `sendto delay [1]` -> `~oscofo.listen("delay", ...)` |
| Python / JavaScript / C++ | exposes score actions through API data | inspect the returned action object |

## Releases

Release are available using [Github](https://github.com/charlesneimog/OpenScofo/releases).

* Installer automatic install all the enviroments (Pd, Max, Csound, etc...); 
* Emscripten is the binary for Web;
* Python is the wheel (is better to install using `pip`);


<release latex="false" interface="All"><i>Loading Releases</i></release>

See also: [Your First Interactive Patch](../getting-started/first-interactive-patch/), [Computer Actions](../score/actions/).
