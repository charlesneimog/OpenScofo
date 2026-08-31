---
icon: custom/max
tags:
  - Host Integration
  - Max
---

# Max

## Overview

Use `openscofo~` to follow a score from a mono signal, report event/BPM, deliver actions to Max receivers, and optionally output descriptors.

| Object | Use |
| --- | --- |
| `openscofo~` | score following |
| `openscofo~ rms centroid mfcc` | score following plus descriptor output |

## Installation

Install the Max external from the OpenScofo release matching your platform, then restart Max if needed.

## Minimal Example

Relative score paths are resolved from the current patch directory.

## Reference Table

### Messages / API

| Message | Purpose |
| --- | --- |
| `score <path>` | load a `.scofo` file |
| `start` | reset to event `0` and start following |
| `follow 1`, `follow 0` | enable or disable following |
| `section <name>` | silently reset following to a section; only the BPM outlet is updated |
| `set event <index>` | force the current score event |
| `set verbosity <0..3>` | set console logging level |
| `set description <0 or 1>` | enable descriptor-only processing |
| `set onnxmodel <path> <descriptors...>` | load an ONNX model with descriptor order |
| `get descriptors <buffer> [start-frame]` | analyze one FFT window from a Max `buffer~` |

Outlets:

| Outlet | Output |
| --- | --- |
| left | current score event index |
| middle | current BPM |
| right | descriptor messages, only when descriptors are requested |

## Score Actions

```openscofo
NOTE C4 1
    sendto delay [1 C4]
    delay 1 tempo sendto delay [0 C4]
```

Receive in Max:

```text
r delay
```

No-argument `sendto` sends a bang; arguments are sent as a Max list. For shared action syntax, see [Computer Actions](../score/actions/).

## Descriptors

Request descriptors as creation arguments:

```text
openscofo~ rms db centroid mfcc chroma
```

Vector descriptors such as `mfcc`, `chroma`, `logmel`, and `magnitude` output lists. Scalar descriptors output one float.

## Complete Example

See [Your First Interactive Patch](../getting-started/first-interactive-patch/) for the smallest complete patch.

## Remarks

Lua actions run when the external is built with Lua support. Inside `LUA { ... }`, use the `max` module; see [Lua](../score/lua/).
