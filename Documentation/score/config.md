---
icon: octicons/gear-16
tags:
  - Score Configuration
---

# Configuring a Score

Configuration lines set piece-level behavior. Put global settings before the first event; tempo and temporal-model settings may also be changed at a section boundary. `BPM` is the only required configuration.

## Reference Table

| Keyword | Purpose | Syntax | Default | Example |
| --- | --- | --- | --- | --- |
| `BPM` | Tempo for following events | `BPM <NUMBER>` | `60` | `BPM 72` |
| `TRANSPOSE` | Transpose subsequent pitches in semitones | `TRANSPOSE <NUMBER>` | `0` | `TRANSPOSE -12` |
| `TUNINGA4` | Reference tuning for A4 in Hz | `TUNINGA4 <NUMBER>` | `440.0` | `TUNINGA4 442` |
| `SR` | Expected sample rate | `SR <NUMBER>` | host sample rate | `SR 48000` |
| `SECTIONRESTRICT` | Restrict inference to the selected section | `SECTIONRESTRICT ON\|OFF` | `OFF` | `SECTIONRESTRICT ON` |

## Sections

Declare a section before its events. Names are stored as strings; numeric and quoted forms are supported.

```openscofo
SECTIONRESTRICT ON

SECTION 1
BPM 60
NOTE C4 1

SECTION "A"
BPM 90
NOTE D4 1
```

With `SECTIONRESTRICT ON`, the forward window is the selected section's `FIRSTEVENT` boundary through its last event. Every section receives this boundary, even when it inherits its BPM. When `BPM` precedes the first `SECTION`, its existing leading `FIRSTEVENT` is attached to that first section. Selecting a section lands silently on `FIRSTEVENT`, resets its BPM and temporal model, and allows the first musical event to trigger only when following advances normally. Host integrations can select one with `SetCurrentSection("A")` or the `section A` Pd/Max message.

## Example

```openscofo hl_lines="1 2 3"
SR 48000
BPM 72
TUNINGA4 442

NOTE C4 1
    sendto section [A]
```

## Remarks

For low-level listening, model, and follower settings, see [Advanced Configuration](advanced-config/).
