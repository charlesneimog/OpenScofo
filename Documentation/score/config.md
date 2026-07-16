---
icon: octicons/gear-16
tags:
  - Score Configuration
---

# Configuring a Score

Configuration lines set piece-level behavior. Put them before the first event. `BPM` is the only required configuration.

## Reference Table

| Keyword | Purpose | Syntax | Default | Example |
| --- | --- | --- | --- | --- |
| `BPM` | Tempo for following events | `BPM <NUMBER>` | `60` | `BPM 72` |
| `TRANSPOSE` | Transpose subsequent pitches in semitones | `TRANSPOSE <NUMBER>` | `0` | `TRANSPOSE -12` |
| `TUNINGA4` | Reference tuning for A4 in Hz | `TUNINGA4 <NUMBER>` | `440.0` | `TUNINGA4 442` |
| `SR` | Expected sample rate | `SR <NUMBER>` | host sample rate | `SR 48000` |

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
