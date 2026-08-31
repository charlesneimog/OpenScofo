---
icon: simple/lua
tags:
  - Lua Actions
---

# Lua for Actions

Lua can be called from score actions when actions need custom logic.

```openscofo hl_lines="8"
LUA {
    function cue(name)
        pd.post(name)
    end
}

NOTE C4 1
    luacall(cue("section A"))
```

For that, you can check the actions below.

## Reference Table

### `OpenScofo` Module

```lua
local oscofo = require("OpenScofo")
```

| Function | Description |
| --- | --- |
| `oscofo.set_db_threshold(value)` | Sets audio threshold. |
| `oscofo.set_tuning(value)` | Sets tuning reference. |
| `oscofo.set_current_event(event)` | Forces score position. |
| `oscofo.set_current_section(section)` | Resets to the first event of a named section. |
| `oscofo.set_harmonics(value)` | Sets pitch-template harmonics. |
| `oscofo.set_pitch_template_sigma(value)` | Sets pitch tolerance. |
| `oscofo.get_live_bpm()` | Returns estimated BPM. |
| `oscofo.get_event_index()` | Returns current event index. |
| `oscofo.get_states()` | Returns current score states. |
| `oscofo.get_pitch_template(freq)` | Returns pitch template for a frequency. |
| `oscofo.get_audio_description()` | Returns current audio descriptors. |

### Host Modules

| Module | Functions |
| --- | --- |
| `pd` | `post`, `error`, `send_bang`, `send_float`, `send_symbol`, `send_list` |
| `max` | `print`, `error`, `send_bang`, `send_float`, `send_symbol`, `send_list` |

```lua
pd.send_list("section", {"A", 1})
max.send_float("reverb", 0.8)
```

## Remarks

See [Computer Actions](actions/) for `luacall` syntax.
