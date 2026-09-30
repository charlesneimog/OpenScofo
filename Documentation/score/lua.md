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
| `oscofo.activate_all_descriptors()` | Enables all descriptors before audio processing; ONNX requires a loaded model. |
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
| `oscofo.schedule(delay_ms, callback, data)` | Schedules a Lua callback and returns a unique timer ID. |
| `oscofo.cancel(timer_id)` | Cancels a pending timer; returns `true` if found and cancelled, or `false` otherwise. |

#### Scheduled Callbacks

Use `openscofo.schedule(delay_ms, callback, data)` to execute a Lua function after a delay:

| Parameter | Description |
| --- | --- |
| `delay_ms` | Non-negative, finite delay in milliseconds. |
| `callback` | Lua function to execute. |
| `data` | Optional Lua value passed as `callback(data)`. May be any Lua value, including a table; omitted data becomes `nil`. |

`schedule()` returns a unique timer ID. Multiple callbacks may be scheduled for exactly the same time; they execute in scheduling order. A zero delay runs at the next timer-processing opportunity, never directly inside `schedule()`.

Timer execution has audio-block precision, determined by the actual audio block size rather than a separate high-resolution timer or MIR analysis hop size. Callbacks execute at the first block boundary at or after their requested time. With the usual 64-sample block size at 48 kHz, the resolution is approximately 1.33 ms; other block sizes change this resolution. Time advances only while audio blocks are processed.

!!! warning "Keep scheduled callbacks lightweight"
    This scheduler is intended for small control/event operations, not heavy computation.

    Scheduled Lua callbacks run as part of OpenScofo's audio processing schedule. Avoid expensive or long-running operations inside callbacks, as they can delay audio processing and potentially cause audio dropouts.

    Keep callbacks lightweight: they should generally perform small, fast operations. Avoid heavy computation, blocking operations, sleeps, large file operations, or other work that may take an unpredictable amount of time.

This example passes a table to a callback scheduled after 500 milliseconds:

```lua
local openscofo = require("OpenScofo")

local id = openscofo.schedule(500, function(data)
    print(data.message)
end, {
    message = "hello"
})
```

Use `openscofo.cancel(timer_id)` to cancel a pending callback. It returns `true` if the pending timer was found and cancelled, or `false` otherwise:

```lua
local openscofo = require("OpenScofo")

local id = openscofo.schedule(5000, function(data)
    print(data)
end, "hello")

openscofo.cancel(id)
```

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
