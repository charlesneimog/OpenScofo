---
icon: fontawesome/solid/computer
tags:
  - Actions
---

# Actions

Actions run when an event is detected. Define them one line below the event.

## Reference Table

| Action | Syntax | Arguments | Example | Remarks |
| --- | --- | --- | --- | --- |
| `sendto` | `sendto <DESTINATION> [<MESSAGE>]` | receiver name, message list | `sendto delay [1]` | Host-specific delivery; see [Platform Integrations](../integrations/#sendto-behavior). |
| `delay` | `delay <VALUE> <UNIT> <ACTION>` | `ms`, `sec`, or `tempo` | `delay 1 tempo sendto echo [1]` | Schedules an action. `tempo` is performer beat-relative. |
| `luacall` | `luacall(<FUNCTION_CALL>)` | Lua function call | `luacall(cue("A"))` | Useful for small logic before sending. |

## Audio state changes

Declare a receiver once at score level to receive every followed audio-state change:

```openscofo
ONAUDIOSTATECHANGE audiostates
```

The forward model stores the best audio observation for every candidate event. After choosing the best event, it
compares that event and its best observation with the previous result. A changed event always sends a notification;
within the same event, a different best observation also sends one. Identical results do not send duplicates.

The receiver uses the same host-specific delivery mechanism as `sendto`. Each message contains the current score
position followed by the selected observation: its MIDI pitch, ONNX label, `silence`, or `onset`. Chords are reported
as `chord` followed by every MIDI pitch in score order; their notification remains stable when the spectral balance
between constituent notes changes.

In Csound the listener is an instrument: notifications use `p2 = 0` and `p3 = 0`, and the payload shown below
starts at `p4`. Other hosts receive the payload as a normal list.

```text
1 60
2 chord 60 64 67
3 pizz
3 silence
```

## Example

<div class="grid cards" markdown>

- __Immediate action__
    ```openscofo hl_lines="2"
    NOTE C4 1
        sendto freeze [1]
    ```

- __Tempo-relative delay__
    ```openscofo hl_lines="2"
    NOTE D4 1
        delay 1 tempo sendto granular [open]
    ```

- __Absolute delay__
    ```openscofo hl_lines="2"
    NOTE E4 1 
        delay 500 ms sendto video [fade_in]
    ```

- __Lua action__
    ```openscofo hl_lines="2"
    NOTE C4 1
        luacall(cue("section A"))
    ```

</div>

## Remarks

Use `tempo` delays for performer-relative timing and `ms` or `sec` for fixed technical timing. See also: [Core Language Concepts](../concepts/core-language-concepts/), [Lua](lua/), [Platform Integrations](../integrations/).
