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
