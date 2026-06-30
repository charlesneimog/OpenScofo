`OpenScofo` associates **actions** with **events**, allowing the host application to react when an event is detected.

Two action types are currently available:

- `sendto`, which communicates with the host environment. Its behavior depends on the integration (Pure Data, Max, Csound, or SuperCollider).
- `luacall`, which executes a Lua function. It is primarily intended for **Pure Data** and **Max**, where Lua provides a convenient way to implement logic that would otherwise be cumbersome using visual patching alone.

All actions can be executed immediately or scheduled using the `delay` keyword.

---

## Scheduling Actions

By default, actions are executed as soon as the associated event is detected. The optional `delay` keyword schedules an action to occur later.

Two kinds of delays are supported:

- **Absolute time**, specified in seconds (`sec`) or milliseconds (`ms`).
- **Relative time**, specified in beats (`tempo`), which automatically follows the performer's estimated tempo.

### Absolute Time

Absolute delays use wall-clock time and are unaffected by tempo changes during the performance.

Examples:

- `delay 2000 ms`
- `delay 2 sec`

```openscofo
NOTE C4 2
    delay 2000 ms sendto e1 [1 2 3 4]
```

In this example, the action is executed exactly two seconds after the `C4` event is detected, regardless of whether the performer is playing faster or slower than the notated tempo.

### Relative Time

Relative delays are expressed in beats and are continuously adjusted using `OpenScofo`'s tempo estimation. This allows delayed actions to remain synchronized with the performer even if the tempo changes.

Examples:

- `delay 1 tempo`
- `delay 0.5 tempo`

```openscofo
NOTE C4 2
    delay 1 tempo sendto e1 [1 2 3 4]
```

In this example, the action is executed one beat after the `C4` event is detected. Because the delay is beat-relative, its duration changes with the performer's tempo: it becomes shorter when the performer plays faster and longer when the performer plays slower.

!!! tip "Use syntax highlighting"
    OpenScofo scores are significantly easier to read and write with syntax highlighting enabled. You can use the **[VS Code extension](https://marketplace.visualstudio.com/items?itemName=charlesneimog.openscofo-language-parse){:target="_blank"}**, the **[Neovim parser](https://github.com/charlesneimog/OpenScofo/tree/main/Sources/Language/nvim){:target="_blank"}**, or the **[OpenScofo Online Editor](https://charlesneimog.github.io/OpenScofo/Editor/){:target="_blank"}** directly in your web browser.

---

## `sendto`

`sendto` communicates with the host application. The destination and interpretation of the transmitted data depend on the integration.

### Pure Data and Max

`sendto` sends data to objects with matching receiver names.

For example,

```openscofo
sendto e1 [1 2 3]
```

sends the list `1 2 3` to a `receive e1` object.

See the [Pure Data](../integrations/puredata.md) and [Max](../integrations/max.md) integration guides for complete examples.

### SuperCollider

`sendto` transmits an OSC message that is received through the `listen` method of an `OpenScofo` instance.

See the [`listen`](../integrations/supercollider.md#2-requesting-audio-descriptors-eg-mfcc) documentation for details.

### Csound

`sendto` schedules an instrument instead of sending a message.

For example,

```openscofo
sendto e1 [0 2 3]
```

executes instrument `e1`, where:

- `0` → `p2` (start delay, in seconds)
- `2` → `p3` (duration, in seconds)
- `3` → `p4` (first user-defined parameter)

This is equivalent to scheduling a normal Csound instrument event.

---

## `luacall`

`luacall` executes a Lua function. Like every other action, it may execute immediately or be scheduled with the `delay` keyword.

Examples:

```openscofo
luacall(myluafunc("hello world"))

delay 1 tempo luacall(myluafunc("hello after one beat"))
```

Global Lua functions are defined inside a `LUA` block:

```lua
LUA {
    function myluafunc(s)
        pd.post(s)
    end
}
```

!!! tip "Using the `pd` and `max` modules"
    When running in **Pure Data** or **Max**, the `pd` and `max` Lua modules provide access to host-specific functionality. See the [Pure Data Lua module](../score/lua/#puredata-lua-module) and the [Max Lua module](../score/lua/#max-lua-module) documentation for details.
