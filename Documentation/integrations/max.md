# Max

**OpenScofo** is available in Max as the `openscofo~` external. It analyzes a mono audio signal for real-time score following, reports the current score event and tempo, triggers score actions through Max receivers, and can output audio descriptors for analysis or machine-learning workflows.

## Syntax

```max
openscofo~
openscofo~ rms centroid mfcc
```

Create `openscofo~` with no arguments for score following only. Add descriptor ids as creation arguments when you also want descriptor messages from the analyzed signal.

## Inlets and Outlets

`openscofo~` has one signal inlet.

* **Signal inlet**: mono audio stream analyzed by OpenScofo.
* **Left outlet**: current score event index.
* **Middle outlet**: current tempo in BPM.
* **Right outlet**: descriptor messages, created only when descriptor names are passed as object arguments.

Descriptor messages use the descriptor id as the selector:

```text
mfcc <float...>
centroid <float>
onnx <label> <float>
```

## Basic Use

1. Put an `openscofo~` object in a Max patch.
2. Connect a mono signal to its inlet.
3. Load a score with `score <score-path>`.
4. Start the follower with `start`.
5. Use the left and middle outlets to track the current score event and tempo.

Relative score paths are resolved from the current patch directory.

## Messages

### `score`

Loads a score file using `score <score_path>`.

After a successful load, `openscofo~` resets the current event to `0`, outputs the current BPM and event index, and reads score configuration values such as FFT size and hop size.

### `start`

Starts following the loaded score from event `0`.

`start` clears pending actions, resets the follower to the beginning of the score, outputs the current BPM and event index, and enables following.

### `follow`

Enables or disables score following without reloading the score:

```max
follow 1
follow 0
```

### `set event`

Sets the current score event index:

```max
set event 12
```

### `set verbosity`

Sets the amount of OpenScofo logging printed in the Max console:

* `0`: warning
* `1`: info
* `2`: debug
* `3`: trace

### `set justdescription`

Switches descriptor-only processing on or off:

```max
set justdescription 1
set justdescription 0
```

This is useful when the object is created with descriptor arguments and you want analysis output without score-following behavior.

### `set onnxmodel`

Loads an ONNX model and tells OpenScofo which descriptors were used to train it:

```max
set onnxmodel /path/to/model.onnx rms centroid mfcc
```

The descriptor order must match the order used during model training.

### `get descriptors`

Processes audio from a Max `buffer~` and outputs the requested descriptors once:

```max
get descriptors <buffer-name>
get descriptors <buffer-name> <start-frame-index>
```

`openscofo~` reads one FFT window from the buffer, beginning at the optional start frame. If the buffer ends before the FFT window is full, the remaining samples are padded with zeroes. Descriptor results are sent through the descriptor outlet.

## Score Actions

OpenScofo score actions are delivered to Max receivers. For example, this score fragment:

```text
NOTE C4 1
    sendto event-hit [1 C4]
    delay 1 tempo sendto event-hit [0 C4]
```

can be received in Max with an object named:

```max
receive event-hit
```

If a `sendto` action has no arguments, `openscofo~` sends a bang. If it has arguments, it sends a Max list. Action arguments can be numbers or symbols.

Lua actions in the score are executed by `openscofo~` when the external is built with Lua support. The Max Lua module is available as `max` inside score `LUA { ... }` blocks:

```text
LUA {
    function hello()
        max.post("hello from OpenScofo")
        max.send_bang("event-hit")
    end
}

NOTE C4 1
    luacall(hello())
```

The `max` Lua module provides:

* `max.post(...)`
* `max.error(...)`
* `max.send_bang(receiver)`
* `max.send_float(receiver, value)`
* `max.send_symbol(receiver, value)`
* `max.send_list(receiver, values)`

See [Score Actions](../score/actions.md) and [Lua for Interactive Actions](../score/lua.md) for the score-language side.

## Descriptors

Request descriptors by passing their ids as creation arguments:

```max
openscofo~ rms db centroid mfcc chroma
```

When DSP is running, descriptors are output after each analyzed block. Vector descriptors such as `mfcc`, `chroma`, `logmel`, and `magnitude` output a list of floats. Scalar descriptors output one float.

The full descriptor reference is available in [Descriptors](../descriptors/index.md).
