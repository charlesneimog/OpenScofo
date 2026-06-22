# Lua for Interactive Actions 

Lua is an amazing language to create interactive events. You can read its 88 pages manual on English [here](https://www.lua.org/manual/5.4){target="_blank"} or in Portuguese [here](https://www.lua.org/manual/5.2/pt){target="_blank"}.

I am not present the language here, there are a lot of tutorials on the internet and YouTube.

* Very Complete [here](https://youtu.be/zi-svfdcUc8?si=-JpbobE-0T-zS4_z&t=1393){target="_blank"}.
* Very Fast, Complex but yet Complete tutorial if you already know how to program in another language [here](https://www.youtube.com/watch?v=jUuqBZwwkQw){target="_blank"}.
* More complet tutorial [here](https://www.youtube.com/watch?v=kgiEF1frHQ8){target="_blank"}.
* Text [here](https://tylerneylon.com/a/learn-lua/){target="_blank"}.


Inside `OpenScofo` you can use `Lua` to create interactive events. `OpenScofo` provide a module where you can create deep interactions with the listening and follower module. Also `OpenScofo` Pure Data and Max object provide `pd` and `max` module to interact with the patches. Below I present the module and its methods.

Also check the [LUAEVENT](events/#luaevent){ data-preview } event type.

!!! danger "Under developement yet"
    Lua Module for Pure Data, Max, and OpenScofo are in developement yet.

## `OpenScofo`: Lua Module

The `OpenScofo` module exposes classes and functions to interact with the `OpenScofo` follower and listening module.

#### `OpenScofo` Module

```lua
local oscofo = require("OpenScofo")
```

These are all functions exposed by `oscofo`.

| Function                                 | Parameters                                           | Returns                         | Description                                                |
| ---------------------------------------- | ---------------------------------------------------- | ------------------------------- | ---------------------------------------------------------- |
| `oscofo.set_db_threshold(value)`         | `value`: float value in dB (default: `-60`)          | None                            | Sets the audio threshold used by the tracker.              |
| `oscofo.set_tuning(value)`               | `value`: tuning reference frequency (default: `440`) | None                            | Sets the tuning used for score parsing and tracking.       |
| `oscofo.set_current_event(event)`        | `event`: integer event index (default: `2`)          | None                            | Forces the current score position.                         |
| `oscofo.set_harmonics(value)`            | `value`: integer number of harmonics (default: `10`) | None                            | Sets the number of harmonics used by the pitch template.   |
| `oscofo.set_pitch_template_sigma(value)` | `value`: float sigma value (default: `0.5`)          | None                            | Sets the pitch template sigma.                             |
| `oscofo.get_live_bpm()`                  | None                                                 | `float`                         | Returns the current estimated BPM.                         |
| `oscofo.get_event_index()`               | None                                                 | `integer`                       | Returns the current tracked score event.                   |
| `oscofo.get_states()`                    | None                                                 | Collection of `OpenScofo.State` | Returns all score states currently loaded.                 |
| `oscofo.get_pitch_template(freq)`        | `freq`: frequency in Hz (`float`)                    | Numeric array                   | Returns the internal pitch template for a given frequency. |
| `oscofo.get_audio_description()`         | None                                                 | `OpenScofo.Description`         | Computes MIR/audio features for the current input buffer.  |


### `OpenScofo` Types

#### `OpenScofo.Description`

| Field                | Description                        |
| -------------------- | ---------------------------------- |
| `mfcc`               | MFCC feature vector.               |
| `onset`              | Onset detection value.             |
| `silence_prob`       | Estimated probability of silence.  |
| `spectral_magnitude` | Spectral magnitude representation. |
| `loudness`           | Estimated loudness.                |
| `spectral_flux`      | Spectral flux measure.             |
| `spectral_flatness`  | Spectral flatness measure.         |
| `harmonicity`        | Harmonicity estimate.              |
| `db`                 | Signal level in dB.                |
| `rms`                | Root mean square energy.           |
| `power`              | Signal power estimate.             |

#### `OpenScofo.State`

| Field            | Description                     |
| ---------------- | ------------------------------- |
| `position`       | Score position.                 |
| `type`           | State type.                     |
| `markov`         | Markov state information.       |
| `forward`        | Forward probability value.      |
| `bpm_expected`   | Expected tempo (BPM).           |
| `bpm_observed`   | Observed tempo (BPM).           |
| `onset_expected` | Expected onset time.            |
| `onset_observed` | Observed onset time.            |
| `phase_expected` | Expected phase.                 |
| `phase_observed` | Observed phase.                 |
| `ioi_phi_n`      | Expected inter-onset interval.  |
| `ioi_hat_phi_n`  | Estimated inter-onset interval. |
| `audiostates`    | Associated audio states.        |
| `duration`       | State duration.                 |
| `line`           | Source score line number.       |


## `PureData` Lua Module

The `pd` module provides access to Pure Data messaging and logging functions.

| Function                              | Parameters                              | Returns | Description                                                              |
| ------------------------------------- | --------------------------------------- | ------- | ------------------------------------------------------------------------ |
| `pd.post(message)`                    | `message`: string                       | None    | Posts a message to the Pure Data console at the default verbosity level. |
| `pd.error(message)`                   | `message`: string                       | None    | Posts an error message to the Pure Data console.                         |
| `pd.send_bang(destination)`           | `destination`: string                   | None    | Sends a bang message to a Pure Data receiver.                            |
| `pd.send_float(destination, value)`   | `destination`: string, `value`: float   | None    | Sends a floating-point value to a Pure Data receiver.                    |
| `pd.send_symbol(destination, symbol)` | `destination`: string, `symbol`: string | None    | Sends a symbol to a Pure Data receiver.                                  |
| `pd.send_list(destination, values)`   | `destination`: string, `values`: list   | None    | Sends a list of values to a Pure Data receiver.                          |

### Examples

```lua
pd.post("hello world")
pd.error("this is an error")

pd.send_bang("mybang")
pd.send_float("myfloat", 20)
pd.send_symbol("mysymbol", "hello")
pd.send_list("mylist", {"hello", 1, 2})
```

## `Max` Lua Module

The `max` module provides access to Max messaging and logging functions.

| Function                               | Parameters                              | Returns | Description                                     |
| -------------------------------------- | --------------------------------------- | ------- | ----------------------------------------------- |
| `max.print(message)`                   | `message`: string                       | None    | Prints a message to the Max console.            |
| `max.error(message)`                   | `message`: string                       | None    | Prints an error message to the Max console.     |
| `max.send_bang(destination)`           | `destination`: string                   | None    | Sends a bang message to a Max receiver.         |
| `max.send_float(destination, value)`   | `destination`: string, `value`: float   | None    | Sends a floating-point value to a Max receiver. |
| `max.send_symbol(destination, symbol)` | `destination`: string, `symbol`: string | None    | Sends a symbol to a Max receiver.               |
| `max.send_list(destination, values)`   | `destination`: string, `values`: list   | None    | Sends a list of values to a Max receiver.       |

### Examples

```lua
max.print("hello world")
max.error("this is an error")

max.send_bang("mybang")
max.send_float("myfloat", 20)
max.send_symbol("mysymbol", "hello")
max.send_list("mylist", {1, 2, "hello"})
```

