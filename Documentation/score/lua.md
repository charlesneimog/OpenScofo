# :simple-lua: Lua for Interactive Actions 

Lua is an amazing language to create interactive events. You can read its 88 pages manual on English [here](https://www.lua.org/manual/5.4){target="_blank"} or in Portuguese [here](https://www.lua.org/manual/5.2/pt){target="_blank"}.

I am not present the language here, there are a lot of tutorials on the internet and YouTube.

* Very Complete [here](https://youtu.be/zi-svfdcUc8?si=-JpbobE-0T-zS4_z&t=1393){target="_blank"}.
* Very Fast, Complex but yet Complete tutorial if you already know how to program in another language [here](https://www.youtube.com/watch?v=jUuqBZwwkQw){target="_blank"}.
* More complet tutorial [here](https://www.youtube.com/watch?v=kgiEF1frHQ8){target="_blank"}.
* Text [here](https://tylerneylon.com/a/learn-lua/){target="_blank"}.


All these modules are for the `OpenScofo` language; What I will present is how to create actions using `lua` and `oscofo`, `pd` and `max` lua module.

!!! danger "Under developement yet"
    Lua Module for Pure Data, Max, and OpenScofo are in developement yet.


### <h2 align="center">:simple-lua: `OpenScofo` Lua Module</h2>

The `OpenScofo` module exposes classes and functions.

### Creating an `OpenScofo` object

```lua
local oscofo = require("OpenScofo")
```

### `OpenScofo` methods

These are all functions exposed using `oscofo`.

- **`oscofo.set_db_threshold(value)`**  
    - `input`: float value in dB (default: -60).  
    - `output`: no output.  
    - `description`: Sets audio threshold used by the tracker.

- **`oscofo.set_tuning(value)`**  
    - `input`: tuning reference value (default: 440).  
    - `output`: no output.  
    - `description`: Sets tuning used for score parsing/tracking.

- **`oscofo.set_current_event(event)`**  
    - `input`: integer event index/position (default: 2).  
    - `output`: no output.  
    - `description`: Forces the current score position.

- **`oscofo.set_harmonics(value)`**  
    - `input`: integer number of harmonics (default: 10).  
    - `output`: no output.  
    - `description`: Sets number of harmonics used by the pitch template.

- **`oscofo.set_pitch_template_sigma(value)`**  
    - `input`: float sigma value (default: 0.5).  
    - `output`: no output.  
    - `description`: Sets pitch template sigma.

- **`oscofo.get_live_bpm()`**  
    - `input`: no input.  
    - `output`: current estimated BPM (float).  
    - `description`: Returns live tempo estimation.

- **`oscofo.get_event_index()`**  
    - `input`: no input.  
    - `output`: current event index (integer).  
    - `description`: Returns the current tracked score event.

- **`oscofo.get_states()`**  
    - `input`: no input.  
    - `output`: return collection of `OpenScofo.State` of the current score.  
    - `description`: Returns all score states currently loaded.

- **`oscofo.get_pitch_template(freq)`**  
    - `input`: frequency in Hz (float).  
    - `output`: numeric array with the pitch template.  
    - `description`: Returns internal pitch template for a given frequency.

- **`oscofo.get_audio_description()`**  
    - `input`: no input.
    - `output`: `OpenScofo.Description` object of the current audio buffer.  
    - `description`: Computes MIR/audio features for the input block.

### Exposed types in `oscofo`

- **`OpenScofo.Description`**
    - `mfcc`
    - `onset`
    - `silence_prob`
    - `spectral_power`
    - `norm_spectral_power`
    - `pseudo_cqt`
    - `loudness`
    - `spectral_flux`
    - `spectral_flatness`
    - `harmonicity`
    - `db`
    - `rms`
    - `power`

- **`OpenScofo.State`**
    - `position`
    - `type`
    - `markov`
    - `forward`
    - `bpm_expected`
    - `bpm_observed`
    - `onset_expected`
    - `onset_observed`
    - `phase_expected`
    - `phase_observed`
    - `ioi_phi_n`
    - `ioi_hat_phi_n`
    - `audiostates`
    - `duration`
    - `line`

### <h2 align="center">:simple-lua: `PureData` Lua Module</h2>

The `pd` module inside Lua allows interaction with Pure Data functionalities, exposing the following functions:

---
- **`pd.post("hello world")`**  
    - **`input`**: A string message to post.  
    - **`output`**: no output.
    - **`description`**: Posts a message at the default verbosity level, used for general logging.

---
- **`pd.error("this is a message error")`**  
    - **`input`**: A string error message.  
    - **`output`**: no output.
    - **`description`**: Logs an error message in the Pure Data console, similar to a `print` with error severity.

---
- **`pd.send_bang("mybang")`**  
    - **`input`**: A string representing the destination symbol in Pure Data.  
    - **`output`**: no output.
    - **`description`**: Sends a bang message to a specified destination in Pure Data.

---
- **`pd.send_float("myfloat", 20)`**  
    - **`input`**: A float value and a string representing the destination symbol in Pure Data.  
    - **`output`**: no output.
    - **`description`**: Sends a floating-point number to a specified destination in Pure Data.

---
- **`pd.send_symbol("mysymbol", "ola")`**  
    - **`input`**: A string symbol and a string representing the destination symbol in Pure Data.  
    - **`output`**: no output.
    - **`description`**: Sends a symbol to a specified destination in Pure Data.

---
- **`pd.send_list("mylist", {"oi", 1, 2})`**  
    - **`input`**: A list of values (mixed types) and a string representing the destination symbol in Pure Data.  
    - **`output`**: no output.
    - **`description`**: Sends a list of values to a specified destination in Pure Data.

### <h2 align="center">:simple-lua: `Max` Lua Module</h2>

The `max` module inside Lua allows interaction with Max functionalities, exposing the following functions:

- **`max.print("hello world")`**  
    - **`input`**: A string message to print.  
    - **`output`**: no output.
    - **`description`**: Logs a message to the console, similar to Max's `print` object.

---
- **`max.error("this is an error")`**  
    - **`input`**: A string error message.  
    - **`output`**: no output.
    - **`description`**: Logs an error message in the Max console, similar to a `print` with error severity.

---
- **`max.send_bang("mybang")`**  
    - **`input`**: A string representing the receiver symbol in Max.  
    - **`output`**: no output.
    - **`description`**: Sends a bang message to a specified receiver in Max.

---
- **`max.send_float("myfloat", 20)`**  
    - **`input`**: A float value and a string representing the destination symbol in Max.  
    - **`output`**: no output.
    - **`description`**: Sends a floating-point number to a specified destination in Max.

---
- **`max.send_symbol("mysymbol", "oi")`**  
    - **`input`**: A string symbol and a string representing the destination symbol in Max.  
    - **`output`**: no output.
    - **`description`**: Sends a symbol to a specified destination in Pure Data.

---
- **`max.send_list("mylist", {1, 2, "oi"})`**  
    - **`input`**: A list of values (mixed types) and a string representing the destination symbol in Max.  
    - **`output`**: no output.
    - **`description`**: Sends a list of values to a specified destination in Max.
