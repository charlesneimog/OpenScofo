# Super Collider

**OpenScofo** is a SuperCollider UGen for real-time score following and audio descriptor extraction. It analyzes an incoming audio stream and communicates with the SuperCollider language (sclang) via OSC messages using `SendReply`.

!!! danger "This is a pre-alpha version!"

## Syntax

```supercollider
OpenScofo.ar(in, sampleRate, fftSize, hopSize)
```

## Description

The `OpenScofo` UGen runs entirely on the server (scsynth) and accepts commands from the client (sclang) to load scores, extract specific audio descriptors, or load ONNX models. Because descriptor arrays (like MFCCs) and string paths cannot be easily output as standard audio/control signals, this UGen heavily utilizes SuperCollider's `UnitCmd` system to receive instructions and `SendReply` to send data back to the client asynchronously.

## Initialization (Inputs)

* **`in`** *(Audio-rate)*: The input audio signal to be analyzed.
* **`sampleRate`** *(Float, Default: 48000.0)*: The expected sample rate.
* **`fftSize`** *(Float, Default: 2048.0)*: The size of the FFT window.
* **`hopSize`** *(Float, Default: 512.0)*: The hop size for analysis blocks.

! warning "Warning"
    The UGen currently outputs a silent audio signal (`0.0`) and is used purely for its side-effects (analysis and OSC replies).

---

## Commands (Client to Server)

You communicate with the `OpenScofo` UGen using `Node.set` or `Synth.cmd`.

### `parseScore`
Loads a text-based score file for the follower.
* **Args:** `path` (String)

### `setFollowScore`
Enables or disables the score following engine.
* **Args:** `follow` (Integer: 1 for true, 0 for false)

### `setEventNotifications`
Enables or disables automatic OSC replies whenever the current score event changes.
* **Args:** `enabled` (Integer: 1 for true, 0 for false)

### `getCurrentEvent`
Requests the current event index from the score follower. The server will reply with an OSC message to `/oscofo/currentEvent`.

### `getDescriptor`
Requests the value(s) of a specific audio descriptor for the current audio block.
* **Args:** `descriptorId` (String, e.g., `"rms"`, `"mfcc"`, `"pitch"`),  check complet list [Descriptors](../descriptors/index.md). 
* **Reply:** The server sends an OSC message to `/oscofo/descriptor/<descriptorId>` containing the float value(s). (Scalar descriptors return 1 float; array descriptors like MFCC return multiple floats).

### `loadOnnxModel`
Loads a custom ONNX machine learning model for advanced descriptor with trained AI model.
* **Args:** 1. `modelPath` (String)
    2. `descriptorIdsCsv` (String: A comma-separated list of descriptor IDs expected by the model).

---

## OSC Replies (Server to Client)

You must set up `OSCFunc` or `OSCdef` listeners in sclang to receive data from the UGen.

* **`/oscofo/currentEvent`**: 
    * Triggered automatically if `setEventNotifications` is 1, OR manually requested via `getCurrentEvent`.
    * **Arguments:** `[ nodeID, replyID, eventIndex ]`
* **`/oscofo/descriptor/<descriptorId>`**: 
    * Triggered by requesting `getDescriptor`.
    * **Arguments:** `[ nodeID, replyID, val1, val2, ... ]`

---

## Examples

### 1. Score Following & Event Notifications

```supercollider
(
// 1. Setup the OSC listener for event updates
OSCdef(\eventTracker, { |msg|
    var eventIndex = msg[3];
    "Current Event Index: %".format(eventIndex).postln;
}, '/oscofo/currentEvent');

// 2. Boot the Synth
x = {
    var sig = SoundIn.ar(0);
    OpenScofo.ar(sig, 48000, 2048, 512);
}.play;
)

// 3. Send commands to the UGen
x.command("parseScore", "/path/to/my/score.txt");
x.command("setFollowScore", 1);
x.command("setEventNotifications", 1); // Start receiving automatic updates

```

### 2. Requesting Audio Descriptors (e.g., MFCC)

```supercollider
(
// 1. Setup the OSC listener for the specific descriptor
OSCdef(\mfccTracker, { |msg|
    // msg[3...] contains the array of MFCC values
    var mfccValues = msg[3..]; 
    "MFCCs: %".format(mfccValues).postln;
}, '/oscofo/descriptor/mfcc');

// 2. Boot the Synth
x = {
    var sig = SoundIn.ar(0);
    OpenScofo.ar(sig, 48000, 2048, 512);
}.play;
)

// 3. Request the descriptor (You would typically put this in a Routine/Task)
Routine({
    loop {
        x.command("getDescriptor", "mfcc");
        0.1.wait; // Request 10 times a second
    }
}).play;
```
