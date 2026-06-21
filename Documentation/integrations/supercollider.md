# Super Collider

<release interface="SuperCollider">Loading Releases</release>

**OpenScofo** is a SuperCollider UGen for real-time score following and audio descriptor extraction. It analyzes an incoming audio stream and communicates with the SuperCollider language (sclang) via OSC messages using `SendReply`.

!!! danger "This is a pre-alpha version!"

## Syntax

```supercollider
~oscofo = OpenScofo.new(
    scorePath: "/path/to/score.txt",
    inBus: ~analysisBus,
    sampleRate: s.sampleRate
);
```

## Description

`OpenScofo` is a client-side SuperCollider object that owns one server-side `OpenScofoUGen` synth. Each `OpenScofo.new(...)` call creates an independent native OpenScofo instance, so different scores and analysis settings can run at the same time. Descriptor arrays (like MFCCs) and string paths cannot be regular audio/control signals, so the wrapper uses SuperCollider's `UnitCmd` system internally and receives asynchronous `SendReply` OSC messages from the server.

## Initialization (Inputs)

* **`scorePath`** *(String or nil)*: The score file loaded during initialization. Pass `nil` for descriptor-only use.
* **`inBus`** *(Audio bus, Default: 0)*: The mono audio bus to analyze.
* **`outBus`** *(Audio bus, Default: 0)*: A silent output bus used by the internal synth.
* **`sampleRate`** *(Float, Default: 48000.0)*: The expected sample rate.
* **`eventNotifications`** *(Boolean, Default: true)*: Enables automatic current-event replies.
* **`eventAction`** *(Function or nil)*: Optional callback receiving `(eventIndex, msg)` whenever `/openscofo/currentEvent` is emitted.

`FFTSIZE` and `HOPSIZE` are read from the score file. They are not SuperCollider constructor arguments.

!!! warning "Warning"
    The UGen currently outputs a silent audio signal (`0.0`) and is used purely for its side-effects (analysis and OSC replies).

---

## Commands (Client to Server)

You communicate with an `OpenScofo` instance using its methods.

### `loadScore` / `parseScore`
Loads a text-based score file for the follower.

* **Args:** `path` (String)

### `setFollowScore`
Enables or disables the score following engine.

* **Args:** `follow` (Boolean)

### `setEventNotifications`
Enables or disables automatic OSC replies whenever the current score event changes.

* **Args:** `enabled` (Boolean)

### `getCurrentEvent`
Requests the current event index from the score follower. The server will reply with an OSC message to `/openscofo/currentEvent`.

### `getDescriptor`
Requests the value(s) of a specific audio descriptor for the current audio block.

* **Args:** `descriptorId` (String, e.g., `"rms"`, `"mfcc"`, `"pitch"`),  check complet list [Descriptors](../descriptors/index.md). 

* **Reply:** The server sends an OSC message to `/openscofo/descriptor/<descriptorId>` containing the float value(s). (Scalar descriptors return 1 float; array descriptors like MFCC return multiple floats).

### `loadOnnxModel`

Loads a custom ONNX machine learning model for advanced descriptor with trained AI model.

* **Args:** 
    1. `modelPath` (String)
    2. `descriptorIdsCsv` (String: A comma-separated list of descriptor IDs expected by the model).

### `listen`

Registers a SuperCollider listener for score `sendto` actions. A score action such as:

```text
sendto delay [1 0.5 250 0.4]
```

is delivered to `/OpenScofo/delay` and can be handled with:

```supercollider
~oscofo.listen("delay", { |args, msg|
    args.postln;
});
```

SuperCollider score actions can contain float, int, and string arguments when received through `~openscofo.listen(...)`. Internally, string arguments are encoded as floats for SuperCollider's public `SendNodeReply` API and decoded by the `OpenScofo` class before your callback is called. If a receiver has no registered listener, the wrapper prints a warning.

---

## OSC Replies (Server to Client)

You must set up `OSCFunc` or `OSCdef` listeners in sclang to receive data from the UGen.

* **`/openscofo/currentEvent`**: 

    * Triggered automatically if `setEventNotifications` is 1, OR manually requested via `getCurrentEvent`.
    * **Arguments:** `[ nodeID, replyID, eventIndex ]`

* **`/openscofo/descriptor/<descriptorId>`**: 

    * Triggered by requesting `getDescriptor`.
    * **Arguments:** `[ nodeID, replyID, val1, val2, ... ]`

---

## Examples

### 1. Score Following & Event Notifications

```supercollider
(
~analysisBus = Bus.audio(s, 1);

SynthDef(\analysisInput, { |analysisBus|
    var sig = SoundIn.ar(0);
    Out.ar(analysisBus, sig);
}).add;

s.sync;

~input = Synth(\analysisInput, [\analysisBus, ~analysisBus]);
~oscofo = OpenScofo.new(
    scorePath: "/path/to/my/score.txt",
    inBus: ~analysisBus,
    sampleRate: s.sampleRate,
    eventAction: { |eventIndex| "Current Event Index: %".format(eventIndex).postln; }
);
)
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

~analysisBus = Bus.audio(s, 1);

SynthDef(\analysisInput, { |analysisBus|
    var sig = SoundIn.ar(0);
    Out.ar(analysisBus, sig);
}).add;

s.sync;

~input = Synth(\analysisInput, [\analysisBus, ~analysisBus]);
~oscofo = OpenScofo.new(nil, inBus: ~analysisBus, sampleRate: s.sampleRate, eventNotifications: false);
~oscofo.setFollowScore(false);

Routine({
    loop {
        ~oscofo.getDescriptor("mfcc");
        0.1.wait; // Request 10 times a second
    }
}).play;
```
