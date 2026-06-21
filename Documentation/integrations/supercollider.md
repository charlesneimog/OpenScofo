# SuperCollider

<release interface="SuperCollider">Loading Releases</release>

**OpenScofo** is a SuperCollider UGen for real-time score following and audio descriptor extraction. It analyzes an incoming audio stream and communicates with the SuperCollider language (`sclang`) through OSC replies from the server.

!!! danger "This is a pre-alpha version!"

## Syntax

```supercollider
~oscofo = OpenScofo.new(
    scorePath: "/path/to/score.txt",
    inBus: ~analysisBus,
    sampleRate: s.sampleRate,
    namespace: "flute"
);
```

## Description

Create one `OpenScofo` object for each score follower you want to run. Each instance can load its own score, analyze its own input bus, and use its own OSC namespace for score actions. Use methods like `loadScore`, `getDescriptor`, and `listen` from the SuperCollider language; the wrapper handles the server communication for you.

## Initialization (Inputs)

* **`scorePath`** *(String or nil)*: The score file loaded during initialization. Pass `nil` for descriptor-only use.
* **`inBus`** *(Audio bus, Default: 0)*: The mono audio bus to analyze.
* **`sampleRate`** *(Float, Default: 48000.0)*: The expected sample rate.
* **`eventNotifications`** *(Boolean, Default: true)*: Enables automatic current-event replies.
* **`eventAction`** *(Function or nil)*: Optional callback receiving `(eventIndex, msg)` whenever `/<namespace>/currentEvent` is emitted.
* **`namespace`** *(String, Default: `"openscofo"`)*: OSC namespace used for current-event replies and score `sendto` action replies. Paths are generated as `/<namespace>/<receiver>`.

`FFTSIZE` and `HOPSIZE` are read from the score file. They are not SuperCollider constructor arguments.

!!! warning "Warning"
    The UGen is used for analysis and OSC replies. It does not produce an audio signal for listening.

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
Requests the current event index from the score follower. The server will reply with an OSC message to `/<namespace>/currentEvent`.

### `getDescriptor`
Requests the value(s) of a specific audio descriptor for the current audio block.

* **Args:** `descriptorId` (String, e.g., `"rms"`, `"mfcc"`, `"pitch"`),  check complet list [Descriptors](../descriptors/index.md). 

* **Reply:** The server sends an OSC message to `/<namespace>/descriptor/<descriptorId>` containing the float value(s). (Scalar descriptors return 1 float; array descriptors like MFCC return multiple floats).

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

is delivered to `/openscofo/delay` by default and can be handled with:

```supercollider
~oscofo.listen("delay", { |args, msg|
    args.postln;
});
```

`listen` is a convenience method for receiving score actions from the instance namespace plus the score receiver name. Pass only the receiver name, not the full OSC address. For example, use `"delay"` or `"buffer-record"`, not `"/<namespace>/delay"`. The wrapper automatically generates action OSC paths as:

```text
/<namespace>/<receiver>
```

The default namespace is `"openscofo"`. You can set another namespace per instance to avoid action-address collisions:

```supercollider
~score1 = OpenScofo.new(
    scorePath: "/path/to/score1.scofo",
    namespace: "openscofo"
);

~score2 = OpenScofo.new(
    scorePath: "/path/to/score2.scofo",
    namespace: "piece2"
);
```

With the same score receiver `buffer-record`, `~score1` sends to `/openscofo/buffer-record` and `~score2` sends to `/piece2/buffer-record`. In both cases, the public API stays the same:

```supercollider
~score2.listen("buffer-record", { |args|
    args.postln;
});
```

Score action arguments can be floats, ints, or strings when received through `~oscofo.listen(...)`. Internally, string arguments are encoded as floats for SuperCollider's public `SendNodeReply` API and decoded by the `OpenScofo` class before your callback is called. If a receiver has no registered listener, the wrapper prints a warning.

For example, this score fragment:

```text
sendto buffer-record [1 bb4-a4 1]
sendto buffer-record [2 bb4 1]
delay 1 tempo sendto buffer-record [2 bb4 0]
```

can be received in SuperCollider with:

```supercollider
~oscofo.listen("buffer-record", { |args, msg|
    var voice = args[0].asInteger;
    var bufferName = args[1].asString;
    var recording = args[2] != 0;

    ["buffer-record", voice, bufferName, recording].postln;
});
```

The callback receives `args` without the `nodeID` and `replyID` fields. The raw OSC message is still available as `msg`.

---

## OSC Replies (Server to Client)

You must set up `OSCFunc` or `OSCdef` listeners in sclang to receive data from the UGen.

* **`/<namespace>/currentEvent`**: 

    * Triggered automatically if `setEventNotifications` is 1, OR manually requested via `getCurrentEvent`.
    * **Arguments:** `[ nodeID, replyID, eventIndex ]`

* **`/<namespace>/descriptor/<descriptorId>`**: 

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
    namespace: "openscofo",
    eventAction: { |eventIndex| "Current Event Index: %".format(eventIndex).postln; }
);
)
```

### 2. Requesting Audio Descriptors (e.g., MFCC)

```supercollider
(
~namespace = "openscofo";

// 1. Setup the OSC listener for the specific descriptor
OSCdef(\mfccTracker, { |msg|
    // msg[3...] contains the array of MFCC values
    var mfccValues = msg[3..]; 
    "MFCCs: %".format(mfccValues).postln;
}, ("/" ++ ~namespace ++ "/descriptor/mfcc").asSymbol);

~analysisBus = Bus.audio(s, 1);

SynthDef(\analysisInput, { |analysisBus|
    var sig = SoundIn.ar(0);
    Out.ar(analysisBus, sig);
}).add;

s.sync;

~input = Synth(\analysisInput, [\analysisBus, ~analysisBus]);
~oscofo = OpenScofo.new(
    nil,
    inBus: ~analysisBus,
    sampleRate: s.sampleRate,
    namespace: ~namespace,
    eventNotifications: false
);
~oscofo.setFollowScore(false);

Routine({
    loop {
        ~oscofo.getDescriptor("mfcc");
        0.1.wait; // Request 10 times a second
    }
}).play;
```
