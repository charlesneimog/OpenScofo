---
icon: custom/supercollider
tags:
  - Host Integration
  - SuperCollider
---

# SuperCollider

## Overview

Use `OpenScofo` in SuperCollider for real-time score following, OSC-style actions, and descriptor requests from `sclang`.

## Installation

Install the SuperCollider wrapper from the OpenScofo release matching your platform and make sure the class is available to SuperCollider.

## Minimal Example

```supercollider
~oscofo = OpenScofo.new(
    scorePath: "/path/to/score.scofo",
    inBus: ~analysisBus,
    sampleRate: s.sampleRate,
    namespace: "flute"
);

~oscofo.listen("delay", { |args, msg|
    args.postln;
});
```

## Reference Table

### Messages / API

| Method / option | Purpose |
| --- | --- |
| `OpenScofo.new(...)` | create a score follower |
| `scorePath` | score loaded at initialization; use `nil` for descriptor-only use |
| `inBus` | mono audio bus to analyze |
| `sampleRate` | expected sample rate |
| `namespace` | OSC namespace for event and action replies |
| `eventNotifications` | enable current-event replies |
| `eventAction` | callback for current-event replies |
| `loadScore(path)` | load a score |
| `setFollowScore(bool)` | enable or disable following |
| `getCurrentEvent()` | request current event |
| `getDescriptor(id)` | request a descriptor |
| `loadOnnxModel(path, descriptorIdsCsv)` | load an ONNX model |
| `listen(receiver, function)` | register a score-action receiver |

## Score Actions

`sendto receiver [...]` is delivered to `/<namespace>/receiver`.

```openscofo
NOTE C4 1
    sendto delay [1 0.5 250 0.4]
```

```supercollider
~oscofo.listen("delay", { |args, msg|
    args.postln;
});
```

Use unique namespaces for multiple followers:

```supercollider
~o_flute = OpenScofo.new(scorePath: "flute.scofo", namespace: "flute");
~o_clarinet = OpenScofo.new(scorePath: "clarinet.scofo", namespace: "clarinet");
```

For shared action syntax, see [Computer Actions](../score/actions/).

## Descriptors

```supercollider
OSCdef(\mfccTracker, { |msg|
    var mfccValues = msg[3..];
    mfccValues.postln;
}, "/openscofo/descriptor/mfcc");

~oscofo.getDescriptor("mfcc");
```

Scalar descriptors return one float. Vector descriptors such as `mfcc` return multiple floats.

## Complete Example

```supercollider
(
s.options.sampleRate = 48000;
s.waitForBoot {
    ~bus = Bus.audio(s, 1);
    ~namespace = "flute";
    ~scriptDir = PathName(thisProcess.nowExecutingPath).pathOnly;
    ~buf = Buffer.read(s, ~scriptDir +/+ "miniatura1.mp3");

    SynthDef(\play, { |buf, bus, amp = 0.5|
        var sig = PlayBuf.ar(1, buf, BufRateScale.kr(buf), doneAction: 2);
        Out.ar(bus, sig);
        Out.ar(0, (sig * amp) ! 2);
    }).add;

    s.sync;

    ~oscofo = OpenScofo.new(
        scorePath: ~scriptDir +/+ "miniatura1.scofo",
        inBus: ~bus,
        sampleRate: s.sampleRate,
        namespace: ~namespace,
        eventNotifications: true,
        eventAction: { |eventIndex| "Event: %".format(eventIndex).postln; }
    );

    ~oscofo.listen("buffer-record", { |args, msg|
        args.postln;
    });

    ~player = Synth(\play, [\buf, ~buf, \bus, ~bus, \amp, 0.7]);
};
)
```

## Remarks

Pass only the receiver name to `listen`, such as `"delay"`, not the full OSC path.
