(
s.waitForBoot {
    var repoRoot = PathName(thisProcess.nowExecutingPath).pathOnly;

    ~audioPath = "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.wav";
    ~scorePath = "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.txt";
    ~lastEvent = -1;

    ~sndBuf = Buffer.read(s, ~audioPath);
    s.sync;

    // Listen for OpenScofo currentEvent messages
    OSCdef(\oscofoCurrentEvent, { |msg|
        var eventIndex = msg[3].asInteger;
        if(eventIndex != ~lastEvent) {
            ~lastEvent = eventIndex;
            ("Current Event Index: " ++ eventIndex).postln;

            switch(eventIndex,
                10, { "FX ON (smooth)".warn; ~oscofoSynth.set(\fxMixTarget, 0.1, \gainTarget, 2.5); },
                25, { "FX OFF (smooth)".warn; ~oscofoSynth.set(\fxMixTarget, 0, \gainTarget, 1.0); },
                45, { "FX MOD (smooth)".warn; ~oscofoSynth.set(\fxMixTarget, 0.2, \gainTarget, 0.4); },
                64, { "FX OFF (smooth)".warn; ~oscofoSynth.set(\fxMixTarget, 0, \gainTarget, 1.0); }
            );
        };
    }, '/oscofo/currentEvent');

    // SynthDef
    SynthDef(\oscofo_test, { |bufnum, fxMixTarget = 0, gainTarget = 1|
        var in, scofo, fx, out;
        var fxMix, gain;

        in = PlayBuf.ar(1, bufnum, BufRateScale.kr(bufnum));

        scofo = OpenScofo.kr(
            in: in,
            sampleRate: 48000,
            fftSize: 2048,
            hopSize: 512
        );

        // Example FX
        fx = in * SinOsc.ar(200);
        fx = (fx * 10).tanh;
        fx = LPF.ar(fx, 2000);

        // Smooth crossfade control
        fxMix = Lag.kr(fxMixTarget, 5.0);
        out = XFade2.ar(in ! 2, fx ! 2, (fxMix * 2) - 1);

        // Smooth gain
        gain = Lag.kr(gainTarget, 5.0);
        out = Limiter.ar(out * gain, 0.9);

        Out.ar(0, out);
    }).add;

    s.sync;

    ~oscofoSynth = Synth(\oscofo_test, [\bufnum, ~sndBuf]);
    s.sync;

    // Send score to OpenScofo
    scofo.parseScore(~oscofoSynth, ~scorePath);

    // Enable event notifications
    if(OpenScofo.respondsTo(\setEventNotifications)) {
        OpenScofo.setEventNotifications(~oscofoSynth, true);
    };
};
)
