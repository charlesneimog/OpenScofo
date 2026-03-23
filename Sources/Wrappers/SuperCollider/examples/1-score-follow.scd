(
s.waitForBoot {
    var repoRoot = PathName(thisProcess.nowExecutingPath).pathOnly;

    ~audioPath = (repoRoot +/+ "../assets/bwv-1013.wav").standardizePath;
    ~scorePath = (repoRoot +/+ "../assets/bwv-1013.txt").standardizePath;
    ~lastEvent = -1;

    ~sndBuf = Buffer.read(s, ~audioPath);
    s.sync;

    OSCdef(\oscofoCurrentEvent, { |msg|
        var eventIndex = msg[3].asInteger;

        if(eventIndex != ~lastEvent) {
            ~lastEvent = eventIndex;

            ("Current Event Index is: " ++ eventIndex).postln;

            if (eventIndex == 10) {
                "FX ON (smooth)".warn;
                ~oscofoSynth.set(\fxMixTarget, 0.1, \gainTarget, 2.5);
            };

            if (eventIndex == 25) {
                "FX OFF (smooth)".warn;
                ~oscofoSynth.set(\fxMixTarget, 0, \gainTarget, 1.0);
            };

			if (eventIndex == 45) {
                "FX OFF (smooth)".warn;
                ~oscofoSynth.set(\fxMixTarget, 0.2, \gainTarget, 0.4);
            };

			if (eventIndex == 64) {
                "FX OFF (smooth)".warn;
                ~oscofoSynth.set(\fxMixTarget, 0, \gainTarget, 1.0);
            };

        };
    }, '/oscofo/currentEvent');


    SynthDef(\oscofo_test, { |bufnum, fxMixTarget = 0, gainTarget = 1|
        var in, scofo, fx, out;
        var fxMix, gain;

        // dry
        in = PlayBuf.ar(1, bufnum, BufRateScale.kr(bufnum));

        // OpenScofo uses dry
        scofo = OpenScofo.kr(
            sampleRate: 48000,
            fftSize: 2048,
            hopSize: 512,
            in: in
        );

        // effect
        fx = in * SinOsc.ar(200);
        fx = (fx * 10).tanh;
        fx = LPF.ar(fx, 2000);

        // smooth crossfade control (IMPORTANT)
        fxMix = Lag.kr(fxMixTarget, 5.0);   // 5 sec fade

        // crossfade (no jump anymore)
        out = XFade2.ar(in ! 2, fx ! 2, (fxMix * 2) - 1);

        // smooth gain
        gain = Lag.kr(gainTarget, 5.0);

        out = Limiter.ar(out * gain, 0.9);

        Out.ar(0, out);
    }).add;

    s.sync;

    ~oscofoSynth = Synth(\oscofo_test, [\bufnum, ~sndBuf]);
    s.sync;

    OpenScofo.parseScore(~oscofoSynth, ~scorePath);

    if(OpenScofo.respondsTo(\setEventNotifications)) {
        OpenScofo.setEventNotifications(~oscofoSynth, true);
    };

    CmdPeriod.doOnce({
        OSCdef(\oscofoCurrentEvent).free;

        if(OpenScofo.respondsTo(\setEventNotifications)) {
            OpenScofo.setEventNotifications(~oscofoSynth, false);
        };

        if(~oscofoSynth.notNil) { ~oscofoSynth.free; };
        if(~sndBuf.notNil) { ~sndBuf.free; };
    });
};
)
