(
s.waitForBoot {
    var exampleDir = PathName(thisProcess.nowExecutingPath).pathOnly;
    var localAudioPath = exampleDir +/+ "bwv-1013.wav";
    var localScorePath = exampleDir +/+ "bwv-1013.txt";

    ~audioPath = if(File.exists(localAudioPath), {
        localAudioPath
    }, {
        (exampleDir +/+ "../../../../Tests/assets/bwv-1013.wav").standardizePath
    });
    ~scorePath = if(File.exists(localScorePath), {
        localScorePath
    }, {
        (exampleDir +/+ "../../../../Tests/assets/bwv-1013.txt").standardizePath
    });
    ~analysisBus = Bus.audio(s, 1);
    ~sndBuf = Buffer.read(s, ~audioPath);
    ~lastEvent = -1;

    SynthDef(\oscofo_player, { |bufnum, analysisBus, fxMixTarget = 0, gainTarget = 1|
        var in, fx, out, fxMix, gain;

        in = PlayBuf.ar(1, bufnum, BufRateScale.kr(bufnum), doneAction: Done.freeSelf);
        Out.ar(analysisBus, in);

        fx = in * SinOsc.ar(200);
        fx = (fx * 10).tanh;
        fx = LPF.ar(fx, 2000);

        fxMix = Lag.kr(fxMixTarget, 5.0);
        out = XFade2.ar(in ! 2, fx ! 2, (fxMix * 2) - 1);

        gain = Lag.kr(gainTarget, 5.0);
        out = Limiter.ar(out * gain, 0.9);

        Out.ar(0, out);
    }).add;

    s.sync;

    ~playerSynth = Synth(\oscofo_player, [
        \bufnum, ~sndBuf,
        \analysisBus, ~analysisBus
    ]);

    ~oscofo = OpenScofo.new(
        scorePath: ~scorePath,
        inBus: ~analysisBus,
        sampleRate: s.sampleRate,
        eventNotifications: true,
        eventAction: { |eventIndex|
            if(eventIndex != ~lastEvent) {
                ~lastEvent = eventIndex;
                ("Current Event Index: " ++ eventIndex).postln;

                switch(eventIndex,
                    10, { "FX ON (smooth)".warn; ~playerSynth.set(\fxMixTarget, 0.1, \gainTarget, 2.5); },
                    25, { "FX OFF (smooth)".warn; ~playerSynth.set(\fxMixTarget, 0, \gainTarget, 1.0); },
                    45, { "FX MOD (smooth)".warn; ~playerSynth.set(\fxMixTarget, 0.2, \gainTarget, 0.4); },
                    64, { "FX OFF (smooth)".warn; ~playerSynth.set(\fxMixTarget, 0, \gainTarget, 1.0); }
                );
            };
        }
    );

    CmdPeriod.doOnce({
        if(~oscofo.notNil) { ~oscofo.free; };
        if(~playerSynth.notNil) { ~playerSynth.free; };
        if(~sndBuf.notNil) { ~sndBuf.free; };
        if(~analysisBus.notNil) { ~analysisBus.free; };
    });
};
)
