(
s.waitForBoot {
    var threshold = 0.45;
    ~lastHfr = -1.0;

    OSCdef(\oscofoHfr, { |msg|
        var hfr = msg[3].asFloat;
        if((hfr > threshold) && ((hfr - ~lastHfr).abs > 0.02)) {
            ~lastHfr = hfr;
            ("Bright sound detected. hfr=" ++ hfr.round(0.001)).postln;
        };
    }, '/oscofo/descriptor/hfr');

    SynthDef(\oscofo_descriptors_in, {
        var in = SoundIn.ar(0);
        var scofo = OpenScofo.kr(
            sampleRate: SampleRate.ir,
            fftSize: 2048,
            hopSize: 512,
            in: in
        );
        Out.ar(0, Silent.ar(2));
    }).add;
    s.sync;

    ~oscofoSynth = Synth(\oscofo_descriptors_in);
    OpenScofo.setFollowScore(~oscofoSynth, false);

    ~descriptorPoll = Routine {
        loop {
            OpenScofo.getDescriptor(~oscofoSynth, "hfr");
            OpenScofo.getDescriptor(~oscofoSynth, "zcr");
            0.05.wait;
        }
    }.play(SystemClock);

    CmdPeriod.doOnce({
        if(~descriptorPoll.notNil) { ~descriptorPoll.stop; };
        OSCdef(\oscofoHfr).free;
        if(~oscofoSynth.notNil) { ~oscofoSynth.free; };
    });
};
)
