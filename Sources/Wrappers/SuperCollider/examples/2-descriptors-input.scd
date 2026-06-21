(
s.waitForBoot {
    var threshold = 0.45;
    ~lastHfr = -1.0;
    ~namespace = "openscofo";
    ~analysisBus = Bus.audio(s, 1);

    OSCdef(\oscofoHfr, { |msg|
        var hfr = msg[3].asFloat;
        if((hfr > threshold) && ((hfr - ~lastHfr).abs > 0.02)) {
            ~lastHfr = hfr;
            ("Bright sound detected. hfr=" ++ hfr.round(0.001)).postln;
        };
    }, ("/" ++ ~namespace ++ "/descriptor/hfr").asSymbol, s.addr);

    SynthDef(\oscofo_mic_input, { |analysisBus|
        Out.ar(analysisBus, SoundIn.ar(0));
        Out.ar(0, Silent.ar(2));
    }).add;

    s.sync;

    ~micSynth = Synth(\oscofo_mic_input, [\analysisBus, ~analysisBus]);
    ~oscofo = OpenScofo.new(
        scorePath: nil,
        inBus: ~analysisBus,
        sampleRate: s.sampleRate,
        namespace: ~namespace,
        eventNotifications: false
    );
    ~oscofo.setFollowScore(false);

    ~descriptorPoll = Routine {
        loop {
            ~oscofo.getDescriptor("hfr");
            ~oscofo.getDescriptor("zcr");
            0.05.wait;
        }
    }.play(SystemClock);

    CmdPeriod.doOnce({
        if(~descriptorPoll.notNil) { ~descriptorPoll.stop; };
        OSCdef(\oscofoHfr).free;
        if(~oscofo.notNil) { ~oscofo.free; };
        if(~micSynth.notNil) { ~micSynth.free; };
        if(~analysisBus.notNil) { ~analysisBus.free; };
    });
};
)
