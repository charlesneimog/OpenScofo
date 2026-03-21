(
s.waitForBoot {
    var repoRoot = PathName(thisProcess.nowExecutingPath).pathOnly;
    4.do { repoRoot = PathName(repoRoot).parentPath };

	// Replace with your audio and score input
    ~audioPath = (repoRoot +/+ "Tests/assets/bwv-1013.wav").standardizePath;
    ~scorePath = (repoRoot +/+ "Tests/assets/bwv-1013.txt").standardizePath;

    ~lastEvent = -1;

    ~sndBuf = Buffer.read(s, ~audioPath);
    s.sync;

    OSCdef(\oscofoCurrentEvent, { |msg|
        var eventIndex = msg[3].asInteger;
        if(eventIndex != ~lastEvent) {
            ~lastEvent = eventIndex;
            ("OpenScofo event: " ++ eventIndex).postln;
        };
    }, '/oscofo/currentEvent');

	// OpenScofo creation
    SynthDef(\oscofo_test, { |bufnum|
        var in, scofo;
        in = PlayBuf.ar(1, bufnum, BufRateScale.kr(bufnum));
        scofo = OpenScofo.kr(
            sampleRate: 48000,
            fftSize: 2048,
            hopSize: 512,
            in: in
        );
        Out.ar(0, in ! 1);
    }).add;
    s.sync;

	//
    ~oscofoSynth = Synth(\oscofo_test, [\bufnum, ~sndBuf]);
    s.sync;

    OpenScofo.parseScore(~oscofoSynth, ~scorePath);
    OpenScofo.setEventNotifications(~oscofoSynth, true);

    CmdPeriod.doOnce({
        OSCdef(\oscofoCurrentEvent).free;
        OpenScofo.setEventNotifications(~oscofoSynth, false);
        if(~oscofoSynth.notNil) { ~oscofoSynth.free; };
        if(~sndBuf.notNil) { ~sndBuf.free; };
    });
};
)
