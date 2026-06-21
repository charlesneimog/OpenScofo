(
s.waitForBoot {
    ~bus = Bus.audio(s, 1);
    ~buf = Buffer.read(s, "/home/neimog/Documents/Git/OpenScofo/Tests/miniaturas/Audios/miniatura1.mp3");

    SynthDef(\play, { |buf, bus, amp = 0.5|
        var sig = PlayBuf.ar(1, buf, BufRateScale.kr(buf), doneAction: 2);

        Out.ar(bus, sig);              // to OpenScofo
        Out.ar(0, (sig * amp) ! 2);    // to speakers
    }).add;

    s.sync;

    ~oscofo = OpenScofo.new(
        scorePath: "/home/neimog/Documents/Git/OpenScofo/Tests/miniaturas/Extras/miniatura1.scofo",
        inBus: ~bus,
        sampleRate: s.sampleRate,
        eventNotifications: true
    );

	~oscofo.listen("buffer-record", { |args, msg|
		// This must execute something
		var voice = args[0].asInteger;
		var bufferName = args[1].asString;
		var recording = args[2] != 0;

		["buffer-record", voice, bufferName, recording].postln;
	});

	~oscofo.listen("buffer-play", { |args, msg|
		// This must execute something
		["buffer-play"].postln;
	});

    ~player = Synth(\play, [\buf, ~buf, \bus, ~bus, \amp, 0.7]);
};
)
