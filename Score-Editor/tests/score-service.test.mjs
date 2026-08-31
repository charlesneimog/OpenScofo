import assert from "node:assert/strict";
import test from "node:test";

import {
    formatOpenScofoNumber,
    generateOpenScofoScore,
    musicXmlDurationToQuarterNotes,
} from "../src/services/score-service.js";

test("converts MusicXML duration units to quarter-note beats", () => {
    assert.equal(musicXmlDurationToQuarterNotes(24, 12), 2);
    assert.equal(musicXmlDurationToQuarterNotes(12, 12), 1);
    assert.equal(musicXmlDurationToQuarterNotes(6, 12), 0.5);
});

test("uses the encoded duration for tuplets without applying their ratio twice", () => {
    assert.equal(musicXmlDurationToQuarterNotes(8, 12), 2 / 3);
});

test("truncates generated numbers to at most three decimal places", () => {
    assert.equal(formatOpenScofoNumber(2 / 3), "0.666");
    assert.equal(formatOpenScofoNumber(1.2349), "1.234");
    assert.equal(formatOpenScofoNumber(2), "2");
    assert.equal(formatOpenScofoNumber(0.5), "0.5");
});

test("generates corrected and truncated MusicXML durations", () => {
    let score = "";
    const context = {
        musicxmlScore: {
            useAIModel: false,
            measures: {
                1: [
                    { bpm: 80, duration: 2, measureNumber: 1, step: "C", octave: 5 },
                    { bpm: 80, duration: 2 / 3, measureNumber: 1, step: "A", octave: 5 },
                    { bpm: 80, duration: 2 / 3, isRest: true, measureNumber: 1 },
                ],
            },
        },
        getPitch(note) {
            return `${note.step}${note.octave}`;
        },
        codeEditor: {
            setValue(value) {
                score = value;
            },
        },
    };

    generateOpenScofoScore.call(context);

    assert.match(score, /BPM 80/);
    assert.match(score, /NOTE C5 2/);
    assert.match(score, /NOTE A5 0\.666/);
    assert.match(score, /REST 0\.666/);
});

test("rejects invalid MusicXML divisions", () => {
    assert.equal(musicXmlDurationToQuarterNotes(24, 0), 0);
    assert.equal(musicXmlDurationToQuarterNotes(24, Number.NaN), 0);
});
