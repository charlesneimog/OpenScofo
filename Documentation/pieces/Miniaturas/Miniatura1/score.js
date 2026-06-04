import * as VexFlow from "https://cdn.jsdelivr.net/npm/vexflow@5.0.0/+esm";

let eventCount = 0;

// ─────────────────────────────────────
function xmlPitchToVexStringPitch(type) {
    const types = {
        whole: "w",
        half: "h",
        quarter: "q",
        eighth: "8",
        "16th": "16",
        "32nd": "32",
        "64th": "64",
        "128th": "128",
    };

    return types[type] ?? alert("Note rhythm value not implemented");
}

// ─────────────────────────────────────
function renderNotehead(staveNote, note, type) {
    var notehead = note.querySelector("notehead")?.textContent;
    if (notehead) {
        notehead = notehead.charAt(0).toUpperCase() + notehead.slice(1);
        if (type === "whole") {
            type = "Whole";
        } else if (type === "half") {
            type = "Half";
        } else {
            type = "Black";
        }

        if (notehead == "Cross") {
            notehead = "Plus";
        } else if (notehead == "Triangle") {
            notehead = "TriangleUp";
        }

        let noteHeadString = `notehead${notehead}${type}`;
        let glyph = VexFlow.VexFlow.Glyphs[noteHeadString];
        staveNote.myNotehead = glyph;
    }
}

// ─────────────────────────────────────
function xmlToVexflowNote(note) {
    const type = note.querySelector("type")?.textContent;

    // Measure rests often have no <type>
    if (!type && note.querySelector("rest")) {
        const duration = Number(note.querySelector("duration")?.textContent ?? 0);

        let vexType = "w";
        if (duration === 6) vexType = "h";
        if (duration === 3) vexType = "q";
        const staveNote = new VexFlow.StaveNote({
            keys: ["b/4"],
            duration: `${vexType}r`,
        });

        return {
            staveNote,
            tieTypes: [],
            beams: [],
            tuplet: null,
            key: null,
        };
    }

    if (!type) return null;
    const vexType = xmlPitchToVexStringPitch(type);
    const tieTypes = [...note.querySelectorAll("tie")].map((tie) => tie.getAttribute("type"));
    const beams = [...note.querySelectorAll("beam")].map((beam) => ({
        number: Number(beam.getAttribute("number") ?? 1),
        value: beam.textContent.trim().toLowerCase(),
    }));

    const dots = note.querySelectorAll("dot").length;

    const timeModification = note.querySelector("time-modification");

    const tuplet = timeModification
        ? {
              actualNotes: Number(timeModification.querySelector("actual-notes")?.textContent),
              normalNotes: Number(timeModification.querySelector("normal-notes")?.textContent),
              type: note.querySelector("notations tuplet")?.getAttribute("type") ?? null,
          }
        : null;

    let durationString = vexType;

    // IMPORTANT:
    // Tell VexFlow that this note/rest is a tuplet note.
    // Otherwise Voice sees normal eighths and tick counts become wrong.
    if (tuplet) {
        // durationString += "t";
    }

    // -----------------------------
    // Rest
    // -----------------------------
    if (note.querySelector("rest")) {
        const staveNote = new VexFlow.StaveNote({
            keys: ["b/4"],
            duration: `${durationString}r`,
        });

        for (let i = 0; i < dots; i++) {
            staveNote.addDotToAll();
        }

        return {
            staveNote,
            tieTypes,
            beams,
            tuplet,
            key: null,
        };
    }

    // -----------------------------
    // Pitched note
    // -----------------------------
    const pitch = note.querySelector("pitch");
    if (!pitch) return null;

    const step = pitch.querySelector("step")?.textContent;
    const octave = pitch.querySelector("octave")?.textContent;
    const alter = Number(pitch.querySelector("alter")?.textContent ?? 0);

    let accidental = "";

    switch (alter) {
        case 2:
            accidental = "##";
            break;
        case 1:
            accidental = "#";
            break;
        case -1:
            accidental = "b";
            break;
        case -2:
            accidental = "bb";
            break;
    }

    const key = `${step.toLowerCase()}${accidental}/${octave}`;

    const stem = note.querySelector("stem")?.textContent;

    let stemDirection;
    if (stem === "up") {
        stemDirection = VexFlow.Stem.UP;
    } else if (stem === "down") {
        stemDirection = VexFlow.Stem.DOWN;
    }

    const staveNote = new VexFlow.StaveNote({
        keys: [key],
        duration: durationString,
        stemDirection,
    });

    if (accidental && !tieTypes.includes("stop")) {
        staveNote.addModifier(new VexFlow.Accidental(accidental), 0);
    }

    for (let i = 0; i < dots; i++) {
        staveNote.addDotToAll();
    }

    renderNotehead(staveNote, note, type);

    return {
        staveNote,
        tieTypes,
        beams,
        tuplet,
        key,
    };
}

// ╭─────────────────────────────────────╮
// │              Parse Xml              │
// ╰─────────────────────────────────────╯
function createFactory() {
    return new VexFlow.Factory({
        renderer: {
            elementId: "vexscore",
            width: 1000,
            height: 700,
        },
    });
}

// ─────────────────────────────────────
async function loadMusicXml(url) {
    const response = await fetch(url);

    if (!response.ok) {
        throw new Error(`HTTP ${response.status}`);
    }

    const xmlText = await response.text();
    const parser = new DOMParser();
    const xmlDoc = parser.parseFromString(xmlText, "application/xml");
    if (xmlDoc.querySelector("parsererror")) {
        throw new Error("Invalid XML");
    }

    return xmlDoc;
}

// ─────────────────────────────────────
function addTechniqueText(staveNote, text) {
    if (!text) return;

    const annotation = new VexFlow.Annotation(text)
        .setVerticalJustification(VexFlow.Annotation.VerticalJustify.TOP)
        .setFont("Bravura", 12, "");

    staveNote.addModifier(annotation, 0);
}

// ─────────────────────────────────────
function renderCustomNoteheads(voice) {
    for (const tick of voice.tickables) {
        if (tick.myNotehead && tick._noteHeads?.length) {
            tick._noteHeads[0]._text = tick.myNotehead;
        }
    }
}

// ─────────────────────────────────────
function processTie(result, openTies, ties) {
    if (!result.key) return;
    if (result.tieTypes.includes("start")) {
        openTies.set(result.key, result.staveNote);
    }
    if (result.tieTypes.includes("stop")) {
        const firstNote = openTies.get(result.key);
        if (!firstNote) return;
        ties.push(
            new VexFlow.StaveTie({
                firstNote,
                lastNote: result.staveNote,
                firstIndexes: [0],
                lastIndexes: [0],
            }),
        );
        openTies.delete(result.key);
    }
}

// ─────────────────────────────────────
function processBeam(result, openBeams, measureBeams) {
    for (const beam of result.beams) {
        if (beam.number !== 1) {
            continue;
        }

        switch (beam.value) {
            case "begin":
                openBeams.set(1, [result.staveNote]);
                break;

            case "continue":
                if (openBeams.has(1)) {
                    openBeams.get(1).push(result.staveNote);
                }
                break;

            case "end":
                if (!openBeams.has(1)) {
                    break;
                }

                const notes = openBeams.get(1);

                notes.push(result.staveNote);

                if (notes.length > 1) {
                    const beam = new VexFlow.Beam(notes);

                    beam.renderOptions.flat_beams = true;

                    measureBeams.push(beam);
                }

                openBeams.delete(1);
                break;
        }
    }
}

// ─────────────────────────────────────
function processTuplet(result, openTuplets, tuplets) {
    if (!result.tuplet) return;
    const id = `${result.tuplet.actualNotes}:${result.tuplet.normalNotes}`;
    switch (result.tuplet.type) {
        case "start":
            openTuplets.set(id, {
                notes: [result.staveNote],
                actualNotes: result.tuplet.actualNotes,
                normalNotes: result.tuplet.normalNotes,
            });
            break;
        case null:
            if (openTuplets.has(id)) {
                openTuplets.get(id).notes.push(result.staveNote);
            }
            break;
        case "stop":
            if (openTuplets.has(id)) {
                const group = openTuplets.get(id);

                group.notes.push(result.staveNote);

                tuplets.push(
                    new VexFlow.Tuplet(group.notes, {
                        numNotes: group.actualNotes,
                        notesOccupied: group.normalNotes,
                    }),
                );

                openTuplets.delete(id);
            }
            break;
    }
}

// ─────────────────────────────────────
function parseMeasureNotes(measure, openTies, ties, scoreEventMap) {
    const vexNotes = [];
    const measureBeams = [];
    const tuplets = [];

    const openBeams = new Map();
    const openTuplets = new Map();
    const children = [...measure.children];

    for (const note of measure.querySelectorAll("note")) {
        if (note.querySelector("chord")) {
            continue;
        }

        const index = children.indexOf(note);
        const previousNode = index > 0 ? children[index - 1] : null;

        const techniqueText =
            previousNode?.localName === "direction" ? (previousNode.querySelector("words")?.textContent ?? "") : "";

        const result = xmlToVexflowNote(note);

        if (!result) {
            continue;
        }

        vexNotes.push(result.staveNote);
        addTechniqueText(result.staveNote, techniqueText);
        processTie(result, openTies, ties);
        processBeam(result, openBeams, measureBeams);
        processTuplet(result, openTuplets, tuplets);

        if (result.key !== null) {
            eventCount++;
        }

        if (result.tieTypes.length > 0) {
            if (result.tieTypes[0] != "start") {
                eventCount--;
            }
        }

        if (result.key !== null) {
            const currentEvent = scoreEventMap.get(eventCount);
            if (currentEvent === undefined) {
                scoreEventMap.set(eventCount, [result.staveNote]);
            } else {
                currentEvent.push(result.staveNote);
            }
        }
    }

    return {
        vexNotes,
        measureBeams,
        tuplets,
    };
}

// ─────────────────────────────────────
function createVoice(factory, beats, beatType, notes) {
    const voice = factory.Voice({
        time: `${beats}/${beatType}`,
    });

    voice.addTickables(notes);

    return voice;
}

// ─────────────────────────────────────
function createMeasureSystem(factory, x, y, width) {
    return factory.System({
        x,
        y,
        width,
        spaceBetweenStaves: 20,
    });
}

// ─────────────────────────────────────
function configureStave(system, voice, currentX, currentMeasure, previousMeasure) {
    const stave = system.addStave({
        voices: [voice],
    });

    if (currentX === 10) {
        stave.addClef("treble");
    }

    if (currentMeasure !== previousMeasure) {
        stave.addTimeSignature(currentMeasure);
    }

    return stave;
}

// ─────────────────────────────────────
function drawBeams(context, beams) {
    for (const beam of beams) {
        beam.setContext(context);
        beam.draw();
    }
}

// ─────────────────────────────────────
function drawTies(context, ties) {
    for (const tie of ties) {
        tie.setContext(context);
        tie.draw();
    }
}

// ─────────────────────────────────────
function drawTuplets(context, tuplets) {
    for (const tuplet of tuplets) {
        tuplet.setContext(context);
        tuplet.draw();
    }
}

// ─────────────────────────────────────
async function loadXml(scorefile) {
    let svgs = document.querySelectorAll("svg");
    if (svgs.length > 0) {
        for (let svg of svgs) {
            console.log(svg);
            svg.remove();
        }
    }

    const xmlDoc = await loadMusicXml(scorefile);
    const factory = createFactory();
    const ties = [];
    const allBeams = [];
    const allTuplets = [];
    const openTies = new Map();
    eventCount = 0;

    let previousMeasure = "";

    let currentX = 10;
    let currentY = 40;

    const measureWidth = 300;
    const maxLineX = 1000;
    const scoreEventMap = new Map();

    for (const part of xmlDoc.querySelectorAll("part")) {
        let beats = 4;
        let beatType = 4;

        for (const measure of part.querySelectorAll("measure")) {
            const time = measure.querySelector("time");
            if (time) {
                beats = Number(time.querySelector("beats")?.textContent);
                beatType = Number(time.querySelector("beat-type")?.textContent);
            }
            const { vexNotes, measureBeams, tuplets } = parseMeasureNotes(measure, openTies, ties, scoreEventMap);
            if (vexNotes.length === 0) {
                continue;
            }

            allTuplets.push(...tuplets);
            allBeams.push(...measureBeams);

            const voice = createVoice(factory, beats, beatType, vexNotes);
            const system = createMeasureSystem(factory, currentX, currentY, measureWidth);
            const currentMeasure = `${beats}/${beatType}`;
            configureStave(system, voice, currentX, currentMeasure, previousMeasure);
            renderCustomNoteheads(voice);
            previousMeasure = currentMeasure;
            currentX += measureWidth;

            if (currentX + measureWidth > maxLineX) {
                currentX = 10;
                currentY += 140;
            }
        }
    }

    factory.draw();
    const context = factory.getContext();
    drawBeams(context, allBeams);
    drawTuplets(context, allTuplets);
    drawTies(context, ties);
    return scoreEventMap;
}

//╭─────────────────────────────────────╮
//│                Init                 │
//╰─────────────────────────────────────╯
window.loadXml = loadXml;
