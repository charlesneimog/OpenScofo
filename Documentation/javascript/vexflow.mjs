import {
    Stave,
    StaveNote,
    Beam,
    Formatter,
    Renderer,
    Accidental,
    Tremolo,
    Barline,
    StaveText,
} from "https://cdn.jsdelivr.net/npm/vexflow@5.0.0/+esm";

// ─────────────────────────────────────
function parseNoteString(noteStr) {
    noteStr = noteStr.trim();

    let keysPart = noteStr;
    let duration = "q"; // default duration
    if (noteStr.includes("/")) {
        const parts = noteStr.split("/");
        keysPart = parts[0].trim();
        duration = parts[1] ? parts[1].trim() : "q";
        // sanitize duration
        const validDurations = ["w", "h", "q", "8", "16", "32", "64"];
        if (!validDurations.includes(duration.replace(/[rd]/, ""))) duration = "q";
    }

    if (!keysPart) keysPart = "b4"; // fallback

    if (keysPart.toLowerCase().endsWith("r")) {
        return new StaveNote({ keys: ["b/4"], duration: duration + "r" });
    }

    let rawKeys = [];
    if (keysPart.startsWith("(") && keysPart.endsWith(")")) {
        rawKeys = keysPart.slice(1, -1).split(" ");
    } else {
        rawKeys = [keysPart];
    }

    const keys = [];
    const accidentals = [];
    rawKeys.forEach((k) => {
        if (!k) return;
        const octave = k.slice(-1);
        let note = k.slice(0, -1);
        let acc = null;
        if (note.endsWith("#") || note.endsWith("b")) {
            acc = note.slice(-1);
            note = note.slice(0, -1);
        }
        keys.push(note + "/" + octave);
        accidentals.push(acc);
    });

    const staveNote = new StaveNote({ keys, duration });
    accidentals.forEach((acc, i) => {
        if (acc) staveNote.addModifier(new Accidental(acc), i);
    });

    staveNote._myOwnSymbols = noteStr;
    return staveNote;
}

// ─────────────────────────────────────
export function renderScores() {
    const scores = document.querySelectorAll("score");
    const scoreWidth = 400;

    scores.forEach((scoreEl) => {
        const raw = scoreEl.getAttribute("notes");
        const timeSig = scoreEl.getAttribute("timeSig") || "4/4";
        const id = scoreEl.getAttribute("id");
        const el = document.getElementById(id);

        scoreEl.innerHTML = "";
        el.style.boxShadow = "0 2px 6px rgba(0,0,0,0.25)";
        el.style.borderRadius = "6px";

        const renderer = new Renderer(el, Renderer.Backends.SVG);
        renderer.resize(scoreWidth * 1.3, 150);
        const context = renderer.getContext();
        context.scale(1.3, 1.3);

        const tokens = raw.split(",");
        const notesArray = [];

        var tremolo = false;

        var tremoloNotes = [];
        for (let i = 0; i < tokens.length; i++) {
            let t = tokens[i].trim();
            const tremMatch = t.match(/:trem(\d+)/);
            if (tremMatch && i + 1 < tokens.length) {
                t = t.replace(/:trem\d+/, "").trim();
                const note1 = parseNoteString(t);
                note1.children[0].hide = true;

                const note2 = parseNoteString(tokens[i + 1].trim().replace(/:trem\d+/, ""));
                note2.children[0].hide = true;

                notesArray.push(note1, note2);
                tremoloNotes.push([note1, note2]);
                tremolo = true;
                i++;
            } else {
                notesArray.push(parseNoteString(t));
            }
        }

        const beams = Beam.generateBeams(notesArray);
        if (tremolo) {
            for (const note of tremoloNotes) {
                if (note.length == 2) {
                    note[0]._noteHeads[0]._text = "\ue0a2";
                    note[1]._noteHeads[0]._text = "\ue0a2";
                }
            }
        }

        for (const note of notesArray) {
            if (note._myOwnSymbols.includes(":diamond")) {
                note._noteHeads[0]._text = "\ue0C6";
            }
            if (note._myOwnSymbols.includes(":x")) {
                note._noteHeads[0]._text = "\ue0A8";
            }
        }

        const stave = new Stave(4, 0, scoreWidth);
        stave.addClef("treble").addTimeSignature(timeSig);
        stave.setContext(context).draw();

        Formatter.FormatAndDraw(context, stave, notesArray);
        beams.forEach((b) => b.setContext(context).draw());

        // After creating and drawing stave:
        if (tremolo) {
            const glyph = "\uE222"; // SMuFL codepoint
            context.save();
            context.setFont("Bravura", 24, "normal"); // font family + size
            context.fillText(glyph, 170, 60); // absolute SVG coords
            context.restore();
        }

        // add final bar
        stave.setEndBarType(Barline.type.END);
        stave.setContext(context).draw();
    });
}

//╭─────────────────────────────────────╮
//│            THEME VALUES             │
//╰─────────────────────────────────────╯
var THEME_FG = null;
const themeQuery = window.matchMedia("(prefers-color-scheme: dark)");

// ─────────────────────────────────────
function updateTheme(isDark) {
    THEME_FG = isDark ? "white" : "black";
}

// ─────────────────────────────────────
function handleThemeChange(e) {
    updateTheme(e.matches);
    renderScores();
}

themeQuery.addEventListener("change", handleThemeChange);

window.onload = function () {
    updateTheme(themeQuery.matches);
    renderScores();
};
