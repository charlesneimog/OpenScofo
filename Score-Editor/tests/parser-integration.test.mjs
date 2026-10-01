import assert from "node:assert/strict";
import test from "node:test";
import { readFileSync } from "node:fs";

import { runFormatterAfterParse } from "../src/parser/parser-integration.js";
import { OPEN_SCOFO_HIGHLIGHT_QUERY } from "../src/config/constants.js";
import { Language, Parser, Query } from "../tree-sitter/web-tree-sitter.js";

const parserWasm = new URL("../tree-sitter/tree-sitter-openscofo.wasm", import.meta.url).pathname;

test("accepts and highlights MIDI pitches while keeping durations separate", async () => {
    await Parser.init();
    const language = await Language.load(parserWasm);
    const parser = new Parser();
    parser.setLanguage(language);
    const fixture = readFileSync(new URL("../../Tests/02-score/midi-pitches.scofo", import.meta.url), "utf8");
    const tree = parser.parse(fixture);
    assert.equal(tree.rootNode.hasError, false, tree.rootNode.toString());
    tree.delete();

    const note = parser.parse("NOTE 60.5 0.25");
    const definition = note.rootNode.namedChild(0).childForFieldName("definition");
    assert.equal(definition.childForFieldName("pitch").childForFieldName("midi").text, "60.5");
    assert.equal(definition.childForFieldName("duration").text, "0.25");
    const query = new Query(language, OPEN_SCOFO_HIGHLIGHT_QUERY);
    const captures = query.captures(note.rootNode);
    assert(captures.some(({ name, node }) => name === "string" && node.text === "60.5"));
    assert(captures.some(({ name, node }) => name === "tempo" && node.text === "0.25"));
    query.delete();
    note.delete();

    const incomplete = parser.parse("NOTE 60.5");
    assert.equal(incomplete.rootNode.hasError, true, "A MIDI pitch still requires a duration");
    incomplete.delete();
    parser.delete();
});

function positionAt(source, index) {
    const lines = source.slice(0, index).split("\n");
    return { line: lines.length - 1, ch: lines.at(-1).length };
}

function createEditor(initialValue) {
    let value = initialValue;

    return {
        getValue() {
            return value;
        },
        indexFromPos(position) {
            const lines = value.split("\n");
            let index = 0;
            for (let line = 0; line < position.line; line++) {
                index += lines[line].length + 1;
            }
            return index + position.ch;
        },
        posFromIndex(index) {
            return positionAt(value, index);
        },
        replaceRange(text, from, to) {
            const start = this.indexFromPos(from);
            const end = this.indexFromPos(to);
            value = value.slice(0, start) + text + value.slice(end);
        },
    };
}

test("indents score events inside sections with a tab", async () => {
    await Parser.init();
    const parser = new Parser();
    const language = await Language.load(parserWasm);
    parser.setLanguage(language);

    const source = ["NOTE B4 1", "SECTION A", "NOTE C#5 2", "sendto a [1 2 3 4]", "NOTE D6 2"].join("\n");
    const editor = createEditor(source);
    const tree = parser.parse(source);

    const changed = runFormatterAfterParse.call({ codeEditor: editor }, tree.rootNode);

    assert.equal(changed, true);
    assert.equal(
        editor.getValue(),
        ["NOTE B4 1", "SECTION A", "\tNOTE C#5 2", "\tsendto a [1 2 3 4]", "\tNOTE D6 2"].join("\n"),
    );

    tree.delete();

    const formattedTree = parser.parse(editor.getValue());
    assert.equal(runFormatterAfterParse.call({ codeEditor: editor }, formattedTree.rootNode), false);

    formattedTree.delete();
    parser.delete();
});
