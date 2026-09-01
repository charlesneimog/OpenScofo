import assert from "node:assert/strict";
import test from "node:test";

import { runFormatterAfterParse } from "../src/parser/parser-integration.js";
import { Language, Parser } from "../tree-sitter/web-tree-sitter.js";

const parserWasm = new URL("../tree-sitter/tree-sitter-openscofo.wasm", import.meta.url).pathname;

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
