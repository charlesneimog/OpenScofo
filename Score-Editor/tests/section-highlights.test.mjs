import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";

import { OPEN_SCOFO_HIGHLIGHT_QUERY } from "../src/config/constants.js";
import { Language, Parser, Query } from "../tree-sitter/web-tree-sitter.js";

const editorRoot = new URL("..", import.meta.url);
const parserWasm = new URL("tree-sitter/tree-sitter-openscofo.wasm", editorRoot).pathname;

test("website and language queries highlight section keywords and names", async () => {
    await Parser.init();
    const language = await Language.load(parserWasm);
    const parser = new Parser();
    parser.setLanguage(language);
    const tree = parser.parse('SECTION "Coda 2"');

    const highlightQueries = [
        ["website", OPEN_SCOFO_HIGHLIGHT_QUERY],
        ["language", readFileSync(new URL("../Sources/Language/highlights.scm", editorRoot), "utf8")],
        [
            "Neovim",
            readFileSync(new URL("../Sources/Language/nvim/queries/openscofo/highlights.scm", editorRoot), "utf8"),
        ],
    ];

    for (const [consumer, source] of highlightQueries) {
        const query = new Query(language, source);
        const sectionCaptures = query.captures(tree.rootNode).map(({ name, node }) => [name, node.text]);
        assert.deepEqual(
            sectionCaptures,
            [
                ["keyword.directive", "SECTION"],
                ["label", '"Coda 2"'],
            ],
            `${consumer} section captures`,
        );
        query.delete();
    }

    tree.delete();
    parser.delete();
});

test("VS Code grammar highlights section keywords and names", () => {
    const grammar = JSON.parse(
        readFileSync(new URL("../Sources/Language/vscode/syntaxes/openscofo.tmLanguage.json", editorRoot), "utf8"),
    );
    const sectionPattern = grammar.repository.sections.patterns[0];
    const match = new RegExp(sectionPattern.match).exec('SECTION "Coda 2"');

    assert.equal(match?.[1], "SECTION");
    assert.equal(match?.[2], '"Coda 2"');
});
