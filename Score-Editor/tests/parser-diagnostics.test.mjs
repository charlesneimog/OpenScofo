import assert from "node:assert/strict";
import test from "node:test";
import { checkErrors, getMissing } from "../src/parser/parser-integration.js";
import { Language, Parser } from "../tree-sitter/web-tree-sitter.js";

await Parser.init();
const language = await Language.load(new URL("../tree-sitter/tree-sitter-openscofo.wasm", import.meta.url).pathname);

function diagnostics(source) {
    const parser = new Parser();
    parser.setLanguage(language);
    const tree = parser.parse(source);
    const consoleElement = {
        style: {},
        messages: [],
        set innerHTML(_) { this.messages = []; },
        appendChild(element) { this.messages.push(element.textContent); },
    };
    const previousDocument = globalThis.document;
    globalThis.document = {
        getElementById: () => consoleElement,
        createElement: () => ({ style: {} }),
    };
    try {
        checkErrors(tree);
        const messages = [...consoleElement.messages];
        const missing = getMissing(tree.rootNode).map((node) => node.type);
        // A subsequent valid parse clears stale errors.
        const valid = parser.parse("BPM 60\nNOTE C4 1");
        try {
            checkErrors(valid);
            assert.deepEqual(consoleElement.messages, []);
        } finally {
            valid.delete();
        }
        return { messages, missing };
    } finally {
        globalThis.document = previousDocument;
        tree.delete();
        parser.delete();
    }
}

test("reports missing duration and anonymous closing brace with precise positions", () => {
    assert.deepEqual(diagnostics("BPM 60\nNOTE C4\nLUA { print(1)"), {
        messages: ["Missing duration at line 2, column 8", 'Missing "}" at line 3, column 15'],
        missing: ["number", "}"],
    });
});

test("reports missing action punctuation", () => {
    assert.deepEqual(diagnostics("BPM 60\nNOTE C4 1 sendto foo [1").messages, [
        'Missing "]" at line 2, column 24',
    ]);
});

test("reports all unexpected regions without duplicating nested errors", () => {
    assert.deepEqual(diagnostics("BPM 60\n@\nNOTE C4 1\n$\nREST 1").messages, [
        'Unexpected text "@" at line 2, column 1',
        'Unexpected text "$" at line 4, column 1',
    ]);
});

test("reports unrecognized identifiers and bounds multiline previews", () => {
    assert.deepEqual(diagnostics("garbage").messages, ['Unexpected text "garbage" at line 1, column 1']);
    const { messages } = diagnostics("invalid ".repeat(30) + "\ninvalid");
    assert.equal(messages.length, 1);
    assert.match(messages[0], /…" at line 1, column 1$/);
    assert(!messages[0].includes("\n"));
});
