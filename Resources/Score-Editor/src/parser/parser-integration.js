import { Parser, Language, Query } from "../../tree-sitter/web-tree-sitter.js";

export function bootstrapParserInitialization(instance) {
    Parser.init().then(() => {
        instance.initParser().then(() => {
            instance.handleCodeChange();
        });
    });
}

export async function initParser() {
    await Parser.init();

    this.ScofoParser = new Parser();
    const scoreScofo = await Language.load(this.parserOpenScofoWasm);
    this.ScofoParser.setLanguage(scoreScofo);

    this.LuaParser = new Parser();
    const luaParser = await Language.load(this.parserLuaWasm);
    this.LuaParser.setLanguage(luaParser);

    this.OpenScofoQuery = new Query(scoreScofo, this.OpenScofoHighlightQuery);

    if (!this.LuaStringQuery) {
        await this.fetchTextFile("highlight/lua.scm");
    }

    if (this.LuaStringQuery) {
        this.LuaQuery = new Query(luaParser, this.LuaStringQuery);
    } else {
        this.LuaQuery = null;
        console.warn("Lua highlight query was not loaded; Lua syntax highlighting is disabled.");
    }
}

export function handleCodeChange(_, changes) {
    const newText = this.codeEditor.getValue();
    const edits = this.tree && changes && changes.map(this.treeEditForEditorChange);
    if (edits) {
        for (const edit of edits) {
            this.tree.edit(edit);
        }
    }
    const newTree = this.ScofoParser.parse(newText, this.tree);

    this.checkErrors(newTree);
    if (this.tree) this.tree.delete();
    this.tree = newTree;

    if (this.debug) {
        console.log(newTree.rootNode.toString());
    }

    this.runTreeQueryOnChange();
    this.saveStateOnChange();
}

export function treeEditForEditorChange(change) {
    const oldLineCount = change.removed.length;
    const newLineCount = change.text.length;
    const lastLineLength = change.text[newLineCount - 1].length;
    const startPosition = { row: change.from.line, column: change.from.ch };
    const oldEndPosition = { row: change.to.line, column: change.to.ch };
    const newEndPosition = {
        row: startPosition.row + newLineCount - 1,
        column: newLineCount === 1 ? startPosition.column + lastLineLength : lastLineLength,
    };
    const startIndex = this.codeEditor.indexFromPos(change.from);
    let newEndIndex = startIndex + newLineCount - 1;
    let oldEndIndex = startIndex + oldLineCount - 1;
    for (let i = 0; i < newLineCount; i++) newEndIndex += change.text[i].length;
    for (let i = 0; i < oldLineCount; i++) oldEndIndex += change.removed[i].length;
    return {
        startIndex,
        oldEndIndex,
        newEndIndex,
        startPosition,
        oldEndPosition,
        newEndPosition,
    };
}

export function luaFormatter(_, __) {
    console.log("here");
    return;
}

export function runFormatterBeforeParser() {
    return;
}

export function formatScore() {
    if (!this.ScofoParser) {
        return;
    }

    const tree = this.ScofoParser.parse(this.codeEditor.getValue());
    const changed = this.codeEditor.operation(() => this.runFormatterAfterParse(tree.rootNode));
    tree.delete();

    if (!changed) {
        this.handleCodeChange();
    }
}

export function runFormatterAfterParse(rootNode) {
    const source = this.codeEditor.getValue();
    const edits = [];
    const shouldFormatStructure = rootNode.hasError;
    let isInsideSection = false;

    function ensureLineStart(instance, node, indentText) {
        const nodeStart = instance.codeEditor.indexFromPos({
            line: node.startPosition.row,
            ch: node.startPosition.column,
        });
        const lineStart = instance.codeEditor.indexFromPos({ line: node.startPosition.row, ch: 0 });

        let wsStart = nodeStart;
        while (wsStart > 0) {
            const ch = source[wsStart - 1];
            if (ch === " " || ch === "\t") {
                wsStart--;
                continue;
            }
            break;
        }

        const hasNewlineBefore = wsStart > 0 && source[wsStart - 1] === "\n";
        if (!hasNewlineBefore) {
            if (wsStart > 0) {
                edits.push({
                    start: wsStart,
                    end: nodeStart,
                    text: `\n${indentText}`,
                });
            }
            return;
        }

        const currentIndent = source.slice(lineStart, nodeStart);
        if (currentIndent !== indentText) {
            edits.push({
                start: lineStart,
                end: nodeStart,
                text: indentText,
            });
        }
    }

    function walk(node, visit) {
        visit(node);
        for (let i = 0; i < node.namedChildCount; i++) {
            walk(node.namedChild(i), visit);
        }
    }

    walk(rootNode, (node) => {
        if (node.type === "SECTION") {
            isInsideSection = true;
        }

        if (node.type === "EVENT" && (shouldFormatStructure || isInsideSection)) {
            ensureLineStart(this, node, isInsideSection ? "\t" : "");
        }

        if (node.type === "action") {
            ensureLineStart(this, node, isInsideSection ? "\t\t" : "\t");
        }

        if (node.type === "lua_body") {
            const luaParent = node.parent;
            if (!luaParent || luaParent.type !== "LUA" || luaParent.hasError) {
                return;
            }

            const isInlineLua = luaParent.startPosition.row === luaParent.endPosition.row;
            if (!isInlineLua) {
                return;
            }

            const parentStart = this.codeEditor.indexFromPos({
                line: luaParent.startPosition.row,
                ch: luaParent.startPosition.column,
            });
            const parentEnd = this.codeEditor.indexFromPos({
                line: luaParent.endPosition.row,
                ch: luaParent.endPosition.column,
            });
            const bodyText = node.text.trim();
            if (bodyText === "") {
                return;
            }

            const parentText = source.slice(parentStart, parentEnd);
            const openBraceIndex = parentText.indexOf("{");
            if (openBraceIndex === -1) {
                return;
            }

            const header = parentText.slice(0, openBraceIndex).trimEnd();
            const formattedLua = `${header} {\n    ${bodyText}\n}`;

            edits.push({
                start: parentStart,
                end: parentEnd,
                text: formattedLua,
            });
        }
    });

    if (edits.length === 0) {
        return false;
    }

    edits
        .sort((a, b) => b.start - a.start)
        .forEach((edit) => {
            this.codeEditor.replaceRange(
                edit.text,
                this.codeEditor.posFromIndex(edit.start),
                this.codeEditor.posFromIndex(edit.end),
            );
        });

    return true;
}

// Missing punctuation is anonymous, so diagnostics must visit all children.
function walkSyntaxNodes(node, visit, insideError = false) {
    visit(node, insideError);
    for (let i = 0; i < node.childCount; i++) {
        walkSyntaxNodes(node.child(i), visit, insideError || node.isError);
    }
}

function missingLabel(node) {
    if (node.type === "number" && node.parent) {
        for (const field of ["duration", "amount"]) {
            if (node.parent.childForFieldName(field)?.id === node.id) {
                return field;
            }
        }
    }
    return node.isNamed ? node.type : JSON.stringify(node.type);
}

export function getMissing(node, list = []) {
    walkSyntaxNodes(node, (child) => {
        if (child.isMissing) list.push(child);
    });
    return list;
}

export function checkErrors(tree) {
    const errorContainer = document.getElementById("editor-console");
    errorContainer.style.color = "var(--red)";
    errorContainer.innerHTML = "";

    walkSyntaxNodes(tree.rootNode, (node, insideError) => {
        let message;
        const position = `at line ${node.startPosition.row + 1}, column ${node.startPosition.column + 1}`;
        if (node.isMissing) {
            message = `Missing ${missingLabel(node)} ${position}`;
        } else if (node.isError && !insideError) {
            const text = node.text.replace(/\s+/g, " ").trim();
            const preview = text.length > 80 ? `${text.slice(0, 80)}…` : text;
            message = preview ? `Unexpected text ${JSON.stringify(preview)} ${position}` : `Syntax error ${position}`;
        } else {
            return;
        }

        const element = document.createElement("p");
        element.style.color = "var(--red)";
        element.textContent = message;
        errorContainer.appendChild(element);
    });
}
