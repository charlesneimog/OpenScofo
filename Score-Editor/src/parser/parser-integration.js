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

    this.OScofoQuery = new Query(scoreScofo, this.OpenScofoHighlightQuery);
    this.LuaQuery = new Query(luaParser, this.LuaStringQuery);
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
    let needParse = this.runFormatterAfterParse(newTree.rootNode);
    if (needParse) {
        this.handleCodeChange(_, changes);
        return;
    }

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

export function runFormatterAfterParse(_) {
    const source = this.codeEditor.getValue();
    const edits = [];
    const shouldFormatStructure = _.hasError;

    function ensureLineStart(instance, node, indentText) {
        const nodeStart = instance.codeEditor.indexFromPos({ line: node.startPosition.row, ch: node.startPosition.column });
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

    walk(_, (node) => {
        if (shouldFormatStructure && node.type === "EVENT") {
            ensureLineStart(this, node, "");
        }

        if (shouldFormatStructure && node.type === "ACTION") {
            ensureLineStart(this, node, "    ");
        }

        if (shouldFormatStructure && node.type === "timedAction") {
            if (!(node.parent && node.parent.type === "ACTION" && node.parent.startPosition.row === node.startPosition.row && node.parent.startPosition.column === node.startPosition.column)) {
                ensureLineStart(this, node, "    ");
            }
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

            const parentStart = this.codeEditor.indexFromPos({ line: luaParent.startPosition.row, ch: luaParent.startPosition.column });
            const parentEnd = this.codeEditor.indexFromPos({ line: luaParent.endPosition.row, ch: luaParent.endPosition.column });
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

export function getMissing(node, list) {
    for (let i = 0; i < node.namedChildCount; i++) {
        let child = node.namedChild(i);
        if (child.isMissing) {
            list.push(child);
        }
        this.getMissing(child);
    }
}

export function checkErrors(node) {
    node = node.rootNode;
    var lineWithErrors = [];
    var lineWithUnexpected = [];
    const errorContainer = document.getElementById("editor-console");
    errorContainer.innerHTML = "";

    function checkNode(node) {
        for (let i = 0; i < node.namedChildCount; i++) {
            let child = node.namedChild(i);

            if (child.hasError) {
                let message = "";
                lineWithErrors.push(child.startPosition.row + 1);

                if (child.type === "number" && child.text === "") {
                    message =
                        "Missing " +
                        child.parent.type +
                        " at line " +
                        (child.parent.endPosition.row + 1) +
                        ", column " +
                        (child.parent.endPosition.column + 1);
                } else if (child.type === "pitch" && child.text === "") {
                    message = `Missing pitch at line ${child.endPosition.row + 1}`;
                }

                if (message !== "") {
                    const errorElement = document.createElement("p");
                    errorElement.textContent = message;
                    errorContainer.appendChild(errorElement);
                }
            }

            let treeString = child.toString();
            if (treeString.startsWith("(UNEXPECTED") && !lineWithUnexpected.includes(child.startPosition.row + 1)) {
                lineWithUnexpected.push(child.startPosition.row + 1);
                const unexpectedElement = document.createElement("p");
                unexpectedElement.textContent = "UNEXPECTED keyword at line " + (child.startPosition.row + 1);
                errorContainer.appendChild(unexpectedElement);
            }

            checkNode(child);
        }
    }

    checkNode(node);
}