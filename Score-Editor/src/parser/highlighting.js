import { SUGGESTIONS } from "../config/constants.js";

export function showCodeSuggestions(cm, startPos) {
    const line = startPos.line;
    const ch = startPos.ch;
    const currentWord = cm.getLine(line).slice(0, ch).split(/\s+/).pop();

    if (SUGGESTIONS[currentWord]) {
        this.addGhostText(line, SUGGESTIONS[currentWord]);
    }
}

export async function fetchTextFile(filePath) {
    try {
        const response = await fetch(filePath);
        if (!response.ok) {
            throw new Error(`Error fetching file: ${response.statusText}`);
        }
        const text = await response.text();
        this.LuaStringQuery = text;
    } catch (error) {
        console.error("Error:", error);
    }
}

export function addGhostText(_, __) {
    return;
}

export function getHighlights(language, capture) {
    return this.highlights[language][capture] || "";
}

export function luaHighlight(luaNode, luaPositionStart) {
    let luaTxt = luaNode.text;
    let luaTree = this.LuaParser.parse(luaTxt);

    if (luaTree.hasError) {
        return;
    }

    const captures = this.LuaQuery.captures(luaTree.rootNode, {
        startPosition: { row: 0, column: 0 },
        endPosition: { row: Infinity, column: Infinity },
    });

    for (const { name, node } of captures) {
        const { startPosition, endPosition } = node;
        let luaStart = {
            row: startPosition.row + luaPositionStart.row,
            column: startPosition.column + (startPosition.row === 0 ? luaPositionStart.column : 0),
        };
        let luaEnd = {
            row: endPosition.row + luaPositionStart.row,
            column: endPosition.column + (endPosition.row === 0 ? luaPositionStart.column : 0),
        };

        this.codeEditor.markText(
            { line: luaStart.row, ch: luaStart.column },
            { line: luaEnd.row, ch: luaEnd.column },
            {
                inclusiveLeft: true,
                inclusiveRight: true,
                css: this.getHighlights("lua", name),
            },
        );
    }
}

export function luaChildOfLuaBody(node) {
    let parent = node.parent;
    while (parent) {
        if (parent.type === "lua_body") {
            return true;
        }
        parent = parent.parent;
    }
    return false;
}

export function luaIndentBody(node) {
    if (node.parent.hasError) {
        return false;
    }
    let luaTxt = node.text;
    if (luaTxt.trim(" ") === "") {
        return false;
    }
    let lines = luaTxt.split("\n");

    const cursorPos = this.codeEditor.getCursor();
    const cursorLine = cursorPos.line;
    const cursorCh = cursorPos.ch;
    let needIndent = false;
    let indentedLines = lines.map((line, index) => {
        let modifiedLine = line;
        let addedSpaces = 0;
        if (!(line.startsWith("    ") || line.startsWith("\t")) && index !== lines.length - 1) {
            needIndent = true;
            const match = line.match(/^(\s*)/);
            let neededSpaces = 4;
            if (match) {
                neededSpaces = 4 - (match[1].length % 4);
            }
            modifiedLine = " ".repeat(neededSpaces) + line;
            addedSpaces = neededSpaces;
        }
        if (index === cursorLine && addedSpaces > 0 && cursorCh <= line.length) {
            cursorPos.ch += addedSpaces;
        }
        return modifiedLine;
    });

    if (needIndent) {
        const indentedText = indentedLines.join("\n");
        const startPosition = { line: node.startPosition.row, ch: node.startPosition.column };
        const endPosition = { line: node.endPosition.row, ch: node.endPosition.column };
        this.codeEditor.replaceRange(indentedText, startPosition, endPosition);
        this.codeEditor.setCursor(cursorPos);
        return true;
    }

    return false;
}

export function runTreeQuery(_, startRow, endRow) {
    if (endRow == null) {
        const viewport = this.codeEditor.getViewport();
        startRow = viewport.from;
        endRow = viewport.to;
    }

    this.codeEditor.operation(() => {
        const marks = this.codeEditor.getAllMarks();
        marks.forEach((m) => m.clear());

        if (this.tree && this.OScofoQuery) {
            const captures = this.OScofoQuery.captures(this.tree.rootNode, {
                startPosition: { row: startRow, column: 0 },
                endPosition: { row: endRow, column: 0 },
            });

            let lastNodeId;
            for (const { name, node } of captures) {
                if (node.id === lastNodeId) continue;
                lastNodeId = node.id;
                const { startPosition, endPosition } = node;
                if (!node.hasError) {
                    if (name === "lua_body" && !this.luaChildOfLuaBody(node)) {
                        this.luaHighlight(node, startPosition);
                    } else {
                        this.codeEditor.markText(
                            { line: startPosition.row, ch: startPosition.column },
                            { line: endPosition.row, ch: endPosition.column },
                            {
                                inclusiveLeft: true,
                                inclusiveRight: true,
                                css: this.getHighlights("oscofo", name),
                            },
                        );
                    }
                }
            }
        }
    });
}