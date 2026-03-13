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
        return text;
    } catch (error) {
        console.error("Error:", error);
        return null;
    }
}

export function addGhostText(_, __) {
    return;
}

export function getHighlights(language, capture) {
    return this.highlights[language][capture] || "";
}

export function luaHighlight(luaNode, luaPositionStart) {
    if (!this.LuaParser || !this.LuaQuery) {
        return;
    }

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

export function highlightTechFallback(startRow, endRow) {
    const techPattern = /(\bTECH\b)\s+([a-zA-Z][a-zA-Z0-9_-]*)(?:\s+(-?[0-9]+(?:\.[0-9]+)?))?/g;

    for (let line = startRow; line < endRow; line++) {
        const text = this.codeEditor.getLine(line);
        if (!text || !text.includes("TECH")) {
            continue;
        }

        techPattern.lastIndex = 0;
        let match;
        while ((match = techPattern.exec(text)) !== null) {
            const keywordStart = match.index + match[0].indexOf(match[1]);
            const techniqueStart = match.index + match[0].indexOf(match[2]);
            const durationStart = match[3] ? match.index + match[0].lastIndexOf(match[3]) : -1;

            this.codeEditor.markText(
                { line, ch: keywordStart },
                { line, ch: keywordStart + match[1].length },
                {
                    inclusiveLeft: true,
                    inclusiveRight: true,
                    css: this.getHighlights("oscofo", "eventKeyword"),
                },
            );

            this.codeEditor.markText(
                { line, ch: techniqueStart },
                { line, ch: techniqueStart + match[2].length },
                {
                    inclusiveLeft: true,
                    inclusiveRight: true,
                    css: this.getHighlights("oscofo", "techniqueId"),
                },
            );

            if (durationStart >= 0) {
                this.codeEditor.markText(
                    { line, ch: durationStart },
                    { line, ch: durationStart + match[3].length },
                    {
                        inclusiveLeft: true,
                        inclusiveRight: true,
                        css: this.getHighlights("oscofo", "duration"),
                    },
                );
            }
        }
    }
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

            for (const { name, node } of captures) {
                const { startPosition, endPosition } = node;
                const isErrorCapture = name === "error";
                if (node.hasError && !isErrorCapture) {
                    continue;
                }

                if (name === "lua_body" && !this.luaChildOfLuaBody(node)) {
                    this.luaHighlight(node, startPosition);
                } else {
                    const css = this.getHighlights("oscofo", name);
                    if (!css) {
                        continue;
                    }

                    this.codeEditor.markText(
                        { line: startPosition.row, ch: startPosition.column },
                        { line: endPosition.row, ch: endPosition.column },
                        {
                            inclusiveLeft: true,
                            inclusiveRight: true,
                            css,
                        },
                    );
                }
            }

            this.highlightTechFallback(startRow, endRow);

            const markMissingNodes = (node) => {
                for (let i = 0; i < node.namedChildCount; i++) {
                    const child = node.namedChild(i);
                    if (child.isMissing) {
                        const row = child.startPosition.row;
                        if (row >= startRow && row < endRow) {
                            const lineText = this.codeEditor.getLine(row) || "";
                            if (lineText.length > 0) {
                                let fromCh = child.startPosition.column;
                                let toCh = fromCh + 1;

                                if (fromCh >= lineText.length) {
                                    fromCh = Math.max(0, lineText.length - 1);
                                    toCh = lineText.length;
                                }

                                if (toCh > fromCh) {
                                    this.codeEditor.markText(
                                        { line: row, ch: fromCh },
                                        { line: row, ch: toCh },
                                        {
                                            inclusiveLeft: true,
                                            inclusiveRight: true,
                                            css: this.getHighlights("oscofo", "error"),
                                        },
                                    );
                                }
                            }
                        }
                    }
                    markMissingNodes(child);
                }
            };

            markMissingNodes(this.tree.rootNode);
        }
    });
}
