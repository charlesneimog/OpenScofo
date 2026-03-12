import {
    HIGHLIGHTS,
    OPEN_SCOFO_HIGHLIGHT_QUERY,
    PARSER_LUA_WASM,
    PARSER_OPEN_SCOFO_WASM,
} from "../config/constants.js";
import { myCustomAutocomplete } from "../utils/autocomplete.js";
import { debounce } from "../utils/debounce.js";
import {
    addGhostText,
    fetchTextFile,
    getHighlights,
    highlightTechFallback,
    luaChildOfLuaBody,
    luaHighlight,
    luaIndentBody,
    runTreeQuery,
    showCodeSuggestions,
} from "../parser/highlighting.js";
import {
    bootstrapParserInitialization,
    checkErrors,
    getMissing,
    handleCodeChange,
    initParser,
    luaFormatter,
    runFormatterAfterParse,
    runFormatterBeforeParser,
    treeEditForEditorChange,
} from "../parser/parser-integration.js";
import {
    downloadScore,
    generateOpenScofoScore,
    getPitch,
    loadTxtScore,
    parseScore,
    uploadScore,
    uploadWav,
} from "../services/score-service.js";

export class OpenScofoOnlineEditor {
    constructor() {
        bootstrapParserInitialization(this);

        this.debug = false;

        this.ScofoParser = null;
        this.parserOpenScofoWasm = PARSER_OPEN_SCOFO_WASM;
        this.parserLuaWasm = PARSER_LUA_WASM;

        this.htmlFormatterCode = "";
        this.linesWithErrors = {};
        this.musicxmlScore = {};

        this.htmlFormatter = {};
        this.bigTextWarning = false;

        this.codeInput = document.getElementById("code-input");
        this.codeEditor = CodeMirror.fromTextArea(this.codeInput, {
            lineNumbers: true,
            showCursorWhenSelecting: true,
            placeholder: "// Editor your score here",
            extraKeys: {
                "Ctrl-Z": function (cm) {
                    console.log("Ctrl-Z");
                    cm.undo();
                },
                "Ctrl-Y": function (cm) {
                    cm.redo();
                },
                "Cmd-Z": function (cm) {
                    cm.undo();
                },
                "Cmd-Shift-Z": function (cm) {
                    cm.redo();
                },
                "Ctrl-Space": function (cm) {
                    myCustomAutocomplete(cm);
                },
                "Ctrl-;": function (cm) {
                    const selection = cm.getSelection();
                    if (selection) {
                        cm.replaceSelection(`/* ${selection} */`);
                    } else {
                        const cursor = cm.getCursor();
                        const line = cm.getLine(cursor.line);
                        cm.replaceRange(
                            "// " + line,
                            { line: cursor.line, ch: 0 },
                            { line: cursor.line, ch: line.length },
                        );
                    }
                },
                "Ctrl-S": function (_) {
                    scofoEditor.downloadScore();
                },
            },
        });

        this.codeEditor.on("cursorActivity", () => {
            const cursorPos = this.codeEditor.getCursor();
            this.showCodeSuggestions(this.codeEditor, cursorPos);
        });

        this.codeContainer = document.getElementById("code-container");

        this.loadState();
        this.saveStateOnChange = this.debounce(this.saveState, 2000);
        this.runTreeQueryOnChange = this.debounce(this.runTreeQuery, 100);
        this.runTreeQueryOnViewportChange = this.debounce(this.runTreeQuery, 50);

        this.handleCodeChange = this.handleCodeChange.bind(this);
        this.treeEditForEditorChange = this.treeEditForEditorChange.bind(this);
        this.codeEditor.on("changes", this.handleCodeChange);
        this.codeEditor.on("viewportChange", (cm, from, to) => {
            this.runTreeQueryOnViewportChange(cm, from, to);
        });

        this.highlights = HIGHLIGHTS;

        this.OpenScofoHighlightQuery = OPEN_SCOFO_HIGHLIGHT_QUERY;

        this.LuaStringQuery = null;
        this.fetchTextFile("highlight/lua.scm");

        const downloadButton = document.getElementById("download-score");
        const uploadButtom = document.getElementById("upload-score");
        const loadButtom = document.getElementById("load-score");
        const audioButtom = document.getElementById("upload-wav");

        if (downloadButton && uploadButtom && loadButtom && audioButtom) {
            downloadButton.addEventListener("click", () => this.downloadScore());
            uploadButtom.addEventListener("click", () => this.uploadScore());
            loadButtom.addEventListener("click", () => this.loadTxtScore());
            audioButtom.addEventListener("click", () => this.uploadWav());
        } else {
            alert("Buttons not found");
        }
    }

    loadState() {
        const sourceCode = localStorage.getItem("sourceCode");
        if (sourceCode) {
            this.codeInput.value = sourceCode;
        } else {
            this.codeInput.value = `// Edit you score here`;
        }
    }

    saveState() {
        localStorage.setItem("sourceCode", this.codeEditor.getValue());
    }

    debounce(func, wait, immediate) {
        return debounce(func, wait, immediate);
    }
}

Object.assign(OpenScofoOnlineEditor.prototype, {
    showCodeSuggestions,
    fetchTextFile,
    addGhostText,
    getHighlights,
    handleCodeChange,
    treeEditForEditorChange,
    luaFormatter,
    luaHighlight,
    luaChildOfLuaBody,
    luaIndentBody,
    highlightTechFallback,
    runFormatterBeforeParser,
    runFormatterAfterParse,
    runTreeQuery,
    initParser,
    getMissing,
    checkErrors,
    downloadScore,
    getPitch,
    generateOpenScofoScore,
    parseScore,
    uploadWav,
    uploadScore,
    loadTxtScore,
});