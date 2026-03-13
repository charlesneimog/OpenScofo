export const PARSER_OPEN_SCOFO_WASM = "tree-sitter/tree-sitter-openscofo.wasm";
export const PARSER_LUA_WASM = "tree-sitter/tree-sitter-lua.wasm";

export const AUTOCOMPLETE_KEYWORDS = ["NOTE", "REST", "TRILL", "CHORD", "TECH", "BPM"];

export const SUGGESTIONS = {
    NOTE: "NOTE",
    TRILL: "TRILL",
    CHORD: "CHORD",
    TECH: "TECH",
};

export const HIGHLIGHTS = {
    oscofo: {
        keyword: "color: var(--fg); font-weight: normal;",
        eventKeyword: "color: var(--red); font-weight: bold;",
        luaKeyword: "color: var(--purple); font-weight: bold;",
        pitch: "color: var(--blue); font-weight: bold;",
        duration: "color: var(--yellow); font-weight: bold;",
        techniqueId: "color: var(--green); font-style: italic; font-weight: bold;",
        attributeId: "color: var(--cyan); font-weight: bold;",
        comment: "color: var(--comment); opacity: 0.4; font-style: italic;",
        config: "color: var(--purple); font-weight: bold;",
        error: "text-decoration: underline; text-decoration-style: wavy; text-decoration-color: var(--red);",
        action: "color: var(--red); font-weight: normal;",
        actionKey: "color: var(--red); font-weight: normal;",
        timeUnit: "color: var(--yellow); font-weight: bold;",
        number: "color: var(--fg);",
        exec: "color: var(--fg); font-weight: normal;",
        receiver: "color: var(--green); font-weight: normal;",
        pdargs: "color: var(--fg); font-weight: normal;",
        pdarg: "color: var(--cyan); font-style: italic; font-weight: normal;",
        timedAction: "color: var(--fg); font-weight: normal;",
    },
    lua: {
        comment: "color: var(--comment); opacity: 0.4; font-style: italic; opacity: var(--lua-opacity);",
        string: "color: var(--green); font-style: italic; font-weight: 100; opacity: var(--lua-opacity);",
        "function.bracket": "color: var(--fg); font-weight: bold; opacity: var(--lua-opacity);",
        "function.call.lua": "color: var(--cyan); font-weight: 100; opacity: var(--lua-opacity);",
        "function.name": "color: var(--cyan); font-weight: 100; opacity: var(--lua-opacity);",
        "keyword.repeat": "color: var(--cyan); font-weight: 100; opacity: var(--lua-opacity);",
        "keyword.conditional": "color: var(--cyan); font-weight: 100; opacity: var(--lua-opacity);",
        "keyword.return": "color: var(--purple); opacity: var(--lua-opacity);",
        "keyword.function": "color: var(--red); opacity: var(--lua-opacity);",
        variable: "color: var(--blue); opacity: var(--lua-opacity);",
    },
};

export const OPEN_SCOFO_HIGHLIGHT_QUERY = `
    ; LUA
    (LUA) @luaKeyword
    (lua_body) @lua_body
    (lua_call) @lua_body
    (lua_comment) @comment

    ; Event keywords
    (note_event) @eventKeyword
    (rest_event) @eventKeyword
    (chord_event) @eventKeyword
    (trill_event) @eventKeyword
    (tech_event) @eventKeyword

    ; Config
    (config_key) @config
    (number) @number

    ; Musical data
    (pitch) @pitch
    (pitch_name) @pitch
    (note_event duration: (number) @duration)
    (rest_event duration: (number) @duration)
    (chord_event duration: (number) @duration)
    (trill_event duration: (number) @duration)
    (tech_event duration: (number) @duration)
    (tech_event technique: (identifier) @techniqueId)
    (attribute) @attributeId

    ; Actions
    (action) @timedAction
    (delay) @actionKey
    (delay amount: (number) @duration)
    (delay unit: (time_unit) @timeUnit)
    (exec) @exec
    (exec "sendto" @actionKey)
    (exec "luacall" @actionKey)
    (exec receiver: (identifier) @receiver)
    (pdargs) @pdargs
    (pdarg (identifier) @pdarg)
    (pdarg (number) @pdarg)

    ; Numbers and comments
    (comment) @comment

    ; Errors
    (ERROR) @error
`;
