export const PARSER_OPEN_SCOFO_WASM = "tree-sitter/tree-sitter-openscofo.wasm";
export const PARSER_LUA_WASM = "tree-sitter/tree-sitter-lua.wasm";

export const AUTOCOMPLETE_KEYWORDS = ["NOTE", "REST", "TRILL", "CHORD", "BPM"];

export const SUGGESTIONS = {
    NOTE: "NOTE",
    TRILL: "TRILL",
    CHORD: "CHORD",
};

export const HIGHLIGHTS = {
    oscofo: {
        keyword: "color: var(--fg); font-weight: normal;",
        eventKeyword: "color: var(--red); font-weight: bold;",
        luaKeyword: "color: var(--purple); font-weight: bold;",
        pitch: "color: var(--blue); font-weight: bold;",
        duration: "color: var(--yellow); font-weight: bold;",
        comment: "color: var(--comment); opacity: 0.4; font-style: italic;",
        config: "color: var(--purple); font-weight: bold;",
        error: "text-decoration: underline; text-decoration-style: wavy; text-decoration-color: var(--red);",
        action: "color: var(--red); font-weight: normal;",
        actionKey: "color: var(--red); font-weight: normal;",
        timeUnit: "color: var(--purple); font-weight: normal;",
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
        ; LUA keyword
        (LUA (identifier) @luaKeyword)

    ; Primary event keywords
        ((keyword) @eventKeyword
            (#match? @eventKeyword "^(NOTE|REST|TRILL|CHORD|EVENT)$"))

        ; Primary config keywords
        ((keyword) @config
            (#match? @config "^(BPM|PHASECOUPLING|SYNCSTRENGTH|TRANSPOSE|ENTROPY|PITCHTEMPLATESIGMA|FFTSIZE|HOPSIZE|DUMMY|TIMBREMODEL)$"))

        ; Secondary action keywords
        (ACTION (keyword) @action)
        "sendto" @action
        "luacall" @action
        "delay" @actionKey
        ((keyword) @action
            (#match? @action "^(sendto|luacall)$"))
        ((keyword) @actionKey
            (#match? @actionKey "^delay$"))

        ; Generic identifiers
        (identifier) @keyword

        ; Configuration nodes
    (numberConfigId) @config
    (symbolConfigId) @config
    (pathConfigId) @config

    ; Musical Attributes
    (pitch) @pitch
    (duration (number) @duration)

    ; Actions and Executions
    (timedAction) @timedAction
    (timedAction value: (number) @duration)
    (timeUnit) @timeUnit
    (exec) @exec
    (receiver) @receiver

    ; Literals and Comments
    (number) @number
    (pdargs) @pdargs
    (pdarg) @pdarg
    (comment) @comment
    (lua_comment) @comment
    (lua_body) @lua_body
    (lua_call) @lua_body

    ; Errors
    (ERROR) @error
`;