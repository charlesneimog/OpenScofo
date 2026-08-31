export const PARSER_OPEN_SCOFO_WASM = "tree-sitter/tree-sitter-openscofo.wasm";
export const PARSER_LUA_WASM = "tree-sitter/tree-sitter-lua.wasm";

export const AUTOCOMPLETE_KEYWORDS = [
    "NOTE",
    "REST",
    "TRILL",
    "CHORD",
    "TECH",
    "BPM",
    "SECTION",
    "SECTIONRESTRICT",
];

export const SUGGESTIONS = {
    NOTE: "NOTE",
    TRILL: "TRILL",
    CHORD: "CHORD",
    TECH: "TECH",
};

export const HIGHLIGHTS = {
    oscofo: {
        keyword: "color: var(--purple); font-weight: bold;",
        "keyword.directive": "color: var(--purple); font-weight: bold;",
        "type.builtin": "color: var(--purple); font-weight: bold;",
        string: "color: var(--green); font-weight: bold;",
        "variable.parameter": "color: var(--pink); font-weight: bold;",
        tempo: "color: var(--red); font-weight: bold;",
        function: "color: var(--blue); font-weight: normal;",
        type: "color: var(--yellow); font-weight: normal;",
        number: "color: var(--orange);",
        comment: "color: var(--comment); opacity: 0.9; font-style: italic;",
        "comment.documentation": "color: var(--comment); opacity: 0.9; font-style: italic;",
        error: "text-decoration: underline; text-decoration-style: wavy; text-decoration-color: var(--red);",
    },
    lua: {
        comment: "color: var(--comment); font-style: italic; opacity: var(--lua-opacity);",
        string: "color: var(--green); font-weight: 100; opacity: var(--lua-opacity);",
        "function.bracket": "color: var(--fg); font-weight: bold; opacity: var(--lua-opacity);",
        "function.call.lua": "color: var(--blue); font-weight: 100; opacity: var(--lua-opacity);",
        "function.name": "color: var(--blue); font-weight: 100; opacity: var(--lua-opacity);",
        "keyword.repeat": "color: var(--purple); font-weight: 100; opacity: var(--lua-opacity);",
        "keyword.conditional": "color: var(--purple); font-weight: 100; opacity: var(--lua-opacity);",
        "keyword.return": "color: var(--purple); opacity: var(--lua-opacity);",
        "keyword.function": "color: var(--purple); opacity: var(--lua-opacity);",
        variable: "color: var(--pink); opacity: var(--lua-opacity);",
    },
};

export const OPEN_SCOFO_HIGHLIGHT_QUERY = `
    ; ================================
    ; COMMENTS
    ; ================================
    (comment) @comment

    ((comment) @comment.documentation
     (#match? @comment "^/\\\\*\\\\*"))

    ; ================================
    ; STRUCTURE
    ; ================================
    (CONFIG) @keyword.directive
    (LUA) @keyword.directive
    (EVENT) @keyword

    ; Lua body is handled by the Lua parser/highlighter.
    (lua_body) @lua_body
    (lua_call) @lua_body
    (lua_comment) @comment

    ; ================================
    ; Config
    ; ================================
    (CONFIG
      value: (_) @number)

    ; ================================
    ; EVENTS (semantic types)
    ; ================================
    (note_event) @type.builtin
    (rest_event) @type.builtin
    (ptech_event) @type.builtin
    (utech_event) @type.builtin
    (chord_event) @type.builtin
    (trill_event) @type.builtin
    (lua_event) @type.builtin

    ; ================================
    ; Musical data
    ; ================================
    (pitch) @string
    [
      (ptech_event
        technique: (identifier) @string)

      (ptech_event
        techniques: (technique_group
          technique: (identifier) @string))

      (utech_event
        technique: (identifier) @string)

      (utech_event
        techniques: (technique_group
          technique: (identifier) @string))
    ]

    [
      (note_event
        duration: (number) @tempo)
      (rest_event
        duration: (number) @tempo)
      (chord_event
        duration: (number) @tempo)
      (trill_event
        duration: (number) @tempo)
      (ptech_event
        duration: (number) @tempo)
      (utech_event
        duration: (number) @tempo)
    ]

    ; ================================
    ; ACTIONS (runtime layer)
    ; ================================
    (exec receiver: (identifier) @type)
    (exec args: (pdargs) @number)

    ; actions keywords
    (exec "sendto" @function)
    (exec "luacall" @function)
    (delay "delay" @function)

    ; delay
    (delay amount: (number) @tempo)
    (delay unit: (time_unit) @tempo)

    ; Errors
    (ERROR) @error
`;
