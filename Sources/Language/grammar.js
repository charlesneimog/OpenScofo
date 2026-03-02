/**
 * @file The OScofo language is a specialized scripting language for OpenScofo, an open-source score follower tailored for contemporary live-electronic music. Inspired by the research behind IRCAM's Antescofo, it leverages Lua for flexible scripting and Tree-sitter to parse events and musical actions. It is built to lower the barrier to entry and setup complexity of live-electronics by enabling interactive scores that integrate seamlessly across PureData (including web deployment via pd4web), Max, Python, and C/C++.
 * @author Charles K. Neimog <charlesneimog@outlook.com>
 * @license GPL3
 */

///
//

module.exports = grammar({
    name: "openscofo",
    precedence: "left",
    rules: {
        score: ($) => repeat($._statement),
        _statement: ($) => choice($.CONFIG, $.EVENT, $.LUA),
        EVENT: ($) => choice($.noteEvent, $.multiPitchEvent, $.restEvent, $.freeEvent),
        keyword: (_) => /[a-zA-Z_][a-zA-Z0-9_]*/,

        //╭─────────────────────────────────────╮
        //│                 LUA                 │
        //╰─────────────────────────────────────╯
        LUA: ($) => seq(alias("LUA", $.identifier), "{", optional($.lua_body), "}"),
        lua_body: ($) => repeat1(choice(/[^{}`]+/, seq("{", optional($.lua_body), "}"), $.lua_comment)),
        lua_call: ($) => field("lua_call", repeat1(choice(/[^()`]+/, seq("(", optional($.lua_call), ")")))),
        lua_function: ($) => $.keyword,
        lua_comment: (_) => /--[^\n]*/,

        //╭─────────────────────────────────────╮
        //│                Config               │
        //╰─────────────────────────────────────╯
        CONFIG: ($) => choice($.numberConfig, $.symbolConfig, $.pathConfig),

        numberConfigId: ($) =>
            field(
                "configId",
                choice(
                    alias("BPM", $.keyword),
                    alias("PHASECOUPLING", $.keyword),
                    alias("SYNCSTRENGTH", $.keyword),
                    alias("TRANSPOSE", $.keyword),
                    alias("ENTROPY", $.keyword),
                    alias("PITCHTEMPLATESIGMA", $.keyword),
                    alias("FFTSIZE", $.keyword),
                    alias("HOPSIZE", $.keyword),
                ),
            ),

        symbolConfigId: ($) => field("configId", choice(alias("DUMMY", $.keyword))),

        pathConfigId: ($) => field("configId", choice(alias("TIMBREMODEL", $.keyword))),

        numberConfig: ($) => seq($.numberConfigId, $.number),
        symbolConfig: ($) => seq($.symbolConfigId, $.symbol),
        pathConfig: ($) => seq($.pathConfigId, $.path),

        //╭─────────────────────────────────────╮
        //│                Events               │
        //╰─────────────────────────────────────╯
        EVENT: ($) => choice($.noteEvent, $.multiPitchEvent, $.restEvent, $.freeEvent),

        noteEvent: ($) => seq(alias("NOTE", $.keyword), $.pitch, $.duration, repeat($.ACTION)),
        multiPitchEvent: ($) =>
            seq(choice(alias("TRILL", $.keyword), alias("CHORD", $.keyword)), $.pitches, $.duration, repeat($.ACTION)),

        restEvent: ($) => seq(alias("REST", $.keyword), $.duration, repeat($.ACTION)),
        freeEvent: ($) => seq(alias("EVENT", $.keyword), $.eventId, $.duration, repeat1($.ACTION)),

        eventId: (_) => choice("timed", "internal"),

        // Pitch
        pitches: ($) => seq("(", repeat1($.pitch), ")"),
        pitch: ($) => choice($.noteName, $.midi),
        midi: (_) => token(/(0|[1-9][0-9]?|1[0-1][0-9]|12[0-7])/),
        noteName: ($) => seq($.pitchname, optional($.alteration), $.octave),

        // Pitch Atoms
        pitchname: (_) => /[A-Ga-g]/,
        alteration: (_) => choice("#", "b"),
        octave: (_) => /[0-9]|1[0-2]/,

        duration: ($) => $.number,

        //╭─────────────────────────────────────╮
        //│               ACTIONS               │
        //╰─────────────────────────────────────╯
        ACTION: ($) =>
            seq(
                optional(alias("ACTION", $.keyword)),
                optional(field("timedAction", $.timedAction)),
                field("exec", $.exec),
                repeat(seq(",", $.exec)),
            ),

        timedAction: ($) => seq(field("actionKey", "delay"), field("value", $.number), field("timeUnit", $.timeUnit)),

        actionKeyword: ($) => choice($.lua_function, $.keyword),

        exec: ($) =>
            choice(
                seq(field("keyword", "sendto"), field("receiver", $.receiver), optional(field("pdargs", $.pdargs))),
                seq(field("keyword", "luacall"), field("luabody", seq("(", $.lua_call, ")"))),
            ),

        pdargs: ($) => seq("[", repeat1($.pdarg), "]"),
        pdarg: ($) => field("pdarg", choice($.symbol, $.number)),
        receiver: ($) => $.keyword,
        actionKey: (_) => /[a-zA-Z][a-zA-Z0-9_]*/,

        //╭─────────────────────────────────────╮
        //│                ATOMS                │
        //╰─────────────────────────────────────╯
        number: (_) => token(choice(/-?[0-9]+/, /-?[0-9]+\.[0-9]+/)),
        symbol: (_) => token(/[a-zA-Z][a-zA-Z0-9_]*/),
        path: (_) => token(choice(seq('"', /[^"]+/, '"'), /[a-zA-Z0-9_.\-\\/]+/)),

        comment: (_) => token(choice(seq("//", /(\\+(.|\r?\n)|[^\\\n])*/), seq("/*", /[^*]*\*+([^/*][^*]*\*+)*/, "/"))),
        timeUnit: (_) => field("timeUnit", choice("tempo", "sec", "ms")),
    },
    extras: ($) => [/\s|\\\r?\n/, $.comment],
});
