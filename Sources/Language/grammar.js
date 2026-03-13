/**
 * @file The OScofo language is a specialized scripting language for OpenScofo, an open-source score follower tailored for contemporary live-electronic music. Inspired by the research behind IRCAM's Antescofo, it leverages Lua for flexible scripting and Tree-sitter to parse EVENTs and musical actions. It is built to lower the barrier to entry and setup complexity of live-electronics by enabling interactive scores that integrate seamlessly across PureData (including web deployment via pd4web), Max, Python, and C/C++.
 * @author Charles K. Neimog <charlesneimog@outlook.com>
 * @license GPL3
 */

///
//

module.exports = grammar({
    name: "openscofo",
    word: ($) => $.identifier,
    rules: {
        source_file: ($) => repeat(choice($.CONFIG, $.LUA, $.EVENT)),

        //╭─────────────────────────────────────╮
        //│             IDENTIFIERS             │
        //╰─────────────────────────────────────╯
        identifier: (_) => token(/[a-zA-Z_][a-zA-Z0-9_-]*/),

        //╭─────────────────────────────────────╮
        //│                 LUA                 │
        //╰─────────────────────────────────────╯
        LUA: ($) => seq("LUA", "{", optional($.lua_body), "}"),
        lua_body: ($) => repeat1(choice(/[^{}`]+/, seq("{", optional($.lua_body), "}"), $.lua_comment)),
        lua_call: ($) => repeat1(choice(/[^()`]+/, seq("(", optional($.lua_call), ")"))),
        lua_comment: (_) => /--[^\n]*/,

        descriptor_list: ($) => seq("[", repeat1(field("descriptor", $.descriptor)), "]"),
        descriptor: (_) =>
            choice(
                "MFCC",
                "LOUDNESS",
                "RMS",
                "POWER",
                "CHROMA",
                "ZCR",
                "HFR",
                "CENTROID",
                "SPREAD",
                "FLATNESS",
                "FLUX",
                "IRREGULARITY",
                "HARMONICITY",
                "YIN",
            ),

        //╭─────────────────────────────────────╮
        //│                Config               │
        //╰─────────────────────────────────────╯
        CONFIG: ($) =>
            seq(field("key", $.config_key), field("value", choice($.number, $.identifier, $.path, $.descriptor_list))),

        config_key: (_) =>
            token(
                choice(
                    "BPM",
                    "PHASECOUPLING",
                    "SYNCSTRENGTH",
                    "TRANSPOSE",
                    "ENTROPY",
                    "PITCHTEMPLATESIGMA",
                    "FFTSIZE",
                    "HOPSIZE",
                    "DUMMY",
                    "ONNXMODEL",
                    "ONNXTENSOR",
                ),
            ),

        //╭─────────────────────────────────────╮
        //│                Events               │
        //╰─────────────────────────────────────╯
        EVENT: ($) =>
            seq(
                field(
                    "definition",
                    choice($.note_event, $.rest_event, $.chord_event, $.trill_event, $.tech_event, $.lua_event),
                ),
                repeat(field("action", $.action)),
            ),

        note_event: ($) =>
            seq(
                "NOTE",
                field("pitch", $.pitch),
                field("duration", $.number),
                optional(field("attribute", $.attribute)),
            ),

        rest_event: ($) => seq("REST", field("duration", $.number)),
        chord_event: ($) => seq("CHORD", field("pitches", $.pitch_group), field("duration", $.number)),
        trill_event: ($) => seq("TRILL", field("pitches", $.pitch_group), field("duration", $.number)),
        tech_event: ($) =>
            seq(
                "TECH",
                field("technique", $.identifier),
                optional(field("pitch", $.pitch)),
                field("duration", $.number),
            ),
        lua_event: ($) => seq("LUAEVENT", field("luacall", $.lua_call, field("duration", $.number))),

        // TODO: Add events attribute
        attribute: (_) => seq("@", field("type", choice("percussive", "other"))),

        // Pitch
        pitch_group: ($) => seq("(", repeat1(field("pitch", $.pitch)), ")"),
        pitch: ($) =>
            seq(
                field("pitch_name", $.pitch_name),
                optional(field("alteration", $.alteration)),
                field("octave", $.octave),
            ),
        pitch_name: (_) => token(/[A-Ga-g]/),
        alteration: (_) => token(choice("#", "b")),
        octave: (_) => token(/(1[0-2]|[0-9])/),

        //╭─────────────────────────────────────╮
        //│               ACTIONS               │
        //╰─────────────────────────────────────╯
        action: ($) =>
            seq(
                optional("ACTION"),
                optional(field("timing", $.delay)),
                field("command", $.exec),
                repeat(seq(",", field("command", $.exec))),
            ),

        delay: ($) => seq("delay", field("amount", $.number), field("unit", $.time_unit)),

        exec: ($) =>
            choice(
                seq("sendto", field("receiver", $.identifier), optional(field("args", $.pdargs))),
                seq("luacall", "(", field("lua", $.lua_call), ")"),
            ),

        pdargs: ($) => seq("[", repeat1($.pdarg), "]"),
        pdarg: ($) => choice($.identifier, $.number),

        //╭─────────────────────────────────────╮
        //│                ATOMS                │
        //╰─────────────────────────────────────╯
        number: (_) => token(choice(/-?[0-9]+/, /-?[0-9]+\.[0-9]+/)),
        path: (_) => token(choice(seq('"', /[^"\n]+/, '"'), /[a-zA-Z0-9_.\-\\/]+/)),

        comment: (_) => token(choice(seq("//", /(\\+(.|\r?\n)|[^\\\n])*/), seq("/*", /[^*]*\*+([^/*][^*]*\*+)*/, "/"))),
        time_unit: (_) => token(choice("tempo", "sec", "ms")),
    },
    extras: ($) => [/\s|\\\r?\n|'/, $.comment],
});
