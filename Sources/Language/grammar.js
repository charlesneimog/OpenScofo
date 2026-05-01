/**
 * @file The OpenScofo language is a specialized scripting language for OpenScofo, an open-source score follower tailored for contemporary live-electronic music. Inspired by the research behind IRCAM's Antescofo, it leverages Lua for flexible scripting and Tree-sitter to parse EVENTs and musical actions. It is built to lower the barrier to entry and setup complexity of live-electronics by enabling interactive scores that integrate seamlessly across PureData (including web deployment via pd4web), Max, Python, and C/C++.
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

        descriptor_list: ($) => repeat1(field("descriptor", $.descriptor)),
        descriptor: (_) =>
            choice(
                "mfcc",
                "logmel",
                "loudness",
                "rms",
                "power",
                "chroma",
                "zcr",
                "hfr",
                "centroid",
                "spread",
                "flatness",
                "flux",
                "irregularity",
                "kurtosis",
                "harmonicity",
                "yin",
            ),

        onset_function: (_) => choice("pow", "pd", "wpd", "sf", "cd", "rcd", "hfc", "mkl"),

        //╭─────────────────────────────────────╮
        //│                Config               │
        //╰─────────────────────────────────────╯
        CONFIG: ($) =>
            seq(
                field("key", $.config_key),
                field(
                    "value",
                    choice(
                        prec(2, $.descriptor_list),
                        $.number,
                        $.identifier,
                        $.path,
                        $.descriptor_list,
                        $.onset_function,
                    ),
                ),
            ),

        config_key: (_) =>
            token(
                choice(
                    "BPM",

                    // Audio
                    "FFTSIZE",
                    "HOPSIZE",

                    // Tempo
                    "PHASECOUPLING",
                    "SYNCSTRENGTH",
                    "TRANSPOSE",

                    // Listening model
                    "PITCHTEMPLATESIGMA",
                    "ONNXMODEL",
                    "ONNXDESCRIPTORS",
                    "ONSETFUNCTION",
                ),
            ),

        //╭─────────────────────────────────────╮
        //│                Events               │
        //╰─────────────────────────────────────╯
        EVENT: ($) =>
            seq(
                field(
                    "definition",
                    choice(
                        $.note_event,
                        $.rest_event,
                        $.chord_event,
                        $.trill_event,
                        $.ptech_event,
                        $.utech_event,
                        $.lua_event,
                    ),
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
        ptech_event: ($) =>
            seq(
                "PTECH",
                field("technique", $.identifier),
                field("pitch", $.pitch),
                field("duration", $.number),
                //
            ),
        utech_event: ($) =>
            seq(
                "UTECH",
                field("technique", $.identifier),
                field("duration", $.number),
                //
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
