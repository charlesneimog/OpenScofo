from pygments.lexer import ExtendedRegexLexer, default, include, words
from pygments.lexers import LuaLexer
from pygments.token import Comment, Keyword, Name, Number, Punctuation, String, Text


NUMBER = r"-?(?:[0-9]+(?:\.[0-9]+)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
IDENTIFIER = r"[a-zA-Z_][a-zA-Z0-9_-]*"
PITCH = r"[A-Ga-g][#b]?(?:1[0-2]|[0-9])(?![\w-])"
DESCRIPTORS = (
    "mfcc", "logmel", "loudness", "rms", "power", "chroma", "zcr", "hfr",
    "centroid", "spread", "flatness", "flux", "irregularity", "kurtosis",
    "harmonicity", "yin",
)
CONFIG_KEYS = (
    "SR", "FFTSIZE", "HOPSIZE", "TUNINGA4", "PHASECOUPLING", "SYNCSTRENGTH",
    "BPM", "TIMETOLERANCE", "PITCHTEMPLATESIGMA", "PITCHTEMPLATEHARMONICS",
    "TRANSPOSE", "ONNXMODEL", "TIMBREMODEL", "ONNXDESCRIPTORS", "ONSETFUNCTION",
    "ONAUDIOSTATECHANGE", "MFCCMELS", "MFCCCOUNT", "MEDSPAN", "DBTHRESHOLD",
    "SPECTRALROLLOFFCUTOFF", "YINTHRESHOLD", "YINMINFREQUENCY", "YINMAXFREQUENCY",
    "CHROMASIZE", "CHROMACENTEROCTAVE", "CHROMAOCTAVEWIDTH", "ZCRCENTER",
    "ZCRPAD", "ZCRZEROPOS", "ZCRTHRESHOLD", "SECTIONRESTRICT",
)


def lua_body(lexer, match, ctx):
    """Delegate balanced Lua bodies/calls without counting braces in strings."""
    opening = match.group()
    closing = "}" if opening == "{" else ")"
    delimiter_token = Keyword.Namespace if opening == "{" else Punctuation
    yield match.start(), delimiter_token, opening
    start = match.end()
    depth = 1
    for offset, token, value in LuaLexer().get_tokens_unprocessed(ctx.text[start:ctx.end]):
        if token in Punctuation:
            for index, char in enumerate(value):
                if char == opening:
                    depth += 1
                elif char == closing:
                    depth -= 1
                    if depth == 0:
                        if index:
                            yield start + offset, token, value[:index]
                        yield start + offset + index, delimiter_token, char
                        ctx.pos = start + offset + index + 1
                        return
        yield start + offset, token, value
    ctx.pos = ctx.end


class OpenScofoLexer(ExtendedRegexLexer):
    """Token categories mirror Sources/Language/highlights.scm.

    directive -> Keyword.Namespace; type.builtin -> Keyword.Type;
    string -> String; variable.parameter -> Name.Variable;
    label -> Name.Label; function -> Name.Function; type -> Name.Class.
    Lua bodies follow the language's Tree-sitter injection query.
    """

    name = "OpenScofo"
    aliases = ["openscofo", "scofo"]
    filenames = ["*.scofo", "*.openscofo"]

    tokens = {
        "root": [
            include("extras"),
            (words(CONFIG_KEYS, suffix=r"\b(?!-)"), Keyword.Namespace, "config-value"),
            (r"SECTION\b", Keyword.Namespace, "section-name"),
            (r"LUA\b", Keyword.Namespace, "lua-block"),
            (r"NOTE\b", Keyword.Type, ("duration", "pitch")),
            (r"REST\b", Keyword.Type, "duration"),
            (r"PTECH\b", Keyword.Type, ("duration", "pitch", "technique")),
            (r"UTECH\b", Keyword.Type, ("duration", "technique")),
            (r"(?:CHORD|TRILL|GLISS)\b", Keyword.Type, ("group-duration", "pitch-group")),
            (r"LUAEVENT\b", Keyword.Type, ("group-duration", "lua-event")),
            (r"ACTION\b", Keyword),
            (r"sendto\b", Name.Function, "receiver"),
            (r"luacall\b", Name.Function, "lua-call"),
            (r"delay\b", Name.Function, ("delay-unit", "duration")),
            (r"\[", Number, "arguments"),
            (r"@(?:percussive|other)\b", Keyword.Type),
            (NUMBER, Number),
            (r'"[^"\n]*"', String),
            (r"[()\[\]{},]", Punctuation),
            (IDENTIFIER, Name),
            (r".", Text),
        ],
        "extras": [
            (r"\s+|\\\r?\n|'", Text.Whitespace),
            (r"//(?:\\(?:\r?\n|.)|[^\\\n])*", Comment.Single),
            (r"/\*\*(?!/)", Comment.Special, "doc-comment"),
            (r"/\*", Comment.Multiline, "comment"),
        ],
        "comment": [
            (r"\*/", Comment.Multiline, "#pop"),
            (r"[^*]+|\*", Comment.Multiline),
        ],
        "doc-comment": [
            (r"\*/", Comment.Special, "#pop"),
            (r"[^*]+|\*", Comment.Special),
        ],
        "config-value": [
            include("extras"),
            # The query captures every config value as @number, including paths
            # and descriptor lists, rather than highlighting their words globally.
            (words(DESCRIPTORS, suffix=r"\b(?!-)"), Number, ("#pop", "descriptors")),
            (NUMBER, Number, "#pop"),
            (r'"[^"\n]*"|[a-zA-Z0-9_.\-\\/]+', Number, "#pop"),
            default("#pop"),
        ],
        "descriptors": [
            include("extras"),
            (words(DESCRIPTORS, suffix=r"\b(?!-)"), Number),
            default("#pop"),
        ],
        "section-name": [
            include("extras"),
            (r'"[^"\n]*"|' + IDENTIFIER + "|" + NUMBER, Name.Label, "#pop"),
            default("#pop"),
        ],
        "pitch": [
            include("extras"),
            (PITCH, String, "#pop"),
            default("#pop"),
        ],
        "duration": [
            include("extras"),
            (NUMBER, Name.Variable, "#pop"),
            default("#pop"),
        ],
        "pitch-group": [
            include("extras"),
            (r"\(", Keyword.Type),
            (PITCH, String),
            (r"\)", Keyword.Type, "#pop"),
            default("#pop"),
        ],
        "group-duration": [
            include("extras"),
            # These durations inherit @type.builtin in highlights.scm.
            (NUMBER, Keyword.Type, "#pop"),
            default("#pop"),
        ],
        "technique": [
            include("extras"),
            (r"\[", Keyword.Type, ("#pop", "technique-group")),
            (IDENTIFIER, String, "#pop"),
            default("#pop"),
        ],
        "technique-group": [
            include("extras"),
            (IDENTIFIER, String),
            (r",", Keyword.Type),
            (r"\]", Keyword.Type, "#pop"),
            default("#pop"),
        ],
        "receiver": [
            include("extras"),
            (IDENTIFIER, Name.Class, "#pop"),
            (NUMBER, Keyword, "#pop"),
            default("#pop"),
        ],
        "arguments": [
            include("extras"),
            (r"\]", Number, "#pop"),
            (IDENTIFIER + "|" + NUMBER, Number),
            default("#pop"),
        ],
        "delay-unit": [
            include("extras"),
            (r"(?:tempo|sec|ms)\b", Name.Variable, "#pop"),
            default("#pop"),
        ],
        "lua-block": [
            include("extras"),
            (r"\{", lua_body, "#pop"),
            default("#pop"),
        ],
        "lua-call": [
            include("extras"),
            (r"\(", lua_body, "#pop"),
            default("#pop"),
        ],
        "lua-event": [
            include("extras"),
            (IDENTIFIER, Name.Function),
            (r"\(", lua_body, "#pop"),
            default("#pop"),
        ],
    }
