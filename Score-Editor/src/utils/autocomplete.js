import { AUTOCOMPLETE_KEYWORDS } from "../config/constants.js";

export function myCustomAutocomplete(cm) {
    cm.showHint({
        hint: customHint,
    });
}

export function customHint(editor) {
    var cur = editor.getCursor();
    var token = editor.getTokenAt(cur);
    var word = token.string;

    var matches = AUTOCOMPLETE_KEYWORDS.filter(function (kw) {
        return kw.startsWith(word);
    });

    return {
        list: matches,
        from: { line: cur.line, ch: token.start },
        to: { line: cur.line, ch: token.end },
    };
}