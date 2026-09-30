import unittest

from pygments.token import Comment, Error, Keyword, Name, Number, String

from openscofo_lexer import OpenScofoLexer


class HighlightTests(unittest.TestCase):
    def tokens(self, source):
        tokens = list(OpenScofoLexer().get_tokens_unprocessed(source))
        self.assertEqual("".join(value for _, _, value in tokens), source)
        self.assertFalse(any(token in Error for _, token, _ in tokens))
        return [(token, value) for _, token, value in tokens if value.strip()]

    def test_action_words_are_receivers_and_arguments_in_context(self):
        tokens = self.tokens("sendto delay [sendto 1] delay 2 tempo sendto 3")
        self.assertEqual(tokens, [
            (Name.Function, "sendto"), (Name.Class, "delay"),
            (Number, "["), (Number, "sendto"), (Number, "1"), (Number, "]"),
            (Name.Function, "delay"), (Name.Variable, "2"),
            (Name.Variable, "tempo"), (Name.Function, "sendto"), (Keyword, "3"),
        ])

    def test_config_values_and_section_names(self):
        tokens = self.tokens(
            'BPM 1e+2 ONNXMODEL "model.onnx" ONNXDESCRIPTORS mfcc rms\n'
            'SECTION 1 SECTION "A" SECTION delay\nNOTE C#4 .5'
        )
        for value in ("1e+2", '"model.onnx"', "mfcc", "rms"):
            self.assertIn((Number, value), tokens)
        for value in ("1", '"A"', "delay"):
            self.assertIn((Name.Label, value), tokens)
        self.assertIn((Keyword.Namespace, "BPM"), tokens)
        self.assertIn((String, "C#4"), tokens)
        self.assertIn((Name.Variable, ".5"), tokens)

    def test_event_fields_and_technique_groups(self):
        tokens = self.tokens(
            "PTECH [tongue-ram, flutter] Bb3 .33\nUTECH delay 2\n"
            "REST 1\nCHORD (C4 E4) 4\nTRILL (F#4 G4) 5\nGLISS (C4 C5) 6"
        )
        for value in ("tongue-ram", "flutter", "delay", "Bb3", "C4", "E4", "F#4", "G4", "C5"):
            self.assertIn((String, value), tokens)
        for value in (".33", "2", "1"):
            self.assertIn((Name.Variable, value), tokens)
        # Group durations inherit the event capture in highlights.scm.
        for value in ("4", "5", "6"):
            self.assertIn((Keyword.Type, value), tokens)

    def test_comments_preserve_field_context(self):
        tokens = self.tokens("NOTE /** pitch */ C4 /* duration */ 1\n// NOTE C5 2\nREST 3")
        self.assertIn((Comment.Special, "/**"), tokens)
        self.assertIn((String, "C4"), tokens)
        self.assertIn((Name.Variable, "1"), tokens)
        self.assertIn((Name.Variable, "3"), tokens)
        self.assertNotIn((String, "C5"), tokens)

    def test_lua_nesting_strings_and_comments_do_not_end_body(self):
        tokens = self.tokens(
            'LUA { local t = {s = "}", long = [=[}]=]} -- }\n'
            'print(t) } NOTE C4 1 luacall(cue(")", nested(1))) REST 2'
        )
        self.assertIn((Keyword.Declaration, "local"), tokens)
        self.assertIn((String, "C4"), tokens)
        self.assertIn((Name.Function, "luacall"), tokens)
        self.assertIn((Name.Variable, "2"), tokens)
        self.assertEqual(sum(token == Keyword.Namespace and value == "}" for token, value in tokens), 1)


if __name__ == "__main__":
    unittest.main()
