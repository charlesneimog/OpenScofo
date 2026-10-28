import contextlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET


spec = importlib.util.spec_from_file_location(
    "svg_score_attr", Path(__file__).with_name("svg-score-attr.py")
)
score = importlib.util.module_from_spec(spec)
spec.loader.exec_module(score)


def xml_note(step="E", octave=4, extra="", attributes=""):
    return (
        f"<note {attributes}>{extra}<pitch><step>{step}</step>"
        f"<octave>{octave}</octave></pitch></note>"
    )


def musicxml(notes, attributes=""):
    return ET.fromstring(
        '<score-partwise><part id="P1"><measure number="1">'
        f"{attributes}{notes}</measure></part></score-partwise>"
    )


def svg(heads, bottoms=(400,), spacing=80):
    root = ET.Element("svg", xmlns="http://www.w3.org/2000/svg")
    for bottom in bottoms:
        for line in range(5):
            y = bottom - line * spacing
            ET.SubElement(root, "polyline", {
                "class": "StaffLines", "points": f"0,{y} 10000,{y}",
            })
    for name, x, y in heads:
        ET.SubElement(root, "path", {
            "class": "Note", "id": name,
            "d": f"M{x - 10},{y - 10} h20 v20 h-20 z",
        })
    return root


def match(xml, drawing):
    notes = score.get_musicxml_notes(xml)
    score.assign_event_numbers(notes)
    return score.match_noteheads(
        [note for note in notes if note["visible"]],
        score.sort_svg_notes(score.get_svg_notes(drawing)),
    )


class NoteheadTests(unittest.TestCase):
    def test_ledger_notes_do_not_split_a_staff(self):
        drawing = svg([
            ("high", 300, -120), ("low", 100, 400),
            ("diamond-position", 200, 480), ("next-system", 100, 2400),
        ], bottoms=(400, 2400))
        xml = musicxml(
            xml_note() + xml_note("C") + xml_note("D", 6) + xml_note()
        )
        self.assertEqual(
            [head["element"].get("id") for _, head in match(xml, drawing)],
            ["low", "diamond-position", "high", "next-system"],
        )

    def test_fallback_uses_gaps_between_notes(self):
        notes = [{"x": 30, "y": 100}, {"x": 10, "y": 450}, {"x": 20, "y": 700}]
        self.assertEqual([n["x"] for n in score.sort_svg_notes(notes)], [10, 20, 30])

    def test_real_musescore_glyph_centers(self):
        # The diamond starts at its right tip; the oval starts near its bottom.
        drawing = svg([], bottoms=(2659.56,))
        ET.SubElement(drawing, "path", {
            "class": "Note", "id": "oval",
            "d": "M1442.26,2705.9 C1475.03,2705.9 1515.4,2675.45 1515.4,2644 "
                 "C1515.4,2625.47 1500.51,2613.22 1479.32,2613.22 "
                 "C1438.94,2613.22 1406.51,2643.34 1406.51,2675.12 "
                 "C1406.51,2694.31 1422.4,2705.9 1442.26,2705.9",
        })
        ET.SubElement(drawing, "path", {
            "class": "Note", "id": "diamond",
            "d": "M6958.81,2743.56 C6959.13,2743.24 6959.46,2742.56 6959.46,2741.91 "
                 "C6959.46,2741.24 6959.13,2740.58 6958.81,2740.25 "
                 "L6920.07,2701.52 C6919.74,2701.19 6919.09,2700.86 6918.09,2700.86 "
                 "C6917.43,2700.86 6916.43,2701.19 6916.1,2701.52 L6877.37,2740.25 "
                 "C6876.71,2740.91 6876.71,2741.24 6876.71,2741.91 "
                 "C6876.71,2742.56 6877.04,2743.24 6877.71,2743.89 L6916.43,2782.95 "
                 "C6917.09,2783.61 6917.43,2783.61 6918.09,2783.61 "
                 "C6918.42,2783.61 6919.09,2783.61 6919.74,2782.95 L6958.81,2743.56",
        })
        heads = score.get_svg_notes(drawing)
        self.assertAlmostEqual(heads[0]["y"], 2659.56)
        self.assertAlmostEqual(heads[1]["y"], 2742.235)
        self.assertAlmostEqual(heads[1]["x"], 6918.085)

    def test_compact_relative_paths_and_exponents(self):
        path = ET.Element("path", d="m1e2,2e2 20-10 h10 v30 l-30,0z")
        self.assertEqual(score.extract_svg_position(path), (100, 200))
        self.assertEqual(score.svg_path_points(path)[-1], (100, 220))

    def test_offset_chord_heads_match_pitch_in_any_xml_order(self):
        xml = musicxml(xml_note("F", 3) + xml_note(extra="<chord/>"))
        drawing = svg([("cross", 100, 400), ("diamond", 107, 640)])
        pairs = match(xml, drawing)
        self.assertEqual(
            [(note["pitch"][0], head["element"].get("id")) for note, head in pairs],
            [("F", "diamond"), ("E", "cross")],
        )
        self.assertEqual([note["event_number"] for note, _ in pairs], [1, 1])
        reversed_xml = musicxml(xml_note() + xml_note("F", 3, "<chord/>"))
        self.assertEqual(
            [head["element"].get("id") for _, head in match(reversed_xml, drawing)],
            ["diamond", "cross"],
        )

    def test_hidden_chord_base_does_not_join_previous_attack(self):
        xml = musicxml(
            xml_note() + xml_note("C", attributes='print-object="no"')
            + xml_note("G", extra="<chord/>")
        )
        pairs = match(xml, svg([("first", 100, 400), ("chord", 200, 320)]))
        self.assertEqual([note["event_number"] for note, _ in pairs], [1, 2])

    def test_tied_member_does_not_change_other_chord_events(self):
        xml = musicxml(
            xml_note("C", extra='<tie type="start"/>') + xml_note("D")
            + xml_note("C", extra='<chord/><tie type="stop"/>')
            + xml_note("F", extra="<chord/>") + xml_note()
        )
        notes = score.assign_event_numbers(score.get_musicxml_notes(xml))
        self.assertEqual([n["event_number"] for n in notes], [1, 2, 1, 2, 3])

    def test_tie_chain_reuses_event_and_start_event_is_preserved(self):
        xml = musicxml(
            xml_note(extra='<notations><tied type="start"/></notations>')
            + xml_note(extra='<tie type="stop"/><tie type="start"/>')
            + xml_note(extra='<tie type="stop"/>') + xml_note("F")
        )
        notes = score.assign_event_numbers(score.get_musicxml_notes(xml), start_event=7)
        self.assertEqual([n["event_number"] for n in notes], [7, 7, 7, 8])

    def test_rests_and_invisible_heads_do_not_consume_svg_heads(self):
        xml = musicxml(
            '<note><rest/></note>' + xml_note(extra="<grace/>")
            + xml_note("F", extra="<notehead>none</notehead>") + xml_note("G")
        )
        pairs = match(xml, svg([("grace", 100, 400), ("normal", 200, 320)]))
        self.assertEqual([note["event_number"] for note, _ in pairs], [1, 3])

    def test_clef_changes_and_octave_clefs(self):
        xml = musicxml(
            xml_note("F", 3)
            + '<attributes><clef><sign>C</sign><line>3</line></clef></attributes>'
            + xml_note("C")
            + '<attributes><clef><sign>G</sign><line>2</line>'
              '<clef-octave-change>-1</clef-octave-change></clef></attributes>'
            + xml_note("G", 3),
            '<attributes><clef><sign>F</sign><line>4</line></clef></attributes>',
        )
        pairs = match(xml, svg([("bass", 100, 160), ("alto", 200, 240), ("octave", 300, 320)]))
        self.assertEqual([note["event_number"] for note, _ in pairs], [1, 2, 3])

    def test_equal_counts_with_wrong_pitch_leave_output_untouched(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            xml_path, svg_path, output = [directory / name for name in ("score.xml", "score.svg", "out.svg")]
            ET.ElementTree(musicxml(xml_note())).write(xml_path)
            ET.ElementTree(svg([("wrong", 100, 360)])).write(svg_path)
            output.write_text("keep existing output")
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(RuntimeError, "pitch mismatch"):
                    score.process(xml_path, svg_path, output)
            self.assertEqual(output.read_text(), "keep existing output")

    def test_count_mismatch_leaves_output_untouched(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            xml_path, svg_path, output = [directory / name for name in ("score.xml", "score.svg", "out.svg")]
            ET.ElementTree(musicxml(xml_note() + xml_note("F"))).write(xml_path)
            ET.ElementTree(svg([("only", 100, 400)])).write(svg_path)
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(RuntimeError, "count mismatch"):
                    score.process(xml_path, svg_path, output)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
