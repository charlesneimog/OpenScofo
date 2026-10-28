#!/usr/bin/env python3

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from statistics import median

SYSTEM_THRESHOLD = 500.0
SVG_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
DIATONIC_STEPS = {step: index for index, step in enumerate("CDEFGAB")}


def local_name(tag):
    """Remove namespace from an XML tag."""
    return tag.split("}")[-1]


def children(element, name):
    return [child for child in element if local_name(child.tag) == name]


def first_child(element, name):
    for child in element:
        if local_name(child.tag) == name:
            return child
    return None


def get_text(element, name, default=None):
    child = first_child(element, name)
    if child is None:
        return default
    return child.text


def get_pitch(note):
    """
    Return a hashable pitch representation.

    Example:
        C#4 -> ('C', 1, 4)
        Eb5 -> ('E', -1, 5)
    """
    pitch = first_child(note, "pitch")

    if pitch is None:
        return None

    step = get_text(pitch, "step")
    alter = float(get_text(pitch, "alter", "0"))
    octave = int(get_text(pitch, "octave", "0"))

    return step, alter, octave


def get_voice(note):
    return get_text(note, "voice", "1")


def get_tie_types(note):
    """
    Get all tie types associated with a note.

    Handles both:

        <tie type="start"/>

    and:

        <notations>
            <tied type="start"/>
        </notations>
    """
    types = set()

    for element in note.iter():
        tag = local_name(element.tag)

        if tag in ("tie", "tied"):
            tie_type = element.get("type")

            if tie_type in ("start", "stop"):
                types.add(tie_type)

    return types


def is_rest(note):
    return first_child(note, "rest") is not None


def is_chord(note):
    return first_child(note, "chord") is not None


def get_musicxml_notes(root):
    """
    Return MusicXML notes, retaining hidden notes for event/tie numbering.

    Rests are ignored because there is no SVG notehead corresponding
    to them.
    """
    result = []

    part_index = -1

    for element in root.iter():
        if local_name(element.tag) == "part":
            part_index += 1

            clefs = {"1": ("G", 2, 0)}
            attacks_by_voice = {}
            attack_number = 0
            for note in element.iter():
                if local_name(note.tag) == "clef":
                    clefs[note.get("number", "1")] = (
                        get_text(note, "sign", "G"),
                        int(get_text(note, "line", "2")),
                        int(get_text(note, "clef-octave-change", "0")),
                    )
                    continue
                if local_name(note.tag) != "note":
                    continue

                voice = get_voice(note)
                chord = is_chord(note)
                if not chord or voice not in attacks_by_voice:
                    attack_number += 1
                    attacks_by_voice[voice] = attack_number
                if is_rest(note):
                    continue

                # Hidden notes still take part in event/tie numbering, but do
                # not consume an SVG notehead. Grace notes do have noteheads.
                staff = get_text(note, "staff", "1")
                result.append(
                    {
                        "element": note,
                        "part": part_index,
                        "voice": voice,
                        "pitch": get_pitch(note),
                        "chord": chord,
                        "attack": attacks_by_voice[voice],
                        "ties": get_tie_types(note),
                        "staff": staff,
                        "clef": clefs.get(staff, ("G", 2, 0)),
                        "visible": (
                            note.get("print-object", "yes") != "no"
                            and get_text(note, "notehead") != "none"
                        ),
                    }
                )

    return result


def extract_svg_position(path):
    """
    Extract the first moveto coordinate from an SVG path.

    Examples:

        M 123,456
        m 123.4,456.7
        m123.4 456.7
    """
    d = path.get("d", "")

    match = re.search(
        rf"[Mm]\s*({SVG_NUMBER})[\s,]*({SVG_NUMBER})",
        d,
    )

    if not match:
        return None

    return float(match.group(1)), float(match.group(2))


def svg_path_points(path):
    """Read endpoints and control points of the outline of a notehead.

    MuseScore's first moveto is on the outline, and its offset differs for
    oval, diamond and cross noteheads. The outline's bounds give a consistent
    center for those glyphs. Relative commands and compact numbers are allowed.
    """
    tokens = re.findall(rf"[A-Za-z]|{SVG_NUMBER}", path.get("d", ""))
    sizes = {"M": 2, "L": 2, "H": 1, "V": 1, "C": 6, "S": 4, "Q": 4, "T": 2}
    points = []
    cursor = (0.0, 0.0)
    start = cursor
    command = None
    index = 0
    while index < len(tokens):
        if tokens[index].isalpha():
            command = tokens[index]
            index += 1
            if command.upper() == "Z":
                cursor = start
                command = None
                continue
        if command is None or command.upper() not in sizes:
            raise RuntimeError(f"Unsupported notehead path command: {command}")
        size = sizes[command.upper()]
        values = tokens[index:index + size]
        if len(values) != size or any(value.isalpha() for value in values):
            raise RuntimeError("Incomplete SVG notehead path command")
        values = [float(value) for value in values]
        index += size
        relative = command.islower()
        if command.upper() == "H":
            segment = [(values[0] + (cursor[0] if relative else 0), cursor[1])]
        elif command.upper() == "V":
            segment = [(cursor[0], values[0] + (cursor[1] if relative else 0))]
        else:
            segment = [
                (values[i] + (cursor[0] if relative else 0),
                 values[i + 1] + (cursor[1] if relative else 0))
                for i in range(0, size, 2)
            ]
        points.extend(segment)
        cursor = segment[-1]
        if command.upper() == "M":
            start = cursor
            command = "l" if relative else "L"
    return points


def get_svg_staffs(root):
    """Locate five-line staves independently of ledger notes and glyph size."""
    lines = set()
    for element in root.iter():
        if "StaffLines" not in element.get("class", "").split():
            continue
        if local_name(element.tag) == "polyline":
            values = [float(value) for value in re.findall(SVG_NUMBER, element.get("points", ""))]
        elif local_name(element.tag) == "line":
            values = [float(element.get(name, "0")) for name in ("x1", "y1", "x2", "y2")]
        else:
            continue
        for index in range(0, len(values) - 3, 2):
            x1, y1, x2, y2 = values[index:index + 4]
            if abs(y2 - y1) < 0.1 and x1 != x2:
                lines.add(round((y1 + y2) / 2, 3))

    lines = sorted(lines)
    staffs = []
    index = 0
    while index + 5 <= len(lines):
        group = lines[index:index + 5]
        gaps = [b - a for a, b in zip(group, group[1:])]
        spacing = median(gaps)
        if spacing > 0 and all(abs(gap - spacing) < spacing * 0.02 for gap in gaps):
            staffs.append({"top": group[0], "bottom": group[-1], "spacing": spacing})
            index += 5
        else:
            index += 1
    return staffs


def get_svg_notes(root):
    """
    Get SVG paths whose class contains 'Note'.
    """
    notes = []
    staffs = get_svg_staffs(root)

    for element in root.iter():
        if local_name(element.tag) != "path":
            continue

        classes = element.get("class", "").split()

        if "Note" not in classes:
            continue

        points = svg_path_points(element)
        if not points:
            raise RuntimeError("Could not determine SVG notehead position")
        xs, ys = zip(*points)
        x = (min(xs) + max(xs)) / 2
        y = (min(ys) + max(ys)) / 2
        staff_index = None
        if staffs:
            staff_index = min(
                range(len(staffs)),
                key=lambda index: abs(y - (staffs[index]["top"] + staffs[index]["bottom"]) / 2),
            )

        notes.append(
            {
                "element": element,
                "x": x,
                "y": y,
                "width": max(xs) - min(xs),
                "staff_index": staff_index,
                "staff": staffs[staff_index] if staff_index is not None else None,
            }
        )

    return notes


def sort_svg_notes(notes, system_threshold=SYSTEM_THRESHOLD):
    """
    Sort SVG noteheads in musical reading order.

    Use detected staves, then read each staff left-to-right. If staff lines
    are unavailable, split rows at vertical gaps between adjacent notes.
    """

    if notes and all(note.get("staff_index") is not None for note in notes):
        return sorted(notes, key=lambda note: (note["staff_index"], note["x"], note["y"]))

    notes = sorted(notes, key=lambda n: n["y"])

    systems = []

    for note in notes:
        placed = False

        for system in systems:
            reference_y = system[-1]["y"]

            if abs(note["y"] - reference_y) < system_threshold:
                system.append(note)
                placed = True
                break

        if not placed:
            systems.append([note])

    # Systems top -> bottom
    systems.sort(key=lambda system: min(note["y"] for note in system))

    result = []

    for system in systems:
        # Notes left -> right.
        #
        # y is only used to make ordering deterministic when several
        # noteheads have the same horizontal position (e.g. chords).
        system.sort(key=lambda n: (n["x"], n["y"]))
        result.extend(system)

    return result


def staff_pitch_position(note):
    """Diatonic steps above the bottom staff line (accidentals do not move it)."""
    pitch = note["pitch"]
    sign, line, octave_change = note["clef"]
    reference = {"G": ("G", 4), "F": ("F", 3), "C": ("C", 4)}.get(sign)
    if pitch is None or reference is None:
        return None
    step, _, octave = pitch
    ref_step, ref_octave = reference
    return (7 * octave + DIATONIC_STEPS[step]
            - 7 * (ref_octave + octave_change) - DIATONIC_STEPS[ref_step]
            + 2 * (line - 1))


def match_noteheads(xml_notes, svg_notes):
    """Pair attacks in reading order, then pair chord members by written pitch.

    Note counts alone cannot detect a misplaced ledger note. Validate the
    written pitch against the staff before accepting any attributes.
    """
    pairs = []
    index = 0
    while index < len(xml_notes):
        end = index + 1
        while (end < len(xml_notes)
               and xml_notes[end]["attack"] == xml_notes[index]["attack"]
               and xml_notes[end]["part"] == xml_notes[index]["part"]
               and xml_notes[end]["voice"] == xml_notes[index]["voice"]):
            end += 1
        group = xml_notes[index:end]
        heads = svg_notes[index:end]
        if len(group) > 1:
            if (len({head["staff_index"] for head in heads}) != 1
                    or max(head["x"] for head in heads) - min(head["x"] for head in heads)
                    > 1.5 * max(head["width"] for head in heads)):
                raise RuntimeError(f"Could not locate chord at MusicXML note {index + 1}")
            group = sorted(group, key=lambda note: (
                7 * note["pitch"][2] + DIATONIC_STEPS[note["pitch"][0]]
                if note["pitch"] else 0
            ))
            heads = sorted(heads, key=lambda head: -head["y"])
        for note, head in zip(group, heads):
            position = staff_pitch_position(note)
            staff = head["staff"]
            if staff is not None and position is not None:
                expected_y = staff["bottom"] - position * staff["spacing"] / 2
                if abs(expected_y - head["y"]) > staff["spacing"] * 0.3:
                    raise RuntimeError(
                        f"Notehead pitch mismatch at MusicXML note {index + 1}: "
                        f"{note['pitch']} expects y={expected_y:.2f}, "
                        f"SVG has y={head['y']:.2f}. "
                        "Aborting instead of assigning incorrect event numbers."
                    )
            pairs.append((note, head))
        index = end
    return pairs


def assign_event_numbers(notes, start_event=1):
    """
    Assign OpenScofo event numbers.

    Important rules:

    1. Normal note:
           creates a new event.

    2. Chord note (<chord/>):
           belongs to the same event as the preceding note
           in the same part/voice.

    3. Tie start:
           creates an event and remembers its event number.

    4. Tie stop:
           reuses the event number from the tie start.

    5. Tie stop + start:
           reuses the previous event AND keeps the tie open.

       Example:

           C4 start          -> 1
           C4 stop + start   -> 1
           C4 stop           -> 1
           D4                -> 2
    """

    next_event = start_event

    # (part, voice, pitch) -> event number
    active_ties = {}

    # Used by <chord/> notes.
    last_event_by_voice = {}

    for item in notes:
        part = item["part"]
        voice = item["voice"]
        pitch = item["pitch"]
        chord = item["chord"]
        ties = item["ties"]

        voice_key = (part, voice)
        tie_key = (part, voice, pitch)

        has_start = "start" in ties
        has_stop = "stop" in ties

        event_number = None

        # ----------------------------------------------------------
        # TIED CONTINUATION
        # ----------------------------------------------------------

        if has_stop and tie_key in active_ties:
            event_number = active_ties[tie_key]

        # ----------------------------------------------------------
        # CHORD MEMBER
        # ----------------------------------------------------------

        elif chord and voice_key in last_event_by_voice:
            event_number = last_event_by_voice[voice_key]

        # ----------------------------------------------------------
        # NEW EVENT
        # ----------------------------------------------------------

        else:
            event_number = next_event
            next_event += 1

        item["event_number"] = event_number

        if not chord:
            last_event_by_voice[voice_key] = event_number

        # ----------------------------------------------------------
        # UPDATE TIE STATE
        # ----------------------------------------------------------

        if has_start:
            # This also correctly handles stop+start:
            # the same event is simply kept active.
            active_ties[tie_key] = event_number

        elif has_stop:
            # Final note of the tie chain.
            active_ties.pop(tie_key, None)

    return notes


def register_svg_namespace(root):
    """
    Prevent ElementTree from writing ns0:path, ns0:svg, etc.
    """
    if root.tag.startswith("{"):
        namespace = root.tag.split("}")[0][1:]
        ET.register_namespace("", namespace)


def process(musicxml_path, svg_path, output_path, start_event=1):
    print(f"MusicXML: {musicxml_path}")
    print(f"SVG:      {svg_path}")

    # MUSICXML
    musicxml_tree = ET.parse(musicxml_path)
    musicxml_root = musicxml_tree.getroot()

    xml_notes = get_musicxml_notes(musicxml_root)
    assign_event_numbers(xml_notes, start_event=start_event)
    xml_notes = [note for note in xml_notes if note["visible"]]

    # SVG
    svg_tree = ET.parse(svg_path)
    svg_root = svg_tree.getroot()
    register_svg_namespace(svg_root)
    svg_notes = get_svg_notes(svg_root)
    svg_notes = sort_svg_notes(svg_notes)

    # SANITY CHECK
    print()
    print(f"MusicXML noteheads: {len(xml_notes)}")
    print(f"SVG noteheads:      {len(svg_notes)}")

    if len(xml_notes) != len(svg_notes):
        raise RuntimeError(
            "\nMusicXML/SVG note count mismatch.\n"
            f"MusicXML: {len(xml_notes)}\n"
            f"SVG:      {len(svg_notes)}\n\n"
            "Aborting instead of assigning incorrect event numbers."
        )

    # APPLY EVENT NUMBERS
    pairs = match_noteheads(xml_notes, svg_notes)
    for xml_note, svg_note in pairs:
        event_number = xml_note["event_number"]
        svg_note["element"].set(
            "openscofo-event-number",
            str(event_number),
        )

    # OUTPUT
    svg_tree.write(
        output_path,
        encoding="utf-8",
        xml_declaration=True,
    )

    unique_events = {note["event_number"] for note in xml_notes}

    print()
    print(f"OpenScofo events:   {len(unique_events)}")
    print(f"Output:             {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Add openscofo-event-number attributes to SVG noteheads "
            "using a MusicXML score."
        )
    )

    parser.add_argument(
        "musicxml",
        type=Path,
        help="Input MusicXML file",
    )

    parser.add_argument(
        "svg",
        type=Path,
        help="Input SVG file",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output SVG file",
    )

    parser.add_argument(
        "--start-event",
        type=int,
        default=1,
        help="First OpenScofo event number (default: 1)",
    )

    args = parser.parse_args()

    output = args.output

    if output is None:
        output = args.svg.with_name(args.svg.stem + "-openscofo.svg")

    process(
        args.musicxml,
        args.svg,
        output,
        start_event=args.start_event,
    )


if __name__ == "__main__":
    main()
