from lxml import etree
import re

svg_file = "/home/neimog/Documents/Git/OpenScofo/Resources/Score/miniatura2.svg"
xml_file = "/home/neimog/Documents/Git/OpenScofo/Resources/Score/miniatura2.musicxml"

# ╭──────────────────────────────────────╮
# │              SVG MATCH               │
# ╰──────────────────────────────────────╯
target_class = "Note"
svg_tree = etree.parse(svg_file)
svg_root = svg_tree.getroot()

note_svg_elems = []
for elem in svg_root.iter():
    class_attr = elem.get("class")
    if class_attr:
        classes = class_attr.split()
        if target_class in classes:
            note_svg_elems.append(elem)


def note_position(note):
    d = note.get("d", "")

    m = re.match(r"[Mm]\s*([-\d.]+)\s*,\s*([-\d.]+)", d)
    if m is None:
        return (float("inf"), float("inf"))

    x = float(m.group(1))
    y = float(m.group(2))
    return (x, y)


# Extract coordinates
notes_with_pos = []

for note in note_svg_elems:
    x, y = note_position(note)
    notes_with_pos.append((note, x, y))

# Sort by vertical position first
notes_with_pos.sort(key=lambda t: t[2])

# Group notes into systems
SYSTEM_THRESHOLD = 700  # adjust if necessary

systems = []

for note, x, y in notes_with_pos:
    if not systems:
        systems.append([(note, x, y)])
        continue

    first_y = systems[-1][0][2]

    if abs(y - first_y) < SYSTEM_THRESHOLD:
        systems[-1].append((note, x, y))
    else:
        systems.append([(note, x, y)])

# Sort each system left-to-right
note_svg_elems = []

for system in systems:
    system.sort(key=lambda t: t[1])  # x coordinate
    for note, x, y in system:
        note_svg_elems.append(note)


# ╭──────────────────────────────────────╮
# │            MUSICXML MATCH            │
# ╰──────────────────────────────────────╯
xml_tree = etree.parse(xml_file)
xml_root = xml_tree.getroot()

note_xml_elems = []
for elem in xml_root.iter("note"):
    if elem.find("rest") is not None:
        continue
    note_xml_elems.append(elem)


event_number = 1
for svg, xml in zip(note_svg_elems, note_xml_elems):
    is_tied = False
    if xml.find("tie") is not None:
        tie = xml.find("tie")
        tie_type = tie.get("type")  # "start" or "stop"
        if tie_type == "start":
            is_tied = True

    svg.set("openscofo-event-number", str(event_number))

    if not is_tied:
        event_number += 1


svg_tree.write("/home/neimog/Documents/Git/OpenScofo/Resources/Score/score.svg")
