#!/usr/bin/env python3

from pathlib import Path
import xml.etree.ElementTree as ET
import html
import re


ROOT = Path(__file__).resolve().parents[1]
XML_DIR = ROOT / "build" / "doxygen" / "xml"
OUTPUT_DIR = ROOT / "Documentation" / "api"


# ----------------------------------------------------------------------
# XML helpers
# ----------------------------------------------------------------------

def text(node):
    if node is None:
        return ""

    value = "".join(node.itertext())
    value = html.unescape(value)
    value = re.sub(r"\s+", " ", value)

    return value.strip()


def paragraphs(node):
    if node is None:
        return ""

    result = []

    for para in node.findall(".//para"):
        value = text(para)

        if value:
            result.append(value)

    return "\n\n".join(result)


def get_return_description(member):
    for section in member.findall(".//simplesect"):
        if section.get("kind") == "return":
            return paragraphs(section)

    return ""


def get_notes(member):
    notes = []

    for section in member.findall(".//simplesect"):
        kind = section.get("kind")

        if kind in ("note", "warning", "attention"):
            value = paragraphs(section)

            if value:
                notes.append((kind, value))

    return notes


def get_parameters(member):
    parameters = []

    descriptions = {}

    for item in member.findall(".//parameteritem"):
        name = text(item.find(".//parametername"))
        description = paragraphs(item.find("parameterdescription"))

        if name:
            descriptions[name] = description

    for param in member.findall("param"):
        type_name = text(param.find("type"))
        name = text(param.find("declname"))

        if not name:
            continue

        parameters.append({
            "name": name,
            "type": type_name,
            "description": descriptions.get(name, ""),
        })

    return parameters


# ----------------------------------------------------------------------
# Classification
# ----------------------------------------------------------------------

def category(name):
    if name in ("OpenScofo", "~OpenScofo"):
        return "Construction"

    if name in (
        "LoadScore",
        "ScoreIsLoaded",
        "SetCurrentEvent",
        "SetCurrentSection",
        "GetCurrentScorePosition",
        "GetCurrentStateIndex",
        "GetCurrentBPM",
        "GetCurrentEventActions",
        "GetAudioStateChangeActions",
    ):
        return "Score Following"

    if name in (
        "ProcessBlock",
        "GetCurrentBufferIndex",
        "GetBlockDuration",
    ):
        return "Audio Processing"

    if (
        "Descriptor" in name
        or "Description" in name
        or name in (
            "GetPitchProb",
            "GetPitchTemplate",
            "ActivateAllDescriptors",
        )
    ):
        return "Audio Descriptors"

    if "ONNX" in name:
        return "ONNX"

    if name.startswith("Lua") or name == "GetLuaCode":
        return "Lua"

    if name in (
        "SetConfiguration",
        "UpdateConfiguration",
        "GetConfiguration",
        "GetSr",
        "GetFFTSize",
        "GetHopSize",
    ):
        return "Configuration"

    if "Error" in name or "Log" in name:
        return "Logging and Errors"

    return "Other"


CATEGORY_ORDER = [
    "Construction",
    "Configuration",
    "Score Following",
    "Audio Processing",
    "Audio Descriptors",
    "ONNX",
    "Lua",
    "Logging and Errors",
    "Other",
]


# ----------------------------------------------------------------------
# Doxygen
# ----------------------------------------------------------------------

def find_openscofo_class():
    index_file = XML_DIR / "index.xml"

    if not index_file.exists():
        raise RuntimeError(
            f"{index_file} does not exist.\n"
            "Run `doxygen Doxyfile.api` first."
        )

    root = ET.parse(index_file).getroot()

    for compound in root.findall("compound"):
        if compound.get("kind") != "class":
            continue

        name = text(compound.find("name"))

        if name.endswith("::OpenScofo") or name == "OpenScofo":
            return compound.get("refid")

    raise RuntimeError("Could not find the OpenScofo class in Doxygen XML.")


def read_functions(refid):
    xml_file = XML_DIR / f"{refid}.xml"

    root = ET.parse(xml_file).getroot()

    functions = []

    for member in root.findall(".//memberdef[@kind='function']"):
        protection = member.get("prot")

        if protection != "public":
            continue

        name = text(member.find("name"))

        brief = paragraphs(member.find("briefdescription"))
        detail = paragraphs(member.find("detaileddescription"))

        definition = text(member.find("definition"))
        args = text(member.find("argsstring"))

        return_type = text(member.find("type"))

        functions.append({
            "name": name,
            "brief": brief,
            "detail": detail,
            "definition": definition,
            "args": args,
            "return_type": return_type,
            "parameters": get_parameters(member),
            "returns": get_return_description(member),
            "notes": get_notes(member),
            "category": category(name),
        })

    return functions


# ----------------------------------------------------------------------
# Markdown
# ----------------------------------------------------------------------

def signature(function):
    definition = function["definition"]
    args = function["args"]

    if definition:
        return f"{definition}{args}"

    return f'{function["name"]}{args}'


def render_function(function):
    lines = []

    name = function["name"]

    lines.append(f"### `{name}`")
    lines.append("")
    lines.append("```cpp")
    lines.append(signature(function))
    lines.append("```")
    lines.append("")

    if function["brief"]:
        lines.append(function["brief"])
        lines.append("")

    parameters = function["parameters"]

    if parameters:
        lines.append("#### Parameters")
        lines.append("")
        lines.append("| Parameter | Type | Description |")
        lines.append("| --- | --- | --- |")

        for parameter in parameters:
            description = parameter["description"].replace("|", "\\|")

            lines.append(
                f'| `{parameter["name"]}` | '
                f'`{parameter["type"]}` | '
                f'{description} |'
            )

        lines.append("")

    if function["returns"]:
        lines.append("#### Returns")
        lines.append("")
        lines.append(function["returns"])
        lines.append("")

    for kind, value in function["notes"]:
        zensical_kind = {
            "note": "note",
            "warning": "warning",
            "attention": "warning",
        }.get(kind, "note")

        lines.append(f'!!! {zensical_kind}')
        lines.append("")

        for paragraph in value.splitlines():
            lines.append(f"    {paragraph}")

        lines.append("")

    return "\n".join(lines)


def generate(functions):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    destination = OUTPUT_DIR / "openscofo.md"

    lines = [
        "---",
        "title: OpenScofo API",
        "---",
        "",
        "# OpenScofo API",
        "",
        "Reference for the main `OpenScofo` C++ class.",
        "",
        "!!! info",
        "",
        "    This page is generated automatically from the OpenScofo source code.",
        "    Do not edit it manually.",
        "",
    ]

    grouped = {}

    for function in functions:
        grouped.setdefault(function["category"], []).append(function)

    for section in CATEGORY_ORDER:
        members = grouped.get(section)

        if not members:
            continue

        lines.append(f"## {section}")
        lines.append("")

        for function in members:
            lines.append(render_function(function))

    destination.write_text("\n".join(lines), encoding="utf-8")

    print(f"Generated {destination}")


def main():
    refid = find_openscofo_class()
    functions = read_functions(refid)

    if not functions:
        raise RuntimeError("No public OpenScofo functions found.")

    generate(functions)

    print(f"{len(functions)} public functions documented.")


if __name__ == "__main__":
    main()
