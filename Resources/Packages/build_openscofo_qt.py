#!/usr/bin/env python3
"""Build an offline Qt Installer Framework installer for OpenScofo on Linux.

Usage:
  python3 build_openscofo_qt.py --repo-root . --version 0.3.1
  python3 build_openscofo_qt.py --repo-root . --strict
  python3 build_openscofo_qt.py --repo-root . --binarycreator /path/to/binarycreator

Compiles nothing. Qt IFW binarycreator must be installed to produce the .run.
"""

import argparse
import datetime
import json
import shutil
import subprocess
import sys
from pathlib import Path
from xml.sax.saxutils import escape

COMPONENTS = {
    "csound": {
        "name": "Csound",
        "dest": ".local/lib/csound/6.0/plugins64",
        "binaries": ["build/Sources/Wrappers/CSound/OpenScofo.so"],
        "required": [],
        "examples": ["Sources/Wrappers/CSound/examples/1-score-follow.csd"],
    },
    "puredata": {
        "name": "Pure Data",
        "dest": "Documents/Pd/externals/openscofo~",
        "binaries": [
            "build/Sources/Wrappers/PureData/openscofo~.linux-amd64-32.so",
            "build/Sources/Wrappers/PureData/openscofo~.so",
        ],
        "required": ["Sources/Wrappers/PureData/openscofo~-help.pd"],
        "examples": [],
    },
    "supercollider": {
        "name": "SuperCollider",
        "dest": ".local/share/SuperCollider/Extensions/OpenScofo",
        "binaries": ["build/Sources/Wrappers/SuperCollider/OpenScofo.so"],
        "required": ["Sources/Wrappers/SuperCollider/OpenScofo.sc"],
        "examples": [
            "Sources/Wrappers/SuperCollider/examples/1-score-follow.scd",
            "Sources/Wrappers/SuperCollider/examples/2-descriptors-input.scd",
        ],
    },
    "vamp": {
        "name": "Vamp",
        "dest": ".vamp",
        "binaries": ["build/Sources/Wrappers/Vamp/OpenScofo.so"],
        "required": [],
        "examples": [],
    },
}

ASSETS = [
    "Tests/01-benchmark/real/bwv-1013.wav",
    "Tests/01-benchmark/real/canticos.wav",
    "Tests/01-benchmark/real/canticos.scofo",
    "Tests/01-benchmark/real/bwv-1013.scofo",
    "Tests/04-miniaturas/Audios/miniatura2.mp3",
    "Tests/04-miniaturas/Extras/miniatura2.scofo",
    "Tests/flute.onnx",
]


def copy_file(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def make_script(key, destination, filenames, version):
    """Generate Qt IFW operations: copy staged files to user-specific locations.

    Data layout: packages/<id>/data/<key>/<filename>.
    Qt IFW extracts data/ contents into @TargetDir@, so source is @TargetDir@/<key>/...
    """
    lines = [
        "function Component() {}",
        "",
        "Component.prototype.createOperations = function() {",
        '    var target = installer.value("TargetDir");',
        '    var archive = "installer://io.github.openscofo.' + key + '/" + '
        + json.dumps(version + key + ".7z") + ";",
        '    component.addOperation("Extract", archive, target);',
        '    var home = installer.value("HomeDir");',
        '    var base = target + "/" + ' + json.dumps(key) + ' + "/";',
        '    var destination = home + "/" + ' + json.dumps(destination) + ' + "/";',
        '    component.addOperation("Mkdir", destination);',
    ]
    for name in filenames:
        quoted = json.dumps(name)
        lines.extend([
            "    if (installer.fileExists(destination + " + quoted + ")) {",
            '        console.log("Skipping existing user file: " + destination + '
            + quoted + ");",
            "    } else {",
            '        component.addOperation("Copy", base + '
            + quoted + ", destination + " + quoted + ");",
            "    }",
        ])
    lines.extend(["};", ""])
    return "\n".join(lines)


def package_xml(name, version, default=True):
    return f'''<?xml version="1.0" encoding="UTF-8"?>
<Package>
    <DisplayName>{escape(name)}</DisplayName>
    <Description>OpenScofo {escape(name)}</Description>
    <Version>{escape(version)}</Version>
    <ReleaseDate>{datetime.date.today().isoformat()}</ReleaseDate>
    <Default>{str(default).lower()}</Default>
    <Script>installscript.qs</Script>
</Package>
'''


def create_package(packages, key, name, dest, version, sources, default=True):
    package = packages / ("io.github.openscofo." + key)
    data = package / "data" / key
    meta = package / "meta"
    data.mkdir(parents=True)
    meta.mkdir(parents=True)
    staged_names = []
    for source in sources:
        if source.name in staged_names:
            raise RuntimeError(f"{key}: duplicate filename {source.name}")
        if not source.is_file():
            print(f"[skip] {key}: source disappeared before staging: {source}")
            continue
        staged = data / source.name
        copy_file(source, staged)
        if not staged.is_file():
            raise RuntimeError(f"{key}: failed to stage: {source}")
        staged_names.append(staged.name)
    if not staged_names:
        shutil.rmtree(package)
        print(f"[skip] {key}: no files were staged")
        return False
    (meta / "package.xml").write_text(package_xml(name, version, default), encoding="utf-8")
    (meta / "installscript.qs").write_text(
        make_script(key, dest, sorted(staged_names), version), encoding="utf-8"
    )
    print(f"[ok] {key}: {len(staged_names)} existing files staged")
    return True


def find_binarycreator(requested):
    if requested:
        path = Path(requested).expanduser().resolve()
        if not path.is_file():
            raise RuntimeError(f"binarycreator not found: {path}")
        return str(path)
    found = shutil.which("binarycreator")
    if found:
        return found
    candidates = sorted((Path.home() / "Qt").glob("**/binarycreator"))
    return str(candidates[0]) if candidates else None


def build(args):
    repo = args.repo_root.expanduser().resolve()
    if not (repo / "Sources/Wrappers").is_dir():
        raise RuntimeError(f"Not an OpenScofo repository root: {repo}")

    root = repo / "Resources/Packages/qt"
    config = root / "config"
    packages = root / "packages"
    config.mkdir(parents=True, exist_ok=True)
    # Only remove generated package directories; never remove other project files.
    for key in [*COMPONENTS, "examples"]:
        path = packages / ("io.github.openscofo." + key)
        if path.exists():
            shutil.rmtree(path)
    packages.mkdir(parents=True, exist_ok=True)

    config_xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<Installer>
    <Name>OpenScofo</Name>
    <Version>{escape(args.version)}</Version>
    <Title>OpenScofo Installer</Title>
    <Publisher>OpenScofo</Publisher>
    <TargetDir>@HomeDir@/.local/share/openscofo-installer</TargetDir>
    <MaintenanceToolName>OpenScofoMaintenanceTool</MaintenanceToolName>
</Installer>
'''
    (config / "config.xml").write_text(config_xml, encoding="utf-8")

    selected_binaries = {}
    for key, spec in COMPONENTS.items():
        # 'binaries' are alternatives: one existing binary is sufficient.
        selected_binaries[key] = next(
            (repo / path for path in spec["binaries"] if (repo / path).is_file()),
            None,
        )

    missing = [key for key, binary in selected_binaries.items() if binary is None]
    if args.strict and missing:
        details = "; ".join(
            f"{key}: {', '.join(COMPONENTS[key]['binaries'])}" for key in missing
        )
        raise RuntimeError(f"missing required component binaries: {details}")

    available = []
    for key, spec in COMPONENTS.items():
        binary = selected_binaries[key]
        if binary is None:
            message = f"{key}: no compiled binary found"
            print("[skip]", message)
            continue
        sources = [binary]
        for rel in spec["required"]:
            path = repo / rel
            if path.is_file():
                sources.append(path)
            else:
                print(f"[skip] {key}: missing auxiliary file: {rel}")
        # Examples stay separate from the plugin directories.
        if create_package(packages, key, spec["name"], spec["dest"],
                          args.version, sources):
            available.append(key)

    if not available:
        raise RuntimeError("No compiled plugins found. Build the plugins first.")

    if not args.no_examples:
        assets = []
        for rel in ASSETS:
            path = repo / rel
            if path.is_file():
                assets.append(path)
            else:
                print(f"[optional] missing asset: {rel}")
        for key, spec in COMPONENTS.items():
            for rel in spec["examples"]:
                path = repo / rel
                if path.is_file():
                    assets.append(path)
                else:
                    print(f"[optional] missing example: {rel}")
        if assets:
            create_package(packages, "examples", "Examples and test assets",
                           ".local/share/openscofo/examples", args.version,
                           assets, default=False)

    output = args.output.expanduser()
    if not output.is_absolute():
        output = repo / output
    output.parent.mkdir(parents=True, exist_ok=True)
    binarycreator = find_binarycreator(args.binarycreator)
    if not binarycreator:
        print("[prepared] Packages are ready at", root)
        print("[next] Install Qt Installer Framework, then run:")
        print(f"  binarycreator --offline-only -c '{config / 'config.xml'}' -p '{packages}' '{output}'")
        return
    command = [binarycreator, "--offline-only", "-c", str(config / "config.xml"),
               "-p", str(packages), str(output)]
    print("[build]", " ".join(command))
    subprocess.run(command, check=True)
    output.chmod(output.stat().st_mode | 0o111)
    print("[done]", output)
    print("[note] Existing destination files are preserved and skipped.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--version", default="0.3.1")
    parser.add_argument("--output", type=Path, default=Path("OpenScofo-Linux.run"))
    parser.add_argument("--binarycreator")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--no-examples", action="store_true",
                        help="Exclude optional examples and test assets")
    args = parser.parse_args()
    try:
        build(args)
    except (OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
