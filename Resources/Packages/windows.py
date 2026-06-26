#!/usr/bin/env python3

"""
OpenScofo Windows component installer builder.

Creates an MSI installer with selectable features so users can install any
subset of integrations (Max, Pure Data, CSound, SuperCollider, Vamp).
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from xml.sax.saxutils import escape

os.chdir(Path(__file__).resolve().parent.parent.parent)


@dataclass
class IntegrationSpec:
    key: str
    title: str
    identifier: str
    destination_root: str
    destination_subdir: str
    destination_subdir_x64: str | None
    required: list[Path]
    optional: list[Path]


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(cwd) if cwd else None)


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def sanitize_id(text: str) -> str:
    out = []
    for ch in text:
        if ch.isalnum() or ch == "_":
            out.append(ch)
        else:
            out.append("_")
    sanitized = "".join(out).strip("_")
    if not sanitized:
        sanitized = "X"
    if sanitized[0].isdigit():
        sanitized = f"X_{sanitized}"
    return sanitized


def short_hash(text: str) -> str:
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:10]


def escape_attr(value: str) -> str:
    return escape(value, {'"': "&quot;", "'": "&apos;"})


def normalize_msi_version(version: str) -> str:
    parts = re.findall(r"\d+", version)
    if not parts:
        return "0.0.0"

    numeric = [str(min(int(part), 65534)) for part in parts[:3]]
    while len(numeric) < 3:
        numeric.append("0")

    return ".".join(numeric)


def to_wix_path(path: Path) -> str:
    return str(path).replace("/", "\\")


def ensure_max_binary(repo_root: Path, build_dir: Path) -> Path | None:
    candidates = [
        repo_root / "max" / "openscofo~.mxe64",
        repo_root / "max" / "openscofo~.mxe",
        build_dir / "Sources" / "Wrappers" / "Max" / "openscofo~.mxe64",
        build_dir / "Sources" / "Wrappers" / "Max" / "openscofo~.mxe",
    ]
    existing = first_existing(candidates)
    if existing is not None:
        return existing

    if build_dir.exists():
        for target in ["Max", "max_o.scofo"]:
            try:
                run(["cmake", "--build", str(build_dir), "--target", target])
                existing = first_existing(candidates)
                if existing is not None:
                    return existing
            except subprocess.CalledProcessError:
                continue

    return None


def copy_item(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(src, dest)
    else:
        shutil.copy2(src, dest)


def stage_max_package_payload(repo_root: Path, build_dir: Path, payload_root: Path) -> int:
    max_root = repo_root / "Sources" / "Wrappers" / "Max"
    copied = 0

    package_root = max_root / "package"
    for filename in ["License.md", "Readme.md", "icon.png", "package-info.json"]:
        src = package_root / filename
        if src.exists():
            copy_item(src, payload_root / filename)
            copied += 1

    docs = [
        max_root / "docs" / "openscofo~.maxref.xml",
    ]
    for src in docs:
        if src.exists():
            copy_item(src, payload_root / "docs" / src.name)
            copied += 1

    help_files = [
        max_root / "openscofo~.maxhelp",
        max_root / "score-max.svg",
        max_root / "score_renderer.js",
        repo_root / "Tests" / "miniaturas" / "Extras" / "ai-flute-model.onnx",
        repo_root / "Tests" / "assets" / "canticos.txt",
        repo_root / "Tests" / "assets" / "canticos.wav",
        repo_root / "Tests" / "miniaturas" / "Audios" / "miniatura1.mp3",
    ]
    for src in help_files:
        if src.exists():
            copy_item(src, payload_root / "help" / src.name)
            copied += 1

    max_score = max_root / "miniatura1-max.scofo"
    if max_score.exists():
        copy_item(max_score, payload_root / "help" / "miniatura1.scofo")
        copied += 1

    binary = ensure_max_binary(repo_root, build_dir)
    if binary is not None:
        copy_item(binary, payload_root / "externals" / binary.name)
        copied += 1

    return copied


def build_destination(spec: IntegrationSpec, arch: str) -> tuple[str, list[str]]:
    sub = spec.destination_subdir_x64 if arch == "x64" and spec.destination_subdir_x64 else spec.destination_subdir
    parts = [p for p in Path(sub).parts if p not in ("", ".")]
    root = spec.destination_root
    if arch == "x64" and root == "ProgramFilesFolder":
        root = "ProgramFiles64Folder"
    return root, parts


def component_specs() -> list[IntegrationSpec]:

    # check files ending with dll inside build/Sources/Wrappers/PureData/
    pd_candidates = list((Path("build") / "Sources" / "Wrappers" / "PureData").glob("*.dll"))
    if not pd_candidates:
        print("Warning: no Pure Data wrapper binary found in build directory. The Pure Data component will be skipped.")


    return [
        IntegrationSpec(
            key="max",
            title="Max (Cycling 74)",
            identifier="org.openscofo.pkg.max",
            destination_root="PersonalFolder",
            destination_subdir="Max 9/Packages/OpenScofo",
            destination_subdir_x64=None,
            required=[
                Path("Sources/Wrappers/Max/package/License.md"),
                Path("Sources/Wrappers/Max/package/Readme.md"),
                Path("Sources/Wrappers/Max/package/icon.png"),
                Path("Sources/Wrappers/Max/package/package-info.json"),
                Path("Sources/Wrappers/Max/docs/openscofo~.maxref.xml"),
                Path("Sources/Wrappers/Max/openscofo~.maxhelp"),
                Path("Sources/Wrappers/Max/score-max.svg"),
                Path("Sources/Wrappers/Max/score_renderer.js"),
                Path("Tests/miniaturas/Extras/ai-flute-model.onnx"),
                Path("Tests/assets/canticos.wav"),
                Path("Tests/assets/canticos.txt"),
                Path("Tests/miniaturas/Audios/miniatura1.mp3"),
                Path("Sources/Wrappers/Max/miniatura1-max.scofo"),
            ],
            optional=[
                Path("max/openscofo~.mxe64"),
                Path("max/openscofo~.mxe"),
                Path("build/Sources/Wrappers/Max/openscofo~.mxe64"),
                Path("build/Sources/Wrappers/Max/openscofo~.mxe"),
            ],
        ),
        IntegrationSpec(
            key="puredata",
            title="Pure Data",
            identifier="org.openscofo.pkg.puredata",
            destination_root="PersonalFolder",
            destination_subdir="Pd/externals/openscofo~",
            destination_subdir_x64=None,
            required=[
                Path("Sources/Wrappers/PureData/openscofo~-help.pd"),
                Path("Tests/assets/canticos.wav"),
                Path("Tests/assets/canticos.txt"),
                Path("Tests/miniaturas/Audios/miniatura1.mp3"),
                Path("Tests/miniaturas/Extras/miniatura1.scofo"),
            ],
            optional=pd_candidates,
        ),
        IntegrationSpec(
            key="csound",
            title="CSound",
            identifier="org.openscofo.pkg.csound",
            destination_root="LocalAppDataFolder",
            destination_subdir="csound/6.0/plugins64",
            destination_subdir_x64=None,
            required=[
                Path("Sources/Wrappers/CSound/examples/1-score-follow.csd"),
                Path("Tests/assets/canticos.wav"),
                Path("Tests/assets/canticos.txt"),
                Path("Tests/miniaturas/Audios/miniatura1.mp3"),
                Path("Tests/miniaturas/Extras/miniatura1.scofo"),
            ],
            optional=[
                Path("build/Sources/Wrappers/CSound/OpenScofo.dll"),
            ],
        ),
        IntegrationSpec(
            key="supercollider",
            title="SuperCollider",
            identifier="org.openscofo.pkg.supercollider",
            destination_root="AppDataFolder",
            destination_subdir="SuperCollider/Extensions/OpenScofo",
            destination_subdir_x64=None,
            required=[
                Path("Sources/Wrappers/SuperCollider/OpenScofo.sc"),
                Path("Sources/Wrappers/SuperCollider/examples/1-score-follow.scd"),
                Path("Sources/Wrappers/SuperCollider/examples/2-descriptors-input.scd"),
                Path("Tests/assets/canticos.wav"),
                Path("Tests/assets/canticos.txt"),
                Path("Tests/miniaturas/Audios/miniatura1.mp3"),
                Path("Tests/miniaturas/Extras/miniatura1.scofo"),
            ],
            optional=[Path("build/Sources/Wrappers/SuperCollider/OpenScofo.scx")],
        ),
        IntegrationSpec(
            key="vamp",
            title="Vamp and Sonic Visualizer",
            identifier="org.openscofo.pkg.vamp",
            destination_root="ProgramFilesFolder",
            destination_subdir="Vamp Plugins",
            destination_subdir_x64=None,
            required=[
                Path("Tests/assets/canticos.wav"),
                Path("Tests/assets/canticos.txt"),
                Path("Tests/miniaturas/Audios/miniatura1.mp3"),
                Path("Tests/miniaturas/Extras/miniatura1.scofo"),
            ],
            optional=[Path("build/Sources/Wrappers/Vamp/OpenScofo.dll")],
        ),
    ]


def stage_component_payload(repo_root: Path, build_dir: Path, work_dir: Path, spec: IntegrationSpec) -> Path:
    payload_root = work_dir / "payload" / spec.key
    payload_root.mkdir(parents=True, exist_ok=True)

    if spec.key == "max":
        stage_max_package_payload(repo_root, build_dir, payload_root)
        return payload_root

    binary = first_existing(repo_root / p for p in spec.optional)
    if binary is not None:
        copy_item(binary, payload_root / binary.name)

    for rel in spec.required:
        src = repo_root / rel
        if src.exists():
            copy_item(src, payload_root / src.name)

    return payload_root


def iter_files_under(path: Path) -> Iterable[Path]:
    for candidate in sorted(path.rglob("*")):
        if candidate.is_file():
            yield candidate


def generate_wix_source(
    path: Path,
    package_name: str,
    manufacturer: str,
    version: str,
    upgrade_code: str,
    arch: str,
    built_components: list[tuple[IntegrationSpec, Path, tuple[str, list[str]]]],
) -> None:
    user_profile_roots = {"PersonalFolder", "AppDataFolder", "LocalAppDataFolder"}
    root_paths: dict[str, set[tuple[str, ...]]] = {}
    directory_ids: dict[tuple[str, tuple[str, ...]], str] = {}
    file_groups: dict[tuple[str, str, tuple[str, ...]], list[Path]] = {}
    destination_leaf_ids: dict[str, str] = {}
    integration_by_key: dict[str, IntegrationSpec] = {spec.key: spec for spec, _, _ in built_components}
    integration_dir_paths: dict[str, set[tuple[str, tuple[str, ...]]]] = {}

    for spec, payload_root, (root, destination_parts) in built_components:
        integration_dir_paths.setdefault(spec.key, set())
        paths = root_paths.setdefault(root, set())
        for depth in range(1, len(destination_parts) + 1):
            path_parts = tuple(destination_parts[:depth])
            paths.add(path_parts)
            integration_dir_paths[spec.key].add((root, path_parts))

        dest_tuple = tuple(destination_parts)
        destination_leaf_ids[spec.key] = (
            directory_ids.setdefault((root, dest_tuple), sanitize_id(f"DIR_{short_hash(root + ':' + '/'.join(dest_tuple))}"))
            if dest_tuple
            else root
        )

        for file_path in iter_files_under(payload_root):
            rel_parent = file_path.relative_to(payload_root).parent
            rel_parts = tuple(rel_parent.parts) if rel_parent != Path(".") else ()
            full_parts = tuple(destination_parts) + rel_parts

            for depth in range(1, len(full_parts) + 1):
                path_parts = full_parts[:depth]
                paths.add(path_parts)
                integration_dir_paths[spec.key].add((root, path_parts))

            group_key = (spec.key, root, full_parts)
            file_groups.setdefault(group_key, []).append(file_path)

    for root, paths in root_paths.items():
        for parts in paths:
            directory_ids.setdefault(
                (root, parts), sanitize_id(f"DIR_{short_hash(root + ':' + '/'.join(parts))}")
            )

    lines: list[str] = []

    lines.append('<?xml version="1.0" encoding="utf-8"?>')
    lines.append('<Wix xmlns="http://schemas.microsoft.com/wix/2006/wi">')
    lines.append(
        f'  <Product Id="*" Name="{escape_attr(package_name)}" Language="1033" '
        f'Version="{escape_attr(version)}" Manufacturer="{escape_attr(manufacturer)}" '
        f'UpgradeCode="{escape_attr(upgrade_code)}">'
    )
    lines.append('    <Package InstallerVersion="500" Compressed="yes" InstallScope="perMachine"/>')
    lines.append('    <MediaTemplate EmbedCab="yes" CompressionLevel="high"/>')
    lines.append(
        '    <MajorUpgrade DowngradeErrorMessage="A newer version of [ProductName] is already installed."/>'
    )
    lines.append('    <UIRef Id="WixUI_FeatureTree"/>')
    lines.append('    <Directory Id="TARGETDIR" Name="SourceDir">')

    standard_roots = [
            "ProgramFilesFolder",
            "ProgramFiles64Folder",
            "PersonalFolder",
            "AppDataFolder",
            "LocalAppDataFolder",
            "CommonAppDataFolder",
    ]

    for root in standard_roots:
        paths = sorted(root_paths.get(root, set()), key=lambda p: (len(p), p))
        if not paths:
            continue

        lines.append(f'      <Directory Id="{root}">')

        children: dict[tuple[str, ...], list[tuple[str, ...]]] = {}
        for parts in paths:
            parent = parts[:-1]
            children.setdefault(parent, []).append(parts)

        def emit_nodes(parent: tuple[str, ...], indent: str) -> None:
            for parts in sorted(children.get(parent, [])):
                dir_id = directory_ids[(root, parts)]
                name = parts[-1]
                lines.append(f'{indent}<Directory Id="{escape_attr(dir_id)}" Name="{escape_attr(name)}">')
                emit_nodes(parts, indent + "  ")
                lines.append(f"{indent}</Directory>")

        emit_nodes((), "        ")
        lines.append("      </Directory>")

    lines.append("    </Directory>")

    all_component_ids: list[str] = []
    integration_components: dict[str, list[str]] = {}

    for spec, _, _ in built_components:
        integration_components.setdefault(spec.key, [])

    for (spec_key, root, full_parts), files in sorted(file_groups.items(), key=lambda x: x[0]):
        spec = integration_by_key[spec_key]
        wix_dir_id = directory_ids[(root, full_parts)] if full_parts else destination_leaf_ids[spec_key]
        component_key = f"{spec.key}:{root}:{'/'.join(full_parts)}"
        component_id = sanitize_id(f"CMP_{spec.key}_{short_hash(component_key)}")
        component_guid = str(
            uuid.uuid5(
                uuid.NAMESPACE_URL,
                f"openscofo/{spec.identifier}/{root}/{'/'.join(full_parts)}",
            )
        )

        lines.append(f'    <DirectoryRef Id="{escape_attr(wix_dir_id)}">')
        lines.append(
            f'      <Component Id="{escape_attr(component_id)}" Guid="{escape_attr(component_guid)}" '
            f'Win64="{"yes" if arch == "x64" else "no"}">'
        )

        per_user_component = root in user_profile_roots
        first_file = True
        for file_path in sorted(files):
            rel_for_id = file_path.name
            file_id = sanitize_id(
                f"FIL_{spec.key}_{short_hash(str(file_path.parent) + ':' + rel_for_id)}"
            )
            source_path = to_wix_path(file_path.resolve())
            keypath_attr = ' KeyPath="yes"' if first_file and not per_user_component else ""
            lines.append(
                f'        <File Id="{escape_attr(file_id)}" Source="{escape_attr(source_path)}"{keypath_attr}/>'
            )
            first_file = False

        if per_user_component:
            reg_id = sanitize_id(f"REG_{spec.key}_{short_hash(component_key)}")
            reg_name = sanitize_id(f"K_{short_hash(component_key)}")
            lines.append(
                '        <RegistryValue Root="HKCU" '
                'Key="Software\\OpenScofo\\Installer\\Components" '
                f'Name="{escape_attr(reg_name)}" Type="integer" Value="1" '
                f'Id="{escape_attr(reg_id)}" KeyPath="yes"/>'
            )

        lines.append("      </Component>")
        lines.append("    </DirectoryRef>")

        all_component_ids.append(component_id)
        integration_components[spec_key].append(component_id)

    for spec, _, _ in built_components:
        for root, parts in sorted(integration_dir_paths.get(spec.key, set()), key=lambda x: (x[0], len(x[1]), x[1])):
            if root not in user_profile_roots:
                continue
            if not parts:
                continue

            dir_id = directory_ids[(root, parts)]
            cleanup_key = f"{spec.key}:{root}:{'/'.join(parts)}:cleanup"
            cleanup_id = sanitize_id(f"CLN_{spec.key}_{short_hash(cleanup_key)}")
            cleanup_guid = str(
                uuid.uuid5(
                    uuid.NAMESPACE_URL,
                    f"openscofo/{spec.identifier}/{root}/{'/'.join(parts)}/cleanup",
                )
            )
            remove_file_id = sanitize_id(f"RMF_{short_hash(cleanup_key)}")
            remove_folder_id = sanitize_id(f"RMD_{short_hash(cleanup_key)}")
            reg_id = sanitize_id(f"RGC_{short_hash(cleanup_key)}")
            reg_name = sanitize_id(f"C_{short_hash(cleanup_key)}")

            lines.append(f'    <DirectoryRef Id="{escape_attr(dir_id)}">')
            lines.append(
                f'      <Component Id="{escape_attr(cleanup_id)}" Guid="{escape_attr(cleanup_guid)}" '
                f'Win64="{"yes" if arch == "x64" else "no"}">'
            )
            lines.append(
                f'        <RemoveFile Id="{escape_attr(remove_file_id)}" Name="*" On="uninstall"/>'
            )
            lines.append(
                f'        <RemoveFolder Id="{escape_attr(remove_folder_id)}" On="uninstall"/>'
            )
            lines.append(
                '        <RegistryValue Root="HKCU" '
                'Key="Software\\OpenScofo\\Installer\\Cleanup" '
                f'Name="{escape_attr(reg_name)}" Type="integer" Value="1" '
                f'Id="{escape_attr(reg_id)}" KeyPath="yes"/>'
            )
            lines.append("      </Component>")
            lines.append("    </DirectoryRef>")

            all_component_ids.append(cleanup_id)
            integration_components[spec.key].append(cleanup_id)

    lines.append('    <Feature Id="MainFeature" Title="OpenScofo" Level="1" Display="expand">')
    for spec, _, _ in built_components:
        feature_id = sanitize_id(f"Feature_{spec.key}")
        lines.append(
            f'      <Feature Id="{feature_id}" Title="{escape_attr(spec.title)}" Level="1" AllowAdvertise="no" Absent="disallow">'
        )
        for component_id in integration_components.get(spec.key, []):
            lines.append(f'        <ComponentRef Id="{escape_attr(component_id)}"/>')
        lines.append("      </Feature>")
    lines.append("    </Feature>")

    lines.append("  </Product>")
    lines.append("</Wix>")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_wix_msi(
    candle_bin: str,
    light_bin: str,
    work_dir: Path,
    wxs_path: Path,
    output_msi: Path,
    arch: str,
) -> None:
    wixobj = wxs_path.with_suffix(".wixobj")

    run([candle_bin, "-arch", arch, str(wxs_path)], cwd=work_dir)
    run(
        [
            light_bin,
            str(wixobj.name),
            "-ext",
            "WixUIExtension",
            "-out",
            str(output_msi.name),
        ],
        cwd=work_dir,
    )


def resolve_executable(bin_name: str) -> str | None:
    candidate = Path(bin_name)
    if candidate.is_file():
        return str(candidate.resolve())
    return shutil.which(bin_name)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build OpenScofo Windows MSI component installer")
    parser.add_argument("--version", default="0.1.0", help="Package version")
    parser.add_argument(
        "--repo-root",
        default=str(Path(__file__).resolve().parents[2]),
        help="Path to OpenScofo repository root",
    )
    parser.add_argument("--work-dir", default="Resources/Packages/out-win", help="Temporary build workspace")
    parser.add_argument("--output", default="OpenScofo-Windows.msi", help="Final MSI output path")
    parser.add_argument("--package-name", default="OpenScofo", help="Installer product name")
    parser.add_argument("--manufacturer", default="OpenScofo", help="Manufacturer string")
    parser.add_argument(
        "--upgrade-code",
        default="54a05500-6d7a-479d-b884-ba844ccaf4ce",
        help="Stable WiX UpgradeCode GUID",
    )
    parser.add_argument(
        "--arch",
        default="x64",
        choices=["x86", "x64"],
        help="MSI architecture",
    )
    parser.add_argument(
        "--components",
        default="max,puredata,csound,supercollider,vamp",
        help="Comma-separated integrations to include",
    )
    parser.add_argument(
        "--build-dir",
        default="build",
        help="CMake build directory to search for wrapper binaries",
    )
    parser.add_argument(
        "--candle",
        default="candle",
        help="Path to WiX candle executable",
    )
    parser.add_argument(
        "--light",
        default="light",
        help="Path to WiX light executable",
    )

    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    build_dir = (repo_root / args.build_dir).resolve()
    work_dir = (repo_root / args.work_dir).resolve()
    output_msi = (repo_root / args.output).resolve()

    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    wanted = {c.strip() for c in args.components.split(",") if c.strip()}
    built_components: list[tuple[IntegrationSpec, Path, tuple[str, list[str]]]] = []

    for spec in component_specs():
        if spec.key not in wanted:
            continue

        payload = stage_component_payload(repo_root, build_dir, work_dir, spec)
        if not any(iter_files_under(payload)):
            print(f"Skipping {spec.key}: no files found")
            continue

        destination = build_destination(spec, args.arch)
        built_components.append((spec, payload, destination))

    if not built_components:
        print("Error: no component payloads were staged.", file=sys.stderr)
        return 2

    wxs_path = work_dir / "OpenScofo.wxs"
    generate_wix_source(
        path=wxs_path,
        package_name=args.package_name,
        manufacturer=args.manufacturer,
        version=normalize_msi_version(args.version),
        upgrade_code=args.upgrade_code,
        arch=args.arch,
        built_components=built_components,
    )

    candle_bin = resolve_executable(args.candle)
    light_bin = resolve_executable(args.light)
    if not candle_bin or not light_bin:
        print("Error: WiX tools candle and light are required.", file=sys.stderr)
        return 1

    output_msi.parent.mkdir(parents=True, exist_ok=True)
    local_out = work_dir / output_msi.name

    build_wix_msi(
        candle_bin=candle_bin,
        light_bin=light_bin,
        work_dir=work_dir,
        wxs_path=wxs_path,
        output_msi=local_out,
        arch=args.arch,
    )

    shutil.copy2(local_out, output_msi)
    print(f"Built installer: {output_msi}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
