#!/usr/bin/env python3

"""
OpenScofo Linux self-extracting installer builder.

Creates a single .sh file containing:
- a bash installer stub
- an appended tar.gz payload with component folders
"""

from __future__ import annotations

import argparse
import os
import shutil
import stat
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.chdir(Path(__file__).resolve().parent.parent.parent)


@dataclass
class ComponentSpec:
    key: str
    required: list[Path]
    optional: list[Path]


def copy_item(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(src, dest)
    else:
        shutil.copy2(src, dest)


def first_existing(paths: Iterable[Path]) -> Path | None:
    for candidate in paths:
        if candidate.exists():
            return candidate
    return None


def copy_if_exists(repo_root: Path, rel: Path, dest_dir: Path) -> bool:
    src = repo_root / rel
    if not src.exists():
        return False
    copy_item(src, dest_dir / src.name)
    return True


def component_specs() -> list[ComponentSpec]:
    shared_assets = [
        Path("Tests/assets/canticos.wav"),
        Path("Tests/assets/canticos.txt"),
        Path("Tests/miniaturas/Audios/miniatura1.mp3"),
        Path("Tests/miniaturas/Extras/miniatura1.scofo"),
    ]

    return [
        ComponentSpec(
            key="csound",
            required=[Path("Sources/Wrappers/CSound/examples/1-score-follow.csd"), *shared_assets],
            optional=[Path("build/Sources/Wrappers/CSound/OpenScofo.so")],
        ),
        ComponentSpec(
            key="puredata",
            required=[Path("Sources/Wrappers/PureData/openscofo~-help.pd"), *shared_assets],
            optional=[
                Path("build/Sources/Wrappers/PureData/openscofo~.linux-amd64-32.so"),
                Path("build/Sources/Wrappers/PureData/openscofo~.so"),
            ],
        ),
        ComponentSpec(
            key="supercollider",
            required=[
                Path("Sources/Wrappers/SuperCollider/OpenScofo.sc"),
                Path("Sources/Wrappers/SuperCollider/examples/1-score-follow.scd"),
                Path("Sources/Wrappers/SuperCollider/examples/2-descriptors-input.scd"),
                *shared_assets,
            ],
            optional=[Path("build/Sources/Wrappers/SuperCollider/OpenScofo.so")],
        ),
        ComponentSpec(
            key="vamp",
            required=[*shared_assets],
            optional=[Path("build/Sources/Wrappers/Vamp/OpenScofo.so")],
        ),
    ]


def stage_component_payload(repo_root: Path, payload_root: Path, spec: ComponentSpec) -> int:
    dest_dir = payload_root / spec.key
    dest_dir.mkdir(parents=True, exist_ok=True)

    copied = 0

    binary = first_existing(repo_root / rel for rel in spec.optional)
    if binary is not None:
        copy_item(binary, dest_dir / binary.name)
        copied += 1

    for rel in spec.required:
        if copy_if_exists(repo_root, rel, dest_dir):
            copied += 1

    # Include extra optional assets if present (for naming differences across branches/platforms).
    for rel in spec.optional:
        src = repo_root / rel
        if src.exists() and not (dest_dir / src.name).exists():
            copy_item(src, dest_dir / src.name)
            copied += 1

    return copied


def create_payload_archive(payload_root: Path, out_tgz: Path) -> None:
    out_tgz.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(out_tgz, "w:gz") as tf:
        tf.add(payload_root, arcname="payload")


def assemble_self_extracting_installer(stub_path: Path, payload_tgz: Path, output_path: Path) -> None:
    marker = b"__ARCHIVE_BELOW__\n"
    stub_bytes = stub_path.read_bytes()

    if marker not in stub_bytes:
        raise RuntimeError(f"Installer stub is missing marker: {stub_path}")

    marker_end = stub_bytes.index(marker) + len(marker)
    final_stub = stub_bytes[:marker_end]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as out:
        out.write(final_stub)
        with payload_tgz.open("rb") as payload:
            shutil.copyfileobj(payload, out)

    mode = output_path.stat().st_mode
    output_path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build OpenScofo Linux self-extracting installer")
    parser.add_argument("--version", default="0.1.0", help="Installer version label")
    parser.add_argument(
        "--repo-root",
        default=str(Path(__file__).resolve().parents[2]),
        help="Path to OpenScofo repository root",
    )
    parser.add_argument(
        "--work-dir",
        default="Resources/Packages/out-linux",
        help="Temporary build workspace",
    )
    parser.add_argument(
        "--stub",
        default="Resources/Packages/linux-installer.stub.sh",
        help="Path to installer shell stub",
    )
    parser.add_argument(
        "--output",
        default="OpenScofo-Linux.sh",
        help="Final self-contained installer path",
    )

    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    work_dir = (repo_root / args.work_dir).resolve()
    payload_root = work_dir / "payload"
    payload_tgz = work_dir / "payload.tar.gz"
    stub_path = (repo_root / args.stub).resolve()
    output_path = (repo_root / args.output).resolve()

    if not stub_path.exists():
        print(f"Error: installer stub not found: {stub_path}", file=sys.stderr)
        return 1

    if work_dir.exists():
        shutil.rmtree(work_dir)
    payload_root.mkdir(parents=True, exist_ok=True)

    print("[linux-builder] note: Max plugin packaging is skipped on Linux (no Linux Max external target).")

    print(f"[linux-builder] staging payload in: {payload_root}")
    total_files = 0
    for spec in component_specs():
        copied = stage_component_payload(repo_root, payload_root, spec)
        total_files += copied
        print(f"[linux-builder] {spec.key}: copied {copied} files")

    if total_files == 0:
        print("Error: payload is empty; no files staged.", file=sys.stderr)
        return 2

    create_payload_archive(payload_root, payload_tgz)
    print(f"[linux-builder] created archive: {payload_tgz}")

    assemble_self_extracting_installer(stub_path, payload_tgz, output_path)
    print(f"[linux-builder] built installer: {output_path}")
    print(f"[linux-builder] version: {args.version}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
