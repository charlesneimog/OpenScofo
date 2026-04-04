#!/usr/bin/env python3

"""
OpenScofo macOS component installer builder.

Creates a product installer with GUI choices so users can install any subset of
integrations (Max, Pure Data, CSound, SuperCollider, Vamp).
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

os.chdir(Path(__file__).resolve().parent.parent.parent)

@dataclass
class IntegrationSpec:
	key: str
	title: str
	identifier: str
	user_destination: str
	required: list[Path]
	optional: list[Path]


def run(cmd: list[str]) -> None:
	print("+", " ".join(cmd))
	subprocess.run(cmd, check=True)


def first_existing(paths: Iterable[Path]) -> Path | None:
	for path in paths:
		if path.exists():
			return path
	return None


def ensure_max_bundle(repo_root: Path) -> Path | None:
	for candidate in [
		repo_root / "max" / "o.scofo~.mxo",
		repo_root / "max" / "o.scofo~.mxo" / "Contents" / "MacOS" / "o.scofo~",
	]:
		if candidate.exists():
			return candidate if candidate.is_dir() else candidate.parent.parent

	build_dir = repo_root / "build"
	if build_dir.exists():
		run(["cmake", "--build", str(build_dir), "--target", "Max"])
		for candidate in [
			repo_root / "max" / "o.scofo~.mxo",
			repo_root / "max" / "o.scofo~.mxo" / "Contents" / "MacOS" / "o.scofo~",
		]:
			if candidate.exists():
				return candidate if candidate.is_dir() else candidate.parent.parent

	return None


def copy_item(src: Path, dest: Path) -> None:
	dest.parent.mkdir(parents=True, exist_ok=True)
	if src.is_dir():
		if dest.exists():
			shutil.rmtree(dest)
		shutil.copytree(src, dest)
	else:
		shutil.copy2(src, dest)


def copy_into_directory(src: Path, dest_dir: Path) -> None:
	copy_item(src, dest_dir / src.name)


def create_postinstall(script_path: Path, install_to: str, stage_root: str) -> None:
	script_path.write_text(
		"""#!/bin/bash
set -euo pipefail

logged_user="$(stat -f%Su /dev/console)"
if [[ -z "$logged_user" || "$logged_user" == "root" ]]; then
  echo "Could not detect a non-root console user."
  exit 1
fi

user_home="$(dscl . -read /Users/$logged_user NFSHomeDirectory | awk '{print $2}')"
if [[ -z "$user_home" ]]; then
  user_home="/Users/$logged_user"
fi

src_dir="%STAGE_ROOT%"

if [[ "%INSTALL_TO%" == /* ]]; then
	dst_dir="%INSTALL_TO%"
else
	dst_dir="$user_home/%INSTALL_TO%"
fi

mkdir -p "$dst_dir"
/usr/bin/ditto "$src_dir" "$dst_dir"

rm -rf "%STAGE_BASE%"
exit 0
"""
		.replace("%INSTALL_TO%", install_to)
		.replace("%STAGE_ROOT%", stage_root)
		.replace("%STAGE_BASE%", str(Path(stage_root).parent)),
		encoding="utf-8",
	)
	script_path.chmod(0o755)


def build_component_pkg(
	pkgbuild_bin: str,
	pkg_dir: Path,
	payload_src: Path,
	work_dir: Path,
	component: IntegrationSpec,
	version: str,
) -> Path:
	stage_base = f"/tmp/OpenScofoInstaller/{component.key}"
	stage_root = f"{stage_base}/payload"
	scripts_dir = work_dir / "scripts" / component.key
	scripts_dir.mkdir(parents=True, exist_ok=True)
	postinstall = scripts_dir / "postinstall"
	create_postinstall(postinstall, component.user_destination, stage_root)

	pkg_path = pkg_dir / f"OpenScofo-{component.key}.pkg"
	run(
		[
			pkgbuild_bin,
			"--root",
			str(payload_src),
			"--identifier",
			component.identifier,
			"--version",
			version,
			"--install-location",
			stage_root,
			"--scripts",
			str(scripts_dir),
			str(pkg_path),
		]
	)
	return pkg_path


def generate_distribution_xml(
	path: Path,
	package_name: str,
	version: str,
	built_components: list[tuple[IntegrationSpec, Path]],
) -> None:
	lines = [
		'<?xml version="1.0" encoding="utf-8"?>',
		'<installer-gui-script minSpecVersion="1">',
		f"    <title>{package_name} Installer</title>",
		'    <license file="License.txt"/>',
		'    <options require-scripts="false" customize="always"/>',
		"    <choices-outline>",
	]

	for spec, _ in built_components:
		lines.append(f'        <line choice="{spec.key}Choice"/>')
	lines.append("    </choices-outline>")

	for spec, pkg_path in built_components:
		lines.append(
			f'    <choice id="{spec.key}Choice" visible="true" start_selected="true" '
			f'title="{spec.title}">'
			f'<pkg-ref id="{spec.identifier}"/></choice>'
		)
		lines.append(
			f'    <pkg-ref id="{spec.identifier}" version="{version}" '
			f'onConclusion="none">{pkg_path.name}</pkg-ref>'
		)

	lines.append("</installer-gui-script>")
	path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def component_specs(repo_root: Path) -> list[IntegrationSpec]:
	return [
		IntegrationSpec(
			key="max",
			title="Max (Cycling 74)",
			identifier="org.openscofo.pkg.max",
			user_destination="Documents/Max 9/Packages/OpenScofo",
			required=[
				Path("Sources/Wrappers/Max/o.scofo~.maxhelp"),
				Path("Tests/assets/bwv-1013.wav"),
				Path("Tests/assets/bwv-1013.txt"),
				Path("Tests/assets/canticos.wav"),
				Path("Tests/assets/canticos.txt"),
			],
			optional=[
				Path("max/o.scofo~.mxo"),
				Path("max/o.scofo~.mxo/Contents/MacOS/o.scofo~"),
				Path("build/Sources/Wrappers/Max/o.scofo~.mxo"),
				Path("build/Sources/Wrappers/Max/o.scofo~.mxo/Contents/MacOS/o.scofo~"),
			],
		),
		IntegrationSpec(
			key="puredata",
			title="Pure Data",
			identifier="org.openscofo.pkg.puredata",
			user_destination="Documents/Pd/externals/o.scofo~",
			required=[
				Path("Sources/Wrappers/PureData/o.scofo~-help.pd"),
				Path("Tests/assets/bwv-1013.wav"),
				Path("Tests/assets/bwv-1013.txt"),
				Path("Tests/assets/canticos.wav"),
				Path("Tests/assets/canticos.txt"),
			],
			optional=[
				Path("build/Sources/Wrappers/PureData/o.scofo~.pd_darwin"),
				Path("build/Sources/Wrappers/PureData/o.scofo~.darwin-fat-32.so"),
				Path("build/Sources/Wrappers/PureData/o.scofo~.so"),
			],
		),
		IntegrationSpec(
			key="csound",
			title="CSound",
			identifier="org.openscofo.pkg.csound",
			user_destination="/Library/Frameworks/CsoundLib64.framework/Versions/6.0/Resources/Opcodes64",
			required=[
				Path("Sources/Wrappers/CSound/examples/1-score-follow.csd"),
				Path("Tests/assets/bwv-1013.wav"),
				Path("Tests/assets/bwv-1013.txt"),
				Path("Tests/assets/canticos.wav"),
				Path("Tests/assets/canticos.txt"),
			],
			optional=[
				Path("build/Sources/Wrappers/CSound/OpenScofo.dylib"),
				Path("build/Sources/Wrappers/CSound/OpenScofo.so"),
			],
		),
		IntegrationSpec(
			key="supercollider",
			title="SuperCollider",
			identifier="org.openscofo.pkg.supercollider",
			user_destination="Library/Application Support/SuperCollider/Extensions/OpenScofo",
			required=[
				Path("Sources/Wrappers/SuperCollider/OpenScofo.sc"),
				Path("Sources/Wrappers/SuperCollider/examples/1-score-follow.scd"),
				Path("Sources/Wrappers/SuperCollider/examples/2-descriptors-input.scd"),
				Path("Tests/assets/bwv-1013.wav"),
				Path("Tests/assets/bwv-1013.txt"),
				Path("Tests/assets/canticos.wav"),
				Path("Tests/assets/canticos.txt"),
			],
			optional=[Path("build/Sources/Wrappers/SuperCollider/OpenScofo.scx")],
		),
		IntegrationSpec(
			key="vamp",
			title="Vamp and Sonic Visualizer",
			identifier="org.openscofo.pkg.vamp",
			user_destination="/Library/Audio/Plug-Ins/Vamp",
			required=[
				Path("Tests/assets/bwv-1013.wav"),
				Path("Tests/assets/bwv-1013.txt"),
				Path("Tests/assets/canticos.wav"),
				Path("Tests/assets/canticos.txt"),
			],
			optional=[Path("build/Sources/Wrappers/Vamp/OpenScofo.dylib")],
		),
	]


def stage_component_payload(repo_root: Path, work_dir: Path, spec: IntegrationSpec) -> Path:
	payload_root = work_dir / "payload" / spec.key
	payload_root.mkdir(parents=True, exist_ok=True)

	if spec.key == "max":
		externals_dir = payload_root / "externals"
		help_dir = payload_root / "help"
		externals_dir.mkdir(parents=True, exist_ok=True)
		help_dir.mkdir(parents=True, exist_ok=True)

		bundle = ensure_max_bundle(repo_root)
		if bundle is not None:
			copy_into_directory(bundle, externals_dir)

		for rel in spec.required:
			src = repo_root / rel
			if src.exists():
				copy_into_directory(src, help_dir)
		return payload_root

	binary = first_existing(repo_root / p for p in spec.optional)
	if binary is not None:
		copy_item(binary, payload_root / binary.name)

	for rel in spec.required:
		src = repo_root / rel
		if src.exists():
			copy_item(src, payload_root / src.name)

	return payload_root


def copy_license(repo_root: Path, resources_dir: Path) -> None:
	src = repo_root / "LICENSE"
	resources_dir.mkdir(parents=True, exist_ok=True)
	if src.exists():
		shutil.copy2(src, resources_dir / "License.txt")
	else:
		(resources_dir / "License.txt").write_text(
			"OpenScofo license file was not found in source tree.\n",
			encoding="utf-8",
		)


def main() -> int:
	parser = argparse.ArgumentParser(description="Build OpenScofo macOS component installer")
	parser.add_argument("--version", default="0.1.0", help="Package version")
	parser.add_argument(
		"--repo-root",
		default=str(Path(__file__).resolve().parents[2]),
		help="Path to OpenScofo repository root",
	)
	parser.add_argument(
		"--work-dir",
		default="Resources/Packages/out",
		help="Temporary build workspace",
	)
	parser.add_argument(
		"--output",
		default="OpenScofo.pkg",
		help="Final productbuild package path",
	)
	parser.add_argument(
		"--package-name",
		default="OpenScofo",
		help="Installer title",
	)
	parser.add_argument(
		"--components",
		default="max,puredata,csound,supercollider,vamp",
		help="Comma-separated integrations to include",
	)

	args = parser.parse_args()

	repo_root = Path(args.repo_root).resolve()
	work_dir = (repo_root / args.work_dir).resolve()
	pkg_dir = work_dir / "pkgs"
	resources_dir = work_dir / "resources"
	scripts_dir = work_dir / "scripts"
	distribution_xml = work_dir / "distribution.xml"
	output_pkg = (repo_root / args.output).resolve()

	if work_dir.exists():
		shutil.rmtree(work_dir)
	pkg_dir.mkdir(parents=True, exist_ok=True)
	scripts_dir.mkdir(parents=True, exist_ok=True)

	pkgbuild_bin = shutil.which("pkgbuild")
	productbuild_bin = shutil.which("productbuild")
	if not pkgbuild_bin or not productbuild_bin:
		print("Error: pkgbuild and productbuild are required on macOS.", file=sys.stderr)
		return 1

	wanted = {c.strip() for c in args.components.split(",") if c.strip()}
	built_components: list[tuple[IntegrationSpec, Path]] = []

	for spec in component_specs(repo_root):
		if spec.key not in wanted:
			continue

		payload = stage_component_payload(repo_root, work_dir, spec)
		if not any(payload.iterdir()):
			print(f"Skipping {spec.key}: no files found")
			continue

		pkg_path = build_component_pkg(
			pkgbuild_bin=pkgbuild_bin,
			pkg_dir=pkg_dir,
			payload_src=payload,
			work_dir=work_dir,
			component=spec,
			version=args.version,
		)
		built_components.append((spec, pkg_path))

	if not built_components:
		print("Error: no component packages were built.", file=sys.stderr)
		return 2

	copy_license(repo_root, resources_dir)
	generate_distribution_xml(distribution_xml, args.package_name, args.version, built_components)

	output_pkg.parent.mkdir(parents=True, exist_ok=True)
	run(
		[
			productbuild_bin,
			"--resources",
			str(resources_dir),
			"--distribution",
			str(distribution_xml),
			"--package-path",
			str(pkg_dir),
			str(output_pkg),
		]
	)

	print(f"Built installer: {output_pkg}")
	return 0


if __name__ == "__main__":
	raise SystemExit(main())
