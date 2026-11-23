#!/usr/bin/env python3
"""OpenScofo CI entry point.

The workflow deliberately contains only orchestration.  Platform-specific build,
validation, installer, and staging logic lives here and is selected from the host
operating system.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import platform
import re
import shlex
import shutil
import ssl
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Callable, Iterable, Sequence

ROOT = Path(__file__).resolve().parents[2]
RUNNER_TEMP = Path(os.environ.get("RUNNER_TEMP", ROOT / ".ci-tmp"))
LIBNAME = os.environ.get("LIBNAME", "openscofo~")
SYSTEM = {"Darwin": "macos", "Linux": "linux", "Windows": "windows"}.get(
    platform.system(), platform.system().lower()
)


def run(
    *command: str | Path,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    capture: bool = False,
) -> str:
    """Run a command, echoing it in a copy/paste-friendly form."""
    args = [str(part) for part in command]
    display_command = (
        subprocess.list2cmdline(args) if SYSTEM == "windows" else shlex.join(args)
    )
    print(f"+ {display_command}", flush=True)
    result = subprocess.run(
        args,
        cwd=cwd or ROOT,
        env=env,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
    )
    return result.stdout.strip() if capture else ""


def require_platform(*supported: str) -> None:
    if SYSTEM not in supported:
        names = ", ".join(supported)
        raise RuntimeError(f"This action supports {names}; detected {SYSTEM}.")


def require_file(path: Path, *, executable: bool = False) -> Path:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Required non-empty file was not found: {path}")
    if executable and not os.access(path, os.X_OK):
        raise PermissionError(f"Required file is not executable: {path}")
    return path


def first_file(
    patterns: Iterable[str], *, exclude_suffix: str | None = None
) -> Path | None:
    matches: list[Path] = []
    for pattern in patterns:
        matches.extend(ROOT.glob(pattern))
    for path in sorted(set(matches)):
        if path.is_file() and path.stat().st_size > 0:
            if exclude_suffix and path.name.endswith(exclude_suffix):
                continue
            return path
    return None


def copy_file(source: str | Path, destination: str | Path) -> None:
    source_path = ROOT / source
    destination_path = ROOT / destination
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, destination_path)


def copy_tree(source: str | Path, destination: str | Path) -> None:
    source_path = ROOT / source
    destination_path = ROOT / destination
    shutil.copytree(source_path, destination_path, dirs_exist_ok=True)


def append_github_file(variable: str, lines: Sequence[str]) -> None:
    destination = os.environ.get(variable)
    if not destination:
        print(f"{variable} is not set; values: {', '.join(lines)}")
        return
    with Path(destination).open("a", encoding="utf-8") as stream:
        for line in lines:
            stream.write(f"{line}\n")


def download(
    url: str,
    destination: Path,
    attempts: int = 3,
    *,
    sha256: str | None = None,
    verify_tls: bool = True,
) -> None:
    if not verify_tls and not sha256:
        raise ValueError(
            "TLS verification may only be disabled for a hash-pinned download"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    context = None if verify_tls else ssl._create_unverified_context()
    for attempt in range(1, attempts + 1):
        try:
            print(f"Downloading {url} -> {destination}", flush=True)
            with (
                urllib.request.urlopen(url, timeout=120, context=context) as response,
                partial.open("wb") as output,
            ):
                shutil.copyfileobj(response, output)
            if sha256:
                digest = hashlib.sha256(partial.read_bytes()).hexdigest()
                if digest.casefold() != sha256.casefold():
                    raise RuntimeError(
                        f"SHA-256 mismatch for {destination.name}: expected {sha256}, got {digest}"
                    )
            partial.replace(destination)
            return
        except Exception:
            partial.unlink(missing_ok=True)
            if attempt == attempts:
                raise
            time.sleep(2 * attempt)


def pip_install(*packages: str, upgrade: bool = False) -> None:
    command = [sys.executable, "-m", "pip", "install"]
    if upgrade:
        command.append("--upgrade")
    run(*command, *packages)


def npm_install_global(*packages: str) -> None:
    if SYSTEM == "windows":
        npm_command = shutil.which("npm.cmd")
        node = shutil.which("node.exe")
        if not npm_command:
            raise FileNotFoundError("npm.cmd was not found on PATH")
        if not node:
            raise FileNotFoundError("node.exe was not found on PATH")
        npm_cli_candidates = (
            Path(npm_command).parent / "node_modules" / "npm" / "bin" / "npm-cli.js",
            Path(node).parent / "node_modules" / "npm" / "bin" / "npm-cli.js",
        )
        npm_cli = next((path for path in npm_cli_candidates if path.is_file()), None)
        if npm_cli is None:
            raise FileNotFoundError(
                "npm-cli.js was not found beside the Node.js installation"
            )
        require_file(npm_cli)
        run(node, npm_cli, "install", "-g", *packages)
        return
    run("npm", "install", "-g", *packages)


def reset_sccache() -> None:
    run("sccache", "--zero-stats")


def show_sccache_stats() -> None:
    run("sccache", "--show-stats")


def install_dependencies() -> None:
    """Install build dependencies for the detected host OS."""
    if SYSTEM == "linux":
        run("sudo", "apt-get", "update")
        run(
            "sudo",
            "apt-get",
            "install",
            "-y",
            "ffmpeg",
            "puredata",
            "lld",
            "libboost-all-dev",
            "libgtest-dev",
            "libssl-dev",
            "ragel",
            "swig",
            "yasm",
            "protobuf-compiler",
            "protobuf-c-compiler",
            "libxcb-image0",
            "libxcb-shape0",
            "libxcb-keysyms1",
            "libxcb-render-util0",
            "libxcb-icccm4",
            "libxcb-xinerama0",
            "libxcb-xkb1",
            "libxcb-xfixes0",
            "libxcb-cursor0",
            "libxcb-randr0",
            "libxkbcommon-x11-0",
            "libgl1",
            "libopengl0",
        )
        llvm_installer = RUNNER_TEMP / "llvm.sh"
        download("https://apt.llvm.org/llvm.sh", llvm_installer)
        llvm_installer.chmod(0o755)
        run("sudo", llvm_installer, "20")
        run("cargo", "install", "tree-sitter-cli", "--locked")
        pip_install(
            "build", "scikit-build-core", "cibuildwheel", "twine", "nanobind", "conan"
        )
        return

    if SYSTEM == "macos":
        run(
            "brew",
            "install",
            "ffmpeg",
            "pd",
            "nanobind",
            "tree-sitter",
            "boost",
            "protobuf",
            "protobuf-c",
            "openssl@3",
            "ragel",
            "swig",
            "yasm",
        )
        npm_install_global("tree-sitter-cli")
        pip_install("build", "scikit-build-core", "cibuildwheel", "twine", "conan")
        return

    if SYSTEM == "windows":
        npm_install_global("tree-sitter-cli")
        winget_common = (
            "--silent",
            "--disable-interactivity",
            "--accept-source-agreements",
            "--accept-package-agreements",
        )

        # run("winget", "install", "--exact", "--id", "LLVM.LLVM", *winget_common)
        # run(
        #     "winget",
        #     "install",
        #     "--exact",
        #     "--id",
        #     "LLVM.LLVM",
        #     "--version",
        #     "21.1.8",
        #     "--force",
        #     *winget_common,
        # )
        run(
            "choco",
            "install",
            "llvm",
            "--version=21.1.8",
            "--allow-downgrade",
            "--force",
            "-y",
        )

        pd_installer_url = (
            "https://msp.ucsd.edu/Software/pd-0.57-0.windows-installer.exe"
        )
        pd_installer_sha256 = (
            "6C1DC1DE251A9E130B8229165EF5932B54F89672376A78ED5FD650C6945439A8"
        )
        pd_installer = RUNNER_TEMP / "pd-0.57-0.windows-installer.exe"
        try:
            download(
                pd_installer_url, pd_installer, attempts=1, sha256=pd_installer_sha256
            )
        except Exception as error:
            print(
                f"Verified TLS download failed ({error}); retrying the hash-pinned installer "
                "without certificate validation.",
                flush=True,
            )
            download(
                pd_installer_url,
                pd_installer,
                sha256=pd_installer_sha256,
                verify_tls=False,
            )
        run(pd_installer, "/S", "/norestart")

        pd_root = Path(os.environ.get("ProgramFiles", "C:/Program Files")) / "Pd"
        require_file(pd_root / "src" / "m_pd.h")
        require_file(pd_root / "bin" / "pd.lib")
        append_github_file(
            "GITHUB_ENV",
            [f"PD_SOURCES_PATH={pd_root / 'src'}", f"PDBINDIR={pd_root / 'bin'}"],
        )
        pip_install("pip", upgrade=True)
        pip_install(
            "build", "scikit-build-core", "cibuildwheel", "twine", "nanobind", "conan"
        )
        return

    raise RuntimeError(f"Unsupported operating system: {SYSTEM}")


def install_qt_framework() -> None:
    require_platform("linux")
    installer = RUNNER_TEMP / "QtInstallerFramework-linux-x64-4.8.1.run"
    destination = RUNNER_TEMP / "QtIFW-4.8.1"
    download(
        "https://download.qt.io/official_releases/qt-installer-framework/4.8.1/"
        "QtInstallerFramework-linux-x64-4.8.1.run",
        installer,
    )
    installer.chmod(0o755)
    run(
        installer,
        "--root",
        destination,
        "--accept-licenses",
        "--default-answer",
        "--confirm-command",
        "install",
    )
    require_file(destination / "bin" / "binarycreator", executable=True)


def build(precision: int) -> None:
    command = [
        sys.executable,
        "-m",
        "build",
        "--wheel",
        "-Cbuild-dir=build",
        f"-Ccmake.define.PD_FLOATSIZE={precision}",
        "-Ccmake.define.CMAKE_BUILD_TYPE=Release",
        "-Ccmake.define.CMAKE_C_COMPILER_LAUNCHER=sccache",
        "-Ccmake.define.CMAKE_CXX_COMPILER_LAUNCHER=sccache",
        "-Ccmake.define.OPENSCOFO_BUILD_TESTS=ON",
        "-Ccmake.define.OPENSCOFO_ENABLE_NATIVE_OPTIMIZATIONS=OFF",
        f"-Ccmake.define.PD_LIB_DIR={ROOT / LIBNAME}",
        "-Ccmake.define.OPENSCOFO_BUILD_ALL=ON",
        "-Ccmake.define.OPENSCOFO_TRAINING=ON",
        "-Ccmake.args=-G Ninja",
    ]
    env = os.environ.copy()

    if SYSTEM == "linux":
        command[5:5] = [
            "-Ccmake.define.CMAKE_C_COMPILER=/usr/bin/clang-20",
            "-Ccmake.define.CMAKE_CXX_COMPILER=/usr/bin/clang++-20",
        ]
        command.insert(-3, f"-Ccmake.define.CMAKE_INSTALL_PREFIX={ROOT / 'install'}")
    elif SYSTEM == "macos":
        openssl = run("brew", "--prefix", "openssl@3", capture=True)
        machine = platform.machine().lower()
        macos_arch = "arm64" if machine in {"arm64", "aarch64"} else "x86_64"
        deployment_target = "11.0" if macos_arch == "arm64" else "10.15"
        command[5:5] = [f"-Ccmake.define.OPENSSL_ROOT_DIR={openssl}"]
        command.insert(-3, f"-Ccmake.define.CMAKE_OSX_ARCHITECTURES={macos_arch}")
        command.insert(
            -3, f"-Ccmake.define.CMAKE_OSX_DEPLOYMENT_TARGET={deployment_target}"
        )
        command.insert(-3, f"-Ccmake.define.CMAKE_INSTALL_PREFIX={ROOT / 'install'}")
        env["ARCHFLAGS"] = f"-arch {macos_arch}"
        env["MACOSX_DEPLOYMENT_TARGET"] = deployment_target
    elif SYSTEM == "windows":
        tree_sitter = shutil.which("tree-sitter.cmd") or shutil.which("tree-sitter")
        if not tree_sitter:
            raise FileNotFoundError("tree-sitter executable was not found on PATH")
        git_config_count = int(env.get("GIT_CONFIG_COUNT", "0"))
        env["GIT_CONFIG_COUNT"] = str(git_config_count + 1)
        env[f"GIT_CONFIG_KEY_{git_config_count}"] = "core.longpaths"
        env[f"GIT_CONFIG_VALUE_{git_config_count}"] = "true"
        command[5:5] = [
            f"-Ccmake.define.TREE_SITTER_EXECUTABLE={tree_sitter}",
            "-Ccmake.define.CMAKE_C_COMPILER=clang-cl",
            "-Ccmake.define.CMAKE_CXX_COMPILER=clang-cl",
            f"-Ccmake.define.PD_SOURCES_PATH={os.environ.get('PD_SOURCES_PATH', '')}",
            f"-Ccmake.define.PDBINDIR={os.environ.get('PDBINDIR', '')}",
        ]
    else:
        raise RuntimeError(f"Unsupported operating system: {SYSTEM}")

    run(*command, env=env)


def run_native_tests() -> None:
    run(
        "ctest",
        "--test-dir",
        "build",
        "-C",
        "Release",
        "--output-on-failure",
        "--no-tests=error",
        "-L",
        "openscofo",
    )


def trainer_binary() -> Path | None:
    suffix = {"linux": ".so", "windows": ".dll"}.get(SYSTEM)
    patterns = [
        "build/Sources/Trainer/openscofo.train*",
        "build/OpenScofo/openscofo.train*",
    ]
    for path in sorted(path for pattern in patterns for path in ROOT.glob(pattern)):
        if not path.is_file() or path.suffix == ".pd" or path.stat().st_size == 0:
            continue
        if suffix and not path.name.endswith(suffix):
            continue
        return path
    return None


def verify_native_trainer() -> None:
    trainer = trainer_binary()
    if not trainer:
        candidates = sorted(
            path
            for directory in (ROOT / "build/Sources/Trainer", ROOT / "build/OpenScofo")
            if directory.exists()
            for path in directory.iterdir()
            if path.is_file()
        )
        raise FileNotFoundError(
            "The openscofo.train Pure Data external was not generated. Candidates: "
            + ", ".join(str(path.relative_to(ROOT)) for path in candidates)
        )
    print(f"Verified native trainer: {trainer.relative_to(ROOT)}")


def repair_linux_wheel() -> None:
    require_platform("linux")
    pip_install("auditwheel", "patchelf>=0.14.5", upgrade=True)
    run("patchelf", "--version")
    wheelhouse = ROOT / "wheelhouse"
    shutil.rmtree(wheelhouse, ignore_errors=True)
    wheelhouse.mkdir()
    wheels = sorted((ROOT / "dist").glob("*-linux_x86_64.whl"))
    if not wheels:
        raise FileNotFoundError("No dist/*-linux_x86_64.whl wheel was found")
    run(
        "auditwheel",
        "repair",
        "--plat",
        "manylinux_2_35_x86_64",
        *wheels,
        "-w",
        wheelhouse,
    )
    for wheel in wheels:
        wheel.unlink()
    repaired = sorted(wheelhouse.glob("*.whl"))
    if not repaired:
        raise FileNotFoundError("auditwheel did not produce a repaired wheel")
    for wheel in repaired:
        shutil.copy2(wheel, ROOT / "dist" / wheel.name)
        print(f"Repaired wheel: {wheel.name} ({wheel.stat().st_size} bytes)")


def build_emscripten() -> None:
    require_platform("linux")
    run(
        "emcmake",
        "cmake",
        ".",
        "-B",
        "build-wasm",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_C_COMPILER_LAUNCHER=sccache",
        "-DCMAKE_CXX_COMPILER_LAUNCHER=sccache",
        "-DOPENSCOFO_BUILD_TESTS=ON",
        "-DOPENSCOFO_ENABLE_NATIVE_OPTIMIZATIONS=OFF",
        "-DOPENSCOFO_BUILD_ALL=ON",
    )
    run("cmake", "--build", "build-wasm")


def extract_wheel_version() -> None:
    wheels = sorted((ROOT / "dist").glob("*.whl"))
    if not wheels:
        raise FileNotFoundError("No wheel found in dist/")
    parts = wheels[0].stem.split("-")
    if len(parts) < 2:
        raise ValueError(f"Invalid wheel filename: {wheels[0].name}")
    name, version = parts[:2]
    print(
        f"Extracted WHEEL_NAME={name} and WHEEL_VERSION={version} from {wheels[0].name}"
    )
    append_github_file("GITHUB_ENV", [f"WHEEL_NAME={name}", f"WHEEL_VERSION={version}"])
    append_github_file("GITHUB_OUTPUT", [f"name={name}", f"version={version}"])


def wheel_version() -> str:
    value = os.environ.get("WHEEL_VERSION")
    if not value:
        raise RuntimeError(
            "WHEEL_VERSION is required; run --extract-wheel-version first"
        )
    return value


def files_test() -> None:
    """Validate platform-specific files used by the component installer."""
    if SYSTEM == "linux":
        for relative in (
            "build/Sources/Wrappers/CSound/OpenScofo.so",
            "build/Sources/Wrappers/SuperCollider/OpenScofo.so",
            "build/Sources/Wrappers/Vamp/OpenScofo.so",
        ):
            require_file(ROOT / relative)
        pd_binary = first_file(
            (
                "build/Sources/Wrappers/PureData/openscofo~.linux-amd64-32.so",
                "build/Sources/Wrappers/PureData/openscofo~.so",
            )
        )
        if not pd_binary:
            raise FileNotFoundError("No Pure Data component binary was generated")
        verify_native_trainer()
        return

    payload = (
        ROOT
        / "Resources/Packages"
        / ("out" if SYSTEM == "macos" else "out-win")
        / "payload/puredata"
    )
    for filename in (
        "flute.onnx",
        "miniatura1.scofo",
        "miniatura1.wav",
        "openscofo.train-help.pd",
    ):
        require_file(payload / filename)
    expected_suffix = ".dll" if SYSTEM == "windows" else None
    trainers = [
        path
        for path in payload.glob("openscofo.train*")
        if path.is_file()
        and path.suffix != ".pd"
        and (expected_suffix is None or path.name.endswith(expected_suffix))
    ]
    if not trainers:
        raise FileNotFoundError(
            f"The {SYSTEM} installer payload is missing openscofo.train"
        )
    score = (payload / "miniatura1.scofo").read_text(encoding="utf-8")
    if not re.search(r"^ONNXMODEL\s+flute\.onnx(?:\s|$)", score, re.MULTILINE):
        raise ValueError(
            "miniatura1.scofo does not reference the adjacent flute.onnx model"
        )
    print(f"Validated installer payload: {payload.relative_to(ROOT)}")


def generate_installer() -> None:
    version = wheel_version()
    if SYSTEM == "macos":
        run(sys.executable, "Resources/Packages/macos.py", "--version", version)
        return
    if SYSTEM == "linux":
        files_test()
        binarycreator = RUNNER_TEMP / "QtIFW-4.8.1/bin/binarycreator"
        require_file(binarycreator, executable=True)
        run(
            sys.executable,
            "Resources/Packages/build_openscofo_qt.py",
            "--repo-root",
            ".",
            "--version",
            version,
            "--binarycreator",
            binarycreator,
            "--output",
            "OpenScofo-Linux.run",
            "--strict",
        )
        require_file(ROOT / "OpenScofo-Linux.run", executable=True)
        return
    if SYSTEM == "windows":
        run("cmake", "--install", "build", "--prefix", "install")
        wix = Path(r"C:\Program Files (x86)\WiX Toolset v3.11\bin")
        candle = shutil.which("candle.exe") or str(wix / "candle.exe")
        light = shutil.which("light.exe") or str(wix / "light.exe")
        require_file(Path(candle))
        require_file(Path(light))
        run(
            sys.executable,
            "Resources/Packages/windows.py",
            "--version",
            version,
            "--arch",
            "x64",
            "--output",
            f"OpenScofo-{version}-Win64.msi",
            "--candle",
            candle,
            "--light",
            light,
        )
        return
    raise RuntimeError(f"Unsupported operating system: {SYSTEM}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def smoke_test_installer() -> None:
    require_platform("linux")
    run_id = os.environ.get("GITHUB_RUN_ID", str(os.getpid()))
    smoke_home = RUNNER_TEMP / f"openscofo-installer-smoke-{run_id}"
    install_root = smoke_home / ".local/share/openscofo-installer"
    examples = smoke_home / ".local/share/openscofo/examples"
    vamp = smoke_home / ".vamp"
    examples.mkdir(parents=True, exist_ok=True)
    vamp.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / "LICENSE", examples / "canticos.scofo")
    shutil.copy2(ROOT / "README.md", vamp / "unrelated-user-file.txt")
    sentinel_hash = sha256(examples / "canticos.scofo")
    env = os.environ.copy()
    env.update(HOME=str(smoke_home), QT_QPA_PLATFORM="offscreen")
    components = (
        "io.github.openscofo.csound",
        "io.github.openscofo.puredata",
        "io.github.openscofo.supercollider",
        "io.github.openscofo.vamp",
        "io.github.openscofo.examples",
    )
    run(
        ROOT / "OpenScofo-Linux.run",
        "--root",
        install_root,
        "--accept-licenses",
        "--default-answer",
        "--confirm-command",
        "install",
        *components,
        env=env,
    )
    pd_root = smoke_home / "Documents/Pd/externals/openscofo~"
    required = (
        smoke_home / ".local/lib/csound/6.0/plugins64/OpenScofo.so",
        pd_root / "flute.onnx",
        pd_root / "miniatura1.scofo",
        pd_root / "miniatura1.wav",
        pd_root / "openscofo.train-help.pd",
        smoke_home / ".local/share/SuperCollider/Extensions/OpenScofo/OpenScofo.so",
        vamp / "OpenScofo.so",
        examples / "bwv-1013.scofo",
    )
    for path in required:
        require_file(path)
    if not any(require_file(path) for path in pd_root.glob("openscofo~*.so")):
        raise FileNotFoundError("Installed Pure Data external was not found")
    if not any(path.is_file() for path in pd_root.glob("openscofo.train*.so")):
        raise FileNotFoundError("Installed openscofo.train external was not found")
    score = (pd_root / "miniatura1.scofo").read_text(encoding="utf-8")
    if not re.search(r"^ONNXMODEL\s+flute\.onnx(?:\s|$)", score, re.MULTILINE):
        raise ValueError("Installed score does not reference flute.onnx")
    if sentinel_hash != sha256(examples / "canticos.scofo"):
        raise RuntimeError("Installer modified an existing example file")
    maintenance_tool = require_file(
        install_root / "OpenScofoMaintenanceTool", executable=True
    )
    run(
        maintenance_tool,
        "--accept-licenses",
        "--default-answer",
        "--confirm-command",
        "remove",
        *components,
        env=env,
    )
    removed = (
        smoke_home / ".local/lib/csound/6.0/plugins64/OpenScofo.so",
        smoke_home / ".local/share/SuperCollider/Extensions/OpenScofo/OpenScofo.so",
        vamp / "OpenScofo.so",
    )
    if any(path.exists() for path in removed) or any(
        pd_root.glob("openscofo.train*.so")
    ):
        raise RuntimeError("The component uninstaller left installed binaries behind")
    require_file(examples / "canticos.scofo")
    require_file(vamp / "unrelated-user-file.txt")


def stage_max_package() -> None:
    require_platform("macos", "windows")
    package_root = ROOT / "max-package/OpenScofo"
    (package_root / "externals").mkdir(parents=True, exist_ok=True)
    if SYSTEM == "windows":
        copy_file(
            "max/openscofo~.mxe64", "max-package/OpenScofo/externals/openscofo~.mxe64"
        )
        copy_file(
            "max/openscofo.train.mxe64",
            "max-package/OpenScofo/externals/openscofo.train.mxe64",
        )
        copy_tree("Sources/Wrappers/Max/package/extras", "max-package/OpenScofo/extras")
        require_file(package_root / "externals/openscofo~.mxe64")
        require_file(package_root / "externals/openscofo.train.mxe64")
        require_file(package_root / "extras/OpenScofo.maxpat")
        return

    for directory in ("docs", "help"):
        (package_root / directory).mkdir(parents=True, exist_ok=True)
    for filename in ("License.md", "Readme.md", "icon.png", "package-info.json"):
        copy_file(f"Sources/Wrappers/Max/package/{filename}", package_root / filename)
    copy_tree("Sources/Wrappers/Max/package/extras", "max-package/OpenScofo/extras")
    files = {
        "Sources/Wrappers/Max/docs/openscofo~.maxref.xml": "docs/openscofo~.maxref.xml",
        "Sources/Wrappers/Max/openscofo~.maxhelp": "help/openscofo~.maxhelp",
        "Sources/Trainer/Max/openscofo.train.maxhelp": "help/openscofo.train.maxhelp",
        "Sources/Wrappers/Max/score-max.svg": "help/score-max.svg",
        "Sources/Wrappers/Max/score_renderer.js": "help/score_renderer.js",
        "Tests/flute.onnx": "help/flute.onnx",
        "Tests/01-benchmark/real/canticos.scofo": "help/canticos.scofo",
        "Tests/01-benchmark/real/canticos.mp3": "help/canticos.mp3",
        "build/OpenScofoPackageAssets/miniatura1.wav": "help/miniatura1.wav",
        "Sources/Wrappers/Max/miniatura1-max.scofo": "help/miniatura1.scofo",
    }
    for source, destination in files.items():
        copy_file(source, package_root / destination)
    copy_tree("max/openscofo~.mxo", "max-package/OpenScofo/externals/openscofo~.mxo")
    copy_tree(
        "max/openscofo.train.mxo",
        "max-package/OpenScofo/externals/openscofo.train.mxo",
    )
    for destination in (
        *files.values(),
        "extras/OpenScofo.maxpat",
        "externals/openscofo~.mxo/Contents/MacOS/openscofo~",
        "externals/openscofo.train.mxo/Contents/MacOS/openscofo.train",
    ):
        require_file(package_root / destination)


def stage_puredata_package() -> None:
    package_root = ROOT / "pd-package" / LIBNAME
    package_root.mkdir(parents=True, exist_ok=True)
    existing = ROOT / "build" / LIBNAME
    if not existing.is_dir():
        existing = ROOT / LIBNAME
    if existing.is_dir():
        shutil.copytree(existing, package_root, dirs_exist_ok=True)
    for stale_cut in package_root.glob("openscofo.cut*"):
        if stale_cut.is_dir():
            shutil.rmtree(stale_cut)
        else:
            stale_cut.unlink()
    for source in (
        "Sources/Wrappers/PureData/openscofo~-help.pd",
        "Sources/Trainer/PureData/openscofo.train-help.pd",
        "README.md",
        "LICENSE",
        "Documentation/assets/icon.png",
        "Tests/flute.onnx",
        "build/OpenScofoPackageAssets/miniatura1.scofo",
        "build/OpenScofoPackageAssets/miniatura1.wav",
    ):
        copy_file(source, package_root / Path(source).name)

    pd_binaries = [
        path
        for path in (ROOT / "build/Sources/Wrappers/PureData").glob("openscofo~*")
        if path.is_file() and path.suffix != ".pd"
    ]
    trainer_patterns = (
        "build/Sources/Trainer/openscofo.train*",
        "build/OpenScofo/openscofo.train*",
    )
    trainers = [
        path
        for pattern in trainer_patterns
        for path in ROOT.glob(pattern)
        if path.is_file()
        and path.suffix != ".pd"
        and (SYSTEM != "windows" or path.name.endswith(".dll"))
    ]
    for source in sorted(set(pd_binaries + trainers)):
        shutil.copy2(source, package_root / source.name)
    staged_pd = [
        path for path in package_root.glob("openscofo~*") if path.suffix != ".pd"
    ]
    staged_trainers = [
        path for path in package_root.glob("openscofo.train*") if path.suffix != ".pd"
    ]
    if not staged_pd:
        raise FileNotFoundError("No Pure Data external binary was staged")
    if not staged_trainers:
        raise FileNotFoundError("No openscofo.train Pure Data external was staged")
    for path in sorted(package_root.iterdir()):
        if path.is_file():
            print(path.relative_to(ROOT))


def publish_pypi() -> None:
    wheels = sorted((ROOT / "dist").glob("*.whl"))
    if not wheels:
        raise FileNotFoundError("No wheel found in dist/")
    run(sys.executable, "-m", "twine", "upload", "--skip-existing", *wheels)


Action = Callable[[argparse.Namespace], None]


def no_args(function: Callable[[], None]) -> Action:
    return lambda _arguments: function()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--precision",
        type=int,
        choices=(32, 64),
        default=32,
        help="Pure Data floating-point precision used by --build (default: 32)",
    )
    actions = parser.add_mutually_exclusive_group(required=True)
    operations: tuple[tuple[tuple[str, ...], Action, str], ...] = (
        (("--reset-sccache",), no_args(reset_sccache), "Reset sccache statistics"),
        (
            ("--show-sccache-stats",),
            no_args(show_sccache_stats),
            "Print sccache statistics",
        ),
        (
            ("--install-dependencies",),
            no_args(install_dependencies),
            "Install dependencies for this OS",
        ),
        (
            ("--install-qt-framework", "--install-qtinstallframwork"),
            no_args(install_qt_framework),
            "Install Qt Installer Framework on Linux",
        ),
        (
            ("--build",),
            lambda args: build(args.precision),
            "Build the native wheel for this OS",
        ),
        (("--run-native-tests",), no_args(run_native_tests), "Run native CTest tests"),
        (
            ("--verify-native-trainer",),
            no_args(verify_native_trainer),
            "Verify the trainer binary",
        ),
        (
            ("--repair-linux-wheel",),
            no_args(repair_linux_wheel),
            "Repair the Linux wheel",
        ),
        (
            ("--build-emscripten",),
            no_args(build_emscripten),
            "Build the Emscripten target",
        ),
        (
            ("--extract-wheel-version",),
            no_args(extract_wheel_version),
            "Export the wheel name and version",
        ),
        (
            ("--generate-installer",),
            no_args(generate_installer),
            "Generate the installer for this OS",
        ),
        (
            ("--files-test",),
            no_args(files_test),
            "Validate installer input or payload files",
        ),
        (
            ("--smoke-test-installer",),
            no_args(smoke_test_installer),
            "Install and remove the Linux installer",
        ),
        (
            ("--stage-max-package",),
            no_args(stage_max_package),
            "Stage the Max package for this OS",
        ),
        (
            ("--stage-puredata-package",),
            no_args(stage_puredata_package),
            "Stage the Pure Data package",
        ),
        (("--publish-pypi",), no_args(publish_pypi), "Upload wheels to PyPI"),
    )
    for flags, function, help_text in operations:
        actions.add_argument(
            *flags, dest="action", action="store_const", const=function, help=help_text
        )
    return parser.parse_args()


def main() -> int:
    os.chdir(ROOT)
    arguments = parse_args()
    print(f"OpenScofo CI: action={arguments.action.__name__}, os={SYSTEM}", flush=True)
    arguments.action(arguments)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
