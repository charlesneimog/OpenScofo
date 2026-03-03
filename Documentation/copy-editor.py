import shutil
import tarfile
import io
from pathlib import Path
import requests

TREE_SITTER_VERSION = "v0.26.6"


def on_serve(server, config, builder):
    original_app = server.application

    def cors_middleware(environ, start_response):
        def custom_start_response(status, headers, exc_info=None):
            headers.append(("Access-Control-Allow-Origin", "*"))
            headers.append(("Access-Control-Allow-Methods", "GET, OPTIONS"))
            headers.append(("Access-Control-Allow-Headers", "*"))
            return start_response(status, headers, exc_info)

        return original_app(environ, custom_start_response)

    server.set_app(cors_middleware)
    return server


def on_post_build(config):
    root = Path(config.config_file_path).parent.resolve()
    site_dir = Path(config["site_dir"]).resolve()

    origem = root / "Score-Editor"
    destino = site_dir / "Editor"

    if destino.exists():
        shutil.rmtree(destino)

    shutil.copytree(origem, destino)

    # ─────────────────────────────────────
    # Download and extract web-tree-sitter files
    # ─────────────────────────────────────

    url = f"https://github.com/tree-sitter/tree-sitter/releases/download/{TREE_SITTER_VERSION}/web-tree-sitter.tar.gz"

    response = requests.get(url)
    response.raise_for_status()

    required_files = {
        "web-tree-sitter.wasm",
        "web-tree-sitter.js",
    }

    with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz") as tar:
        members = tar.getmembers()

        for filename in required_files:
            member = next((m for m in members if m.name.endswith(filename)), None)

            if member is None:
                raise FileNotFoundError(f"{filename} not found in archive")

            extracted = tar.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Failed to extract {filename}")

            output_path = destino / filename

            with open(output_path, "wb") as f:
                f.write(extracted.read())
