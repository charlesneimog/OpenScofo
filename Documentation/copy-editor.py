import shutil
from pathlib import Path

def on_post_build(config):
    # Root do projeto (onde está o mkdocs.yml)
    root = Path(config.config_file_path).parent.resolve()

    # Diretório de saída (site ou tmp no serve)
    site_dir = Path(config["site_dir"]).resolve()

    origem = root / "Resources" / "Online-Editor"
    destino = site_dir / "Editor"

    if destino.exists():
        shutil.rmtree(destino)

    shutil.copytree(origem, destino)
