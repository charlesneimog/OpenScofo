import shutil
from pathlib import Path


def on_pre_build(config):
    # Root do projeto (onde está o mkdocs.yml)
    root = Path(config.config_file_path).parent.resolve()

    # docs_dir definido no mkdocs.yml (no seu caso: Documentation)
    docs_dir = root / config["docs_dir"]

    origem = root / "Resources" / "Online-Editor"
    destino = docs_dir / "Editor"

    if destino.exists():
        shutil.rmtree(destino)

    shutil.copytree(origem, destino)
