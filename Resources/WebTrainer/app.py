"""Train OpenScofo models from musician-friendly ZIP datasets."""

import hashlib
import logging
import multiprocessing
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import tempfile
import zipfile

import gradio as gr

WORK_DIR = Path(tempfile.gettempdir()) / "openscofo-web-trainer"
AUDIO_EXTENSIONS = {".wav", ".aif", ".aiff"}
DESCRIPTORS = [
    "mfcc",
    "logmel",
    "centroid",
    "flatness",
    "hfr",
    "flux",
    "zcr",
    "irregularity",
    "kurtosis",
]
MAX_BYTES = 1024 * 1024 * 1024
MAX_ENTRIES = 10000


class DatasetError(ValueError):
    """A dataset problem that can be explained to the user."""


def extract_dataset(archive_path, destination):
    """Extract only audio and visible directories, then detect one optional wrapper."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(archive_path) as archive:
            entries = archive.infolist()
            if (
                len(entries) > MAX_ENTRIES
                or sum(e.file_size for e in entries) > MAX_BYTES
            ):
                raise DatasetError(
                    "The ZIP is too large. Use at most 10,000 entries and 1 GiB of uncompressed files."
                )
            targets = set()
            for entry in entries:
                name = entry.filename
                path = PurePosixPath(name)
                mode = entry.external_attr >> 16
                if (
                    path.is_absolute()
                    or ".." in path.parts
                    or "\\" in name
                    or ":" in name
                    or any(ord(c) < 32 for c in name)
                    or stat.S_ISLNK(mode)
                    or (
                        stat.S_IFMT(mode)
                        and not (stat.S_ISREG(mode) or stat.S_ISDIR(mode))
                    )
                ):
                    raise DatasetError(
                        "The ZIP contains an unsafe path or special file. Please create a new ZIP from your class folders."
                    )
                if not path.parts or any(
                    p.startswith(".") or p == "__MACOSX" for p in path.parts
                ):
                    continue
                target = destination.joinpath(*path.parts)
                if not target.resolve().is_relative_to(destination.resolve()):
                    raise DatasetError("The ZIP contains an unsafe path.")
                if entry.is_dir():
                    target.mkdir(parents=True, exist_ok=True)
                    continue
                # Preserve folders even when they contain only unsupported files:
                # these must be reported as empty classes rather than disappearing.
                target.parent.mkdir(parents=True, exist_ok=True)
                if path.suffix.lower() not in AUDIO_EXTENSIONS:
                    continue
                # OpenScofo currently selects lowercase extensions. No audio conversion
                # is needed; normalize just the suffix of uppercase upload filenames.
                target = target.with_suffix(path.suffix.lower())
                if target in targets:
                    raise DatasetError(
                        "The ZIP contains duplicate audio filenames. Give each recording a unique filename within its class."
                    )
                targets.add(target)
                with archive.open(entry) as source, target.open("xb") as output:
                    shutil.copyfileobj(source, output)
    except DatasetError:
        raise
    except (
        zipfile.BadZipFile,
        zipfile.LargeZipFile,
        OSError,
        RuntimeError,
        EOFError,
        ValueError,
        NotImplementedError,
    ) as exc:
        raise DatasetError(
            "Cannot read this ZIP. Upload a valid, unencrypted ZIP containing your recordings."
        ) from exc

    root = destination
    children = list(root.iterdir())
    if len(children) == 1 and children[0].is_dir():
        candidate = children[0]
        if any(p.is_dir() for p in candidate.iterdir()):
            root = candidate
    classes = sorted(p for p in root.iterdir() if p.is_dir())
    if any(p.is_file() for p in root.iterdir()) or any(
        p.is_dir() for folder in classes for p in folder.iterdir()
    ):
        raise DatasetError(
            "Use one folder per technique, with recordings directly inside it. One outer instrument folder is optional; nested folders and audio outside class folders are not supported."
        )
    counts = [
        (folder.name, sum(p.is_file() for p in folder.iterdir())) for folder in classes
    ]
    if not any(count for _, count in counts):
        raise DatasetError(
            "No supported audio files found. Put WAV, AIF, or AIFF recordings inside technique folders."
        )
    if len(counts) < 2:
        raise DatasetError(
            "At least two technique folders are required to train a classifier."
        )
    empty = [label for label, count in counts if not count]
    if empty:
        raise DatasetError(
            "Every technique folder needs at least one WAV, AIF, or AIFF recording. Remove or fill the empty folders shown in your ZIP."
        )
    return root, counts


def archive_digest(path):
    with open(path, "rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def preview_dataset(archive_path):
    """The state records exactly which archive and labels the user reviewed."""
    if not archive_path:
        return (
            [],
            "Upload a ZIP to see its technique labels.",
            None,
            gr.Button(interactive=False),
            None,
            "",
        )
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(prefix="preview-", dir=WORK_DIR) as work:
            _, counts = extract_dataset(archive_path, Path(work) / "dataset")
            reviewed = {"digest": archive_digest(archive_path), "counts": counts}
        total = sum(count for _, count in counts)
        return (
            counts,
            f"{total} audio files total. Check these labels before training.",
            reviewed,
            gr.Button(interactive=True),
            None,
            "",
        )
    except DatasetError as exc:
        return [], str(exc), None, gr.Button(interactive=False), None, ""
    except Exception:
        logging.exception("Dataset preview failed")
        return (
            [],
            "Could not inspect the dataset. Please upload the ZIP again.",
            None,
            gr.Button(interactive=False),
            None,
            "",
        )


def training_worker(connection, work, dataset_root, descriptors, config):
    """Keep library RNG state and CatBoost working files local to each process."""
    try:
        import OpenScofo

        os.chdir(work)
        trainer = OpenScofo.ExtendedTechniqueClassifier(
            sample_rate=48000,
            fft_size=2048,
            hop_size=512,
            base_path=work,
            model_type="catboost",
        )
        trainer.set_print_callback(lambda message: None)
        trainer.set_onprogress_callback(
            lambda count: (
                connection.send(("progress", count)) if count % 500 == 0 else None
            )
        )
        trainer.set_descriptors(descriptors)
        trainer.set_catboost_config(**config)
        trainer.set_train_folder(dataset_root)
        connection.send(
            ("status", "Analyzing recordings and preparing training samples…")
        )
        trainer.analyze()
        # A file may load successfully but contain only silence or be too short.
        if set(map(str, trainer.y_np_train)) != set(trainer.folders):
            raise DatasetError(
                "Some techniques have no usable audio samples. Use audible recordings longer than 0.05 seconds in every folder."
            )
        connection.send(
            ("status", "Training the model… This may take several minutes.")
        )
        trainer.train()
        connection.send(("status", "Exporting the OpenScofo model…"))
        trainer.export_model("extended-techniques.onnx")
        connection.send(("done", None))
    except DatasetError as exc:
        connection.send(("error", str(exc)))
    except Exception:
        logging.exception("OpenScofo training failed")
        connection.send(
            (
                "error",
                "Training could not finish. Check that every recording is readable and contains audible sound, then try again. If it continues, ask the server administrator to check the training log.",
            )
        )
    finally:
        connection.close()


def train_model(
    archive_path,
    reviewed,
    descriptors,
    iterations,
    learning_rate,
    depth,
    early_stopping,
    progress=gr.Progress(),
):
    """Revalidate a reviewed archive and stream progress from an isolated job."""
    yield "Preparing your dataset…", None
    if not archive_path or not reviewed:
        yield "Upload a valid ZIP and inspect its detected classes first.", None
        return
    if not descriptors or any(d not in DESCRIPTORS for d in descriptors):
        yield "Select at least one descriptor in Advanced Settings.", None
        return
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    try:
        config = {
            "iterations": int(iterations),
            "learning_rate": float(learning_rate),
            "depth": int(depth),
            "early_stopping_rounds": int(early_stopping),
            "thread_count": max(1, min(4, os.cpu_count() or 1)),
        }
        if not (
            1 <= config["iterations"] <= 10000
            and 0 < config["learning_rate"] <= 1
            and 1 <= config["depth"] <= 10
            and 1 <= config["early_stopping_rounds"] <= 1000
        ):
            raise DatasetError(
                "Choose training settings within the ranges shown in Advanced Settings."
            )
        with tempfile.TemporaryDirectory(prefix="job-", dir=WORK_DIR) as work:
            # A private snapshot makes the review check and extraction use the same bytes.
            snapshot = Path(work) / "dataset.zip"
            shutil.copyfile(archive_path, snapshot)
            if archive_digest(snapshot) != reviewed["digest"]:
                raise DatasetError(
                    "The upload changed. Upload it again and check the detected labels before training."
                )
            root, counts = extract_dataset(snapshot, Path(work) / "dataset")
            if counts != reviewed["counts"]:
                raise DatasetError(
                    "The dataset changed. Upload it again to review the labels."
                )
            context = multiprocessing.get_context("spawn")
            receiver, sender = context.Pipe(duplex=False)
            worker = context.Process(
                target=training_worker,
                args=(sender, work, str(root), descriptors, config),
            )
            try:
                worker.start()
                sender.close()
                while True:
                    if receiver.poll(0.2):
                        try:
                            kind, value = receiver.recv()
                        except EOFError:
                            raise DatasetError(
                                "The training worker stopped unexpectedly. The server may need more memory; try a smaller dataset."
                            )
                        if kind == "error":
                            raise DatasetError(value)
                        if kind == "done":
                            progress(1, desc="Model ready")
                            # Gradio copies yielded files into its managed download cache
                            # before resuming this generator and removing the job directory.
                            yield "Model ready. Download your OpenScofo model below.", str(
                                Path(work) / "extended-techniques.onnx"
                            )
                            break
                        message = (
                            f"Preparing training samples… {value:,} samples extracted."
                            if kind == "progress"
                            else value
                        )
                        progress(None, desc=message)
                        yield message, None
                    elif not worker.is_alive():
                        raise DatasetError(
                            "The training worker stopped unexpectedly. Try a smaller dataset or ask the server administrator to check available memory."
                        )
            finally:
                if worker.pid is not None:
                    if worker.is_alive():
                        worker.terminate()
                    worker.join()
                receiver.close()
                sender.close()
    except DatasetError as exc:
        yield str(exc), None
    except Exception:
        logging.exception("Training job failed")
        yield "Could not start training. Please upload the ZIP again or ask the server administrator to check the training log.", None


with gr.Blocks(title="OpenScofo Model Trainer", delete_cache=(3600, 86400)) as app:
    gr.Markdown(
        "# OpenScofo Model Trainer\nCreate an extended-technique classifier from your recordings, then download a model to use in OpenScofo."
    )
    reviewed = gr.State()
    gr.Markdown(
        "## 1. Upload Dataset\nUpload a ZIP with one folder per technique. WAV, AIF, and AIFF recordings are supported."
    )
    upload = gr.File(label="Dataset ZIP", file_types=[".zip"], type="filepath")
    gr.Markdown("## 2. Dataset")
    dataset = gr.Dataframe(
        headers=["Technique", "Recordings"],
        datatype=["str", "number"],
        interactive=False,
    )
    summary = gr.Textbox(
        value="Upload a ZIP to see its technique labels.",
        show_label=False,
        interactive=False,
    )
    gr.Markdown("## 3. Train")
    with gr.Accordion("Advanced Settings", open=False):
        descriptors = gr.CheckboxGroup(
            DESCRIPTORS, value=DESCRIPTORS, label="Audio descriptors"
        )
        iterations = gr.Slider(1, 10000, value=1000, step=1, label="Iterations")
        learning_rate = gr.Slider(
            0.001, 1, value=0.05, step=0.001, label="Learning rate"
        )
        depth = gr.Slider(1, 10, value=6, step=1, label="Tree depth")
        early_stopping = gr.Slider(
            1, 1000, value=70, step=1, label="Early stopping rounds"
        )
    train = gr.Button("Train Model", variant="primary", interactive=False)
    status = gr.Textbox(label="Training progress", interactive=False)
    download = gr.File(label="Download OpenScofo Model", interactive=False)
    upload.change(
        preview_dataset, upload, [dataset, summary, reviewed, train, download, status]
    )
    train.click(
        train_model,
        [
            upload,
            reviewed,
            descriptors,
            iterations,
            learning_rate,
            depth,
            early_stopping,
        ],
        [status, download],
        concurrency_limit=2,
    )


if __name__ == "__main__":
    app.queue()
    app.launch(
        server_name="0.0.0.0",
        server_port=int(os.environ.get("PORT", 7860)),
    )
