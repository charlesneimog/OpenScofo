# OpenScofo Model Trainer

A Gradio web interface for musicians to upload recordings, review technique labels,
train an extended-technique classifier, and download `extended-techniques.onnx`.
Extraction, augmentation, CatBoost training, and ONNX export use OpenScofo's
`OpenScofo.ExtendedTechniqueClassifier`; no separate training implementation is used.

## Dataset

Upload an unencrypted ZIP containing at least two technique folders:

```text
Flute.zip
└── Flute/                 # optional outer folder
    ├── normal/
    │   ├── normal01.wav
    │   └── normal02.aiff
    └── key-click/
        ├── click01.aif
        └── click02.wav
```

The class folders can also be directly at the ZIP root. Folder names become model
labels. Recordings must be directly inside each class folder; nested folders are
rejected. WAV, AIF, and AIFF extensions are accepted, including uppercase extensions.
Hidden files, `__MACOSX`, README files, and other unsupported files are ignored.
Every visible class folder must contain supported audio. Use audible recordings
longer than 0.05 seconds, preferably several independent recordings per technique.
Archives are limited to 1 GiB uncompressed and 10,000 entries; unsafe paths, links,
and duplicate audio paths are rejected.

After upload, check the technique table, then click **Train Model**. Defaults use
all nine descriptors and the library's standard CatBoost configuration.
**Advanced Settings** exposes descriptors, iterations, learning rate, tree depth,
and early stopping. Analysis uses the library's default training sample rates
(44.1 and 48 kHz), with the standard 48 kHz runtime configuration, FFT size 2048,
and hop size 512. Early stopping depends on usable validation samples; one-file
classes cannot provide an independent recording for test evaluation.

## Run locally

Use Python 3.12–3.14 (3.13 is recommended) and a virtual environment. Install
OpenScofo and Gradio with pip. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r Resources/WebTrainer/requirements.txt
python Resources/WebTrainer/app.py
```

Open `http://localhost:7860`. The requirements install the published OpenScofo
package and Gradio 6. OpenScofo supplies its audio and training dependencies.
On platforms with a compatible wheel, no native compilation is needed.

For development against the current repository instead of the published package,
install it into the same environment:

```bash
export CMAKE_ARGS="-DOPENSCOFO_BUILD_ALL=OFF -DOPENSCOFO_BUILD_PY_MODULE=ON -DOPENSCOFO_BUILD_TESTS=OFF -DOPENSCOFO_UPDATE_LANGUAGE=OFF"
python -m pip install -e .
```

Only the source installation needs a C/C++ compiler, CMake 3.30+, Git, and network
access for native build dependencies. The published Linux wheel tested with this
app is OpenScofo 0.2.10b1. Its ONNX exporter does not add the OpenScofo descriptor
and label metadata added by the current repository; use the source installation
when you need that metadata.

## Railway

Connect this repository and keep the **repository root** as the build root.
Select Python 3.13. Set the following build and start commands explicitly so that
Railway installs the nested web requirements into the runtime environment.
The pip wheel installation does not need `CMAKE_ARGS`.

Build command:

```bash
python -m pip install -r Resources/WebTrainer/requirements.txt
```

Start command:

```bash
python Resources/WebTrainer/app.py
```

The server binds to `0.0.0.0` and Railway's `PORT` (locally, 7860).
No Gradio public sharing tunnel is enabled.

If the runtime reports `ModuleNotFoundError: No module named 'gradio'`, check that
the build log shows Gradio being installed from the requirements file above,
then rebuild and redeploy. Installing OpenScofo alone does not install Gradio.
Use `python -m pip` to install into the same Python environment used to start
the application.

Two jobs can train at once; each has a separate process and temporary directory,
including dataset caches and CatBoost working files. Temporary datasets are
removed after preview/training. Gradio caches the downloadable model before
job cleanup; cached downloads are removed after 24 hours, checked hourly.
Downloads are temporary and should be saved before restarting the server.
Training time and memory use depend on recording duration and class balance;
the ZIP size limit is not a bound on descriptor memory use.

Run dataset validation tests from the repository root:

```bash
python -m unittest discover -s Resources/WebTrainer/tests -v
```
