"""
Train and smoke-test an OpenScofo extended-technique model.

The test uses the local flute extended-technique dataset when it is available:

    /home/neimog/Nextcloud/Flute

It trains through OpenScofo.ExtendedTechniqueClassifier, exports CatBoost to
ONNX, then loads that ONNX model back into the OpenScofo engine and performs
window-level inference on real audio.
"""

from pathlib import Path
import argparse
import os
import tempfile
import unittest

os.environ.setdefault("NUMBA_CACHE_DIR", tempfile.gettempdir())
os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

import numpy as np
from catboost import CatBoostClassifier
from scipy.io import wavfile
from scipy.signal import resample_poly

import OpenScofo


TRAIN_FOLDER = Path("/home/neimog/Nextcloud/Flute")
OUTPUT_MODEL = Path("Tests/extended-techniques.onnx")
SAMPLE_RATE = 48000
FFT_SIZE = 2048
HOP_SIZE = 512

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


def _descriptor_enums():
    return [getattr(OpenScofo.Descriptors, name.upper()) for name in DESCRIPTORS]


def _audio_files_by_class(dataset):
    return {
        folder.name: sorted(folder.glob("*.wav"))
        for folder in sorted(dataset.iterdir())
        if folder.is_dir()
    }


def _load_audio_file(path, sample_rate=SAMPLE_RATE):
    sr, data = wavfile.read(path)
    data = np.asarray(data)

    if data.ndim > 1:
        data = np.mean(data, axis=1)

    if np.issubdtype(data.dtype, np.integer):
        data = data.astype(np.float32) / np.iinfo(data.dtype).max
    else:
        data = data.astype(np.float32)

    if sr != sample_rate:
        gcd = np.gcd(sr, sample_rate)
        data = resample_poly(data, sample_rate // gcd, sr // gcd).astype(np.float32)

    return data


def _build_trainer(base_path, quiet=True):
    trainer = OpenScofo.ExtendedTechniqueClassifier(
        sample_rate=SAMPLE_RATE,
        fft_size=FFT_SIZE,
        hop_size=HOP_SIZE,
        base_path=str(base_path),
        model_type="catboost",
    )
    if quiet:
        trainer.set_print_callback(lambda *_args, **_kwargs: None)
    trainer.set_random_seed(42)
    trainer.set_catboost_config(
        iterations=500,
        learning_rate=0.05,
        depth=6,
        early_stopping_rounds=60,
        thread_count=1,
    )
    trainer.test_fraction = 0.2
    trainer.min_db = -60.0
    trainer.max_augmentation_loops = 128
    trainer.set_balanced_augmentation_factor(1.0)
    trainer.single_file_class_window_split = True
    trainer.set_descriptors(DESCRIPTORS)
    trainer.set_train_folder(str(TRAIN_FOLDER))
    trainer._load_audio = _load_audio_file
    trainer._build_classifier = lambda: CatBoostClassifier(
        iterations=trainer.iterations,
        depth=trainer.depth,
        learning_rate=trainer.learning_rate,
        loss_function="MultiClass",
        thread_count=trainer.thread_count,
        early_stopping_rounds=trainer.early_stopping_rounds,
        random_seed=trainer.random_state,
        auto_class_weights="Balanced",
        allow_writing_files=False,
        verbose=not quiet,
    )
    return trainer


def _average_catboost_probabilities(trainer, wav_path):
    y = _load_audio_file(wav_path)
    features = []

    for idx in range(0, len(y) - FFT_SIZE + 1, HOP_SIZE):
        window = np.asarray(y[idx : idx + FFT_SIZE], dtype=np.float32)
        trainer.oscofo.process_block(window)
        desc = trainer.oscofo.get_description()
        if desc.db < trainer.min_db:
            continue
        features.append(trainer._descriptor_values(desc))

    if not features:
        return {}

    probabilities = trainer.clf.predict_proba(np.asarray(features))
    mean_probabilities = np.mean(probabilities, axis=0)
    return {
        str(label): float(value)
        for label, value in zip(trainer.clf.classes_, mean_probabilities)
    }


def _average_onnx_probabilities(model_path, wav_path):
    scofo = OpenScofo.OpenScofo(SAMPLE_RATE, FFT_SIZE, HOP_SIZE)
    scofo.request_descriptor(OpenScofo.Descriptors.ONNX)
    scofo.load_onnx_model(str(model_path), _descriptor_enums())

    y = _load_audio_file(wav_path)
    probabilities = []

    for idx in range(0, len(y) - FFT_SIZE + 1, HOP_SIZE):
        window = np.asarray(y[idx : idx + FFT_SIZE], dtype=np.float32)
        scofo.process_block(window)
        desc = scofo.get_description()
        if desc.db < -60.0:
            continue
        frame_probs = dict(desc.onnx)
        if frame_probs:
            probabilities.append(frame_probs)

    if not probabilities:
        return {}

    labels = sorted(probabilities[0])
    return {
        label: float(np.mean([frame[label] for frame in probabilities]))
        for label in labels
    }


def _evaluate_files(files_by_class, predict_probabilities):
    total = 0
    correct = 0
    per_class = {}

    for label, files in files_by_class.items():
        class_total = 0
        class_correct = 0
        examples = []

        for wav_path in files:
            probs = predict_probabilities(wav_path)
            if not probs:
                continue

            predicted = max(probs, key=probs.get)
            is_correct = predicted == label
            class_total += 1
            class_correct += int(is_correct)
            examples.append((wav_path, predicted, probs))

        total += class_total
        correct += class_correct
        per_class[label] = {
            "total": class_total,
            "correct": class_correct,
            "examples": examples,
        }

    accuracy = correct / total if total else 0.0
    return correct, total, accuracy, per_class


def _ensure_training_dataset():
    if not TRAIN_FOLDER.exists():
        raise FileNotFoundError(f"Training dataset not found: {TRAIN_FOLDER}")

    files_by_class = {
        label: files
        for label, files in _audio_files_by_class(TRAIN_FOLDER).items()
        if files
    }
    if len(files_by_class) < 2:
        raise RuntimeError("Need at least two non-empty class folders")
    return files_by_class


def train_and_show(output_model=OUTPUT_MODEL):
    files_by_class = _ensure_training_dataset()
    output_model = Path(output_model)
    output_model.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="openscofo-train-") as tmp:
        trainer = _build_trainer(Path(tmp), quiet=False)

        print(f"Dataset: {TRAIN_FOLDER}")
        for label, files in files_by_class.items():
            print(f"  {label}: {len(files)} wav files")

        print("\nExtracting OpenScofo descriptors...")
        trainer.analyze()
        trainer.print_data()

        print("\nTraining CatBoost classifier...")
        trainer.train()

        exported = Path(tmp) / output_model.name
        trainer.export_model(exported.name)
        output_model.write_bytes(exported.read_bytes())

        print(f"\nModel written to: {output_model}")
        print("\nCatBoost inference on every wav file:")
        _print_evaluation(
            *_evaluate_files(
                files_by_class,
                lambda wav_path: _average_catboost_probabilities(trainer, wav_path),
            )
        )

        print("\nOpenScofo ONNX inference on every wav file:")
        _print_evaluation(
            *_evaluate_files(
                files_by_class,
                lambda wav_path: _average_onnx_probabilities(output_model, wav_path),
            )
        )


def _format_probabilities(probs):
    return ", ".join(
        f"{label}={value:.3f}"
        for label, value in sorted(
            probs.items(), key=lambda item: item[1], reverse=True
        )
    )


def _print_evaluation(correct, total, accuracy, per_class):
    print(f"  Accuracy: {accuracy:.3f} ({correct}/{total})")
    for label, result in per_class.items():
        class_total = result["total"]
        class_correct = result["correct"]
        class_accuracy = class_correct / class_total if class_total else 0.0
        print(f"  {label}: {class_accuracy:.3f} ({class_correct}/{class_total})")

        for wav_path, predicted, probs in result["examples"][:3]:
            print(f"    {wav_path.name} -> {predicted}")
            print(f"      {_format_probabilities(probs)}")


def main():
    parser = argparse.ArgumentParser(
        description="Train an OpenScofo extended-technique model and show inference."
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="run the unittest smoke test instead of the visible training demo",
    )
    parser.add_argument(
        "--output",
        default=str(OUTPUT_MODEL),
        help=f"ONNX model output path, default: {OUTPUT_MODEL}",
    )
    args = parser.parse_args()

    if args.test:
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(
            ExtendedTechniqueModelTest
        )
        raise SystemExit(unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() is False)

    train_and_show(args.output)


class ExtendedTechniqueModelTest(unittest.TestCase):
    def setUp(self):
        if not TRAIN_FOLDER.exists():
            self.skipTest(f"Training dataset not found: {TRAIN_FOLDER}")

        self.files_by_class = _audio_files_by_class(TRAIN_FOLDER)
        self.files_by_class = {
            label: files for label, files in self.files_by_class.items() if files
        }
        if len(self.files_by_class) < 2:
            self.skipTest("Need at least two non-empty class folders")

    def test_train_export_and_infer_with_openscofo_model(self):
        with tempfile.TemporaryDirectory(prefix="openscofo-train-test-") as tmp:
            tmp_path = Path(tmp)
            trainer = _build_trainer(tmp_path)

            trainer.analyze()
            self.assertGreater(len(trainer.x_np_train), 0)
            self.assertGreaterEqual(len(set(trainer.y_np_train)), 2)

            trainer.train()
            model_path = tmp_path / "extended-techniques-test.onnx"
            trainer.export_model(model_path.name)

            self.assertTrue(model_path.exists())
            self.assertGreater(model_path.stat().st_size, 0)

            expected_labels = set(self.files_by_class)
            correct = 0

            for label, files in self.files_by_class.items():
                probs = _average_onnx_probabilities(model_path, files[0])
                self.assertEqual(set(probs), expected_labels)
                self.assertTrue(all(np.isfinite(value) for value in probs.values()))

                predicted_label = max(probs, key=probs.get)
                correct += int(predicted_label == label)

            self.assertGreaterEqual(correct, 1)


if __name__ == "__main__":
    main()
