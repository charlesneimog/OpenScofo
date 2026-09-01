"""
Copyright (c) 2024-2026 Charles K. Neimog
Website: charlesneimog.github.io

This file is part of a project licensed under the
GNU General Public License v3.0 or later (GPL-3.0-or-later).
See the LICENSE file for details.
"""

import math
import os
import random
import json
import librosa
from collections.abc import Iterable
from importlib import metadata as importlib_metadata
from typing import List

import numpy as np
from catboost import CatBoostClassifier
from scipy.signal import fftconvolve, resample
from ._OpenScofo import OpenScofo


class ExtendedTechniqueClassifier:
    SUPPORTED_MODEL_TYPES = ("catboost",)
    DATASET_CACHE_VERSION = 13

    def __init__(
        self,
        sample_rate=48000,
        training_sample_rates=(44100, 48000),
        fft_size=2048,
        hop_size=512,
        base_path=None,
        model_type="catboost",
    ):
        self.print = print

        self.sample_rate = sample_rate
        self.training_sample_rates = self._normalize_sample_rates(training_sample_rates)
        self.fft_size = fft_size
        self.hop_size = hop_size
        self.base_path = base_path or os.getcwd()

        self.tabname = "train"
        self.redraw_tab_after = 120
        self.processed_count = 0
        self.descriptors = []
        self.redraw = 0

        self.iterations = 1000
        self.random_state = 80
        self.test_fraction = 0.1
        self.min_db = -60.0
        self.oscofo: OpenScofo = OpenScofo(sample_rate, fft_size, hop_size)
        self._oscofos = {sample_rate: self.oscofo}
        self.max_augmentation_loops = 128
        self.balanced_augmentation_factor = 1.0
        self.single_file_class_window_split = False
        self.learning_rate = 0.05
        self.depth = 6
        self.early_stopping_rounds = 70
        self.thread_count = -1

        self.trainfolder = ""
        self.folders = {}
        self.currtraindata = []
        self.currtestdata = []

        self.x_np_train, self.y_np_train = [], []
        self.x_np_val, self.y_np_val = [], []
        self.x_np_test, self.y_np_test = [], []
        self.x_train, self.y_train = [], []
        self.x_val, self.y_val = [], []
        self.x_test, self.y_test = [], []
        self.dataset_counts = {}

        self.set_model_type(model_type)
        self.clf = None

        self.ir_files = None
        self.ir_folders = []
        self.on_progress = None
        self._single_file_split_cache = {}
        self.set_random_seed(self.random_state)

    def _normalize_sample_rates(self, sample_rates):
        if sample_rates is None:
            sample_rates = (self.sample_rate,)
        if isinstance(sample_rates, (int, np.integer)):
            sample_rates = (int(sample_rates),)

        normalized = []
        for sample_rate in sample_rates:
            sample_rate = int(sample_rate)
            if sample_rate <= 0:
                raise ValueError("training sample rates must be positive integers")
            if sample_rate not in normalized:
                normalized.append(sample_rate)

        if not normalized:
            raise ValueError("at least one training sample rate is required")
        return tuple(normalized)

    def _oscofo_for_sample_rate(self, sample_rate):
        sample_rate = int(sample_rate)
        if sample_rate not in self._oscofos:
            self._oscofos[sample_rate] = OpenScofo(
                sample_rate, self.fft_size, self.hop_size
            )
        return self._oscofos[sample_rate]

    def resolve_trainfolder(self, path):
        """Resolve a training-related folder against base_path when needed."""
        if os.path.exists(path):
            return os.path.abspath(path)
        candidate = os.path.join(self.base_path, path)
        if not os.path.exists(candidate):
            raise FileNotFoundError(
                f"Folder not found: '{path}' (also tried '{candidate}')"
            )
        return os.path.abspath(candidate)

    def apply_reverb(self, y, sr):
        """Apply a randomly selected impulse response to a signal."""
        if not self.ir_files:
            return y

        ir_file = random.choice(self.ir_files)
        ir, _ = librosa.load(ir_file, sr=sr)
        ir = ir / (np.max(np.abs(ir)) + 1e-8)
        y_conv = fftconvolve(y, ir, mode="full")
        y_conv = y_conv[: len(y)]
        y_conv = y_conv / (np.max(np.abs(y_conv)) + 1e-8)
        return y_conv.astype(np.float32)

    def _load_audio(self, filepath, sample_rate=None):
        sample_rate = self.sample_rate if sample_rate is None else int(sample_rate)
        y, _ = librosa.load(filepath, sr=sample_rate)
        return y

    def _seed_random_generators(self):
        random.seed(self.random_state)
        np.random.seed(self.random_state)

    def set_random_seed(self, seed):
        """Set the seed used by Python, NumPy, and CatBoost randomness."""
        self.random_state = int(seed)
        self._seed_random_generators()

    def set_balanced_augmentation_factor(self, factor):
        """Set the class target multiplier used for train augmentation."""
        factor = float(factor)
        if factor < 1.0:
            raise ValueError("balanced_augmentation_factor must be >= 1")
        self.balanced_augmentation_factor = factor

    def _split_train_test(self, files):
        if not files:
            return [], [], False

        if len(files) == 1:
            if self.single_file_class_window_split:
                # A single recording cannot give independent file-level splits.
                # The duplicated path is safe because _process_file crops it into
                # non-overlapping temporal regions before any train augmentation.
                return files, files, True
            # With one recording, there is no independent test file. By default
            # we train on the whole file and omit this class from test metrics
            # rather than reporting a leaky or unrepresentative score.
            return files, [], False

        n_test = int(len(files) * self.test_fraction)
        n_test = max(1, n_test)
        n_test = min(n_test, len(files) - 1)

        test_files = random.sample(files, n_test)
        train_files = [f for f in files if f not in test_files]
        return train_files, test_files, False

    def _valid_window_starts(self, y, sample_rate):
        oscofo = self._oscofo_for_sample_rate(sample_rate)
        starts = []
        idx = 0
        while idx <= len(y) - self.fft_size:
            window = np.asarray(y[idx : idx + self.fft_size], dtype=np.float32)
            oscofo.process_block(window)
            desc = oscofo.get_description()
            if desc.db >= self.min_db:
                starts.append(idx)
            idx += self.hop_size
        return starts

    def _single_file_split_regions(self, y, filepath, label, sample_rate):
        cache_key = (os.path.abspath(filepath), label, int(sample_rate), len(y))
        if cache_key in self._single_file_split_cache:
            return self._single_file_split_cache[cache_key]

        empty = np.asarray([], dtype=np.float32)
        if len(y) < self.fft_size * 2:
            self.print(
                f"Class {label}: skipping temporal test split for "
                f"{os.path.basename(filepath)} because it is shorter than "
                f"two FFT windows."
            )
            regions = (y, empty)
            self._single_file_split_cache[cache_key] = regions
            return regions

        valid_starts = self._valid_window_starts(y, sample_rate)
        if len(valid_starts) < 2:
            self.print(
                f"Class {label}: skipping temporal test split for "
                f"{os.path.basename(filepath)} because it has fewer than two "
                f"non-silent windows."
            )
            regions = (y, empty)
            self._single_file_split_cache[cache_key] = regions
            return regions

        n_test = int(round(len(valid_starts) * self.test_fraction))
        n_test = max(1, min(n_test, len(valid_starts) - 1))
        split_at = valid_starts[-n_test]

        # Split at the first held-out non-silent window. This preserves temporal
        # order, avoids silent-only test tails, and prevents train augmentations
        # from being generated from any segment that can appear in test data.
        regions = (y[:split_at], y[split_at:])
        self._single_file_split_cache[cache_key] = regions
        return regions

    def _single_file_temporal_region(self, y, role, filepath, label, sample_rate):
        if role not in ("train", "test"):
            return y

        train_region, test_region = self._single_file_split_regions(
            y, filepath, label, sample_rate
        )
        if role == "train":
            return train_region
        return test_region

    def _count_valid_windows(self, y, sample_rate):
        oscofo = self._oscofo_for_sample_rate(sample_rate)
        count = 0
        idx = 0
        while idx <= len(y) - self.fft_size:
            window = np.asarray(y[idx : idx + self.fft_size], dtype=np.float32)
            oscofo.process_block(window)
            desc = oscofo.get_description()
            if desc.db >= self.min_db:
                count += 1
            idx += self.hop_size
        return count

    def _augmentation_variant(self, y, sr):
        """Return one deterministic-seed-controlled augmented variant."""
        choices = ["time_stretch", "pitch_shift", "noise"]
        if self.ir_files:
            choices.append("reverb")

        choice = random.choice(choices)
        if choice == "time_stretch":
            rate = random.uniform(0.7, 1.2)
            return self._time_stretch(y, rate)

        if choice == "pitch_shift":
            steps = random.uniform(-2, 2)
            return self._pitch_shift(y, steps)

        if choice == "noise":
            noise_val = random.uniform(0.005, 0.01)
            noise = np.random.normal(0, noise_val, len(y))
            return y + noise

        return self.apply_reverb(y, sr)

    def _time_stretch(self, y, rate):
        """Stretch/compress audio with SciPy to avoid librosa.effects JIT paths."""
        if len(y) == 0:
            return y
        target_length = max(1, int(round(len(y) / rate)))
        return resample(y, target_length).astype(np.float32)

    def _pitch_shift(self, y, steps):
        """Approximate pitch-shift augmentation without importing librosa.effects."""
        if len(y) == 0:
            return y
        factor = 2.0 ** (steps / 12.0)
        shifted_length = max(1, int(round(len(y) / factor)))
        shifted = resample(y, shifted_length)
        return resample(shifted, len(y)).astype(np.float32)

    def _augmentation_target(self, original_counts_by_class):
        counts = [count for count in original_counts_by_class.values() if count > 0]
        if not counts:
            return 0
        largest_count = max(counts)
        return int(math.ceil(largest_count * self.balanced_augmentation_factor))

    def _label_counts(self, data):
        counts = {}
        for _, label in data:
            counts[label] = counts.get(label, 0) + 1
        return counts

    def _split_original_train_validation(self, original_train_data):
        """Split only original train samples into fit and validation sets."""
        if not original_train_data:
            return [], []

        rng = np.random.default_rng(self.random_state)
        indices_by_label = {}
        for idx, (_, label) in enumerate(original_train_data):
            indices_by_label.setdefault(label, []).append(idx)

        validation_indices = []
        for label in sorted(indices_by_label):
            indices = np.asarray(indices_by_label[label], dtype=int)
            if len(indices) < 2:
                self.print(
                    f"Class {label}: skipping validation split because it has "
                    "fewer than two original train samples."
                )
                continue

            indices = rng.permutation(indices)
            n_val = int(round(len(indices) * self.test_fraction))
            n_val = max(1, min(n_val, len(indices) - 1))
            validation_indices.extend(indices[:n_val])

        validation_set = set(int(i) for i in validation_indices)
        validation_data = [
            sample
            for idx, sample in enumerate(original_train_data)
            if idx in validation_set
        ]

        if len(set(label for _, label in validation_data)) < 2:
            if validation_data:
                self.print(
                    "Validation disabled because original validation samples "
                    "contain fewer than two classes."
                )
            return original_train_data, []

        original_fit_data = [
            sample
            for idx, sample in enumerate(original_train_data)
            if idx not in validation_set
        ]
        return original_fit_data, validation_data

    def _assert_dataset_partitions(
        self,
        original_fit_data,
        augmented_train_data,
        validation_data,
        test_data,
    ):
        fit_data = original_fit_data + augmented_train_data
        if len(fit_data) != len([label for _, label in fit_data]):
            raise AssertionError("Fit feature/label length mismatch")
        if len(validation_data) != len([label for _, label in validation_data]):
            raise AssertionError("Validation feature/label length mismatch")
        if len(test_data) != len([label for _, label in test_data]):
            raise AssertionError("Test feature/label length mismatch")

        augmented_ids = {id(sample) for sample in augmented_train_data}
        validation_ids = {id(sample) for sample in validation_data}
        test_ids = {id(sample) for sample in test_data}
        fit_ids = {id(sample) for sample in fit_data}

        if augmented_ids & validation_ids:
            raise AssertionError("Augmented sample entered validation")
        if test_ids & fit_ids:
            raise AssertionError("Test sample entered fit data")
        if test_ids & validation_ids:
            raise AssertionError("Test sample entered validation")

        validation_classes = set(label for _, label in validation_data)
        if validation_data and len(validation_classes) < 2:
            raise AssertionError(
                "Validation must contain at least two classes or be disabled"
            )

    def _descriptor_values(self, desc):
        descriptors_values = []

        for d in self.descriptors:
            value = getattr(desc, d)
            if isinstance(value, (int, float, np.number)):
                descriptors_values.append(float(value))
            elif isinstance(value, np.ndarray):
                descriptors_values.extend(np.asarray(value, dtype=np.float64).ravel())
            elif isinstance(value, Iterable) and not isinstance(
                value, (str, bytes, dict)
            ):
                try:
                    descriptors_values.extend(
                        np.asarray(list(value), dtype=np.float64).ravel()
                    )
                except (TypeError, ValueError) as exc:
                    raise RuntimeError(
                        f"Descriptor '{d}' is not a numeric 1D container"
                    ) from exc
            else:
                raise RuntimeError(f"Descriptor '{d}' has unsupported type")

        return np.asarray(descriptors_values, dtype=np.float64)

    def _append_sample_from_description(self, desc, label, target_list):
        target_list.append((self._descriptor_values(desc), label))
        if self.on_progress:
            self.on_progress(self.processed_count)
        return True

    def _first_non_silent_description(self, y, sample_rate):
        oscofo = self._oscofo_for_sample_rate(sample_rate)
        idx = 0
        while idx <= len(y) - self.fft_size:
            window = np.asarray(y[idx : idx + self.fft_size], dtype=np.float32)
            oscofo.process_block(window)
            desc = oscofo.get_description()
            if desc.db >= self.min_db:
                return desc
            idx += self.hop_size
        return None

    def random_window(self, signal, size=2048):
        """Return a deterministic-seed-controlled random window from a signal."""
        if len(signal) < size:
            raise ValueError("Signal shorter than window size")

        start = random.randint(0, len(signal) - size)
        return signal[start : start + size]

    def _process_file(
        self,
        filepath,
        label,
        target_list,
        mode,
        window_split_role=None,
        augment=False,
        max_samples=None,
        sample_rate=None,
    ):
        sr = self.sample_rate if sample_rate is None else int(sample_rate)
        oscofo = self._oscofo_for_sample_rate(sr)
        y = self._load_audio(filepath, sr)
        y = self._single_file_temporal_region(y, window_split_role, filepath, label, sr)
        signal = self._augmentation_variant(y, sr) if augment else y
        appended_count = 0

        if len(signal) < self.fft_size:
            return appended_count

        idx = 0
        while idx <= len(signal) - self.fft_size:
            self.processed_count += 1
            window = np.asarray(signal[idx : idx + self.fft_size], dtype=np.float32)
            oscofo.process_block(window)
            desc = oscofo.get_description()

            if desc.db < self.min_db:
                idx += self.hop_size
                continue

            if self._append_sample_from_description(desc, label, target_list):
                appended_count += 1
                if max_samples is not None and appended_count >= max_samples:
                    return appended_count
            idx += self.hop_size

        if (
            mode == "traindata"
            and not augment
            and appended_count == 0
            and len(y) >= self.fft_size
        ):
            fallback_desc = self._first_non_silent_description(y, sr)
            if fallback_desc is not None and self._append_sample_from_description(
                fallback_desc, label, target_list
            ):
                appended_count += 1
                self.print(
                    f"Class {label}: fallback train sample added from {os.path.basename(filepath)}"
                )
            else:
                self.print(
                    f"Class {label}: no non-silent fallback window found in "
                    f"{os.path.basename(filepath)}"
                )

        return appended_count

    def get_train_mir(self):
        """Extract descriptor datasets from the configured training folder."""
        if not self.trainfolder or not self.folders:
            raise RuntimeError("You need to define the train folder first")

        self._seed_random_generators()
        self._single_file_split_cache = {}
        original_train_data, augmented_train_data, test_data = [], [], []

        class_files = {
            label: [
                os.path.join(folder, f)
                for f in os.listdir(folder)
                if f.endswith((".aif", ".aiff", ".wav"))
            ]
            for label, folder in self.folders.items()
        }

        self.min_y = None
        self.max_y = None
        file_raw_frames = {}

        for label, all_files in class_files.items():
            self.print(f"Class {label} raw files: {len(all_files)}")
            class_raw_frames = 0
            for file in all_files:
                file_raw_frames[file] = {}
                for sample_rate in self.training_sample_rates:
                    y = self._load_audio(file, sample_rate)
                    length = len(y)
                    if self.min_y is None or length < self.min_y:
                        self.min_y = length
                    if self.max_y is None or length > self.max_y:
                        self.max_y = length

                    valid_frames = self._count_valid_windows(y, sample_rate)
                    file_raw_frames[file][sample_rate] = valid_frames
                    class_raw_frames += valid_frames

            self.print(f"Class {label} raw non-silent frames: {class_raw_frames}")

        self.print(f"Minimum length: {self.min_y}")
        self.print(f"Maximum length: {self.max_y}")

        train_files_by_class = {}
        test_files_by_class = {}
        raw_train_frames_by_class = {}
        raw_test_frames_by_class = {}
        shared_single_file_split_by_class = {}

        for label, all_files in class_files.items():
            train_files, test_files, shared_single_file_split = self._split_train_test(
                all_files
            )
            train_files_by_class[label] = train_files
            test_files_by_class[label] = test_files
            shared_single_file_split_by_class[label] = shared_single_file_split
            if shared_single_file_split and train_files:
                train_frames = 0
                test_frames = 0
                for sample_rate in self.training_sample_rates:
                    y = self._load_audio(train_files[0], sample_rate)
                    train_region, test_region = self._single_file_split_regions(
                        y, train_files[0], label, sample_rate
                    )
                    train_frames += self._count_valid_windows(train_region, sample_rate)
                    test_frames += self._count_valid_windows(test_region, sample_rate)
                raw_train_frames_by_class[label] = train_frames
                raw_test_frames_by_class[label] = test_frames
            else:
                raw_train_frames_by_class[label] = sum(
                    sum(file_raw_frames.get(f, {}).values()) for f in train_files
                )
                raw_test_frames_by_class[label] = sum(
                    sum(file_raw_frames.get(f, {}).values()) for f in test_files
                )

        for label in self.folders.keys():
            raw_train_files = len(train_files_by_class.get(label, []))
            raw_test_files = len(test_files_by_class.get(label, []))
            raw_train_frames = raw_train_frames_by_class.get(label, 0)
            raw_test_frames = raw_test_frames_by_class.get(label, 0)
            self.print(
                f"Class {label}: "
                f"raw train files={raw_train_files}, raw test files={raw_test_files}, "
                f"raw train non-silent frames={raw_train_frames}, "
                f"raw test non-silent frames={raw_test_frames}, "
                f"single-file-window-split={shared_single_file_split_by_class.get(label, False)}"
            )

        generated_train_by_class = {label: 0 for label in self.folders.keys()}
        augmented_train_by_class = {label: 0 for label in self.folders.keys()}
        generated_test_by_class = {label: 0 for label in self.folders.keys()}

        for label, _ in self.folders.items():
            self.print(f"Processing {label} class originals")
            train_files = train_files_by_class.get(label, [])
            test_files = test_files_by_class.get(label, [])
            shared_split = shared_single_file_split_by_class.get(label, False)

            for f in test_files:
                for sample_rate in self.training_sample_rates:
                    generated_test_by_class[label] += self._process_file(
                        f,
                        label,
                        test_data,
                        "testdata",
                        window_split_role=("test" if shared_split else None),
                        sample_rate=sample_rate,
                    )

            for f in train_files:
                for sample_rate in self.training_sample_rates:
                    generated_train_by_class[label] += self._process_file(
                        f,
                        label,
                        original_train_data,
                        "traindata",
                        window_split_role=("train" if shared_split else None),
                        sample_rate=sample_rate,
                    )

        train_samples_target = self._augmentation_target(generated_train_by_class)
        if train_samples_target > 0:
            self.print(f"Target train samples/class: {train_samples_target}")

        for label, _ in self.folders.items():
            train_files = train_files_by_class.get(label, [])
            shared_split = shared_single_file_split_by_class.get(label, False)
            original_count = generated_train_by_class[label]
            needed = max(0, train_samples_target - original_count)

            if needed == 0:
                continue

            if not train_files:
                self.print(
                    f"Class {label}: warning, cannot augment because there are "
                    "no training files."
                )
                continue

            max_attempts = int(self.max_augmentation_loops) * len(train_files)
            if max_attempts <= 0:
                self.print(
                    f"Class {label}: warning, cannot augment because "
                    "max_augmentation_loops is not positive."
                )
                continue

            self.print(f"Augmenting {label} class: target additional samples={needed}")
            attempts = 0
            file_index = 0
            while (
                original_count + augmented_train_by_class[label] < train_samples_target
                and attempts < max_attempts
            ):
                f = train_files[file_index % len(train_files)]
                sample_rate = self.training_sample_rates[
                    attempts % len(self.training_sample_rates)
                ]
                remaining = (
                    train_samples_target
                    - original_count
                    - augmented_train_by_class[label]
                )
                added = self._process_file(
                    f,
                    label,
                    augmented_train_data,
                    "traindata",
                    window_split_role=("train" if shared_split else None),
                    augment=True,
                    max_samples=remaining,
                    sample_rate=sample_rate,
                )
                augmented_train_by_class[label] += added
                attempts += 1
                file_index += 1

            final_count = original_count + augmented_train_by_class[label]
            if final_count < train_samples_target:
                self.print(
                    f"Class {label}: warning, augmentation stopped before target. "
                    f"Reached {final_count}/{train_samples_target} samples after "
                    f"{attempts}/{max_attempts} augmentation attempts. "
                    "Possible reasons: source files are too short, augmented "
                    "variants are below min_db, or max_augmentation_loops is too low."
                )

        original_fit_data, validation_data = self._split_original_train_validation(
            original_train_data
        )
        fit_data = original_fit_data + augmented_train_data
        self._assert_dataset_partitions(
            original_fit_data,
            augmented_train_data,
            validation_data,
            test_data,
        )

        original_fit_counts = self._label_counts(original_fit_data)
        validation_counts = self._label_counts(validation_data)

        for label in self.folders.keys():
            train_count = generated_train_by_class.get(label, 0)
            augmented_count = augmented_train_by_class.get(label, 0)
            original_fit_count = original_fit_counts.get(label, 0)
            validation_count = validation_counts.get(label, 0)
            test_count = generated_test_by_class.get(label, 0)
            self.print(
                f"Class {label}:\n"
                f"    original train samples total: {train_count}\n"
                f"    original fit samples: {original_fit_count}\n"
                f"    original validation samples: {validation_count}\n"
                f"    target train samples: {train_samples_target}\n"
                f"    needed augmented samples: "
                f"{max(0, train_samples_target - train_count)}\n"
                f"    augmented fit samples: {augmented_count}\n"
                f"    final generated train samples: {train_count + augmented_count}\n"
                f"    final CatBoost fit samples: {original_fit_count + augmented_count}\n"
                f"    independent test samples: {test_count}"
            )
            if train_count < train_samples_target and augmented_count == 0:
                raise RuntimeError(
                    f"Augmentation failed for class {label}: "
                    f"original_count={train_count}, "
                    f"target_per_class={train_samples_target}, "
                    "augmented_count=0"
                )
            if train_count == 0:
                self.print(
                    f"Class {label}: warning, no training samples survived filtering."
                )
            if test_count == 0:
                self.print(
                    f"Class {label}: warning, no independent test samples "
                    f"survived filtering."
                )

        self.dataset_counts = {
            "original_train": generated_train_by_class,
            "original_fit": original_fit_counts,
            "validation": validation_counts,
            "augmented_fit": augmented_train_by_class,
            "test": generated_test_by_class,
            "target_per_class": train_samples_target,
        }

        self.x_train, self.y_train = zip(*fit_data) if fit_data else ([], [])
        self.x_val, self.y_val = zip(*validation_data) if validation_data else ([], [])
        self.x_test, self.y_test = zip(*test_data) if test_data else ([], [])

        self.x_np_train = np.array(self.x_train)
        self.y_np_train = np.array(self.y_train)
        self.x_np_val = np.array(self.x_val)
        self.y_np_val = np.array(self.y_val)
        self.x_np_test = np.array(self.x_test)
        self.y_np_test = np.array(self.y_test)
        assert len(self.x_np_train) == len(self.y_np_train)
        assert len(self.x_np_val) == len(self.y_np_val)
        assert len(self.x_np_test) == len(self.y_np_test)
        assert len(self.y_np_val) == 0 or len(np.unique(self.y_np_val)) >= 2
        self.print("Done!")

        self.save_data(os.path.join(self.base_path, "dataset.npz"))

    def save_data(self, path="dataset.npz"):
        """Save the extracted dataset and cache metadata."""
        np.savez_compressed(
            path,
            x_train=self.x_np_train,
            y_train=self.y_np_train,
            x_val=self.x_np_val,
            y_val=self.y_np_val,
            x_test=self.x_np_test,
            y_test=self.y_np_test,
            min_y=self.min_y,
            max_y=self.max_y,
            dataset_counts=np.asarray(self.dataset_counts, dtype=object),
            metadata=np.asarray(self._dataset_metadata(), dtype=object),
        )

    def load_data(self, path="dataset.npz"):
        """Load a previously extracted dataset."""
        data = np.load(path, allow_pickle=True)
        self.x_np_train = data["x_train"]
        self.y_np_train = data["y_train"]
        self.x_np_val = data["x_val"]
        self.y_np_val = data["y_val"]
        self.x_np_test = data["x_test"]
        self.y_np_test = data["y_test"]
        self.min_y = data["min_y"]
        self.max_y = data["max_y"]
        self.x_train = list(self.x_np_train)
        self.y_train = list(self.y_np_train)
        self.x_val = list(self.x_np_val)
        self.y_val = list(self.y_np_val)
        self.x_test = list(self.x_np_test)
        self.y_test = list(self.y_np_test)
        self.dataset_counts = (
            data["dataset_counts"].item() if "dataset_counts" in data.files else {}
        )

    def _audio_file_signature(self):
        files = []
        for label, folder in sorted(self.folders.items()):
            for name in sorted(os.listdir(folder)):
                if not name.endswith((".aif", ".aiff", ".wav")):
                    continue
                path = os.path.join(folder, name)
                stat = os.stat(path)
                files.append(
                    {
                        "label": label,
                        "path": os.path.relpath(path, self.trainfolder),
                        "size": stat.st_size,
                        "mtime_ns": stat.st_mtime_ns,
                    }
                )
        return files

    def _dataset_metadata(self):
        """Return configuration that affects descriptor dataset generation."""
        return {
            "cache_version": self.DATASET_CACHE_VERSION,
            "descriptors": list(self.descriptors),
            "sample_rate": int(self.sample_rate),
            "training_sample_rates": [int(sr) for sr in self.training_sample_rates],
            "fft_size": int(self.fft_size),
            "hop_size": int(self.hop_size),
            "trainfolder": (
                os.path.abspath(self.trainfolder) if self.trainfolder else ""
            ),
            "audio_files": self._audio_file_signature() if self.folders else [],
            "min_db": float(self.min_db),
            "test_fraction": float(self.test_fraction),
            "single_file_class_window_split": bool(self.single_file_class_window_split),
            "max_augmentation_loops": int(self.max_augmentation_loops),
            "balanced_augmentation_factor": float(self.balanced_augmentation_factor),
            "random_state": int(self.random_state),
            "ir_files": [os.path.abspath(f) for f in self.ir_files or []],
        }

    def _dataset_matches_metadata(self, path="dataset.npz"):
        """Return True when a cached dataset was built with current settings."""
        try:
            data = np.load(path, allow_pickle=True)
            if "metadata" not in data.files:
                return False
            cached_metadata = data["metadata"].item()
        except (OSError, ValueError, KeyError, AttributeError):
            return False

        return cached_metadata == self._dataset_metadata()

    def _print_evaluation(self):
        y_true = np.asarray(self.y_np_test).ravel()
        y_pred = np.asarray(self.clf.predict(self.x_np_test)).ravel()

        total = len(y_true)
        correct = int(np.sum(y_true == y_pred))
        accuracy = correct / total if total else 0.0
        self.print(f"Test accuracy: {accuracy:.3f} ({correct}/{total})")

        per_class_accuracy = []
        labels = sorted(set(np.asarray(self.y_np_train).ravel()) | set(y_true))
        for label in labels:
            mask = y_true == label
            label_total = int(np.sum(mask))
            if label_total == 0:
                self.print(f"Class {label}: no independent test samples")
                continue
            label_correct = int(np.sum(y_pred[mask] == label))
            label_accuracy = label_correct / label_total
            per_class_accuracy.append(label_accuracy)
            self.print(
                f"Class {label}: {label_correct}/{label_total} correct "
                f"({label_accuracy:.3f})"
            )

        if per_class_accuracy:
            balanced_accuracy = float(np.mean(per_class_accuracy))
            self.print(f"Balanced test accuracy: {balanced_accuracy:.3f}")

        prediction_counts = {
            label: int(np.sum(y_pred == label)) for label in sorted(set(y_pred))
        }
        self.print(f"Test prediction counts: {prediction_counts}")

    def _validation_eval_set(self):
        if len(self.x_np_val) == 0 or len(self.y_np_val) == 0:
            return None

        if len(self.x_np_val) != len(self.y_np_val):
            raise AssertionError("Validation feature/label length mismatch")

        if len(np.unique(self.y_np_val)) < 2:
            raise AssertionError(
                "Validation must contain at least two classes or be disabled"
            )

        return (self.x_np_val, self.y_np_val)

    def _train(self):
        if len(self.x_np_train) == 0 or len(self.y_np_train) == 0:
            raise RuntimeError("No training samples found.")

        unique_train = np.unique(self.y_np_train)
        if len(unique_train) < 2:
            raise RuntimeError("Training target has only one class.")

        # Early stopping uses only original train validation samples prepared
        # during dataset generation. Augmented and test samples never enter it.
        eval_set = self._validation_eval_set()

        fit_kwargs = {"use_best_model": eval_set is not None}
        self.clf.fit(
            self.x_np_train,
            self.y_np_train,
            eval_set=eval_set,
            **fit_kwargs,
        )

        if eval_set:
            self.print(
                f"Model shrunk to {self.clf.tree_count_} trees via early stopping."
            )

        if len(self.x_np_test) > 0 and len(self.y_np_test) > 0:
            self._print_evaluation()
        else:
            self.print(
                "No independent test samples available; skipped test evaluation."
            )

        self.print("Training finished!")

    def set_descriptors(self, args: List[str]):
        """Select OpenScofo descriptor names used as classifier features."""
        arr = np.zeros(self.fft_size)
        self.oscofo.process_block(arr)
        random_des = self.oscofo.get_description()
        for a in args:
            if not hasattr(random_des, a):
                raise RuntimeError(f"The descriptor '{a}' does not exist on OpenScofo")
        self.descriptors = args

    def print_data(self):
        """Print dataset sample and label counts."""
        self.print(f"Train samples: {len(self.x_train)}")
        self.print(f"Train labels: {len(self.y_train)}")
        assert len(self.x_train) == len(self.y_train)

        self.print(f"Validation samples: {len(self.x_val)}")
        self.print(f"Validation labels: {len(self.y_val)}")
        assert len(self.x_val) == len(self.y_val)

        self.print(f"Test samples: {len(self.x_test)}")
        self.print(f"Test labels: {len(self.y_test)}")
        assert len(self.x_test) == len(self.y_test)

    def analyze(self):
        """Build or load the descriptor dataset for the configured training set."""
        if len(self.descriptors) == 0:
            raise RuntimeError("You need to define the descriptors first")

        path = os.path.join(self.base_path, "dataset.npz")
        if os.path.exists(path) and self._dataset_matches_metadata(path):
            self.load_data(path)
        else:
            if os.path.exists(path):
                self.print("Cached dataset metadata mismatch; rebuilding dataset.")
            self.get_train_mir()

    def clear_data(self):
        """Delete the cached dataset, if it exists."""
        path = os.path.join(self.base_path, "dataset.npz")
        if os.path.exists(path):
            os.remove(path)
            self.print(f"Deleted cached dataset: {path}")
        else:
            self.print(f"No cached dataset found at: {path}")

    def set_ir_folders(self, args: List[str]):
        """Set absolute or base_path-relative impulse response folders."""
        self.ir_folders = []
        self.ir_files = []
        for folder in args:
            resolved = self.resolve_trainfolder(folder)
            if not os.path.isdir(resolved):
                raise NotADirectoryError(f"IR folder is not a directory: {resolved}")
            self.ir_folders.append(resolved)
            self.ir_files.extend(
                os.path.join(resolved, f)
                for f in sorted(os.listdir(resolved))
                if f.endswith((".wav", ".aif", ".aiff"))
            )

        if self.ir_files:
            self.print(f"Found {len(self.ir_files)} IR files")

    def set_train_folder(self, folder: str):
        """Set the dataset root; each subfolder is treated as one class label."""
        self.trainfolder = self.resolve_trainfolder(folder)
        self.folders = {
            f: os.path.join(self.trainfolder, f)
            for f in sorted(os.listdir(self.trainfolder))
            if os.path.isdir(os.path.join(self.trainfolder, f))
        }
        self.print("Train folder: " + self.trainfolder)

    def _write_onnx_metadata(self, model_path: str):
        """Add OpenScofo training metadata to an exported ONNX model."""
        try:
            import onnx
        except ImportError as exc:
            raise RuntimeError(
                "The Python 'onnx' package is required to write ONNX metadata."
            ) from exc

        model = onnx.load(model_path)
        metadata = {
            "openscofo.descriptors": json.dumps(list(self.descriptors)),
            "openscofo.labels": json.dumps(
                [str(label) for label in getattr(self.clf, "classes_", [])]
            ),
            "openscofo.sample_rate": str(int(self.sample_rate)),
            "openscofo.fft_size": str(int(self.fft_size)),
            "openscofo.hop_size": str(int(self.hop_size)),
            "openscofo.version": self._openscofo_version(),
        }

        existing = {prop.key: prop for prop in model.metadata_props}
        for key, value in metadata.items():
            prop = existing.get(key)
            if prop is None:
                prop = model.metadata_props.add()
                prop.key = key
            prop.value = value

        onnx.save(model, model_path)

    def _openscofo_version(self):
        try:
            from . import _OpenScofo

            version = getattr(_OpenScofo, "__version__", None)
            if version:
                return str(version)
        except ImportError:
            pass

        try:
            return importlib_metadata.version("OpenScofo")
        except importlib_metadata.PackageNotFoundError:
            return "unknown"

    def export_model(self, model_path: str):
        """Export the trained CatBoost model to CatBoost or ONNX format."""
        if self.clf is None:
            raise RuntimeError("No trained model to export. Call train() first.")
        is_fitted = getattr(self.clf, "is_fitted", None)
        if callable(is_fitted):
            is_fitted = is_fitted()
        if is_fitted is False:
            raise RuntimeError("No trained model to export. Call train() first.")

        path = os.path.join(self.base_path, model_path)
        if os.path.exists(path):
            self.print("Model already exists, will replace it!")

        if path.endswith(".onnx"):
            self.clf.save_model(path, format="onnx")
            self._write_onnx_metadata(path)
        else:
            self.clf.save_model(path)
        self.print(f"Model exported to {path}")

    def set_print_callback(self, func):
        """Set the function used for status messages."""
        self.print = func

    def set_onprogress_callback(self, func):
        """Set a callback invoked after each appended training sample."""
        self.on_progress = func

    def set_model_type(self, model_type: str):
        """Select the classifier backend."""
        normalized = model_type.lower().replace("-", "").replace("_", "")
        aliases = {
            "cat": "catboost",
            "catboost": "catboost",
        }

        if normalized not in aliases:
            supported = ", ".join(self.SUPPORTED_MODEL_TYPES)
            raise ValueError(
                f"Unsupported model_type '{model_type}'. Use one of: {supported}"
            )

        self.model_type = aliases[normalized]

    def set_catboost_config(
        self,
        iterations=1000,
        learning_rate=0.05,
        depth=6,
        early_stopping_rounds=70,
        thread_count=-1,
    ):
        """Configure CatBoost training parameters."""
        self.iterations = iterations
        self.learning_rate = learning_rate
        self.depth = depth
        self.early_stopping_rounds = early_stopping_rounds
        self.thread_count = thread_count

    def _build_classifier(self):
        return CatBoostClassifier(
            iterations=self.iterations,
            depth=self.depth,
            learning_rate=self.learning_rate,
            loss_function="MultiClass",
            thread_count=self.thread_count,
            early_stopping_rounds=self.early_stopping_rounds,
            random_seed=self.random_state,
        )

    def train(self):
        """Train the classifier from the extracted or cached dataset."""
        self.clf = self._build_classifier()

        self.print(f"Training {self.model_type}, wait...")

        if len(self.descriptors) == 0:
            raise RuntimeError("You need to define the descriptors used to train")

        if len(self.x_np_test) == 0 or len(self.x_np_train) == 0:
            self.analyze()

        self._train()
