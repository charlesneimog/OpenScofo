try:
    from ._OpenScofo import (
        OpenScofo,
        EventType,
        HMMType,
        AudioState,
        MarkovState,
        Description,
        Configuration,
    )
except:
    from _OpenScofo import (
        OpenScofo,
        EventType,
        HMMType,
        AudioState,
        MarkovState,
        Description,
        Configuration,
    )

__all__ = [
    "OpenScofo",
    "EventType",
    "HMMType",
    "AudioState",
    "MarkovState",
    "Description",
    "Configuration",
]


import os
import random
import json
import math
from typing import Optional, Callable, List
import numpy as np
import librosa
from scipy.signal import fftconvolve

from sklearn.metrics import classification_report
from catboost import CatBoostClassifier
import joblib

from skl2onnx import convert_sklearn
from skl2onnx.common.data_types import FloatTensorType

import onnx


class ExtendedTechniqueClassifier:

    def __init__(self, sample_rate=48000, fft_size=2048, hop_size=512, base_path=None):
        self.print = print

        # configuration
        self.sample_rate = sample_rate
        self.fft_size = fft_size
        self.hop_size = hop_size
        self.base_path = base_path or os.getcwd()

        # internal state
        self.tabname = "train"
        self.redraw_tab_after = 120
        self.processed_count = 0
        self.descriptors = []
        self.redraw = 0

        # train parameters
        self.iterations = 1000
        self.random_state = 42
        self.test_fraction = 0.2
        self.min_db = -60.0
        self.oscofo: OpenScofo = OpenScofo(sample_rate, fft_size, hop_size)
        self.base_augmentation_loops = 1
        self.max_augmentation_loops = 64
        self.single_file_class_window_split = True

        # folders
        self.trainfolder = ""
        self.folders = {}
        self.currtraindata = []
        self.currtestdata = []

        # datasets
        self.x_np_train, self.y_np_train = [], []
        self.x_np_test, self.y_np_test = [], []
        self.x_train, self.y_train = [], []
        self.x_test, self.y_test = [], []

        # model
        self.model_type = "catboost"

        # ir
        self.ir_files = None
        self.on_progress = None

    # Folder / Dataset Management
    def resolve_trainfolder(self, path):
        if os.path.exists(path):
            return path
        candidate = os.path.join(self.base_path, path)
        if not os.path.exists(candidate):
            raise FileNotFoundError(f"{path} folder not found")
        return candidate

    def apply_reverb(self, y, sr):
        if not self.ir_files:
            return y

        ir_file = random.choice(self.ir_files)
        ir, _ = librosa.load(ir_file, sr=sr)
        ir = ir / (np.max(np.abs(ir)) + 1e-8)
        y_conv = fftconvolve(y, ir, mode="full")
        y_conv = y_conv[: len(y)]
        y_conv = y_conv / (np.max(np.abs(y_conv)) + 1e-8)
        return y_conv.astype(np.float32)

    # Audio processing
    def _load_audio(self, filepath):
        y, _ = librosa.load(filepath, sr=self.sample_rate)
        return y

    def _split_train_test(self, files):
        if not files:
            return [], [], False

        if len(files) == 1:
            if self.single_file_class_window_split:
                return files, files, True
            return files, [], False

        n_test = int(len(files) * self.test_fraction)
        n_test = max(1, n_test)
        n_test = min(n_test, len(files) - 1)

        test_files = random.sample(files, n_test)
        train_files = [f for f in files if f not in test_files]
        return train_files, test_files, False

    def _count_valid_windows(self, y):
        count = 0
        idx = 0
        while idx <= len(y) - self.fft_size:
            window = np.asarray(y[idx : idx + self.fft_size], dtype=np.float32)
            self.oscofo.process_block(window)
            desc = self.oscofo.get_description()
            if desc.db >= self.min_db:
                count += 1
            idx += self.hop_size
        return count

    def _augmented_data(self, y, sr, mode="traindata", loop_iterations=None):
        variants = [y]
        if mode != "traindata":
            return variants

        if loop_iterations is None:
            target_frames = 500
            base_frames = max(1, self._count_valid_windows(y))
            total_variants_needed = target_frames / base_frames
            loop_iterations = int(total_variants_needed / 4)
            loop_iterations = max(1, min(loop_iterations, self.max_augmentation_loops))

        for _ in range(loop_iterations):
            rate = random.uniform(0.7, 1.2)
            variants.append(librosa.effects.time_stretch(y=y, rate=rate))

            steps = random.uniform(-2, 2)
            variants.append(librosa.effects.pitch_shift(y=y, sr=sr, n_steps=steps))

            noise_val = random.uniform(0.005, 0.01)
            noise = np.random.normal(0, noise_val, len(y))
            variants.append(y + noise)

            variants.append(self.apply_reverb(y, sr))

        return variants

    def _compute_class_augmentation_loops(self, raw_train_frames_by_class):
        raw_frames = {
            label: int(frames)
            for label, frames in raw_train_frames_by_class.items()
            if int(frames) > 0
        }

        if not raw_frames:
            return {}

        max_raw_frames = max(raw_frames.values())
        loops_by_class = {}

        for label, frame_count in raw_frames.items():
            ratio = max_raw_frames / max(1, frame_count)
            loops = int(math.ceil(self.base_augmentation_loops * ratio))
            loops = max(1, min(loops, self.max_augmentation_loops))
            loops_by_class[label] = loops

        return loops_by_class

    def _append_sample_from_window(self, window, label, target_list):
        self.oscofo.process_block(window)
        desc = self.oscofo.get_description()

        descriptors_values = []
        for d in self.descriptors:
            value = getattr(desc, d)
            if isinstance(value, (int, float)):
                descriptors_values.append(value)
            elif isinstance(value, list):
                descriptors_values.extend(value)
            else:
                raise RuntimeError("Wrong descriptor")

        desc_arr = np.asarray(descriptors_values, dtype=np.float64)
        target_list.append((desc_arr, label))
        if self.on_progress:
            self.on_progress(self.processed_count)
        return True

    def random_window(self, signal, size=2048):
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
        augmentation_loops=None,
        window_split_role=None,
    ):
        y = self._load_audio(filepath)
        sr = self.sample_rate
        variants = self._augmented_data(
            y,
            sr,
            mode=mode,
            loop_iterations=augmentation_loops,
        )
        appended_count = 0
        non_silent_seen = 0
        split_period = max(2, int(round(1.0 / max(1e-6, self.test_fraction))))

        for signal in variants:
            if len(signal) < self.fft_size:
                continue

            idx = random.randint(0, min(1024, len(signal) - self.fft_size))
            while idx <= len(signal) - self.fft_size:
                self.processed_count += 1
                window = np.asarray(signal[idx : idx + self.fft_size], dtype=np.float32)
                self.oscofo.process_block(window)
                desc = self.oscofo.get_description()

                if desc.db < self.min_db:
                    idx += self.hop_size
                    continue

                if window_split_role in ("train", "test"):
                    in_test_bucket = (non_silent_seen % split_period) == 0
                    non_silent_seen += 1

                    if window_split_role == "test" and not in_test_bucket:
                        idx += self.hop_size
                        continue

                    if window_split_role == "train" and in_test_bucket:
                        idx += self.hop_size
                        continue

                if self._append_sample_from_window(window, label, target_list):
                    appended_count += 1
                idx += self.hop_size

        if mode == "traindata" and appended_count == 0 and len(y) >= self.fft_size:
            fallback_window = np.asarray(y[: self.fft_size], dtype=np.float32)
            self.processed_count += 1
            if self._append_sample_from_window(fallback_window, label, target_list):
                appended_count += 1
                self.print(
                    f"Class {label}: fallback train sample added from {os.path.basename(filepath)}"
                )

        return appended_count

    # ----------------------------
    # Dataset Build
    # ----------------------------
    def get_train_mir(self):
        train_data, test_data = [], []

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
                y = self._load_audio(file)
                length = len(y)
                if self.min_y is None or length < self.min_y:
                    self.min_y = length
                if self.max_y is None or length > self.max_y:
                    self.max_y = length

                valid_frames = self._count_valid_windows(y)
                file_raw_frames[file] = valid_frames
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
            raw_train_frames_by_class[label] = sum(
                file_raw_frames.get(f, 0) for f in train_files
            )
            raw_test_frames_by_class[label] = sum(
                file_raw_frames.get(f, 0) for f in test_files
            )

        loops_by_class = self._compute_class_augmentation_loops(
            raw_train_frames_by_class
        )

        for label in self.folders.keys():
            raw_train_files = len(train_files_by_class.get(label, []))
            raw_test_files = len(test_files_by_class.get(label, []))
            raw_train_frames = raw_train_frames_by_class.get(label, 0)
            raw_test_frames = raw_test_frames_by_class.get(label, 0)
            loops = loops_by_class.get(label, 1)
            self.print(
                f"Class {label}: "
                f"raw train files={raw_train_files}, raw test files={raw_test_files}, "
                f"raw train non-silent frames={raw_train_frames}, "
                f"raw test non-silent frames={raw_test_frames}, "
                f"augmentation loops/file={loops}, "
                f"single-file-window-split={shared_single_file_split_by_class.get(label, False)}"
            )

        generated_train_by_class = {label: 0 for label in self.folders.keys()}
        generated_test_by_class = {label: 0 for label in self.folders.keys()}

        for label, _ in self.folders.items():
            self.print(f"Processing {label} class")
            train_files = train_files_by_class.get(label, [])
            test_files = test_files_by_class.get(label, [])
            augmentation_loops = loops_by_class.get(label, 1)
            shared_split = shared_single_file_split_by_class.get(label, False)

            for f in test_files:
                generated_test_by_class[label] += self._process_file(
                    f,
                    label,
                    test_data,
                    "testdata",
                    augmentation_loops=0,
                    window_split_role=("test" if shared_split else None),
                )

            for f in train_files:
                generated_train_by_class[label] += self._process_file(
                    f,
                    label,
                    train_data,
                    "traindata",
                    augmentation_loops=augmentation_loops,
                    window_split_role=("train" if shared_split else None),
                )

        for label in self.folders.keys():
            self.print(
                f"Class {label} generated samples: "
                f"train={generated_train_by_class.get(label, 0)}, "
                f"test={generated_test_by_class.get(label, 0)}"
            )

        self.x_train, self.y_train = zip(*train_data) if train_data else ([], [])
        self.x_test, self.y_test = zip(*test_data) if test_data else ([], [])

        self.x_np_train = np.array(self.x_train)
        self.y_np_train = np.array(self.y_train)
        self.x_np_test = np.array(self.x_test)
        self.y_np_test = np.array(self.y_test)
        self.print("Done!")

        p = self.base_path
        self.save_data(p + "/" + "dataset.npz")

    def save_data(self, path="dataset.npz"):
        np.savez_compressed(
            path,
            x_train=self.x_np_train,
            y_train=self.y_np_train,
            x_test=self.x_np_test,
            y_test=self.y_np_test,
            min_y=self.min_y,
            max_y=self.max_y,
        )

    def load_data(self, path="dataset.npz"):
        data = np.load(path, allow_pickle=True)
        self.x_np_train = data["x_train"]
        self.y_np_train = data["y_train"]
        self.x_np_test = data["x_test"]
        self.y_np_test = data["y_test"]
        self.min_y = data["min_y"]
        self.max_y = data["max_y"]

    # ----------------------------
    # Training / Export
    # ----------------------------
    def _train(self):
        if len(self.x_np_train) == 0 or len(self.y_np_train) == 0:
            raise RuntimeError("No training samples found.")

        unique_train = np.unique(self.y_np_train)
        if len(unique_train) < 2:
            raise RuntimeError("Training target has only one class.")

        eval_set = None
        if len(self.x_np_test) > 0 and len(self.y_np_test) > 0:
            eval_set = (self.x_np_test, self.y_np_test)

        self.clf.fit(
            self.x_np_train,
            self.y_np_train,
            eval_set=eval_set,
            use_best_model=True,
        )

        if eval_set:
            y_true = np.asarray(self.y_np_test).ravel()
            y_pred = np.asarray(self.clf.predict(self.x_np_test)).ravel()
            self.print(classification_report(y_true, y_pred))
            self.print(
                f"Model shrunk to {self.clf.tree_count_} trees via early stopping."
            )

        self.print("Training finished!")

    def _build_payload(self, serialized_model):
        return {
            "type": "openscofo.ml_model.v1",
            "model_type": self.model_type,
            "descriptors": list(self.descriptors),
            "sample_rate": self.sample_rate,
            "model": serialized_model,
        }

    def _export_onnx(self, path):
        if self.model_type == "catboost":
            self.clf.save_model(path, format="onnx")
            return

        if convert_sklearn is None or FloatTensorType is None:
            raise RuntimeError("skl2onnx is not available")

        n_features = int(self.x_np_train.shape[1])
        onnx_model = convert_sklearn(
            self.clf,
            initial_types=[("input", FloatTensorType([None, n_features]))],
            target_opset=17,
        )

        metadata = self._build_payload("onnx")
        metadata["descriptors"] = json.dumps(metadata["descriptors"])
        for key, value in metadata.items():
            prop = onnx_model.metadata_props.add()
            prop.key = str(key)
            prop.value = str(value)

        onnx.save(onnx_model, path)

    # ╭──────────────────────────────────────╮
    # │           Public Methods             │
    # ╰──────────────────────────────────────╯

    def set_descriptors(self, args: List[str]):
        arr = np.zeros(self.fft_size)
        self.oscofo.process_block(arr)
        random_des = self.oscofo.get_description()
        for a in args:
            _ = getattr(random_des, a)
            if not (hasattr(random_des, a)):
                raise RuntimeError(f"The descriptor '{a}' does not exists on OpenScofo")
        self.descriptors = args

    def print_data(self):
        self.print(f"Train samples: {len(self.x_train)}")
        self.print(f"Train labels: {len(self.y_train)}")
        assert len(self.x_train) == len(self.y_train)

        self.print(f"Test samples: {len(self.x_test)}")
        self.print(f"Test labels: {len(self.y_test)}")
        assert len(self.x_test) == len(self.y_test)

    def analyze(self):
        if len(self.descriptors) == 0:
            raise RuntimeError("You need to defined the descriptors first")

        p = self.base_path + "/" + "dataset.npz"
        if os.path.exists(p):
            self.load_data(self.base_path + "/" + "dataset.npz")
        else:
            self.get_train_mir()

    def clear_data(self):
        p = self.base_path + "/" + "dataset.npz"
        if os.path.exists(p):
            os.remove(p)

    def set_ir_folders(self, args: List[str]):
        self.ir_folders = args
        for folder in args:
            self.ir_files = [
                os.path.join(folder, f)
                for f in os.listdir(folder)
                if f.endswith((".wav", ".aif", ".aiff"))
            ]

        if self.ir_files:
            self.print(f"Found {len(self.ir_files)} IR files")

    def set_train_folder(self, folder: str):
        self.trainfolder = self.resolve_trainfolder(folder)
        self.folders = {
            f: os.path.join(self.trainfolder, f)
            for f in os.listdir(self.trainfolder)
            if os.path.isdir(os.path.join(self.trainfolder, f))
        }
        self.print("Train folder: " + self.trainfolder)

    def export_model(self, model_path: str):
        file = model_path
        path = os.path.join(self.base_path, file)
        if os.path.exists(path):
            self.print("Model already exists, will replace it!")

        if file.endswith(".onnx"):
            try:
                self._export_onnx(path)
                self.print(f"Model exported to {path}")
                return
            except Exception as exc:
                fallback_path = path[:-5] + ".joblib"
                payload = self._build_payload(self.clf)
                joblib.dump(payload, fallback_path)
                self.print(
                    f"ONNX export failed ({exc}). Fallback to joblib at {fallback_path}"
                )
                return

        payload = self._build_payload(self.clf)
        joblib.dump(payload, path)
        self.print(f"Model exported to {path}")

    def set_print_callback(self, func):
        self.print = func

    def set_onprogress_callback(self, func):
        self.on_progress = func

    def set_catboost_config(
        self,
        iterations=1000,
        learning_rate=0.05,
        depth=6,
        early_stopping_rounds=70,
        thread_count=-1,
    ):
        self.iterations = iterations
        self.learning_rate = learning_rate
        self.depth = depth
        self.early_stopping_rounds = early_stopping_rounds
        self.thread_count = thread_count

    def train(self):
        self.set_catboost_config()

        self.clf = CatBoostClassifier(
            iterations=self.iterations,
            depth=self.depth,
            learning_rate=self.learning_rate,
            loss_function="MultiClass",
            thread_count=self.thread_count,
            early_stopping_rounds=self.early_stopping_rounds,
        )

        self.print("Training, wait...")

        if len(self.descriptors) == 0:
            raise RuntimeError("You need to define the descriptors used to train")

        if len(self.x_np_test) == 0 or len(self.x_np_train) == 0:
            self.load_data(self.base_path + "/" + "dataset.npz")

        self._train()
