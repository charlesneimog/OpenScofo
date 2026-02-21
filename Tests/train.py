#!/usr/bin/env python3

import os
import re
import random
import time
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import librosa

from sklearn.metrics import classification_report
from catboost import CatBoostClassifier

import OScofo


# ----------------------------
# Config
# ----------------------------
ITERATIONS = 150
RANDOM_STATE = 42
TEST_FRACTION = 0.2

SR = 48000
FFT = 4096
HOP = 1024

ENERGY_THRESHOLD = 1e-4

CACHE_NAME = f"dataset_cache_fft{FFT}_hop{HOP}.npz"

PITCH_RE = re.compile(r"-([A-G]#?\d)-")


# ----------------------------
# Progress helper
# ----------------------------
class Progress:
    def __init__(self):
        self.start = time.time()
        self.files = 0
        self.frames = 0
        self.per_pitch = defaultdict(int)

    def print_status(self):
        elapsed = time.time() - self.start
        print(
            f"[dataset] files={self.files:4d} "
            f"frames={self.frames:7d} "
            f"time={elapsed:6.1f}s"
        )


# ----------------------------
# Audio utils
# ----------------------------
def load_audio(filepath):
    y, _ = librosa.load(filepath, sr=SR, mono=True)
    return y


def generate_variants(y, mode):
    if mode != "traindata":
        return [y]

    variants = [y]

    for _ in range(2):
        rate = random.uniform(0.7, 1.2)
        variants.append(librosa.effects.time_stretch(y=y, rate=rate))

    for _ in range(2):
        steps = random.uniform(-2, 2)
        variants.append(librosa.effects.pitch_shift(y=y, sr=SR, n_steps=steps))

    noise_amp = random.uniform(0.005, 0.009)
    noise = np.random.normal(0, noise_amp, len(y))
    variants.append(y + noise)

    return variants


def frame_signal(y):
    for i in range(0, len(y) - FFT + 1, HOP):
        yield y[i : i + FFT]


# ----------------------------
# Pitch + features
# ----------------------------
def pitch_from_filename(filename):
    m = PITCH_RE.search(filename)
    return m.group(1) if m else None


def extract_features(frame):
    scofo = OScofo.OScofo(SR, FFT, HOP)
    desc = scofo.get_audio_description(frame)
    return np.asarray(desc.pseudo_cqt, dtype=np.float32)


# ----------------------------
# Worker (1 file)
# ----------------------------
def process_file(args):
    filepath, mode = args

    pitch = pitch_from_filename(os.path.basename(filepath))
    if pitch is None:
        return None

    y = load_audio(filepath)
    variants = generate_variants(y, mode)

    x_local = []
    y_local = []
    frames = 0

    for signal in variants:
        for frame in frame_signal(signal):
            if np.max(np.abs(frame)) < ENERGY_THRESHOLD:
                continue

            feat = extract_features(frame)
            x_local.append(feat)
            y_local.append(pitch)
            frames += 1

    if not x_local:
        return None

    return (
        np.vstack(x_local),
        np.asarray(y_local),
        pitch,
        frames,
    )


# ----------------------------
# File listing
# ----------------------------
def list_ordinario_files(root):
    for dirpath, _, filenames in os.walk(root):
        if os.path.basename(dirpath) != "ordinario":
            continue
        for f in filenames:
            if f.lower().endswith((".wav", ".aif", ".aiff")):
                yield os.path.join(dirpath, f)


# ----------------------------
# Dataset build (parallel)
# ----------------------------
def build_dataset(root, workers=None):
    # ---------- cache ----------
    if os.path.exists(CACHE_NAME):
        print(f"Loading cached dataset → {CACHE_NAME}")
        data = np.load(CACHE_NAME, allow_pickle=True)

        x_train = data["x_train"]
        y_train = data["y_train"]
        x_test = data["x_test"]
        y_test = data["y_test"]

        MAX_TRAIN = 200_000
        MAX_TEST = 30_000

        if len(x_train) > MAX_TRAIN:
            idx = np.random.choice(len(x_train), MAX_TRAIN, replace=False)
            x_train = x_train[idx]
            y_train = y_train[idx]

        if len(x_test) > MAX_TEST:
            idx = np.random.choice(len(x_test), MAX_TEST, replace=False)
            x_test = x_test[idx]
            y_test = y_test[idx]

        return x_train, y_train, x_test, y_test

    files = list(list_ordinario_files(root))
    if not files:
        raise RuntimeError("No audio files found inside ordinario folders")

    random.shuffle(files)
    split = int(len(files) * (1 - TEST_FRACTION))

    train_files = files[:split]
    test_files = files[split:]

    x_train, y_train = [], []
    x_test, y_test = [], []

    progress = Progress()
    workers = workers or os.cpu_count()

    for mode, filelist in [("traindata", train_files), ("testdata", test_files)]:
        print(f"\n--- Processing {mode} ({len(filelist)} files) ---")

        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(process_file, (fp, mode)) for fp in filelist]

            for fut in as_completed(futures):
                result = fut.result()
                if result is None:
                    continue

                x, y, pitch, frames = result

                if mode == "traindata":
                    x_train.append(x)
                    y_train.append(y)
                else:
                    x_test.append(x)
                    y_test.append(y)

                progress.files += 1
                progress.frames += frames
                progress.per_pitch[pitch] += frames

                if progress.files % 10 == 0:
                    progress.print_status()

    print("\n--- Dataset summary (frames per pitch) ---")
    for p in sorted(progress.per_pitch):
        print(f"{p:4s}: {progress.per_pitch[p]}")

    x_train = np.vstack(x_train)
    y_train = np.concatenate(y_train)
    x_test = np.vstack(x_test)
    y_test = np.concatenate(y_test)

    print(f"\nSaving dataset cache → {CACHE_NAME}")
    np.savez_compressed(
        CACHE_NAME,
        x_train=x_train,
        y_train=y_train,
        x_test=x_test,
        y_test=y_test,
        fft=FFT,
        hop=HOP,
        sr=SR,
    )

    return x_train, y_train, x_test, y_test


# ----------------------------
# Model
# ----------------------------
def init_model():
    return CatBoostClassifier(
        iterations=ITERATIONS,
        depth=6,
        learning_rate=0.1,
        loss_function="MultiClass",
        random_seed=RANDOM_STATE,
        verbose=100,
        early_stopping_rounds=20,
        thread_count=-1,
    )


# ----------------------------
# Main
# ----------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("orchidea_root", help="Root Orchidea folder")
    parser.add_argument("output", help="Output model (.onnx)")
    parser.add_argument("--workers", type=int, default=None)
    args = parser.parse_args()

    print("Building dataset (parallel, cached)...")
    x_train, y_train, x_test, y_test = build_dataset(
        args.orchidea_root,
        workers=args.workers,
    )

    print(f"\nTrain samples: {len(x_train)}")
    print(f"Test samples : {len(x_test)}")
    print(f"Pitch classes: {len(set(y_train))}")

    clf = init_model()

    print("\nTraining...")
    clf.fit(
        x_train,
        y_train,
        eval_set=(x_test, y_test),
    )

    print("\nEvaluation:")
    y_pred = clf.predict(x_test)
    print(classification_report(y_test, y_pred))

    print(f"\nSaving model → {args.output}")
    clf.save_model(args.output, format="onnx")


if __name__ == "__main__":
    main()
