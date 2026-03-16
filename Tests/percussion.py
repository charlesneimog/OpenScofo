import os
import time
import itertools
import math
import numpy as np
import librosa
from concurrent.futures import ProcessPoolExecutor

import OpenScofo
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

# --- Configuration ---
SR = 48000
FFTSIZE = 2048
HOPSIZE = 512
ACTIVE_DB_THRESHOLD = -60.0
CLASSIFICATION_THRESHOLD = 0.5
RANDOM_STATE = 42
MAX_PER_CLASS = 5000
CACHE_FILE = "percussion_descriptor_dataset_first10.npz"
N_BEGIN = 10  # first 10 frames
MIN_FEATURES_IN_EQUATION = 1
MAX_FEATURES_IN_EQUATION = 4
TOP_K_EQUATIONS = 10
EPS = 1e-6
PROGRESS_EVERY = 1000

# --- Descriptor set (max 8 as requested) ---
fields = [
    "onset",
    "silence_prob",
    "percussive_prob",
    "max_amp",
    "loudness",
    "spectral_flux",
    "spectral_irregularity",
    "spectral_crest",
    "spectral_centroid",
    "centroid_velocity",
    "spectral_spread",
    "spectral_flatness",
    "high_freq_ratio",
    "harmonicity",
    "peakiness",
    "zero_crossing_rate",
    "std_dev",
    "pitch",
    "pitch_confidence",
    "db",
    "rms",
]

# --- Audio folders and labels ---
base_folder = "/home/neimog/Nextcloud/MusicData/Samples/Orchidea"
folder_labels = {
    f"{base_folder}/Winds/Flute/pizzicato": 1.0,
    f"{base_folder}/Winds/Flute/key_click": 1.0,
    f"{base_folder}/Winds/Flute/tongue_ram": 1.0,
    f"{base_folder}/Winds/Flute/jet_whistle": 1.0,
    f"{base_folder}/Strings/Violoncello/pizzicato_secco": 1.0,
    f"{base_folder}/Winds/Flute/ordinario": 0.0,
    f"{base_folder}/Strings/Violoncello/ordinario": 0.0,
    f"{base_folder}/Strings/Violin/ordinario": 0.0,
}


# --- Helpers ---
def list_files():
    files = []
    for root, label in folder_labels.items():
        if os.path.isdir(root):
            for f in os.listdir(root):
                if f.lower().endswith(".wav"):
                    files.append((os.path.join(root, f), label))
        elif os.path.isfile(root):
            files.append((root, label))
    if len(files) == 0:
        raise RuntimeError("No audio files found")
    return files


def build_feature_vector(desc):
    return [float(getattr(desc, f)) for f in fields]


def process_file(args):
    file, label = args
    scofo = OpenScofo.OpenScofo(SR, FFTSIZE, HOPSIZE)
    X_temp = []
    y_temp = []
    x_audio, _ = librosa.load(file, sr=SR)
    count = 0
    for i in range(0, len(x_audio) - FFTSIZE, HOPSIZE):
        if count >= N_BEGIN:  # only first N_BEGIN frames
            break
        frame = x_audio[i : i + FFTSIZE]
        desc = scofo.get_audio_description(frame)
        if desc.db <= ACTIVE_DB_THRESHOLD:
            continue
        features = build_feature_vector(desc)
        if not np.all(np.isfinite(features)):
            continue
        X_temp.append(features)
        y_temp.append(label)
        count += 1
    return X_temp, y_temp


def extract_dataset(files):
    X = []
    y = []
    with ProcessPoolExecutor() as executor:
        for X_part, y_part in executor.map(process_file, files):
            X.extend(X_part)
            y.extend(y_part)
    return np.array(X), np.array(y)


def load_cached():
    if not os.path.exists(CACHE_FILE):
        return None
    cache = np.load(CACHE_FILE, allow_pickle=False)
    if cache["feature_names"].tolist() != fields:
        return None
    return cache["X"], cache["y"]


def save_cache(X, y):
    np.savez_compressed(CACHE_FILE, X=X, y=y, feature_names=np.array(fields))


# --- Load dataset ---
files = list_files()
file_labels = [label for _, label in files]
train_files, test_files = train_test_split(
    files,
    test_size=0.2,
    random_state=RANDOM_STATE,
    stratify=file_labels,
)

cached = load_cached()
if cached is None:
    print("Extracting features for first 10 frames...")
    start = time.perf_counter()
    X_train, y_train = extract_dataset(train_files)
    X_test, y_test = extract_dataset(test_files)
    save_cache(np.vstack((X_train, X_test)), np.concatenate((y_train, y_test)))
    end = time.perf_counter()
    print(f"Feature extraction time: {end-start:.2f}s")
else:
    print("Loading cached dataset")
    X_all, y_all = cached
    X_train, X_test, y_train, y_test = train_test_split(
        X_all, y_all, test_size=0.2, random_state=RANDOM_STATE, stratify=y_all
    )

# --- Balance classes ---
rng = np.random.default_rng(RANDOM_STATE)


def balance(X, y):
    perc_idx = np.where(y == 1.0)[0]
    pitch_idx = np.where(y == 0.0)[0]
    perc_sample = rng.choice(perc_idx, min(MAX_PER_CLASS, len(perc_idx)), replace=False)
    pitch_sample = rng.choice(
        pitch_idx, min(MAX_PER_CLASS, len(pitch_idx)), replace=False
    )
    idx = np.concatenate((perc_sample, pitch_sample))
    rng.shuffle(idx)
    return X[idx], y[idx]


X_train, y_train = balance(X_train, y_train)
X_test, y_test = balance(X_test, y_test)

# --- Normalize to [0, 1] for stable multiplication/inversion ---
fmin = X_train.min(axis=0)
fmax = X_train.max(axis=0)
frange = np.maximum(fmax - fmin, EPS)
X_train = np.clip((X_train - fmin) / frange, EPS, 1.0)
X_test = np.clip((X_test - fmin) / frange, EPS, 1.0)


def normalize_scores(scores):
    smin = float(scores.min())
    smax = float(scores.max())
    if smax - smin < EPS:
        return np.zeros_like(scores)
    return (scores - smin) / (smax - smin)


def term_to_string(index, use_complement):
    name = fields[index]
    if use_complement:
        return f"(1 - {name})"
    return name


def equation_to_string(indices, transforms, signs):
    numerator = [
        term_to_string(i, t == 1)
        for i, t, s in zip(indices, transforms, signs)
        if s == 1
    ]
    denominator = [
        term_to_string(i, t == 1)
        for i, t, s in zip(indices, transforms, signs)
        if s == -1
    ]

    if not numerator:
        num_str = "1"
    else:
        num_str = " * ".join(numerator)

    if denominator:
        den_str = " * ".join(denominator)
        return f"({num_str}) / ({den_str})"
    return num_str


def search_best_equations(X_train, y_train, min_terms, max_terms, top_k):
    n_features = X_train.shape[1]
    best = []

    total_candidates = 0
    total_to_test = 0
    for n_terms in range(min_terms, max_terms + 1):
        total_to_test += (
            math.comb(n_features, n_terms)
            * (2**n_terms)  # x or (1 - x) per selected feature
            * ((2**n_terms) - 1)  # numerator/denominator assignment
        )

    start = time.perf_counter()

    for n_terms in range(min_terms, max_terms + 1):
        for combo in itertools.combinations(range(n_features), n_terms):
            combo_values = X_train[:, combo]
            for transforms in itertools.product([0, 1], repeat=n_terms):
                transformed_values = np.empty_like(combo_values)
                for j, use_complement in enumerate(transforms):
                    if use_complement == 1:
                        transformed_values[:, j] = 1.0 - combo_values[:, j]
                    else:
                        transformed_values[:, j] = combo_values[:, j]

                transformed_values = np.clip(transformed_values, EPS, 1.0)
                log_values = np.log(transformed_values)

                for signs in itertools.product([-1, 1], repeat=n_terms):
                    if all(s == -1 for s in signs):
                        continue

                    total_candidates += 1
                    signed = np.array(signs, dtype=np.float64)
                    raw = log_values @ signed
                    pred = normalize_scores(raw)

                    mse = float(np.mean((pred - y_train) ** 2))
                    binary = (pred > CLASSIFICATION_THRESHOLD).astype(int)
                    acc = float(accuracy_score(y_train, binary))

                    best.append((mse, -acc, combo, transforms, signs))
                    if len(best) > (top_k * 8):
                        best.sort(key=lambda item: (item[0], item[1]))
                        best = best[: top_k * 4]

                    if (total_candidates % PROGRESS_EVERY) == 0:
                        elapsed_now = time.perf_counter() - start
                        progress = (100.0 * total_candidates) / max(total_to_test, 1)
                        print(
                            f"Progress: {total_candidates}/{total_to_test} "
                            f"({progress:.1f}%) - elapsed {elapsed_now:.1f}s"
                        )

    best.sort(key=lambda item: (item[0], item[1]))
    elapsed = time.perf_counter() - start
    return best[:top_k], total_candidates, elapsed


print(
    f"\nBrute force search with multiplication + inversion only "
    f"(terms {MIN_FEATURES_IN_EQUATION}..{MAX_FEATURES_IN_EQUATION})"
)
top_equations, tested_count, search_time = search_best_equations(
    X_train,
    y_train,
    MIN_FEATURES_IN_EQUATION,
    MAX_FEATURES_IN_EQUATION,
    TOP_K_EQUATIONS,
)
print(f"Tested {tested_count} equations in {search_time:.2f}s")

best_mse, best_neg_acc, best_combo, best_transforms, best_signs = top_equations[0]
best_acc_train = -best_neg_acc
print("\nBest equation:")
best_eq = equation_to_string(best_combo, best_transforms, best_signs)
print(best_eq)
print(f"Train MSE: {best_mse:.6f}")
print(f"Train accuracy: {best_acc_train:.4f}")

print("\nTop equations:")
for rank, (mse, neg_acc, combo, transforms, signs) in enumerate(top_equations, start=1):
    eq = equation_to_string(combo, transforms, signs)
    print(f"{rank:02d}. mse={mse:.6f} acc={-neg_acc:.4f} :: {eq}")

# --- Evaluate on test set using best train equation ---
best_values_test = X_test[:, best_combo].copy()
for j, use_complement in enumerate(best_transforms):
    if use_complement == 1:
        best_values_test[:, j] = 1.0 - best_values_test[:, j]
best_values_test = np.clip(best_values_test, EPS, 1.0)
raw_test = np.log(best_values_test) @ np.array(best_signs)
pred = normalize_scores(raw_test)
binary = (pred > CLASSIFICATION_THRESHOLD).astype(int)
acc = accuracy_score(y_test, binary)
print("\nValidation accuracy:", acc)

perc_mask = y_test == 1
pitch_mask = y_test == 0
perc_violations = np.sum(pred[perc_mask] < CLASSIFICATION_THRESHOLD)
pitch_violations = np.sum(pred[pitch_mask] > CLASSIFICATION_THRESHOLD)
print("Percussive violations:", perc_violations)
print("Pitch violations:", pitch_violations)
