import os
import time
import numpy as np
import librosa
from concurrent.futures import ProcessPoolExecutor

import OpenScofo
from pysr import PySRRegressor
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

SR = 48000
FFTSIZE = 2048
HOPSIZE = 512

ACTIVE_DB_THRESHOLD = -60.0
CLASSIFICATION_THRESHOLD = 0.5

RANDOM_STATE = 42
MAX_PER_CLASS = 5000

CACHE_FILE = "percussion_descriptor_dataset.npz"

VIOLATION_PENALTY = 500.0
STRICT_MAX_VIOLATIONS = 0


FEATURE_SPECS = [
    ("Harmonicity", "harmonicity"),
    ("SpectralFlatness", "spectral_flatness"),
    ("SpectralFlux", "spectral_flux"),
    ("SpectralCrest", "spectral_crest"),
    ("CentroidVelocity", "centroid_velocity"),
    ("HighFreqRatio", "high_freq_ratio"),
    ("Peakiness", "peakiness"),
    ("PitchConfidence", "pitch_confidence"),
]

feature_names = [name for name, _ in FEATURE_SPECS]
feature_attrs = [attr for _, attr in FEATURE_SPECS]

base_folder = "/home/neimog/Nextcloud/MusicData/Samples/Orchidea"

folder_labels = {
    f"{base_folder}/Winds/Flute/pizzicato": 1.0,
    f"{base_folder}/Winds/Flute/key_click": 1.0,
    f"{base_folder}/Winds/Flute/tongue_ram": 1.0,
    f"{base_folder}/Winds/Flute/jet_whistle": 1.0,
    f"{base_folder}/Winds/Flute/ordinario": 0.0,
    f"{base_folder}/Winds/Flute/aeolian": 0.0,
    f"{base_folder}/Winds/Flute/staccato": 0.0,
    f"{base_folder}/Strings/Violoncello/pizzicato_secco": 1.0,
    f"{base_folder}/Strings/Violoncello/hit_on_body": 1.0,
    f"{base_folder}/Strings/Violin/pizzicato_secco": 1.0,
    f"{base_folder}/Winds/Bassoon/key_click": 1.0,
    f"{base_folder}/Winds/Clarinet_Bb/key_click": 1.0,
    f"{base_folder}/Winds/Oboe/key_click": 1.0,
}


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
    return [float(getattr(desc, attr)) for attr in feature_attrs]


def process_file(args):
    file, label = args
    scofo = OpenScofo.OpenScofo(SR, FFTSIZE, HOPSIZE)
    X_temp = []
    y_temp = []
    x_audio, _ = librosa.load(file, sr=SR)
    for i in range(0, len(x_audio) - FFTSIZE, HOPSIZE):
        frame = x_audio[i : i + FFTSIZE]
        desc = scofo.get_audio_description(frame)
        if desc.db <= ACTIVE_DB_THRESHOLD:
            continue
        features = build_feature_vector(desc)
        if not np.all(np.isfinite(features)):
            continue
        X_temp.append(features)
        y_temp.append(label)
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
    if cache["feature_names"].tolist() != feature_names:
        return None
    return cache["X"], cache["y"]


def save_cache(X, y):
    np.savez_compressed(CACHE_FILE, X=X, y=y, feature_names=np.array(feature_names))


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
    print("Extracting features...")
    start = time.perf_counter()
    X_train, y_train = extract_dataset(train_files)
    X_test, y_test = extract_dataset(test_files)
    save_cache(
        np.vstack((X_train, X_test)),
        np.concatenate((y_train, y_test)),
    )
    end = time.perf_counter()
    print(f"Feature extraction time: {end-start:.2f}s")

else:
    print("Loading cached dataset")
    X_all, y_all = cached
    X_train, X_test, y_train, y_test = train_test_split(
        X_all,
        y_all,
        test_size=0.2,
        random_state=RANDOM_STATE,
        stratify=y_all,
    )


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

mean = X_train.mean(axis=0)
std = X_train.std(axis=0) + 1e-8

X_train = (X_train - mean) / std
X_test = (X_test - mean) / std

print("Training frames:", len(X_train))
print("Validation frames:", len(X_test))

constraints = {"const": lambda x: -0.01 <= x <= 0.01}


model = PySRRegressor(
    niterations=150,
    populations=12,
    maxsize=18,
    binary_operators=["*", "/"],
    unary_operators=["exp", "log"],  # allows some non-linearity
    model_selection="best",
    elementwise_loss="loss(x, y) = sum((x - y)^2)",  # or your ELEMENTWISE_LOSS
    loss_scale="linear",
    parsimony=1e-3,
    parallelism="serial",
    deterministic=True,
    random_state=RANDOM_STATE,
)

model.fit(
    X_train,
    y_train,
    variable_names=feature_names,
)

print("\nBest equation:")
print(model.sympy())

pred = model.predict(X_test)
pred = np.clip(pred, 0, 1)
binary = (pred > CLASSIFICATION_THRESHOLD).astype(int)
acc = accuracy_score(y_test, binary)
print("\nValidation accuracy:", acc)

perc_mask = y_test == 1
pitch_mask = y_test == 0

perc_violations = np.sum(pred[perc_mask] < CLASSIFICATION_THRESHOLD)
pitch_violations = np.sum(pred[pitch_mask] > CLASSIFICATION_THRESHOLD)

print("Percussive violations:", perc_violations)
print("Pitch violations:", pitch_violations)


def validate_files(files, label_name, n_frames=5):
    print(f"=== Random {label_name} Sounds ====")
    scofo = OpenScofo.OpenScofo(SR, FFTSIZE, HOPSIZE)
    max_name_len = 50  # adjust width if needed
    for file, label in rng.choice(files, size=5, replace=False):
        x_audio, _ = librosa.load(file, sr=SR)
        frame_indices = np.arange(0, len(x_audio) - FFTSIZE, HOPSIZE)
        rng.shuffle(frame_indices)
        preds = []
        count = 0
        for idx in frame_indices:
            frame = x_audio[idx : idx + FFTSIZE]
            desc = scofo.get_audio_description(frame)
            if desc.db <= ACTIVE_DB_THRESHOLD:
                continue
            features = build_feature_vector(desc)
            if not np.all(np.isfinite(features)):
                continue
            features = (np.array(features) - mean) / std
            pred = model.predict([features])[0]
            preds.append(pred)
            count += 1
            if count >= n_frames:
                break
        if len(preds) == 0:
            print(f"{os.path.basename(file):{max_name_len}} | no valid frames")
            continue
        avg_pred = np.mean(preds)
        print(
            f"{os.path.basename(file):{max_name_len}} | avg result eq: {avg_pred:.4f}"
        )


# select percussive and pitch files separately
perc_files = [f for f in files if f[1] == 1.0]
pitch_files = [f for f in files if f[1] == 0.0]

validate_files(perc_files, "Percussive")
validate_files(pitch_files, "Pitch")

