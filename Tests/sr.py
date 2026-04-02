import os
import time
import numpy as np
import librosa
from concurrent.futures import ProcessPoolExecutor

import OpenScofo
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from pyoperon.sklearn import SymbolicRegressor

# --- Configuration ---
SR = 48000
FFTSIZE = 2048
HOPSIZE = 512
ACTIVE_DB_THRESHOLD = -60.0
CLASSIFICATION_THRESHOLD = 0.5
RANDOM_STATE = 42
MAX_PER_CLASS = 5000
CACHE_FILE = "percussion_descriptor_dataset_first10.npz"
N_BEGIN = 10  
EPS = 1e-6

fields = [
    "onset", "max_amp", "loudness",
    "flux", "irregularity", "crest", "centroid",
    "centroid_velocity", "flatness", "hfr",
    "harmonicity", "zcr", "pitch",
    "pitch_confidence", "db", "rms"
]

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

# --- Extraction Helpers ---
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
    X_temp, y_temp = [], []
    x_audio, _ = librosa.load(file, sr=SR)
    count = 0
    for i in range(0, len(x_audio) - FFTSIZE, HOPSIZE):
        if count >= N_BEGIN:
            break
        frame = x_audio[i : i + FFTSIZE]
        ok = scofo.process_block(frame)
        desc = scofo.get_description()
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
    X, y = [], []
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

# --- Data Pipeline ---
files = list_files()
file_labels = [label for _, label in files]
train_files, test_files = train_test_split(files, test_size=0.2, random_state=RANDOM_STATE, stratify=file_labels)

cached = load_cached()
if cached is None:
    print("Extracting features for first 10 frames...")
    start = time.perf_counter()
    X_train, y_train = extract_dataset(train_files)
    X_test, y_test = extract_dataset(test_files)
    save_cache(np.vstack((X_train, X_test)), np.concatenate((y_train, y_test)))
    print(f"Feature extraction time: {time.perf_counter()-start:.2f}s")
else:
    print("Loading cached dataset")
    X_all, y_all = cached
    X_train, X_test, y_train, y_test = train_test_split(X_all, y_all, test_size=0.2, random_state=RANDOM_STATE, stratify=y_all)

# Balance classes
rng = np.random.default_rng(RANDOM_STATE)
def balance(X, y):
    perc_idx, pitch_idx = np.where(y == 1.0)[0], np.where(y == 0.0)[0]
    perc_sample = rng.choice(perc_idx, min(MAX_PER_CLASS, len(perc_idx)), replace=False)
    pitch_sample = rng.choice(pitch_idx, min(MAX_PER_CLASS, len(pitch_idx)), replace=False)
    idx = np.concatenate((perc_sample, pitch_sample))
    rng.shuffle(idx)
    return X[idx], y[idx]

X_train, y_train = balance(X_train, y_train)
X_test, y_test = balance(X_test, y_test)

# Normalize
fmin, fmax = X_train.min(axis=0), X_train.max(axis=0)
frange = np.maximum(fmax - fmin, EPS)
X_train = np.clip((X_train - fmin) / frange, EPS, 1.0)
X_test = np.clip((X_test - fmin) / frange, EPS, 1.0)

# --- PyOperon Symbolic Regression ---
print("\nStarting PyOperon Symbolic Regression...")
start = time.perf_counter()

clf = SymbolicRegressor(
    allowed_symbols="add,sub,mul,aq,variable,constant", # 'aq' is safe division
    symbolic_mode=True,                 # Forces integer/simple coefficients
    add_model_scale_term=False,         # Disables the massive outer multiplier
    add_model_intercept_term=False,     # Disables the massive outer offset
    optimizer_iterations=0,             # Stops it from fine-tuning decimals
    objectives=['r2', 'length'],        # Punishes overly complex equations
    generations=250,           
    population_size=1000,       
    max_length=15,             
    random_state=RANDOM_STATE,
    n_threads=16,
)

clf.fit(X_train, y_train)
elapsed = time.perf_counter() - start

# Map X1, X2... back to feature names
best_eq = clf.get_model_string(clf.model_)
for i in sorted(range(len(fields)), reverse=True):
    best_eq = best_eq.replace(f"X{i+1}", fields[i])

print(f"Search completed in {elapsed:.2f}s")
print("\nBest equation found by PyOperon:")
print(best_eq)

# Evaluate thresholding at 0.5
train_preds = (clf.predict(X_train) > CLASSIFICATION_THRESHOLD).astype(int)
test_preds = (clf.predict(X_test) > CLASSIFICATION_THRESHOLD).astype(int)

print(f"\nTrain accuracy: {accuracy_score(y_train, train_preds):.4f}")
print(f"Validation accuracy: {accuracy_score(y_test, test_preds):.4f}")

perc_mask = y_test == 1.0
pitch_mask = y_test == 0.0

perc_violations = np.sum(test_preds[perc_mask] == 0)
pitch_violations = np.sum(test_preds[pitch_mask] == 1)

print("\nExtended Technique (Percussive) violations:", perc_violations)
print("Traditional Note (Pitch) violations:", pitch_violations)
