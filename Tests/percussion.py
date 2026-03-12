import os
import time
import numpy as np
import librosa
import OpenScofo

from pysr import PySRRegressor
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split

SR = 48000
FFTSIZE = 2048
HOPSIZE = 512
CLASSIFICATION_THRESHOLD = 0.3
MARGIN_PENALTY_WEIGHT = 30.0
RANGE_PENALTY_WEIGHT = 0.5

os.chdir(os.path.dirname(__file__))
base_folder = "/home/neimog/Nextcloud/MusicData/Samples/Orchidea"

folder_labels = {
    f"{base_folder}/Winds/Flute/pizzicato": "e",
    f"{base_folder}/Winds/Flute/key_click": "e",
    f"{base_folder}/Winds/Flute/tongue_ram": "e",
    f"{base_folder}/Winds/Flute/jet_whistle": "e",
    f"{base_folder}/Winds/Flute/ordinario": "o",
    f"{base_folder}/Winds/Flute/ordinario_quartertones": "o",
    f"{base_folder}/Winds/Flute/flatterzunge_to_ordinario": "o",
    f"{base_folder}/Winds/Flute/aeolian": "o",
    f"{base_folder}/Winds/Flute/discolored_fingering": "o",
    f"{base_folder}/Winds/Flute/sforzato": "o",
    f"{base_folder}/Winds/Flute/play_and_sing_unison": "o",
    f"{base_folder}/Winds/Flute/whistle_tones": "o",
    f"{base_folder}/Winds/Flute/decrescendo": "o",
    f"{base_folder}/Winds/Flute/ordinario_to_flatterzunge": "o",
    f"{base_folder}/Winds/Flute/note_lasting": "o",
    f"{base_folder}/Winds/Flute/play_and_sing": "o",
    f"{base_folder}/Winds/Flute/crescendo": "o",
    f"{base_folder}/Winds/Flute/crescendo_to_decrescendo": "o",
    f"{base_folder}/Winds/Flute/aeolian_and_ordinario": "o",
    f"{base_folder}/Winds/Flute/aeolian_to_ordinario": "o",
    f"{base_folder}/Winds/Flute/chromatic_scale": "o",
    f"{base_folder}/Winds/Flute/ordinario_to_aeolian": "o",
    f"{base_folder}/Winds/Flute/trill_minor_second_up": "o",
    f"{base_folder}/Winds/Flute/staccato": "o",
    f"{base_folder}/Brass/Bass_Tuba/slap_pitched": "e",
    f"{base_folder}/Brass/Trombone/slap_pitched": "e",
    f"{base_folder}/Keyboards/Accordion/breath/Acc-breath-N-mf-2-N.wav": "e",
    f"{base_folder}/Strings/Violoncello/pizzicato_secco": "e",
    f"{base_folder}/Strings/Violoncello/pizzicato_bartok/Vc-pizz_bartok-B3-ff-1c-N.wav": "e",
    f"{base_folder}/Strings/Violoncello/hit_on_body": "e",
    f"{base_folder}/Strings/Violin/pizzicato_secco": "e",
    f"{base_folder}/Winds/Bassoon/key_click": "e",
    f"{base_folder}/Winds/Clarinet_Bb/key_click": "e",
    f"{base_folder}/Winds/Oboe/key_click": "e",
    f"{base_folder}/Winds/Oboe/kiss": "e",
    f"{base_folder}/Winds/Sax_Alto/exploding_slap_pitched": "e",
    f"{base_folder}/Winds/Sax_Alto/key_click": "e",
    f"{base_folder}/Winds/Sax_Alto/slap_pitched": "e",
    f"{base_folder}/Winds/Sax_Alto/slap_unpitched": "e",
}

files_percussive = []
files_pitch = []

for root, label in folder_labels.items():
    if os.path.isdir(root):
        wav_files = [
            os.path.join(root, f)
            for f in os.listdir(root)
            if f.lower().endswith(".wav")
        ]
    elif os.path.isfile(root) and root.lower().endswith(".wav"):
        wav_files = [root]
    else:
        continue

    if label == "e":
        files_percussive.extend(wav_files)
    elif label == "o":
        files_pitch.extend(wav_files)

if len(files_percussive) == 0 or len(files_pitch) == 0:
    raise Exception("Invalid audio files")

scofo = OpenScofo.OpenScofo(SR, FFTSIZE, HOPSIZE)

feature_names = [
    "Flatness",
    "SpectralFlux",
    "Harmonicity",
    "ZeroCrossingRate",
    "DeltaRMS",
    "SpectralIrregularity",
    "SpectralCrest",
    "CentroidVelocity",
    "HighFreqRatio",
    "Peakiness",
    "SilenceGate",
    "FluxOverHarmonicity",
    "InstPercSeed",
]


def build_feature_vector(desc, delta_rms):
    silence_gate = 1.0 - desc.silence_prob
    flux_over_harmonicity = desc.spectral_flux / (desc.harmonicity + 1e-6)
    inst_perc_seed = 0.21974495 * flux_over_harmonicity * silence_gate

    return [
        desc.spectral_flatness,
        desc.spectral_flux,
        desc.harmonicity,
        desc.zero_crossing_rate,
        delta_rms,
        desc.spectral_irregularity,
        desc.spectral_crest,
        desc.centroid_velocity,
        desc.high_freq_ratio,
        desc.peakiness,
        silence_gate,
        flux_over_harmonicity,
        inst_perc_seed,
    ]


def extract_features(file_list, target_label):
    X_temp, y_temp = [], []
    for file in file_list:
        x_audio, sr = librosa.load(file, sr=SR)
        frames = librosa.util.frame(x_audio, frame_length=FFTSIZE, hop_length=HOPSIZE).T
        prev_rms = 0.0

        for frame in frames:
            desc = scofo.get_audio_description(frame)
            delta_rms = max(0.0, desc.rms - prev_rms)

            if desc.db > -60.0:
                X_temp.append(build_feature_vector(desc, delta_rms))
                y_temp.append(target_label)

            prev_rms = desc.rms

    return X_temp, y_temp


X_file = "X.npy"
y_file = "y.npy"

if os.path.exists(X_file) and os.path.exists(y_file):
    print("Loading cached dataset...")
    X = np.load(X_file)
    y = np.load(y_file)
else:
    print("Extracting features...")
    start = time.perf_counter()

    X_perc, y_perc = extract_features(files_percussive, 1.0)
    X_pitch, y_pitch = extract_features(files_pitch, 0.0)

    end = time.perf_counter()
    print(f"Elapsed time: {end - start:.3f} seconds")

    X = np.array(X_perc + X_pitch)
    y = np.array(y_perc + y_pitch)

    np.save(X_file, X)
    np.save(y_file, y)

print(f"Original dataset size: {len(X)}")

MAX_PER_CLASS = 5000

perc_idx = np.where(y == 1.0)[0]
pitch_idx = np.where(y == 0.0)[0]

perc_sample = np.random.choice(
    perc_idx, min(MAX_PER_CLASS, len(perc_idx)), replace=False
)
pitch_sample = np.random.choice(
    pitch_idx, min(MAX_PER_CLASS, len(pitch_idx)), replace=False
)

idx = np.concatenate([perc_sample, pitch_sample])
np.random.shuffle(idx)

X = X[idx].astype(np.float32)
y = y[idx].astype(np.float32)

print(f"Extracted {len(X)} active frames for training.")

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42,
    stratify=y,
)

print(f"Training frames: {len(X_train)} | Validation frames: {len(X_test)}")

elementwise_loss = (
    "loss(prediction, target) = begin "
    f"threshold = {CLASSIFICATION_THRESHOLD:.8f}f0; "
    "p = clamp(prediction, 0.0f0, 1.0f0); "
    "range_penalty = (prediction - p)^2; "
    "perc_violation = target * max(0.0f0, threshold - p)^2; "
    "pitch_violation = (1.0f0 - target) * max(0.0f0, p - threshold)^2; "
    "target_penalty = (p - target)^2; "
    f"return target_penalty + {MARGIN_PENALTY_WEIGHT:.8f}f0 * (perc_violation + pitch_violation) + "
    f"{RANGE_PENALTY_WEIGHT:.8f}f0 * range_penalty; "
    "end"
)

# Setup Symbolic Regression
model = PySRRegressor(
    niterations=40,
    binary_operators=["+", "*", "-", "/"],
    unary_operators=[
        "exp",
        "log",
    ],
    # Force it to find simple equations, not complex unreadable ones
    # Use logistic loss to map the output equation to probabilities [0, 1]
    elementwise_loss="loss(prediction, target) = (prediction - target)^2",
)

# Run the evolutionary search
model.fit(X_train, y_train, variable_names=feature_names)

# Print the discovered equations
print("\n=== Best Discovered Equations ===")
print(model.sympy())

print("\n=== Classification Accuracy ===")

equations = model.equations_
if equations is None:
    raise Exception("Not found equations")


def summarize_predictions(y_true, y_pred):
    y_pred_clipped = np.clip(y_pred, 0.0, 1.0)
    y_pred_binary = (y_pred_clipped > CLASSIFICATION_THRESHOLD).astype(int)

    acc = accuracy_score(y_true, y_pred_binary)

    perc_mask = y_true == 1
    pitch_mask = y_true == 0
    perc_acc = np.mean(y_pred_binary[perc_mask] == 1)
    pitch_acc = np.mean(y_pred_binary[pitch_mask] == 0)
    perc_mean = float(np.mean(y_pred_clipped[perc_mask]))
    pitch_mean = float(np.mean(y_pred_clipped[pitch_mask]))
    perc_min = float(np.min(y_pred_clipped[perc_mask]))
    pitch_max = float(np.max(y_pred_clipped[pitch_mask]))
    perc_violations = int(np.sum(y_pred_clipped[perc_mask] < CLASSIFICATION_THRESHOLD))
    pitch_violations = int(
        np.sum(y_pred_clipped[pitch_mask] > CLASSIFICATION_THRESHOLD)
    )

    return (
        acc,
        perc_acc,
        pitch_acc,
        perc_mean,
        pitch_mean,
        perc_min,
        pitch_max,
        perc_violations,
        pitch_violations,
    )


for i in range(len(equations)):
    row = equations.iloc[i]

    train_pred = model.predict(X_train, index=i)
    test_pred = model.predict(X_test, index=i)

    (
        train_acc,
        train_perc_acc,
        train_pitch_acc,
        train_perc_mean,
        train_pitch_mean,
        train_perc_min,
        train_pitch_max,
        train_perc_violations,
        train_pitch_violations,
    ) = summarize_predictions(y_train, train_pred)
    (
        test_acc,
        test_perc_acc,
        test_pitch_acc,
        test_perc_mean,
        test_pitch_mean,
        test_perc_min,
        test_pitch_max,
        test_perc_violations,
        test_pitch_violations,
    ) = summarize_predictions(y_test, test_pred)

    hard_pass = test_perc_violations == 0 and test_pitch_violations == 0

    print(
        f"Equation {i:2d} | Complexity {row['complexity']:2d} | "
        f"Train Global: {train_acc*100:.2f}% | "
        f"Train Perc: {train_perc_acc*100:.2f}% | "
        f"Train Pitch: {train_pitch_acc*100:.2f}% | "
        f"Train Bounds => perc min: {train_perc_min:.3f}, pitch max: {train_pitch_max:.3f} | "
        f"Train Violations => perc: {train_perc_violations}, pitch: {train_pitch_violations} | "
        f"Val Global: {test_acc*100:.2f}% | "
        f"Val Perc: {test_perc_acc*100:.2f}% | "
        f"Val Pitch: {test_pitch_acc*100:.2f}% | "
        f"Val Means => perc: {test_perc_mean:.3f}, pitch: {test_pitch_mean:.3f} | "
        f"Val Bounds => perc min: {test_perc_min:.3f}, pitch max: {test_pitch_max:.3f} | "
        f"Val Violations => perc: {test_perc_violations}, pitch: {test_pitch_violations} | "
        f"HardPass: {hard_pass}"
    )
