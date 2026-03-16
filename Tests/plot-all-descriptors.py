import OpenScofo
import librosa
import numpy as np
import itertools

# --- Configuration ---
PERCUSIVE_AUDIO_FILE = [
    "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Winds/Flute/key_click/Fl-key_cl-B3-f-N-N.wav",
    "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Winds/Flute/tongue_ram/Fl-tng_ram-C3-mf-N-N.wav",
]

PITCH_AUDIO_FILE = [
    "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Winds/Flute/ordinario/Fl-ord-A4-mf-N-N.wav",
    "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Winds/Flute/flatterzunge/Fl-flatt-A4-mf-N-N.wav",
]

SR = 48000
FFT_SIZE = 2048
HOP_SIZE = 512
N_BEGIN = 10  # <-- first 10 frames
MAX_COMBO = 5

fields = [
    "silence_prob",
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
    "magnitude",
]


# --- Helpers ---
def extract_field(desc_list, field):
    arr = []
    for d in desc_list[:N_BEGIN]:  # <-- take only first 10 frames
        value = getattr(d, field)
        if isinstance(value, (list, np.ndarray)):
            value = np.mean(value)
        arr.append(value)
    arr = np.array(arr, dtype=np.float64)
    min_val, max_val = arr.min(), arr.max()
    if max_val > min_val:
        arr = (arr - min_val) / (max_val - min_val)
    else:
        arr = np.zeros_like(arr)
    return arr


def process_file(file_path):
    y, sr = librosa.load(file_path, sr=SR)
    oscofo = OpenScofo.OpenScofo(SR, FFT_SIZE, HOP_SIZE)
    num_blocks = (len(y) - FFT_SIZE) // HOP_SIZE + 1
    descriptors = []
    for i in range(min(num_blocks, N_BEGIN)):  # <-- only process first 10 blocks
        block = y[i * HOP_SIZE : i * HOP_SIZE + FFT_SIZE]
        descriptors.append(oscofo.get_audio_description(block))
    return {f: extract_field(descriptors, f) for f in fields}


# --- Load all audio data ---
perc_data = [process_file(f) for f in PERCUSIVE_AUDIO_FILE]
pitch_data = [process_file(f) for f in PITCH_AUDIO_FILE]

# --- Test all combinations ---
results = []

for r in range(1, MAX_COMBO + 1):
    combos = list(itertools.combinations(fields, r))
    total = len(combos) * (2**r)
    count = 0
    print(f"Testing combinations of length {r}, total {total}")

    for combo in combos:
        for inv_flags in itertools.product([False, True], repeat=r):
            count += 1
            if count % 1000 == 0 or count == total:
                print(f"Progress: {count}/{total} ({count/total*100:.1f}%)")

            # --- Compute mean start for percussive ---
            perc_vals = []
            for d in perc_data:
                arrays = []
                for f, inv in zip(combo, inv_flags):
                    arrays.append(1 - d[f] if inv else d[f])
                combined = np.prod(np.vstack(arrays), axis=0)
                perc_vals.append(combined.mean())
            perc_mean = np.mean(perc_vals)

            # --- Compute mean start for pitched ---
            pitch_vals = []
            for d in pitch_data:
                arrays = []
                for f, inv in zip(combo, inv_flags):
                    arrays.append(1 - d[f] if inv else d[f])
                combined = np.prod(np.vstack(arrays), axis=0)
                pitch_vals.append(combined.mean())
            pitch_mean = np.mean(pitch_vals)

            # --- Compute score ---
            score = perc_mean - pitch_mean
            results.append(
                {
                    "combo": combo,
                    "inverted": inv_flags,
                    "score": score,
                    "perc_mean": perc_mean,
                    "pitch_mean": pitch_mean,
                }
            )

# --- Sort by score ---
results_sorted = sorted(results, key=lambda x: x["score"], reverse=True)

print("Top combinations for strong percussive / low pitched (first 10 frames):")
for r in results_sorted[:10]:
    print(r)
