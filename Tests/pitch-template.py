import OScofo
import numpy as np
import librosa
import music21
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from skopt import gp_minimize
from skopt.space import Real, Integer

# ----------------------------
# Configuration
# ----------------------------
SR = 48000
FFT_SIZE = 4096
HOP_SIZE = 1024
N_WORKERS = os.cpu_count()

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# ----------------------------
# Dataset
# ----------------------------
root = Path("/home/neimog/Nextcloud/MusicData/Samples/Orchidea")
files = sorted(p for p in root.rglob("ordinario/**/*.wav"))


# ----------------------------
# Audio analysis helper
# ----------------------------
def get_audioanalisys(path, scofo):
    y, _ = librosa.load(path, sr=SR)
    start = round(SR * 0.15)
    segment = y[start : start + FFT_SIZE]
    return scofo.get_audio_description(segment)


# ----------------------------
# Worker (one file)
# ----------------------------
def process_file(args):
    f, amplitude_decay, harmonics, sigma = args

    scofo = OScofo.OScofo(SR, FFT_SIZE, HOP_SIZE)
    scofo.parse_score("../Resources/tests/canticos.txt")
    scofo.set_amplitude_decay(amplitude_decay)
    scofo.set_harmonics(int(harmonics))
    scofo.set_pitch_template_sigma(sigma)

    desc = get_audioanalisys(f, scofo)
    spectral = np.asarray(desc.norm_spectral_power)

    max_result = -np.inf
    max_note = None

    for midi in range(24, 96):
        note = music21.note.Note(midi)
        freq = note.pitch.freq440
        pitch_template = np.asarray(scofo.get_pitch_template(freq))

        peak_p = np.max(np.abs(spectral))
        peak_q = np.max(np.abs(pitch_template))
        if peak_p > 0 and peak_q > 0:
            pitch_template *= peak_p / peak_q

        eps = 1e-24
        p = np.maximum(spectral, 0) + eps
        q = np.maximum(pitch_template, 0) + eps
        p /= np.sum(p)
        q /= np.sum(q)

        kl = np.sum(p * np.log(p / q))
        score = np.exp(-0.5 * kl)

        if score > max_result:
            max_result = score
            max_note = note

    gt = os.path.basename(f).split("-")[2]
    return int(max_note.nameWithOctave != gt)


# ----------------------------
# Evaluation (parallel)
# ----------------------------
def evaluate_params(amplitude_decay, harmonics, sigma):
    args = [(f, amplitude_decay, harmonics, sigma) for f in files]

    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        misses = sum(ex.map(process_file, args))

    print(
        f"misses={misses:4d} | "
        f"amp_decay={amplitude_decay:.3f} "
        f"harmonics={harmonics:2d} "
        f"sigma={sigma:.3f}"
    )

    return misses


# ----------------------------
# Bayesian optimization
# ----------------------------
space = [
    Real(0.1, 2.0),  # Amplitude decay (allow very slow decay for brass/strings)
    Integer(5, 20),  # Harmonics (Real instruments need >10 harmonics)
    Real(0.5, 2.0),  # Sigma (Width in semitones/log scale needs room to breathe)
]

if __name__ == "__main__":
    result = gp_minimize(
        func=lambda x: evaluate_params(x[0], x[1], x[2]),
        dimensions=space,
        n_calls=40,
        n_initial_points=10,
        random_state=0,
    )

    print("\n===== BEST PARAMETERS =====")
    print("Amplitude decay      :", result.x[0])
    print("Harmonics            :", result.x[1])
    print("Pitch template sigma :", result.x[2])
    print("Misses               :", result.fun)
