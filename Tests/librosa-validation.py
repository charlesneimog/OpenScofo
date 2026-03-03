import numpy as np
import librosa
import OpenScofo
import os
import random

os.chdir(os.path.dirname(__file__))

sr = 48000
n_fft = 2048
hop = 512

scofo = OpenScofo.OpenScofo(sr, n_fft, hop)

y, sr = librosa.load(
    "./assets/bwv-1013.mp3",
    sr=sr,
)


# ---------------- MFCC TEST ----------------
def run_test_mfcc(window, label):

    scofo_desc = scofo.get_audio_description(window)

    mfcc = librosa.feature.mfcc(
        y=window,
        sr=sr,
        n_mfcc=13,
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
        n_mels=40,
        power=2.0,
    )

    librosa_mfcc = mfcc[:, 0].tolist()
    scofo_mfcc = scofo_desc.mfcc

    max_diff = 0.0
    max_idx = -1
    max_vals = (0.0, 0.0)

    for i, (l, s) in enumerate(zip(librosa_mfcc, scofo_mfcc)):
        abs_diff = abs(l - s)
        if abs_diff > max_diff:
            max_diff = abs_diff
            max_idx = i
            max_vals = (l, s)

    print(
        f"{label} | "
        f"MFCC idx {max_idx:02d} | "
        f"L: {max_vals[0]:+010.5f} | "
        f"S: {max_vals[1]:+010.5f} | "
        f"D: {max_diff:+010.5f}"
    )


# ---------------- FLATNESS TEST ----------------
def run_test_flatness(window, label):

    scofo_desc = scofo.get_audio_description(window)

    l_flatness = librosa.feature.spectral_flatness(
        y=window,
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
        power=2.0,
    )[0, 0]

    s_flatness = scofo_desc.spectral_flatness
    diff = abs(l_flatness - s_flatness)

    print(
        f"{label} | "
        f"FLAT | "
        f"L: {l_flatness:+010.5f} | "
        f"S: {s_flatness:+010.5f} | "
        f"D: {diff:+010.5f}"
    )


# ---------------- SILÊNCIO ----------------
silence = np.zeros(n_fft, dtype=np.float32)

run_test_mfcc(silence, "SILENCE")
run_test_flatness(silence, "SILENCE")


# ---------------- TESTES ALEATÓRIOS ----------------
n_tests = 20

for _ in range(n_tests):
    max_start = len(y) - n_fft
    start = random.randint(0, max_start)
    window = y[start : start + n_fft]

    run_test_mfcc(window, f"Start {start:08d}")
    run_test_flatness(window, f"Start {start:08d}")
    print("")
