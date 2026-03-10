import numpy as np
import OpenScofo
import os
import random


import librosa
import pyloudnorm
import essentia
import essentia.standard as es

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


# ---------------- RMS TEST ----------------
def run_test_rms(window, label):
    scofo_desc = scofo.get_audio_description(window)

    # RMS-based loudness (power=2.0 equivalent)
    l_loudness = librosa.feature.rms(
        y=window,
        frame_length=n_fft,
        hop_length=hop,
        center=False,
    )[0, 0]

    s_loudness = scofo_desc.rms
    diff = abs(l_loudness - s_loudness)

    print(
        f"{label} | "
        f"RMS | "
        f"L: {l_loudness:+010.5f} | "
        f"S: {s_loudness:+010.5f} | "
        f"D: {diff:+010.5f}"
    )


# ---------------- ZEROCROSS-RATING TEST --------
def run_test_zcr(window, label):
    scofo_desc = scofo.get_audio_description(window)

    l_zcr = librosa.feature.zero_crossing_rate(
        y=window,
        frame_length=n_fft,
        hop_length=hop,
        center=True,
        threshold=1e-10,
        ref_magnitude=1.0,
        pad=False,
        zero_pos=True,
    )[0, 0]

    s_zcr = scofo_desc.zero_crossing_rate
    diff = abs(l_zcr - s_zcr)

    print(
        f"{label} | "
        f"ZCR | "
        f"L: {l_zcr:+015.10f} | "
        f"S: {s_zcr:+015.10f} | "
        f"D: {diff:+015.10f}"
    )


# TODO: Create loudness validation
# ---------------- LOUDNESS TEST ----------------


# ---------------- TESTES ALEATÓRIOS ----------------
n_tests = 20

for _ in range(n_tests):
    max_start = len(y) - n_fft
    start = random.randint(0, max_start)
    window = y[start : start + n_fft]

    run_test_mfcc(window, f"Start {start:08d}")
    run_test_flatness(window, f"Start {start:08d}")
    run_test_rms(window, f"Start {start:08d}")
    run_test_zcr(window, f"Start {start:08d}")

    print("")
