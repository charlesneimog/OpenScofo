import os
import random

import OpenScofo
import librosa

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

    # print(librosa_mfcc)
    # print(scofo_mfcc)

    for i, (l, s) in enumerate(zip(librosa_mfcc, scofo_mfcc)):
        abs_diff = abs(l - s)
        if abs_diff > max_diff:
            max_diff = abs_diff
            max_idx = i
            max_vals = (l, s)

    print(
        f"{label} | "
        f"MFCC | "
        f"L: {max_vals[0]:+012.5f} | "
        f"S: {max_vals[1]:+012.5f} | "
        f"D: {max_diff:+012.5f}"
    )


# ---------------- Centroid TEST ----------------
def run_test_centroid(window, label):
    scofo_desc = scofo.get_audio_description(window)

    # Librosa spectral centroid
    l_centroid = librosa.feature.spectral_centroid(
        y=window,
        sr=sr,  # make sure sr is defined
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
    )[
        0, 0
    ]  # first frame

    # SCOFO spectral centroid
    s_centroid = scofo_desc.spectral_centroid
    diff = abs(l_centroid - s_centroid)

    print(
        f"{label} | "
        f"CENT | "
        f"L: {l_centroid:+012.5f} | "
        f"S: {s_centroid:+012.5f} | "
        f"D: {diff:+012.5f}"
    )


# ---------------- Spread TEST ----------------
def run_test_spread(window, label):
    scofo_desc = scofo.get_audio_description(window)

    # Librosa spectral spread
    l_spread = librosa.feature.spectral_bandwidth(
        y=window,
        sr=sr,  # make sure sr is defined
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
        p=2,  # second-order → spectral spread
    )[
        0, 0
    ]  # first frame

    # SCOFO spectral spread
    s_spread = scofo_desc.spectral_spread
    diff = abs(l_spread - s_spread)

    print(
        f"{label} | "
        f"SPRE | "
        f"L: {l_spread:+012.5f} | "
        f"S: {s_spread:+012.5f} | "
        f"D: {diff:+012.5f}"
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
        f"L: {l_flatness:+012.5f} | "
        f"S: {s_flatness:+012.5f} | "
        f"D: {diff:+012.5f}"
    )


# ---------------- CHROMA TEST ----------------
def run_test_chroma(window, label):
    scofo_desc = scofo.get_audio_description(window)
    librosa_chroma = librosa.feature.chroma_stft(
        y=window,
        sr=sr,
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
        n_chroma=12,
        tuning=0.0,
        norm=None,
    )[:, 0].tolist()

    scofo_chroma = scofo_desc.chroma

    max_diff = 0.0
    max_vals = (0.0, 0.0)

    for i, (l_val, s_val) in enumerate(zip(librosa_chroma, scofo_chroma)):
        abs_diff = abs(l_val - s_val)
        if abs_diff > max_diff:
            max_diff = abs_diff
            max_vals = (l_val, s_val)

    print(
        f"{label} | "
        f"CHRO | "
        f"L: {max_vals[0]:+012.5f} | "
        f"S: {max_vals[1]:+012.5f} | "
        f"D: {max_diff:+012.5f}"
    )


# ---------------- FLUX TEST ----------------
def run_test_harmonicity(window, label):
    scofo_desc = scofo.get_audio_description(window)
    l_flux = librosa.feature.spectral_flux(
        y=window,
        frame_length=n_fft,
        hop_length=hop,
        center=False,
    )[0, 0]
    pass


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
        f"RMS  | "
        f"L: {l_loudness:+012.5f} | "
        f"S: {s_loudness:+012.5f} | "
        f"D: {diff:+012.5f}"
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
        f"ZCR  | "
        f"L: {l_zcr:+012.5f} | "
        f"S: {s_zcr:+012.5f} | "
        f"D: {diff:+012.5f}"
    )


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
    run_test_spread(window, f"Start {start:08d}")
    run_test_centroid(window, f"Start {start:08d}")
    run_test_chroma(window, f"Start {start:08d}")

    print("")
