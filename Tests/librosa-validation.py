import os
import math
import random
import OpenScofo
import librosa

import essentia
import essentia.standard as es
import numpy as np

os.chdir(os.path.dirname(__file__))

sr = 48000
n_fft = 2048
hop = 512

scofo = OpenScofo.OpenScofo(sr, n_fft, hop)

y, sr = librosa.load(
    "./assets/bwv-1013.mp3",
    sr=sr,
)


# ---------------- Helper function ----------------
def format_diff(diff):
    if diff == 0:
        return "0.0e+00", "0"
    else:
        order = math.floor(math.log10(diff))
        return f"{diff:.5e}", f"10^{order}"


# ---------------- Track max difference globally ----------------
max_diffs = {
    "MFCC": 0.0,
    "CENT": 0.0,
    "SPRE": 0.0,
    "FLAT": 0.0,
    "CHRO": 0.0,
    "RMS": 0.0,
    "ZCR": 0.0,
    "Flux": 0.0,
}


def update_max(desc, diff):
    if diff > max_diffs[desc]:
        max_diffs[desc] = diff


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
    max_vals = (0.0, 0.0)

    for l, s in zip(librosa_mfcc, scofo_mfcc):
        abs_diff = abs(l - s)
        if abs_diff > max_diff:
            max_diff = abs_diff
            max_vals = (l, s)

    update_max("MFCC", max_diff)
    diff_str, order_str = format_diff(max_diff)

    print(
        f"{label} | MFCC  | L: {max_vals[0]:+012.5f} | S: {max_vals[1]:+012.5f} | D: {diff_str} | Order: {order_str}"
    )


# ---------------- Centroid TEST ----------------
def run_test_centroid(window, label):
    scofo_desc = scofo.get_audio_description(window)

    l_centroid = librosa.feature.spectral_centroid(
        y=window,
        sr=sr,
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
    )[0, 0]

    s_centroid = scofo_desc.spectral_centroid
    diff = abs(l_centroid - s_centroid)

    update_max("CENT", diff)
    diff_str, order_str = format_diff(diff)

    print(
        f"{label} | CENT  | L: {l_centroid:+012.5f} | S: {s_centroid:+012.5f} | D: {diff_str} | Order: {order_str}"
    )


# ---------------- Spread TEST ----------------
def run_test_spread(window, label):
    scofo_desc = scofo.get_audio_description(window)

    l_spread = librosa.feature.spectral_bandwidth(
        y=window,
        sr=sr,
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
        p=2,
    )[0, 0]

    s_spread = scofo_desc.spectral_spread
    diff = abs(l_spread - s_spread)

    update_max("SPRE", diff)
    diff_str, order_str = format_diff(diff)

    print(
        f"{label} | SPRE  | L: {l_spread:+012.5f} | S: {s_spread:+012.5f} | D: {diff_str} | Order: {order_str}"
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

    update_max("FLAT", diff)
    diff_str, order_str = format_diff(diff)

    print(
        f"{label} | FLAT  | L: {l_flatness:+012.5f} | S: {s_flatness:+012.5f} | D: {diff_str} | Order: {order_str}"
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

    for l_val, s_val in zip(librosa_chroma, scofo_chroma):
        abs_diff = abs(l_val - s_val)
        if abs_diff > max_diff:
            max_diff = abs_diff
            max_vals = (l_val, s_val)

    update_max("CHRO", max_diff)
    diff_str, order_str = format_diff(max_diff)

    print(
        f"{label} | CHRO  | L: {max_vals[0]:+012.5f} | S: {max_vals[1]:+012.5f} | D: {diff_str} | Order: {order_str}"
    )


# ---------------- RMS TEST ----------------
def run_test_rms(window, label):
    scofo_desc = scofo.get_audio_description(window)

    # librosa
    l_rms = librosa.feature.rms(
        y=window,
        frame_length=n_fft,
        hop_length=hop,
        center=False,
    )[0, 0]

    # scofo
    s_rms = scofo_desc.rms

    # essentia
    e_rms = es.RMS()(window)

    # differences
    diff_ls = abs(l_rms - s_rms)
    diff_le = abs(l_rms - e_rms)
    diff_se = abs(s_rms - e_rms)

    update_max("RMS", diff_ls)
    diff_str, order_str = format_diff(diff_ls)

    print(
        f"{label} | RMS | "
        f"L: {l_rms:+012.5f} | "
        f"S: {s_rms:+012.5f} | "
        f"E: {e_rms:+012.5f} | "
        f"D(L-S): {diff_str} | "
        f"D(L-E): {diff_le:.3e} | "
        f"D(S-E): {diff_se:.3e} | "
        f"Order: {order_str}"
    )


# ---------------- RUN LOUDNESS TEST ----------------
spectrum = es.Spectrum()
flux_alg = es.Flux()


def run_test_flux(prev_window, window, label):

    scofo_desc = scofo.get_audio_description(window)
    s_flux = scofo_desc.spectral_flux

    prev_spec = spectrum(prev_window)
    curr_spec = spectrum(window)

    # initialize internal state
    flux_alg(prev_spec)

    # compute flux between prev_spec and curr_spec
    e_flux = flux_alg(curr_spec)

    diff = abs(s_flux - e_flux)

    update_max("Flux", diff)
    diff_str, order_str = format_diff(diff)

    print(
        f"{label} | FLUX | "
        f"S: {s_flux:+012.5f} | "
        f"E: {e_flux:+012.5f} | "
        f"D(S-E): {diff_str} | "
        f"Order: {order_str}"
    )


# ---------------- ZEROCROSS-RATE TEST ----------------
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

    update_max("ZCR", diff)
    diff_str, order_str = format_diff(diff)

    print(
        f"{label} | ZCR   | L: {l_zcr:+012.5f} | S: {s_zcr:+012.5f} | D: {diff_str} | Order: {order_str}"
    )


# ---------------- RUN TESTS ----------------
n_tests = 100

for _ in range(n_tests):
    max_start = len(y) - 2 * n_fft
    start = random.randint(n_fft, max_start)

    prev_window = y[start - n_fft : start]
    window = y[start : start + n_fft]

    # run_test_mfcc(window, f"Start {start:08d}")
    # run_test_flatness(window, f"Start {start:08d}")
    # run_test_rms(window, f"Start {start:08d}")
    # run_test_zcr(window, f"Start {start:08d}")
    # run_test_spread(window, f"Start {start:08d}")
    # run_test_centroid(window, f"Start {start:08d}")
    # run_test_chroma(window, f"Start {start:08d}")
    run_test_flux(prev_window, window, f"Start {start:08d}")
    print("")


# ---------------- PRINT MAX DIFFERENCE ----------------
print("========== MAX DIFFERENCES (worst-case precision) ==========")
for desc, diff in max_diffs.items():
    _, order_str = format_diff(diff)
    print(f"{desc}: {diff:.5e} | Order: {order_str}")
