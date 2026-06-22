import random
import os

import OpenScofo
import librosa
import essentia

essentia.log.infoActive = False
import essentia.standard as es

import numpy as np
from scipy.stats import kurtosis

sr = 48000
n_fft = 2048
hop = 512
analysis_hop = n_fft

os.chdir(os.path.dirname(__file__))

# Use one analysis per window so flux compares consecutive frames.
scofo = OpenScofo.OpenScofo(sr, n_fft, analysis_hop)
scofo.set_requested_descriptors(
    [
        OpenScofo.Descriptors.LOGMEL,
        OpenScofo.Descriptors.MFCC,
        OpenScofo.Descriptors.RMS,
        OpenScofo.Descriptors.ZCR,
        OpenScofo.Descriptors.SPREAD,
        OpenScofo.Descriptors.CENTROID,
        OpenScofo.Descriptors.CHROMA,
        OpenScofo.Descriptors.ROLLOFF,
        OpenScofo.Descriptors.ENTROPY,
        OpenScofo.Descriptors.CREST,
        OpenScofo.Descriptors.FLUX,
    ]
)

y, sr = librosa.load(
    "./assets/bwv-1013.wav",
    sr=sr,
)

MAX_DIFF = 0
MAX_DIFF_ALLOWED = 0.003

RED = "\033[31m"
RESET = "\033[0m"

ess_window = es.Windowing(type="hann")
ess_spectrum = es.Spectrum()


def essentia_mag_spectrum(frame):
    frame32 = np.ascontiguousarray(frame, dtype=np.float32)
    return ess_spectrum(ess_window(frame32))


def print_with_threshold(diff, message):
    if diff > MAX_DIFF_ALLOWED:
        print(f"{RED}{message}{RESET}")
    else:
        print(message)


def analyze_window(window):
    scofo.process_block(window)
    return scofo.get_description()


def run_test_logmel(window, scofo_desc, label):
    global MAX_DIFF

    # Librosa log-mel
    mel = librosa.feature.melspectrogram(
        y=window,
        sr=sr,
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
        n_mels=40,
        power=2.0,
    )

    librosa_logmel = librosa.power_to_db(mel, ref=1.0)[:, 0]
    scofo_logmel = scofo_desc.logmel

    n = min(len(librosa_logmel), len(scofo_logmel))
    if n == 0:
        print_with_threshold(
            MAX_DIFF_ALLOWED + 1.0,
            "Librosa  "
            f"{label} | "
            f"LOGM | Empty feature vector (librosa={len(librosa_logmel)}, scofo={len(scofo_logmel)})",
        )
        MAX_DIFF = max(MAX_DIFF, MAX_DIFF_ALLOWED + 1.0)
        return

    max_diff = 0.0
    max_idx = 0

    for i, (l, s) in enumerate(zip(librosa_logmel[:n], scofo_logmel[:n])):
        diff = abs(l - s)
        if diff > max_diff:
            max_diff = diff
            max_idx = i

    MAX_DIFF = max(MAX_DIFF, max_diff)

    print_with_threshold(
        max_diff,
        "Librosa  "
        f"{label} | "
        f"LOGM | "
        f"L: {librosa_logmel[max_idx]:+012.5f} | "
        f"S: {scofo_logmel[max_idx]:+012.5f} | "
        f"D: {max_diff:+012.5f} | "
        f"N: {n:03d}",
    )


# ---------------- MFCC TEST ----------------
def run_test_mfcc(window, scofo_desc, label):
    global MAX_DIFF

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

    for i, (l, s) in enumerate(zip(librosa_mfcc, scofo_mfcc)):
        abs_diff = abs(l - s)
        if abs_diff > max_diff:
            max_diff = abs_diff
            max_vals = (l, s)

    MAX_DIFF = max(abs(max_diff), MAX_DIFF)

    print_with_threshold(
        max_diff,
        "Librosa  "
        f"{label} | "
        f"MFCC | "
        f"L: {max_vals[0]:+012.5f} | "
        f"S: {max_vals[1]:+012.5f} | "
        f"D: {max_diff:+012.5f}",
    )


# ──────────────────────────────────────
def run_test_kurtosis(window, scofo_desc, label):
    global MAX_DIFF
    S = np.abs(librosa.stft(window))
    k = kurtosis(S, fisher=True, bias=True)

    s_kurtosis = scofo_desc.kurtosis
    diff = abs(k - s_kurtosis)
    MAX_DIFF = max(abs(diff), MAX_DIFF)

    print_with_threshold(
        diff,
        "Librosa  "
        f"{label} | "
        f"CENT | "
        f"L: {k:+012.5f} | "
        f"S: {s_kurtosis:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- Centroid TEST ----------------
def run_test_centroid(window, scofo_desc, label):
    global MAX_DIFF

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
    s_centroid = scofo_desc.centroid
    diff = abs(l_centroid - s_centroid)
    MAX_DIFF = max(abs(diff), MAX_DIFF)

    print_with_threshold(
        diff,
        "Librosa  "
        f"{label} | "
        f"CENT | "
        f"L: {l_centroid:+012.5f} | "
        f"S: {s_centroid:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- Spread TEST ----------------
def run_test_spread(window, scofo_desc, label):
    global MAX_DIFF

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
    s_spread = scofo_desc.spreadhz
    diff = abs(l_spread - s_spread)
    MAX_DIFF = max(abs(diff), MAX_DIFF)

    print_with_threshold(
        diff,
        "Librosa  "
        f"{label} | "
        f"SPRE | "
        f"L: {l_spread:+012.5f} | "
        f"S: {s_spread:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- FLATNESS TEST ----------------
def run_test_flatness(window, scofo_desc, label):
    global MAX_DIFF

    l_flatness = librosa.feature.spectral_flatness(
        y=window,
        n_fft=n_fft,
        hop_length=hop,
        win_length=n_fft,
        window="hann",
        center=False,
        power=2.0,
    )[0, 0]

    s_flatness = scofo_desc.flatness
    diff = abs(l_flatness - s_flatness)
    MAX_DIFF = max(abs(diff), MAX_DIFF)

    print_with_threshold(
        diff,
        "Librosa  "
        f"{label} | "
        f"FLAT | "
        f"L: {l_flatness:+012.5f} | "
        f"S: {s_flatness:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- CREST TEST ----------------
def run_test_crest(scofo_desc, label):
    global MAX_DIFF

    o_crest = scofo_desc.crest

    crest = es.Crest()

    # Compare algorithm only (same spectrum values for both sides).
    spectrum = np.ascontiguousarray(scofo_desc.magnitude, dtype=np.float32)
    essentia_crest = crest(spectrum)

    diff = abs(essentia_crest - o_crest)
    MAX_DIFF = max(diff, MAX_DIFF)

    print_with_threshold(
        diff,
        "Essentia "
        f"{label} | "
        f"CRES | "
        f"E: {essentia_crest:+012.5f} | "
        f"S: {o_crest:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- FLUX TEST ----------------
def run_test_flux(prev_desc, curr_desc, label):
    global MAX_DIFF
    o_flux = curr_desc.flux

    # Essentia (Flux memory behavior: first call seeds memory, second call computes prev->curr)
    flux_algo = es.Flux(norm="L2", halfRectify=False)

    # 1. Use OpenScofo spectra so the comparison is algorithm-only.
    spec_prev = np.ascontiguousarray(prev_desc.magnitude, dtype=np.float32)
    spec_curr = np.ascontiguousarray(curr_desc.magnitude, dtype=np.float32)

    # 2. Calculate Flux
    _ = flux_algo(spec_prev)
    essentia_flux = flux_algo(spec_curr)

    diff = abs(essentia_flux - o_flux)
    MAX_DIFF = max(diff, MAX_DIFF)

    print_with_threshold(
        diff,
        "Essentia "
        f"{label} | "
        f"FLUX | "
        f"E: {essentia_flux:+012.5f} | "
        f"S: {o_flux:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- ROLLOFF TEST ----------------
def run_test_rolloff(scofo_desc, label):
    global MAX_DIFF
    o_rolloff = scofo_desc.rolloff

    # Essentia rolloff (uses magnitude spectrum)
    rolloff_algo = es.RollOff(cutoff=0.85)  # match typical definition

    spectrum = np.ascontiguousarray(scofo_desc.magnitude, dtype=np.float32)
    e_rolloff = rolloff_algo(spectrum)

    diff = abs(e_rolloff - o_rolloff)
    MAX_DIFF = max(diff, MAX_DIFF)

    print_with_threshold(
        diff,
        "Essentia "
        f"{label} | "
        f"ROLL | "
        f"E: {e_rolloff:+012.5f} | "
        f"S: {o_rolloff:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- ENTROPY TEST ----------------
def run_test_entropy(scofo_desc, label):
    global MAX_DIFF
    o_entropy = scofo_desc.entropy

    # Essentia spectral entropy
    entropy_algo = es.Entropy()

    spectrum = np.ascontiguousarray(scofo_desc.magnitude, dtype=np.float32)

    # Normalize to probability distribution (required for entropy)
    s = spectrum.sum()
    if s > 0:
        spectrum = spectrum / s

    e_entropy = entropy_algo(spectrum)

    diff = abs(e_entropy - o_entropy)
    MAX_DIFF = max(diff, MAX_DIFF)

    print_with_threshold(
        diff,
        "Essentia "
        f"{label} | "
        f"ENTR | "
        f"E: {e_entropy:+012.5f} | "
        f"S: {o_entropy:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- FLUX TEST ----------------
def run_test_spread_skew_kurt(scofo_desc, label):
    global MAX_DIFF

    # OpenScofo targets
    o_spread = scofo_desc.spread_variance
    o_skew = scofo_desc.skewness
    o_kurt = scofo_desc.kurtosis

    # Essentia setup
    central_moments = es.CentralMoments()
    dist_shape = es.DistributionShape()

    # Use OpenScofo magnitude spectrum
    spec = np.ascontiguousarray(scofo_desc.magnitude, dtype=np.float32)

    # Normalize to probability distribution
    s = spec.sum()
    if s > 0:
        spec = spec / s

    # Compute moments → shape
    moments = central_moments(spec)
    e_spread, e_skew, e_kurt = dist_shape(moments)

    # Compute diffs
    spread_diff = abs(e_spread - o_spread)
    skew_diff = abs(e_skew - o_skew)
    kurt_diff = abs(e_kurt - o_kurt)

    MAX_DIFF = max(MAX_DIFF, spread_diff, skew_diff, kurt_diff)

    print_with_threshold(
        spread_diff,
        "Essentia "
        f"{label} | "
        f"SPRE | "
        f"E: {e_spread:+012.5f} | "
        f"S: {o_spread:+012.5f} | "
        f"D: {spread_diff:+012.5f}",
    )

    print_with_threshold(
        skew_diff,
        "Essentia "
        f"{label} | "
        f"SKEW | "
        f"E: {e_skew:+012.5f} | "
        f"S: {o_skew:+012.5f} | "
        f"D: {skew_diff:+012.5f}",
    )

    print_with_threshold(
        kurt_diff,
        "Essentia "
        f"{label} | "
        f"KURT | "
        f"E: {e_kurt:+012.5f} | "
        f"S: {o_kurt:+012.5f} | "
        f"D: {kurt_diff:+012.5f}",
    )


# ---------------- CHROMA TEST ----------------
def run_test_chroma(window, scofo_desc, label):
    global MAX_DIFF

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

    MAX_DIFF = max(abs(max_diff), MAX_DIFF)

    print_with_threshold(
        max_diff,
        "Librosa  "
        f"{label} | "
        f"CHRO | "
        f"L: {max_vals[0]:+012.5f} | "
        f"S: {max_vals[1]:+012.5f} | "
        f"D: {max_diff:+012.5f}",
    )


# ---------------- RMS TEST ----------------
def run_test_rms(window, scofo_desc, label):
    global MAX_DIFF

    # RMS-based loudness (power=2.0 equivalent)
    l_loudness = librosa.feature.rms(
        y=window,
        frame_length=n_fft,
        hop_length=hop,
        center=False,
    )[0, 0]

    s_loudness = scofo_desc.rms
    diff = abs(l_loudness - s_loudness)
    MAX_DIFF = max(abs(diff), MAX_DIFF)

    print_with_threshold(
        diff,
        "Librosa  "
        f"{label} | "
        f"RMS  | "
        f"L: {l_loudness:+012.5f} | "
        f"S: {s_loudness:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- ZEROCROSS-RATING TEST --------
def run_test_zcr(window, scofo_desc, label):
    global MAX_DIFF

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

    s_zcr = scofo_desc.zcr
    diff = abs(l_zcr - s_zcr)
    MAX_DIFF = max(abs(diff), MAX_DIFF)

    print_with_threshold(
        diff,
        "Librosa  "
        f"{label} | "
        f"ZCR  | "
        f"L: {l_zcr:+012.5f} | "
        f"S: {s_zcr:+012.5f} | "
        f"D: {diff:+012.5f}",
    )


# ---------------- TESTES ALEATÓRIOS ----------------
n_tests = 20

for _ in range(n_tests):
    max_start = len(y) - n_fft
    start = random.randint(n_fft, max_start)
    prev_window = y[start - hop : start - hop + n_fft]
    window = y[start : start + n_fft]

    prev_desc = analyze_window(prev_window)
    curr_desc = analyze_window(window)

    run_test_logmel(window, curr_desc, f"Start {start:08d}")
    run_test_mfcc(window, curr_desc, f"Start {start:08d}")
    run_test_flatness(window, curr_desc, f"Start {start:08d}")
    run_test_rms(window, curr_desc, f"Start {start:08d}")
    run_test_zcr(window, curr_desc, f"Start {start:08d}")
    run_test_spread(window, curr_desc, f"Start {start:08d}")
    run_test_centroid(window, curr_desc, f"Start {start:08d}")
    run_test_chroma(window, curr_desc, f"start {start:08d}")

    # Essentia
    # run_test_rolloff(curr_desc, f"Start {start:08d}")
    # run_test_entropy(curr_desc, f"Start {start:08d}")
    run_test_crest(curr_desc, f"start {start:08d}")
    run_test_flux(prev_desc, curr_desc, f"Start {start:08d}")
    run_test_spread_skew_kurt(curr_desc, f"Start {start:08d}")

    print("----")


if MAX_DIFF > MAX_DIFF_ALLOWED:
    raise Exception(f"MAX_DIFF {MAX_DIFF} | MAX_DIFF_ALLOWED {MAX_DIFF_ALLOWED}")


print(f"MAX_DIFF {MAX_DIFF} | MAX_DIFF_ALLOWED {MAX_DIFF_ALLOWED}")
