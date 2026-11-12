import random
import os
import json

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
scofo.activate_all_descriptors()

y, sr = librosa.load(
    # "./../01-benchmark/real/bwv-1013.wav",
    "./../04-miniaturas/miniatura2-larissa.wav",
    sr=sr,
)

MAX_DIFF = 0
MAX_DIFF_ALLOWED = 0.003
RANDOM_SEED = 0
N_TESTS = 20
RESULTS_FILE = "descriptor_results.json"

# Save the exact inputs so a C/C++ benchmark can reproduce each test without Python.
RESULTS = {
    "schema_version": 1,
    "configuration": {
        "audio_file": "../01-benchmark/real/bwv-1013.wav",
        "sample_rate": sr,
        "fft_size": n_fft,
        "reference_hop": hop,
        "openscofo_hop": analysis_hop,
        "threshold": MAX_DIFF_ALLOWED,
        "random_seed": RANDOM_SEED,
        "n_tests": N_TESTS,
        "essentia_spectrum_source": "openscofo.magnitude",
    },
    "cases": [],
    "summary": {},
}
CURRENT_CASE = None

RED = "\033[31m"
RESET = "\033[0m"

ess_window = es.Windowing(type="hann")
ess_spectrum = es.Spectrum()


# ─────────────────────────────────────
def json_array(values):
    array = np.asarray(values, dtype=np.float64).reshape(-1)
    return [float(v) if np.isfinite(v) else None for v in array]


# ─────────────────────────────────────
def record_result(name, library, reference, openscofo):
    """Save complete reference/OpenScofo values and their absolute differences."""
    global MAX_DIFF

    reference_array = np.asarray(reference, dtype=np.float64)
    openscofo_array = np.asarray(openscofo, dtype=np.float64)
    scalar = reference_array.ndim == 0 and openscofo_array.ndim == 0
    reference_array = reference_array.reshape(-1)
    openscofo_array = openscofo_array.reshape(-1)

    differences = []
    max_difference = None
    max_index = None
    for i in range(min(len(reference_array), len(openscofo_array))):
        a = float(reference_array[i])
        b = float(openscofo_array[i])
        difference = abs(a - b) if np.isfinite(a) and np.isfinite(b) else None
        differences.append(difference)
        if difference is not None and (
            max_difference is None or difference > max_difference
        ):
            max_difference = difference
            max_index = i

    valid = (
        len(reference_array) > 0
        and len(reference_array) == len(openscofo_array)
        and bool(np.isfinite(reference_array).all())
        and bool(np.isfinite(openscofo_array).all())
    )
    passed = bool(
        valid and max_difference is not None and max_difference <= MAX_DIFF_ALLOWED
    )

    # JSON does not support NaN/Infinity. Represent invalid measurements as null.
    reference_values = [float(v) if np.isfinite(v) else None for v in reference_array]
    openscofo_values = [float(v) if np.isfinite(v) else None for v in openscofo_array]
    CURRENT_CASE["descriptors"][name] = {
        "library": library,
        "reference": reference_values[0] if scalar else reference_values,
        "openscofo": openscofo_values[0] if scalar else openscofo_values,
        "absolute_difference": (
            differences[0] if scalar and differences else differences
        ),
        "max_absolute_difference": max_difference,
        "max_error_index": max_index if not scalar else None,
        "reference_length": len(reference_array),
        "openscofo_length": len(openscofo_array),
        "passed": passed,
    }

    # A length mismatch or non-finite value must fail, even if the shared prefix matches.
    MAX_DIFF = max(MAX_DIFF, max_difference if valid else MAX_DIFF_ALLOWED + 1.0)


# ─────────────────────────────────────
def save_results():
    """Save all cases, with a per-descriptor summary and strict, portable JSON."""
    summary = {}
    for case in RESULTS["cases"]:
        for name, record in case["descriptors"].items():
            if name not in summary:
                summary[name] = {
                    "library": record["library"],
                    "n_tests": 0,
                    "max_absolute_difference": None,
                    "passed": True,
                }
            item = summary[name]
            item["n_tests"] += 1
            item["passed"] = item["passed"] and record["passed"]
            value = record["max_absolute_difference"]
            if value is not None and (
                item["max_absolute_difference"] is None
                or value > item["max_absolute_difference"]
            ):
                item["max_absolute_difference"] = value

    RESULTS["summary"] = summary
    RESULTS["max_absolute_difference"] = MAX_DIFF
    RESULTS["passed"] = all(item["passed"] for item in summary.values()) and bool(
        summary
    )

    # Write atomically to avoid leaving an incomplete JSON on interruption.
    temporary_file = RESULTS_FILE + ".tmp"
    with open(temporary_file, "w", encoding="utf-8") as file:
        json.dump(RESULTS, file, indent=2, sort_keys=True, allow_nan=False)
        file.write("\n")
    os.replace(temporary_file, RESULTS_FILE)
    print(f"Saved {RESULTS_FILE} ({len(RESULTS['cases'])} test cases)")


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

    record_result("logmel", "librosa", librosa_logmel, scofo_logmel)
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
    record_result("mfcc", "librosa", librosa_mfcc, scofo_mfcc)

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
    record_result("kurtosis_scipy", "scipy", k, s_kurtosis)
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
    record_result("centroid", "librosa", l_centroid, s_centroid)
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
    record_result("spreadhz", "librosa", l_spread, s_spread)
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
    record_result("flatness", "librosa", l_flatness, s_flatness)
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
    record_result("crest", "essentia", essentia_crest, o_crest)

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
    record_result("flux", "essentia", essentia_flux, o_flux)

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
    record_result("rolloff", "essentia", e_rolloff, o_rolloff)

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
    record_result("entropy", "essentia", e_entropy, o_entropy)

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
    record_result("spread_v", "essentia", e_spread, o_spread)
    record_result("skewness", "essentia", e_skew, o_skew)
    record_result("kurtosis", "essentia", e_kurt, o_kurt)

    # Compute diffs
    spread_diff = abs(e_spread - o_spread)
    skew_diff = abs(e_skew - o_skew)
    kurt_diff = abs(e_kurt - o_kurt)

    MAX_DIFF = max(MAX_DIFF, spread_diff, skew_diff, kurt_diff)

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
    record_result("chroma", "librosa", librosa_chroma, scofo_chroma)

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
    record_result("rms", "librosa", l_loudness, s_loudness)
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
    record_result("zcr", "librosa", l_zcr, s_zcr)
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
random.seed(RANDOM_SEED)

for _ in range(N_TESTS):
    max_start = len(y) - n_fft
    start = random.randint(n_fft, max_start)
    prev_window = y[start - hop : start - hop + n_fft]
    window = y[start : start + n_fft]

    # Inputs are saved as float32-derived JSON numbers, matching process_block().
    CURRENT_CASE = {
        "start_sample": start,
        "input": {
            "previous": json_array(prev_window),
            "current": json_array(window),
        },
        "descriptors": {},
    }
    RESULTS["cases"].append(CURRENT_CASE)

    prev_desc = analyze_window(prev_window)
    curr_desc = analyze_window(window)

    # Essentia's spectrum-based comparisons use these exact OpenScofo spectra.
    CURRENT_CASE["input"]["previous_magnitude"] = json_array(prev_desc.magnitude)
    CURRENT_CASE["input"]["current_magnitude"] = json_array(curr_desc.magnitude)

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


save_results()  # Save even when the tolerance test fails.

if MAX_DIFF > MAX_DIFF_ALLOWED:
    raise Exception(f"MAX_DIFF {MAX_DIFF} | MAX_DIFF_ALLOWED {MAX_DIFF_ALLOWED}")


print(f"MAX_DIFF {MAX_DIFF} | MAX_DIFF_ALLOWED {MAX_DIFF_ALLOWED}")
