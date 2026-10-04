#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import resample_poly

# ============================================================
# Configuration
# ============================================================

EPS = 1e-24
SPECTRUM_EPS = 1e-12

MIN_F0 = 32.70
MAX_F0 = 4186.0

MIN_HARMONIC_DECAY = 0.2
MAX_HARMONIC_DECAY = 1.8

INHARMONICITY_B = 0.0001


# ============================================================
# Pitch utilities
# ============================================================

PITCH_RE = re.compile(r"(?<![A-Za-z0-9])([A-Ga-g])([#b]?)(-?\d+)(?=-|_|\.|$)")

NOTE_TO_SEMITONE = {
    "C": 0,
    "D": 2,
    "E": 4,
    "F": 5,
    "G": 7,
    "A": 9,
    "B": 11,
}


def extract_pitch_from_filename(path: Path) -> str | None:
    """
    Example:
        Acc-ord-A1-mf-alt1-N.wav -> A1
        Vln-ord-C#4-mf.wav       -> C#4
        Vc-ord-Bb2-p.wav         -> Bb2
    """
    match = PITCH_RE.search(path.stem)

    if match is None:
        return None

    note = match.group(1).upper()
    accidental = match.group(2)
    octave = match.group(3)

    return f"{note}{accidental}{octave}"


def pitch_to_midi(pitch: str) -> float:
    match = re.fullmatch(r"([A-G])([#b]?)(-?\d+)", pitch)

    if match is None:
        raise ValueError(f"Invalid pitch: {pitch}")

    note = match.group(1)
    accidental = match.group(2)
    octave = int(match.group(3))

    semitone = NOTE_TO_SEMITONE[note]

    if accidental == "#":
        semitone += 1
    elif accidental == "b":
        semitone -= 1

    return 12 * (octave + 1) + semitone


def midi_to_freq(midi: float, tuning: float = 440.0) -> float:
    return tuning * (2.0 ** ((midi - 69.0) / 12.0))


def pitch_to_freq(pitch: str, tuning: float = 440.0) -> float:
    return midi_to_freq(pitch_to_midi(pitch), tuning)


def shift_frequency_cents(freq: float, cents: float) -> float:
    return freq * (2.0 ** (cents / 1200.0))


# ============================================================
# Audio
# ============================================================


def pcm_to_float(data: np.ndarray) -> np.ndarray:
    if np.issubdtype(data.dtype, np.floating):
        return data.astype(np.float64)

    if np.issubdtype(data.dtype, np.signedinteger):
        bits = np.iinfo(data.dtype).bits
        scale = float(2 ** (bits - 1))
        return data.astype(np.float64) / scale

    if np.issubdtype(data.dtype, np.unsignedinteger):
        info = np.iinfo(data.dtype)
        midpoint = (info.max + 1) / 2.0
        return (data.astype(np.float64) - midpoint) / midpoint

    raise TypeError(f"Unsupported WAV dtype: {data.dtype}")


def load_audio(path: Path, target_sr: int) -> np.ndarray:
    sr, data = wavfile.read(path)

    data = pcm_to_float(data)

    if data.ndim == 2:
        data = np.mean(data, axis=1)

    if sr != target_sr:
        gcd = math.gcd(sr, target_sr)

        up = target_sr // gcd
        down = sr // gcd

        data = resample_poly(data, up, down)

    return np.asarray(data, dtype=np.float64)


# ============================================================
# OpenScofo spectrum
# ============================================================


def periodic_hann(size: int) -> np.ndarray:
    """
    Matches OpenScofo:

        0.5 * (1 - cos(2*pi*i/N))

    This is the periodic Hann, not np.hanning(N).
    """
    n = np.arange(size, dtype=np.float64)

    return 0.5 * (1.0 - np.cos(2.0 * np.pi * n / float(size)))


@dataclass
class FrameSpectrum:
    spectrum: np.ndarray
    stddev: float
    rms: float


def get_spectrum(
    frame: np.ndarray,
    fft_size: int,
    window: np.ndarray,
) -> FrameSpectrum:

    windowed = frame * window

    fft = np.fft.rfft(windowed, n=fft_size)

    magnitude = np.abs(fft)

    # OpenScofo:
    #
    # SpectralMagnitudeNorm[i] = mag / FFTSize
    #
    magnitude_norm = magnitude / float(fft_size)

    sum_power = np.sum(magnitude_norm)

    # Match:
    #
    # normSp =
    #   (SpectralMagnitudeNorm[i] + 1e-12)
    #   / (SumPower + 1e-12)
    #
    spectrum = (magnitude_norm + SPECTRUM_EPS) / (sum_power + SPECTRUM_EPS)

    n_half = spectrum.size

    mean = 1.0 / float(n_half)

    variance = np.sum((spectrum - mean) ** 2)

    # OpenScofo:
    #
    # Desc.StdDev = sqrt(Variance * Mean)
    #
    stddev = math.sqrt(variance * mean)

    rms = math.sqrt(np.mean(frame * frame) + 1e-30)

    return FrameSpectrum(
        spectrum=spectrum,
        stddev=stddev,
        rms=rms,
    )


# ============================================================
# Shared harmonic envelope
# ============================================================


def harmonic_beta(
    freq: float,
    amplitude_decay: float,
) -> float:

    f0_norm = math.log2(freq / MIN_F0) / math.log2(MAX_F0 / MIN_F0)

    f0_norm = max(
        0.0,
        min(1.0, f0_norm),
    )

    shaped = (1.0 - f0_norm) ** amplitude_decay

    return MAX_HARMONIC_DECAY - shaped * (MAX_HARMONIC_DECAY - MIN_HARMONIC_DECAY)


def harmonic_envelope(
    harmonic: int,
    beta: float,
) -> float:

    envelope = math.exp(-beta * (harmonic - 1))

    # Current OpenScofo behavior
    if harmonic > 1:
        envelope *= 1.25

    return envelope


# ============================================================
# CURRENT OpenScofo template
# ============================================================


def build_current_template(
    freq: float,
    sr: int,
    fft_size: int,
    harmonics: int,
    pitch_template_sigma: float,
    amplitude_decay: float,
) -> np.ndarray:
    """
    Reproduces the current OpenScofo BuildPitchTemplate()
    as closely as possible.

    Important current behavior:

        sigmaLog = PitchTemplateSigma / 12
        sigmaConst = 2^sigmaLog - 1

        sigmaHz = harmonicFreq * sigmaConst

        sigmaHz >= binWidth * 0.75
    """

    bins_count = fft_size // 2

    template = np.full(
        bins_count,
        EPS,
        dtype=np.float64,
    )

    bin_width = sr / float(fft_size)

    beta = harmonic_beta(
        freq,
        amplitude_decay,
    )

    sigma_log = pitch_template_sigma / 12.0

    sigma_const = (2.0**sigma_log) - 1.0

    for k in range(1, harmonics + 1):

        stretch = math.sqrt(1.0 + INHARMONICITY_B * k * k)

        harmonic_freq = freq * k * stretch

        if harmonic_freq >= sr / 2.0:
            break

        sigma_hz = harmonic_freq * sigma_const

        # Important limitation in current implementation
        sigma_hz = max(
            sigma_hz,
            bin_width * 0.75,
        )

        envelope = harmonic_envelope(
            k,
            beta,
        )

        if envelope < 1e-5:
            break

        range_hz = 4.0 * sigma_hz

        min_bin = int(math.floor((harmonic_freq - range_hz) / bin_width))

        max_bin = int(math.ceil((harmonic_freq + range_hz) / bin_width))

        min_bin = max(
            0,
            min_bin,
        )

        max_bin = min(
            bins_count - 1,
            max_bin,
        )

        two_sigma_sq = 2.0 * sigma_hz * sigma_hz

        normalization_factor = 1.0 / (sigma_hz * math.sqrt(2.0 * math.pi))

        indices = np.arange(
            min_bin,
            max_bin + 1,
        )

        bin_freqs = indices * bin_width

        diff = bin_freqs - harmonic_freq

        gaussian = normalization_factor * np.exp(-(diff * diff) / two_sigma_sq)

        template[indices] += envelope * gaussian

    total = np.sum(template)

    if total > 1e-12:
        template /= total

    return template


# ============================================================
# PROPOSED high-resolution template
# ============================================================


def build_proposed_template(
    freq: float,
    sr: int,
    fft_size: int,
    harmonics: int,
    sigma_cents: float,
    amplitude_decay: float,
) -> np.ndarray:
    """
    Proposed template.

    Differences from current OpenScofo:

    1. No binWidth * 0.75 minimum sigma.
    2. Width is expressed directly in cents.
    3. Harmonic Gaussian is evaluated using cents distance.
    4. Each harmonic Gaussian is normalized BEFORE applying
       its harmonic envelope.

    Everything still ends as a full FFT-bin probability
    distribution and is compared with the SAME KL divergence.
    """

    bins_count = fft_size // 2

    template = np.full(
        bins_count,
        EPS,
        dtype=np.float64,
    )

    bin_width = sr / float(fft_size)

    frequencies = (
        np.arange(
            bins_count,
            dtype=np.float64,
        )
        * bin_width
    )

    beta = harmonic_beta(
        freq,
        amplitude_decay,
    )

    for k in range(1, harmonics + 1):

        stretch = math.sqrt(1.0 + INHARMONICITY_B * k * k)

        harmonic_freq = freq * k * stretch

        if harmonic_freq >= sr / 2.0:
            break

        envelope = harmonic_envelope(
            k,
            beta,
        )

        if envelope < 1e-5:
            break

        # +/- 4 sigma in cents
        lower_freq = harmonic_freq * 2.0 ** (-4.0 * sigma_cents / 1200.0)

        upper_freq = harmonic_freq * 2.0 ** (4.0 * sigma_cents / 1200.0)

        min_bin = max(
            1,
            int(math.floor(lower_freq / bin_width)),
        )

        max_bin = min(
            bins_count - 1,
            int(math.ceil(upper_freq / bin_width)),
        )

        #
        # For very low fundamentals, a narrow cents window
        # can fall entirely between FFT bins.
        #
        # In that case include the nearest bin. Higher
        # harmonics will provide most of the useful resolution.
        #
        nearest_bin = int(round(harmonic_freq / bin_width))

        nearest_bin = max(
            1,
            min(
                bins_count - 1,
                nearest_bin,
            ),
        )

        min_bin = min(
            min_bin,
            nearest_bin,
        )

        max_bin = max(
            max_bin,
            nearest_bin,
        )

        indices = np.arange(
            min_bin,
            max_bin + 1,
        )

        bin_freqs = frequencies[indices]

        valid = bin_freqs > 0.0

        cents = np.zeros_like(bin_freqs)

        cents[valid] = 1200.0 * np.log2(bin_freqs[valid] / harmonic_freq)

        gaussian = np.exp(-0.5 * (cents / sigma_cents) ** 2)

        #
        # Normalize each harmonic independently.
        #
        # This prevents high harmonics from receiving
        # more total template mass merely because the
        # same cents width covers more Hz / FFT bins.
        #
        gaussian_sum = np.sum(gaussian)

        if gaussian_sum <= 1e-30:
            gaussian = np.zeros_like(gaussian)

            local_index = int(np.argmin(np.abs(bin_freqs - harmonic_freq)))

            gaussian[local_index] = 1.0

        else:
            gaussian /= gaussian_sum

        template[indices] += envelope * gaussian

    total = np.sum(template)

    if total > 1e-12:
        template /= total

    return template


# ============================================================
# KL probability
# ============================================================


def kl_probability(
    template: np.ndarray,
    spectrum: np.ndarray,
    stddev: float,
    pitch_scaling_factor: float,
) -> float:
    """
    Matches current OpenScofo direction:

        P = template
        Q = observed spectrum

        KL = sum(P * log(P / Q))

        KL /= (1 + StdDev)

        probability =
            exp(-PitchScalingFactor * KL)
    """

    bins = min(
        len(template),
        len(spectrum),
    )

    p = template[:bins]
    q = spectrum[:bins]

    q = np.maximum(
        q,
        1e-300,
    )

    mask = p > 0.0

    kl = np.sum(p[mask] * np.log(p[mask] / q[mask]))

    kl /= 1.0 + stddev

    probability = math.exp(-pitch_scaling_factor * kl)

    # Both vectors should produce KL >= 0.
    # Clamp tiny numerical excursions.
    return float(
        np.clip(
            probability,
            0.0,
            1.0,
        )
    )


# ============================================================
# Frame analysis
# ============================================================


@dataclass
class MethodResult:
    median: float
    mean: float
    maximum: float


@dataclass
class FileResult:
    path: Path
    pitch: str
    frequency: float

    current: MethodResult
    proposed: MethodResult

    current_minus_50: float
    current_plus_50: float

    proposed_minus_50: float
    proposed_plus_50: float

    current_margin_50: float
    proposed_margin_50: float

    current_best_offset: int
    proposed_best_offset: int

    frames: int


def aggregate(values: list[float]) -> MethodResult:
    x = np.asarray(
        values,
        dtype=np.float64,
    )

    return MethodResult(
        median=float(np.median(x)),
        mean=float(np.mean(x)),
        maximum=float(np.max(x)),
    )


def analyze_file(
    path: Path,
    pitch: str,
    sr: int,
    fft_size: int,
    hop_size: int,
    harmonics: int,
    current_sigma: float,
    proposal_sigma_cents: float,
    amplitude_decay: float,
    pitch_scaling_factor: float,
    skip_attack_ms: float,
    rms_threshold_db: float,
) -> FileResult | None:

    audio = load_audio(
        path,
        sr,
    )

    if len(audio) < fft_size:
        return None

    freq = pitch_to_freq(pitch)

    window = periodic_hann(fft_size)

    offsets = (
        -100,
        -50,
        0,
        50,
        100,
    )

    current_templates = {}
    proposed_templates = {}

    for cents in offsets:

        candidate_freq = shift_frequency_cents(
            freq,
            cents,
        )

        current_templates[cents] = build_current_template(
            candidate_freq,
            sr,
            fft_size,
            harmonics,
            current_sigma,
            amplitude_decay,
        )

        proposed_templates[cents] = build_proposed_template(
            candidate_freq,
            sr,
            fft_size,
            harmonics,
            proposal_sigma_cents,
            amplitude_decay,
        )

    frames: list[FrameSpectrum] = []

    for start in range(
        0,
        len(audio) - fft_size + 1,
        hop_size,
    ):
        frame = audio[start : start + fft_size]

        frames.append(
            get_spectrum(
                frame,
                fft_size,
                window,
            )
        )

    if not frames:
        return None

    rms_values = np.asarray([x.rms for x in frames])

    max_rms = float(np.max(rms_values))

    if max_rms <= 1e-15:
        return None

    rms_threshold = max_rms * 10.0 ** (rms_threshold_db / 20.0)

    skip_frames = int(round((skip_attack_ms / 1000.0) * sr / hop_size))

    active_indices = np.where(rms_values >= rms_threshold)[0]

    if len(active_indices) == 0:
        return None

    first_active = int(active_indices[0])

    minimum_frame = first_active + skip_frames

    selected = [frames[i] for i in active_indices if i >= minimum_frame]

    # Very short sample: use active frames rather than fail.
    if not selected:
        selected = [frames[i] for i in active_indices]

    current_scores: dict[
        int,
        list[float],
    ] = {offset: [] for offset in offsets}

    proposed_scores: dict[
        int,
        list[float],
    ] = {offset: [] for offset in offsets}

    for frame in selected:

        for offset in offsets:

            current_scores[offset].append(
                kl_probability(
                    current_templates[offset],
                    frame.spectrum,
                    frame.stddev,
                    pitch_scaling_factor,
                )
            )

            proposed_scores[offset].append(
                kl_probability(
                    proposed_templates[offset],
                    frame.spectrum,
                    frame.stddev,
                    pitch_scaling_factor,
                )
            )

    current_medians = {
        offset: float(np.median(values)) for offset, values in current_scores.items()
    }

    proposed_medians = {
        offset: float(np.median(values)) for offset, values in proposed_scores.items()
    }

    current_best_offset = max(
        offsets,
        key=lambda x: current_medians[x],
    )

    proposed_best_offset = max(
        offsets,
        key=lambda x: proposed_medians[x],
    )

    current_neighbor = max(
        current_medians[-50],
        current_medians[50],
    )

    proposed_neighbor = max(
        proposed_medians[-50],
        proposed_medians[50],
    )

    return FileResult(
        path=path,
        pitch=pitch,
        frequency=freq,
        current=aggregate(current_scores[0]),
        proposed=aggregate(proposed_scores[0]),
        current_minus_50=(current_medians[-50]),
        current_plus_50=(current_medians[50]),
        proposed_minus_50=(proposed_medians[-50]),
        proposed_plus_50=(proposed_medians[50]),
        current_margin_50=(current_medians[0] - current_neighbor),
        proposed_margin_50=(proposed_medians[0] - proposed_neighbor),
        current_best_offset=(current_best_offset),
        proposed_best_offset=(proposed_best_offset),
        frames=len(selected),
    )


# ============================================================
# Dataset
# ============================================================


def find_ordinario_files(
    root: Path,
) -> list[Path]:

    files = []

    for path in root.rglob("*.wav"):

        parents = {parent.name.lower() for parent in path.parents}

        if "ordinario" in parents:
            files.append(path)

    return sorted(files)


# ============================================================
# CSV
# ============================================================


def write_csv(
    results: list[FileResult],
    output: Path,
) -> None:

    with output.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as file:

        writer = csv.writer(file)

        writer.writerow(
            [
                "file",
                "pitch",
                "frequency_hz",
                "frames",
                "current_median",
                "current_mean",
                "current_max",
                "proposed_median",
                "proposed_mean",
                "proposed_max",
                "current_minus_50",
                "current_plus_50",
                "proposed_minus_50",
                "proposed_plus_50",
                "current_margin_50",
                "proposed_margin_50",
                "current_best_offset_cents",
                "proposed_best_offset_cents",
            ]
        )

        for r in results:

            writer.writerow(
                [
                    str(r.path),
                    r.pitch,
                    f"{r.frequency:.6f}",
                    r.frames,
                    f"{r.current.median:.8f}",
                    f"{r.current.mean:.8f}",
                    f"{r.current.maximum:.8f}",
                    f"{r.proposed.median:.8f}",
                    f"{r.proposed.mean:.8f}",
                    f"{r.proposed.maximum:.8f}",
                    f"{r.current_minus_50:.8f}",
                    f"{r.current_plus_50:.8f}",
                    f"{r.proposed_minus_50:.8f}",
                    f"{r.proposed_plus_50:.8f}",
                    f"{r.current_margin_50:.8f}",
                    f"{r.proposed_margin_50:.8f}",
                    r.current_best_offset,
                    r.proposed_best_offset,
                ]
            )


# ============================================================
# Summary
# ============================================================


def print_summary(
    results: list[FileResult],
) -> None:

    if not results:
        print("No valid files analyzed.")
        return

    current_prob = np.asarray([x.current.median for x in results])

    proposed_prob = np.asarray([x.proposed.median for x in results])

    current_margin = np.asarray([x.current_margin_50 for x in results])

    proposed_margin = np.asarray([x.proposed_margin_50 for x in results])

    current_correct = np.asarray([x.current_best_offset == 0 for x in results])

    proposed_correct = np.asarray([x.proposed_best_offset == 0 for x in results])

    print()
    print("=" * 72)
    print("SUMMARY")
    print("=" * 72)

    print(f"Files analyzed: " f"{len(results)}")

    print()

    print("True-pitch median probability")

    print(f"  Current KL : " f"{np.mean(current_prob):.6f}")

    print(f"  Proposed   : " f"{np.mean(proposed_prob):.6f}")

    print()

    print("Mean 50-cent discrimination margin")

    print("  probability(true) - " "max(probability(-50c), probability(+50c))")

    print(f"  Current KL : " f"{np.mean(current_margin):+.6f}")

    print(f"  Proposed   : " f"{np.mean(proposed_margin):+.6f}")

    print()

    print("Local pitch accuracy")

    print("  Candidates: " "-100, -50, 0, +50, +100 cents")

    print(f"  Current KL : " f"{100.0 * np.mean(current_correct):.2f}%")

    print(f"  Proposed   : " f"{100.0 * np.mean(proposed_correct):.2f}%")

    print()

    improvement = proposed_margin - current_margin

    print(
        "Proposal improves 50-cent "
        f"margin in "
        f"{100.0 * np.mean(improvement > 0):.2f}% "
        "of files"
    )


# ============================================================
# Main
# ============================================================


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Compare current OpenScofo "
            "pitch templates against a "
            "higher-resolution cents-domain "
            "KL template on Orchidea "
            "ordinario samples."
        )
    )

    parser.add_argument(
        "root",
        type=Path,
        help=("Root of the Orchidea " "sample library"),
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=Path("pitch_template_results.csv"),
    )

    parser.add_argument(
        "--sample-rate",
        type=int,
        default=48000,
    )

    parser.add_argument(
        "--fft-size",
        type=int,
        default=2048,
    )

    parser.add_argument(
        "--hop-size",
        type=int,
        default=512,
    )

    parser.add_argument(
        "--harmonics",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--current-sigma",
        type=float,
        default=0.5,
        help=("Current OpenScofo " "PitchTemplateSigma " "in semitones"),
    )

    parser.add_argument(
        "--proposal-sigma-cents",
        type=float,
        default=25.0,
        help=("Gaussian sigma in cents " "for proposed template"),
    )

    parser.add_argument(
        "--amplitude-decay",
        type=float,
        default=0.5,
    )

    parser.add_argument(
        "--pitch-scaling",
        type=float,
        default=0.5,
    )

    parser.add_argument(
        "--skip-attack-ms",
        type=float,
        default=50.0,
        help=(
            "Ignore beginning of note " "when evaluating stationary " "pitch template"
        ),
    )

    parser.add_argument(
        "--rms-threshold-db",
        type=float,
        default=-30.0,
        help=(
            "Only analyze frames within "
            "this many dB of the sample's "
            "maximum frame RMS"
        ),
    )

    args = parser.parse_args()

    wav_files = find_ordinario_files(args.root)

    print(f"Found {len(wav_files)} " "WAV files inside " "'ordinario' directories.")

    results: list[FileResult] = []

    skipped = 0

    for index, path in enumerate(
        wav_files,
        start=1,
    ):

        pitch = extract_pitch_from_filename(path)

        if pitch is None:
            print(f"[skip] Could not parse " f"pitch: {path.name}")
            skipped += 1
            continue

        try:
            result = analyze_file(
                path=path,
                pitch=pitch,
                sr=args.sample_rate,
                fft_size=args.fft_size,
                hop_size=args.hop_size,
                harmonics=args.harmonics,
                current_sigma=args.current_sigma,
                proposal_sigma_cents=(args.proposal_sigma_cents),
                amplitude_decay=(args.amplitude_decay),
                pitch_scaling_factor=(args.pitch_scaling),
                skip_attack_ms=(args.skip_attack_ms),
                rms_threshold_db=(args.rms_threshold_db),
            )

        except Exception as exc:
            print(f"[error] {path}: {exc}")
            skipped += 1
            continue

        if result is None:
            print(f"[skip] No usable audio: " f"{path}")
            skipped += 1
            continue

        results.append(result)

        if result.proposed_margin_50 > result.current_margin_50:
            winner = "NEW"
        elif result.current_margin_50 > result.proposed_margin_50:
            winner = "KL"
        else:
            winner = "TIE"

        print(
            f"[{index:5d}/{len(wav_files):5d}] "
            f"{pitch:4s} "
            f"{result.frequency:9.3f} Hz | "
            f"KL={result.current.median:.5f} "
            f"(50c margin {result.current_margin_50:+.5f}) | "
            f"NEW={result.proposed.median:.5f} "
            f"(50c margin {result.proposed_margin_50:+.5f}) | "
            f"WINNER={winner}"
        )

    write_csv(
        results,
        args.output,
    )

    print_summary(results)

    print()
    print(f"CSV: {args.output}")

    if skipped:
        print(f"Skipped: {skipped}")


if __name__ == "__main__":
    main()
