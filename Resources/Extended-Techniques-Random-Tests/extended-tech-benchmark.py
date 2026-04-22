#!/usr/bin/env python3
"""
Extended Techniques Benchmark for OpenScofo.

Evaluation follows the protocol described in Cont (Section 8.2) and MIREX
score-following campaigns:

- Missed note:
  - reference event not detected, or
  - detected but |offset| > tolerance.
- False positive:
  - misaligned event (|offset| > tolerance), and counted as part of misses.
- Offset statistics:
  - mean absolute offset, mean signed offset, standard deviation.
- Precision metrics:
  - piecewise precision = mean_i((N_ref_i - N_miss_i) / N_ref_i)
  - overall precision = (sum_i N_ref_i - sum_i N_miss_i) / sum_i N_ref_i

Latency is fixed to zero because this benchmark runs in strict online mode.
"""

import argparse
import hashlib
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import librosa
import numpy as np
import OpenScofo

os.chdir(os.path.dirname(__file__))

# ============================================================================
# GLOBAL CONFIGURATION
# ============================================================================
SR = 48000
FFT = 2048
HOP = 256

TEST_FILES = [
    {"audio": "./score-1.wav", "score": "./score-1.txt"},
    {"audio": "./score-2.wav", "score": "./score-2.txt"},
    {"audio": "./score-3.wav", "score": "./score-3.txt"},
    {"audio": "./score-4.wav", "score": "./score-4.txt"},
    {"audio": "./score-5.wav", "score": "./score-5.txt"},
    {"audio": "./score-6.wav", "score": "./score-6.txt"},
    {"audio": "./score-7.wav", "score": "./score-7.txt"},
    {"audio": "./score-8.wav", "score": "./score-8.txt"},
    {"audio": "./score-9.wav", "score": "./score-9.txt"},
    {"audio": "./score-10.wav", "score": "./score-10.txt"},
    {"audio": "./score-11.wav", "score": "./score-11.txt"},
    {"audio": "./score-12.wav", "score": "./score-12.txt"},
    {"audio": "./score-13.wav", "score": "./score-13.txt"},
    {"audio": "./score-14.wav", "score": "./score-14.txt"},
    {"audio": "./score-15.wav", "score": "./score-15.txt"},
    {"audio": "./score-16.wav", "score": "./score-16.txt"},
    {"audio": "./score-17.wav", "score": "./score-17.txt"},
    {"audio": "./score-18.wav", "score": "./score-18.txt"},
    {"audio": "./score-19.wav", "score": "./score-19.txt"},
    {"audio": "./score-20.wav", "score": "./score-20.txt"},
]

RESULTS_PATH = "follower_validation.json"
RESULT_FILES_DIR = "mirex_results"

PROTOCOL_TOLERANCE_MS = {
    "cont-8.2": 250.0,
    "mirex-2006": 2000.0,
}
DEFAULT_PROTOCOL = "cont-8.2"
ONLINE_LATENCY_MS = 0.0


class ScoreFollowerValidator:
    """Compute and persist protocol-based evaluation metrics."""

    def __init__(
        self,
        results_path: str,
        protocol_name: str,
        tolerance_ms: float,
    ):
        self.results_path = Path(results_path)
        self.protocol_name = protocol_name
        self.tolerance_ms = float(tolerance_ms)

    def resolve_implementation_name(self, implementation_name: Optional[str]) -> str:
        """Normalize implementation name used for summary and persistence."""
        if implementation_name:
            return implementation_name
        return f"impl_{self.hash_implementation()}"

    @staticmethod
    def hash_implementation() -> str:
        """Create a short hash for the currently loaded OpenScofo module."""
        try:
            module_path = Path(OpenScofo.__file__)
            if module_path.exists():
                with open(module_path, "rb") as handle:
                    return hashlib.sha256(handle.read()).hexdigest()[:8]
        except Exception:
            pass

        return hashlib.sha256(str(datetime.now()).encode("utf-8")).hexdigest()[:8]

    def compute_piece_metrics(
        self,
        detected_events: List[Tuple[int, float]],
        expected_times: Dict[int, float],
    ) -> Dict:
        """
        Compute per-piece metrics according to the protocol.

        detected_events: list of tuples (score_pos, detected_time_seconds)
        expected_times: map score_pos -> reference_time_seconds
        """
        detected_by_pos: Dict[int, List[float]] = {}
        for pos, detected_time in detected_events:
            detected_by_pos.setdefault(pos, []).append(float(detected_time))

        expected_positions = set(expected_times.keys())

        missed_positions: List[int] = []
        misaligned_positions: List[int] = []
        matched_offsets_ms: List[float] = []

        for pos in sorted(expected_positions):
            reference_time = expected_times[pos]
            candidate_detections = detected_by_pos.get(pos, [])

            if not candidate_detections:
                missed_positions.append(pos)
                continue

            best_detection = min(
                candidate_detections,
                key=lambda value: abs(value - reference_time),
            )
            offset_ms = float((best_detection - reference_time) * 1000.0)

            if abs(offset_ms) > self.tolerance_ms:
                missed_positions.append(pos)
                misaligned_positions.append(pos)
            else:
                matched_offsets_ms.append(offset_ms)

        unexpected_positions = sorted(set(detected_by_pos.keys()) - expected_positions)

        n_reference = len(expected_positions)
        n_missed = len(missed_positions)
        n_false_positives = len(misaligned_positions)
        n_matched = n_reference - n_missed

        precision = (n_matched / n_reference) if n_reference > 0 else 0.0

        if matched_offsets_ms:
            offsets = np.array(matched_offsets_ms, dtype=float)
            mean_abs_offset_ms = float(np.mean(np.abs(offsets)))
            mean_offset_ms = float(np.mean(offsets))
            std_offset_ms = float(np.std(offsets))
        else:
            mean_abs_offset_ms = 0.0
            mean_offset_ms = 0.0
            std_offset_ms = 0.0

        return {
            "reference_events": n_reference,
            "detected_transitions": len(detected_events),
            "detected_event_tags": len(detected_by_pos),
            "matched_events": n_matched,
            "missed_notes": n_missed,
            "missed_notes_pct": (
                float((n_missed / n_reference) * 100.0) if n_reference > 0 else 0.0
            ),
            "false_positives": n_false_positives,
            "unexpected_event_tags": len(unexpected_positions),
            "precision": float(precision),
            "precision_pct": float(precision * 100.0),
            "mean_abs_offset_ms": mean_abs_offset_ms,
            "mean_offset_ms": mean_offset_ms,
            "std_offset_ms": std_offset_ms,
            "latency_ms_mean": ONLINE_LATENCY_MS,
            "protocol": self.protocol_name,
            "tolerance_ms": self.tolerance_ms,
            "missing_positions": missed_positions,
            "misaligned_positions": misaligned_positions,
            "unexpected_positions": unexpected_positions,
            "matched_offsets_ms": matched_offsets_ms,
        }

    def save_run(
        self,
        piece_results: List[Dict],
        global_metrics: Dict,
        implementation_name: Optional[str] = None,
    ) -> None:
        """Persist per-piece history entries and run-level global summary."""
        if implementation_name is None:
            implementation_name = f"impl_{self.hash_implementation()}"

        timestamp = datetime.now().isoformat()

        if self.results_path.exists():
            with open(self.results_path, "r", encoding="utf-8") as handle:
                storage = json.load(handle)
        else:
            storage = {}

        history = storage.setdefault("history", [])
        run_history = storage.setdefault("run_history", [])

        for piece in piece_results:
            history.append(
                {
                    "timestamp": timestamp,
                    "implementation": implementation_name,
                    "protocol": self.protocol_name,
                    "tolerance_ms": self.tolerance_ms,
                    "audio_file": piece["audio_file"],
                    "score_file": piece["score_file"],
                    "metrics": piece["metrics"],
                }
            )

        run_history.append(
            {
                "timestamp": timestamp,
                "implementation": implementation_name,
                "protocol": self.protocol_name,
                "tolerance_ms": self.tolerance_ms,
                "files_processed": len(piece_results),
                "global_metrics": global_metrics,
                "pieces": [
                    {
                        "audio_file": piece["audio_file"],
                        "score_file": piece["score_file"],
                    }
                    for piece in piece_results
                ],
            }
        )

        storage["schema_version"] = "mirex_cont_v1"
        storage["last_updated"] = timestamp

        with open(self.results_path, "w", encoding="utf-8") as handle:
            json.dump(storage, handle, indent=2)

        print(f"\nResults saved to {self.results_path}")


def write_result_file(
    output_path: Path,
    detected_events: List[Tuple[int, float]],
    expected_times: Dict[int, float],
) -> None:
    """
    Write MIREX-style ASCII output with four columns:
    1) estimated onset time in performance [ms]
    2) detection time [ms]
    3) score start time [ms] (used as event identifier)
    4) score position integer
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lines: List[str] = []
    for pos, detected_time in detected_events:
        onset_ms = detected_time * 1000.0
        detection_ms = detected_time * 1000.0
        score_ms = expected_times.get(pos, -1.0) * 1000.0
        lines.append(f"{onset_ms:.3f} {detection_ms:.3f} {score_ms:.3f} {pos}")

    output_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def process_audio_file(
    audio_path: str,
    score_path: str,
) -> Tuple[List[Tuple[int, float]], Dict[int, float]]:
    """
    Process one audio file with OpenScofo.

    Returns:
        detected_events: list[(score_pos, detected_time_seconds)]
        expected_times: dict[score_pos] = reference_time_seconds
    """
    print(f"\n--- Processing {audio_path} ---")

    # Load audio at benchmark sample rate.
    audio, _ = librosa.load(audio_path, sr=SR)

    #audio = audio * 2
    scofo = OpenScofo.OpenScofo(SR, FFT, HOP)
    scofo.parse_score(Path(score_path))

    # Build reference tag -> time map from parser states.
    expected_times: Dict[int, float] = {}
    for state in scofo.get_states():
        if not hasattr(state, "score_pos") or not hasattr(state, "onset_expected"):
            continue
        try:
            pos = int(state.score_pos)
            expected_times[pos] = float(state.onset_expected)
        except (TypeError, ValueError):
            continue

    n_samples = len(audio)
    prev_pos: Optional[int] = None
    detected_events: List[Tuple[int, float]] = []

    for start in range(0, n_samples, 64):
        end = min(start + 64, n_samples)
        frame = np.zeros(64, dtype=audio.dtype)
        valid = end - start
        if valid > 0:
            frame[:valid] = audio[start:end]

        scofo.process_block(frame)
        pos = int(scofo.get_current_score_position())

        if pos < 0 or pos == prev_pos:
            continue

        detected_time = start / SR
        detected_events.append((pos, detected_time))

        reference_time = expected_times.get(pos)
        if reference_time is None:
            print(
                f"pos={pos:03d}  det={detected_time:8.3f}s  ref=   N/A   offset=   N/A"
            )
        else:
            offset_ms = (detected_time - reference_time) * 1000.0
            print(
                f"pos={pos:03d}  det={detected_time:8.3f}s  "
                f"ref={reference_time:8.3f}s  offset={offset_ms:+8.2f} ms"
            )

        prev_pos = pos

    detected_positions = {position for position, _ in detected_events}
    expected_positions = set(expected_times.keys())
    not_reported = expected_positions - detected_positions
    unexpected = detected_positions - expected_positions

    print(f"Detected transitions: {len(detected_events)}")
    print(f"Reference events: {len(expected_positions)}")
    print(
        f"Reference tags with no detection: {len(not_reported)}"
        if not_reported
        else "Reference tags with no detection: 0"
    )
    print(
        f"Unexpected detected tags: {len(unexpected)}"
        if unexpected
        else "Unexpected detected tags: 0"
    )

    return detected_events, expected_times


def compute_global_metrics(piece_results: List[Dict]) -> Dict:
    """Compute global metrics over all processed files."""
    if not piece_results:
        return {
            "files_processed": 0,
            "total_reference_events": 0,
            "total_missed_notes": 0,
            "total_false_positives": 0,
            "overall_precision": 0.0,
            "overall_precision_pct": 0.0,
            "piecewise_precision": 0.0,
            "piecewise_precision_pct": 0.0,
            "global_mean_abs_offset_ms": 0.0,
            "piecewise_mean_abs_offset_ms": 0.0,
            "global_mean_offset_ms": 0.0,
            "global_std_offset_ms": 0.0,
            "average_latency_ms": ONLINE_LATENCY_MS,
        }

    total_reference_events = sum(
        piece["metrics"]["reference_events"] for piece in piece_results
    )
    total_missed_notes = sum(
        piece["metrics"]["missed_notes"] for piece in piece_results
    )
    total_false_positives = sum(
        piece["metrics"]["false_positives"] for piece in piece_results
    )

    overall_precision = (
        (total_reference_events - total_missed_notes) / total_reference_events
        if total_reference_events > 0
        else 0.0
    )

    piecewise_precision_values = [
        piece["metrics"]["precision"] for piece in piece_results
    ]
    piecewise_precision = (
        float(np.mean(piecewise_precision_values))
        if piecewise_precision_values
        else 0.0
    )

    all_offsets: List[float] = []
    piecewise_abs_offsets: List[float] = []

    for piece in piece_results:
        metrics = piece["metrics"]
        all_offsets.extend(metrics.get("matched_offsets_ms", []))
        if metrics.get("matched_events", 0) > 0:
            piecewise_abs_offsets.append(metrics["mean_abs_offset_ms"])

    if all_offsets:
        offsets = np.array(all_offsets, dtype=float)
        global_mean_abs_offset_ms = float(np.mean(np.abs(offsets)))
        global_mean_offset_ms = float(np.mean(offsets))
        global_std_offset_ms = float(np.std(offsets))
    else:
        global_mean_abs_offset_ms = 0.0
        global_mean_offset_ms = 0.0
        global_std_offset_ms = 0.0

    piecewise_mean_abs_offset_ms = (
        float(np.mean(piecewise_abs_offsets)) if piecewise_abs_offsets else 0.0
    )

    return {
        "files_processed": len(piece_results),
        "total_reference_events": int(total_reference_events),
        "total_missed_notes": int(total_missed_notes),
        "total_false_positives": int(total_false_positives),
        "overall_precision": float(overall_precision),
        "overall_precision_pct": float(overall_precision * 100.0),
        "piecewise_precision": float(piecewise_precision),
        "piecewise_precision_pct": float(piecewise_precision * 100.0),
        "global_mean_abs_offset_ms": global_mean_abs_offset_ms,
        "piecewise_mean_abs_offset_ms": piecewise_mean_abs_offset_ms,
        "global_mean_offset_ms": global_mean_offset_ms,
        "global_std_offset_ms": global_std_offset_ms,
        "average_latency_ms": ONLINE_LATENCY_MS,
    }


def load_run_history(results_path: Path) -> List[Dict]:
    """Load previously stored run history entries from results JSON."""
    if not results_path.exists():
        return []

    try:
        with open(results_path, "r", encoding="utf-8") as handle:
            storage = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return []

    run_history = storage.get("run_history", [])
    if not isinstance(run_history, list):
        return []

    return [entry for entry in run_history if isinstance(entry, dict)]


def _as_float(value, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_int(value, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def select_best_run(
    run_history: List[Dict],
    protocol_name: str,
    tolerance_ms: float,
    files_processed: int,
) -> Tuple[Optional[Dict], str]:
    """Select best run, preferring same protocol/tolerance and file coverage."""
    runs_with_metrics = [
        run for run in run_history if isinstance(run.get("global_metrics"), dict)
    ]
    if not runs_with_metrics:
        return None, "no-history"

    def same_tolerance(run: Dict) -> bool:
        return abs(_as_float(run.get("tolerance_ms"), tolerance_ms) - tolerance_ms) < 1e-9

    exact_scope = [
        run
        for run in runs_with_metrics
        if run.get("protocol") == protocol_name
        and same_tolerance(run)
        and _as_int(run.get("files_processed"), -1) == files_processed
    ]

    protocol_scope = [
        run
        for run in runs_with_metrics
        if run.get("protocol") == protocol_name and same_tolerance(run)
    ]

    if exact_scope:
        candidates = exact_scope
        scope_label = "same protocol, same tolerance, same files processed"
    elif protocol_scope:
        candidates = protocol_scope
        scope_label = "same protocol and tolerance"
    else:
        candidates = runs_with_metrics
        scope_label = "all recorded runs"

    def run_rank_key(run: Dict) -> Tuple[float, float, float, int, str]:
        metrics = run.get("global_metrics", {})
        overall_precision = _as_float(metrics.get("overall_precision"), 0.0)
        mean_abs_offset = _as_float(metrics.get("global_mean_abs_offset_ms"), float("inf"))
        missed_notes = _as_float(metrics.get("total_missed_notes"), float("inf"))
        processed = _as_int(run.get("files_processed"), 0)
        timestamp = str(run.get("timestamp", ""))

        # Higher precision is better; lower offset/missed are better.
        return (
            overall_precision,
            -mean_abs_offset,
            -missed_notes,
            processed,
            timestamp,
        )

    return max(candidates, key=run_rank_key), scope_label


def compare_runs(current_metrics: Dict, best_metrics: Dict) -> Dict:
    """Build numeric deltas for current run against best run."""
    return {
        "delta_mean_abs_offset_ms": _as_float(
            current_metrics.get("global_mean_abs_offset_ms"), 0.0
        )
        - _as_float(best_metrics.get("global_mean_abs_offset_ms"), 0.0),
        "delta_mean_offset_ms": _as_float(current_metrics.get("global_mean_offset_ms"), 0.0)
        - _as_float(best_metrics.get("global_mean_offset_ms"), 0.0),
        "delta_std_offset_ms": _as_float(current_metrics.get("global_std_offset_ms"), 0.0)
        - _as_float(best_metrics.get("global_std_offset_ms"), 0.0),
        "delta_missed_notes": _as_int(current_metrics.get("total_missed_notes"), 0)
        - _as_int(best_metrics.get("total_missed_notes"), 0),
        "delta_false_positives": _as_int(current_metrics.get("total_false_positives"), 0)
        - _as_int(best_metrics.get("total_false_positives"), 0),
    }


def format_percent_vs_best(current_pct: float, best_pct: float) -> str:
    """Format precision comparison as relative percent better/worse."""
    if abs(current_pct - best_pct) < 1e-12:
        return (
            f"matches best (current={current_pct:.2f}%, best={best_pct:.2f}%)"
        )

    if best_pct <= 0.0:
        if current_pct > 0.0:
            return (
                "is above a zero best baseline "
                f"(current={current_pct:.2f}%, best={best_pct:.2f}%)"
            )
        return (
            "cannot be compared against zero best baseline "
            f"(current={current_pct:.2f}%, best={best_pct:.2f}%)"
        )

    relative = abs((current_pct - best_pct) / best_pct * 100.0)
    relation = "better" if current_pct > best_pct else "worse"
    of_best = (current_pct / best_pct) * 100.0

    return (
        f"is {relative:.2f}% {relation} than best "
        f"({of_best:.2f}% of best, current={current_pct:.2f}%, best={best_pct:.2f}%)"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Evaluate OpenScofo alignment following Cont/MIREX metrics."
    )
    parser.add_argument(
        "--protocol",
        choices=sorted(PROTOCOL_TOLERANCE_MS.keys()),
        default=DEFAULT_PROTOCOL,
        help="Evaluation protocol. Defines default miss tolerance.",
    )
    parser.add_argument(
        "--tolerance-ms",
        type=float,
        default=None,
        help="Override tolerance (ms). If omitted, protocol default is used.",
    )
    parser.add_argument(
        "--results-path",
        default=RESULTS_PATH,
        help="JSON file where piece history and run summaries are stored.",
    )
    parser.add_argument(
        "--implementation-name",
        default=None,
        help="Optional implementation label. Defaults to OpenScofo hash.",
    )
    parser.add_argument(
        "--write-result-files",
        action="store_true",
        help="Write MIREX-style ASCII output files for each processed piece.",
    )
    parser.add_argument(
        "--result-files-dir",
        default=RESULT_FILES_DIR,
        help="Directory for per-piece MIREX-style ASCII result files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    tolerance_ms = (
        float(args.tolerance_ms)
        if args.tolerance_ms is not None
        else PROTOCOL_TOLERANCE_MS[args.protocol]
    )

    print("=" * 72)
    print("EXTENDED TECHNIQUES BENCHMARK - MIREX/CONT PROTOCOL")
    print("=" * 72)
    print(f"Protocol: {args.protocol}")
    print(f"Missing/misaligned tolerance: {tolerance_ms:.1f} ms")
    print("Missed = no event or |offset| > tolerance")
    print("False positive = misaligned event (subset of misses)")
    print("Latency metric fixed to zero (strict online)")
    print()

    try:
        module_time = Path(OpenScofo.__file__).stat().st_mtime
        print(f"OpenScofo module: {Path(OpenScofo.__file__).name}")
        print(
            "Last modified: "
            f"{datetime.fromtimestamp(module_time).strftime('%Y-%m-%d %H:%M:%S')}"
        )
    except Exception:
        pass
    print()

    validator = ScoreFollowerValidator(
        results_path=args.results_path,
        protocol_name=args.protocol,
        tolerance_ms=tolerance_ms,
    )

    piece_results: List[Dict] = []

    for test_file in TEST_FILES:
        audio_path = test_file["audio"]
        score_path = test_file["score"]

        if not Path(audio_path).exists():
            print(f"[WARN] Missing file: {audio_path}")
            continue
        if not Path(score_path).exists():
            print(f"[WARN] Missing file: {score_path}")
            continue

        detected_events, expected_times = process_audio_file(audio_path, score_path)
        metrics = validator.compute_piece_metrics(detected_events, expected_times)

        piece_results.append(
            {
                "audio_file": audio_path,
                "score_file": score_path,
                "metrics": metrics,
            }
        )

        if args.write_result_files:
            result_name = f"{Path(audio_path).stem}.result.txt"
            output_path = Path(args.result_files_dir) / result_name
            write_result_file(output_path, detected_events, expected_times)
            print(f"Result file written: {output_path}")

        print(
            "Piece metrics: "
            f"ref={metrics['reference_events']}  "
            f"missed={metrics['missed_notes']}  "
            f"fp={metrics['false_positives']}  "
            f"precision={metrics['precision_pct']:.2f}%  "
            f"mean|offset|={metrics['mean_abs_offset_ms']:.2f} ms  "
            f"mean={metrics['mean_offset_ms']:.2f} ms  "
            f"std={metrics['std_offset_ms']:.2f} ms"
        )

    print("\n" + "=" * 72)
    print("GLOBAL SUMMARY")
    print("=" * 72)

    if not piece_results:
        print("No file was processed successfully.")
        return

    global_metrics = compute_global_metrics(piece_results)
    current_implementation = validator.resolve_implementation_name(
        args.implementation_name
    )
    current_run = {
        "timestamp": datetime.now().isoformat(),
        "implementation": current_implementation,
        "protocol": args.protocol,
        "tolerance_ms": tolerance_ms,
        "files_processed": len(piece_results),
        "global_metrics": global_metrics,
    }
    historical_runs = load_run_history(Path(args.results_path))
    best_run, best_scope = select_best_run(
        run_history=historical_runs + [current_run],
        protocol_name=args.protocol,
        tolerance_ms=tolerance_ms,
        files_processed=len(piece_results),
    )

    print(f"Files processed: {global_metrics['files_processed']}")
    print(f"Total reference events: {global_metrics['total_reference_events']}")
    print(f"Total missed notes: {global_metrics['total_missed_notes']}")
    print(f"Total false positives: {global_metrics['total_false_positives']}")
    print(f"Overall precision: {global_metrics['overall_precision_pct']:.2f}%")
    print(f"Piecewise precision: {global_metrics['piecewise_precision_pct']:.2f}%")
    print(
        "Global mean absolute offset: "
        f"{global_metrics['global_mean_abs_offset_ms']:.2f} ms"
    )
    print(
        "Piecewise mean absolute offset: "
        f"{global_metrics['piecewise_mean_abs_offset_ms']:.2f} ms"
    )
    print(f"Global mean offset: {global_metrics['global_mean_offset_ms']:.2f} ms")
    print(f"Global offset std: {global_metrics['global_std_offset_ms']:.2f} ms")
    print(f"Average latency: {global_metrics['average_latency_ms']:.2f} ms")

    if best_run is not None:
        best_metrics = best_run["global_metrics"]
        comparison = compare_runs(global_metrics, best_metrics)
        overall_precision_cmp = format_percent_vs_best(
            _as_float(global_metrics.get("overall_precision_pct"), 0.0),
            _as_float(best_metrics.get("overall_precision_pct"), 0.0),
        )
        piecewise_precision_cmp = format_percent_vs_best(
            _as_float(global_metrics.get("piecewise_precision_pct"), 0.0),
            _as_float(best_metrics.get("piecewise_precision_pct"), 0.0),
        )
        best_impl = str(best_run.get("implementation", "unknown"))
        best_time = str(best_run.get("timestamp", "N/A"))

        print("\n" + "-" * 72)
        print("BEST IMPLEMENTATION")
        print("-" * 72)
        print(f"Comparison scope: {best_scope}")
        print(f"Best implementation: {best_impl}")
        print(f"Best run timestamp: {best_time}")
        print(
            "Best metrics: "
            f"overall_precision={_as_float(best_metrics.get('overall_precision_pct')):.2f}%  "
            f"piecewise_precision={_as_float(best_metrics.get('piecewise_precision_pct')):.2f}%  "
            f"mean|offset|={_as_float(best_metrics.get('global_mean_abs_offset_ms')):.2f} ms  "
            f"missed={_as_int(best_metrics.get('total_missed_notes'))}  "
            f"fp={_as_int(best_metrics.get('total_false_positives'))}"
        )

        print("\nCurrent vs best:")
        print(f"Current implementation: {current_implementation}")
        print(f"Best implementation: {best_impl}")
        print(f"Overall precision: current {overall_precision_cmp}")
        print(f"Piecewise precision: current {piecewise_precision_cmp}")
        print(
            f"Delta mean absolute offset: {comparison['delta_mean_abs_offset_ms']:+.2f} ms"
        )
        print(f"Delta mean offset: {comparison['delta_mean_offset_ms']:+.2f} ms")
        print(f"Delta offset std: {comparison['delta_std_offset_ms']:+.2f} ms")
        print(f"Delta missed notes: {comparison['delta_missed_notes']:+d}")
        print(f"Delta false positives: {comparison['delta_false_positives']:+d}")

    validator.save_run(
        piece_results=piece_results,
        global_metrics=global_metrics,
        implementation_name=current_implementation,
    )


if __name__ == "__main__":
    main()

