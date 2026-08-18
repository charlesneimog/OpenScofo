#!/usr/bin/env python3
"""
Extended Techniques Benchmark for OpenScofo.

Evaluation follows the protocol described in Cont (Section 8.2) and MIREX
score-following campaigns:

- Missed note:
  - reference event not detected, or
  - detected but |offset| > tolerance.
- False positive:
    - spurious detection with no corresponding reference event.
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
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

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
    {"audio": "./audios/score-1.wav", "score": "./audios/score-1.txt"},
    {"audio": "./audios/score-2.wav", "score": "./audios/score-2.txt"},
    {"audio": "./audios/score-3.wav", "score": "./audios/score-3.txt"},
    {"audio": "./audios/score-4.wav", "score": "./audios/score-4.txt"},
    {"audio": "./audios/score-5.wav", "score": "./audios/score-5.txt"},
    {"audio": "./audios/score-6.wav", "score": "./audios/score-6.txt"},
    {"audio": "./audios/score-7.wav", "score": "./audios/score-7.txt"},
    {"audio": "./audios/score-8.wav", "score": "./audios/score-8.txt"},
    {"audio": "./audios/score-9.wav", "score": "./audios/score-9.txt"},
    {"audio": "./audios/score-10.wav", "score": "./audios/score-10.txt"},
    {"audio": "./audios/score-11.wav", "score": "./audios/score-11.txt"},
    {"audio": "./audios/score-12.wav", "score": "./audios/score-12.txt"},
    {"audio": "./audios/score-13.wav", "score": "./audios/score-13.txt"},
    {"audio": "./audios/score-14.wav", "score": "./audios/score-14.txt"},
    {"audio": "./audios/score-15.wav", "score": "./audios/score-15.txt"},
    {"audio": "./audios/score-16.wav", "score": "./audios/score-16.txt"},
    {"audio": "./audios/score-17.wav", "score": "./audios/score-17.txt"},
    {"audio": "./audios/score-18.wav", "score": "./audios/score-18.txt"},
    {"audio": "./audios/score-19.wav", "score": "./audios/score-19.txt"},
    {"audio": "./audios/score-20.wav", "score": "./audios/score-20.txt"},
    {"audio": "./audios/score-21.wav", "score": "./audios/score-21.txt"},
    {"audio": "./audios/score-22.wav", "score": "./audios/score-22.txt"},
    {"audio": "./audios/score-23.wav", "score": "./audios/score-23.txt"},
    {"audio": "./audios/score-24.wav", "score": "./audios/score-24.txt"},
    {"audio": "./audios/score-25.wav", "score": "./audios/score-25.txt"},
    {"audio": "./audios/score-26.wav", "score": "./audios/score-26.txt"},
    {"audio": "./audios/score-27.wav", "score": "./audios/score-27.txt"},
    {"audio": "./audios/score-28.wav", "score": "./audios/score-28.txt"},
    {"audio": "./audios/score-29.wav", "score": "./audios/score-29.txt"},
    {"audio": "./audios/score-30.wav", "score": "./audios/score-30.txt"},
    {"audio": "./audios/score-31.wav", "score": "./audios/score-31.txt"},
    {"audio": "./audios/score-32.wav", "score": "./audios/score-32.txt"},
    {"audio": "./audios/score-33.wav", "score": "./audios/score-33.txt"},
    {"audio": "./audios/score-34.wav", "score": "./audios/score-34.txt"},
    {"audio": "./audios/score-35.wav", "score": "./audios/score-35.txt"},
    {"audio": "./audios/score-36.wav", "score": "./audios/score-36.txt"},
    {"audio": "./audios/score-37.wav", "score": "./audios/score-37.txt"},
    {"audio": "./audios/score-38.wav", "score": "./audios/score-38.txt"},
    {"audio": "./audios/score-39.wav", "score": "./audios/score-39.txt"},
    {"audio": "./audios/score-40.wav", "score": "./audios/score-40.txt"}
]

COLOR_RED = "\033[91m"
COLOR_GREEN = "\033[92m"
COLOR_RESET = "\033[0m"

RESULTS_PATH = "follower_validation.json"
RESULT_FILES_DIR = "mirex_results"

PROTOCOL_TOLERANCE_MS = {
    "cont-8.2": 250.0,
    "mirex-2006": 2000.0,
}

DEFAULT_PROTOCOL = "cont-8.2"
ONLINE_LATENCY_MS = 0.0
SCHEMA_VERSION = "mirex_2006_v2"

# Bulletproof fallback map assuming standard C++ enum 0-indexing.
FALLBACK_EVENT_TYPES = {
    0: "REST",
    1: "NOTE",
    2: "CHORD",
    3: "TRILL",
    4: "MULTI",
    5: "PTECH",
    6: "UTECH",
}


def _get_enum_name(enum_obj: Any) -> str:
    """Extract clean string representation from pybind11/nanobind enums with integer fallback."""
    # 1. Try native nanobind name resolution
    if hasattr(enum_obj, "name") and enum_obj.name:
        return enum_obj.name

    # 2. Try raw integer memory cast (bypasses unexported enum bindings)
    try:
        val = int(enum_obj)
        if val in FALLBACK_EVENT_TYPES:
            return FALLBACK_EVENT_TYPES[val]
    except (TypeError, ValueError):
        pass

    # 3. String parsing for mangled objects
    s = str(enum_obj)
    if "." in s:
        return s.split(".")[-1].replace(">", "").strip()
    return s


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
        base = Path(OpenScofo.__file__).resolve().parent
        h = hashlib.sha256()

        for path in sorted(base.rglob("*")):
            if path.is_file():
                h.update(path.relative_to(base).as_posix().encode())
                with open(path, "rb") as f:
                    for chunk in iter(lambda: f.read(8192), b""):
                        h.update(chunk)
        return h.hexdigest()[:16]

    def compute_piece_metrics_mirex(
        self,
        detected_events: List[Tuple[int, float]],
        expected_times: Dict[int, float],
        tolerance_ms: float,
    ) -> Dict:
        detected_by_pos: Dict[int, List[float]] = {}  # MIREX-COMPLIANT
        for pos, detected_time in detected_events:
            detected_by_pos.setdefault(pos, []).append(
                float(detected_time)
            )  # MIREX-COMPLIANT

        expected_positions = set(expected_times.keys())  # MIREX-COMPLIANT

        missed_positions: List[int] = []  # MIREX-COMPLIANT
        false_positive_positions: List[int] = []  # MIREX-COMPLIANT
        misaligned_positions: List[int] = []  # MIREX-COMPLIANT
        matched_offsets_ms: List[float] = []  # MIREX-COMPLIANT

        for pos in sorted(expected_positions):
            reference_time = expected_times[pos]  # MIREX-COMPLIANT
            candidate_detections = detected_by_pos.get(pos, [])  # MIREX-COMPLIANT

            if not candidate_detections:
                missed_positions.append(pos)  # MIREX-COMPLIANT
                continue

            candidate_detections.sort()  # MIREX-COMPLIANT
            first_detection = candidate_detections[0]  # MIREX-COMPLIANT

            offset_ms = float(
                (first_detection - reference_time) * 1000.0
            )  # MIREX-COMPLIANT

            if abs(offset_ms) > float(tolerance_ms):  # MIREX-COMPLIANT
                missed_positions.append(pos)  # MIREX-COMPLIANT
                misaligned_positions.append(pos)  # MIREX-COMPLIANT
            else:
                matched_offsets_ms.append(offset_ms)  # MIREX-COMPLIANT

        unexpected_positions = sorted(
            set(detected_by_pos.keys()) - expected_positions
        )  # MIREX-COMPLIANT
        false_positive_positions.extend(unexpected_positions)  # MIREX-COMPLIANT

        n_reference = len(expected_positions)  # MIREX-COMPLIANT
        n_missed = len(missed_positions)  # MIREX-COMPLIANT
        n_false_positives = len(false_positive_positions)  # MIREX-COMPLIANT
        n_matched = n_reference - n_missed  # MIREX-COMPLIANT

        precision = (
            ((n_reference - n_missed) / n_reference) if n_reference > 0 else 0.0
        )  # MIREX-COMPLIANT

        if matched_offsets_ms:
            offsets = np.array(matched_offsets_ms, dtype=float)  # MIREX-COMPLIANT
            mean_abs_offset_ms = float(np.mean(np.abs(offsets)))  # MIREX-COMPLIANT
            mean_offset_ms = float(np.mean(offsets))  # MIREX-COMPLIANT
            std_offset_ms = float(np.std(offsets))  # MIREX-COMPLIANT
        else:
            mean_abs_offset_ms = 0.0  # MIREX-COMPLIANT
            mean_offset_ms = 0.0  # MIREX-COMPLIANT
            std_offset_ms = 0.0  # MIREX-COMPLIANT

        return {
            "reference_events": n_reference,  # MIREX-COMPLIANT
            "matched_events": n_matched,  # MIREX-COMPLIANT
            "missed_notes": n_missed,  # MIREX-COMPLIANT
            "missed_notes_pct": (
                float((n_missed / n_reference) * 100.0) if n_reference > 0 else 0.0
            ),  # MIREX-COMPLIANT
            "false_positives": n_false_positives,  # MIREX-COMPLIANT
            "precision": float(precision),  # MIREX-COMPLIANT
            "precision_pct": float(precision * 100.0),  # MIREX-COMPLIANT
            "mean_abs_offset_ms": mean_abs_offset_ms,  # MIREX-COMPLIANT
            "mean_offset_ms": mean_offset_ms,  # MIREX-COMPLIANT
            "std_offset_ms": std_offset_ms,  # MIREX-COMPLIANT
            "latency_ms_mean": ONLINE_LATENCY_MS,  # MIREX-COMPLIANT
            "protocol": self.protocol_name,  # MIREX-COMPLIANT
            "tolerance_ms": float(tolerance_ms),  # MIREX-COMPLIANT
            "missing_positions": missed_positions,  # MIREX-COMPLIANT
            "misaligned_positions": misaligned_positions,  # MIREX-COMPLIANT
            "unexpected_positions": unexpected_positions,  # MIREX-COMPLIANT
            "false_positive_positions": false_positive_positions,  # MIREX-COMPLIANT
            "matched_offsets_ms": matched_offsets_ms,  # MIREX-COMPLIANT
        }  # MIREX-COMPLIANT

    def save_run(
        self,
        piece_results: List[Dict],
        global_metrics: Dict,
        implementation_name: Optional[str] = None,
    ) -> None:
        """Persist per-piece history entries and run-level global summary."""
        if implementation_name is None:
            implementation_name = (
                f"impl_{self.hash_implementation()}"  # MIREX-COMPLIANT
            )

        timestamp = datetime.now().isoformat()  # MIREX-COMPLIANT

        if self.results_path.exists():
            with open(
                self.results_path, "r", encoding="utf-8"
            ) as handle:  # MIREX-COMPLIANT
                storage = json.load(handle)  # MIREX-COMPLIANT
        else:
            storage = {}  # MIREX-COMPLIANT

        history = storage.setdefault("history", [])  # MIREX-COMPLIANT
        run_history = storage.setdefault("run_history", [])  # MIREX-COMPLIANT

        for piece in piece_results:
            history.append(  # MIREX-COMPLIANT
                {
                    "schema_version": SCHEMA_VERSION,  # MIREX-COMPLIANT
                    "timestamp": timestamp,  # MIREX-COMPLIANT
                    "implementation": implementation_name,  # MIREX-COMPLIANT
                    "protocol": self.protocol_name,  # MIREX-COMPLIANT
                    "tolerance_ms": self.tolerance_ms,  # MIREX-COMPLIANT
                    "audio_file": piece["audio_file"],  # MIREX-COMPLIANT
                    "score_file": piece["score_file"],  # MIREX-COMPLIANT
                    "metrics": piece["metrics"],  # MIREX-COMPLIANT
                }
            )

        run_history.append(  # MIREX-COMPLIANT
            {
                "schema_version": SCHEMA_VERSION,  # MIREX-COMPLIANT
                "timestamp": timestamp,  # MIREX-COMPLIANT
                "implementation": implementation_name,  # MIREX-COMPLIANT
                "protocol": self.protocol_name,  # MIREX-COMPLIANT
                "tolerance_ms": self.tolerance_ms,  # MIREX-COMPLIANT
                "files_processed": len(piece_results),  # MIREX-COMPLIANT
                "global_metrics": global_metrics,  # MIREX-COMPLIANT
                "pieces": [  # MIREX-COMPLIANT
                    {
                        "audio_file": piece["audio_file"],  # MIREX-COMPLIANT
                        "score_file": piece["score_file"],  # MIREX-COMPLIANT
                    }
                    for piece in piece_results
                ],
            }
        )

        storage["schema_version"] = SCHEMA_VERSION  # MIREX-COMPLIANT
        storage["last_updated"] = timestamp  # MIREX-COMPLIANT

        with open(
            self.results_path, "w", encoding="utf-8"
        ) as handle:  # MIREX-COMPLIANT
            json.dump(storage, handle, indent=2)  # MIREX-COMPLIANT

        print(f"\nResults saved to {self.results_path}")  # MIREX-COMPLIANT


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
    tolerance_ms: float,
) -> Tuple[List[Tuple[int, float]], Dict[int, float], Dict[int, Dict[str, str]]]:
    """
    Process one audio file with OpenScofo.

    Returns:
        detected_events: list[(score_pos, detected_time_seconds)]
        expected_times: dict[score_pos] = reference_time_seconds
        score_context: dict[score_pos] = {"type": current, "prev_type": previous}
    """
    print(f"\n--- Processing {audio_path} ---")

    # Load audio at benchmark sample rate.
    audio, _ = librosa.load(audio_path, sr=SR)

    scofo = OpenScofo.OpenScofo(SR, FFT, HOP)
    scofo.load_score(Path(score_path))

    expected_times: Dict[int, float] = {}
    score_context: Dict[int, Dict[str, str]] = {}
    prev_type_str = "START_OF_SCORE"

    for state in scofo.get_states():
        if not hasattr(state, "score_pos") or not hasattr(state, "onset_expected"):
            continue
        try:
            pos = int(state.score_pos)

            # Extract robust state type for precise context analysis
            current_type_str = _get_enum_name(getattr(state, "type", "UNKNOWN"))

            # Only register the expected onset of the FIRST state mapped to this score_pos
            if pos not in expected_times:
                expected_times[pos] = float(state.onset_expected)
                score_context[pos] = {
                    "type": current_type_str,
                    "prev_type": prev_type_str,
                }

            # Advance the topological memory
            prev_type_str = current_type_str

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

            # Check absolute offset against the provided tolerance
            if abs(offset_ms) > tolerance_ms:
                color = COLOR_RED
            else:
                color = COLOR_GREEN

            # Print with color and reset at the end of the line
            print(
                f"{color}pos={pos:03d}  det={detected_time:8.3f}s  "
                f"ref={reference_time:8.3f}s  offset={offset_ms:+8.2f} ms{COLOR_RESET}"
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

    return detected_events, expected_times, score_context


def compute_global_metrics_mirex(piece_results: List[Dict]) -> Dict:
    if not piece_results:
        return {
            "files_processed": 0,  # MIREX-COMPLIANT
            "total_reference_events": 0,  # MIREX-COMPLIANT
            "total_missed_notes": 0,  # MIREX-COMPLIANT
            "total_false_positives": 0,  # MIREX-COMPLIANT
            "overall_precision": 0.0,  # MIREX-COMPLIANT
            "overall_precision_pct": 0.0,  # MIREX-COMPLIANT
            "piecewise_precision": 0.0,  # MIREX-COMPLIANT
            "piecewise_precision_pct": 0.0,  # MIREX-COMPLIANT
            "global_mean_abs_offset_ms": 0.0,  # MIREX-COMPLIANT
            "piecewise_mean_abs_offset_ms": 0.0,  # MIREX-COMPLIANT
            "global_mean_offset_ms": 0.0,  # MIREX-COMPLIANT
            "global_std_offset_ms": 0.0,  # MIREX-COMPLIANT
            "average_latency_ms": ONLINE_LATENCY_MS,  # MIREX-COMPLIANT
        }  # MIREX-COMPLIANT

    total_reference_events = sum(
        piece["metrics"]["reference_events"] for piece in piece_results
    )  # MIREX-COMPLIANT
    total_missed_notes = sum(
        piece["metrics"]["missed_notes"] for piece in piece_results
    )  # MIREX-COMPLIANT
    total_false_positives = sum(
        piece["metrics"]["false_positives"] for piece in piece_results
    )  # MIREX-COMPLIANT

    overall_precision = (
        ((total_reference_events - total_missed_notes) / total_reference_events)
        if total_reference_events > 0
        else 0.0
    )  # MIREX-COMPLIANT

    piecewise_precision_values = [
        piece["metrics"]["precision"] for piece in piece_results
    ]  # MIREX-COMPLIANT
    piecewise_precision = (
        float(np.mean(piecewise_precision_values))
        if piecewise_precision_values
        else 0.0
    )  # MIREX-COMPLIANT

    all_offsets: List[float] = []  # MIREX-COMPLIANT

    for piece in piece_results:
        metrics = piece["metrics"]  # MIREX-COMPLIANT
        all_offsets.extend(metrics.get("matched_offsets_ms", []))  # MIREX-COMPLIANT

    if all_offsets:
        offsets = np.array(all_offsets, dtype=float)  # MIREX-COMPLIANT
        global_mean_abs_offset_ms = float(np.mean(np.abs(offsets)))  # MIREX-COMPLIANT
        global_mean_offset_ms = float(np.mean(offsets))  # MIREX-COMPLIANT
        global_std_offset_ms = float(np.std(offsets))  # MIREX-COMPLIANT
    else:
        global_mean_abs_offset_ms = 0.0  # MIREX-COMPLIANT
        global_mean_offset_ms = 0.0  # MIREX-COMPLIANT
        global_std_offset_ms = 0.0  # MIREX-COMPLIANT

    piecewise_mean_abs_offset_ms = global_mean_abs_offset_ms  # MIREX-COMPLIANT

    return {
        "files_processed": len(piece_results),  # MIREX-COMPLIANT
        "total_reference_events": int(total_reference_events),  # MIREX-COMPLIANT
        "total_missed_notes": int(total_missed_notes),  # MIREX-COMPLIANT
        "total_false_positives": int(total_false_positives),  # MIREX-COMPLIANT
        "overall_precision": float(overall_precision),  # MIREX-COMPLIANT
        "overall_precision_pct": float(overall_precision * 100.0),  # MIREX-COMPLIANT
        "piecewise_precision": float(piecewise_precision),  # MIREX-COMPLIANT
        "piecewise_precision_pct": float(
            piecewise_precision * 100.0
        ),  # MIREX-COMPLIANT
        "global_mean_abs_offset_ms": global_mean_abs_offset_ms,  # MIREX-COMPLIANT
        "piecewise_mean_abs_offset_ms": piecewise_mean_abs_offset_ms,  # MIREX-COMPLIANT
        "global_mean_offset_ms": global_mean_offset_ms,  # MIREX-COMPLIANT
        "global_std_offset_ms": global_std_offset_ms,  # MIREX-COMPLIANT
        "average_latency_ms": ONLINE_LATENCY_MS,  # MIREX-COMPLIANT
    }  # MIREX-COMPLIANT


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

    return [  # MIREX-COMPLIANT
        entry
        for entry in run_history
        if isinstance(entry, dict) and entry.get("schema_version") == SCHEMA_VERSION
    ]  # MIREX-COMPLIANT


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
    runs_with_metrics = [  # MIREX-COMPLIANT
        run
        for run in run_history
        if isinstance(run.get("global_metrics"), dict)
        and run.get("schema_version") == SCHEMA_VERSION
    ]  # MIREX-COMPLIANT
    if not runs_with_metrics:
        return None, "no-history"

    def same_tolerance(run: Dict) -> bool:
        return (
            abs(_as_float(run.get("tolerance_ms"), tolerance_ms) - tolerance_ms) < 1e-9
        )

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
        mean_abs_offset = _as_float(
            metrics.get("global_mean_abs_offset_ms"), float("inf")
        )
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
        "delta_mean_offset_ms": _as_float(
            current_metrics.get("global_mean_offset_ms"), 0.0
        )
        - _as_float(best_metrics.get("global_mean_offset_ms"), 0.0),
        "delta_std_offset_ms": _as_float(
            current_metrics.get("global_std_offset_ms"), 0.0
        )
        - _as_float(best_metrics.get("global_std_offset_ms"), 0.0),
        "delta_missed_notes": _as_int(current_metrics.get("total_missed_notes"), 0)
        - _as_int(best_metrics.get("total_missed_notes"), 0),
        "delta_false_positives": _as_int(
            current_metrics.get("total_false_positives"), 0
        )
        - _as_int(best_metrics.get("total_false_positives"), 0),
    }


def format_percent_vs_best(current_pct: float, best_pct: float) -> str:
    """Format precision comparison as relative percent better/worse."""
    if abs(current_pct - best_pct) < 1e-12:
        return f"matches best (current={current_pct:.2f}%, best={best_pct:.2f}%)"

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


def process_audio_file_worker(
    test_file: Dict,
    tolerance_ms: float,
    protocol_name: str,
) -> Tuple[Optional[Dict], str]:
    """
    Worker function for multiprocessing pool.
    Processes one audio file and returns results and status.

    Returns:
        Tuple of (piece_result_dict, status_message)
    """
    audio_path = test_file["audio"]
    score_path = test_file["score"]

    try:
        if not Path(audio_path).exists():
            return None, f"[WARN] Missing file: {audio_path}"
        if not Path(score_path).exists():
            return None, f"[WARN] Missing file: {score_path}"

        # Suppress worker process output to avoid interleaving
        import io

        with redirect_stdout(io.StringIO()):
            detected_events, expected_times, score_context = process_audio_file(
                audio_path, score_path, tolerance_ms
            )

        validator = ScoreFollowerValidator(  # MIREX-COMPLIANT
            results_path=RESULTS_PATH,  # MIREX-COMPLIANT
            protocol_name=protocol_name,  # MIREX-COMPLIANT
            tolerance_ms=tolerance_ms,  # MIREX-COMPLIANT
        )  # MIREX-COMPLIANT

        metrics = validator.compute_piece_metrics_mirex(  # MIREX-COMPLIANT
            detected_events, expected_times, tolerance_ms  # MIREX-COMPLIANT
        )  # MIREX-COMPLIANT

        piece_result = {
            "audio_file": audio_path,
            "score_file": score_path,
            "metrics": metrics,
            "detected_events": detected_events,
            "expected_times": expected_times,
        }

        status_msg = (
            f"✓ {Path(audio_path).name}: "
            f"ref={metrics['reference_events']}  "
            f"missed={metrics['missed_notes']}  "
            f"fp={metrics['false_positives']}  "
            f"precision={metrics['precision_pct']:.2f}%  "
            f"mean|offset|={metrics['mean_abs_offset_ms']:.2f}ms  "
            f"std={metrics['std_offset_ms']:.2f}ms"
        )

        return piece_result, status_msg

    except Exception as e:
        return None, f"✗ {Path(audio_path).name}: {str(e)}"


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
    parser.add_argument(  # MIREX-COMPLIANT
        "--mirex-strict",  # MIREX-COMPLIANT
        dest="mirex_strict",  # MIREX-COMPLIANT
        action="store_true",  # MIREX-COMPLIANT
        default=True,  # MIREX-COMPLIANT
        help="Force 250 ms tolerance and MIREX-only metrics.",  # MIREX-COMPLIANT
    )  # MIREX-COMPLIANT
    parser.add_argument(  # MIREX-COMPLIANT
        "--no-mirex-strict",  # MIREX-COMPLIANT
        dest="mirex_strict",  # MIREX-COMPLIANT
        action="store_false",  # MIREX-COMPLIANT
        help="Allow protocol-default tolerance and extra diagnostics.",  # MIREX-COMPLIANT
    )  # MIREX-COMPLIANT
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers. Defaults to CPU count.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    tolerance_ms = (
        250.0
        if args.mirex_strict
        else (
            float(args.tolerance_ms)
            if args.tolerance_ms is not None
            else PROTOCOL_TOLERANCE_MS[args.protocol]
        )
    )  # MIREX-COMPLIANT

    # Determine number of workers
    import os as os_module  # MIREX-COMPLIANT

    num_workers = (
        args.workers if args.workers is not None else os_module.cpu_count()
    )  # MIREX-COMPLIANT

    print("=" * 72)
    print("EXTENDED TECHNIQUES BENCHMARK - MIREX/CONT PROTOCOL")
    print("=" * 72)
    print(f"Protocol: {args.protocol}")
    print(f"Missing/misaligned tolerance: {tolerance_ms:.1f} ms")
    print(f"MIREX strict mode: {'enabled' if args.mirex_strict else 'disabled'}")
    print("Missed = no event or |offset| > tolerance")
    print(
        "False positive = spurious detection with no reference event"
    )  # MIREX-COMPLIANT
    print("Latency metric fixed to zero (strict online)")
    print(f"Parallel workers: {num_workers}")
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

    piece_results: List[Dict] = []

    print("Processing audio files with multiprocessing...")
    print("-" * 72)

    # Process files in parallel using ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit all tasks
        futures = {
            executor.submit(
                process_audio_file_worker, test_file, tolerance_ms, args.protocol
            ): test_file  # MIREX-COMPLIANT
            for test_file in TEST_FILES
        }

        # Collect results as they complete
        for future in as_completed(futures):
            test_file = futures[future]
            try:
                piece_result, status_msg = future.result()
                print(status_msg)

                if piece_result is not None:
                    piece_results.append(piece_result)

                    # Write result file if requested
                    if args.write_result_files:
                        result_name = (
                            f"{Path(piece_result['audio_file']).stem}.result.txt"
                        )
                        output_path = Path(args.result_files_dir) / result_name
                        write_result_file(
                            output_path,
                            piece_result["detected_events"],
                            piece_result["expected_times"],
                        )
                        print(f"Result file written: {output_path}")

            except Exception as e:
                print(f"[ERROR] Task failed for {test_file}: {str(e)}")

    print()
    print("=" * 72)
    print("TABLE 3 - MIREX-COMPLIANT SUMMARY")
    print("=" * 72)

    if not piece_results:
        print("No file was processed successfully.")
        return

    global_metrics = compute_global_metrics_mirex(piece_results)  # MIREX-COMPLIANT
    validator = ScoreFollowerValidator(
        results_path=args.results_path,
        protocol_name=args.protocol,
        tolerance_ms=tolerance_ms,
    )
    current_implementation = validator.resolve_implementation_name(
        args.implementation_name
    )
    current_run = {
        "schema_version": SCHEMA_VERSION,
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

    rows = sorted(
        piece_results, key=lambda item: Path(item["audio_file"]).name
    )  # MIREX-COMPLIANT
    header = f"{'File':<22} {'Ref':>6} {'Miss':>6} {'FP':>6} {'Mean|Off|(ms)':>15} {'MeanOff(ms)':>12} {'StdOff(ms)':>11} {'Precision(%)':>13}"  # MIREX-COMPLIANT
    print(header)  # MIREX-COMPLIANT
    print("-" * len(header))  # MIREX-COMPLIANT
    for piece in rows:  # MIREX-COMPLIANT
        metrics = piece["metrics"]  # MIREX-COMPLIANT
        print(  # MIREX-COMPLIANT
            f"{Path(piece['audio_file']).name:<22} "  # MIREX-COMPLIANT
            f"{metrics['reference_events']:>6d} "  # MIREX-COMPLIANT
            f"{metrics['missed_notes']:>6d} "  # MIREX-COMPLIANT
            f"{metrics['false_positives']:>6d} "  # MIREX-COMPLIANT
            f"{metrics['mean_abs_offset_ms']:>15.1f} "  # MIREX-COMPLIANT
            f"{metrics['mean_offset_ms']:>12.1f} "  # MIREX-COMPLIANT
            f"{metrics['std_offset_ms']:>11.1f} "  # MIREX-COMPLIANT
            f"{metrics['precision_pct']:>13.2f}"  # MIREX-COMPLIANT
        )  # MIREX-COMPLIANT

    print("=" * 72)  # MIREX-COMPLIANT
    print("GLOBAL:")  # MIREX-COMPLIANT
    print(
        f"Total Reference Events: {global_metrics['total_reference_events']}"
    )  # MIREX-COMPLIANT
    print(
        f"Total Missed Notes:     {global_metrics['total_missed_notes']}"
    )  # MIREX-COMPLIANT
    print(
        f"Total False Positives:   {global_metrics['total_false_positives']}"
    )  # MIREX-COMPLIANT
    print(
        f"Overall Precision:      {global_metrics['overall_precision_pct']:.2f}%"
    )  # MIREX-COMPLIANT
    print(
        f"Piecewise Precision:    {global_metrics['piecewise_precision_pct']:.2f}%"
    )  # MIREX-COMPLIANT
    print(
        f"Global Mean |Offset|:   {global_metrics['global_mean_abs_offset_ms']:.1f} ms"
    )  # MIREX-COMPLIANT
    print(
        f"Global Mean Offset:     {global_metrics['global_mean_offset_ms']:+.1f} ms"
    )  # MIREX-COMPLIANT
    print(
        f"Global Std Offset:      {global_metrics['global_std_offset_ms']:.1f} ms"
    )  # MIREX-COMPLIANT
    print(
        f"Average Latency:        {global_metrics['average_latency_ms']:.2f} ms"
    )  # MIREX-COMPLIANT

    if best_run is not None:  # MIREX-COMPLIANT
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

        print("\n" + "-" * 72)  # MIREX-COMPLIANT
        print("BEST IMPLEMENTATION")  # MIREX-COMPLIANT
        print("-" * 72)  # MIREX-COMPLIANT
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
