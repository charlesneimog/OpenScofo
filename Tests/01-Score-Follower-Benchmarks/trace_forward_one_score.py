#!/usr/bin/env python3
"""Capture OpenScofo forward-state probabilities around one score jump.

Run from anywhere. By default this inspects score-32 near its 3 -> 10 jump.
The output contains one snapshot per analysis hop, including every state's
Forward, ExitProb, BestObs, and InitProb. It does not change the follower.
"""

import argparse
import json
import math
from pathlib import Path

import librosa
import numpy as np
import OpenScofo


SR = 48000
FFT = 2048
HOP = 256
BLOCK = 64
DEFAULT_CENTERS = {15: 16.959, 30: 2.2813333333333334, 32: 1.2253333333333334}


def finite_number(value):
    value = float(value)
    return value if math.isfinite(value) else None


def event_type(state):
    # FIRSTEVENT (0) is not exported by the current nanobind enum.
    if int(state.score_pos) == 0:
        return "FIRSTEVENT"
    enum = state.type
    return str(enum.name) if hasattr(enum, "name") else str(enum)


def description_snapshot(description):
    fields = ("silence", "ext", "pitch_confidence", "pitch", "onset")
    result = {}
    for field in fields:
        if hasattr(description, field):
            result[field] = finite_number(getattr(description, field))
    try:
        result["onnx"] = {
            str(label): finite_number(probability)
            for label, probability in description.onnx.items()
        }
    except (AttributeError, TypeError, ValueError):
        pass
    return result


def read_state_snapshot(scofo):
    # GetCurrentBufferIndex() returns the slot for the NEXT analysis frame:
    # GetEvent() increments Tau after writing and normalizing the current one.
    next_slot = int(scofo.get_current_buffer_index())
    states = []
    for state in scofo.get_states():
        forward = state.forward
        if not forward:
            continue
        last_slot = (next_slot - 1) % len(forward)
        states.append({
            "index": int(state.index),
            "score_pos": int(state.score_pos),
            "type": event_type(state),
            "forward": finite_number(forward[last_slot]),
            "exit_prob": finite_number(state.exit_prob[last_slot]),
            "best_obs": finite_number(state.best_obs[last_slot]),
            "init_prob": finite_number(state.init_prob),
            "onset_expected_s": finite_number(state.onset_expected),
        })
    ranked = sorted(
        states, key=lambda item: item["forward"] if item["forward"] is not None else -1.0,
        reverse=True,
    )
    return next_slot, last_slot, states, [state["index"] for state in ranked[:12]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--score", type=int, choices=sorted(DEFAULT_CENTERS), default=32)
    parser.add_argument("--audio", type=Path, help="Audio path; defaults to audios/score-N.wav")
    parser.add_argument("--score-file", type=Path, help="Score path; defaults to audios/score-N.txt")
    parser.add_argument("--center-s", type=float, help="Center of the diagnostic window")
    parser.add_argument("--half-window-s", type=float, default=0.18)
    parser.add_argument("--output", type=Path, help="Output JSON path")
    args = parser.parse_args()

    if args.half_window_s <= 0:
        parser.error("--half-window-s must be positive")

    base = Path(__file__).resolve().parent
    audio_path = args.audio or base / "audios" / f"score-{args.score}.wav"
    score_path = args.score_file or base / "audios" / f"score-{args.score}.txt"
    output_path = args.output or base / f"forward_trace_score{args.score}.json"
    center = args.center_s if args.center_s is not None else DEFAULT_CENTERS[args.score]
    window_start = max(0.0, center - args.half_window_s)
    window_end = center + args.half_window_s

    if not audio_path.is_file() or not score_path.is_file():
        parser.error(f"Missing audio or score file: {audio_path}, {score_path}")

    audio, _ = librosa.load(str(audio_path), sr=SR)
    scofo = OpenScofo.OpenScofo(SR, FFT, HOP)
    scofo.load_score(score_path)

    # This static map includes REST states, which share the preceding ScorePos.
    score_states = []
    reference_times = {}
    for state in scofo.get_states():
        pos = int(state.score_pos)
        record = {
            "index": int(state.index),
            "score_pos": pos,
            "type": event_type(state),
            "onset_expected_s": finite_number(state.onset_expected),
            "duration_beats": finite_number(state.duration),
            "score_line": int(state.line) if pos > 0 else None,
        }
        score_states.append(record)
        if pos > 0 and pos not in reference_times:
            reference_times[pos] = record["onset_expected_s"]

    frames = []
    transitions = []
    previous_pos = int(scofo.get_current_score_position())
    previous_slot = int(scofo.get_current_buffer_index())
    for start in range(0, len(audio), BLOCK):
        frame = np.zeros(BLOCK, dtype=audio.dtype)
        n_valid = min(BLOCK, len(audio) - start)
        frame[:n_valid] = audio[start:start + n_valid]
        scofo.process_block(frame)
        position = int(scofo.get_current_score_position())
        time_s = start / SR  # Same timestamp convention as the benchmark.

        if position != previous_pos:
            transitions.append({
                "time_s": time_s,
                "from_score_pos": previous_pos,
                "to_score_pos": position,
                "skipped_score_positions": (
                    list(range(previous_pos + 1, position))
                    if position > previous_pos + 1 else []
                ),
                "reference_time_s": reference_times.get(position),
            })
            previous_pos = position

        next_slot = int(scofo.get_current_buffer_index())
        if next_slot == previous_slot:
            continue  # No analysis hop occurred in this 64-sample block.
        previous_slot = next_slot
        if not window_start <= time_s <= window_end:
            continue

        next_slot, observed_slot, states, top_indices = read_state_snapshot(scofo)
        frames.append({
            "time_s": time_s,
            "score_pos": position,
            "next_buffer_slot": next_slot,
            "observed_buffer_slot": observed_slot,
            "description": description_snapshot(scofo.get_description()),
            "top_forward_indices": top_indices,
            "states": states,
        })

    if not frames:
        parser.error("No analysis frames found in the requested time window")

    result = {
        "score": args.score,
        "audio_file": str(audio_path.resolve()),
        "score_file": str(score_path.resolve()),
        "sample_rate": SR,
        "fft_size": FFT,
        "hop_size": HOP,
        "benchmark_block_size": BLOCK,
        "center_time_s": center,
        "window_start_s": window_start,
        "window_end_s": window_end,
        "score_states": score_states,
        "reference_times_s": reference_times,
        "transitions": transitions,
        "frames": frames,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)

    print(f"Saved {len(frames)} analysis frames to {output_path.resolve()}")
    for event in transitions:
        if window_start <= event["time_s"] <= window_end:
            print(
                f"{event['time_s']:.3f}s: {event['from_score_pos']} -> "
                f"{event['to_score_pos']} skipped={event['skipped_score_positions']}"
            )


if __name__ == "__main__":
    main()
