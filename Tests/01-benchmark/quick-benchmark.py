#!/usr/bin/env python3

import os
import re
import json
from pathlib import Path

import numpy as np
import soundfile as sf
import OpenScofo

SR = 48000
FFT = 2048
HOP = 256
BLOCK = 64
TOLERANCE = 0.250  # seconds

os.chdir(os.path.dirname(__file__))


def natural_key(path):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", path.name)]


def list_files():
    tests = []

    # Miniaturas
    miniaturas = Path("../04-miniaturas")

    for audio in sorted(
        (miniaturas / "Audios").glob("miniatura*.mp3"),
        key=natural_key,
    ):
        name = audio.stem

        score = miniaturas / "Extras" / f"{name}.scofo"
        annotations = miniaturas / "Extras" / f"{name}.json"

        if not score.exists():
            print(f"Skipping {audio.name}: missing score")
            continue

        if not annotations.exists():
            print(f"Skipping {audio.name}: missing annotations")
            continue

        tests.append(
            {
                "audio": audio,
                "score": score,
                "annotations": annotations,
            }
        )

    # Real recordings
    for audio in sorted(Path("./real").glob("*.wav"), key=natural_key):
        name = audio.stem
        score = Path(f"./real/{name}.scofo")
        annotations = Path(f"./real/{name}.json")

        if score.exists() and annotations.exists():
            tests.append(
                {
                    "audio": audio,
                    "score": score,
                    "annotations": annotations,
                }
            )

    # Synthetic benchmark
    for audio in sorted(Path("./audios").glob("score-*.wav"), key=natural_key):
        score = audio.with_suffix(".txt")

        if score.exists():
            tests.append(
                {
                    "audio": audio,
                    "score": score,
                    "annotations": None,
                }
            )

    return tests


def test(item):
    audio_path = item["audio"]
    score_path = item["score"]
    annotations_path = item.get("annotations")

    audio, sr = sf.read(audio_path, dtype="float32")
    if sr != SR:
        raise RuntimeError(f"{audio_path}: expected {SR} Hz, got {sr} Hz")

    # soundfile stereo layout = (samples, channels)
    if audio.ndim == 2:
        audio = audio[:, 0]

    audio = np.ascontiguousarray(audio, dtype=np.float32)
    scofo = OpenScofo.OpenScofo(SR, FFT, HOP)
    if not scofo.load_score(score_path):
        raise RuntimeError(f"Could not load score: {score_path}")

    scofo.set_current_event(0)

    if annotations_path:
        with open(annotations_path, encoding="utf-8") as f:
            data = json.load(f)

        reference = {
            int(e["event"]): float(e["timestamp_seconds"]) for e in data["events"]
        }

    else:
        reference = {}

        for state in scofo.get_states():
            pos = int(state.score_pos)

            if pos > 0 and pos not in reference:
                reference[pos] = float(state.onset_expected)

    # Run follower
    detected = {}
    previous = 0

    for start in range(0, len(audio), BLOCK):
        frame = np.zeros(BLOCK, dtype=np.float32)
        chunk = audio[start : start + BLOCK]
        frame[: len(chunk)] = chunk
        scofo.process_block(frame)
        pos = int(scofo.get_current_score_position())
        if pos > 0 and pos != previous:
            detected.setdefault(pos, start / SR)
            previous = pos

    # Evaluate
    offsets = []
    missed = []
    for pos, ref in reference.items():
        if pos not in detected:
            missed.append(pos)
            continue

        offset = detected[pos] - ref
        if abs(offset) > TOLERANCE:
            missed.append(pos)
        else:
            offsets.append(offset * 1000.0)

    false_positives = sorted(set(detected) - set(reference))
    total = len(reference)
    matched = total - len(missed)
    precision = 100.0 * matched / total if total else 0.0
    mean_offset = np.mean(np.abs(offsets)) if offsets else float("nan")

    print(
        f"{audio_path.name:<20} "
        f"{matched:>3}/{total:<3} "
        f"{precision:6.2f}%  "
        f"miss={len(missed):>3}  "
        f"fp={len(false_positives):>2}  "
        f"|off|={mean_offset:6.1f} ms"
    )

    return total, matched, len(false_positives), offsets


def main():
    total_ref = 0
    total_matched = 0
    total_fp = 0
    all_offsets = []
    for item in list_files():
        ref, matched, fp, offsets = test(item)
        total_ref += ref
        total_matched += matched
        total_fp += fp
        all_offsets.extend(offsets)

    print("-" * 75)
    precision = 100.0 * total_matched / total_ref if total_ref else 0.0
    mean_offset = np.mean(np.abs(all_offsets)) if all_offsets else float("nan")

    print(
        f"TOTAL                "
        f"{total_matched}/{total_ref} "
        f"{precision:.2f}%  "
        f"miss={total_ref - total_matched}  "
        f"fp={total_fp}  "
        f"|off|={mean_offset:.1f} ms"
    )


if __name__ == "__main__":
    main()
