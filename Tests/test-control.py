import json
import os
import time
import librosa
import OpenScofo
import numpy as np

TIME_TOL = 1e-5
BPM_TOL = 1

SR = 48000
FFTSIZE = 2048
HOPSIZE = 512

os.chdir(os.path.dirname(__file__))


def run_scofo(audio_path, score_path):
    start = time.perf_counter()
    oscofo = OpenScofo.OpenScofo(SR, FFTSIZE, HOPSIZE)

    x, sr = librosa.load(audio_path, sr=SR)
    oscofo.parse_score(score_path)

    event_marks = []
    CURRENT_EVENT = 0

    for i in range(0, len(x), HOPSIZE):
        block = x[i : i + FFTSIZE]

        if len(block) < FFTSIZE:
            padded = np.zeros(FFTSIZE, dtype=x.dtype)
            padded[: len(block)] = block
            block = padded

        oscofo.process_block(block)

        event = oscofo.get_event_index()
        bpm = oscofo.get_live_bpm()

        if CURRENT_EVENT != event:
            CURRENT_EVENT = event
            t = i / sr
            event_marks.append((t, event, bpm))

    elapsed_seconds = time.perf_counter() - start
    return event_marks, elapsed_seconds


def compare_events(ref, cur, time_tol=TIME_TOL, bpm_tol=BPM_TOL):
    diffs = []

    if len(ref) != len(cur):
        diffs.append(f"event count mismatch: control={len(ref)} current={len(cur)}")

    n = min(len(ref), len(cur))

    for i in range(n):
        t_ref = ref[i]["time"]
        ev_ref = ref[i]["event"]
        bpm_ref = ref[i]["bpm"]

        t_cur, ev_cur, bpm_cur = cur[i]

        if (
            ev_ref != ev_cur
            or abs(t_ref - t_cur) > time_tol
            or abs(bpm_ref - bpm_cur) > bpm_tol
        ):
            diffs.append(
                {
                    "index": i,
                    "control": (t_ref, ev_ref, bpm_ref),
                    "current": (t_cur, ev_cur, bpm_cur),
                    "time_diff": abs(t_ref - t_cur),
                    "bpm_diff": abs(bpm_ref - bpm_cur),
                }
            )

    return diffs


def test_against_control(control_path="control_events.json"):
    tests = {
        "bwv-1013": (
            "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.mp3",
            "/home/neimog/Documents/Git/OpenScofo/Tests/assets/bwv-1013.txt",
        ),
        "canticos": (
            "/home/neimog/Documents/Git/OpenScofo/Tests/assets/canticos.mp3",
            "/home/neimog/Documents/Git/OpenScofo/Tests/assets/canticos.txt",
        ),
    }

    with open(control_path) as f:
        control = json.load(f)

    all_ok = True

    for name, (audio, score) in tests.items():
        current_events, elapsed_seconds = run_scofo(audio, score)

        control_entry = control[name]
        if isinstance(control_entry, dict):
            ref_events = control_entry.get("events", [])
            ref_elapsed_seconds = control_entry.get("elapsed_seconds")
        else:
            # Backward compatibility with old control format (list of events only).
            ref_events = control_entry
            ref_elapsed_seconds = None

        diffs = compare_events(ref_events, current_events)

        if diffs:
            all_ok = False
            print(f"\nFAIL: {name}")
            for d in diffs:
                print(d)
        else:
            print(f"PASS: {name}")

        if ref_elapsed_seconds is None:
            print(f"Timing [{name}] current={elapsed_seconds:.6f}s | control=not recorded")
        else:
            delta = elapsed_seconds - ref_elapsed_seconds
            print(
                f"Timing [{name}] current={elapsed_seconds:.6f}s | "
                f"control={ref_elapsed_seconds:.6f}s | delta={delta:+.6f}s"
            )

    return all_ok


test_against_control()
