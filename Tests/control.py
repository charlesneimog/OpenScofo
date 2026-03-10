import OpenScofo
import librosa
import numpy as np
import json
import os
import time

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


def create_control(control_path="control_events.json"):
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

    results = {}

    for name, (audio, score) in tests.items():
        events, elapsed_seconds = run_scofo(audio, score)
        results[name] = {
            "elapsed_seconds": elapsed_seconds,
            "events": [{"time": t, "event": ev, "bpm": bpm} for t, ev, bpm in events],
        }
        print(f"Saved control for {name}: {elapsed_seconds:.6f}s")

    with open(control_path, "w") as f:
        json.dump(results, f, indent=4)


create_control()
