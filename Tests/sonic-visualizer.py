import OpenScofo
import librosa
import numpy as np
import os

os.chdir(os.path.dirname(__file__))


SR = 48000
FFTSIZE = 2048
HOPSIZE = 512

oscofo = OpenScofo.OpenScofo(SR, FFTSIZE, HOPSIZE)
x, sr = librosa.load("./assets/canticos.mp3", sr=SR)
oscofo.parse_score("./assets/canticos.txt")


event_marks = []
CURRENT_EVENT = 0

for i in range(0, len(x), HOPSIZE):
    block = x[i : i + FFTSIZE]

    if len(block) < FFTSIZE:
        padded = np.zeros(FFTSIZE, dtype=x.dtype)
        padded[: len(block)] = block
        block = padded

    ok = oscofo.process_block(block)
    event = oscofo.get_event_index()
    bpm = oscofo.get_live_bpm()

    if CURRENT_EVENT != event:
        CURRENT_EVENT = event
        t = i / sr
        event_marks.append((t, event))


# Allows to load these events on sonic visualizer
with open("events.lab", "w") as f:
    for t, ev in event_marks:
        f.write(f"{t:.6f}\tevent_{ev}\n")
