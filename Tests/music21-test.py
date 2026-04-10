import os
import random
import string

from music21 import stream, note, meter, instrument, articulations, tempo
import subprocess
import librosa
import numpy as np
import OpenScofo
import soundfile as sf
from copy import deepcopy

curr_dir = os.path.dirname(__file__) + "/test"
os.makedirs(curr_dir, exist_ok=True)
s = "".join(random.choices(string.ascii_letters + string.digits, k=8))

measure_counter_max = 4
bpm = 60

score = stream.Score()
part = stream.Part()
part.insert(0, instrument.Flute())  # define como flauta

part.append(meter.TimeSignature("4/4"))

# rhythmic groups (each sums to 1 beat)
groups = [
    [1.0, 1.0],
    [1.0, 1.0],
    [2.0],
    [1.0],
    [0.5, 0.5],
    [1 / 3, 1 / 3, 1 / 3],
    [0.75, 0.25],
]

noteheads = ["x", "diamond", "normal", "normal", "normal"]

measure_counter = 0
time = 0.0
measure = stream.Measure(number=1)
pitch = 0

while measure_counter < measure_counter_max:
    group = random.choice(groups)
    for k in group:
        # prevent overflow of the measure
        if time + k > 4.0 + 1e-9:
            break

        prev_pitch = pitch
        while True:
            pitch = random.randint(60, 72)
            if pitch != 0 and prev_pitch != pitch:
                break

        n = note.Note(pitch, quarterLength=k)
        # random
        if (pitch > 59) and (pitch < 72):
            n.notehead = random.choice(noteheads)

        if k < 0.5 and n.notehead == "normal":
            n.articulations.append(articulations.Staccato())

        if n.notehead != "normal":
            n.volume.velocity = 0

        measure.append(n)
        time += k

    # close measure when full
    if abs(time - 4.0) < 1e-6:
        part.append(measure)
        measure_counter += 1

        # reset
        time = 0.0
        measure = stream.Measure(number=measure_counter + 1)

# append to score and export
score.append(part)
part.makeBeams(inPlace=True)

# Put tempo directly in measure 1 so MusicXML exporters keep the BPM marking.
first_measure = part.measure(1)
if first_measure is not None:
    first_measure.insert(0, tempo.MetronomeMark(number=bpm))


# ╭──────────────────────────────────────╮
# │        RENDER AUDIO AND SCORE        │
# ╰──────────────────────────────────────╯
score_original = deepcopy(score)  # your original score object

# --- generate a random base name ---
base_name = "".join(random.choices(string.ascii_letters + string.digits, k=8))

# --- 1. Prepare audio score ---
score_audio = deepcopy(score_original)
for n in score_audio.flatten().notes:
    if n.notehead != "normal":
        r = note.Rest(quarterLength=n.quarterLength)
        measure = n.getContextByClass(stream.Measure)
        if measure:
            measure.replace(n, r)

# --- 2. Write MusicXML for audio and render MP3 ---
score_audio.write("musicxml", f"{curr_dir}/{base_name}-audio.musicxml")
subprocess.run(
    [
        "mscore",
        f"{curr_dir}/{base_name}-audio.musicxml",
        "-o",
        f"{curr_dir}/{base_name}-audio.mp3",
        "--sound-profile",
        "MuseSounds",
        "--style",
        f"{curr_dir}/../flute.mss",
    ],
    check=True,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)

# --- 3. Write MusicXML for PDF display ---
score_original.write("musicxml", f"{curr_dir}/{base_name}.musicxml")
subprocess.run(
    ["mscore", f"{curr_dir}/{base_name}.musicxml", "-o", f"{curr_dir}/{base_name}.pdf"],
    check=True,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)

# --- 4. Load audio and overlay key-click samples ---
x, sr = librosa.load(f"{curr_dir}/{base_name}-audio.mp3", sr=48000)
y = x.copy()
events = score_original.flatten()
onset = 0

for n in events.notes:
    dur_sec = float(n.seconds)
    start = int(onset * sr)
    length = int(dur_sec * sr)
    end = start + length
    if n.notehead == "normal":
        onset += dur_sec
        continue

    else:
        if n.notehead == "x":
            pitch = n.pitch.nameWithOctave
            if "-" in pitch:
                pitch = n.pitch.getEnharmonic().nameWithOctave
            path = f"/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Winds/Flute/key_click/Fl-key_cl-{pitch}-f-N-N.wav"

        elif n.notehead == "diamond":
            n.transpose(-12, inPlace=True)
            pitch = n.pitch.nameWithOctave
            if "-" in pitch:
                pitch = n.pitch.getEnharmonic().nameWithOctave
            path = f"/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Winds/Flute/tongue_ram/Fl-tng_ram-{pitch}-mf-N-N.wav"
        else:
            raise Exception("Wrong notehead")

        sound, _ = librosa.load(path, sr=sr)

        if len(sound) < length:
            sound = np.pad(sound, (0, length - len(sound)))
        else:
            sound = sound[:length]

        if start >= len(y):
            onset += dur_sec
            continue
        if end > len(y):
            sound = sound[: len(y) - start]
            end = len(y)

        y[start:end] += sound

    onset += dur_sec

# --- 5. Normalize and save WAV ---
y = np.clip(y, -1.0, 1.0)
sf.write(f"{curr_dir}/{base_name}-audio.wav", y * 2, sr)
