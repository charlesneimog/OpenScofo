import OScofo
import numpy as np
import matplotlib.pyplot as plt
import music21
import os
from pathlib import Path


for midi in range(30, 86):
    p = music21.pitch.Pitch(midi)
    p.accidental.alter = round(p.accidental.alter)  # ensure accidental exists
    p.preferSharps = True
    print(p.nameWithOctave)

# root = Path("/home/neimog/Nextcloud/MusicData/Samples/Orchidea")
# files = sorted(
#     p for p in root.rglob("ordinario/**/*.wav") if "-" + NOTE + "-" in p.name
# )
