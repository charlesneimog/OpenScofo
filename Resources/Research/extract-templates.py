import OpenScofo
import numpy as np
import librosa
import music21
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from skopt import gp_minimize
from skopt.space import Real, Integer

# ----------------------------
# Configuration
# ----------------------------
SR = 48000
FFT_SIZE = 4096
HOP_SIZE = 1024
N_WORKERS = os.cpu_count()

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# ----------------------------
# Dataset
# ----------------------------
root = Path("/home/neimog/Nextcloud/MusicData/Samples/Orchidea")
files = sorted(p for p in root.rglob("ordinario/**/*.wav"))


def get_audioanalisys(path, scofo):
    y, _ = librosa.load(path, sr=SR)
    start = round(SR * 0.15)
    segment = y[start : start + FFT_SIZE]
    return scofo.get_audio_description(segment)


allnotes = {}
files = sorted(p for p in root.rglob("ordinario/**/*.wav"))

scofo = OpenScofo.OpenScofo(SR, FFT_SIZE, HOP_SIZE)
scofo.parse_score("../Resources/tests/canticos.txt")

# for f in files:
#     note = os.path.basename(f).split("-")[2]
#     desc = get_audioanalisys(f, scofo)
#     allnotes.setdefault(note, []).append(desc.norm_spectral_power)

# import pickle
#
# with open("allnotes.bin", "wb") as f:
# pickle.dump(allnotes, f, protocol=pickle.HIGHEST_PROTOCOL)

import pickle

with open("allnotes.bin", "rb") as f:
    allnotes = pickle.load(f)

import numpy as np
import matplotlib.pyplot as plt

templates = allnotes["C2"]

plt.figure()

# all templates
# for tpl in templates:
#     plt.plot(tpl, alpha=0.4)

# mean template
c4_mean = np.mean(templates, axis=0)
plt.plot(c4_mean, linewidth=1, label="Mean")

plt.title("C4 – norm spectral power")
plt.xlabel("Bin / Frame")
plt.ylabel("Power")
plt.legend()
plt.show()
