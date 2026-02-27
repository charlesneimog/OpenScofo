import os
import OpenScofo
import numpy as np
import matplotlib.pyplot as plt

os.chdir(os.path.dirname(__file__))

scofo = OpenScofo.OpenScofo(48000, 2048, 1024)
ok = scofo.parse_score("../../Tests/canticos.txt")
if not ok:
    raise Exception("Error to parse score")

states = scofo.get_states()
onset_expected = np.array([s.onset_expected for s in states], dtype=float)

block_duration = scofo.get_block_duration()
time = np.arange(0.0, onset_expected[-1], block_duration)

plt.figure()

# sigma controls width (in blocks)
sigma = 150  # increase for wider peaks

block_indices = np.arange(len(time))

for onset in onset_expected:
    center_idx = int(onset / block_duration)

    distance = block_indices - center_idx

    gaussian = np.exp(-0.5 * (distance**2) / (sigma**2)) / (sigma * np.sqrt(2 * np.pi))

    plt.plot(time, gaussian)

plt.xlim(0, onset_expected[-1])
plt.xlabel("Time (s)")
plt.title("Gaussian Peaks Around Onsets")
plt.show()
