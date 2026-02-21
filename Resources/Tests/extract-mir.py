import os
import librosa
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt


from OScofo import OScofo


this_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(this_dir)

piece_file = "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Strings/Viola/ordinario/Va-ord-A4-ff-4c-N.wav"
piece_file = "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Strings/Viola/pizzicato_l_vib/Va-pizz_lv-A3-ff-4c-N.wav"
piece_file = "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Percussion/hard-plastic-mallet/middle/tt-plast_mid-mf-N.wav"
# piece_file = "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/PluckedStrings/Harp/pizzicato_bartok/Hp-pizz_bartok-G4-ff-N-N.wav"
# piece_file = "/home/neimog/Nextcloud/MusicData/Samples/Orchidea/Winds/Flute/jet_whistle/Fl-jet_wh-N-N-N-N.wav"

FFT_SIZE = 4096
HOP_SIZE = 1024
y, sr = librosa.load(piece_file, sr=None)
score_file = piece_file.replace(".mp3", ".txt")
follower = OScofo(sr, FFT_SIZE, HOP_SIZE)

ok = follower.parse_score("bwv-1013.txt")


flux = []
flat = []
harm = []
p_tonal = []

prev_p = 0.0  # for smoothing / hysteresis

for w in range(0, len(y), 1024):
    audio = y[w : w + 4096]
    if len(audio) < 4096:
        break

    desc = follower.get_audio_description(audio)

    # raw descriptors (optional logging)
    f = desc.spectral_flux
    b = desc.spectral_flatness
    h = desc.harmonicity

    flux.append(f)
    flat.append(b)
    harm.append(h)

    # ---- tonal probability (in-place) ----
    F = 1.0 - (1.0 - np.exp(-f * 8.0))  # low flux → tonal
    B = 1.0 - np.clip(b, 0.0, 1.0)  # low flatness → tonal
    H = np.clip(h * 4.0, 0.0, 1.0)  # strong peaks → tonal

    p = 0.4 * H + 0.3 * B + 0.3 * F
    p = np.clip(p, 0.0, 1.0)

    # temporal smoothing (important)
    p = 0.8 * prev_p + 0.2 * p
    prev_p = p

    p_tonal.append(p)


x = np.arange(len(flux))  # frame index (or time if you scale it)

plt.figure()
plt.plot(x, flux, label="flux")
plt.plot(x, flat, label="flat")
plt.plot(x, harm, label="harm")
plt.plot(x, p_tonal, label="tonal")

plt.xlabel("Frame")
plt.ylabel("Value")
plt.ylim(0.0, 1.0)  # ← force 0–1 scale
plt.legend()
plt.tight_layout()
plt.show()
