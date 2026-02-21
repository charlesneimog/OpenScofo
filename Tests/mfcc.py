import numpy as np
import librosa

import OScofo

# Parameters
sr = 48000
fft_size = 4096
hop = 1024

# sinusoid
freq = 440.0
duration = 1.0

# OScofo internal config
n_mfcc = 13
n_mels = 40
fmin = 0.0
fmax = sr * 0.5

# ╭──────────────────────────────────────╮
# │                OScofo                │
# ╰──────────────────────────────────────╯
scofo = OScofo.OScofo(sr, fft_size, hop)
t = np.arange(int(fft_size * duration)) / sr
audio = np.sin(2 * np.pi * freq * t).astype(np.float64)
oscofo_desc = scofo.get_audio_description(audio)

for _ in range(60):
    audio_loop = np.zeros(int(fft_size * duration))
    for _ in range(15):
        freq = np.random.uniform(100, 2000)  # 100 Hz to 2 kHz
        audio_loop += np.sin(
            2 * np.pi * freq * np.arange(int(fft_size * duration)) / sr
        )
    audio_loop = (audio / np.max(np.abs(audio))).astype(np.float64)  # normalize
    _ = scofo.get_audio_description(audio)


# ╭──────────────────────────────────────╮
# │               Librosa                │
# ╰──────────────────────────────────────╯
stft = librosa.stft(audio, n_fft=fft_size, hop_length=hop, center=False)
power_spec = np.abs(stft) ** 2
mel_basis = librosa.filters.mel(
    sr=sr, n_fft=fft_size, n_mels=n_mels, fmin=fmin, fmax=fmax, htk=False
)
mel_energy = np.dot(mel_basis, power_spec)
mel_energy = np.maximum(mel_energy, 1e-30)
mel_db = 10.0 * np.log10(mel_energy)
mel_db = np.maximum(mel_db, mel_db.max() - 80.0)
mfcc_librosa = librosa.feature.mfcc(S=mel_db, n_mfcc=n_mfcc, dct_type=2)
mfcc_list = mfcc_librosa.tolist()  # gives a list of lists

librosa_mfcc = mfcc_librosa.flatten().tolist()
librosa_loudness = np.sum(mel_energy, axis=0)


# ╭──────────────────────────────────────╮
# │              Comparação              │
# ╰──────────────────────────────────────╯
oscofo = np.array(oscofo_desc.mfcc)
librosa_vals = np.array(librosa_mfcc)
print("Difference is of", np.max(np.abs(oscofo - librosa_vals)))
assert (
    np.max(np.abs(oscofo - librosa_vals)) < 0.0095
), "MFCC is too different from librosa"
