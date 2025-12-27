import numpy as np
import librosa

import OScofo

# Parameters
sr = 48000
fft_size = 2048
hop = 512

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
print(oscofo_desc.loudness)
print(oscofo_desc.silence_prob)

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
print(librosa_loudness)


# ╭──────────────────────────────────────╮
# │              Comparação              │
# ╰──────────────────────────────────────╯
oscofo = np.array(oscofo_desc.mfcc)
librosa_vals = np.array(librosa_mfcc)
