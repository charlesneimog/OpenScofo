import numpy as np
import librosa
import matplotlib.pyplot as plt
import OpenScofo

def plot_audio_db(filepath):
    # 1. Initialize OpenScofo with your training parameters
    sr = 48000
    window_size = 2048
    hop_size = 512
    oscofo = OpenScofo.OpenScofo(sr, window_size, hop_size)

    # 2. Load the audio file
    print(f"Loading {filepath}...")
    y, _ = librosa.load(filepath, sr=sr)

    # 3. Extract dB values window by window
    db_values = []
    idx = 0
    
    while idx <= len(y) - window_size:
        window = np.asarray(y[idx : idx + window_size], dtype=np.float32)
        desc = oscofo.get_audio_description(window)
        db_values.append(desc.db)
        idx += hop_size

    # 4. Plot the results
    plt.figure(figsize=(10, 5))
    
    # Calculate time axis in seconds
    time_axis = np.arange(len(db_values)) * hop_size / sr

    # plt.plot(time_axis, db_values, label='Audio dB Level', color='black', linewidth=1.5)
    plt.scatter(time_axis, db_values, label='Frame dB Level', color='black', s=15)
    
    # Add your threshold lines for debugging
    plt.axhline(y=-60.0, color='red', linestyle='--', label='-60 dB (Old Threshold)')
    plt.axhline(y=-100.0, color='blue', linestyle='--', label='-100 dB (New Threshold)')

    file_name = filepath.split('/')[-1]
    plt.title(f"Energy (dB) per Frame: {file_name}")
    plt.xlabel("Time (seconds)")
    plt.ylabel("Decibels (dB)")
    plt.legend()
    plt.grid(True, alpha=0.4)
    plt.tight_layout()
    plt.show()

# Run it on one of your files:
plot_audio_db("/home/neimog/Downloads/Flute/key_click/Fl-key_cl-A4-f-N-N.wav")
plot_audio_db("/home/neimog/Downloads/Flute/tongue_ram/Fl-tng_ram-C3-mf-N-N.wav")
