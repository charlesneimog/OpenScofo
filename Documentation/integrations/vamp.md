# Vamp Plugin

**OpenScofo** is also available as a Vamp plugin for offline audio feature extraction. It is designed to be used in host applications like Sonic Visualiser or Audacity.

!!! warning "Note on Score Following"
    **Due to inherent limitations in the architecture of Vamp plugins**, which are strictly designed for offline analysis and do not support dynamic string-based file parsing from the host GUI during processing, **this plugin does NOT perform score following.** Instead, it acts as a comprehensive **audio descriptor extractor**, exposing all of OpenScofo's internal feature extraction capabilities directly to your Vamp host.

---

## Features

The plugin analyzes incoming audio blocks (default size: 2048 samples, step size: 512 samples) and outputs a wide array of both vector and scalar descriptors. 

### Vector Descriptors (Arrays)
These output multiple values per step (e.g., bins or bands):

* **MFCC** (13 bins)

* **LogMelSpectrum / Melogram** (40 bins)

* **Chroma** (12 bins)

### Scalar Descriptors (Single Values)
These output a single continuous value per step:

* **Amplitude & Energy:** Max Amplitude, Loudness, RMS, dB, Amplitude Standard Deviation.

* **Spectral Features:** Spectral Flux, Spectral Irregularity, Spectral Crest, Spectral Centroid, Spectral Spread (Hz), Spectral Spread Variance, Spectral Flatness, High Frequency Ratio.

* **Pitch & Tone:** Pitch, Pitch Confidence, Harmonicity, Zero Crossing Rate.

* **Probabilities & Classifications:** Onset Probability, Silence Probability, Extended Technique Probability.

* **Motion:** Centroid Velocity.

---

## Usage (e.g., in Sonic Visualiser)

1.  Ensure the compiled Vamp plugin library (`.so`, `.dylib`, or `.dll`) is placed in your system's standard Vamp plugin directory.
2.  Open your audio file in Sonic Visualiser.
3.  Go to **Transform > Analysis by Plugin Name > Open Score Follower**.
4.  Select the specific descriptor you wish to visualize (e.g., `MFCC`, `Pitch`, `Spectral Centroid`, etc).
5.  The host will process the audio offline and draw the requested descriptor data over your waveform/spectrogram.
