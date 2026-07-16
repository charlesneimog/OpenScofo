---
icon: custom/vamp
tags:
  - Host Integration
  - Vamp
---

# Vamp

## Overview

Use the Vamp plugin for offline descriptor extraction in hosts such as Sonic Visualiser or Audacity.

## Installation

Place the compiled Vamp plugin library (`.so`, `.dylib`, or `.dll`) in your system Vamp plugin directory. Use the installer `OpenScofo` [installer](https://github.com/charlesneimog/OpenScofo/releases/) for this.

## Minimal Descriptors Example 

1. Open an audio file in Sonic Visualiser after install `OpenScofo`.
2. Choose **Transform > Analysis by Plugin Name > Open Score Follower**.
3. Select a descriptor, such as `MFCC`, `Pitch`, or `Spectral Centroid`.

## Minimal Score Follower Example 

1. Put score files (`.scofo`) inside your `Documents` folder. 
2. Open an audio recording of the piece in Sonic Visualiser after install `OpenScofo`.
3. Choose **Transform > Analysis by Plugin Name > Open Score Follower -> Score Marks**.
4. Select the piece score file, click on `OK`.

## Reference Table

### Messages / API

Vamp is host-driven. Select descriptors from the host interface.

## Score Actions

Vamp does not supported actions. **Vamp is offline analysis only**.

## Descriptors

| Type | Descriptors |
| --- | --- |
| Vector | `MFCC`, `LogMelSpectrum`, `Melogram`, `Chroma` |
| Amplitude | max amplitude, loudness, RMS, dB, amplitude standard deviation |
| Spectral | flux, irregularity, crest, centroid, spread, flatness, HFR, rolloff |
| Pitch / time | pitch, pitch confidence, harmonicity, zero crossing rate |
| Probability | onset, silence, extended technique |
| Motion | centroid velocity |

## Remarks

Vamp plugins cannot dynamically parse score files from the host GUI during processing, so this integration exposes OpenScofo's descriptor extractor only.
