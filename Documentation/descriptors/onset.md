---
tags:
  - Onset Detection
---

# Onset Descriptor

Use this page to tune how OpenScofo detects attacks and new sound events.

## Overview

**ID**: `onset`

Boolean value (`0` or `1`) indicating whether an onset was detected. OpenScofo uses the `onsetsds` library; see _Adaptive Whitening for Improved Real-Time Audio Onset Detection_ for the method family.

The default method uses modified Kullback-Leibler divergence:

$$MKL_n = \sum_{k=0}^{K} \log \left ( 1 + \frac{|S_{n,k}|}{|S_{n-1,k}|} \right )$$

All methods below can be selected with [`ONSETFUNCTION`](../score/advanced-config/). Here, $S_{n,k}$ is an FFT bin.

## Reference Table

| ID | Name | Use |
| ----- | ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `pow` | Power | Energy changes; good for clear attacks. |
| `pd` | Phase deviation | Phase discontinuities at note attacks. |
| `wpd` | Weighted phase deviation | Phase deviation weighted by stronger spectral bins. |
| `sf` | Spectral flux | Positive magnitude change between frames. |
| `cd` | Complex domain | Magnitude and phase prediction; strong general-purpose choice. |
| `rcd` | Rectified complex domain | Complex-domain increases only; suppresses some false detections. |
| `hfc` | High-frequency content | Percussive attacks with strong high-frequency energy. |
| `mkl` | Modified Kullback-Leibler | Spectral-distribution changes; useful for softer pitched attacks. |

## Example

```openscofo
ONSETFUNCTION cd
BPM 72
```

## Remarks

No onset function wins for every signal. Start with `cd` for mixed material, `pow` or `hfc` for percussion, and `mkl` or `cd` for softer monophonic pitched signals.
