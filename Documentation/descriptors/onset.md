# Onset Descriptor

## `Onset Detection`
**ID**: `onset` 

Boolean value (`0` or `1`) indicating whether a sound event onset has been detected. Evaluated using the `onsetsds` library's. We highly recommend the reading of the paper _Adaptive whitening for improved real-time audio onset detection_ to better understand how the onset detection works and how to choose the best method for the music you will use (Table 1 and Table 2). 

Default method uses: 

$$MKL_n = \sum_{k=0}^{K} \log \left ( 1 + \frac{|S_{n,k}|}{|S_{n-1,k}|} \right )$$

But all methods descripted in the article can be used using the [`ONSETFUNCTION`](./../score/config.md#onsetfunction){ data-preview } on config. Where $S_{n,k}$ is a FFT bin.

## Different methods




| sigla | name                      | description                                                                                                                                                                             |
| ----- | ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **pow**   | Power                     | Measures frame-to-frame changes in signal energy. Works well for signals with clear amplitude attacks (e.g., percussion) because onsets often correspond to sudden increases in energy. |
| **pd**    | Phase Deviation           | Detects onsets by measuring irregularities in the phase progression of spectral bins. Effective when phase continuity breaks at note attacks.                                           |
| **wpd**   | Weighted Phase Deviation  | Similar to phase deviation but weights bins by magnitude, emphasizing stronger spectral components. Improves robustness when relevant partials dominate the spectrum.                   |
| **sf**    | Spectral Flux             | Measures the positive change in magnitude spectrum between consecutive frames. Commonly used for onset detection since note attacks typically introduce new spectral energy.            |
| **cd**    | Complex Domain            | Uses both magnitude and phase to estimate expected spectral evolution and measures the deviation from this prediction. Captures subtle onset cues in complex signals.                   |
| **rcd**   | Rectified Complex Domain  | A rectified version of the complex-domain method that counts only increases in deviation. Helps suppress false detections from decreases or cancellations.                              |
| **hfc** | High Frequency Content    | Emphasizes changes in high-frequency bins. Effective for percussive sounds where attacks introduce strong high-frequency components.                                                    |
| **mkl**   | Modified Kullback–Leibler | Measures divergence between spectral distributions of consecutive frames. Sensitive to structural changes in the spectrum, which often correspond to note onsets.                       |


## Considerations of use

In a review of Tables 1 and 2 from the article *Adaptive Whitening for Improved Real‑Time Audio Onset Detection*, it is possible to conclude that no single onset detection function consistently dominates across all types of audio material. Performance varies according to signal characteristics such as percussiveness, harmonic content, and polyphonic complexity. Energy-based and high-frequency methods tend to perform well for percussive signals, while approaches that incorporate phase or spectral distribution information—particularly the complex-domain method—show stronger and more stable performance across a wider range of datasets. The results also indicate that adaptive whitening generally improves detection accuracy for complex mixtures and pitched material, suggesting that preprocessing the spectrum can significantly enhance the robustness of onset detection functions.


| situation                  | best odf     | reason                                                                                                                      |
| -------------------------- | ------------ | --------------------------------------------------------------------------------------------------------------------------- |
| percussive / drums         | cd, pow, hfc | Strong transient attacks produce large energy increases and high-frequency bursts, which these methods capture effectively. |
| polyphonic pitched music   | cd           | Combines magnitude and phase prediction, allowing it to detect subtle spectral changes in harmonic textures.                |
| monophonic pitched signals | mkl, cd      | Sensitive to distribution changes in the spectrum, which helps when attacks are softer and energy change is smaller.        |
| complex mixtures           | cd, pow      | More robust to heterogeneous signals where both energy changes and spectral deviations occur.                               |
| general-purpose            | cd           | Typically the most robust overall because it models expected spectral evolution using both magnitude and phase information. |
