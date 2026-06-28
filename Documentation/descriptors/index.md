# Listening Module

One of the primary objectives of `OpenScofo` is to extend the range of sonic capabilities supported by score-following systems. In particular, the framework is designed to handle complex timbral events, including extended and percussive instrumental techniques, which are often not well captured in traditional approaches.

To support this goal, `OpenScofo` integrates a machine learning module based on the `ONNX` runtime and `libonnx`, along with a set of audio descriptors. These components enable both the training and deployment of models for the recognition and classification of diverse sound events. Through this approach, the system can be adapted to a wide variety of sonic materials and performance contexts.

!!! tip "You can test the descriptors online"
    Check the link [testing OpenScofo descriptors](https://charlesneimog.github.io/OpenScofo/descriptors/testing/){target="_blank"}.

The following section presents the audio descriptors used during the training process, together with their corresponding mathematical formulations. These descriptors form the feature representation used by the machine learning models within the `OpenScofo` framework, and can be accessed through the `py.o.train` (check [AI Model](ai.md)) object or the `train` model available in the `OpenScofo` Python module.

For easier integration with AI-oriented environments such as Python, `OpenScofo` descriptors aim to be compatible with both `librosa` and `essentia`. `librosa` is the primary reference due to its strong support for spectral descriptors, although it lacks coverage of some other descriptor types. `essentia` is more comprehensive, but less widely used and somewhat less streamlined.

* Descriptors compatible with `librosa` are marked with the icon: :custom-librosa:
* Descriptors compatible with `essentia` are marked with the icon: :custom-essentia:
