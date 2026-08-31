---
icon: material/language-cpp
tags:
  - Host Integration
  - C++
---

# C++

## Overview

Use the C++ API to embed score following or descriptor extraction in an audio application, plugin, or offline analysis tool.

## Installation

Clone with submodules and build the library:

```bash
git clone --recursive https://github.com/charlesneimog/OpenScofo.git
cd OpenScofo
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target OpenScofo
```

The main header and CMake target are:

```cpp
#include <OpenScofo.hpp>
```

```cmake
add_subdirectory(path/to/OpenScofo)
target_link_libraries(my_app PRIVATE OpenScofo)
```

The build expects the Tree-sitter CLI. If needed:

```bash
npm install -g tree-sitter-cli
```

## Minimal Example

```cpp
OpenScofo::OpenScofo scofo(48000.0f, 2048.0f, 512.0f);
scofo.LoadScore("score.scofo");
scofo.ProcessBlock(audio.data(), audio.size());
```

## Reference Table

### Construction and Processing

| API | Purpose |
| --- | --- |
| `OpenScofo(float sampleRate, float fftSize, float hopSize)` | create engine |
| `LoadScore(path)` | load score |
| `ScoreIsLoaded()` | check score state |
| `ProcessBlock(buffer, n)` | process one audio block |
| `SetCurrentEvent(event)` | force score position |
| `SetCurrentSection(section)` | reset to the first event of a named section |
| `SetConfiguration(config)` | update configuration |

### State

| API | Purpose |
| --- | --- |
| `GetCurrentBPM()` | current estimated BPM |
| `GetCurrentScorePosition()` | current score event index |
| `GetCurrentStateIndex()` | internal state index |
| `GetStates()` | Markov states |
| `GetCurrentEventActions()` | actions attached to the current event |

### Descriptors and Models

| API | Purpose |
| --- | --- |
| `SetRequestedDescriptors(list)` | compute selected descriptors |
| `RequestDescriptor(descriptor)` | request one descriptor |
| `GetDescription()` | current descriptor container |
| `GetDescriptionFloat(desc, descriptor)` | scalar descriptor |
| `GetDescriptionArray(desc, descriptor)` | vector descriptor |
| `GetDescriptorsEnum(name)` | descriptor id lookup |
| `LoadONNXModel(path, descriptors)` | load custom model |

## Score Actions

Score actions are returned as `OpenScofo::ScoreAction`.

```cpp
struct ScoreAction {
    bool isLua;
    std::string Lua;
    std::string Receiver;
    std::vector<std::variant<float, int, std::string>> Args;
    bool AbsoluteTime;
    double Time;
};
```

For `sendto`, use `Receiver` and `Args`; for `luacall`, use `Lua`. If `AbsoluteTime` is false, `Time` is in beats:

```cpp
double delayMs = 60.0 / scofo.GetCurrentBPM() * action.Time * 1000.0;
```

If `AbsoluteTime` is true, `Time` is already in milliseconds.

## Descriptors

```cpp
scofo.SetRequestedDescriptors({
    OpenScofo::RMS,
    OpenScofo::CENTROID,
    OpenScofo::MFCC,
});

OpenScofo::Description desc = scofo.GetDescription();
double rms = scofo.GetDescriptionFloat(desc, OpenScofo::RMS);
std::vector<double> &mfcc = scofo.GetDescriptionArray(desc, OpenScofo::MFCC);
```

Common descriptor ids:

```text
onset loudness db maxamp rms stddev
magnitude power silence mfcc chroma logmel
zcr hfr centroid spreadhz spread_variance
crest flatness entropy rolloff flux
skewness slope kurtosis irregularity
harmonicity yin yin_confidence ext onnx
```

## Complete Example

```cpp
#include <OpenScofo.hpp>

#include <iostream>
#include <vector>

int main() {
    OpenScofo::OpenScofo scofo(48000.0f, 2048.0f, 512.0f);

    if (!scofo.LoadScore("score.scofo")) {
        std::cerr << "Could not load score\n";
        return 1;
    }

    std::vector<float> block(512);
    int lastState = -1;

    while (get_next_audio_block(block)) {
        if (!scofo.ProcessBlock(block.data(), block.size())) {
            continue;
        }

        int state = scofo.GetCurrentStateIndex();
        if (state == lastState) {
            continue;
        }

        lastState = state;
        std::cout << "event " << scofo.GetCurrentScorePosition()
                  << ", bpm " << scofo.GetCurrentBPM() << "\n";

        for (const OpenScofo::ScoreAction &action : scofo.GetCurrentEventActions()) {
            if (action.isLua) {
                std::cout << "lua: " << action.Lua << "\n";
            } else {
                std::cout << "sendto " << action.Receiver << "\n";
            }
        }
    }
}
```

Replace `get_next_audio_block(...)` with your audio input code.

## Remarks

* Keep `ProcessBlock(...)` calls in chronological order.
* Do not call `LoadScore(...)`, `SetConfiguration(...)`, or `LoadONNXModel(...)` from a real-time audio callback.
* `GetDescriptionArray(...)` returns data owned by the `Description` object you pass in.
* Lua API functions are available when OpenScofo is built with `OPENSCOFO_BUILD_WITH_LUA=ON`.
* The API is still changing while OpenScofo is in pre-alpha.
