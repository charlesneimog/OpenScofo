# C++

```
git clone https://github.com/charlesneimog/OpenScofo.git
```

The C++ API provides direct access to the `OpenScofo` engine with minimal abstraction. It is intended for real-time systems, high-performance applications, and integration into custom audio pipelines.

## Overview

The library exposes a single main class:

- `OpenScofo::OpenScofo` — core score-following engine

It operates on blocks of audio samples and maintains an internal state representing score position and audio descriptors.

---

## Version

```cpp 
OSCOFO_VERSION_MAJOR
OSCOFO_VERSION_MINOR
OSCOFO_VERSION_PATCH
OSCOFO_BUILD_TIME
```

---

## Constructor

```cpp
OpenScofo::OpenScofo(float sampleRate, float fftSize, float hopSize);
```

* `sampleRate` — audio sampling rate
* `fftSize` — analysis window size
* `hopSize` — hop size between frames

---

## Core Workflow

### Load Score

```cpp 
bool ParseScore(std::string scorePath);
```

---

### Process Audio

```cpp
bool ProcessBlock(const std::vector<double>& audioBuffer);
```

* `audioBuffer` must match the configured block size

---

### Check State

```cpp 
bool ScoreIsLoaded();
```

---

## Configuration

```cpp
void SetAmplitudeDecay(double decay);
void SetPitchTemplateSigma(double sigma);
void SetHarmonics(int harmonics);
void SetdBTreshold(double dB);
void SetTunning(double tuning);
void SetCurrentEvent(int eventIndex);
void SetMinEntropy(double entropy);
void SetNewAudioParameters(float sr, float fftSize, float hopSize);
```

---

## ONNX Model

```cpp 
void LoadONNXModel(fs::path model, std::vector<Descriptors> descriptors);
```

Loads a neural model for extended feature estimation.

---

## Getters

### Transport / Score

```cpp 
double GetLiveBPM();
int GetEventIndex();
States GetStates();
EventActions GetEventActions(int index);
```

---

### Audio / Signal

```cpp 
double GetdBValue();
double GetPitchProb(double frequency);
std::vector<double> GetSpectrumPower() const;
double GetBlockDuration();
```

---

### Audio Descriptors

```cpp 
Description GetDescription();
Description GetAudioDescription(const std::vector<double>& audioBuffer);
```

---

### Pitch

```cpp 
PitchTemplateArray GetPitchTemplate(double frequency);
```

---

### System Parameters

```cpp 
double GetSr();
double GetFFTSize();
double GetHopSize();
```

---

## Descriptor Access

```cpp 
Descriptors GetDescriptorsEnum(const char* name);
const char* GetDescriptionId(Descriptors d);

double GetDescriptionFloat(Description& desc, Descriptors d);
std::vector<double>& GetDescriptionArray(Description& desc, Descriptors d);
```

Provides dynamic access to descriptor values.

---

## Error Handling

### Callback

```cpp 
void SetErrorCallback(
  std::function<void(const spdlog::details::log_msg&, void* data)> cb,
  void* data = nullptr
);
```

### Logging

```cpp
void SetLogLevel(spdlog::level::level_enum level);
```

### Reset

```cpp 
void ClearErrors();
```

---

## Lua Integration (Optional)

Enabled with `OSCOFO_LUA`.

```cpp 
void InitLuaModule();
bool LuaExecute(std::string code);
std::string LuaGetError();

bool LuaAddModule(std::string name, lua_CFunction func);
bool LuaAddPointer(void* pointer, const char* name);
void LuaAddPath(std::string path);
```

---

## Minimal Example

```cpp 
#include <OpenScofo.hpp>
#include <vector>

int main() {
    OpenScofo::OpenScofo scofo(44100, 1024, 512);

    scofo.ParseScore("score.txt");

    std::vector<double> buffer(512);

    while (true) {
        // fill buffer with audio input

        scofo.ProcessBlock(buffer);

        int event = scofo.GetEventIndex();
        double bpm = scofo.GetLiveBPM();
    }
}
```

---

## Notes

* Designed for real-time use; avoid allocations in the audio thread;
* `ProcessBlock` should be called with consistent timing;
* Descriptor extraction is integrated into processing;
* API is stable but the project is currently in `alpha`;
