# C++

The C++ API gives direct access to the `OpenScofo` engine. Use it when you want to embed score following or descriptor extraction in your own audio application, plugin, or offline analysis tool.

## Build

Clone the repository with submodules:

```bash
git clone --recursive https://github.com/charlesneimog/OpenScofo.git
cd OpenScofo
```

Build the static library with CMake:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target OpenScofo
```

The build expects the Tree-sitter CLI to be available. If CMake cannot find it, install it with:

```bash
npm install -g tree-sitter-cli
```

The main target is:

```cmake
OpenScofo
```

The main header is:

```cpp
#include <OpenScofo.hpp>
```

If you use OpenScofo from another CMake project, add the repository as a subdirectory and link to the target:

```cmake
add_subdirectory(path/to/OpenScofo)
target_link_libraries(my_app PRIVATE OpenScofo)
```

## Basic Workflow

Create an engine with the audio parameters used by your audio stream:

```cpp
OpenScofo::OpenScofo scofo(48000.0f, 2048.0f, 512.0f);
```

Then choose one of two modes:

* **Score following**: call `LoadScore(...)`, then process audio blocks.
* **Descriptor extraction**: request descriptors, then process audio blocks without loading a score.

Audio is processed with a pointer and a sample count:

```cpp
bool ok = scofo.ProcessBlock(audio.data(), audio.size());
```

`ProcessBlock` is explicitly instantiated for `float` and `double` buffers.

## Score Following

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

        const int state = scofo.GetCurrentStateIndex();
        if (state != lastState) {
            lastState = state;

            const int event = scofo.GetCurrentScorePosition();
            const double bpm = scofo.GetCurrentBPM();

            std::cout << "event " << event << ", bpm " << bpm << "\n";

            for (const OpenScofo::ScoreAction &action : scofo.GetCurrentEventActions()) {
                if (action.isLua) {
                    std::cout << "lua: " << action.Lua << "\n";
                } else {
                    std::cout << "sendto " << action.Receiver << "\n";
                }
            }
        }
    }
}
```

Replace `get_next_audio_block(...)` with your own audio input code. For real-time use, avoid allocations and heavy work in the audio callback. A common pattern is to process audio in the callback, then send event/action data to another thread.

## Descriptor Extraction

Request descriptors before processing audio:

```cpp
#include <OpenScofo.hpp>

#include <iostream>
#include <vector>

int main() {
    OpenScofo::OpenScofo scofo(48000.0f, 2048.0f, 512.0f);

    scofo.SetRequestedDescriptors({
        OpenScofo::RMS,
        OpenScofo::CENTROID,
        OpenScofo::MFCC,
    });

    std::vector<double> block(512);

    while (get_next_audio_block(block)) {
        scofo.ProcessBlock(block.data(), block.size());

        OpenScofo::Description desc = scofo.GetDescription();
        double rms = scofo.GetDescriptionFloat(desc, OpenScofo::RMS);
        double centroid = scofo.GetDescriptionFloat(desc, OpenScofo::CENTROID);
        std::vector<double> &mfcc = scofo.GetDescriptionArray(desc, OpenScofo::MFCC);

        std::cout << "rms " << rms << ", centroid " << centroid
                  << ", mfcc count " << mfcc.size() << "\n";
    }
}
```

You can also request one descriptor at a time:

```cpp
scofo.RequestDescriptor(OpenScofo::CHROMA);
```

## Public API

### Construction

```cpp
OpenScofo::OpenScofo(float sampleRate, float fftSize, float hopSize);
```

### Main Methods

```cpp
bool LoadScore(std::filesystem::path scorePath);
bool ScoreIsLoaded();
template <typename T> bool ProcessBlock(const T *audioBuffer, size_t n);

void LoadONNXModel(std::filesystem::path model,
                   std::vector<OpenScofo::Descriptors> descriptors);

void SetCurrentEvent(int event);
void SetConfiguration(OpenScofo::Configuration &config);
void SetRequestedDescriptors(std::vector<OpenScofo::Descriptors> descriptors);
void RequestDescriptor(OpenScofo::Descriptors descriptor);
```

`SetConfiguration(...)` updates audio and analysis settings. Do not call it while another thread is processing audio.

### Score State

```cpp
double GetCurrentBPM();
int GetCurrentScorePosition();
int GetCurrentStateIndex();
int GetCurrentBufferIndex();

OpenScofo::States &GetStates();
OpenScofo::EventActions GetCurrentEventActions();
```

`GetCurrentScorePosition()` returns the current score event index. `GetCurrentStateIndex()` returns the internal state index, which is useful for detecting a new event/state transition.

### Audio and Timing

```cpp
double GetSr();
double GetFFTSize();
double GetHopSize();
double GetBlockDuration();
double GetPitchProb(double frequency);
OpenScofo::PitchTemplateArray GetPitchTemplate(double frequency);
OpenScofo::Configuration GetConfiguration();
```

### Descriptors

```cpp
OpenScofo::Description GetDescription();

OpenScofo::Descriptors GetDescriptorsEnum(const char *name);
const char *GetDescriptionId(OpenScofo::Descriptors descriptor);

double GetDescriptionFloat(OpenScofo::Description &desc,
                           OpenScofo::Descriptors descriptor);

std::vector<double> &GetDescriptionArray(OpenScofo::Description &desc,
                                         OpenScofo::Descriptors descriptor);
```

Scalar descriptors should be read with `GetDescriptionFloat(...)`. Vector descriptors such as `MFCC`, `CHROMA`, `MELOGRAM`, and `MAGNITUDE` should be read with `GetDescriptionArray(...)`.

Common descriptor ids include:

```text
onset loudness db maxamp rms stddev
magnitude power silence mfcc chroma logmel
zcr hfr centroid spreadhz spread_variance
crest flatness entropy rolloff flux
skewness slope kurtosis irregularity
harmonicity yin yin_confidence ext onnx
```

The full descriptor reference is available in [Descriptors](../descriptors/index.md).

### Configuration

```cpp
OpenScofo::Configuration config = scofo.GetConfiguration();
config.SR = 48000;
config.FFTSize = 2048;
config.HOPSize = 512;
config.TunningA4 = 440.0;
config.dBTreshold = -60;
scofo.SetConfiguration(config);
```

Important fields include:

* `SR`, `FFTSize`, `HOPSize`: audio analysis parameters.
* `TunningA4`: reference tuning.
* `PitchTemplateSigma`, `PitchTemplateHarmonics`: pitch template settings.
* `MFCCMels`, `MFCCCount`: MFCC settings.
* `dBTreshold`: silence threshold.
* `YINThreshold`, `YINMinFrequency`, `YINMaxFrequency`: pitch detection settings.
* `RequestedDescriptors`: descriptor list computed during processing.

## Score Actions

Score actions are returned as `OpenScofo::ScoreAction` values:

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

For `sendto` actions, use `Receiver` and `Args`. For `luacall` actions, use `Lua`. If `AbsoluteTime` is false, `Time` is expressed in beats and should be converted using the current BPM:

```cpp
double delayMs = 60.0 / scofo.GetCurrentBPM() * action.Time * 1000.0;
```

If `AbsoluteTime` is true, `Time` is already in milliseconds. In score files, `delay 2 sec` is converted to `2000`, and `delay 2000 ms` stays `2000`.

## ONNX Models

Load a custom ONNX model with the descriptor order used during training:

```cpp
scofo.LoadONNXModel("model.onnx", {
    OpenScofo::RMS,
    OpenScofo::CENTROID,
    OpenScofo::MFCC,
});
```

When ONNX descriptors are requested, results are stored in:

```cpp
OpenScofo::Description desc = scofo.GetDescription();
for (const auto &[label, value] : desc.ONNX) {
    std::cout << label << ": " << value << "\n";
}
```

## Errors and Logging

Register a callback to receive OpenScofo log messages:

```cpp
scofo.SetErrorCallback(
    [](const spdlog::details::log_msg &msg, void *data) {
        std::string text(msg.payload.data(), msg.payload.size());
        std::cerr << "[OpenScofo] " << text << "\n";
    });
```

Set the logging level:

```cpp
scofo.SetLogLevel(spdlog::level::info);
```

Clear stored error state:

```cpp
scofo.ClearErrors();
```

## Lua

When OpenScofo is built with `OPENSCOFO_BUILD_WITH_LUA=ON`, CMake defines `OPENSCOFO_LUA` and the C++ API exposes the embedded Lua runtime:

```cpp
bool LuaExecute(std::string code);
std::string LuaGetError();

bool LuaAddModule(std::string name, lua_CFunction func);
bool LuaAddPointer(void *pointer, const char *name);
void LuaAddPath(std::string path);
```

Score-level `LUA { ... }` code can be retrieved with:

```cpp
std::string code = scofo.GetLuaCode();
```

See [Lua for Interactive Actions](../score/lua.md) for the Lua-side API.

## Notes

* Keep `ProcessBlock(...)` calls in chronological order and use a stable sample rate, FFT size, and hop size.
* `LoadScore(...)`, `SetConfiguration(...)`, and `LoadONNXModel(...)` can allocate and should not be called from a real-time audio callback.
* `GetDescriptionArray(...)` returns a reference to data owned by the `Description` object you pass in, so keep that object alive while using the reference.
* The API is still changing while OpenScofo is in pre-alpha.
