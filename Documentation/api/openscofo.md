---
title: OpenScofo API
---

# OpenScofo API

Reference for the main `OpenScofo` C++ class.

!!! info

    This page is generated automatically from the OpenScofo source code.
    Do not edit it manually.

## Construction

### `OpenScofo`

```cpp
OpenScofo::OpenScofo::OpenScofo(float Sr, float WindowSize, float HopSize)
```

Initialize OpenScofo processing pipeline.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Sr` | `float` | Sampling rate (Hz) |
| `WindowSize` | `float` |  |
| `HopSize` | `float` | Hop size (samples) |

!!! note

    Initializes Forward model, MIR extractor, and score handler.

!!! note

    Configures global spdlog logger (overwrites default).

## Configuration

### `SetConfiguration`

```cpp
void OpenScofo::OpenScofo::SetConfiguration(Configuration &Config)
```

Apply configuration settings to the system.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Config` | `Configuration &` | Configuration object |

!!! note

    Updates internal modules according to provided configuration.

### `GetSr`

```cpp
double OpenScofo::OpenScofo::GetSr()
```

Get current sampling rate.

#### Returns

Sampling rate in Hz

### `GetFFTSize`

```cpp
double OpenScofo::OpenScofo::GetFFTSize()
```

Get FFT window size.

#### Returns

FFT size in samples

### `GetHopSize`

```cpp
double OpenScofo::OpenScofo::GetHopSize()
```

Get hop size used for frame processing.

#### Returns

Hop size in samples

### `GetConfiguration`

```cpp
Configuration OpenScofo::OpenScofo::GetConfiguration()
```

Retrieve current system configuration.

#### Returns

Copy of the current configuration object

!!! note

    Returned by value (snapshot, not a live reference).

## Score Following

### `LoadScore`

```cpp
bool OpenScofo::OpenScofo::LoadScore(fs::path ScorePath)
```

Load and initialize a score from file.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `ScorePath` | `fs::path` | Path to score file |

#### Returns

true if loading and initialization succeeded, false otherwise

!!! note

    Resets internal state and reinitializes processing pipeline.

!!! note

    May load and validate an ONNX model if present in the score.

### `ScoreIsLoaded`

```cpp
bool OpenScofo::OpenScofo::ScoreIsLoaded()
```

Check if a score is currently loaded.

#### Returns

true if score data is available, false otherwise

### `SetCurrentEvent`

```cpp
void OpenScofo::OpenScofo::SetCurrentEvent(int Event)
```

Set active score event and reset decoding state.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Event` | `int` | Event index in the loaded score (0 = reset) |

!!! note

    Resets forward model state, buffers, and descriptors.

!!! note

    Updates current score position based on event mapping if valid.

### `SetCurrentSection`

```cpp
bool OpenScofo::OpenScofo::SetCurrentSection(const std::string &Section)
```

Select a score section and reset decoding at its first state.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Section` | `const std::string &` | Section identifier without quotes |

#### Returns

true when the section exists, false otherwise

### `GetCurrentBPM`

```cpp
double OpenScofo::OpenScofo::GetCurrentBPM()
```

Get estimated current tempo.

#### Returns

Current BPM estimate from the forward model.

### `GetCurrentScorePosition`

```cpp
int OpenScofo::OpenScofo::GetCurrentScorePosition()
```

Get current position in the score (following Antescofo, Rest does not count for this).

#### Returns

Score position index computed by the forward model.

### `GetCurrentStateIndex`

```cpp
int OpenScofo::OpenScofo::GetCurrentStateIndex()
```

Get current state index from the forward model.

#### Returns

Index of the active internal state.

### `GetCurrentEventActions`

```cpp
EventActions OpenScofo::OpenScofo::GetCurrentEventActions()
```

Get actions associated with the current score event.

#### Returns

Event action list (empty if no score is loaded).

!!! note

    Delegates to the forward model when a score is active.

### `GetAudioStateChangeActions`

```cpp
EventActions OpenScofo::OpenScofo::GetAudioStateChangeActions()
```

## Audio Processing

### `ProcessBlock`

```cpp
template bool OpenScofo::OpenScofo::ProcessBlock< double >(const T *AudioBuffer, size_t n)
```

Process an incoming audio block.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `AudioBuffer` | `const T *` | Input audio buffer |
| `n` | `size_t` | Number of samples in buffer |

#### Returns

true if processing succeeded

!!! note

    Maintains internal circular buffer state.

!!! note

    Triggers analysis every hop size.

!!! note

    Updates descriptors and score position depending on mode.

### `GetBlockDuration`

```cpp
double OpenScofo::OpenScofo::GetBlockDuration()
```

Get processing block duration in seconds.

#### Returns

Block duration in seconds (derived from hop size and sampling rate)

### `GetCurrentBufferIndex`

```cpp
int OpenScofo::OpenScofo::GetCurrentBufferIndex()
```

Get current processing buffer index.

#### Returns

Current index within the analysis buffer (forward model state)

## Audio Descriptors

### `SetRequestedDescriptors`

```cpp
void OpenScofo::OpenScofo::SetRequestedDescriptors(std::vector< Descriptors > Descriptors)
```

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Descriptors` | `std::vector< Descriptors >` |  |

### `RequestDescriptor`

```cpp
void OpenScofo::OpenScofo::RequestDescriptor(Descriptors Descriptor)
```

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Descriptor` | `Descriptors` |  |

### `ActivateAllDescriptors`

```cpp
void OpenScofo::OpenScofo::ActivateAllDescriptors()
```

### `GetPitchProb`

```cpp
double OpenScofo::OpenScofo::GetPitchProb(double Freq)
```

Compute pitch probability for a given frequency.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Freq` | `double` | Frequency in Hz |

#### Returns

Probability score from forward model

!!! note

    Uses current description frame as input to the model.

### `GetPitchTemplate`

```cpp
std::vector< double > OpenScofo::OpenScofo::GetPitchTemplate(double Freq)
```

Generate pitch template for a given frequency.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Freq` | `double` | Target frequency in Hz |

#### Returns

Pitch template vector computed by the forward model

### `GetDescription`

```cpp
Description OpenScofo::OpenScofo::GetDescription()
```

Get current audio description frame.

#### Returns

Copy of the current descriptor structure

!!! note

    Returns by value (snapshot, not live reference).

### `GetDescriptorsEnum`

```cpp
Descriptors OpenScofo::OpenScofo::GetDescriptorsEnum(const char *s)
```

Convert string identifier to descriptor enum.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `s` | `const char *` | Descriptor name (e.g. "mfcc", "chroma") |

#### Returns

Corresponding Descriptors enum value, or INVALID on failure

!!! note

    Logs an error if the string is not recognized.

### `GetDescriptionId`

```cpp
const char * OpenScofo::OpenScofo::GetDescriptionId(Descriptors d)
```

Convert descriptor enum to string identifier.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `d` | `Descriptors` | Descriptor enum value |

#### Returns

Human-readable identifier string (e.g. "mfcc", "chroma")

### `GetDescriptionFloat`

```cpp
double OpenScofo::OpenScofo::GetDescriptionFloat(Description &Desc, Descriptors d)
```

Extract scalar descriptor value from a Description.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Desc` | `Description &` | Audio description container |
| `d` | `Descriptors` | Descriptor type |

#### Returns

Scalar value for the requested descriptor, or -1.0 on error / invalid type

!!! note

    Logs an error if the descriptor is vector-valued or invalid.

### `GetDescriptionArray`

```cpp
std::vector< double > & OpenScofo::OpenScofo::GetDescriptionArray(Description &Desc, Descriptors d)
```

Access vector-valued descriptor data.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Desc` | `Description &` | Audio description container |
| `d` | `Descriptors` | Descriptor type (must be vector-valued) |

#### Returns

Reference to internal descriptor array

!!! note

    Logs an error before throwing for invalid descriptor types.

## ONNX

### `LoadONNXModel`

```cpp
void OpenScofo::OpenScofo::LoadONNXModel(fs::path Model, std::vector< Descriptors > Descriptors)
```

Load an ONNX model for descriptor inference.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `Model` | `fs::path` | Path to .onnx model file |
| `Descriptors` | `std::vector< Descriptors >` | List of descriptors expected by the model |

!!! note

    Only .onnx models are supported.

!!! note

    Delegates initialization to the MIR module.

!!! warning

    Invalid file extension or incompatible descriptors will result in an error log.

## Lua

### `GetLuaCode`

```cpp
std::string OpenScofo::OpenScofo::GetLuaCode()
```

Return Lua code string defined in global events using LUA {}.

#### Returns

Lua script as a string

## Logging and Errors

### `ClearErrors`

```cpp
void OpenScofo::OpenScofo::ClearErrors()
```

Reset internal error state.

!!! note

    Clears m_HasErrors unless a critical error was previously set.

!!! warning

    If a critical error occurred, the state is not reset and recovery is not possible without reinitializing the instance.

### `SetErrorCallback`

```cpp
void OpenScofo::OpenScofo::SetErrorCallback(std::function< void(const spdlog::details::log_msg &, void *data)> cb, void *data=nullptr)
```

Set callback for log/error messages.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `cb` | `std::function< void(const spdlog::details::log_msg &, void *data)>` | Callback invoked on each log message |
| `data` | `void *` | User-defined pointer passed to the callback |

!!! note

    The callback is triggered by the internal logging sink.

!!! note

    Updates the internal error flag (m_HasErrors) automatically.

!!! note

    Log level is set to debug (debug builds) or info (release builds).

!!! warning

    Overwrites any previously registered callback.

### `SetLogLevel`

```cpp
void OpenScofo::OpenScofo::SetLogLevel(spdlog::level::level_enum level)
```

Set logging verbosity level.

#### Parameters

| Parameter | Type | Description |
| --- | --- | --- |
| `level` | `spdlog::level::level_enum` | spdlog log level |

!!! note

    Affects the global default spdlog logger.

## Other

### `GetStates`

```cpp
States & OpenScofo::OpenScofo::GetStates()
```

Access internal score state machine.

#### Returns

Reference to forward model states container.

!!! note

    Exposes internal mutable state (no copy is made).
