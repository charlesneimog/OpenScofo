# Running the tests

Tests are enabled by default (`OPENSCOFO_BUILD_TESTS=ON`). Normal builds also run
the `check` target (`OPENSCOFO_RUN_TESTS_ON_BUILD=ON`); a failure fails the build.
CTest is enabled at the project root on Linux, macOS, Windows, and Emscripten.

For a core-only build, without the host application wrappers:

```sh
cmake -S . -B build-tests -DOPENSCOFO_BUILD_ALL=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build-tests --config Release --parallel
```

Use the project's usual compiler and dependency setup. GoogleTest is found on
native Linux or built from pinned source when needed. macOS universal binaries,
Windows, and WebAssembly build GoogleTest with the same toolchain as OpenScofo.
The normal grammar-generation setup requires the Tree-sitter CLI.

To build and run just the tests, or rerun already-built tests:

```sh
cmake --build build-tests --target check --config Release --parallel
ctest --test-dir build-tests -C Release --output-on-failure --no-tests=error -L openscofo
```

To run only the configuration checks:

```sh
ctest --test-dir build-tests -C Release -R '^Configuration$' --output-on-failure --no-tests=error
```

`--config` and `-C` select the configuration for multi-configuration generators
such as Visual Studio and Xcode. Explicitly building an individual library or
wrapper target does not run the `ALL` target; use `check` afterward in that case.

To compile tests without running them during every build, configure with
`-DOPENSCOFO_RUN_TESTS_ON_BUILD=OFF`. The explicit `check` target still works.
`-DOPENSCOFO_BUILD_TESTS=OFF` removes the tests entirely. Performance executables
are separately enabled with `-DOPENSCOFO_BUILD_BENCHMARKS=ON` and require AudioFile.

## Recorded event timing benchmarks

The `benchmark_tests` executable reports the percentage of annotated
score events detected within **±250 ms** of their reference timestamps:

**match percentage = 100 × matched events / annotated events**

Missing events and detections outside tolerance reduce the percentage. Normal
native builds run the benchmarks through the existing `check` target, even when
no files changed, when `OPENSCOFO_RUN_TESTS_ON_BUILD=ON`.
Each piece has its own GoogleTest `TEST(Inference, ...)` case registered directly
with CTest under the `openscofo` and `benchmark` labels. Each case prints one
percentage alongside normal GoogleTest output. Use `ctest -V` to see percentages
for passing tests.

Build failure based on match quality is optional and disabled by default:

```cmake
OPENSCOFO_BENCHMARK_FAIL_BELOW_MINIMUM=OFF
OPENSCOFO_BENCHMARK_MIN_MINIATURA1=14.81
OPENSCOFO_BENCHMARK_MIN_MINIATURA2=100.00
OPENSCOFO_BENCHMARK_MIN_MINIATURA3=91.94
OPENSCOFO_BENCHMARK_MIN_CANTICOS=97.56
```

These are CMake cache options. Set them with `-D`, without a `%` suffix:

```sh
cmake -S . -B build-tests -DOPENSCOFO_BENCHMARK_FAIL_BELOW_MINIMUM=ON -DOPENSCOFO_BENCHMARK_MIN_MINIATURA1=90
cmake --build build-tests --parallel
```

With the toggle on, a percentage below that performance's minimum fails the
corresponding GoogleTest case and the build. Equality passes. The comparison uses
the displayed percentage rounded to two decimal places, so a displayed 91.94%
meets a 91.94 minimum. With the toggle off, low percentages never fail the build.

Install the native development packages for **nlohmann_json** and **libsndfile**
(1.1 or newer with MP3 support), then build and run:

```sh
cmake -S . -B build-tests -DOPENSCOFO_BUILD_ALL=OFF -DOPENSCOFO_BUILD_BENCHMARK_TESTS=ON
cmake --build build-tests --target benchmark_tests --parallel
ctest --test-dir build-tests -L benchmark -V
```

To report just one performance:

```sh
ctest --test-dir build-tests -R '^Benchmark.Miniatura1$' -V
```

The GoogleTest cases are `Inference.Miniatura1`, `Inference.Miniatura2`,
`Inference.Miniatura3`, and `Inference.Canticos`. For full GoogleTest diagnostics,
run the executable directly with `--gtest_filter=Inference.Miniatura1`.

Detailed results remain available in
`build-tests/Tests/benchmark-results/<piece>.json`. Each report contains
`summary.match_percentage`, expected and detected times, signed errors in
milliseconds, missing events, and unexpected event IDs. Negative errors mean
early detection; positive errors mean late detection. Missing detections have
null times and errors. Unexpected IDs are reported separately and do not count
as matches. The mean absolute error includes all detected reference events,
including those outside tolerance.

The runner uses `event` and `timestamp_seconds` from each annotation JSON and
reads the recording named by its `audio` field. Miniatura recordings belong in
`Tests/04-miniaturas/Audios/`; Cânticos audio belongs in
`Tests/01-benchmark/real/`. Recordings are not bundled in the repository. The
Miniatura scores also require the checked-in `Tests/flute.onnx`. Missing assets,
invalid annotations, decoding failures, and processing errors produce an error
instead of a misleading percentage; only these execution problems cause a
nonzero exit status when the percentage gate is disabled.

Audio is decoded from WAV or MP3 and processed in 64-sample blocks using the
first channel. The first detection of each event is timestamped at the end of
the consumed block, with no latency compensation. The last partial block is
processed without adding silence. Comparisons use the recording timeline, not
wall-clock processing time or the score's ideal tempo. The score controls
sample rate, FFT size, and hop size; sample-rate mismatches are reported as
errors. These settings are included in each report.

The runner is enabled by default for native builds with `OPENSCOFO_BUILD_TESTS`
enabled. Existing CMake caches retain their previous setting; use
`-DOPENSCOFO_BUILD_BENCHMARK_TESTS=ON` to enable it in an older build directory.
Use `OFF` to opt out. Emscripten remains unsupported by this audio-file runner
and defaults to disabled. Its option is separate from
`OPENSCOFO_BUILD_BENCHMARKS`, which builds the performance/profiling tools.

## Configuration coverage

[`02-score/configuration.scofo`](02-score/configuration.scofo) overrides every
score-configurable default. The parser checks include BPM, transposition, and
time tolerance stored on score states, as well as all configuration fields.
Aliases such as `TIMBREMODEL` refer to the same setting; the fixture uses the
canonical `ONNXMODEL`. `RequestedDescriptors` is an API setting, not a score
keyword, and is activated explicitly by the tests.

[`05-tests/test_configuration.cpp`](05-tests/test_configuration.cpp) also loads
the score through `OpenScofo::LoadScore` and checks runtime behavior: analysis
dimensions and update timing, tuning/transposition, harmonic templates, phase
and tempo adaptation, onset initialization, MFCC/chroma, YIN, rolloff, ZCR,
section selection, notification receivers, and ONNX predictions. Existing section
tests check that section restriction limits inference.

The small checked-in `configuration.onnx` classifies RMS into `quiet` and `loud`.
It uses the score's descriptor list, so tests need no training data, Python, or
external model download. To regenerate it, install Python's `onnx` package and
run `python Tests/05-tests/generate_configuration_model.py` from the repository root.

Two passing **characterization tests document unimplemented runtime wiring**:

- `KnownLimitationDbThresholdDoesNotAffectSilenceProbability`: `DBTHRESHOLD` is
  parsed, but silence probability still uses the fixed loudness midpoint.
- `KnownLimitationTimeToleranceDoesNotAffectDurationDistribution`:
  `TIMETOLERANCE` is stored on events, but the duration calculation ignores it.

These tests do not certify that those settings work. When implementing them,
replace the characterization assertions with tests of their intended effects.

## CI and WebAssembly

The binary workflow enables tests for each Linux, macOS, Windows, and Emscripten
compilation and explicitly runs CTest before packaging/upload steps. A test
failure stops that build job. Its existing tag/manual triggers are unchanged.

Emscripten tests run through CMake's cross-compiling emulator (Node.js), with
`NODERAWFS` enabled only for test executables so they can read the same fixtures:

```sh
emcmake cmake -S . -B build-wasm-tests -DOPENSCOFO_BUILD_ALL=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build-wasm-tests --target check --parallel
```

For other cross-compilers, set `CMAKE_CROSSCOMPILING_EMULATOR` to a compatible
runner, or disable automatic execution and run the binaries on the target.
