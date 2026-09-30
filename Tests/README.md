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
