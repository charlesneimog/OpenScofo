/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#pragma once

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32) && !defined(OPENSCOFO_CATBOOST_STATIC)
#if defined(OPENSCOFO_CATBOOST_BRIDGE_BUILD)
#define OPENSCOFO_CATBOOST_API __declspec(dllexport)
#else
#define OPENSCOFO_CATBOOST_API __declspec(dllimport)
#endif
#elif defined(_WIN32)
#define OPENSCOFO_CATBOOST_API
#else
#define OPENSCOFO_CATBOOST_API __attribute__((visibility("default")))
#endif

#if defined(__cplusplus)
extern "C" {
#endif

typedef struct OpenScofoCatBoostHandle OpenScofoCatBoostHandle;

typedef struct OpenScofoCatBoostData {
    const float *Features;
    const int32_t *Labels;
    size_t SampleCount;
    size_t FeatureCount;
} OpenScofoCatBoostData;

typedef int (*OpenScofoCatBoostProgress)(void *userData, size_t iteration, size_t totalIterations);

OPENSCOFO_CATBOOST_API int OpenScofoCatBoostTrain(const OpenScofoCatBoostData *train,
                                                  const OpenScofoCatBoostData *validation,
                                                  const char *const *labelNames, size_t labelCount,
                                                  const char *parameters, size_t iterations,
                                                  OpenScofoCatBoostProgress progress, void *userData,
                                                  OpenScofoCatBoostHandle **result, char *error, size_t errorSize);

OPENSCOFO_CATBOOST_API int OpenScofoCatBoostPredict(const OpenScofoCatBoostHandle *handle, const float *features,
                                                    size_t sampleCount, size_t featureCount, int32_t *predictions,
                                                    size_t *outputDimension, char *error, size_t errorSize);

OPENSCOFO_CATBOOST_API int OpenScofoCatBoostExport(OpenScofoCatBoostHandle *handle, const char *path, int onnx,
                                                   const char *const *metadataKeys, const char *const *metadataValues,
                                                   size_t metadataCount, char *error, size_t errorSize);

OPENSCOFO_CATBOOST_API void OpenScofoCatBoostDestroy(OpenScofoCatBoostHandle *handle);

#if defined(__cplusplus)
}
#endif
