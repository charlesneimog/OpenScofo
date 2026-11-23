/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "CatBoostBridge.h"

#include <algorithm>
#include <cstdio>

#include <catboost/libs/data/data_provider_builders.h>
#include <catboost/libs/helpers/json_helpers.h>
#include <catboost/libs/model/model.h>
#include <catboost/libs/model/model_export/model_exporter.h>
#include <catboost/libs/train_lib/train_model.h>
#include <library/cpp/json/json_reader.h>
#include <library/cpp/threading/local_executor/local_executor.h>
#include <util/generic/array_ref.h>
#include <util/generic/ptr.h>
#include <util/generic/string.h>
#include <util/generic/vector.h>
#include <util/generic/yexception.h>

// ─────────────────────────────────────
struct OpenScofoCatBoostHandle {
    TFullModel *Model = nullptr;
};

// ─────────────────────────────────────
struct CallbackState {
    OpenScofoCatBoostProgress Progress = nullptr;
    void *UserData = nullptr;
    std::size_t Iterations = 1;
};

// ─────────────────────────────────────
void SetError(char *error, std::size_t errorSize, const TString &message) {
    if (error == nullptr || errorSize == 0) {
        return;
    }
    std::snprintf(error, errorSize, "%s", message.c_str());
}

// ─────────────────────────────────────
void ClearError(char *error, std::size_t errorSize) {
    if (error != nullptr && errorSize != 0) {
        error[0] = '\0';
    }
}

// ─────────────────────────────────────
bool AfterIteration(const TMetricsAndTimeLeftHistory &history, void *data) {
    CallbackState *state = static_cast<CallbackState *>(data);
    if (state->Progress == nullptr) {
        return true;
    }
    return state->Progress(state->UserData, history.LearnMetricsHistory.size(), state->Iterations) != 0;
}

// ─────────────────────────────────────
NCB::TDataProviderPtr MakeProvider(const OpenScofoCatBoostData &data, const char *const *labelNames,
                                   std::size_t labelCount, TVector<float> &features, TVector<TString> &targets) {
    if (data.Features == nullptr || data.Labels == nullptr || data.SampleCount == 0 || data.FeatureCount == 0) {
        ythrow yexception() << "CatBoost received an empty or invalid dataset";
    }

    features.resize(data.SampleCount * data.FeatureCount);
    for (std::size_t sample = 0; sample < data.SampleCount; ++sample) {
        for (std::size_t feature = 0; feature < data.FeatureCount; ++feature) {
            features[feature * data.SampleCount + sample] = data.Features[sample * data.FeatureCount + feature];
        }
    }

    targets.resize(data.SampleCount);
    for (std::size_t sample = 0; sample < data.SampleCount; ++sample) {
        const int32_t label = data.Labels[sample];
        if (label < 0 || static_cast<std::size_t>(label) >= labelCount || labelNames[label] == nullptr) {
            ythrow yexception() << "CatBoost received an invalid class index";
        }
        targets[sample] = labelNames[label];
    }

    NCB::TDataProviderBuilderOptions options;
    THolder<NCB::IDataProviderBuilder> builder;
    NCB::IRawFeaturesOrderDataVisitor *visitor = nullptr;
    NCB::CreateDataProviderBuilderAndVisitor(options, &NPar::LocalExecutor(), &builder, &visitor);

    NCB::TDataMetaInfo metaInfo;
    metaInfo.TargetType = NCB::ERawTargetType::String;
    metaInfo.TargetCount = 1;
    metaInfo.FeaturesLayout = new NCB::TFeaturesLayout(static_cast<ui32>(data.FeatureCount));
    visitor->Start(metaInfo, static_cast<ui32>(data.SampleCount), NCB::EObjectsOrder::Ordered, {});

    for (std::size_t feature = 0; feature < data.FeatureCount; ++feature) {
        const TConstArrayRef<float> values(features.data() + feature * data.SampleCount, data.SampleCount);
        auto holder = MakeIntrusive<NCB::TTypeCastArrayHolder<float, float>>(
            NCB::TMaybeOwningConstArrayHolder<float>::CreateNonOwning(values));
        visitor->AddFloatFeature(static_cast<ui32>(feature), holder);
    }
    visitor->AddTarget(TConstArrayRef<TString>(targets.data(), targets.size()));
    visitor->Finish();
    return builder->GetResult();
}

// ─────────────────────────────────────
int OpenScofoCatBoostTrain(const OpenScofoCatBoostData *train, const OpenScofoCatBoostData *validation,
                           const char *const *labelNames, std::size_t labelCount, const char *parameters,
                           std::size_t iterations, OpenScofoCatBoostProgress progress, void *userData,
                           OpenScofoCatBoostHandle **result, char *error, std::size_t errorSize) {
    ClearError(error, errorSize);
    if (train == nullptr || labelNames == nullptr || labelCount < 2 || parameters == nullptr || result == nullptr) {
        SetError(error, errorSize, "CatBoost training arguments are invalid");
        return 0;
    }
    *result = nullptr;

    try {
        TVector<float> trainFeatures;
        TVector<float> validationFeatures;
        TVector<TString> trainTargets;
        TVector<TString> validationTargets;

        NCB::TDataProviders providers;
        providers.Learn = MakeProvider(*train, labelNames, labelCount, trainFeatures, trainTargets);
        if (validation != nullptr && validation->SampleCount != 0) {
            if (validation->FeatureCount != train->FeatureCount) {
                ythrow yexception() << "Training and validation feature dimensions do not match";
            }
            providers.Test.push_back(
                MakeProvider(*validation, labelNames, labelCount, validationFeatures, validationTargets));
        }

        NJson::TJsonValue jsonParameters;
        NJson::ReadJsonTree(TString(parameters), &jsonParameters);

        CallbackState callbackState;
        callbackState.Progress = progress;
        callbackState.UserData = userData;
        callbackState.Iterations = std::max<std::size_t>(1, iterations);
        TCustomCallbackDescriptor callbackDescriptor;
        callbackDescriptor.CustomData = &callbackState;
        callbackDescriptor.AfterIterationFunc = AfterIteration;

        THolder<TFullModel> model = MakeHolder<TFullModel>();
        TMetricsAndTimeLeftHistory history;
        TEvalResult validationResult;
        TVector<TEvalResult *> evaluationResults;
        if (validation != nullptr && validation->SampleCount != 0) {
            evaluationResults.push_back(&validationResult);
        }

        NCB::TQuantizedFeaturesInfoPtr quantizedFeaturesInfo;
        TMaybe<TCustomObjectiveDescriptor> objectiveDescriptor;
        TMaybe<TCustomMetricDescriptor> metricDescriptor;
        TrainModel(jsonParameters, quantizedFeaturesInfo, objectiveDescriptor, metricDescriptor, callbackDescriptor,
                   std::move(providers), Nothing(), nullptr, TString(), model.Get(), evaluationResults, &history);

        THolder<OpenScofoCatBoostHandle> output = MakeHolder<OpenScofoCatBoostHandle>();
        output->Model = model.Release();
        *result = output.Release();
        return 1;
    } catch (...) {
        SetError(error, errorSize, CurrentExceptionMessage());
        return 0;
    }
}

// ─────────────────────────────────────
int OpenScofoCatBoostPredict(const OpenScofoCatBoostHandle *handle, const float *features, std::size_t sampleCount,
                             std::size_t featureCount, int32_t *predictions, std::size_t *outputDimension, char *error,
                             std::size_t errorSize) {
    ClearError(error, errorSize);
    if (handle == nullptr || handle->Model == nullptr || features == nullptr || predictions == nullptr ||
        outputDimension == nullptr || sampleCount == 0 || featureCount == 0) {
        SetError(error, errorSize, "CatBoost prediction arguments are invalid");
        return 0;
    }

    try {
        TVector<TVector<float>> rows(sampleCount);
        for (std::size_t sample = 0; sample < sampleCount; ++sample) {
            rows[sample].assign(features + sample * featureCount, features + (sample + 1) * featureCount);
        }

        const std::size_t dimensions = handle->Model->GetDimensionsCount();
        TVector<double> values(sampleCount * dimensions);
        handle->Model->CalcFlat(rows, values);
        for (std::size_t sample = 0; sample < sampleCount; ++sample) {
            const double *begin = values.data() + sample * dimensions;
            predictions[sample] = static_cast<int32_t>(std::max_element(begin, begin + dimensions) - begin);
        }
        *outputDimension = dimensions;
        return 1;
    } catch (...) {
        SetError(error, errorSize, CurrentExceptionMessage());
        return 0;
    }
}

// ─────────────────────────────────────
int OpenScofoCatBoostExport(OpenScofoCatBoostHandle *handle, const char *path, int onnx,
                            const char *const *metadataKeys, const char *const *metadataValues,
                            std::size_t metadataCount, char *error, std::size_t errorSize) {
    ClearError(error, errorSize);
    if (handle == nullptr || handle->Model == nullptr || path == nullptr || path[0] == '\0') {
        SetError(error, errorSize, "CatBoost export arguments are invalid");
        return 0;
    }

    try {
        for (std::size_t i = 0; i < metadataCount; ++i) {
            if (metadataKeys == nullptr || metadataValues == nullptr || metadataKeys[i] == nullptr ||
                metadataValues[i] == nullptr) {
                ythrow yexception() << "CatBoost metadata arguments are invalid";
            }
            handle->Model->ModelInfo[metadataKeys[i]] = metadataValues[i];
        }
        NCB::ExportModel(*handle->Model, TString(path), onnx != 0 ? EModelType::Onnx : EModelType::CatboostBinary);
        return 1;
    } catch (...) {
        SetError(error, errorSize, CurrentExceptionMessage());
        return 0;
    }
}

// ─────────────────────────────────────
void OpenScofoCatBoostDestroy(OpenScofoCatBoostHandle *handle) {
    if (handle == nullptr) {
        return;
    }
    delete handle->Model;
    delete handle;
}
