/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "Trainer.hpp"

#include "CatBoostBridge.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <stdexcept>
#include <string>
#include <vector>

namespace OpenScofo {

std::string Trainer::TrainingParameters(bool haveValidation) const {
    std::string json = "{";
    json += "\"iterations\":" + std::to_string(m_Config.Iterations);
    json += ",\"learning_rate\":" + std::to_string(m_Config.LearningRate);
    json += ",\"depth\":" + std::to_string(m_Config.Depth);
    json += ",\"loss_function\":\"MultiClass\"";
    json += ",\"random_seed\":" + std::to_string(m_Config.RandomSeed);
    // CatBoost's native C++ option parser only accepts positive thread counts,
    // unlike the Python API where -1 means automatic/all available threads.
    if (m_Config.ThreadCount > 0) {
        json += ",\"thread_count\":" + std::to_string(m_Config.ThreadCount);
    }
    json += ",\"allow_writing_files\":false";
    json += ",\"logging_level\":\"Silent\"";
    if (haveValidation) {
        json += ",\"use_best_model\":true";
        json += ",\"od_type\":\"Iter\"";
        json += ",\"od_wait\":" + std::to_string(m_Config.EarlyStoppingRounds);
    } else {
        json += ",\"use_best_model\":false";
    }
    json += "}";
    return json;
}

std::string Trainer::JsonString(const std::string &value) const {
    std::string output = "\"";
    for (unsigned char character : value) {
        if (character == '"' || character == '\\') {
            output.push_back('\\');
            output.push_back(static_cast<char>(character));
        } else if (character == '\n') {
            output += "\\n";
        } else if (character == '\r') {
            output += "\\r";
        } else if (character == '\t') {
            output += "\\t";
        } else if (character >= 0x20) {
            output.push_back(static_cast<char>(character));
        }
    }
    output.push_back('"');
    return output;
}

std::string Trainer::JsonArray(const std::vector<std::string> &values) const {
    std::string output = "[";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            output.push_back(',');
        }
        output += JsonString(values[i]);
    }
    output.push_back(']');
    return output;
}

void Trainer::TrainModel() {
    if (m_Data.Train.Empty()) {
        throw std::runtime_error("No training samples found");
    }
    if (m_Data.Labels.size() < 2) {
        throw std::runtime_error("Training requires at least two classes");
    }

    const OpenScofoCatBoostData train = {m_Data.Train.Features.data(), m_Data.Train.Labels.data(),
                                         m_Data.Train.SampleCount(), m_Data.Train.FeatureCount};
    const OpenScofoCatBoostData validation = {m_Data.Validation.Features.data(), m_Data.Validation.Labels.data(),
                                              m_Data.Validation.SampleCount(), m_Data.Validation.FeatureCount};
    const OpenScofoCatBoostData *validationPointer = m_Data.Validation.Empty() ? nullptr : &validation;
    const std::string parameters = TrainingParameters(validationPointer != nullptr);

    struct CallbackState {
        std::atomic<bool> *Cancelled = nullptr;
        const ProgressCallback *Progress = nullptr;
    };
    CallbackState callbackState{&m_Cancelled, &m_Progress};
    const auto afterIteration = +[](void *data, std::size_t iteration, std::size_t iterations) -> int {
        CallbackState *state = static_cast<CallbackState *>(data);
        if (state->Cancelled->load()) {
            return 0;
        }
        try {
            if (*state->Progress) {
                const double value = std::min(1.0, static_cast<double>(iteration) / static_cast<double>(iterations));
                (*state->Progress)("training", value);
            }
        } catch (...) {
            state->Cancelled->store(true);
            return 0;
        }
        return state->Cancelled->load() ? 0 : 1;
    };

    std::vector<const char *> labelNames;
    labelNames.reserve(m_Data.Labels.size());
    for (const std::string &label : m_Data.Labels) {
        labelNames.push_back(label.c_str());
    }

    std::array<char, 4096> error{};
    OpenScofoCatBoostHandle *model = nullptr;
    if (!OpenScofoCatBoostTrain(&train, validationPointer, labelNames.data(), labelNames.size(), parameters.c_str(),
                                static_cast<std::size_t>(std::max(1, m_Config.Iterations)), afterIteration,
                                &callbackState, &model, error.data(), error.size())) {
        throw std::runtime_error(error[0] == '\0' ? "CatBoost training failed" : error.data());
    }
    if (m_Cancelled.load()) {
        OpenScofoCatBoostDestroy(model);
        throw std::runtime_error("Training cancelled");
    }

    ClearModel();
    m_Model = model;
    if (m_Progress) {
        m_Progress("training", 1.0);
    }
}

Trainer::Evaluation Trainer::EvaluateModel() const {
    Evaluation evaluation;
    if (m_Model == nullptr || m_Data.Test.Empty()) {
        return evaluation;
    }

    std::vector<std::int32_t> predictions(m_Data.Test.SampleCount());
    std::size_t dimensions = 0;
    std::array<char, 4096> error{};
    if (!OpenScofoCatBoostPredict(m_Model, m_Data.Test.Features.data(), m_Data.Test.SampleCount(),
                                  m_Data.Test.FeatureCount, predictions.data(), &dimensions, error.data(),
                                  error.size())) {
        throw std::runtime_error(error[0] == '\0' ? "CatBoost prediction failed" : error.data());
    }
    if (dimensions != m_Data.Labels.size()) {
        throw std::runtime_error("CatBoost output dimension does not match the OpenScofo class list");
    }

    std::vector<std::size_t> totals(m_Data.Labels.size(), 0);
    std::vector<std::size_t> correct(m_Data.Labels.size(), 0);
    std::size_t allCorrect = 0;
    for (std::size_t sample = 0; sample < predictions.size(); ++sample) {
        const std::int32_t expectedValue = m_Data.Test.Labels[sample];
        if (expectedValue < 0 || static_cast<std::size_t>(expectedValue) >= m_Data.Labels.size()) {
            throw std::runtime_error("Test partition contains an invalid class index");
        }
        const std::size_t expected = static_cast<std::size_t>(expectedValue);
        ++totals[expected];
        if (predictions[sample] == expectedValue) {
            ++correct[expected];
            ++allCorrect;
        }
    }

    evaluation.TestedSamples = predictions.size();
    evaluation.Accuracy = static_cast<double>(allCorrect) / static_cast<double>(predictions.size());
    double balancedTotal = 0.0;
    std::size_t balancedClasses = 0;
    for (std::size_t label = 0; label < m_Data.Labels.size(); ++label) {
        if (totals[label] == 0) {
            continue;
        }
        const double accuracy = static_cast<double>(correct[label]) / static_cast<double>(totals[label]);
        evaluation.ClassAccuracy[m_Data.Labels[label]] = accuracy;
        balancedTotal += accuracy;
        ++balancedClasses;
    }
    if (balancedClasses != 0) {
        evaluation.BalancedAccuracy = balancedTotal / static_cast<double>(balancedClasses);
    }
    return evaluation;
}

void Trainer::ExportTrainedModel(const std::filesystem::path &path) const {
    if (m_Model == nullptr) {
        throw std::runtime_error("No trained model to export. Call Train() first");
    }
    if (path.empty()) {
        throw std::invalid_argument("Model export path is empty");
    }
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    const std::array<const char *, 6> keys = {"openscofo.descriptors", "openscofo.labels",   "openscofo.sample_rate",
                                              "openscofo.fft_size",    "openscofo.hop_size", "openscofo.version"};
    const std::array<std::string, 6> values = {
        JsonArray(m_Data.Descriptors),       JsonArray(m_Data.Labels),
        std::to_string(m_Config.SampleRate), std::to_string(m_Config.FFTSize),
        std::to_string(m_Config.HopSize),    OPENSCOFO_VERSION};
    std::array<const char *, 6> valuePointers{};
    for (std::size_t i = 0; i < values.size(); ++i) {
        valuePointers[i] = values[i].c_str();
    }

    std::string extension = path.extension().string();
    for (char &character : extension) {
        character = static_cast<char>(std::tolower(static_cast<unsigned char>(character)));
    }
    const int onnx = extension == ".onnx" ? 1 : 0;
    const std::string pathString = path.string();
    std::array<char, 4096> error{};
    if (!OpenScofoCatBoostExport(m_Model, pathString.c_str(), onnx, keys.data(), valuePointers.data(), keys.size(),
                                 error.data(), error.size())) {
        throw std::runtime_error(error[0] == '\0' ? "CatBoost model export failed" : error.data());
    }
}

void Trainer::ClearModel() noexcept {
    OpenScofoCatBoostDestroy(m_Model);
    m_Model = nullptr;
}

} // namespace OpenScofo
