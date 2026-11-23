/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "Trainer.hpp"

#include "CatBoostBridge.h"

#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <system_error>
#include <vector>

namespace OpenScofo {

// ─────────────────────────────────────
std::string Trainer::TrainingParameters(bool haveValidation) const {
    nlohmann::json parameters = {
        {"iterations", m_Config.Iterations}, {"learning_rate", m_Config.LearningRate}, {"depth", m_Config.Depth},
        {"loss_function", "MultiClass"},     {"random_seed", m_Config.RandomSeed},     {"allow_writing_files", false},
        {"logging_level", "Silent"},         {"use_best_model", haveValidation},
    };
    // CatBoost's native C++ option parser only accepts positive thread counts,
    // unlike the Python API where -1 means automatic/all available threads.
    if (m_Config.ThreadCount > 0) {
        parameters["thread_count"] = m_Config.ThreadCount;
    }
    if (haveValidation) {
        parameters["od_type"] = "Iter";
        parameters["od_wait"] = m_Config.EarlyStoppingRounds;
    }
    return parameters.dump();
}

// ─────────────────────────────────────
// Called by the CatBoost bridge after every boosting iteration. Exceptions must not cross the C bridge, so a failing
// progress callback cancels the training instead of propagating.
int Trainer::TrainingProgress(void *data, std::size_t iteration, std::size_t iterations) {
    Trainer *trainer = static_cast<Trainer *>(data);
    if (trainer->m_Cancelled.load()) {
        return 0;
    }
    try {
        trainer->Report("training", static_cast<double>(iteration) / static_cast<double>(iterations));
    } catch (...) {
        trainer->m_Logger->error("Progress callback failed; cancelling training");
        trainer->m_Cancelled.store(true);
        return 0;
    }
    return trainer->m_Cancelled.load() ? 0 : 1;
}

// ─────────────────────────────────────
bool Trainer::TrainModel() {
    const OpenScofoCatBoostData train = {m_Data.Train.Features.data(), m_Data.Train.Labels.data(),
                                         m_Data.Train.SampleCount(), m_Data.Train.FeatureCount};
    const OpenScofoCatBoostData validation = {m_Data.Validation.Features.data(), m_Data.Validation.Labels.data(),
                                              m_Data.Validation.SampleCount(), m_Data.Validation.FeatureCount};
    const OpenScofoCatBoostData *validationPointer = m_Data.Validation.Empty() ? nullptr : &validation;
    const std::string parameters = TrainingParameters(validationPointer != nullptr);

    std::vector<const char *> labelNames;
    labelNames.reserve(m_Data.Labels.size());
    for (const std::string &label : m_Data.Labels) {
        labelNames.push_back(label.c_str());
    }

    std::array<char, 4096> error{};
    OpenScofoCatBoostHandle *model = nullptr;
    if (!OpenScofoCatBoostTrain(&train, validationPointer, labelNames.data(), labelNames.size(), parameters.c_str(),
                                static_cast<std::size_t>(m_Config.Iterations), TrainingProgress, this, &model,
                                error.data(), error.size())) {
        m_Logger->error("{}", error[0] == '\0' ? "CatBoost training failed" : error.data());
        return false;
    }
    if (m_Cancelled.load()) {
        OpenScofoCatBoostDestroy(model);
        return false;
    }

    ClearModel();
    m_Model = model;
    Report("training", 1.0);
    return true;
}

// ─────────────────────────────────────
bool Trainer::EvaluateModel(Evaluation &evaluation) const {
    evaluation = {};
    if (m_Data.Test.Empty()) {
        return true;
    }

    std::vector<std::int32_t> predictions(m_Data.Test.SampleCount());
    std::size_t dimensions = 0;
    std::array<char, 4096> error{};
    if (!OpenScofoCatBoostPredict(m_Model, m_Data.Test.Features.data(), m_Data.Test.SampleCount(),
                                  m_Data.Test.FeatureCount, predictions.data(), &dimensions, error.data(),
                                  error.size())) {
        m_Logger->error("{}", error[0] == '\0' ? "CatBoost prediction failed" : error.data());
        return false;
    }
    if (dimensions != m_Data.Labels.size()) {
        m_Logger->error("CatBoost output dimension does not match the OpenScofo class list");
        return false;
    }

    std::vector<std::size_t> totals(m_Data.Labels.size(), 0);
    std::vector<std::size_t> correct(m_Data.Labels.size(), 0);
    std::size_t allCorrect = 0;
    for (std::size_t sample = 0; sample < predictions.size(); ++sample) {
        const std::int32_t expected = m_Data.Test.Labels[sample];
        ++totals[static_cast<std::size_t>(expected)];
        if (predictions[sample] == expected) {
            ++correct[static_cast<std::size_t>(expected)];
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
    return true;
}

// ─────────────────────────────────────
bool Trainer::ExportTrainedModel(const std::filesystem::path &path) const {
    if (m_Model == nullptr) {
        m_Logger->error("No trained model to export. Call Train() first");
        return false;
    }
    if (path.empty()) {
        m_Logger->error("Model export path is empty");
        return false;
    }
    if (path.has_parent_path()) {
        std::error_code errorCode;
        std::filesystem::create_directories(path.parent_path(), errorCode);
        if (errorCode) {
            m_Logger->error("Could not create model export directory '{}': {}", path.parent_path().string(),
                            errorCode.message());
            return false;
        }
    }

    // CatBoost stores model metadata as string pairs. Descriptors and labels are JSON arrays kept inside those strings.
    // Class folder names are not guaranteed to be UTF-8 (e.g. Windows code pages), so invalid bytes are replaced rather
    // instead of rejecting the metadata after training has finished.
    constexpr nlohmann::json::error_handler_t replaceInvalidUtf8 = nlohmann::json::error_handler_t::replace;
    const std::map<std::string, std::string> metadata = {
        {"openscofo.descriptors", nlohmann::json(m_Data.Descriptors).dump(-1, ' ', false, replaceInvalidUtf8)},
        {"openscofo.labels", nlohmann::json(m_Data.Labels).dump(-1, ' ', false, replaceInvalidUtf8)},
        {"openscofo.sample_rate", std::to_string(m_Config.SampleRate)},
        {"openscofo.fft_size", std::to_string(m_Config.FFTSize)},
        {"openscofo.hop_size", std::to_string(m_Config.HopSize)},
        {"openscofo.version", OPENSCOFO_VERSION},
    };
    std::vector<const char *> keys;
    std::vector<const char *> values;
    for (const auto &[key, value] : metadata) {
        keys.push_back(key.c_str());
        values.push_back(value.c_str());
    }

    const int onnx = Lowercase(path.extension().string()) == ".onnx" ? 1 : 0;
    const std::string pathString = path.string();
    std::array<char, 4096> error{};
    if (!OpenScofoCatBoostExport(m_Model, pathString.c_str(), onnx, keys.data(), values.data(), keys.size(),
                                 error.data(), error.size())) {
        m_Logger->error("{}", error[0] == '\0' ? "CatBoost model export failed" : error.data());
        return false;
    }
    return true;
}

// ─────────────────────────────────────
void Trainer::ClearModel() noexcept {
    OpenScofoCatBoostDestroy(m_Model);
    m_Model = nullptr;
}

} // namespace OpenScofo
