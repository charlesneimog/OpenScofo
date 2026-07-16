/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#pragma once

#include "states.hpp"

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

struct onnx_context_t;
struct onnx_tensor_t;

namespace OpenScofo {

class ONNXModel {
  public:
    ONNXModel() = default;
    ~ONNXModel();

    ONNXModel(const ONNXModel &) = delete;
    ONNXModel &operator=(const ONNXModel &) = delete;
    ONNXModel(ONNXModel &&) = delete;
    ONNXModel &operator=(ONNXModel &&) = delete;

    bool Load(const std::filesystem::path &path, std::vector<Descriptors> descriptors,
              const Configuration &configuration);
    void Execute(Description &description);
    void Reset();

    bool IsLoaded() const;
    const std::vector<std::string> &GetLabels() const;
    const std::vector<Descriptors> &GetDescriptors() const;

  private:
    void ReadMetadata(std::vector<Descriptors> &descriptors, const Configuration &configuration,
                      std::vector<std::string> &metadataLabels);
    bool ValidateOpsets() const;
    bool ReadLabelsFromModel(const std::vector<std::string> &metadataLabels);
    bool PrepareDescriptorBuffer(const Configuration &configuration);
    bool FindTensors();
    void WriteDescriptor(Descriptors descriptor, const Description &description, float *&output) const;

    static void SkipWhitespace(std::string_view text, std::size_t &position);
    static bool ParseJsonString(std::string_view text, std::size_t &position, std::string &output);
    static bool ParseJsonStringArray(const char *metadata, std::vector<std::string> &output);
    static bool ParseIntMetadata(const char *metadata, int &output);
    static Descriptors DescriptorFromMetadataName(std::string_view name);
    static bool ParseDescriptorMetadata(const char *metadata, std::vector<Descriptors> &output);
    static bool IsFloatTensor(const struct onnx_tensor_t *tensor);

    static const int CurrentOpset;

    bool m_ModelLoaded = false;
    std::filesystem::path m_ModelPath;
    struct onnx_context_t *m_Context = nullptr;
    struct onnx_tensor_t *m_InputTensor = nullptr;
    struct onnx_tensor_t *m_OutputTensor = nullptr;

    int m_MFCCCount = 0;
    int m_MFCCMels = 0;
    int m_ChromaSize = 0;
    int m_FFTSize = 0;
    std::size_t m_DescriptorSize = 0;

    std::vector<std::string> m_Labels;
    std::vector<Descriptors> m_Descriptors;
    std::vector<float> m_DescriptorValues;
};

} // namespace OpenScofo
