/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file onnx.cpp
 * @brief ONNX classifier loading, metadata parsing, and descriptor inference.
 *
 * @note Uses ordered descriptor buffers and string labels to map classifier output into audio descriptions.
 * @warning Model loading, reset, and inference share owned buffers and must be serialized.
 */

#include <onnx.h>

#include "onnx.hpp"
#include "log.hpp"

#include <charconv>
#include <cctype>
#include <cstring>
#include <system_error>
#include <utility>

namespace OpenScofo {

const int ONNXModel::CurrentOpset = 24;

// ─────────────────────────────────────
/**
 * @brief Release the owned ONNX context.
 *
 * @note Calls Reset() to clear model state and inference buffers.
 */
ONNXModel::~ONNXModel() {
    Reset();
}

// ─────────────────────────────────────
/**
 * @brief Load and validate a descriptor-based ONNX classifier.
 *
 * @param path Path to the ONNX model file.
 * @param descriptors Ordered model input descriptors; metadata may replace them.
 * @param configuration Audio and descriptor dimensions used to prepare and validate model inputs.
 *
 * @return True if the model is ready or the same path is already loaded; false on validation or loading failure.
 *
 * @note Reads metadata, checks opsets and labels, sizes input storage, and selects compatible tensors.
 * @warning A different path resets the previous model before validation. Reloading the same path skips
 * configuration changes.
 */
bool ONNXModel::Load(const std::filesystem::path &path, std::vector<Descriptors> descriptors,
                     const Configuration &configuration) {
    if (m_ModelLoaded && m_ModelPath == path) {
        spdlog::warn("ONNX model already loaded");
        return true;
    }

    Reset();
    m_ModelPath = path;

    const std::u8string utf8Path = path.u8string();
    const std::string modelPath(utf8Path.begin(), utf8Path.end());
    m_Context = onnx_context_alloc_from_file(modelPath.c_str(), nullptr, 0);
    if (m_Context == nullptr) {
        spdlog::error("Failed to load ONNX model: {}.", path.string());
        return false;
    }

    std::vector<std::string> metadataLabels;
    ReadMetadata(descriptors, configuration, metadataLabels);

    if (!ValidateOpsets() || !ReadLabelsFromModel(metadataLabels)) {
        return false;
    }

    m_Descriptors = std::move(descriptors);
    if (!PrepareDescriptorBuffer(configuration) || !FindTensors()) {
        return false;
    }

    m_ModelLoaded = true;
    spdlog::debug("ONNX model loaded");
    return true;
}

// ─────────────────────────────────────
/**
 * @brief Run classifier inference on an audio description.
 *
 * @param description Audio frame containing input descriptors and receiving label evidence.
 *
 * @note Flattens descriptors in model order and maps numeric output values to labels; clears ONNX results if
 * unavailable.
 * @warning Descriptor arrays must match the dimensions recorded at load time; serialize with Load() and Reset().
 */
void ONNXModel::Execute(Description &description) {
    if (!m_ModelLoaded || m_InputTensor == nullptr || m_OutputTensor == nullptr) {
        description.ONNX.clear();
        return;
    }

    float *output = m_DescriptorValues.data();
    for (const Descriptors descriptor : m_Descriptors) {
        WriteDescriptor(descriptor, description, output);
    }

    if (output != m_DescriptorValues.data() + m_DescriptorSize) {
        spdlog::error("ONNX descriptor writer produced {} values, expected {}", output - m_DescriptorValues.data(),
                      m_DescriptorSize);
        description.ONNX.clear();
        return;
    }

    onnx_tensor_apply(m_InputTensor, m_DescriptorValues.data(),
                      m_DescriptorValues.size() * sizeof(m_DescriptorValues[0]));
    onnx_run(m_Context);

    description.ONNX.clear();
    if (m_OutputTensor->type == ONNX_TENSOR_TYPE_FLOAT32) {
        const float *data = static_cast<const float *>(m_OutputTensor->datas);
        for (std::size_t i = 0; i < m_Labels.size(); ++i) {
            description.ONNX[m_Labels[i]] = data[i];
        }
        return;
    }

    if (m_OutputTensor->type == ONNX_TENSOR_TYPE_FLOAT64) {
        const double *data = static_cast<const double *>(m_OutputTensor->datas);
        for (std::size_t i = 0; i < m_Labels.size(); ++i) {
            description.ONNX[m_Labels[i]] = static_cast<float>(data[i]);
        }
        return;
    }

    spdlog::error("Tensor output type not supported {}", static_cast<int>(m_OutputTensor->type));
}

// ─────────────────────────────────────
/**
 * @brief Release the ONNX context and clear model configuration.
 *
 * @note Clears tensor pointers, descriptors, labels, model path, and flattened input storage.
 * @warning Invalidates tensor pointers and elements referenced through the model descriptor or label containers.
 */
void ONNXModel::Reset() {
    if (m_Context != nullptr) {
        onnx_context_free(m_Context);
    }

    m_ModelLoaded = false;
    m_ModelPath.clear();
    m_Context = nullptr;
    m_InputTensor = nullptr;
    m_OutputTensor = nullptr;
    m_DescriptorSize = 0;
    m_Labels.clear();
    m_Descriptors.clear();
    m_DescriptorValues.clear();
}

// ─────────────────────────────────────
/**
 * @brief Check whether model loading and validation completed.
 *
 * @return True if the model is ready; false otherwise.
 *
 * @note The flag is set only after descriptor storage and tensors are prepared.
 */
bool ONNXModel::IsLoaded() const {
    return m_ModelLoaded;
}

// ─────────────────────────────────────
/**
 * @brief Access the loaded classifier label list.
 *
 * @return Const reference to the class labels.
 *
 * @note Returns model-owned storage without copying.
 * @warning Reset() or loading a different model clears the contents and invalidates references to its elements.
 */
const std::vector<std::string> &ONNXModel::GetLabels() const {
    return m_Labels;
}

// ─────────────────────────────────────
/**
 * @brief Access the ordered model input descriptor list.
 *
 * @return Const reference to the model input descriptors.
 *
 * @note The order determines how Execute() flattens each frame.
 * @warning Reset() or loading a different model clears the contents and invalidates references to its elements.
 */
const std::vector<Descriptors> &ONNXModel::GetDescriptors() const {
    return m_Descriptors;
}

// ─────────────────────────────────────
/**
 * @brief Read OpenScofo training configuration, descriptor, and label metadata.
 *
 * @param descriptors Ordered model input descriptors; metadata may replace them.
 * @param configuration Audio and descriptor dimensions used to prepare and validate model inputs.
 * @param metadataLabels Output label list parsed from metadata when available.
 *
 * @note Valid descriptor metadata overrides caller inputs; mismatched training dimensions are logged as warnings.
 * @warning Requires an allocated model context; malformed metadata falls back to caller descriptors or classifier
 * labels.
 */
void ONNXModel::ReadMetadata(std::vector<Descriptors> &descriptors, const Configuration &configuration,
                             std::vector<std::string> &metadataLabels) {

    const char *sampleRateMetadata = onnx_metadata_get(m_Context, "openscofo.sample_rate");
    const char *fftSizeMetadata = onnx_metadata_get(m_Context, "openscofo.fft_size");
    const char *hopSizeMetadata = onnx_metadata_get(m_Context, "openscofo.hop_size");
    const char *descriptorMetadata = onnx_metadata_get(m_Context, "openscofo.descriptors");
    const char *labelMetadata = onnx_metadata_get(m_Context, "openscofo.labels");

    int metadataSampleRate = 0;
    if (ParseIntMetadata(sampleRateMetadata, metadataSampleRate) && metadataSampleRate != configuration.SR) {
        spdlog::warn("ONNX model was trained at {} Hz, but OpenScofo is running at {} Hz.", metadataSampleRate,
                     configuration.SR);
    }

    int metadataFFTSize = 0;
    if (ParseIntMetadata(fftSizeMetadata, metadataFFTSize) && metadataFFTSize != configuration.FFTSize) {
        spdlog::warn("ONNX model was trained at FFT Size of {}, but OpenScofo is running FFT Size {}.", metadataFFTSize,
                     configuration.FFTSize);
    }

    int metadataHopSize = 0;
    if (ParseIntMetadata(hopSizeMetadata, metadataHopSize) && metadataHopSize != configuration.HOPSize) {
        spdlog::warn("ONNX model was trained at Hop Size of {}, but OpenScofo is running Hop Size {}.", metadataHopSize,
                     configuration.HOPSize);
    }

    std::vector<Descriptors> metadataDescriptors;
    if (ParseDescriptorMetadata(descriptorMetadata, metadataDescriptors) && !metadataDescriptors.empty()) {
        descriptors = std::move(metadataDescriptors);
    } else if (descriptorMetadata != nullptr) {
        spdlog::warn("Could not parse ONNX descriptor metadata, using descriptors supplied by the score/API.");
    }

    if (!ParseJsonStringArray(labelMetadata, metadataLabels) && labelMetadata != nullptr) {
        spdlog::warn("Could not parse ONNX label metadata, falling back to TreeEnsembleClassifier labels.");
    }
}

// ─────────────────────────────────────
/**
 * @brief Check that every graph node uses a supported opset.
 *
 * @return True when the graph exists and all opsets are supported; false otherwise.
 *
 * @note Rejects missing graphs and opsets newer than CurrentOpset with an error log.
 */
bool ONNXModel::ValidateOpsets() const {
    if (m_Context == nullptr || m_Context->g == nullptr) {
        spdlog::error("ONNX graph not found");
        return false;
    }

    const struct onnx_graph_t *graph = m_Context->g;
    for (int i = 0; i < graph->nlen; ++i) {
        const struct onnx_node_t *node = &graph->nodes[i];
        if (node->opset > CurrentOpset) {
            spdlog::error("Unsupported opset => {} {}.", node->proto->op_type, node->opset);
            return false;
        }
    }

    return true;
}

// ─────────────────────────────────────
/**
 * @brief Collect classifier labels from metadata or graph attributes.
 *
 * @param metadataLabels Class labels parsed from model metadata.
 *
 * @return True if a TreeEnsembleClassifier and a nonempty label list are found; false otherwise.
 *
 * @note Metadata labels take precedence over classlabels_strings attributes.
 * @warning Requires an allocated graph; classifiers without string labels or supported metadata are rejected.
 */
bool ONNXModel::ReadLabelsFromModel(const std::vector<std::string> &metadataLabels) {
    bool classifierFound = false;
    struct onnx_graph_t *graph = m_Context->g;

    for (int i = 0; i < graph->nlen; ++i) {
        struct onnx_node_t *node = &graph->nodes[i];
        if (node->proto == nullptr || std::strcmp(node->proto->op_type, "TreeEnsembleClassifier") != 0) {
            continue;
        }

        classifierFound = true;
        if (!metadataLabels.empty()) {
            continue;
        }

        for (std::size_t attributeIndex = 0; attributeIndex < node->proto->n_attribute; ++attributeIndex) {
            Onnx__AttributeProto *attribute = node->proto->attribute[attributeIndex];
            if (std::strcmp(attribute->name, "classlabels_strings") != 0) {
                continue;
            }

            for (std::size_t labelIndex = 0; labelIndex < attribute->n_strings; ++labelIndex) {
                const ProtobufCBinaryData &label = attribute->strings[labelIndex];
                m_Labels.emplace_back(reinterpret_cast<const char *>(label.data), label.len);
            }
        }
    }

    if (!metadataLabels.empty()) {
        m_Labels = metadataLabels;
    }

    if (!classifierFound) {
        spdlog::error(
            "TreeEnsembleClassifier not found in model, please use PureData py.train to train models for OpenScofo");
        return false;
    }
    if (m_Labels.empty()) {
        spdlog::error("TreeEnsembleClassifier labels not found in ONNX model");
        return false;
    }

    return true;
}

// ─────────────────────────────────────
/**
 * @brief Calculate model input width and allocate the flattened descriptor buffer.
 *
 * @param configuration Audio and descriptor dimensions used to prepare and validate model inputs.
 *
 * @return True when at least one input value is selected; false otherwise.
 *
 * @note Expands vector descriptors using configuration dimensions and ignores ONNX and INVALID entries.
 */
bool ONNXModel::PrepareDescriptorBuffer(const Configuration &configuration) {
    m_MFCCCount = configuration.MFCCCount;
    m_MFCCMels = configuration.MFCCMels;
    m_ChromaSize = configuration.ChromaSize;
    m_FFTSize = static_cast<int>(configuration.FFTSize);
    m_DescriptorSize = 0;

    for (const Descriptors descriptor : m_Descriptors) {
        if (descriptor == MFCC) {
            m_DescriptorSize += static_cast<std::size_t>(m_MFCCCount);
        } else if (descriptor == CHROMA) {
            m_DescriptorSize += static_cast<std::size_t>(m_ChromaSize);
        } else if (descriptor == POWERARRAY || descriptor == MAGNITUDE) {
            m_DescriptorSize += static_cast<std::size_t>(m_FFTSize / 2 + 1);
        } else if (descriptor == LOGMEL) {
            m_DescriptorSize += static_cast<std::size_t>(m_MFCCMels);
        } else if (descriptor != ONNX && descriptor != INVALID) {
            ++m_DescriptorSize;
        }
    }

    if (m_DescriptorSize == 0) {
        spdlog::error("No valid descriptors were selected for the ONNX model");
        return false;
    }

    m_DescriptorValues.resize(m_DescriptorSize);
    spdlog::debug("ONNX model input array has {} values", m_DescriptorSize);
    return true;
}

// ─────────────────────────────────────
/**
 * @brief Select numeric model input and probability-output tensors.
 *
 * @return True if both matching tensors are found; false otherwise.
 *
 * @note Matches input width to descriptor count and output width to label count using graph tensor names.
 * @warning Requires model context, labels, and descriptor dimensions.
 */
bool ONNXModel::FindTensors() {
    Onnx__GraphProto *graph = m_Context->model == nullptr ? nullptr : m_Context->model->graph;
    if (graph == nullptr) {
        spdlog::error("ONNX graph not found");
        return false;
    }

    for (std::size_t i = 0; i < graph->n_input; ++i) {
        const char *name = graph->input[i] == nullptr ? nullptr : graph->input[i]->name;
        struct onnx_tensor_t *tensor = onnx_tensor_search(m_Context, name);
        if (IsFloatTensor(tensor) && tensor->ndata == m_DescriptorSize) {
            m_InputTensor = tensor;
            break;
        }
    }
    if (m_InputTensor == nullptr) {
        spdlog::error("ONNX input tensor with {} float values not found", m_DescriptorSize);
        return false;
    }

    for (std::size_t i = 0; i < graph->n_output; ++i) {
        const char *name = graph->output[i] == nullptr ? nullptr : graph->output[i]->name;
        struct onnx_tensor_t *tensor = onnx_tensor_search(m_Context, name);
        if (IsFloatTensor(tensor) && tensor->ndata == m_Labels.size()) {
            m_OutputTensor = tensor;
            break;
        }
    }
    if (m_OutputTensor == nullptr) {
        spdlog::error("ONNX numeric probability output with {} values not found", m_Labels.size());
        return false;
    }

    return true;
}

// ─────────────────────────────────────
/**
 * @brief Append a scalar or vector descriptor to the flattened float input.
 *
 * @param descriptor Descriptor to append to the input buffer.
 * @param description Audio frame supplying descriptor values.
 * @param output Writable float cursor, advanced past the appended descriptor values.
 *
 * @note Advances the caller cursor for every written value; ONNX and invalid descriptors add no values.
 * @warning The cursor must have enough writable space and source arrays must match the recorded model dimensions.
 */
void ONNXModel::WriteDescriptor(Descriptors descriptor, const Description &description, float *&output) const {
    if (descriptor == MFCC) {
        for (int i = 0; i < m_MFCCCount; ++i) {
            *output++ = static_cast<float>(description.MFCC[static_cast<std::size_t>(i)]);
        }
    } else if (descriptor == CHROMA) {
        for (int i = 0; i < m_ChromaSize; ++i) {
            *output++ = static_cast<float>(description.Chroma[static_cast<std::size_t>(i)]);
        }
    } else if (descriptor == POWERARRAY) {
        for (const double value : description.Power) {
            *output++ = static_cast<float>(value);
        }
    } else if (descriptor == LOGMEL) {
        for (int i = 0; i < m_MFCCMels; ++i) {
            *output++ = static_cast<float>(description.LogMelSpectrum[static_cast<std::size_t>(i)]);
        }
    } else if (descriptor == MAGNITUDE) {
        for (const double value : description.Magnitude) {
            *output++ = static_cast<float>(value);
        }
    } else if (descriptor == LOUDNESS) {
        *output++ = static_cast<float>(description.Loudness);
    } else if (descriptor == RMS) {
        *output++ = static_cast<float>(description.RMS);
    } else if (descriptor == ZCR) {
        *output++ = static_cast<float>(description.ZeroCrossingRate);
    } else if (descriptor == HFR) {
        *output++ = static_cast<float>(description.HighFreqRatio);
    } else if (descriptor == CENTROID) {
        *output++ = static_cast<float>(description.SpectralCentroid);
    } else if (descriptor == SPREADHZ) {
        *output++ = static_cast<float>(description.SpectralSpreadHz);
    } else if (descriptor == FLATNESS) {
        *output++ = static_cast<float>(description.SpectralFlatness);
    } else if (descriptor == ENTROPY) {
        *output++ = static_cast<float>(description.SpectralEntropy);
    } else if (descriptor == ROLLOFF) {
        *output++ = static_cast<float>(description.SpectralRolloff);
    } else if (descriptor == FLUX) {
        *output++ = static_cast<float>(description.SpectralFlux);
    } else if (descriptor == IRREGULARITY) {
        *output++ = static_cast<float>(description.SpectralIrregularity);
    } else if (descriptor == HARMONICITY) {
        *output++ = static_cast<float>(description.Harmonicity);
    } else if (descriptor == YIN) {
        *output++ = static_cast<float>(description.Pitch);
    } else if (descriptor == SILENCEPROB) {
        *output++ = static_cast<float>(description.SilenceProb);
    } else if (descriptor == EXTENDEDTECHNIQUE) {
        *output++ = static_cast<float>(description.ExtendedTechProb);
    } else if (descriptor == ODSONSET) {
        *output++ = description.Onset ? 1.0F : 0.0F;
    } else if (descriptor == STDDEV) {
        *output++ = static_cast<float>(description.StdDev);
    } else if (descriptor == SKEWNESS) {
        *output++ = static_cast<float>(description.SpectralSkewness);
    } else if (descriptor == SLOPE) {
        *output++ = static_cast<float>(description.SpectralSlope);
    } else if (descriptor == DB) {
        *output++ = static_cast<float>(description.dB);
    } else if (descriptor == MAXAMP) {
        *output++ = static_cast<float>(description.MaxAmp);
    } else if (descriptor == SPREADVARIANCE) {
        *output++ = static_cast<float>(description.SpectralSpreadVariance);
    } else if (descriptor == CREST) {
        *output++ = static_cast<float>(description.SpectralCrest);
    } else if (descriptor == KURTOSIS) {
        *output++ = static_cast<float>(description.SpectralKurtosis);
    } else if (descriptor == CENTROIDVEL) {
        *output++ = static_cast<float>(description.CentroidVelocity);
    } else if (descriptor == YINCONFIDENCE) {
        *output++ = static_cast<float>(description.PitchConfidence);
    } else if (descriptor != ONNX && descriptor != INVALID) {
        spdlog::error("ONNX descriptor {} has no value writer", static_cast<int>(descriptor));
        *output++ = 0.0F;
    }
}

// ─────────────────────────────────────
/**
 * @brief Advance a text cursor past whitespace.
 *
 * @param text Source text to inspect.
 * @param position Current text offset, advanced as characters are consumed.
 *
 * @note Stops at the first non-whitespace character or the end of the view.
 */
void ONNXModel::SkipWhitespace(std::string_view text, std::size_t &position) {
    while (position < text.size() && std::isspace(static_cast<unsigned char>(text[position]))) {
        ++position;
    }
}

// ─────────────────────────────────────
/**
 * @brief Parse a quoted JSON string at the current cursor.
 *
 * @param text Source text to inspect.
 * @param position Current text offset, advanced as characters are consumed.
 * @param output Decoded string output.
 *
 * @return True if a closing quote is found and all escapes are supported; false otherwise.
 *
 * @note Supports standard single-character escapes and skips leading whitespace.
 * @warning Unicode escape sequences are unsupported; failed parsing may advance the cursor and leave partial
 * output.
 */
bool ONNXModel::ParseJsonString(std::string_view text, std::size_t &position, std::string &output) {
    SkipWhitespace(text, position);
    if (position >= text.size() || text[position] != '"') {
        return false;
    }

    ++position;
    output.clear();
    while (position < text.size()) {
        const char character = text[position++];
        if (character == '"') {
            return true;
        }
        if (character != '\\') {
            output.push_back(character);
            continue;
        }
        if (position >= text.size()) {
            return false;
        }

        const char escaped = text[position++];
        if (escaped == '"' || escaped == '\\' || escaped == '/') {
            output.push_back(escaped);
        } else if (escaped == 'b') {
            output.push_back('\b');
        } else if (escaped == 'f') {
            output.push_back('\f');
        } else if (escaped == 'n') {
            output.push_back('\n');
        } else if (escaped == 'r') {
            output.push_back('\r');
        } else if (escaped == 't') {
            output.push_back('\t');
        } else {
            return false;
        }
    }

    return false;
}

// ─────────────────────────────────────
/**
 * @brief Parse metadata as an array of quoted strings.
 *
 * @param metadata Null-terminated metadata text, or null when the entry is absent.
 * @param output Output list of decoded string values.
 *
 * @return True if the entire array is accepted; false for missing or malformed metadata.
 *
 * @note Clears the output initially and requires the complete text, apart from whitespace, to form an array.
 * @warning Failed parsing can leave parsed elements when trailing text follows a closing bracket.
 */
bool ONNXModel::ParseJsonStringArray(const char *metadata, std::vector<std::string> &output) {
    output.clear();
    if (metadata == nullptr) {
        return false;
    }

    const std::string_view text(metadata);
    std::size_t position = 0;
    SkipWhitespace(text, position);
    if (position >= text.size() || text[position] != '[') {
        return false;
    }

    ++position;
    SkipWhitespace(text, position);
    if (position < text.size() && text[position] == ']') {
        ++position;
        SkipWhitespace(text, position);
        return position == text.size();
    }

    while (position < text.size()) {
        std::string value;
        if (!ParseJsonString(text, position, value)) {
            output.clear();
            return false;
        }
        output.push_back(std::move(value));

        SkipWhitespace(text, position);
        if (position >= text.size()) {
            output.clear();
            return false;
        }
        if (text[position] == ']') {
            ++position;
            SkipWhitespace(text, position);
            return position == text.size();
        }
        if (text[position] != ',') {
            output.clear();
            return false;
        }
        ++position;
    }

    output.clear();
    return false;
}

// ─────────────────────────────────────
/**
 * @brief Parse an integer metadata value after trimming whitespace.
 *
 * @param metadata Null-terminated metadata text, or null when the entry is absent.
 * @param output Integer output, unchanged on failure.
 *
 * @return True for a valid in-range integer; false for missing or malformed metadata.
 *
 * @note Uses std::from_chars() and updates output only when the complete trimmed value parses.
 */
bool ONNXModel::ParseIntMetadata(const char *metadata, int &output) {
    if (metadata == nullptr) {
        return false;
    }

    const std::string_view text(metadata);
    std::size_t begin = 0;
    std::size_t end = text.size();
    while (begin < end && std::isspace(static_cast<unsigned char>(text[begin]))) {
        ++begin;
    }
    while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) {
        --end;
    }

    int value = 0;
    const char *first = text.data() + begin;
    const char *last = text.data() + end;
    const std::from_chars_result result = std::from_chars(first, last, value);
    if (result.ec != std::errc() || result.ptr != last) {
        return false;
    }

    output = value;
    return true;
}

// ─────────────────────────────────────
/**
 * @brief Map a metadata descriptor name to its enum.
 *
 * @param name Descriptor identifier stored in model metadata.
 *
 * @return Matching descriptor enum, or INVALID for an unknown name.
 *
 * @note Accepts supported canonical names and aliases without case conversion.
 */
Descriptors ONNXModel::DescriptorFromMetadataName(std::string_view name) {
    if (name == "mfcc")
        return MFCC;
    if (name == "logmel")
        return LOGMEL;
    if (name == "rms")
        return RMS;
    if (name == "loudness")
        return LOUDNESS;
    if (name == "db")
        return DB;
    if (name == "maxamp" || name == "max_amp")
        return MAXAMP;
    if (name == "magnitude")
        return MAGNITUDE;
    if (name == "power" || name == "powerarray")
        return POWERARRAY;
    if (name == "stddev")
        return STDDEV;
    if (name == "chroma")
        return CHROMA;
    if (name == "silence")
        return SILENCEPROB;
    if (name == "harmonicity")
        return HARMONICITY;
    if (name == "centroid")
        return CENTROID;
    if (name == "zcr")
        return ZCR;
    if (name == "hfr")
        return HFR;
    if (name == "spread")
        return SPREADHZ;
    if (name == "spread_variance")
        return SPREADVARIANCE;
    if (name == "crest")
        return CREST;
    if (name == "flatness")
        return FLATNESS;
    if (name == "entropy")
        return ENTROPY;
    if (name == "rolloff")
        return ROLLOFF;
    if (name == "centroid_velocity")
        return CENTROIDVEL;
    if (name == "flux")
        return FLUX;
    if (name == "skewness")
        return SKEWNESS;
    if (name == "slope")
        return SLOPE;
    if (name == "kurtosis")
        return KURTOSIS;
    if (name == "ext")
        return EXTENDEDTECHNIQUE;
    if (name == "onset")
        return ODSONSET;
    if (name == "irregularity")
        return IRREGULARITY;
    if (name == "yin")
        return YIN;
    if (name == "yin_confidence" || name == "pitch_confidence")
        return YINCONFIDENCE;
    if (name == "onnx")
        return ONNX;
    return INVALID;
}

// ─────────────────────────────────────
/**
 * @brief Parse an ordered descriptor-name array from metadata.
 *
 * @param metadata Null-terminated metadata text, or null when the entry is absent.
 * @param output Output descriptors in metadata order; unchanged if array parsing fails, cleared for unknown names.
 *
 * @return True when the string array parses and every name is recognized; false otherwise.
 *
 * @note Converts every name and rejects the entire selection when a name is unsupported.
 */
bool ONNXModel::ParseDescriptorMetadata(const char *metadata, std::vector<Descriptors> &output) {
    std::vector<std::string> names;
    if (!ParseJsonStringArray(metadata, names)) {
        return false;
    }

    output.clear();
    for (const std::string &name : names) {
        const Descriptors descriptor = DescriptorFromMetadataName(name);
        if (descriptor == INVALID) {
            spdlog::error("Invalid ONNX metadata descriptor '{}'", name);
            output.clear();
            return false;
        }
        output.push_back(descriptor);
    }
    return true;
}

// ─────────────────────────────────────
/**
 * @brief Check whether a tensor contains float32 or float64 values.
 *
 * @param tensor Tensor pointer to inspect; may be null.
 *
 * @return True for a floating-point tensor; false otherwise.
 *
 * @note Null pointers and non-floating tensor types are rejected.
 */
bool ONNXModel::IsFloatTensor(const struct onnx_tensor_t *tensor) {
    return tensor != nullptr && (tensor->type == ONNX_TENSOR_TYPE_FLOAT32 || tensor->type == ONNX_TENSOR_TYPE_FLOAT64);
}

} // namespace OpenScofo
