/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file config.cpp
 * @brief Score configuration parsing, validation, and assignment.
 */

#include "OpenScofo.hpp"
#include <algorithm>
#include <cctype>
#include <tree_sitter/api.h>

namespace OpenScofo {

// ─────────────────────────────────────
/**
 * @brief Parse a numeric configuration value.
 *
 * @param id Field name, fallback node type, or configuration key.
 * @param valueType Tree-sitter type of the configuration value.
 * @param value Source text of the configuration value.
 * @param pos Source position used for diagnostics.
 * @param out Numeric output assigned on success.
 *
 * @return True when the value is accepted and assigned; false for the wrong syntax type.
 *
 * @note Logs and rejects values whose syntax type is not number.
 * @warning std::stod() can throw for malformed or out-of-range text.
 */
bool Score::GetConfigNumber(const std::string &id, const std::string &valueType, const std::string &value, TSPoint pos,
                            double &out) {
    if (valueType != "number") {
        spdlog::error("Invalid numeric value for {} on line {}.", id, pos.row + 1);
        return false;
    }

    out = std::stod(value);
    return true;
}

// ─────────────────────────────────────
/**
 * @brief Parse a boolean configuration value and supported aliases.
 *
 * @param id Field name, fallback node type, or configuration key.
 * @param valueType Tree-sitter type of the configuration value.
 * @param value Source text of the configuration value.
 * @param pos Source position used for diagnostics.
 * @param out Boolean output, unchanged on failure.
 *
 * @return True if a supported value was assigned; false otherwise.
 *
 * @note Accepts case-insensitive true/false, on/off, yes/no, and numeric 1/0.
 */
bool Score::GetConfigBool(const std::string &id, const std::string &valueType, std::string value, TSPoint pos,
                          bool &out) {
    if (valueType != "identifier" && valueType != "number") {
        spdlog::error("Invalid boolean value for {} on line {}.", id, pos.row + 1);
        return false;
    }

    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (value == "true" || value == "on" || value == "yes" || value == "1") {
        out = true;
        return true;
    }
    if (value == "false" || value == "off" || value == "no" || value == "0") {
        out = false;
        return true;
    }

    spdlog::error("Invalid boolean value for {} on line {}.", id, pos.row + 1);
    return false;
}

// ─────────────────────────────────────
/**
 * @brief Apply a parsed score configuration declaration.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param node Tree-sitter node to inspect in the current score tree.
 * @param Config Score configuration to read or update.
 *
 * @note Validates supported keys, updates parser timing state, and resolves model paths relative to the score file.
 * @warning Some settings mutate parser state or append boundaries; numeric conversions can throw.
 */
void Score::NewConfig(const std::string &ScoreStr, TSNode node, Configuration &Config) {
    TSNode keyNode = GetField(node, "key");
    TSNode valueNode = GetField(node, "value");
    TSPoint pos = ts_node_start_point(node);

    if (ts_node_is_null(keyNode) || ts_node_is_null(valueNode)) {
        spdlog::error("Invalid CONFIG on line {}.", pos.row + 1);
        return;
    }

    std::string id = GetCodeStr(ScoreStr, keyNode);
    std::string valueType = ts_node_type(valueNode);
    std::string value = GetCodeStr(ScoreStr, valueNode);

    if (id == "BPM" || id == "TRANSPOSE" || id == "SR" || id == "FFTSIZE" || id == "HOPSIZE" || id == "TUNINGA4" ||
        id == "PHASECOUPLING" || id == "SYNCSTRENGTH" || id == "PITCHTEMPLATESIGMA" || id == "PITCHTEMPLATEHARMONICS" ||
        id == "MFCCMELS" || id == "MFCCCOUNT" || id == "MEDSPAN" || id == "DBTRESHOLD" || id == "DBTHRESHOLD" ||
        id == "SPECTRALROLLOFFCUTOFF" || id == "YINTHRESHOLD" || id == "YINMINFREQUENCY" || id == "YINMAXFREQUENCY" ||
        id == "CHROMASIZE" || id == "CHROMACENTEROCTAVE" || id == "CHROMAOCTAVEWIDTH" || id == "ZCRTHRESHOLD" ||
        id == "TIMETOLERANCE") {
        double v = 0.0;
        if (!GetConfigNumber(id, valueType, value, pos, v)) {
            return;
        }

        if (id == "BPM") {
            m_CurrentBPM = v;
            if (m_CurrentBPM < 1) {
                spdlog::error("BPM must be bigger than 0");
                m_CurrentBPM = 1;
            }

            ScoreState Begin = GetFirstEvent();
            ProcessEventTime(Begin);
            Begin.BPMExpected = v;
            m_ScoreStates.emplace_back(Begin);
        } else if (id == "TRANSPOSE") {
            if (v < -36 || v > 36) {
                spdlog::warn("Weird transpose value on line {}.", pos.row + 1);
            }
            m_Transpose = v;
        } else if (id == "PHASECOUPLING") {
            if (v < 0 || v > 2) {
                spdlog::error("Invalid value for PHASECOUPLING on line {}.", pos.row + 1);
            } else {
                Config.PhaseCoupling = v;
            }
        } else if (id == "SYNCSTRENGTH") {
            if (v < 0 || v > 1) {
                spdlog::error("Invalid value for SYNCSTRENGTH on line {}.", pos.row + 1);
            } else {
                Config.SyncStrength = v;
            }
        } else if (id == "TIMETOLERANCE") {
            double tolerance = std::clamp(v, 0.0, 1.0);
            double r = 64.0 - tolerance * 60.0;
            m_TimeTolerance = r;
        } else if (id == "PITCHTEMPLATESIGMA") {
            if (v < 0 || v > 1) {
                spdlog::error("Invalid value for PITCHTEMPLATESIGMA on line {}.", pos.row + 1);
            } else {
                Config.PitchTemplateSigma = v;
            }
        } else if (id == "PITCHTEMPLATEHARMONICS") {
            int harmonics = static_cast<int>(v);
            if (harmonics > 0) {
                Config.PitchTemplateHarmonics = harmonics;
            } else {
                spdlog::error("PITCHTEMPLATEHARMONICS must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "SR") {
            int sr = static_cast<int>(v);
            if (sr > 0) {
                Config.SR = sr;
            } else {
                spdlog::error("SR must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "FFTSIZE") {
            int fft = static_cast<int>(v);
            if (fft > 0 && (fft & (fft - 1)) == 0) {
                Config.FFTSize = fft;
            } else {
                spdlog::error("FFTSIZE must be a power of two.");
            }
        } else if (id == "HOPSIZE") {
            int hop = static_cast<int>(v);
            if (hop > 0 && (hop & (hop - 1)) == 0) {
                Config.HOPSize = hop;
            } else {
                spdlog::error("HOPSIZE must be a power of two.");
            }
        } else if (id == "TUNINGA4") {
            if (v > 0) {
                Config.TuningA4 = v;
                m_Tunning = v;
            } else {
                spdlog::error("{} must be bigger than 0 on line {}.", id, pos.row + 1);
            }
        } else if (id == "MFCCMELS") {
            int mels = static_cast<int>(v);
            if (mels > 0) {
                Config.MFCCMels = mels;
            } else {
                spdlog::error("MFCCMELS must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "MFCCCOUNT") {
            int count = static_cast<int>(v);
            if (count > 0) {
                Config.MFCCCount = count;
            } else {
                spdlog::error("MFCCCOUNT must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "MEDSPAN") {
            int medSpan = static_cast<int>(v);
            if (medSpan > 0) {
                Config.MedSpan = medSpan;
            } else {
                spdlog::error("MEDSPAN must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "DBTRESHOLD" || id == "DBTHRESHOLD") {
            Config.dBTreshold = v;
        } else if (id == "SPECTRALROLLOFFCUTOFF") {
            if (v >= 0 && v <= 1) {
                Config.SpectralRolloffCutoff = v;
            } else {
                spdlog::error("SPECTRALROLLOFFCUTOFF must be between 0 and 1 on line {}.", pos.row + 1);
            }
        } else if (id == "YINTHRESHOLD") {
            if (v >= 0 && v <= 1) {
                Config.YINThreshold = v;
            } else {
                spdlog::error("YINTHRESHOLD must be between 0 and 1 on line {}.", pos.row + 1);
            }
        } else if (id == "YINMINFREQUENCY") {
            if (v > 0) {
                Config.YINMinFrequency = v;
            } else {
                spdlog::error("YINMINFREQUENCY must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "YINMAXFREQUENCY") {
            if (v > 0) {
                Config.YINMaxFrequency = v;
            } else {
                spdlog::error("YINMAXFREQUENCY must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "CHROMASIZE") {
            int size = static_cast<int>(v);
            if (size > 0) {
                Config.ChromaSize = size;
            } else {
                spdlog::error("CHROMASIZE must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "CHROMACENTEROCTAVE") {
            Config.ChromaCenterOctave = v;
        } else if (id == "CHROMAOCTAVEWIDTH") {
            if (v > 0) {
                Config.ChromaOctaveWidth = v;
            } else {
                spdlog::error("CHROMAOCTAVEWIDTH must be bigger than 0 on line {}.", pos.row + 1);
            }
        } else if (id == "ZCRTHRESHOLD") {
            if (v >= 0) {
                Config.ZCRThreshold = v;
            } else {
                spdlog::error("ZCRTHRESHOLD must be non-negative on line {}.", pos.row + 1);
            }
        } else {
            spdlog::error("Unrecognized variable " + id);
            return;
        }
        return;
    }

    if (id == "SECTIONRESTRICT") {
        if (valueType != "identifier") {
            spdlog::error("{} must be ON or OFF on line {}.", id, pos.row + 1);
            return;
        }

        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char Character) { return static_cast<char>(std::toupper(Character)); });
        if (value != "ON" && value != "OFF") {
            spdlog::error("{} must be ON or OFF on line {}.", id, pos.row + 1);
            return;
        }
        Config.SectionRestrict = value == "ON";
        return;
    }

    if (id == "ZCRCENTER" || id == "ZCRPAD" || id == "ZCRZEROPOS") {
        bool v = false;
        if (!GetConfigBool(id, valueType, value, pos, v)) {
            return;
        }

        if (id == "ZCRCENTER") {
            Config.ZCRCenter = v;
        } else if (id == "ZCRPAD") {
            Config.ZCRPad = v;
        } else if (id == "ZCRZEROPOS") {
            Config.ZCRZeroPos = v;
        }
        return;
    }

    if (id == "ONSETFUNCTION") {
        if (valueType != "onset_function" && valueType != "identifier") {
            spdlog::error("Invalid onset function value for {} on line {}.", id, pos.row + 1);
            return;
        }

        if (value == "pow") {
            Config.OnsetType = ODS_ODF_POWER;
        } else if (value == "sf" || value == "hfc") {
            Config.OnsetType = ODS_ODF_MAGSUM;
        } else if (value == "cd") {
            Config.OnsetType = ODS_ODF_COMPLEX;
        } else if (value == "rcd") {
            Config.OnsetType = ODS_ODF_RCOMPLEX;
        } else if (value == "pd") {
            Config.OnsetType = ODS_ODF_PHASE;
        } else if (value == "wpd") {
            Config.OnsetType = ODS_ODF_WPHASE;
        } else if (value == "mkl") {
            Config.OnsetType = ODS_ODF_MKL;
        } else {
            spdlog::error("Invalid onset function {} on line {}.", value, pos.row + 1);
        }
        return;
    }

    if (id == "ONNXMODEL" || id == "TIMBREMODEL") {
        if (valueType != "path") {
            spdlog::error("Invalid path value for {} on line {}.", id, pos.row + 1);
            return;
        }

        std::string path = value;
        if (!path.empty() && path.front() == '"' && path.size() >= 2 && path.back() == '"') {
            path = path.substr(1, path.size() - 2);
        }

        Config.TimbreONNXModel = m_ScoreRootPath / fs::path(path);
        if (!fs::exists(Config.TimbreONNXModel)) {
            spdlog::error("Model path not found: {}", Config.TimbreONNXModel.string());
        }
        return;
    }

    if (id == "ONAUDIOSTATECHANGE") {
        if (valueType != "identifier" && valueType != "number") {
            spdlog::error("Invalid receiver for {} on line {}.", id, pos.row + 1);
            return;
        }
        Config.AudioStateChangeReceiver = value;
        return;
    }

    if (id == "ONNXDESCRIPTORS") {
        Config.ONNXDescriptors.clear();
        uint32_t count = ts_node_named_child_count(valueNode);
        for (uint32_t i = 0; i < count; ++i) {
            TSNode child = ts_node_named_child(valueNode, i);
            std::string d = GetCodeStr(ScoreStr, child);
            Config.ONNXDescriptors.push_back(d);
        }
        return;
    }
}

} // namespace OpenScofo
