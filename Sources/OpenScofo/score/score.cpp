/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file score.cpp
 * @brief Score parsing coordination, state management, and cleanup.
 *
 * @note Builds internal score states, section boundaries, timing metadata, and global Lua source.
 * @warning Parsing mutates the score object and logs recoverable syntax errors; serialize access to parser state.
 */

#include "OpenScofo.hpp"
#include <fstream>
#include <sstream>
#include <utility>
#include <tree_sitter/api.h>

namespace OpenScofo {

extern "C" TSLanguage *tree_sitter_openscofo();

// ─────────────────────────────────────
/**
 * @brief Copy the Lua source collected from global score blocks.
 *
 * @return Collected global Lua source.
 *
 * @note The source is accumulated during Parse() and is not executed by the parser.
 */
std::string Score::GetLuaCode() {
    return m_LuaCode;
}

// ─────────────────────────────────────
/**
 * @brief Read the parser loaded flag.
 *
 * @return True if the stored loaded flag is set; false otherwise.
 *
 * @note Set when Parse() reaches its final return; does not independently validate runtime readiness.
 * @warning The flag does not guarantee an error-free score and is not reset at the start of a subsequent parse.
 */
bool Score::ScoreIsLoaded() {
    return m_ScoreLoaded;
}

// ─────────────────────────────────────
/**
 * @brief Read a score file and build its configuration and internal state vector.
 *
 * @param ScoreFilePath Path to the score file to parse.
 *
 * @return Parsed configuration and states; a default pair on missing files, read failure, or events before BPM.
 *
 * @note Collects global Lua blocks and processes configurations, sections, and events in source order.
 * @warning Clears previous states and Lua source before loading; syntax errors can be logged while returning
 * partial results.
 */
std::pair<Configuration, States> Score::Parse(fs::path ScoreFilePath) {
    m_ScoreStates.clear();
    m_LuaCode.clear();
    Configuration Config = Configuration();

    if (fs::exists(ScoreFilePath) == false) {
        spdlog::error("Score File not found");
        return {};
    }
    m_ScoreRootPath = ScoreFilePath.parent_path();

    // Open the score file for reading
    std::ifstream File(ScoreFilePath, std::ios::binary);
    if (File.is_open() == false) {
        spdlog::error("Not possible to open score file");
        return {};
    }

    File.clear(); // Clear error flags

    std::ostringstream Buffer;
    Buffer << File.rdbuf(); // Safely read the entire file
    std::string ScoreStr = Buffer.str();

    // Proceed with parsing ScoreStr...
    m_LineCount = 0;
    m_MarkovIndex = 0;
    m_LastOnset = 0;
    m_PrevDuration = 0;
    m_CurrentBPM = -1;
    m_Transpose = 0;
    m_CurrentSection.clear();
    m_HasSection = false;
    m_SectionStartPending = false;
    m_ScorePosition = 0;
    std::string Line;

    // read and process score
    TSParser *parser = ts_parser_new();
    ts_parser_set_language(parser, tree_sitter_openscofo());
    TSTree *tree = ts_parser_parse_string(parser, nullptr, ScoreStr.c_str(), ScoreStr.size());
    TSNode rootNode = ts_tree_root_node(tree);

    if (ts_node_has_error(rootNode)) {
        FindErrors(rootNode, ScoreStr);
    }

    uint32_t child_count = ts_node_child_count(rootNode);
    for (uint32_t i = 0; i < child_count; i++) {
        TSNode Child = ts_node_child(rootNode, i);
        std::string type = ts_node_type(Child);
        if (type == "EVENT") {
            if (m_CurrentBPM == -1) {
                spdlog::error("BPM is not defined");
                return {};
            }
            NewEvent(ScoreStr, Child, Config);
        } else if (type == "CONFIG") {
            NewConfig(ScoreStr, Child, Config);
        } else if (type == "SECTION") {
            NewSection(ScoreStr, Child);
        } else if (type == "LUA") {
            std::string lua_body = GetChildStringFromField(ScoreStr, Child, "lua_body");
            m_LuaCode += lua_body;
        } else if (type == "comment") {
            // ignore
        } else {
            spdlog::error("Not recognized {}", type);
        }
    }

    // Cleanup
    ts_tree_delete(tree);
    ts_parser_delete(parser);

    m_ScoreLoaded = true;
    return {Config, m_ScoreStates};
}

} // namespace OpenScofo
