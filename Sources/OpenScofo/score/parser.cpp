/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file parser.cpp
 * @brief Tree-sitter utilities, syntax validation, and source-text extraction.
 */

#include "OpenScofo.hpp"
#include <algorithm>
#include <charconv>
#include <cctype>
#include <fstream>
#include <iostream>
#include <string_view>
#include <system_error>
#include <tree_sitter/api.h>

namespace OpenScofo {

// ─────────────────────────────────────
/**
 * @brief Print a syntax subtree recursively for debugging.
 *
 * @param node Tree-sitter node to inspect in the current score tree.
 * @param indent Initial indentation width in spaces.
 *
 * @note Indents child nodes by four additional spaces and omits printing the root when indent is zero.
 */
void Score::PrintTreeSitterNode(TSNode node, int indent) {
    const char *type = ts_node_type(node);
    std::string text = ts_node_string(node);
    if (indent != 0) {
        std::cout << std::string(indent, ' ') << type << ": " << text << std::endl;
    }
    uint32_t child_count = ts_node_child_count(node);
    for (uint32_t i = 0; i < child_count; i++) {
        PrintTreeSitterNode(ts_node_child(node, i), indent + 4);
    }
}

// ─────────────────────────────────────
/**
 * @brief Find a child node by its Tree-sitter field name.
 *
 * @param Node Tree-sitter node to inspect in the current score tree.
 * @param s Tree-sitter field name to look up.
 *
 * @return Matching field node, or a null node if the field is absent.
 *
 * @note Delegates to ts_node_child_by_field_name().
 */
TSNode Score::GetField(TSNode Node, std::string s) {
    int strLen = s.length();
    TSNode field = ts_node_child_by_field_name(Node, s.c_str(), strLen);
    return field;
}

// ─────────────────────────────────────
/**
 * @brief Check whether text is a complete floating-point literal.
 *
 * @param str Text to validate as a complete floating-point number.
 *
 * @return True if the complete text parses as float; false otherwise.
 *
 * @note Uses std::from_chars() and rejects empty text or trailing unparsed characters.
 */
bool Score::isNumber(const std::string &str) {
    if (str.empty()) {
        return false;
    }

    float value = 0.0f;
    const char *begin = str.data();
    const char *end = begin + str.size();
    auto [ptr, ec] = std::from_chars(begin, end, value);
    return ec == std::errc{} && ptr == end;
}

// ╭─────────────────────────────────────╮
// │       Parse File of the Score       │
// ╰─────────────────────────────────────╯
/**
 * @brief Copy the source text covered by a syntax node.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Source substring covered by Node.
 *
 * @note Uses Tree-sitter start and end byte offsets to select the text.
 * @warning Node must be valid and its byte range must lie within ScoreStr.
 */
std::string Score::GetCodeStr(const std::string &ScoreStr, TSNode Node) {
    int start = ts_node_start_byte(Node);
    int end = ts_node_end_byte(Node);
    return std::string(std::string_view(ScoreStr.data() + start, end - start));
}

// ─────────────────────────────────────
/**
 * @brief Read child source text by field name or fallback node type.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param node Tree-sitter node to inspect in the current score tree.
 * @param id Field name, fallback node type, or configuration key.
 *
 * @return Matching child text, or an empty string if no child matches.
 *
 * @note Tries a named field first, then inspects direct children for a matching type.
 */
std::string Score::GetChildStringFromField(const std::string &ScoreStr, TSNode node, std::string id) {
    TSNode field = ts_node_child_by_field_name(node, id.c_str(), id.length());
    if (!ts_node_is_null(field)) {
        return GetCodeStr(ScoreStr, field);
    }

    int child_count = ts_node_child_count(node);
    for (int i = 0; i < child_count; i++) {
        TSNode child = ts_node_child(node, i);
        const char *type = ts_node_type(child);
        if (id == type) {
            return GetCodeStr(ScoreStr, child);
        }
    }
    return "";
}

// ─────────────────────────────────────
/**
 * @brief Read the first direct named child of a specified type.
 *
 * @param source Complete source text corresponding to the parent node byte offsets.
 * @param parent Parent node whose named children are inspected.
 * @param wanted_type Named child node type to find.
 *
 * @return Matching child source text, or an empty string if no child matches.
 *
 * @note Ignores unnamed children and extracts text using the matched child byte range.
 */
std::string Score::GetChildStringFromType(const std::string &source, TSNode parent, const std::string &wanted_type) {
    uint32_t count = ts_node_child_count(parent);
    for (uint32_t i = 0; i < count; ++i) {
        TSNode child = ts_node_child(parent, i);

        if (!ts_node_is_named(child))
            continue;

        if (wanted_type == ts_node_type(child)) {
            uint32_t start = ts_node_start_byte(child);
            uint32_t end = ts_node_end_byte(child);
            return source.substr(start, end - start);
        }
    }
    return {};
}

// ─────────────────────────────────────
/**
 * @brief Log missing tokens and unexpected source text in a syntax tree.
 *
 * @param Node Tree-sitter node to inspect in the current score tree.
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param InsideError Whether an ancestor is already an error node, suppressing duplicate recovery messages.
 *
 * @note Traverses recursively, avoiding duplicate nested error messages and shortening recovery text to 80
 * characters.
 */
void Score::FindErrors(TSNode Node, const std::string &ScoreStr, bool InsideError) {
    if (ts_node_is_null(Node)) {
        return;
    }

    const TSPoint Position = ts_node_start_point(Node);
    if (ts_node_is_missing(Node)) {
        std::string Label = ts_node_type(Node);
        const TSNode Parent = ts_node_parent(Node);
        if (Label == "number" && !ts_node_is_null(Parent)) {
            for (const std::string Field : {"duration", "amount"}) {
                const TSNode Value = ts_node_child_by_field_name(Parent, Field.c_str(), Field.size());
                if (!ts_node_is_null(Value) && ts_node_eq(Value, Node)) {
                    Label = Field;
                    break;
                }
            }
        }
        if (!ts_node_is_named(Node)) {
            Label = "\"" + Label + "\"";
        }
        spdlog::error("Missing {} at line {}, column {}", Label, Position.row + 1, Position.column + 1);
    } else if (ts_node_is_error(Node) && !InsideError) {
        const uint32_t Start = ts_node_start_byte(Node);
        std::string Text = ScoreStr.substr(Start, ts_node_end_byte(Node) - Start);
        // Keep multi-line recovery nodes readable in host consoles.
        std::replace_if(Text.begin(), Text.end(), [](unsigned char C) { return std::isspace(C); }, ' ');
        if (Text.size() > 80) {
            Text = Text.substr(0, 80) + "…";
        }
        spdlog::error("Unexpected text '{}' at line {}, column {}", Text, Position.row + 1, Position.column + 1);
    }

    // has_error also marks ancestors; only the actual ERROR/MISSING nodes get messages.
    for (uint32_t Index = 0; Index < ts_node_child_count(Node); ++Index) {
        FindErrors(ts_node_child(Node, Index), ScoreStr, InsideError || ts_node_is_error(Node));
    }
}

// ─────────────────────────────────────
/**
 * @brief Heuristically classify a file as text from its leading bytes.
 *
 * @param path Path to the file whose leading bytes are inspected.
 *
 * @return True if the sample appears textual; false for unreadable or apparently binary files.
 *
 * @note Samples up to 4096 bytes, accepts empty files, and rejects null bytes or a suspicious-byte ratio of at
 * least 5%.
 * @warning This byte heuristic is not a Unicode-aware validator and can reject non-ASCII text.
 */
bool Score::ScoreIsText(const std::string &path) {
    std::ifstream file(path, std::ios::binary);
    if (!file)
        return false;

    constexpr std::size_t SampleSize = 4096;
    std::vector<unsigned char> buffer(SampleSize);

    file.read(reinterpret_cast<char *>(buffer.data()), SampleSize);
    std::size_t bytesRead = static_cast<std::size_t>(file.gcount());

    if (bytesRead == 0)
        return true; // empty file → treat as text

    std::size_t suspicious = 0;

    for (std::size_t i = 0; i < bytesRead; ++i) {
        unsigned char c = buffer[i];

        // Null byte strongly indicates binary
        if (c == 0)
            return false;

        // Allow printable ASCII and common whitespace
        if (!(std::isprint(c) || c == '\n' || c == '\r' || c == '\t'))
            suspicious++;
    }

    double ratio = static_cast<double>(suspicious) / bytesRead;

    // Threshold: >5% suspicious bytes → likely binary
    return ratio < 0.05;
}

} // namespace OpenScofo
