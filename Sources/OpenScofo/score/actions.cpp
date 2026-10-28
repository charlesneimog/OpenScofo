/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file actions.cpp
 * @brief Score action and Lua-related command extraction.
 */

#include "OpenScofo.hpp"
#include <algorithm>
#include <cctype>
#include <tree_sitter/api.h>

namespace OpenScofo {

// ─────────────────────────────────────
/**
 * @brief Extract timed Lua or receiver commands and append event actions.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 * @param Event Score event to populate or extend in place.
 *
 * @note Stores absolute delays in milliseconds and tempo-relative delays in beats, preserving command arguments.
 * @warning Numeric conversions can throw; a malformed command can leave earlier actions already appended.
 */
void Score::NewEventAction(const std::string &ScoreStr, TSNode Node, ScoreState &Event) {
    ScoreAction BaseAction;
    BaseAction.AbsoluteTime = true;
    BaseAction.Time = 0;

    TSNode timingNode = GetField(Node, "timing");
    if (!ts_node_is_null(timingNode)) {
        TSNode amountNode = GetField(timingNode, "amount");
        TSNode unitNode = GetField(timingNode, "unit");

        if (!ts_node_is_null(amountNode)) {
            BaseAction.Time = std::stof(GetCodeStr(ScoreStr, amountNode));
        }

        if (!ts_node_is_null(unitNode)) {
            std::string unit = GetCodeStr(ScoreStr, unitNode);
            if (unit == "sec") {
                BaseAction.Time *= 1000.0;
                BaseAction.AbsoluteTime = true;
            } else if (unit == "ms") {
                BaseAction.AbsoluteTime = true;
            } else if (unit == "tempo") {
                BaseAction.AbsoluteTime = false;
            }
        }
    }

    uint32_t childCount = ts_node_child_count(Node);
    for (uint32_t i = 0; i < childCount; ++i) {
        const char *fieldName = ts_node_field_name_for_child(Node, i);
        if (fieldName == nullptr || std::string(fieldName) != "command") {
            continue;
        }

        TSNode execNode = ts_node_child(Node, i);
        if (std::string(ts_node_type(execNode)) != "exec") {
            continue;
        }

        ScoreAction NewAction = BaseAction;

        TSNode luaNode = GetField(execNode, "lua");
        TSNode receiverNode = GetField(execNode, "receiver");

        if (!ts_node_is_null(luaNode)) {
            NewAction.Lua = GetCodeStr(ScoreStr, luaNode);
            NewAction.isLua = true;
        } else if (!ts_node_is_null(receiverNode)) {
            NewAction.Receiver = GetCodeStr(ScoreStr, receiverNode);
            NewAction.isLua = false;

            TSNode args = GetField(execNode, "args");
            if (!ts_node_is_null(args)) {
                uint32_t argsCount = ts_node_child_count(args);
                std::string pendingErrorPrefix;
                for (uint32_t j = 0; j < argsCount; j++) {
                    TSNode arg = ts_node_child(args, j);
                    std::string argType = ts_node_type(arg);

                    if (argType == "ERROR") {
                        std::string token = GetCodeStr(ScoreStr, arg);
                        token.erase(
                            std::remove_if(token.begin(), token.end(), [](unsigned char c) { return std::isspace(c); }),
                            token.end());
                        pendingErrorPrefix += token;
                        continue;
                    }

                    if (argType != "pdarg") {
                        continue;
                    }

                    TSNode pdarg = ts_node_child(arg, 0);
                    if (ts_node_is_null(pdarg)) {
                        continue;
                    }

                    std::string pdargType = ts_node_type(pdarg);
                    std::string token = GetCodeStr(ScoreStr, pdarg);
                    if (!pendingErrorPrefix.empty()) {
                        token = pendingErrorPrefix + token;
                        pendingErrorPrefix.clear();
                    }

                    if (pdargType == "number") {
                        if (isNumber(token)) {
                            float f = std::stof(token);
                            NewAction.Args.emplace_back(f);
                        } else {
                            spdlog::error("Invalid number argument on line {}.", ts_node_start_point(pdarg).row + 1);
                            return;
                        }
                    } else if (pdargType == "identifier" || pdargType == "symbol") {
                        std::string s = token;
                        NewAction.Args.emplace_back(s);
                    } else if (!token.empty()) {
                        std::string s = token;
                        NewAction.Args.emplace_back(s);
                    }
                }

                if (!pendingErrorPrefix.empty()) {
                    NewAction.Args.push_back(pendingErrorPrefix);
                }
            }
        } else {
            TSPoint Pos = ts_node_start_point(execNode);
            spdlog::error("Invalid action command on line {}.", Pos.row + 1);
            continue;
        }

        Event.Actions.push_back(NewAction);
    }
}

} // namespace OpenScofo
