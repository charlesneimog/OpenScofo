#include <OpenScofo.hpp>
#include <algorithm>
#include <filesystem>
#include <unordered_set>
#include <vector>

namespace {

std::string NodeText(TSNode node, const std::string &source) {
    uint32_t start = ts_node_start_byte(node);
    uint32_t end = ts_node_end_byte(node);

    if (start >= source.size() || end > source.size() || start >= end)
        return {};

    return source.substr(start, end - start);
}

bool IsPunctuationOnly(const std::string &text) {
    if (text.empty())
        return false;

    static const std::string punctuation = "()[]{}.,:;";
    return std::all_of(text.begin(), text.end(), [&](char c) { return punctuation.find(c) != std::string::npos; });
}

bool IsHiddenLiteral(const std::string &text) {
    static const std::unordered_set<std::string> hidden_literals = {
        "NOTE", "REST", "CHORD", "TRILL", "TECH", "ACTION", "delay", "sendto", "luacall", "@", "LUA",
    };

    return hidden_literals.find(text) != hidden_literals.end();
}

bool ShouldPrintNode(TSNode node, const std::string &source) {
    if (ts_node_is_named(node))
        return true;

    std::string text = NodeText(node, source);
    return !IsPunctuationOnly(text) && !IsHiddenLiteral(text);
}

struct VisibleChild {
    TSNode node;
    uint32_t index;
    const char *field;
};

std::vector<VisibleChild> GetVisibleChildren(TSNode node, const std::string &source) {
    std::vector<VisibleChild> children;
    uint32_t count = ts_node_child_count(node);
    children.reserve(count);

    for (uint32_t i = 0; i < count; ++i) {
        TSNode child = ts_node_child(node, i);
        if (!ShouldPrintNode(child, source))
            continue;

        children.push_back({
            child,
            i,
            ts_node_field_name_for_child(node, i),
        });
    }

    return children;
}

std::string NodeLabel(TSNode node, const std::string &source) {
    const std::string type = ts_node_type(node);
    if (!ts_node_is_named(node))
        return type + " \"" + NodeText(node, source) + "\"";

    auto children = GetVisibleChildren(node, source);

    if (children.size() == 1) {
        const auto &only_child = children.front();
        if (!ts_node_is_named(only_child.node) && only_child.field == nullptr)
            return NodeText(only_child.node, source);
    }

    if (children.empty()) {
        const std::string text = NodeText(node, source);
        if (!text.empty() && text != type)
            return text;
    }

    return type;
}

} // namespace

extern "C" TSLanguage *tree_sitter_openscofo();

// ─────────────────────────────────────
void PrintTreeSitterNode(TSNode node, const std::string &source, const std::string &prefix = "", bool is_last = true,
                         const char *field_name = nullptr, bool is_root = true) {
    if (!ShouldPrintNode(node, source))
        return;

    std::cout << prefix;
    if (!is_root)
        std::cout << (is_last ? "└─ " : "├─ ");

    const std::string node_label = NodeLabel(node, source);
    const std::string node_type = ts_node_type(node);

    if (field_name) {
        const std::string field = field_name;
        const bool hide_redundant =
            field == node_label || field == node_type || (field == "definition" && node_type != "definition") ||
            (field == "action" && node_type == "action") || (field == "timing" && node_type == "delay") ||
            (field == "command" && node_type == "exec");

        if (hide_redundant)
            std::cout << (field == "timing" ? field : node_label);
        else
            std::cout << field << ": " << node_label;
    } else {
        std::cout << node_label;
    }
    std::cout << std::endl;

    const auto children = GetVisibleChildren(node, source);
    const std::string child_prefix = prefix + (is_root ? "" : (is_last ? "   " : "│  "));

    for (size_t i = 0; i < children.size(); ++i) {
        const bool child_is_last = (i + 1 == children.size());
        const auto &child = children[i];
        PrintTreeSitterNode(child.node, source, child_prefix, child_is_last, child.field, false);
    }
}

// ─────────────────────────────────────
int main() {
    std::string m_ScoreRootPath = //"/home/neimog/Documents/Git/OpenScofo/Tests/assets/all-tokens.txt";
        "/home/neimog/Documents/Git/OpenScofo/Tests/assets/canticos.txt";
    std::ifstream File(m_ScoreRootPath, std::ios::binary);

    if (File.is_open() == false) {
        spdlog::error("Not possible to open score file");
        return {};
    }

    File.clear(); // Clear error flags

    std::ostringstream Buffer;
    Buffer << File.rdbuf(); // Safely read the entire file
    std::string ScoreStr = Buffer.str();

    TSParser *parser = ts_parser_new();
    ts_parser_set_language(parser, tree_sitter_openscofo());

    TSTree *tree = ts_parser_parse_string(parser, nullptr, ScoreStr.c_str(), ScoreStr.size());
    TSNode rootNode = ts_tree_root_node(tree);
    PrintTreeSitterNode(rootNode, ScoreStr);

    return 0;
}
