#include "OpenScofo.hpp"
#include <tree_sitter/api.h>

extern "C" TSLanguage *tree_sitter_openscofo();

namespace OpenScofo {

// ─────────────────────────────────────
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
TSNode Score::GetField(TSNode Node, std::string s) {
    int strLen = s.length();
    TSNode field = ts_node_child_by_field_name(Node, s.c_str(), strLen);
    return field;
}

// ─────────────────────────────────────
std::string Score::GetLuaCode() {
    return m_LuaCode;
}

// ─────────────────────────────────────
void Score::UpdateAudioParameters(float FftSize, float HopSize) {
    m_FFTSize = FftSize;
    m_HopSize = HopSize;
}

// ─────────────────────────────────────
void Score::SetTunning(double Tunning) {
    m_Tunning = Tunning;
}

// ─────────────────────────────────────
double Score::GetFFTSize() {
    return m_FFTSize;
}

// ─────────────────────────────────────
double Score::GetHopSize() {
    return m_HopSize;
}

// ─────────────────────────────────────
double Score::GetPitchTemplateSigma() {
    return m_PitchTemplateSigma;
}

// ─────────────────────────────────────
double Score::GetPhaseCoupling() {
    return m_PhaseCoupling;
}

// ─────────────────────────────────────
double Score::GetSyncStrength() {
    return m_SyncStrength;
}

// ─────────────────────────────────────
bool Score::ScoreIsLoaded() {
    return m_ScoreLoaded;
}

// ─────────────────────────────────────
bool Score::isNumber(std::string str) {
    if (str.empty()) {
        return false;
    }

    if (std::isspace(static_cast<unsigned char>(str[0]))) {
        return false;
    }

    const char *start = str.c_str();
    char *endptr;
    errno = 0;

    std::strtof(start, &endptr);

    if (endptr == start || errno == ERANGE || endptr != start + str.size()) {
        return false;
    }

    return true;
}

// ─────────────────────────────────────
void Score::PitchNode2Freq(const std::string ScoreStr, TSNode node, AudioState &State) {
    TSNode pitch = node;
    std::string type = ts_node_type(pitch);
    if (type == "midi") {
        int midi = std::stof(GetCodeStr(ScoreStr, pitch));
        State.Midi = midi;
        State.Freq = m_Tunning * pow(2, (midi - 69.0) / 12);
        State.Type = PITCH;
        return;
    } else if (type != "pitch") {
        TSPoint Pos = ts_node_start_point(pitch);
        spdlog::error("Invalid pitch type on line {}", std::to_string(Pos.row + 1));
        return;
    }

    std::string pitchNameStr = GetChildStringFromField(ScoreStr, pitch, "pitch_name");
    if (pitchNameStr.empty()) {
        pitchNameStr = GetChildStringFromField(ScoreStr, pitch, "pitchname");
    }
    if (pitchNameStr.empty()) {
        pitchNameStr = GetChildStringFromField(ScoreStr, pitch, "noteName");
    }

    std::string octave = GetChildStringFromField(ScoreStr, pitch, "octave");
    if (pitchNameStr.empty() || octave.empty()) {
        TSPoint Pos = ts_node_start_point(pitch);
        spdlog::error("Invalid pitch on line {}", std::to_string(Pos.row + 1));
        return;
    }

    char pitchName = static_cast<char>(std::toupper(static_cast<unsigned char>(pitchNameStr[0])));
    std::string alt = GetChildStringFromField(ScoreStr, pitch, "alteration");

    int classNote = -1;
    switch (pitchName) {
    case 'C':
        classNote = 0;
        break;
    case 'D':
        classNote = 2;
        break;
    case 'E':
        classNote = 4;
        break;
    case 'F':
        classNote = 5;
        break;
    case 'G':
        classNote = 7;
        break;
    case 'A':
        classNote = 9;
        break;
    case 'B':
        classNote = 11;
        break;
    default:
        TSPoint Pos = ts_node_start_point(pitch);
        spdlog::error("Invalid note name on line line {}", std::to_string(Pos.row + 1));
        return;
    }

    if (alt != "") {
        if (alt == "#") {
            classNote++;
        } else if (alt == "b") {
            classNote--;
        } else if (alt == "##") {
            classNote += 2;
        } else if (alt == "bb") {
            classNote -= 2;
        }
    }

    int midi = classNote + 12 + (12 * std::stoi(octave));
    midi = midi + m_Transpose;
    State.Midi = midi;
    State.Freq = m_Tunning * pow(2, (midi - 69.0) / 12);
    State.Type = PITCH;
}

// ─────────────────────────────────────
double Score::ModPhases(double Phase) {
    Phase = std::fmod(Phase + 0.5, 1.0);
    if (Phase < 0.0) {
        Phase += 1.0;
    }
    return Phase - 0.5;
}

// ╭─────────────────────────────────────╮
// │       Parse File of the Score       │
// ╰─────────────────────────────────────╯
std::string Score::GetCodeStr(const std::string &ScoreStr, TSNode Node) {
    int start = ts_node_start_byte(Node);
    int end = ts_node_end_byte(Node);
    return std::string(std::string_view(ScoreStr.data() + start, end - start));
}

// ─────────────────────────────────────
double Score::GetDurationFromNode(const std::string &ScoreStr, TSNode Node) {
    std::string dur_type = ts_node_type(Node);
    if (dur_type == "number") {
        std::string dur_str = GetCodeStr(ScoreStr, Node);
        return std::stof(dur_str);
    }

    uint32_t count = ts_node_child_count(Node);
    if (count == 1) {
        TSNode dur = ts_node_child(Node, 0);
        dur_type = ts_node_type(dur);
        if (dur_type == "number") {
            std::string dur_str = GetCodeStr(ScoreStr, dur);
            return std::stof(dur_str);
        }
    }

    TSPoint Pos = ts_node_start_point(Node);
    spdlog::error("Invalid duration type on line {}", Pos.row + 1);
    return 0;
}

// ─────────────────────────────────────
MarkovState Score::AddDummySilence() {
    MarkovState Event;
    Event.HSMMType = MARKOV;
    Event.Type = REST;
    Event.ScorePos = m_ScorePosition;
    Event.Index = m_ScoreStates.size();

    AudioState Silence;
    Silence.Type = SILENCE;
    Silence.Freq = 0;
    Silence.Midi = 0;
    Silence.Index = 0;

    Event.AudioStates.emplace_back(Silence);
    return Event;
}

// ─────────────────────────────────────
MarkovState Score::GetFirstEvent() {
    MarkovState Event;
    Event.HSMMType = MARKOV;
    Event.Type = REST;
    Event.ScorePos = 0;
    Event.Index = m_ScoreStates.size();

    AudioState Silence;
    Silence.Type = SILENCE;
    Silence.Freq = 0;
    Silence.Midi = 0;
    Silence.Index = 0;

    Event.AudioStates.emplace_back(Silence);

    return Event;
}

// ─────────────────────────────────────
MarkovState Score::NewPitchEvent(const std::string &ScoreStr, TSNode Node) {
    m_ScorePosition++;

    MarkovState Event;
    Event.Line = ts_node_start_point(Node).row + 1;
    Event.HSMMType = SEMIMARKOV;
    Event.Index = m_ScoreStates.size();
    Event.ScorePos = m_ScorePosition;

    if (ts_node_has_error(Node)) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Pitch event with syntax error on line {}", Init.row + 1);
        return {};
    }

    TSNode PitchNode = ts_node_child_by_field_name(Node, "pitch", 5);
    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);
    TSNode AttributeNode = ts_node_child_by_field_name(Node, "attribute", 9);

    bool Percussive = false;
    if (!ts_node_is_null(AttributeNode)) {
        std::string attr_type = GetChildStringFromField(ScoreStr, AttributeNode, "type");
        Percussive = (attr_type == "percussive");
    }

    if (ts_node_is_null(PitchNode) || ts_node_is_null(DurationNode)) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Invalid NOTE event structure on line {}", Init.row + 1);
        return {};
    }

    if (std::string(ts_node_type(PitchNode)) != "pitch" || std::string(ts_node_type(DurationNode)) != "number") {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Unexpected NOTE event tokens on line {}", Init.row + 1);
        return {};
    }

    Event.Type = NOTE;

    // Pitch
    AudioState SubState;
    PitchNode2Freq(ScoreStr, PitchNode, SubState);
    Event.AudioStates.push_back(SubState);

    // Silence
    if (Percussive) {
        AudioState PercussiveDesc;
        PercussiveDesc.Type = SILENCE;
        Event.AudioStates.push_back(PercussiveDesc);
    }

    // Duration
    double duration = GetDurationFromNode(ScoreStr, DurationNode);
    Event.Duration = duration;

    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
MarkovState Score::NewMultiPitchEvent(const std::string &ScoreStr, TSNode Node) {
    m_ScorePosition++;

    MarkovState Event;
    Event.Line = ts_node_start_point(Node).row + 1;
    Event.HSMMType = SEMIMARKOV;
    Event.Index = m_ScoreStates.size();
    Event.ScorePos = m_ScorePosition;

    if (ts_node_has_error(Node)) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Multi pitch event with syntax error on line {}", Init.row + 1);
        return {};
    }

    TSNode PitchesNode = ts_node_child_by_field_name(Node, "pitches", 7);
    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);

    if (ts_node_is_null(PitchesNode) || ts_node_is_null(DurationNode)) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Invalid multi pitch event structure on line {}", Init.row + 1);
        return {};
    }

    std::string nodeType = ts_node_type(Node);
    if (nodeType == "trill_event") {
        Event.Type = TRILL;
    } else if (nodeType == "chord_event") {
        Event.Type = CHORD;
    } else {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Wrong Type on {}", Init.row + 1);
        return {};
    }

    uint32_t pitch_count = ts_node_named_child_count(PitchesNode);
    if (pitch_count == 0) {
        TSPoint Init = ts_node_start_point(PitchesNode);
        spdlog::error("Missing pitches on line {}", Init.row + 1);
        return {};
    }

    for (uint32_t i = 0; i < pitch_count; i++) {
        TSNode PitchNode = ts_node_named_child(PitchesNode, i);
        if (std::string(ts_node_type(PitchNode)) != "pitch") {
            continue;
        }
        AudioState SubState;
        PitchNode2Freq(ScoreStr, PitchNode, SubState);
        Event.AudioStates.push_back(SubState);
    }

    double duration = GetDurationFromNode(ScoreStr, DurationNode);
    Event.Duration = duration;

    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
MarkovState Score::NewRestEvent(const std::string &ScoreStr, TSNode Node) {
    // m_ScorePosition++;
    // Note that Rest do not count for the score position, as they are not considered in the evaluation

    if (m_ScorePosition == 0) {
        spdlog::warn("OpenScofo cannot detect the start of a piece when the first events are REST. "
                     "It cannot distinguish between silence before the piece and the actual start of the piece. "
                     "As a result, the current event (line {}) and its associated actions will not be added.",
                     ts_node_start_point(Node).row + 1);
        return {};
    }

    MarkovState Event;
    Event.Line = ts_node_start_point(Node).row + 1;
    Event.HSMMType = SEMIMARKOV;
    Event.Index = m_ScoreStates.size();
    Event.ScorePos = m_ScorePosition;

    if (ts_node_has_error(Node)) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Rest event with syntax error on line {}", Init.row + 1);
        return {};
    }

    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);

    if (ts_node_is_null(DurationNode)) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Invalid REST event structure on line {}", Init.row + 1);
        return {};
    }

    double duration = GetDurationFromNode(ScoreStr, DurationNode);
    Event.Duration = duration;
    Event.Type = REST;
    m_MinimalDuration = std::min(m_MinimalDuration, duration);

    AudioState Silence;
    Silence.Type = SILENCE;
    Event.AudioStates.push_back(Silence);

    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
void Score::ProcessEventTime(MarkovState &Event) {
    if (Event.Index != 0) {
        int index = Event.Index;
        MarkovState &prev = m_ScoreStates[index - 1];

        double psiPrev = 60.0f / prev.BPMExpected;
        double ioibeats = prev.Duration;

        Event.OnsetExpected = prev.OnsetExpected + ioibeats * psiPrev;
        Event.IOIHatPhiN = ModPhases(prev.IOIHatPhiN + ioibeats);
        Event.PhaseExpected = Event.IOIHatPhiN;
        Event.IOIPhiN = Event.IOIHatPhiN;
    } else {
        Event.PhaseExpected = 0;
        Event.IOIHatPhiN = 0;
        Event.IOIPhiN = 0;
        Event.OnsetExpected = 0;
    }
    Event.PhaseCoupling = m_PhaseCoupling;
    Event.SyncStrength = m_SyncStrength;
    Event.BPMExpected = m_CurrentBPM;

    spdlog::debug("Added Time for Event {}, BPM {}, Phase PhaseCoupling {}, SyncStrength {}, "
                  "PhaseExpected {}, Onset Expected {}",
                  Event.ScorePos, Event.BPMExpected, Event.PhaseCoupling, Event.SyncStrength, Event.PhaseExpected,
                  Event.OnsetExpected);
}

// ─────────────────────────────────────
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
void Score::NewEvent(const std::string &ScoreStr, TSNode Node) {
    MarkovState Event;

    TSNode definition = GetField(Node, "definition");
    if (ts_node_is_null(definition)) {
        TSPoint Pos = ts_node_start_point(Node);
        spdlog::error("Invalid EVENT on line {}", Pos.row + 1);
        return;
    }

    std::string defType = ts_node_type(definition);
    if (defType == "note_event") {
        Event = NewPitchEvent(ScoreStr, definition);
    } else if (defType == "chord_event" || defType == "trill_event") {
        Event = NewMultiPitchEvent(ScoreStr, definition);
    } else if (defType == "rest_event") {
        Event = NewRestEvent(ScoreStr, definition);
    } else {
        TSPoint Pos = ts_node_start_point(definition);
        spdlog::error("Type not implemented {} on line {}.", defType, Pos.row + 1);
        return;
    }

    if (Event.AudioStates.empty()) {
        return;
    }

    uint32_t child_count = ts_node_child_count(Node);
    for (uint32_t i = 0; i < child_count; i++) {
        TSNode child = ts_node_child(Node, i);
        std::string type = ts_node_type(child);
        if (type == "action") {
            NewEventAction(ScoreStr, child, Event);
        }
    }

    m_PrevDuration = Event.Duration;
    m_LastOnset = Event.OnsetExpected;
    m_ScoreStates.emplace_back(Event);
}

// ─────────────────────────────────────
std::string GetChildStringFromType(const std::string &source, TSNode parent, const std::string &wanted_type) {
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
void Score::NewConfig(const std::string &ScoreStr, TSNode node) {
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

    if (id == "BPM" || id == "TRANSPOSE" || id == "PHASECOUPLING" || id == "SYNCSTRENGTH" ||
        id == "PITCHTEMPLATESIGMA" || id == "FFTSIZE" || id == "HOPSIZE") {
        if (valueType != "number") {
            spdlog::error("Invalid numeric value for {} on line {}.", id, pos.row + 1);
            return;
        }

        float v = std::stof(value);
        if (id == "BPM") {
            m_CurrentBPM = v;
            MarkovState Begin = GetFirstEvent();
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
                m_PhaseCoupling = v;
            }
        } else if (id == "SYNCSTRENGTH") {
            if (v < 0 || v > 1) {
                spdlog::error("Invalid value for SYNCSTRENGTH on line {}.", pos.row + 1);
            } else {
                m_SyncStrength = v;
            }
        } else if (id == "PITCHTEMPLATESIGMA") {
            if (v < 0 || v > 1) {
                spdlog::error("Invalid value for PITCHTEMPLATESIGMA on line {}.", pos.row + 1);
            } else {
                m_PitchTemplateSigma = v;
            }
        } else if (id == "FFTSIZE") {
            int fft = static_cast<int>(v);
            if (fft > 0 && (fft & (fft - 1)) == 0) {
                m_FFTSize = fft;
            } else {
                spdlog::error("FFTSIZE must be a power of two.");
            }
        } else if (id == "HOPSIZE") {
            int hop = static_cast<int>(v);
            if (hop > 0 && (hop & (hop - 1)) == 0) {
                m_HopSize = hop;
            } else {
                spdlog::error("HOPSIZE must be a power of two.");
            }
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

        m_TimbreModel = m_ScoreRootPath / fs::path(path);
        if (!fs::exists(m_TimbreModel)) {
            spdlog::error("Model path not found: {}", m_TimbreModel.string());
        }
        return;
    }
}

// ─────────────────────────────────────
void Score::NewEventAction(const std::string &ScoreStr, TSNode Node, MarkovState &Event) {
    Action BaseAction;
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

        Action NewAction = BaseAction;

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
                        token.erase(std::remove_if(token.begin(), token.end(), [](unsigned char c) {
                            return std::isspace(c);
                        }),
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
                            NewAction.Args.push_back(std::stof(token));
                        } else {
                            spdlog::error("Invalid number argument on line {}.", ts_node_start_point(pdarg).row + 1);
                            return;
                        }
                    } else if (pdargType == "identifier" || pdargType == "symbol") {
                        NewAction.Args.push_back(token);
                    } else if (!token.empty()) {
                        NewAction.Args.push_back(token);
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

// ─────────────────────────────────────
void FindErrors(TSNode &root, TSNode &node, const std::string &source_code) {
    if (!ts_node_is_null(node) && !ts_node_eq(root, node) && ts_node_has_error(node)) {
        TSPoint start = ts_node_start_point(node);
        uint32_t row = start.row + 1; // 1-based
        uint32_t column = start.column + 1;
        uint32_t byteStart = ts_node_start_byte(node);
        uint32_t byteEnd = ts_node_end_byte(node);
        std::string tokenText = source_code.substr(byteStart, byteEnd - byteStart);
        spdlog::error("Fond Error {}, line: {}, column: {}, token: '{}'", ts_node_type(node), row, column, tokenText);
    }

    uint32_t childCount = ts_node_child_count(node);
    for (uint32_t i = 0; i < childCount; i++) {
        TSNode child = ts_node_child(node, i);
        FindErrors(root, child, source_code);
    }
}

// ─────────────────────────────────────
bool ScoreIsText(const std::string &path) {
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

// ─────────────────────────────────────
States Score::Parse(std::string ScoreFile) {
    m_ScoreStates.clear();
    m_LuaCode.clear();

    fs::path ScoreFilePath = fs::path(ScoreFile);
    if (fs::exists(ScoreFilePath) == false) {
        spdlog::error("Score File not found");
        return {};
    }
    m_ScoreRootPath = ScoreFilePath.parent_path();

    // Open the score file for reading
    std::ifstream File(ScoreFile, std::ios::binary);
    if (File.is_open() == false) {
        spdlog::error("Not possible to open score file");
        return {};
    }

    File.clear(); // Clear error flags

    std::ostringstream Buffer;
    Buffer << File.rdbuf(); // Safely read the entire file
    std::string ScoreStr = Buffer.str();

    // Proceed with parsing ScoreStr...
    // Config Values
    m_CurrentBPM = -1;
    m_Transpose = 0;
    m_PitchTemplateSigma = 0.5;

    m_LineCount = 0;
    m_MarkovIndex = 0;
    m_ScorePosition = 0;
    m_LastOnset = 0;
    m_PrevDuration = 0;
    std::string Line;

    // read and process score
    TSParser *parser = ts_parser_new();
    ts_parser_set_language(parser, tree_sitter_openscofo());
    TSTree *tree = ts_parser_parse_string(parser, nullptr, ScoreStr.c_str(), ScoreStr.size());
    TSNode rootNode = ts_tree_root_node(tree);

    if (ts_node_has_error(rootNode)) {
        FindErrors(rootNode, rootNode, ScoreStr);
    }

    uint32_t child_count = ts_node_child_count(rootNode);
    for (uint32_t i = 0; i < child_count; i++) {
        TSNode child = ts_node_child(rootNode, i);
        std::string type = ts_node_type(child);
        if (type == "EVENT") {
            if (m_CurrentBPM == -1) {
                spdlog::error("BPM is not defined");
                return {};
            }
            NewEvent(ScoreStr, child);
        } else if (type == "CONFIG") {
            NewConfig(ScoreStr, child);
        } else if (type == "LUA") {
            std::string lua_body = GetChildStringFromField(ScoreStr, child, "lua_body");
            lua_body += "\n\n";
            m_LuaCode += lua_body;
        } else if (type == "comment") {
        } else {
            spdlog::error("Not recognized {}", type);
        }
    }

    // Cleanup
    ts_tree_delete(tree);
    ts_parser_delete(parser);

    m_ScoreLoaded = true;
    return m_ScoreStates;
}
} // namespace OpenScofo
