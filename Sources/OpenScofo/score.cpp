/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#include "OpenScofo.hpp"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <tree_sitter/api.h>

namespace OpenScofo {

extern "C" TSLanguage *tree_sitter_openscofo();

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
bool Score::ScoreIsLoaded() {
    return m_ScoreLoaded;
}

// ─────────────────────────────────────
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

// ─────────────────────────────────────
void Score::PitchNode2Freq(const std::string ScoreStr, TSNode node, Observation &State) {
    TSNode pitch = node;
    std::string type = ts_node_type(pitch);
    TSNode midiNode = ts_node_child_by_field_name(pitch, "midi", 4);

    if (type == "midi" || !ts_node_is_null(midiNode)) {
        const std::string midiText = GetCodeStr(ScoreStr, type == "midi" ? pitch : midiNode);
        double parsedMidi = 0.0;

        if (!ParseDouble(midiText, parsedMidi)) {
            spdlog::error("Invalid MIDI pitch on line {}", ts_node_start_point(pitch).row + 1);
            return;
        }

        const double midi = parsedMidi + m_Transpose;
        const double frequency = m_Tunning * std::pow(2.0, (midi - 69.0) / 12.0);
        if (!std::isfinite(midi) || !std::isfinite(frequency) || frequency <= 0.0) {
            spdlog::error("MIDI pitch out of range on line {}", ts_node_start_point(pitch).row + 1);
            return;
        }

        State.Midi = midi;
        State.Freq = frequency;
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
        } else if (alt == "+") {
            classNote += 0.5;
        } else if (alt == "#+") {
            classNote += 1.5;
        } else if (alt == "-") {
            classNote -= 0.5;
        } else if (alt == "b-") {
            classNote -= 1.5;
        }
    }

    float midi = classNote + 12 + (12 * std::stoi(octave));
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
void Score::AddDummySilence(const ScoreState &Next) {
    if (m_ScoreStates.empty() || Next.Type == REST || Next.Type == FIRSTEVENT) {
        return;
    }
    const ScoreState &Previous = m_ScoreStates.back();
    const bool Sounded =
        Previous.Type == NOTE || Previous.Type == CHORD || Previous.Type == PTECH || Previous.Type == UTECH;
    if (!Sounded || Previous.Section != Next.Section) {
        return;
    }

    ScoreState Event{};
    Event.HSMMType = MARKOV;
    Event.Type = REST;
    Event.IsInterEventSilence = true;
    Event.ScorePos = Previous.ScorePos;
    Event.Index = static_cast<int>(m_ScoreStates.size());
    Event.Duration = 0;
    Event.Section = Previous.Section;
    Event.Line = Previous.Line;
    Event.BPMExpected = Previous.BPMExpected;
    Event.OnsetExpected = Next.OnsetExpected;
    Event.PhaseExpected = Next.PhaseExpected;
    Event.IOIHatPhiN = Next.IOIHatPhiN;
    Event.IOIPhiN = Next.IOIPhiN;
    Event.SyncStrength = Previous.SyncStrength;
    Event.PhaseCoupling = Previous.PhaseCoupling;
    Event.TimeTolerance = Previous.TimeTolerance;
    Event.Observations.push_back({SILENCE});

    m_ScoreStates.emplace_back(std::move(Event));
}

// ─────────────────────────────────────
ScoreState Score::GetFirstEvent() {
    ScoreState Event;
    Event.HSMMType = MARKOV;
    Event.Type = FIRSTEVENT;
    Event.ScorePos = 0;
    Event.Index = m_ScoreStates.size();
    Event.Duration = 0.0;

    Observation Silence;
    Silence.Type = SILENCE;
    Silence.Freq = 0;
    Silence.Midi = 0;
    Silence.Index = 0;

    Event.Observations.emplace_back(Silence);

    return Event;
}

// ─────────────────────────────────────
ScoreState Score::NewPitchEvent(const std::string &ScoreStr, TSNode Node) {
    m_ScorePosition++;

    ScoreState Event;
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
    Observation SubState{};
    PitchNode2Freq(ScoreStr, PitchNode, SubState);
    if (!(SubState.Freq > 0.0)) {
        return {};
    }
    Event.Observations.push_back(SubState);

    if (Percussive) {
        // TODO: need tests
        Observation Onset;
        Onset.Type = ONSET;
        Event.Observations.push_back(Onset);
    }

    // Duration
    double duration = GetDurationFromNode(ScoreStr, DurationNode);
    Event.Duration = duration;

    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
ScoreState Score::NewChordEvent(const std::string &ScoreStr, TSNode Node) {
    m_ScorePosition++;

    ScoreState Event{};
    Event.Line = ts_node_start_point(Node).row + 1;
    Event.HSMMType = SEMIMARKOV;
    Event.Type = CHORD;
    Event.MicroTopologyType = NO_MICROSTATES;
    Event.Index = m_ScoreStates.size();
    Event.ScorePos = m_ScorePosition;

    TSNode PitchesNode = ts_node_child_by_field_name(Node, "pitches", 7);
    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);
    if (ts_node_has_error(Node) || ts_node_is_null(PitchesNode) || ts_node_is_null(DurationNode)) {
        spdlog::error("Invalid CHORD event structure on line {}", Event.Line);
        return {};
    }
    const uint32_t PitchCount = ts_node_named_child_count(PitchesNode);
    if (PitchCount == 0) {
        spdlog::error("Missing CHORD pitches on line {}", Event.Line);
        return {};
    }

    for (uint32_t i = 0; i < PitchCount; ++i) {
        TSNode PitchNode = ts_node_named_child(PitchesNode, i);
        Observation Pitch{};
        PitchNode2Freq(ScoreStr, PitchNode, Pitch);
        if (!(Pitch.Freq > 0.0)) {
            return {};
        }
        Event.Observations.push_back(Pitch);
    }

    Event.Duration = GetDurationFromNode(ScoreStr, DurationNode);
    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
ScoreState Score::NewTrillEvent(const std::string &ScoreStr, TSNode Node) {
    m_ScorePosition++;

    ScoreState Event{};
    Event.Line = ts_node_start_point(Node).row + 1;
    Event.HSMMType = SEMIMARKOV;
    Event.Type = TRILL;
    Event.MicroTopologyType = UNORDERED;
    Event.Index = m_ScoreStates.size();
    Event.ScorePos = m_ScorePosition;

    TSNode PitchesNode = ts_node_child_by_field_name(Node, "pitches", 7);
    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);
    if (ts_node_has_error(Node) || ts_node_is_null(PitchesNode) || ts_node_is_null(DurationNode)) {
        spdlog::error("Invalid TRILL event structure on line {}", Event.Line);
        return {};
    }
    const uint32_t PitchCount = ts_node_named_child_count(PitchesNode);
    if (PitchCount == 0) {
        spdlog::error("Missing TRILL pitches on line {}", Event.Line);
        return {};
    }

    for (uint32_t i = 0; i < PitchCount; ++i) {
        TSNode PitchNode = ts_node_named_child(PitchesNode, i);
        Observation Pitch{};
        PitchNode2Freq(ScoreStr, PitchNode, Pitch);
        if (!(Pitch.Freq > 0.0)) {
            return {};
        }
        MarkovMicroState MicroState;
        MicroState.Observations.push_back(Pitch);
        Event.MicroStates.push_back(std::move(MicroState));
    }

    Event.Duration = GetDurationFromNode(ScoreStr, DurationNode);
    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
ScoreState Score::NewMultiEvent(const std::string &ScoreStr, TSNode Node) {
    m_ScorePosition++;

    ScoreState Event{};
    Event.Line = ts_node_start_point(Node).row + 1;
    Event.HSMMType = SEMIMARKOV;
    Event.Type = GLISS;
    Event.MicroTopologyType = LEFT_RIGHT;
    Event.Index = m_ScoreStates.size();
    Event.ScorePos = m_ScorePosition;

    TSNode PitchesNode = ts_node_child_by_field_name(Node, "pitches", 7);
    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);
    if (ts_node_has_error(Node) || ts_node_is_null(PitchesNode) || ts_node_is_null(DurationNode)) {
        spdlog::error("Invalid GLISS event structure on line {}", Event.Line);
        return {};
    }
    const uint32_t PitchCount = ts_node_named_child_count(PitchesNode);
    if (PitchCount == 0) {
        spdlog::error("Missing GLISS pitches on line {}", Event.Line);
        return {};
    }

    // Build the complete glissando here. The forward model only follows
    // this ordered chain; it never generates intermediate pitches.
    auto AddPitch = [&](const Observation &Pitch) {
        MarkovMicroState MicroState;
        MicroState.Observations.push_back(Pitch);
        Event.MicroStates.push_back(std::move(MicroState));
    };

    for (uint32_t i = 0; i < PitchCount; ++i) {
        TSNode PitchNode = ts_node_named_child(PitchesNode, i);
        Observation Pitch{};
        PitchNode2Freq(ScoreStr, PitchNode, Pitch);
        if (!(Pitch.Freq > 0.0)) {
            return {};
        }
        if (!Event.MicroStates.empty()) {
            const double PreviousMidi = Event.MicroStates.back().Observations[0].Midi;
            if (Pitch.Midi == PreviousMidi) {
                continue;
            }
            const double Step = Pitch.Midi > PreviousMidi ? 0.5 : -0.5;
            for (double Midi = PreviousMidi + Step; Step > 0 ? Midi < Pitch.Midi : Midi > Pitch.Midi; Midi += Step) {
                Observation Intermediate{};
                Intermediate.Type = PITCH;
                Intermediate.Midi = Midi;
                Intermediate.Freq = m_Tunning * std::pow(2.0, (Midi - 69.0) / 12.0);
                AddPitch(Intermediate);
            }
        }
        AddPitch(Pitch);
    }

    Event.Duration = GetDurationFromNode(ScoreStr, DurationNode);
    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
ScoreState Score::NewPTechEvent(const std::string &ScoreStr, TSNode Node) {
    m_ScorePosition++;

    ScoreState Event;
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

    Event.Type = PTECH;
    Event.MicroTopologyType = UNORDERED;

    // Sounded alternatives: technique labels or pitch.
    Event.MicroStates.resize(2);
    bool HasTechnique = false;
    TSNode TechniquesNode = ts_node_child_by_field_name(Node, "techniques", 10);

    if (!ts_node_is_null(TechniquesNode)) {
        uint32_t count = ts_node_named_child_count(TechniquesNode);
        for (uint32_t i = 0; i < count; ++i) {
            TSNode TechniqueNode = ts_node_named_child(TechniquesNode, i);
            if (std::string(ts_node_type(TechniqueNode)) != "identifier") {
                continue;
            }

            Observation SubState;
            SubState.Label = GetCodeStr(ScoreStr, TechniqueNode);
            SubState.Type = LABEL;
            Event.MicroStates[0].Observations.push_back(SubState);
            HasTechnique = true;
        }
    } else {
        std::string Label = GetChildStringFromField(ScoreStr, Node, "technique");
        if (!Label.empty()) {
            Observation SubState;
            SubState.Label = Label;
            SubState.Type = LABEL;
            Event.MicroStates[0].Observations.push_back(SubState);
            HasTechnique = true;
        }
    }

    if (!HasTechnique) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("PTECH event without technique on line {}", Init.row + 1);
        return {};
    }

    // Pitch
    Observation Pitch{};
    PitchNode2Freq(ScoreStr, PitchNode, Pitch);
    if (!(Pitch.Freq > 0.0)) {
        return {};
    }
    Event.MicroStates[1].Observations.push_back(Pitch);

    // Duration
    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);
    double Duration = GetDurationFromNode(ScoreStr, DurationNode);
    Event.Duration = Duration;

    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
ScoreState Score::NewUTechEvent(const std::string &ScoreStr, TSNode Node) {
    m_ScorePosition++;

    ScoreState Event;
    Event.Line = ts_node_start_point(Node).row + 1;
    Event.HSMMType = SEMIMARKOV;
    Event.Index = m_ScoreStates.size();
    Event.ScorePos = m_ScorePosition;

    if (ts_node_has_error(Node)) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("Pitch event with syntax error on line {}", Init.row + 1);
        return {};
    }

    Event.Type = UTECH;
    Event.MicroTopologyType = UNORDERED;

    // Alternative technique labels share one sounded microstate.
    Event.MicroStates.resize(1);
    bool HasTechnique = false;
    TSNode TechniquesNode = ts_node_child_by_field_name(Node, "techniques", 10);
    if (!ts_node_is_null(TechniquesNode)) {
        uint32_t count = ts_node_named_child_count(TechniquesNode);
        for (uint32_t i = 0; i < count; ++i) {
            TSNode TechniqueNode = ts_node_named_child(TechniquesNode, i);
            if (std::string(ts_node_type(TechniqueNode)) != "identifier") {
                continue;
            }
            Observation SubState;
            SubState.Label = GetCodeStr(ScoreStr, TechniqueNode);
            SubState.Type = LABEL;
            Event.MicroStates[0].Observations.push_back(SubState);
            HasTechnique = true;
        }
    } else {
        std::string Label = GetChildStringFromField(ScoreStr, Node, "technique");
        if (!Label.empty()) {
            Observation SubState;
            SubState.Label = Label;
            SubState.Type = LABEL;
            Event.MicroStates[0].Observations.push_back(SubState);
            HasTechnique = true;
        }
    }

    if (!HasTechnique) {
        TSPoint Init = ts_node_start_point(Node);
        spdlog::error("UTECH event without technique on line {}", Init.row + 1);
        return {};
    }

    // Duration
    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);
    double Duration = GetDurationFromNode(ScoreStr, DurationNode);
    Event.Duration = Duration;

    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
ScoreState Score::NewRestEvent(const std::string &ScoreStr, TSNode Node) {
    if (m_ScorePosition == 0) {
        spdlog::warn("OpenScofo cannot detect the start of a piece when the first events are REST. "
                     "It cannot distinguish between silence before the piece and the actual start of the piece. "
                     "As a result, the first event (line {}) and its associated actions will not be added.",
                     ts_node_start_point(Node).row + 1);
        return {};
    }

    ScoreState Event;
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

    Observation Silence;
    Silence.Type = SILENCE;
    Event.Observations.push_back(Silence);

    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
void Score::ProcessEventTime(ScoreState &Event) {
    Event.Section = m_CurrentSection;

    const bool IsFirstStateInSection =
        Event.Index == 0 || m_ScoreStates[static_cast<size_t>(Event.Index - 1)].Section != Event.Section;
    if (!IsFirstStateInSection) {
        int index = Event.Index;
        ScoreState &prev = m_ScoreStates[index - 1];

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

    Event.BPMExpected = m_CurrentBPM;

    spdlog::debug("Added Time for Event {}, BPM {}, Phase PhaseCoupling {}, SyncStrength {}, "
                  "PhaseExpected {}, Onset Expected {}",
                  Event.ScorePos, Event.BPMExpected, Event.PhaseCoupling, Event.SyncStrength, Event.PhaseExpected,
                  Event.OnsetExpected);
}

// ─────────────────────────────────────
void Score::NewSection(const std::string &ScoreStr, TSNode Node) {
    TSNode NameNode = GetField(Node, "name");
    TSPoint Position = ts_node_start_point(Node);
    if (ts_node_is_null(NameNode)) {
        spdlog::error("Invalid SECTION on line {}.", Position.row + 1);
        return;
    }

    std::string Name = GetCodeStr(ScoreStr, NameNode);
    if (Name.size() >= 2 && Name.front() == '"' && Name.back() == '"') {
        Name = Name.substr(1, Name.size() - 2);
    }
    if (Name.empty()) {
        spdlog::error("SECTION name cannot be empty on line {}.", Position.row + 1);
        return;
    }

    const bool HasMusicalEvent = std::any_of(m_ScoreStates.begin(), m_ScoreStates.end(),
                                             [](const ScoreState &State) { return State.Type != FIRSTEVENT; });

    // BPM is commonly declared before the first SECTION. In that case its
    // leading FIRSTEVENT is the boundary of the first section, not an
    // unsectioned state.
    if (!m_HasSection && !HasMusicalEvent) {
        for (ScoreState &State : m_ScoreStates) {
            if (State.Type == FIRSTEVENT && State.Section.empty()) {
                State.Section = Name;
            }
        }
    }

    // At the first musical event this either refreshes an existing BPM-created
    // FIRSTEVENT or inserts a boundary when the section inherits its BPM.
    m_SectionStartPending = true;
    m_HasSection = true;
    m_CurrentSection = std::move(Name);
}

// ─────────────────────────────────────
void Score::EnsureSectionStart(TSNode EventNode, Configuration &Config) {
    if (!m_SectionStartPending) {
        return;
    }

    ScoreState *SectionStart = nullptr;
    if (!m_ScoreStates.empty() && m_ScoreStates.back().Type == FIRSTEVENT &&
        m_ScoreStates.back().Section == m_CurrentSection) {
        SectionStart = &m_ScoreStates.back();
    } else {
        ScoreState Begin = GetFirstEvent();
        Begin.Line = static_cast<int>(ts_node_start_point(EventNode).row + 1);
        Begin.TimeTolerance = m_TimeTolerance;
        Begin.SyncStrength = Config.SyncStrength;
        Begin.PhaseCoupling = Config.PhaseCoupling;
        ProcessEventTime(Begin);
        m_ScoreStates.emplace_back(std::move(Begin));
        SectionStart = &m_ScoreStates.back();
    }

    SectionStart->BPMExpected = m_CurrentBPM;
    SectionStart->TimeTolerance = m_TimeTolerance;
    SectionStart->SyncStrength = Config.SyncStrength;
    SectionStart->PhaseCoupling = Config.PhaseCoupling;
    m_SectionStartPending = false;
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
void Score::NewEvent(const std::string &ScoreStr, TSNode Node, Configuration &Config) {
    EnsureSectionStart(Node, Config);
    ScoreState Event;

    TSNode definition = GetField(Node, "definition");
    if (ts_node_is_null(definition)) {
        TSPoint Pos = ts_node_start_point(Node);
        spdlog::error("Invalid EVENT on line {}", Pos.row + 1);
        return;
    }

    std::string defType = ts_node_type(definition);
    if (defType == "note_event") {
        Event = NewPitchEvent(ScoreStr, definition);
    } else if (defType == "chord_event") {
        Event = NewChordEvent(ScoreStr, definition);
    } else if (defType == "trill_event") {
        Event = NewTrillEvent(ScoreStr, definition);
    } else if (defType == "multi_event") {
        Event = NewMultiEvent(ScoreStr, definition);
    } else if (defType == "ptech_event") {
        Event = NewPTechEvent(ScoreStr, definition);
    } else if (defType == "utech_event") {
        Event = NewUTechEvent(ScoreStr, definition);
    } else if (defType == "rest_event") {
        Event = NewRestEvent(ScoreStr, definition);
    } else {
        TSPoint Pos = ts_node_start_point(definition);
        spdlog::error("Type not implemented {} on line {}.", defType, Pos.row + 1);
        return;
    }

    // Configuration by Event
    Event.TimeTolerance = m_TimeTolerance;

    if ((Event.Observations.empty() && Event.MicroStates.empty()) ||
        std::any_of(Event.MicroStates.begin(), Event.MicroStates.end(),
                    [](const MarkovMicroState &MicroState) { return MicroState.Observations.empty(); })) {
        return;
    }

    if (Event.Type == TRILL || Event.Type == GLISS || Event.Type == UTECH || Event.Type == PTECH) {
        const bool Technique = Event.Type == UTECH || Event.Type == PTECH;
        const size_t ExpectedMicroStates = Event.Type == PTECH ? 2 : 1;
        if (!Event.Observations.empty() || Event.MicroStates.empty() ||
            (Technique && Event.MicroStates.size() != ExpectedMicroStates)) {
            spdlog::error("Invalid microstate structure on line {}", Event.Line);
            return;
        }

        for (size_t k = 0; k < Event.MicroStates.size(); ++k) {
            AudioDescType Expected = PITCH;

            if (Event.Type == PTECH) {
                // Unordered technique and pitch alternatives.
                Expected = k == 0 ? LABEL : PITCH;
            } else if (Event.Type == UTECH) {
                Expected = LABEL;
            }

            const auto &Observations = Event.MicroStates[k].Observations;

            if (Expected != LABEL && Observations.size() != 1) {
                spdlog::error("Invalid microstate observations on line {}", Event.Line);
                return;
            }

            for (const Observation &Obs : Observations) {
                if (Obs.Type != Expected || (Expected == LABEL && Obs.Label.empty()) ||
                    (Expected == PITCH && !(Obs.Freq > 0.0))) {

                    spdlog::error("Invalid microstate observation on line {}", Event.Line);
                    return;
                }
            }
        }
    }

    uint32_t child_count = ts_node_child_count(Node);
    for (uint32_t i = 0; i < child_count; i++) {
        TSNode child = ts_node_child(Node, i);
        std::string type = ts_node_type(child);
        if (type == "action") {
            NewEventAction(ScoreStr, child, Event);
        }
    }

    Event.SyncStrength = Config.SyncStrength;
    Event.PhaseCoupling = Config.PhaseCoupling;

    m_PrevDuration = Event.Duration;
    m_LastOnset = Event.OnsetExpected;
    // Event timing is already computed from the previous scored event.
    AddDummySilence(Event);

    Event.Index = static_cast<int>(m_ScoreStates.size());
    m_ScoreStates.emplace_back(std::move(Event));
}

// ─────────────────────────────────────
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

// ─────────────────────────────────────
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

// ─────────────────────────────────────
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

// ─────────────────────────────────────
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
