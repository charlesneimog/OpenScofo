/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

/**
 * @file events.cpp
 * @brief Score event construction, timing, sections, and microstates.
 */

#include "OpenScofo.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <utility>
#include <tree_sitter/api.h>

namespace OpenScofo {

// ─────────────────────────────────────
/**
 * @brief Read the parser loaded flag.
 *
 * @return True if the stored loaded flag is set; false otherwise.
 *
 * @note Set when Parse() reaches its final return; does not independently validate runtime readiness.
 * @warning The flag does not guarantee an error-free score and is not reset at the start of a subsequent parse.
 */
bool Score::GetAttribute(const std::string ScoreStr, TSNode Node, const std::string &attr) {
    TSNode AttributeNode = ts_node_child_by_field_name(Node, "attribute", 9);
    if (!ts_node_is_null(AttributeNode)) {
        std::string attr_type = GetChildStringFromField(ScoreStr, AttributeNode, "type");
        return (attr_type == attr);
    }
    return false;
}

// ─────────────────────────────────────
/**
 * @brief Decode a pitch node into frequency and MIDI pitch.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param node Tree-sitter node to inspect in the current score tree.
 * @param State Audio observation receiving the decoded pitch and MIDI value.
 *
 * @note Applies configured A4 tuning and transposition to MIDI values or note names with accidentals.
 * @warning Node offsets must correspond to ScoreStr; invalid pitches are logged and may leave the output unchanged.
 */
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

    double classNote = -1;
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

    double midi = classNote + 12 + (12 * std::stoi(octave));
    midi = midi + m_Transpose;
    State.Midi = midi;
    State.Freq = m_Tunning * pow(2, (midi - 69.0) / 12);
    State.Type = PITCH;
}

// ─────────────────────────────────────
/**
 * @brief Wrap an event phase to the centered unit interval.
 *
 * @param Phase Phase in cycles.
 *
 * @return Wrapped phase in [-0.5, 0.5).
 *
 * @note Preserves phase equivalence modulo one cycle.
 */
double Score::ModPhases(double Phase) {
    Phase = std::fmod(Phase + 0.5, 1.0);
    if (Phase < 0.0) {
        Phase += 1.0;
    }
    return Phase - 0.5;
}

// ─────────────────────────────────────
/**
 * @brief Read an event duration from a numeric syntax node.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Parsed duration in beats, or zero for an unsupported node structure.
 *
 * @note Accepts a number directly or as the sole child; logs unsupported structures.
 * @warning Numeric conversion can throw for malformed or out-of-range text.
 */
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
/**
 * @brief Insert an optional inter-event silence state when applicable.
 *
 * @param Next Next musical event, used to decide and time an optional silence state.
 *
 * @note Adds a zero-beat Markov REST after supported sounded events in the same section, preserving public score
 * position.
 */
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
/**
 * @brief Construct an initial silence boundary state.
 *
 * @return Initial boundary state with one silence observation.
 *
 * @note Creates a zero-duration Markov FIRSTEVENT; timing and section metadata are assigned separately.
 */
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
/**
 * @brief Construct a NOTE event from a syntax node.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Constructed note event, or a default state for an invalid structure or pitch.
 *
 * @note Advances public score numbering, converts pitch, reads duration, and applies event timing.
 * @warning Requires source-matched nodes; failed construction can still advance the parser score position.
 */
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

    bool Percussive = GetAttribute(ScoreStr, Node, "percussive");
    if (Percussive) {
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
/**
 * @brief Construct a CHORD event from its pitch list and duration.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Constructed chord event, or a default state for invalid structure or pitches.
 *
 * @note Advances public score numbering and places simultaneous pitch observations in one event.
 */
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
/**
 * @brief Construct a TRILL event with unordered pitch alternatives.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Constructed trill event, or a default state for invalid structure or pitches.
 *
 * @note Creates one microstate per listed pitch and applies the parent duration and timing.
 */
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
/**
 * @brief Construct a GLISS event as an ordered pitch microstate chain.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Constructed glissando event, or a default state for invalid structure or pitches.
 *
 * @note Interpolates intermediate pitches in half-semitone steps and skips duplicate adjacent endpoints.
 */
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
/**
 * @brief Construct a pitched extended-technique event.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Constructed PTECH event, or a default state for invalid syntax, missing labels, or invalid pitch.
 *
 * @note Creates unordered technique-label and pitch alternatives with a shared parent duration.
 */
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
    Event.MicroTopologyType = LEFT_RIGHT;

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
/**
 * @brief Construct an unpitched extended-technique event.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Constructed UTECH event, or a default state for invalid syntax or missing technique labels.
 *
 * @note Groups alternative technique labels in a single unordered microstate.
 */
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

    bool Percussive = GetAttribute(ScoreStr, Node, "percussive");
    if (Percussive) {
        Observation Onset;
        Onset.Type = ONSET;
        Event.Observations.push_back(Onset);
    }

    // Duration
    TSNode DurationNode = ts_node_child_by_field_name(Node, "duration", 8);
    double Duration = GetDurationFromNode(ScoreStr, DurationNode);
    Event.Duration = Duration;

    ProcessEventTime(Event);
    return Event;
}

// ─────────────────────────────────────
/**
 * @brief Construct a REST event without advancing public score numbering.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @return Constructed rest event, or a default state for an initial rest or invalid syntax.
 *
 * @note Uses a silence observation and semi-Markov duration.
 * @warning Leading rests and their actions are skipped because they cannot be distinguished from pre-performance
 * silence.
 */
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
/**
 * @brief Assign section, expected onset, phase, and tempo to an event.
 *
 * @param Event Score event to populate or extend in place.
 *
 * @note Continues timing from the preceding state within a section and starts a new section at onset zero.
 * @warning Event.Index must identify the next insertion point or another valid predecessor position.
 */
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
/**
 * @brief Record a section declaration and defer creation of its boundary.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 *
 * @note Removes surrounding quotes and can attach the first section name to an existing leading boundary state.
 */
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
/**
 * @brief Create or refresh a pending section boundary before its first event.
 *
 * @param EventNode First event node of a pending section, used for source-line metadata.
 * @param Config Score configuration to read or update.
 *
 * @note Copies the current tempo, tolerance, and synchronization configuration to a FIRSTEVENT state.
 */
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
/**
 * @brief Parse, validate, and append a score event and its actions.
 *
 * @param ScoreStr Complete score source corresponding to the node byte offsets.
 * @param Node Tree-sitter node to inspect in the current score tree.
 * @param Config Score configuration to read or update.
 *
 * @note Dispatches by event type, resolves synchronization settings, and inserts eligible optional silence states.
 * @warning Updates parser state incrementally; failed events can leave numbering or boundary changes in place.
 */
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
        if (!Event.Observations.empty() || Event.MicroStates.empty()) {
            spdlog::error("Invalid microstate structure on line {}", Event.Line);
            return;
        }

        for (size_t k = 0; k < Event.MicroStates.size(); ++k) {
            AudioDescType Expected = PITCH;
            if (Event.Type == PTECH) {
                Expected = k == 0 ? LABEL : PITCH;
            } else if (Event.Type == UTECH) {
                Expected = LABEL;
            }

            const auto &Observations = Event.MicroStates[k].Observations;

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

} // namespace OpenScofo
