/*
    Copyright (c) 2024-2026 Charles K. Neimog
    Website: charlesneimog.github.io

    This file is part of a project licensed under the
    GNU General Public License v3.0 or later (GPL-3.0-or-later).
    See the LICENSE file for details.
*/

#pragma once

#include <string>
#include <vector>
#include <filesystem>

#include "states.hpp"
#include <tree_sitter/api.h>

namespace OpenScofo {

namespace fs = std::filesystem;

// ╭─────────────────────────────────────╮
// │                Score                │
// ╰─────────────────────────────────────╯

class Score {
  public:
    Score() = default;

    std::pair<Configuration, States> Parse(fs::path ScoreFile);
    bool ScoreIsLoaded();
    std::string GetLuaCode();

  private:
    Configuration *m_Config;

    // Helpers
    MarkovState AddDummySilence();
    double ModPhases(double Phase);
    MarkovState AddTransState(MarkovState &State, int ScoreEvent, int BPM);
    double PitchName2Midi(char pitchName, std::string alt, std::string octave);
    void PitchNode2Freq(const std::string Score, TSNode node, AudioState &State);
    void ParseInput(const std::string &Score);
    void PrintTreeSitterNode(TSNode node, int indent = 0);
    TSNode GetField(TSNode Node, std::string s);
    bool isNumber(std::string str);
    static bool ScoreIsText(const std::string &Path);
    static void FindErrors(TSNode &Root, TSNode &Node, const std::string &Score);
    static bool GetConfigNumber(const std::string &Id, const std::string &ValueType, const std::string &Value,
                                TSPoint Position, double &Output);
    static bool GetConfigBool(const std::string &Id, const std::string &ValueType, std::string Value, TSPoint Position,
                              bool &Output);

    void NewEvent(const std::string &Score, TSNode Event, Configuration &Config);
    void NewConfig(const std::string &Score, TSNode Node, Configuration &Config);
    void NewSection(const std::string &Score, TSNode Node);
    void EnsureSectionStart(TSNode Event, Configuration &Config);
    void NewEventAction(const std::string &Score, TSNode Node, MarkovState &Event);

    void ProcessEventTime(MarkovState &Event);
    void ProcessNote(TSNode Note);

    // Get TreeSitter Values
    std::string GetCodeStr(const std::string &Score, TSNode Node);
    double GetDurationFromNode(const std::string &Score, TSNode Node);
    std::string GetChildStringFromField(const std::string &Score, TSNode node, std::string id);
    static std::string GetChildStringFromType(const std::string &Score, TSNode Parent, const std::string &WantedType);

    // Events
    MarkovState NewRestEvent(const std::string &Score, TSNode Node);
    MarkovState NewPitchEvent(const std::string &Score, TSNode Node);
    MarkovState NewMultiPitchEvent(const std::string &Score, TSNode Node);
    MarkovState NewPTechEvent(const std::string &ScoreStr, TSNode Node);
    MarkovState NewUTechEvent(const std::string &ScoreStr, TSNode Node);

    // Add events
    MarkovState GetFirstEvent();
    MarkovState AddDumpSilence();
    void AddAction(std::vector<std::string> Tokens);

  private:
    States m_ScoreStates;
    std::string m_LuaCode;
    double m_CurrentBPM = 60;
    double m_Transpose = 0;
    std::string m_CurrentSection;
    bool m_HasSection = false;
    bool m_SectionStartPending = false;
    int m_ScorePosition;

    // Paths
    fs::path m_ScoreRootPath;

    // Variables
    int m_LineCount = 0;
    int m_MarkovIndex = 0;
    double m_LastOnset = 0;
    double m_LastPhase = 0;
    double m_PrevDuration;
    double m_Tunning = 440;
    bool m_ScoreLoaded = false;

    // State configuration
    double m_TimeTolerance = 16;

    // Errors
    bool m_HasErrors = false;
    std::vector<std::string> m_Errors;
};
} // namespace OpenScofo
