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

    void NewEvent(const std::string &Score, TSNode Event, Configuration &Config);
    void NewConfig(const std::string &Score, TSNode Node, Configuration &Config);
    void NewEventAction(const std::string &Score, TSNode Node, MarkovState &Event);

    void ProcessEventTime(MarkovState &Event);
    void ProcessNote(TSNode Note);

    // Get TreeSitter Values
    std::string GetCodeStr(const std::string &Score, TSNode Node);
    double GetDurationFromNode(const std::string &Score, TSNode Node);
    std::string GetChildStringFromField(const std::string &Score, TSNode node, std::string id);

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

    // Errors
    bool m_HasErrors = false;
    std::vector<std::string> m_Errors;
};
} // namespace OpenScofo
