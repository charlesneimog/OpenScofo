#pragma once

#include <string>
#include <vector>
#include <filesystem>

#include "states.hpp"
#include <tree_sitter/api.h>

namespace OpenScofo {

#ifndef TWO_PI
#define TWO_PI (2 * M_PI)
#endif

namespace fs = std::filesystem;

// ╭─────────────────────────────────────╮
// │                Score                │
// ╰─────────────────────────────────────╯

class Score {
  public:
    Score(float FftSize, float HopSize) {
        m_FFTSize = FftSize;
        m_HopSize = HopSize;
    }

    States Parse(std::string ScoreFile);
    void SetTunning(double Tunning);
    bool ScoreIsLoaded();

    // Get Score
    std::string GetLuaCode();
    double GetFFTSize();
    double GetHopSize();
    double GetPitchTemplateSigma();
    double GetPhaseCoupling();
    double GetSyncStrength();

    // AI
    bool HasTimbreModel() {
        if (m_TimbreModel.empty()) {
            return false;
        } else {
            return true;
        }
    }
    fs::path GetTimbreModel() {
        return m_TimbreModel;
    }

  private:
    States m_ScoreStates;
    std::string m_LuaCode;

    // Helpers
    MarkovState AddDummySilence();
    double ModPhases(double Phase);
    MarkovState AddTransState(MarkovState &State, int ScoreEvent, int BPM);
    double PitchName2Midi(char pitchName, std::string alt, std::string octave);
    void PitchNode2Freq(const std::string Score, TSNode node, AudioState &State);
    // bool SpaceTab(const std::string &line, int numSpaces);
    void PrintTreeSitterNode(TSNode node, int indent = 0);
    TSNode GetField(TSNode Node, std::string s);
    bool isNumber(std::string str);

    void ProcessEvent(const std::string &Score, TSNode Event);
    void ProcessEventTime(MarkovState &Event);
    void ProcessConfig(const std::string &Score, TSNode Node);
    void ProcessAction(const std::string &Score, TSNode Node, MarkovState &Event);

    void ProcessNote(TSNode Note);

    // Get TreeSitter Values
    std::string GetCodeStr(const std::string &Score, TSNode Node);
    // double GetFreqsFromNode(const std::string &Score, TSNode Node);
    double GetDurationFromNode(const std::string &Score, TSNode Node);
    std::string GetChildStringFromField(const std::string &Score, TSNode node, std::string id);

    // Events
    MarkovState NewRestEvent(const std::string &Score, TSNode Node);
    MarkovState NewPitchEvent(const std::string &Score, TSNode Node);
    MarkovState TrillEvent(const std::string &Score, TSNode Node);

    // Add events
    MarkovState GetFirstEvent();
    MarkovState AddNote(std::vector<std::string> Tokens);
    MarkovState AddChord(std::vector<std::string> Tokens);
    MarkovState AddTrill(std::vector<std::string> Tokens);
    MarkovState AddMulti(std::vector<std::string> Tokens);
    MarkovState AddDumpSilence();
    MarkovState AddRest(std::vector<std::string> Tokens);
    void AddAction(std::vector<std::string> Tokens);

    // Some ScoreConfigs
    double m_CurrentBPM = 60;
    float m_Transpose = 0;
    double m_Entropy = 0;
    double m_PitchTemplateSigma = 0.5;
    double m_SyncStrength = 0.5;
    double m_PhaseCoupling = 0.5;
    double m_FFTSize;
    double m_HopSize;

    // Paths
    fs::path m_ScoreRootPath;
    fs::path m_TimbreModel;

    // Variables
    int m_ScorePosition = 1;
    int m_LineCount = 0;
    int m_MarkovIndex = 0;
    double m_LastOnset = 0;
    double m_LastPhase = 0;
    double m_PrevDuration;
    double m_Tunning = 440;
    bool m_ScoreLoaded = false;
    double m_MinimalDuration = std::numeric_limits<double>::infinity();

    // Errors
    bool m_HasErrors = false;
    std::vector<std::string> m_Errors;
};
} // namespace OpenScofo
