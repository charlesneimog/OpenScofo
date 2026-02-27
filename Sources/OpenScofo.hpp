#pragma once

#include <OpenScofo/mdp.hpp>
#include <OpenScofo/mir.hpp>
#include <OpenScofo/score.hpp>
#include <OpenScofo/states.hpp>
#include <OpenScofo/log.hpp>

#if defined(OSCOFO_LUA)
extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}
#endif

#define OSCOFO_VERSION_MAJOR 0
#define OSCOFO_VERSION_MINOR 1
#define OSCOFO_VERSION_PATCH 4

// vector
#include <vector>
#include <functional>

// log
#include <chrono>
#include <iostream>

#include <cmath>
#include <numeric>
#include <spdlog/spdlog.h>

class Timer {
  public:
    Timer(const std::string &name = "") : m_Name(name), m_Start(std::chrono::high_resolution_clock::now()) {
    }

    ~Timer() {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - m_Start).count();
        if (!m_Name.empty())
            std::cout << m_Name << " ";
        std::cout << "Time elapsed: " << duration << " µs" << std::endl;
    }

  private:
    std::string m_Name;
    std::chrono::high_resolution_clock::time_point m_Start;
};

namespace OpenScofo {

class OpenScofo;

class OpenScofo {
  public:
    OpenScofo(float Sr, float WindowSize, float HopSize);
    // Main Functions
    bool ParseScore(std::string ScorePath);
    bool ProcessBlock(std::vector<double> &AudioBuffer);

    bool ScoreIsLoaded();

    // Set Functions
    void SetAmplitudeDecay(double decay);
    void SetPitchTemplateSigma(double Sigma);
    void SetHarmonics(int Harmonics);
    void SetdBTreshold(double dB);
    void SetTunning(double Tunning);
    void SetCurrentEvent(int Event);
    void SetMinEntropy(double EntropyValue);
    void SetNewAudioParameters(float Sr, float FftSize, float HopSize);

    // Get Functions
    double GetLiveBPM();
    int GetEventIndex();
    double GetKappa();
    double GetdBValue();
    ActionVec GetEventActions(int Index);
    std::string GetLuaCode();
    double GetPitchProb(double f);
    States GetStates();
    PitchTemplateArray GetPitchTemplate(double Freq);
    std::vector<double> GetSpectrumPower() const;
    double GetSr();
    double GetFFTSize();
    double GetHopSize();
    double GetBlockDuration();
    Description GetAudioDescription(std::vector<double> &AudioBuffer);
    Description GetDescription();
    std::vector<double> GetCQTTemplate(double Freq);
    std::vector<float> GetTimeCoherenceTemplate(int pos, int timeInEvent = 0);
    double GetTimeCoherenceConfiability(const std::vector<double> &eventValues);

    // Config
    void ClearErrors();

#if defined(OSCOFO_LUA)
    void InitLua();
    bool LuaExecute(std::string code);
    std::string LuaGetError();
    bool LuaAddModule(std::string name, lua_CFunction func);
    bool LuaAddPointer(void *pointer, const char *name);
    void LuaAddPath(std::string path);
#endif

    // Errors
    void SetErrorCallback(std::function<void(const spdlog::details::log_msg &, void *data)> cb, void *data = nullptr);
    void SetLogLevel(spdlog::level::level_enum level);

  private:
    MDP m_MDP;
    MIR m_MIR;
    Score m_Score;
    std::shared_ptr<OpenScofoLog<std::mutex>> m_Log;

#if defined(OSCOFO_LUA)
    lua_State *m_LuaState;
#endif

    States m_States;
    Description m_Desc;
    int m_CurrentScorePosition = -1;

    double m_Sr;
    double m_FFTSize;
    double m_HopSize;
    unsigned m_BlockIndex;

    // Errors
    spdlog::level::level_enum m_HasErrors;
    std::vector<std::string> m_Errors;
    std::function<void(const std::string &)> m_ErrorCallback = nullptr;
    std::vector<double> m_InBuffer;
};

} // namespace OpenScofo
