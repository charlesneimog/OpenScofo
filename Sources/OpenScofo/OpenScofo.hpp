#pragma once

#include "mdp.hpp"
#include "mir.hpp"
#include "score.hpp"
#include "states.hpp"
#include "log.hpp"

#if defined(OSCOFO_LUA)
extern "C" {
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}
#endif

#define OSCOFO_BUILD_TIME (__DATE__ " " __TIME__)

// vector
#include <vector>
#include <functional>

// log
#include <chrono>
#include <iostream>

#include <concepts>
#include <cmath>
#include <numbers>
#include <numeric>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

namespace OpenScofo {

template <typename T>
concept OpenScofoPrecision = std::same_as<T, float> || std::same_as<T, double>;

class OpenScofo {
  public:
    OpenScofo(float Sr, float WindowSize, float HopSize);
    // Main Functions
    bool ParseScore(fs::path ScorePath);
    bool ScoreIsLoaded();

    template <OpenScofoPrecision T> bool ProcessBlock(const T *AudioBuffer, size_t n);

    // ONNX
    void LoadONNXModel(fs::path Model, std::vector<Descriptors> Descriptors);

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
    EventActions GetEventActions(int Index);
    std::string GetLuaCode();
    double GetPitchProb(double f);
    States GetStates();
    PitchTemplateArray GetPitchTemplate(double Freq);
    std::vector<double> GetSpectrumPower() const;
    double GetSr();
    double GetFFTSize();
    double GetHopSize();
    double GetBlockDuration();
    Description GetDescription();

    Descriptors GetDescriptorsEnum(const char *s);
    const char *GetDescriptionId(Descriptors d);
    double GetDescriptionFloat(Description &Desc, Descriptors d);
    std::vector<double> &GetDescriptionArray(Description &Desc, Descriptors d);

    // Config
    void ClearErrors();

#if defined(OSCOFO_LUA)
    void InitLuaModule();
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

    Mode m_Mode = DESCRIPTORS;
    States m_States;
    Description m_Desc;
    int m_CurrentScorePosition = -1;

    int m_Sr;
    int m_FFTSize;
    int m_HopSize;
    int m_BlockIndex = 0;
    std::vector<double> m_InputBuffer;

    // AI
    bool m_LoadingONNX = false;

    // Errors
    spdlog::level::level_enum m_HasErrors;
    std::vector<std::string> m_Errors;
    std::function<void(const std::string &)> m_ErrorCallback = nullptr;
    std::vector<double> m_InBuffer;
};

} // namespace OpenScofo
