#include <OpenScofo.hpp>

// ╭─────────────────────────────────────╮
// │     Construstor and Destructor      │
// ╰─────────────────────────────────────╯
namespace OpenScofo {

#if defined(OSCOFO_LUA)
int luaopen_oscofo(lua_State *L);
#endif

//  ─────────────────────────────────────
OpenScofo::OpenScofo(float Sr, float FftSize, float HopSize) : m_MDP(Sr, FftSize, HopSize), m_MIR(Sr, FftSize, HopSize), m_Score(FftSize, HopSize) {
    m_States = States();
    m_Desc = Description();
    m_Sr = Sr;
    m_FFTSize = FftSize;
    m_HopSize = HopSize;
    m_InBuffer.reserve(FftSize);
    m_BlockIndex = 0;

    if (m_MIR.HasErrors() || m_MDP.HasErrors()) {
        for (auto &error : m_MIR.GetErrorMessage()) {
            SetError(error);
        }
        m_MIR.ClearError();
        for (auto &error : m_MDP.GetErrorMessage()) {
            SetError(error);
        }
        m_MDP.ClearError();
        return;
    }

#if defined(OSCOFO_LUA)
    InitLua();
#endif
}

//  ─────────────────────────────────────
void OpenScofo::SetNewAudioParameters(float Sr, float FftSize, float HopSize) {
    if (m_FFTSize == FftSize && m_HopSize == HopSize && m_Sr == Sr) {
        return;
    }
    m_Sr = Sr;
    m_FFTSize = FftSize;
    m_HopSize = HopSize;
    m_MDP = MDP(Sr, FftSize, HopSize);
    m_MIR = MIR(Sr, FftSize, HopSize);
}

// ╭─────────────────────────────────────╮
// │               Errors                │
// ╰─────────────────────────────────────╯
bool OpenScofo::HasErrors() {
    return m_HasErrors;
}

// ─────────────────────────────────────
std::vector<std::string> OpenScofo::GetErrorMessage() {
    return m_Errors;
}

// ─────────────────────────────────────
void OpenScofo::SetError(const std::string &message) {
    if (m_ErrorCallback) {
        m_HasErrors = true;
        m_ErrorCallback(message);
    } else {
        m_HasErrors = true;
        m_Errors.push_back(message);
    }
}

// ─────────────────────────────────────
void OpenScofo::ClearError() {
    m_HasErrors = false;
    m_Errors.clear();
}

// ╭─────────────────────────────────────╮
// │                 Lua                 │
// ╰─────────────────────────────────────╯
#if defined(OSCOFO_LUA)
void OpenScofo::InitLua() {
    m_LuaState = luaL_newstate();
    luaL_openlibs(m_LuaState); // NOTE: Rethink if I load all functions
    lua_newtable(m_LuaState);
    lua_pushlightuserdata(m_LuaState, this);
    lua_setfield(m_LuaState, -2, "pointer");
    lua_setglobal(m_LuaState, "_OpenScofo");
    luaL_requiref(m_LuaState, "oscofo", luaopen_oscofo, 1);
}

// ─────────────────────────────────────
bool OpenScofo::LuaAddModule(std::string name, lua_CFunction func) {
    if (m_LuaState == nullptr) {
        return false;
    }
    luaL_requiref(m_LuaState, name.c_str(), func, 1);
    if (lua_isnil(m_LuaState, -1)) {
        return false;
    }
    return true;
}

// ─────────────────────────────────────
bool OpenScofo::LuaExecute(std::string code) {
    if (m_LuaState == nullptr) {
        return false;
    }
    int status = luaL_loadstring(m_LuaState, code.c_str());
    if (status == LUA_OK) {
        status = lua_pcall(m_LuaState, 0, LUA_MULTRET, 0);
        if (status != LUA_OK) {
            return false;
        }
        return true;
    } else {
        return false;
    }
}

// ─────────────────────────────────────
bool OpenScofo::LuaAddPointer(void *pointer, const char *name) {
    if (m_LuaState == nullptr) {
        return false;
    }
    lua_pushlightuserdata(m_LuaState, pointer);
    lua_setglobal(m_LuaState, name);
    return true;
}

// ─────────────────────────────────────
void OpenScofo::LuaAddPath(std::string path) {
    if (m_LuaState == nullptr) {
        return;
    }

    lua_getglobal(m_LuaState, "package");
    lua_getfield(m_LuaState, -1, "path");
    const char *current_path = lua_tostring(m_LuaState, -1);
    if (path.back() != '/') {
        lua_pushfstring(m_LuaState, "%s;%s/?.lua", current_path, path.c_str());
    } else {
        lua_pushfstring(m_LuaState, "%s;%s?.lua", current_path, path.c_str());
    }

    lua_setfield(m_LuaState, -3, "path");
    lua_pop(m_LuaState, 1);
}

// ─────────────────────────────────────
std::string OpenScofo::LuaGetError() {
    if (m_LuaState == nullptr) {
        return "m_LuaState is null";
    }
    if (lua_isstring(m_LuaState, -1)) {
        std::string errorMsg = lua_tostring(m_LuaState, -1);
        lua_pop(m_LuaState, 1);
        return errorMsg;
    }
    return "Unknown error";
}
#endif

// ╭─────────────────────────────────────╮
// │            Set Functions            │
// ╰─────────────────────────────────────╯
void OpenScofo::SetPitchTemplateSigma(double Sigma) {
    m_MDP.SetPitchTemplateSigma(Sigma);
    m_MDP.UpdateAudioTemplate();
}

// ─────────────────────────────────────
void OpenScofo::SetAmplitudeDecay(double decay) {
    m_MDP.SetAmplitudeDecay(decay);
}

// ─────────────────────────────────────
void OpenScofo::SetMinEntropy(double EntropyValue) {
    m_MDP.SetMinEntropy(EntropyValue);
}

// ─────────────────────────────────────
void OpenScofo::SetHarmonics(int Harmonics) {
    m_MDP.SetHarmonics(Harmonics);
    m_MDP.UpdateAudioTemplate();
}

// ─────────────────────────────────────
double OpenScofo::GetdBValue() {
    return 0;
}

// ─────────────────────────────────────
void OpenScofo::SetdBTreshold(double dB) {
    m_MDP.SetdBTreshold(dB);
    m_MIR.SetdBTreshold(dB);
}

// ─────────────────────────────────────
void OpenScofo::SetTunning(double Tunning) {
    m_Score.SetTunning(Tunning);
}

// ─────────────────────────────────────
void OpenScofo::SetCurrentEvent(int Event) {
    m_MDP.SetCurrentEvent(Event);
}

// ╭─────────────────────────────────────╮
// │            Get Functions            │
// ╰─────────────────────────────────────╯
int OpenScofo::GetEventIndex() {
    return m_CurrentScorePosition; // TODO: Implement yet
}

// ─────────────────────────────────────
double OpenScofo::GetLiveBPM() {
    return m_MDP.GetLiveBPM();
}

// ─────────────────────────────────────
ActionVec OpenScofo::GetEventActions(int Index) {
    return m_MDP.GetEventActions(Index);
}

// ─────────────────────────────────────
double OpenScofo::GetKappa() {
    return m_MDP.GetKappa();
}

// ─────────────────────────────────────
double OpenScofo::GetPitchProb(double f) {
    return m_MDP.GetPitchSimilarity(f);
}

// ─────────────────────────────────────
std::string OpenScofo::GetLuaCode() {
    return m_Score.GetLuaCode();
}

// ╭─────────────────────────────────────╮
// │          Helpers Functions          │
// ╰─────────────────────────────────────╯
bool OpenScofo::ScoreIsLoaded() {
    return m_Score.ScoreIsLoaded();
}

// ╭─────────────────────────────────────╮
// │ Python Research and Test Functions  │
// ╰─────────────────────────────────────╯
States OpenScofo::GetStates() {
    if (m_States.size() != 0) {
        return m_States;
    }
    SetError("No states found, please use the ScoreParse first");
    return m_States;
}

// ─────────────────────────────────────
std::vector<double> OpenScofo::GetPitchTemplate(double Freq) {
    return m_MDP.GetPitchTemplate(Freq);
}

// ─────────────────────────────────────
std::vector<double> OpenScofo::GetCQTTemplate(double Freq) {
    PitchTemplateArray p = m_MDP.GetPitchTemplate(Freq);
    std::vector<std::pair<int, int>> cqt = m_MIR.GetCQT();

    std::vector<double> cqt_template;
    cqt_template.resize(cqt.size());

    for (size_t k = 0; k < cqt.size(); ++k) {
        auto [b0, b1] = cqt[k];
        double sum = 0.0;
        for (int i = b0; i <= b1; i++) {
            sum += p[i];
        }
        cqt_template[k] = sum / double(b1 - b0 + 1);
    }

    return cqt_template;
}

// ─────────────────────────────────────
std::vector<double> OpenScofo::GetSpectrumPower() const {
    return m_Desc.NormSpectralPower;
}

// ─────────────────────────────────────
double OpenScofo::GetSr() {
    return m_Sr;
}

// ─────────────────────────────────────
double OpenScofo::GetFFTSize() {
    return m_FFTSize;
}

// ─────────────────────────────────────
double OpenScofo::GetHopSize() {
    return m_HopSize;
}

// ─────────────────────────────────────
std::vector<float> OpenScofo::GetTimeCoherenceTemplate(int pos, int timeInEvent) {
    return m_MIR.GetTimeCoherenceTemplate(m_States, pos, timeInEvent);
}

// ─────────────────────────────────────
double OpenScofo::GetTimeCoherenceConfiability(const std::vector<double> &eventValues) {
    return m_MIR.GetTimeCoherenceConfiability(eventValues);
}

// ╭─────────────────────────────────────╮
// │           Main Functions            │
// ╰─────────────────────────────────────╯
bool OpenScofo::ParseScore(std::string ScorePath) {
    m_Score = Score(m_FFTSize, m_HopSize);
    m_States.clear();
    m_States = m_Score.Parse(ScorePath);
    if (m_Score.HasErrors()) {
        for (auto &error : m_Score.GetErrorMessage()) {
            SetError(error);
            m_Score.ClearError();
        }
        return false;
    }

    // Time coherence
    m_MIR.BuildTimeCoherenceTemplate(m_States);

    // Timbre detection
    if (m_Score.HasTimbreModel()) {
        m_MIR.LoadONNXModel(m_Score.GetTimbreModel());
    }

    m_FFTSize = m_Score.GetFFTSize();
    m_HopSize = m_Score.GetHopSize();
    SetNewAudioParameters(m_Sr, m_FFTSize, m_HopSize);

    // Parse Config
    m_MDP.SetPitchTemplateSigma(m_Score.GetPitchTemplateSigma());

    // Add States
    m_MDP.SetScoreStates(m_States);
    return true;
}

// ─────────────────────────────────────
Description OpenScofo::GetDescription() {
    return m_Desc;
}

// ─────────────────────────────────────
Description OpenScofo::GetAudioDescription(std::vector<double> &AudioBuffer) {
    if (m_FFTSize != AudioBuffer.size()) {
        SetError(std::format("AudioBuffer size ({}) differs from FFT size ({})", AudioBuffer.size(), m_FFTSize));
        return {};
    }

    SetNewAudioParameters(m_Sr, m_FFTSize, m_HopSize);
    States GoodStates = m_MDP.GetStatesForProcessing();
    m_MIR.GetDescription(AudioBuffer, m_Desc, GoodStates);
    return m_Desc;
}

// ─────────────────────────────────────
bool OpenScofo::ProcessBlock(std::vector<double> &AudioBuffer) {
    if (!m_Score.ScoreIsLoaded()) {
        return false;
    }

    States GoodStates = m_MDP.GetStatesForProcessing();
    m_MIR.GetDescription(AudioBuffer, m_Desc, GoodStates);
    m_CurrentScorePosition = m_MDP.GetEvent(m_Desc);
    // m_MIR.AddReverb(m_Desc, 0.001);

    if (m_MDP.HasErrors()) {
        for (auto &error : m_MDP.GetErrorMessage()) {
            SetError(error);
        }
        m_MDP.ClearError();
        return false;
    }

    return true;
}

// bool OpenScofo::ProcessBlock(std::vector<double> &AudioBuffer) {
//     // move accumation to here buffer
//     size_t n = AudioBuffer.size();
//     memcpy(m_InBuffer.data() + m_BlockIndex, AudioBuffer.data(), n * sizeof(double));
//     m_BlockIndex += n;
//
//     if (m_BlockIndex > m_FFTSize) {
//         SetError("The configuration of OpenScofo is wrong, please review");
//         return false;
//     }
//
//     if (m_FFTSize != m_BlockIndex) {
//         return true;
//     }
//
//     States GoodStates = m_MDP.GetStatesForProcessing();
//     m_MIR.GetDescription(m_InBuffer, m_Desc, GoodStates);
//     m_CurrentScorePosition = m_MDP.GetEvent(m_Desc);
//
//     if (m_MDP.HasErrors()) {
//         for (auto &error : m_MDP.GetErrorMessage()) {
//             SetError(error);
//         }
//         m_MDP.ClearError();
//         return false;
//     }
//
//     m_BlockIndex = 0;
//     return true;
// }

} // namespace OpenScofo
