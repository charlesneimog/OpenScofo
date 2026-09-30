#include <gtest/gtest.h>
#include <OpenScofo.hpp>

#if defined(OPENSCOFO_LUA)
namespace {

TEST(LuaTimers, UsesActualBlockSizeAndRunsSimultaneousTimersInOrder) {
    OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    ASSERT_TRUE(Scofo.LuaExecute(R"(
        openscofo = require('OpenScofo')
        seen = {}
        local function record(data) table.insert(seen, data.message) end
        later = openscofo.schedule(2, record, {message = 'later'})
        first = openscofo.schedule(1, record, {message = 'A'})
        second = openscofo.schedule(1, record, {message = 'B'})
        assert(later ~= first and first ~= second)
        collectgarbage()
    )")) << Scofo.LuaGetError();
    const float Audio[64] = {};
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 31));
    ASSERT_TRUE(Scofo.LuaExecute("assert(#seen == 0)")) << Scofo.LuaGetError();
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 16));
    ASSERT_TRUE(Scofo.LuaExecute("assert(#seen == 0)")) << Scofo.LuaGetError();
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 17));
    ASSERT_TRUE(Scofo.LuaExecute("assert(table.concat(seen) == 'AB')")) << Scofo.LuaGetError();
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 32));
    ASSERT_TRUE(Scofo.LuaExecute("assert(table.concat(seen) == 'ABlater'); assert(not openscofo.cancel(first))"))
        << Scofo.LuaGetError();
}

TEST(LuaTimers, RoundsFractionalSamplesUpAndUsesConfiguredSampleRate) {
    OpenScofo::OpenScofo Scofo(44100, 2048, 512);
    ASSERT_TRUE(Scofo.LuaExecute(R"(
        fired = false
        require('OpenScofo').schedule(1, function() fired = true end)
    )")) << Scofo.LuaGetError();
    const double Audio[44] = {};
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 44));
    ASSERT_TRUE(Scofo.LuaExecute("assert(not fired)")) << Scofo.LuaGetError();
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 1));
    ASSERT_TRUE(Scofo.LuaExecute("assert(fired)")) << Scofo.LuaGetError();
}

TEST(LuaTimers, AcceptsAnyDataAndDefersZeroDelay) {
    OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    ASSERT_TRUE(Scofo.LuaExecute(R"(
        local openscofo = require('OpenScofo')
        count = 0
        local values = {42, 'hello', true, false, {}, function() end, coroutine.create(function() end)}
        for _, value in ipairs(values) do
            openscofo.schedule(0, function(data)
                assert(data == value)
                count = count + 1
            end, value)
        end
        local function missing(...)
            assert(select('#', ...) == 1 and (...) == nil)
            count = count + 1
        end
        openscofo.schedule(0, missing)
        openscofo.schedule(0, missing, nil)
        assert(count == 0)
    )")) << Scofo.LuaGetError();
    const float Audio[64] = {};
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 64));
    ASSERT_TRUE(Scofo.LuaExecute("assert(count == 9)")) << Scofo.LuaGetError();
}

TEST(LuaTimers, CancelsPendingTimersIncludingFromCallbacks) {
    OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    ASSERT_TRUE(Scofo.LuaExecute(R"(
        local openscofo = require('OpenScofo')
        fired = false
        local function callback() fired = true end
        local id = openscofo.schedule(5000, callback, 'cancel test')
        assert(openscofo.cancel(id))
        assert(not openscofo.cancel(id))
        assert(not openscofo.cancel(9999))
        assert(not openscofo.cancel(-1))
        local pending
        openscofo.schedule(0, function() assert(openscofo.cancel(pending)) end)
        pending = openscofo.schedule(0, callback)
    )")) << Scofo.LuaGetError();
    const float Audio[64] = {};
    for (int Block = 0; Block < 3751; ++Block) {
        ASSERT_TRUE(Scofo.ProcessBlock(Audio, 64));
    }
    ASSERT_TRUE(Scofo.LuaExecute("assert(not fired)")) << Scofo.LuaGetError();
}

TEST(LuaTimers, SchedulesFromCallbacksWithoutExecutingInsideSchedule) {
    OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    ASSERT_TRUE(Scofo.LuaExecute(R"(
        local openscofo = require('OpenScofo')
        seen = {}
        openscofo.schedule(500, function(data)
            table.insert(seen, data)
            openscofo.schedule(500, function(inner) table.insert(seen, inner) end, 'inner')
            local scheduling = true
            openscofo.schedule(0, function() assert(not scheduling); table.insert(seen, 'zero') end)
            scheduling = false
        end, 'outer')
    )")) << Scofo.LuaGetError();
    const float Audio[64] = {};
    for (int Block = 0; Block < 375; ++Block) {
        ASSERT_TRUE(Scofo.ProcessBlock(Audio, 64));
    }
    ASSERT_TRUE(Scofo.LuaExecute("assert(table.concat(seen, ',') == 'outer,zero')")) << Scofo.LuaGetError();
    for (int Block = 0; Block < 374; ++Block) {
        ASSERT_TRUE(Scofo.ProcessBlock(Audio, 64));
    }
    ASSERT_TRUE(Scofo.LuaExecute("assert(#seen == 2)")) << Scofo.LuaGetError();
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 64));
    ASSERT_TRUE(Scofo.LuaExecute("assert(table.concat(seen, ',') == 'outer,zero,inner')")) << Scofo.LuaGetError();
}

TEST(LuaTimers, LogsErrorsContinuesAndReleasesReferences) {
    std::string Errors;
    OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    Scofo.SetErrorCallback([&](const spdlog::details::log_msg &Log, void *) {
        if (Log.level == spdlog::level::err) {
            Errors.append(Log.payload.data(), Log.payload.size());
        }
    });
    lua_State *State = nullptr;
    Scofo.LuaAddPointer(&State, "test_state");
    Scofo.LuaAddModule("capture_state", [](lua_State *L) {
        lua_getglobal(L, "test_state");
        *static_cast<lua_State **>(lua_touserdata(L, -1)) = L;
        lua_pop(L, 1);
        lua_pushboolean(L, true);
        return 1;
    });
    ASSERT_NE(State, nullptr);
    ASSERT_TRUE(Scofo.LuaExecute(R"(
        local openscofo = require('OpenScofo')
        weak = setmetatable({}, {__mode = 'v'})
        local callbacks = {
            function() error('expected timer failure') end,
            function() error({}) end,
            function(data) assert(data.message == 'hello'); succeeded = true end,
            function() error('cancelled callback executed') end
        }
        for i, callback in ipairs(callbacks) do
            local data = {message = 'hello'}
            weak[2*i - 1], weak[2*i] = callback, data
            local id = openscofo.schedule(0, callback, data)
            if i == 4 then assert(openscofo.cancel(id)) end
        end
        callbacks = nil
        collectgarbage()
        assert(weak[1] and weak[2] and weak[5] and weak[6])
        assert(weak[7] == nil and weak[8] == nil)
    )")) << Scofo.LuaGetError();
    const int StackTop = lua_gettop(State);
    const float Audio[64] = {};
    ASSERT_TRUE(Scofo.ProcessBlock(Audio, 64));
    EXPECT_EQ(lua_gettop(State), StackTop);
    EXPECT_NE(Errors.find("expected timer failure"), std::string::npos);
    EXPECT_NE(Errors.find("Unknown error"), std::string::npos);
    ASSERT_TRUE(Scofo.LuaExecute("assert(succeeded); collectgarbage(); assert(next(weak) == nil)"))
        << Scofo.LuaGetError();
    Scofo.SetErrorCallback(nullptr);
}

TEST(LuaTimers, RejectsInvalidArguments) {
    OpenScofo::OpenScofo Scofo(48000, 2048, 512);
    ASSERT_TRUE(Scofo.LuaExecute(R"(
        local openscofo = require('OpenScofo')
        local function callback() end
        for _, delay in ipairs({-1, '500', false, {}, math.huge, -math.huge, 0/0, 1e300}) do
            assert(not pcall(openscofo.schedule, delay, callback))
        end
        assert(not pcall(openscofo.schedule, nil, callback))
        assert(not pcall(openscofo.schedule, 1, {}))
        assert(not pcall(openscofo.schedule, 1))
        assert(not pcall(openscofo.cancel, {}))
        assert(not pcall(openscofo.cancel, 1.5))
    )")) << Scofo.LuaGetError();
}

TEST(LuaTimers, CleansUpOnReinitializationAndDestruction) {
    int Finalized = 0;
    {
        OpenScofo::OpenScofo Scofo(48000, 2048, 512);
        const auto SchedulePending = [&]() {
            Scofo.LuaAddPointer(&Finalized, "finalized");
            Scofo.LuaAddModule("finalize", [](lua_State *L) {
                lua_pushcfunction(L, [](lua_State *State) {
                    lua_getglobal(State, "finalized");
                    ++*static_cast<int *>(lua_touserdata(State, -1));
                    return 0;
                });
                return 1;
            });
            return Scofo.LuaExecute(R"(
                require('OpenScofo').schedule(5000, function() error('stale timer') end,
                    setmetatable({}, {__gc = require('finalize')}))
            )");
        };
        ASSERT_TRUE(SchedulePending()) << Scofo.LuaGetError();
        Scofo.InitLuaModule();
        EXPECT_EQ(Finalized, 1);
        ASSERT_TRUE(Scofo.LuaExecute("assert(not require('OpenScofo').cancel(1))")) << Scofo.LuaGetError();
        const float Audio[64] = {};
        ASSERT_TRUE(Scofo.ProcessBlock(Audio, 64));
        ASSERT_TRUE(SchedulePending()) << Scofo.LuaGetError();
    }
    EXPECT_EQ(Finalized, 2);
}

} // namespace
#endif
