#ifdef OSCOFO_LUA

#include <OpenScofo.hpp>
#include <ext.h>
#include <z_dsp.h>

// ─────────────────────────────────────
static int max_Post(lua_State *L) {
    int num_args = lua_gettop(L);
    for (int i = 1; i <= num_args; i++) {
        if (lua_type(L, i) == LUA_TSTRING) {
            object_post(nullptr, "%s", lua_tostring(L, i));
        } else if (lua_type(L, i) == LUA_TNUMBER) {
            object_post(nullptr, "%f", lua_tonumber(L, i));
        } else if (lua_type(L, i) == LUA_TBOOLEAN) {
            if (lua_toboolean(L, i)) {
                object_post(nullptr, "true");
            } else {
                object_post(nullptr, "false");
            }
        } else {
            object_post(nullptr, "Unsupported type: %s", lua_typename(L, lua_type(L, i)));
        }
    }

    return 0;
}

// ─────────────────────────────────────
static int max_Error(lua_State *L) {
    int num_args = lua_gettop(L);
    for (int i = 1; i <= num_args; i++) {
        if (lua_type(L, i) == LUA_TSTRING) {
            object_error(nullptr, "%s", lua_tostring(L, i));
        } else if (lua_type(L, i) == LUA_TNUMBER) {
            object_error(nullptr, "%f", lua_tonumber(L, i));
        } else if (lua_type(L, i) == LUA_TBOOLEAN) {
            if (lua_toboolean(L, i)) {
                object_error(nullptr, "true");
            } else {
                object_error(nullptr, "false");
            }
        } else {
            object_error(nullptr, "Unsupported type: %s", lua_typename(L, lua_type(L, i)));
        }
    }
    return 0;
}

// ─────────────────────────────────────
static int max_sendBang(lua_State *L) {
    const char *r = luaL_checkstring(L, 1);
    t_symbol *symbol = gensym(r);
    if (symbol->s_thing) {
        object_method_typed(symbol->s_thing, gensym("bang"), 0, nullptr, nullptr);
    } else {
        luaL_error(L, "Receiver not found");
    }
    return 0;
}

// ─────────────────────────────────────
static int max_sendFloat(lua_State *L) {
    const char *r = luaL_checkstring(L, 1);
    const double f = luaL_checknumber(L, 2);
    t_symbol *symbol = gensym(r);
    if (symbol->s_thing) {
        t_atom atom;
        atom_setfloat(&atom, f);
        object_method_typed(symbol->s_thing, gensym("float"), 1, &atom, nullptr);
    } else {
        luaL_error(L, "Receiver not found");
    }
    return 0;
}

// ─────────────────────────────────────
static int max_sendSymbol(lua_State *L) {
    const char *r = luaL_checkstring(L, 1);
    const char *s = luaL_checkstring(L, 2);
    t_symbol *symbol = gensym(r);
    if (symbol->s_thing) {
        t_atom atom;
        atom_setsym(&atom, gensym(s));
        object_method_typed(symbol->s_thing, gensym("symbol"), 1, &atom, nullptr);
    } else {
        luaL_error(L, "Receiver not found");
    }
    return 0;
}

// ─────────────────────────────────────
static int max_sendList(lua_State *L) {
    const char *r = luaL_checkstring(L, 1);
    luaL_checktype(L, 2, LUA_TTABLE);
    t_symbol *symbol = gensym(r);
    if (!symbol->s_thing) {
        return luaL_error(L, "Receiver not found");
    }

    const int listSize = luaL_len(L, 2);
    std::vector<t_atom> atomList(listSize);

    for (int i = 0; i < listSize; i++) {
        lua_rawgeti(L, 2, i + 1);
        if (lua_isnumber(L, -1)) {
            atom_setfloat(&atomList[i], lua_tonumber(L, -1));
        } else if (lua_isstring(L, -1)) {
            atom_setsym(&atomList[i], gensym(lua_tostring(L, -1)));
        } else {
            return luaL_error(L, "Table contains unsupported value type");
        }
        lua_pop(L, 1);
    }
    object_method_typed(symbol->s_thing, gensym("list"), listSize, atomList.data(), nullptr);
    return 0;
}

// ─────────────────────────────────────
static const luaL_Reg max_funcs[] = {

    // Log
    {"post", max_Post},
    {"error", max_Error},

    // Max
    {"send_bang", max_sendBang},
    {"send_float", max_sendFloat},
    {"send_symbol", max_sendSymbol},
    {"send_list", max_sendList},

    // Sentinela
    {NULL, NULL}};

// ─────────────────────────────────────
int luaopen_max(lua_State *L) {
    luaL_newlib(L, max_funcs);
    return 1;
}

#endif
