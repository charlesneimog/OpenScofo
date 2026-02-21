cmake_minimum_required(VERSION 3.25)
project(OScofo)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS OFF)

option(BUILD_PY_MODULE "Build Python Module" ON)
option(BUILD_PD_OBJECT "Build PureData Object" ON)
option(BUILD_MAX_OBJECT "Build Max Object" OFF)
option(BUILD_ALL "Build Python Module and PureData Object" ON)
option(BUILD_WITH_LUA "Build Lua Module embbeded on OScofo" ON)
option(UPDATE_OSCOFO_LANGUAGE "Update OScofo language grammar.js" OFF)

# check if it is EMSCRIPTEN
if(EMSCRIPTEN)
    message(STATUS "Emscripten compilation")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s USE_BOOST_HEADERS=1")
endif()

if(BUILD_ALL)
    set(BUILD_PY_MODULE ON)
    set(BUILD_PD_OBJECT ON)
    if(NOT LINUX)
        set(BUILD_MAX_OBJECT ON)
    endif()
    set(BUILD_TESTS ON)
endif()

if(BUILD_ALL_OBJECTS)
    set(BUILD_PD_OBJECT ON)
    if(NOT LINUX)
        set(BUILD_MAX_OBJECT ON)
    endif()
endif()

if(BUILD_PD_OBJECT)
    set(PDCMAKE_FILE ${CMAKE_BINARY_DIR}/pd.cmake)
    set(PDCMAKE_VERSION "v0.2.12")
    if(NOT EXISTS "${PDCMAKE_FILE}")
        file(
            DOWNLOAD
            https://raw.githubusercontent.com/pure-data/pd.cmake/refs/heads/naming-conventions/pd.cmake
            ${PDCMAKE_FILE})
    endif()
    include(${PDCMAKE_FILE})
endif()

# ╭──────────────────────────────────────╮
# │             CPM Package              │
# ╰──────────────────────────────────────╯
set(CPM_FILE ${CMAKE_BINARY_DIR}/CPM.cmake)
set(CPM_VERSION "0.42.0")
if(NOT EXISTS "${CPM_FILE}")
    file(DOWNLOAD
         "https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_VERSION}/CPM.cmake"
         ${CPM_FILE})
endif()
include(${CPM_FILE})

# ╭──────────────────────────────────────╮
# │             Tree-Sitter              │
# ╰──────────────────────────────────────╯
cpmaddpackage("gh:tree-sitter/tree-sitter#v0.26.3")

find_program(
    TREE_SITTER_EXECUTABLE
    NAMES tree-sitter
    DOC "tree-sitter CLI executable")
if(NOT TREE_SITTER_EXECUTABLE)
    message(
        FATAL_ERROR "tree-sitter CLI not found. Install it with: npm install -g tree-sitter-cli")
endif()

set(GRAMMAR_BUILD_DIR "${CMAKE_BINARY_DIR}/oscofo_grammar")
set(GRAMMAR_JS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/Sources/OScofo/grammar.js")
set(GRAMMAR_JS_BUILD "${GRAMMAR_BUILD_DIR}/grammar.js")
set(PARSER_OUTPUT "${GRAMMAR_BUILD_DIR}/src/parser.c")
set(PARSER_BYPRODUCTS "${GRAMMAR_BUILD_DIR}/src/node-types.json")

add_custom_command(
    OUTPUT "${PARSER_OUTPUT}"
    BYPRODUCTS "${PARSER_BYPRODUCTS}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${GRAMMAR_BUILD_DIR}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${GRAMMAR_BUILD_DIR}/src"
    COMMAND ${CMAKE_COMMAND} -E copy "${GRAMMAR_JS_SRC}" "${GRAMMAR_JS_BUILD}"
    COMMAND ${CMAKE_COMMAND} -E chdir "${GRAMMAR_BUILD_DIR}" ${TREE_SITTER_EXECUTABLE} generate
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    DEPENDS "${GRAMMAR_JS_SRC}"
    COMMENT "Copying grammar.js and generating Tree-sitter parser..."
    VERBATIM)

# Target custom
add_custom_target(generate_tree_sitter_parser DEPENDS "${PARSER_OUTPUT}")

# Biblioteca
add_library(score_parser STATIC "${PARSER_OUTPUT}")
set_target_properties(score_parser PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_dependencies(score_parser generate_tree_sitter_parser)

# ╭──────────────────────────────────────╮
# │                 Lua                  │
# ╰──────────────────────────────────────╯
if(BUILD_WITH_LUA)
    cpmaddpackage("gh:lua/lua@5.4.0")
    if(lua_ADDED)
        file(GLOB lua_sources ${lua_SOURCE_DIR}/*.c)
        add_library(lua STATIC "${lua_SOURCE_DIR}/onelua.c")
        add_definitions("-DMAKE_LIB")
        target_compile_definitions(lua PUBLIC MAKE_LIB)
        target_include_directories(lua PUBLIC $<BUILD_INTERFACE:${lua_SOURCE_DIR}>)
        set_target_properties(lua PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif()
    add_definitions(-DOSCOFO_LUA)
else()
    message(STATUS "Lua not included")
endif()

# ╭──────────────────────────────────────╮
# │                FFTW3                 │
# ╰──────────────────────────────────────╯
set(FFTW3_FILE ${CMAKE_BINARY_DIR}/fftw-3.3.10.tar.gz)
if(NOT EXISTS ${FFTW3_FILE})
    message(STATUS "Downloading FFTW3")
    file(
        DOWNLOAD https://www.fftw.org/fftw-3.3.10.tar.gz ${FFTW3_FILE}
        STATUS DOWNLOAD_STATUS)
endif()

set(BUILD_SHARED_LIBS
    OFF
    CACHE BOOL "Build static libs" FORCE)
set(BUILD_TESTS
    OFF
    CACHE BOOL "Disable FFTW tests" FORCE)
set(ENABLE_FLOAT
    OFF
    CACHE BOOL "Enable single precision" FORCE)

file(ARCHIVE_EXTRACT INPUT ${CMAKE_BINARY_DIR}/fftw-3.3.10.tar.gz DESTINATION ${CMAKE_BINARY_DIR}/)
add_subdirectory(${CMAKE_BINARY_DIR}/fftw-3.3.10 ${CMAKE_BINARY_DIR}/fftw3 EXCLUDE_FROM_ALL)
set_target_properties(fftw3 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
set_target_properties(fftw3 PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ╭──────────────────────────────────────╮
# │               LibONNX                │
# ╰──────────────────────────────────────╯
cpmaddpackage("gh:charlesneimog/libonnx#master")
cpmaddpackage("gh:charlesneimog/onsetsds_tilde#main")
cpmaddpackage("gh:ThePhD/sol2#v3.5.0")

# ╭──────────────────────────────────────╮
# │                OScofo                │
# ╰──────────────────────────────────────╯
list(APPEND ScofoSrc Sources/OScofo/mdp.cpp)
list(APPEND ScofoSrc Sources/OScofo/mir.cpp)
list(APPEND ScofoSrc Sources/OScofo/score.cpp)
list(APPEND ScofoSrc Sources/OScofo/lua.cpp)
list(APPEND ScofoSrc Sources/OScofo/OScofo.cpp)
add_library(OScofo STATIC ${ScofoSrc})

target_include_directories(OScofo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Sources/")
set_target_properties(OScofo PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(
    OScofo
    PUBLIC sol2
           fftw3
           tree-sitter
           score_parser
           lua
           onnx
           onsetsds)

if(LINUX)
    target_link_options(OScofo PRIVATE -static-libstdc++ -Wl,--exclude-libs,ALL -static,)
endif()

if(NOT EMSCRIPTEN AND NOT MSVC)
    target_compile_options(
        OScofo
        PRIVATE -Wall # Enable all common warnings
                -Wextra # Enable extra warnings
                -Wshadow # Warn about shadowed variable declarations
                -Wformat=2 # Enable stricter format string checks
                -Wnull-dereference # Warn on dereferences of null pointers
                -march=native
                # Debug
                $<$<CONFIG:Debug>:
                -O2
                -g
                -fno-omit-frame-pointer
                >
                # Release
                $<$<CONFIG:Release>:
                -DNDEBUG
                -O3
                -march=native
                >)
endif()

# ╭──────────────────────────────────────╮
# │            Python Module             │
# ╰──────────────────────────────────────╯
if(BUILD_PY_MODULE)
    add_subdirectory(Sources/Python)
endif()

# ╭──────────────────────────────────────╮
# │           PureData Object            │
# ╰──────────────────────────────────────╯
if(BUILD_PD_OBJECT)
    add_subdirectory(Sources/PureData)
endif()

# ╭──────────────────────────────────────╮
# │             Max Obj                  │
# ╰──────────────────────────────────────╯
if(BUILD_MAX_OBJECT AND NOT LINUX)
    add_subdirectory(Sources/Max)
endif()

# ╭──────────────────────────────────────╮
# │                Tests                 │
# ╰──────────────────────────────────────╯
if(BUILD_TESTS)
    add_subdirectory(Tests)
endif()
