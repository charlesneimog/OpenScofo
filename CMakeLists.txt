cmake_minimum_required(VERSION 3.25)
project(OpenScofo VERSION 0.1.4)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS OFF)
set(SPDLOG_BUILD_SHARED OFF)
set(BUILD_EXAMPLES OFF)
set(CMAKE_OSX_DEPLOYMENT_TARGET 13.3)
set(TREE_SITTER_VERSION "v0.26.6")

option(BUILD_PY_MODULE "Build Python Module" OFF)
option(BUILD_PD_OBJECT "Build PureData Object" ON)
option(BUILD_MAX_OBJECT "Build Max Object" OFF)
option(BUILD_ALL "Build Python Module and PureData Object" ON)
option(BUILD_WITH_LUA "Build Lua Module embbeded on OpenScofo" ON)
option(UPDATE_OSCOFO_LANGUAGE "Update OpenScofo language grammar.js" OFF)
option(OSCOFO_BUILD_TESTS "Run Tests" OFF)

if(BUILD_ALL)
    set(BUILD_PY_MODULE ON)
    set(BUILD_PD_OBJECT ON)
    if(NOT LINUX)
        set(BUILD_MAX_OBJECT ON)
    endif()
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
        file(DOWNLOAD https://raw.githubusercontent.com/pure-data/pd.cmake/refs/heads/main/pd.cmake
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
cpmaddpackage("gh:tree-sitter/tree-sitter#${TREE_SITTER_VERSION}")

find_program(TREE_SITTER_EXECUTABLE NAMES tree-sitter DOC "tree-sitter CLI executable")
if(NOT TREE_SITTER_EXECUTABLE)
    message(
        FATAL_ERROR "tree-sitter CLI not found. Install it with: npm install -g tree-sitter-cli")
endif()

set(GRAMMAR_BUILD_DIR "${CMAKE_BINARY_DIR}/oscofo_grammar")
set(GRAMMAR_JS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/Sources/Language/grammar.js")
set(GRAMMAR_JS_BUILD "${GRAMMAR_BUILD_DIR}/grammar.js")

set(GRAMMAR_JSON_SRC "${CMAKE_CURRENT_SOURCE_DIR}/Resources/OScofo Language/tree-sitter.json")
set(GRAMMAR_JSON_BUILD "${GRAMMAR_BUILD_DIR}/tree-sitter.json")

set(PARSER_OUTPUT "${GRAMMAR_BUILD_DIR}/src/parser.c")
set(PARSER_BYPRODUCTS "${GRAMMAR_BUILD_DIR}/src/node-types.json")

set(OSCOFO_ONLINE_EDITOR_WASM
    "${CMAKE_CURRENT_SOURCE_DIR}/Score-Editor/tree-sitter/tree-sitter-openscofo.wasm")

add_custom_command(
    OUTPUT "${PARSER_OUTPUT}"
    BYPRODUCTS "${PARSER_BYPRODUCTS}" "${GRAMMAR_BUILD_DIR}/tree-sitter-${GRAMMAR_NAME}.wasm"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${GRAMMAR_BUILD_DIR}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${GRAMMAR_BUILD_DIR}/src"
    COMMAND ${CMAKE_COMMAND} -E copy "${GRAMMAR_JS_SRC}" "${GRAMMAR_JS_BUILD}"
    COMMAND ${CMAKE_COMMAND} -E copy "${GRAMMAR_JSON_SRC}" "${GRAMMAR_JSON_BUILD}"
    COMMAND ${CMAKE_COMMAND} -E chdir "${GRAMMAR_BUILD_DIR}" ${TREE_SITTER_EXECUTABLE} generate
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    DEPENDS "${GRAMMAR_JS_SRC}"
    COMMENT "Generating Tree-sitter parser and WASM..."
    VERBATIM)

# Target custom
add_custom_target(generate_tree_sitter_parser DEPENDS "${PARSER_OUTPUT}")

# Biblioteca
add_library(score_parser STATIC "${PARSER_OUTPUT}")
set_target_properties(score_parser PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_dependencies(score_parser generate_tree_sitter_parser)

# Custom command to generate the WASM file if it doesn't exist
if(NOT EXISTS "${OSCOFO_ONLINE_EDITOR_WASM}")
    add_custom_command(
        OUTPUT "${OSCOFO_ONLINE_EDITOR_WASM}"
        COMMAND ${CMAKE_COMMAND} -E chdir "${GRAMMAR_BUILD_DIR}" ${TREE_SITTER_EXECUTABLE} build
                --wasm
        COMMAND ${CMAKE_COMMAND} -E copy "${GRAMMAR_BUILD_DIR}/tree-sitter-openscofo.wasm"
                "${CMAKE_CURRENT_SOURCE_DIR}/Score-Editor/tree-sitter/tree-sitter-openscofo.wasm"
        DEPENDS "${GRAMMAR_JS_SRC}"
        COMMENT "Generating Tree-sitter parser for OScofo Online Editor and WASM..."
        VERBATIM)

    set(TS_URL
        "https://github.com/tree-sitter/tree-sitter/releases/download/${TREE_SITTER_VERSION}/web-tree-sitter.tar.gz"
    )
    set(TS_ARCHIVE "${CMAKE_BINARY_DIR}/web-tree-sitter.tar.gz")
    set(TS_URL "https://example.com/path/to/web-tree-sitter.tar.gz")
    set(TS_TEMP_DIR "${CMAKE_BINARY_DIR}/web-tree-sitter-temp")
    set(TS_DEST "${CMAKE_CURRENT_SOURCE_DIR}/Score-Editor/tree-sitter")

    # Download the archive
    if(NOT EXISTS "${TS_ARCHIVE}")
        file(DOWNLOAD "${TS_URL}" "${TS_ARCHIVE}" SHOW_PROGRESS STATUS download_status)
        list(GET download_status 0 status_code)
        if(NOT status_code EQUAL 0)
            message(FATAL_ERROR "Failed to download web-tree-sitter archive: ${download_status}")
        endif()
    endif()

    # Extract all files to temp dir
    file(MAKE_DIRECTORY "${TS_TEMP_DIR}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf "${TS_ARCHIVE}"
                    WORKING_DIRECTORY "${TS_TEMP_DIR}")

    # Copy only the required files
    file(COPY "${TS_TEMP_DIR}/web-tree-sitter.wasm" "${TS_TEMP_DIR}/web-tree-sitter.js"
         DESTINATION "${TS_DEST}")

    # Create a target that depends on this output
    add_custom_target(build_online_editor_wasm ALL DEPENDS "${OSCOFO_ONLINE_EDITOR_WASM}")
endif()

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

        if(LINUX)
            target_compile_definitions(lua PUBLIC LUA_USE_LINUX)
        elseif(APPLE)
            target_compile_definitions(lua PUBLIC LUA_USE_MACOSX)
        endif()

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
    file(DOWNLOAD https://www.fftw.org/fftw-3.3.10.tar.gz ${FFTW3_FILE} STATUS DOWNLOAD_STATUS)
endif()

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build static libs" FORCE)
set(ENABLE_FLOAT OFF CACHE BOOL "Enable single precision" FORCE)

file(ARCHIVE_EXTRACT INPUT ${CMAKE_BINARY_DIR}/fftw-3.3.10.tar.gz DESTINATION ${CMAKE_BINARY_DIR}/)
add_subdirectory(${CMAKE_BINARY_DIR}/fftw-3.3.10 ${CMAKE_BINARY_DIR}/fftw3 EXCLUDE_FROM_ALL)
set_target_properties(fftw3 PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
set_target_properties(fftw3 PROPERTIES POSITION_INDEPENDENT_CODE ON)
include_directories("${CMAKE_BINARY_DIR}/fftw-3.3.10/api")

# ╭──────────────────────────────────────╮
# │               LibONNX                │
# ╰──────────────────────────────────────╯
cpmaddpackage("gh:charlesneimog/libonnx#master")
cpmaddpackage("gh:charlesneimog/onsetsds_tilde#main")
cpmaddpackage("gh:gabime/spdlog@1.17.0")

set_target_properties(spdlog PROPERTIES POSITION_INDEPENDENT_CODE ON)

# ╭──────────────────────────────────────╮
# │                OpenScofo             │
# ╰──────────────────────────────────────╯
list(APPEND ScofoSrc Sources/OpenScofo/mdp.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/mir.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/score.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/lua.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/OpenScofo.cpp)
add_library(OpenScofo STATIC ${ScofoSrc})

target_include_directories(OpenScofo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Sources/")
target_compile_definitions(OpenScofo PUBLIC OPENSCOFO_VERSION="${PROJECT_VERSION}")
set_target_properties(OpenScofo PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(OpenScofo PUBLIC fftw3 tree-sitter score_parser onnx spdlog onsetsds)

if(BUILD_WITH_LUA)
    target_link_libraries(OpenScofo PUBLIC lua)
endif()

if(EMSCRIPTEN)
    target_compile_options(
        OpenScofo
        PUBLIC -Wall
               -Wextra
               -Wshadow
               $<$<CONFIG:Debug>:
               -O0
               -g
               >
               $<$<CONFIG:Release>:
               -O3
               -DNDEBUG
               -flto
               >)

elseif(MSVC)
    target_compile_options(
        OpenScofo
        PUBLIC /W4
               /permissive-
               /Zc:__cplusplus
               $<$<CONFIG:Debug>:
               /Od
               /Zi
               >
               $<$<CONFIG:Release>:
               /O2
               /DNDEBUG
               /GL
               >)

    target_link_options(OpenScofo PUBLIC $<$<CONFIG:Release>: /LTCG >)

elseif(WIN32)
    target_compile_options(
        OpenScofo
        PUBLIC -Wall
               -Wextra
               -Wshadow
               -Wformat=2
               -Wnull-dereference
               $<$<CONFIG:Debug>:
               -O0
               -g
               -fno-omit-frame-pointer
               >
               $<$<CONFIG:Release>:
               -O3
               -DNDEBUG
               -flto=auto
               -ffunction-sections
               -fdata-sections
               >)

    target_link_options(OpenScofo PUBLIC $<$<CONFIG:Release>: -flto=auto -Wl,--gc-sections -s >)

elseif(APPLE)
    find_package(Boost REQUIRED)
    target_link_libraries(OpenScofo PRIVATE Boost::headers)
    target_compile_definitions(OpenScofo PRIVATE USE_BOOST_BESSEL)
    target_compile_options(
        OpenScofo
        PUBLIC -Wall
               -Wextra
               -Wshadow
               -Wformat=2
               -Wnull-dereference
               $<$<CONFIG:Debug>:
               -O0
               -g
               -fno-omit-frame-pointer
               >
               $<$<CONFIG:Release>:
               -O3
               -DNDEBUG
               -flto
               -ffunction-sections
               -fdata-sections
               >)

    target_link_options(OpenScofo PUBLIC $<$<CONFIG:Release>: -flto -Wl,-dead_strip >)

else() # Linux (GCC or Clang)

    target_compile_options(
        OpenScofo
        PUBLIC -Wall
               -Wextra
               -Wshadow
               -Wformat=2
               -Wnull-dereference
               -march=native
               $<$<CONFIG:Debug>:
               -O0
               -g
               -fno-omit-frame-pointer
               >
               $<$<CONFIG:RelWithDebInfo>:
               -g
               -fno-omit-frame-pointer
               >
               $<$<CONFIG:Release>:
               -O3
               -DNDEBUG
               -flto=auto
               -ffunction-sections
               -fdata-sections
               >)

    target_link_options(OpenScofo PUBLIC $<$<CONFIG:Release>: -flto=auto -Wl,--gc-sections -s >)

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
if(OSCOFO_BUILD_TESTS)
    add_subdirectory(Tests)
endif()
