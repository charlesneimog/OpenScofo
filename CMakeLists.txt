cmake_minimum_required(VERSION 3.30)
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15")

project(OpenScofo VERSION 0.2.11 LANGUAGES C CXX
        HOMEPAGE_URL "https://charlesneimog.github.io/OpenScofo/"
        DESCRIPTION "OpenScofo is a Open Source Score Follower designed for contemporary music")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
set(SPDLOG_BUILD_SHARED OFF)
set(BUILD_EXAMPLES OFF)
set(TREE_SITTER_VERSION "v0.26.10")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fdiagnostics-color=always")

option(OPENSCOFO_BUILD_PD_OBJECT "Build PureData Object" OFF)
option(OPENSCOFO_BUILD_PY_MODULE "Build Python Module" OFF)
option(OPENSCOFO_BUILD_MAX_OBJECT "Build Max Object" OFF)
option(OPENSCOFO_BUILD_CSOUND_PLUGIN "Build CSound Plugin" OFF)
option(OPENSCOFO_BUILD_SC_PLUGIN "Build SuperCollider Plugin" OFF)
option(OPENSCOFO_BUILD_VAMP_PLUGIN "Build Vamp Plugin" OFF)
option(OPENSCOFO_BUILD_TESTS "Build Tests" ON)
option(OPENSCOFO_BUILD_ALL "Build All OpenScofo Tests and Wrappers" ON)
option(OPENSCOFO_ENABLE_SANITIZERS
       "Build OpenScofo with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
option(OPENSCOFO_BUILD_WITH_LUA "Build Lua Module embbeded on OpenScofo" ON)
option(UPDATE_OPENSCOFO_LANGUAGE "Update OpenScofo language grammar.js" ON)
option(OPENSCOFO_USE_BESSEL "Use Bessel function to calculate Kappa" ON)
option(OPENSCOFO_ENABLE_NATIVE_OPTIMIZATIONS "Build with native CPU optimizations" ON)

add_definitions(-DOPENSCOFO_VERSION="${PROJECT_VERSION}")
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

if(EMSCRIPTEN)
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -O3 -flto -pthread -matomics -mbulk-memory -msimd128 -sUSE_BOOST_HEADERS=1"
    )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -flto -pthread -matomics -mbulk-memory -msimd128")
elseif(WIN32)
    add_compile_definitions(_WIN32_WINNT=0x0601)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
endif()

if(OPENSCOFO_BUILD_ALL)
    set(OPENSCOFO_BUILD_PY_MODULE ON)
    set(OPENSCOFO_BUILD_PD_OBJECT ON)
    set(OPENSCOFO_BUILD_MAX_OBJECT ON)
    set(OPENSCOFO_BUILD_VAMP_PLUGIN ON)
    set(OPENSCOFO_BUILD_SC_PLUGIN ON)
    if(NOT WIN32)
        set(OPENSCOFO_BUILD_CSOUND_PLUGIN ON)
    endif()
endif()

# ╭──────────────────────────────────────╮
# │               PD CMAKE               │
# ╰──────────────────────────────────────╯
set(PDCMAKE_FILE "${CMAKE_BINARY_DIR}/pd.cmake")
set(PDCMAKE_VERSION "v2.0.17")

if(NOT EXISTS "${PDCMAKE_FILE}")
    file(DOWNLOAD
         "https://github.com/pure-data/pd.cmake/releases/download/${PDCMAKE_VERSION}/pd.cmake"
         "${PDCMAKE_FILE}" STATUS PDCMAKE_DOWNLOAD_STATUS SHOW_PROGRESS)
    list(GET PDCMAKE_DOWNLOAD_STATUS 0 PDCMAKE_DOWNLOAD_CODE)
    list(GET PDCMAKE_DOWNLOAD_STATUS 1 PDCMAKE_DOWNLOAD_MESSAGE)
    if(NOT PDCMAKE_DOWNLOAD_CODE EQUAL 0)
        message(FATAL_ERROR "Failed to download pd.cmake: ${PDCMAKE_DOWNLOAD_MESSAGE}")
    endif()
endif()

include("${PDCMAKE_FILE}")
# ╭──────────────────────────────────────╮
# │             CPM Package              │
# ╰──────────────────────────────────────╯
set(CPM_FILE ${CMAKE_BINARY_DIR}/CPM.cmake)
set(CPM_VERSION "0.42.1")
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

set(GRAMMAR_BUILD_DIR "${CMAKE_BINARY_DIR}/openscofo_grammar")
set(GRAMMAR_JS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/Sources/Language/grammar.js")
set(GRAMMAR_JS_BUILD "${GRAMMAR_BUILD_DIR}/grammar.js")

set(GRAMMAR_JSON_SRC "${CMAKE_CURRENT_SOURCE_DIR}/Resources/OpenScofo Language/tree-sitter.json")
set(GRAMMAR_JSON_BUILD "${GRAMMAR_BUILD_DIR}/tree-sitter.json")

set(PARSER_GENERATED_OUTPUT "${GRAMMAR_BUILD_DIR}/src/parser.c")
set(PARSER_COMMITTED_OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/Sources/Language/src/parser.c")
set(PARSER_BYPRODUCTS "${GRAMMAR_BUILD_DIR}/src/node-types.json")

set(OPENSCOFO_ONLINE_EDITOR_WASM
    "${CMAKE_CURRENT_SOURCE_DIR}/Score-Editor/tree-sitter/tree-sitter-openscofo.wasm")

if(UPDATE_OPENSCOFO_LANGUAGE)
    find_program(TREE_SITTER_EXECUTABLE NAMES tree-sitter DOC "tree-sitter CLI executable")
    if(NOT TREE_SITTER_EXECUTABLE)
        message(
            FATAL_ERROR "tree-sitter CLI not found. Install it with: npm install -g tree-sitter-cli"
        )
    endif()

    add_custom_command(
        OUTPUT "${PARSER_GENERATED_OUTPUT}"
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

    add_custom_target(generate_tree_sitter_parser DEPENDS "${PARSER_GENERATED_OUTPUT}")
    set(PARSER_OUTPUT "${PARSER_GENERATED_OUTPUT}")
else()
    set(PARSER_OUTPUT "${PARSER_COMMITTED_OUTPUT}")
endif()

# Biblioteca
add_library(score_parser STATIC "${PARSER_OUTPUT}")
set_target_properties(score_parser PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(UPDATE_OPENSCOFO_LANGUAGE)
    add_dependencies(score_parser generate_tree_sitter_parser)
endif()

# Custom command to generate the WASM file if it doesn't exist
if(NOT EXISTS "${OPENSCOFO_ONLINE_EDITOR_WASM}" OR ${UPDATE_OPENSCOFO_LANGUAGE})
    add_custom_command(
        OUTPUT "${OPENSCOFO_ONLINE_EDITOR_WASM}"
        COMMAND ${CMAKE_COMMAND} -E chdir "${GRAMMAR_BUILD_DIR}" ${TREE_SITTER_EXECUTABLE} build
                --wasm
        COMMAND ${CMAKE_COMMAND} -E copy "${GRAMMAR_BUILD_DIR}/tree-sitter-openscofo.wasm"
                "${CMAKE_CURRENT_SOURCE_DIR}/Score-Editor/tree-sitter/tree-sitter-openscofo.wasm"
        DEPENDS "${GRAMMAR_JS_SRC}"
        COMMENT "Generating Tree-sitter parser for OpenScofo Online Editor and WASM..."
        VERBATIM)

    set(TS_URL
        "https://github.com/tree-sitter/tree-sitter/releases/download/${TREE_SITTER_VERSION}/web-tree-sitter.tar.gz"
    )

    set(TS_ARCHIVE "${CMAKE_BINARY_DIR}/web-tree-sitter.tar.gz")
    set(TS_TEMP_DIR "${CMAKE_BINARY_DIR}/web-tree-sitter-temp")
    set(TS_DEST "${CMAKE_CURRENT_SOURCE_DIR}/Score-Editor/tree-sitter")

    file(MAKE_DIRECTORY "${TS_TEMP_DIR}")
    file(MAKE_DIRECTORY "${TS_DEST}")

    # Download archive if missing
    if(NOT EXISTS "${TS_ARCHIVE}")
        file(DOWNLOAD "${TS_URL}" "${TS_ARCHIVE}" STATUS download_status)
        list(GET download_status 0 status_code)
        if(NOT status_code EQUAL 0)
            message(FATAL_ERROR "Failed to download web-tree-sitter archive")
        endif()
    endif()

    # Extract
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf "${TS_ARCHIVE}"
                    WORKING_DIRECTORY "${TS_TEMP_DIR}")

    # Verify files exist
    if(NOT EXISTS "${TS_TEMP_DIR}/web-tree-sitter.js" OR NOT EXISTS
                                                         "${TS_TEMP_DIR}/web-tree-sitter.wasm")
        message(FATAL_ERROR "web-tree-sitter extraction failed")
    endif()

    # Copy required runtime files
    file(COPY "${TS_TEMP_DIR}/web-tree-sitter.js" "${TS_TEMP_DIR}/web-tree-sitter.wasm"
         DESTINATION "${TS_DEST}")

    add_custom_target(build_online_editor_wasm ALL DEPENDS "${OPENSCOFO_ONLINE_EDITOR_WASM}")

endif()

# ╭──────────────────────────────────────╮
# │                 Lua                  │
# ╰──────────────────────────────────────╯
if(OPENSCOFO_BUILD_WITH_LUA)
    cpmaddpackage("gh:lua/lua@5.5.0")
    if(lua_ADDED)
        file(GLOB lua_sources ${lua_SOURCE_DIR}/*.c)
        add_library(OpenScofo_Lua STATIC "${lua_SOURCE_DIR}/onelua.c")
        add_definitions("-DMAKE_LIB")
        target_compile_definitions(OpenScofo_Lua PUBLIC MAKE_LIB)
        target_include_directories(OpenScofo_Lua PUBLIC $<BUILD_INTERFACE:${lua_SOURCE_DIR}>)
        set_target_properties(OpenScofo_Lua PROPERTIES POSITION_INDEPENDENT_CODE ON)

        if(APPLE)
            target_compile_definitions(OpenScofo_Lua PRIVATE LUA_USE_MACOSX)
        elseif(WIN32)
            target_compile_definitions(OpenScofo_Lua PRIVATE LUA_USE_WINDOWS)
        elseif(UNIX)
            target_compile_definitions(OpenScofo_Lua PRIVATE LUA_USE_LINUX)
        endif()

    endif()
    add_definitions(-DOPENSCOFO_LUA)
else()
    message(STATUS "Lua not included")
endif()

# ╭──────────────────────────────────────╮
# │                FFTW3                 │
# ╰──────────────────────────────────────╯
cpmaddpackage(NAME MIPP GITHUB_REPOSITORY "aff3ct/MIPP" GIT_TAG "master" DOWNLOAD_ONLY YES)
set(MIPP_INCLUDE_DIRS "${MIPP_SOURCE_DIR}/include")

if(WIN32)
    set(PFFFT_DISABLE_LINK_WITH_M ON CACHE BOOL "Do not link libm on Windows" FORCE)
endif()
cpmaddpackage(NAME pffft GITHUB_REPOSITORY "marton78/pffft" GIT_TAG "v1.1.0" OPTIONS
              "CMAKE_POSITION_INDEPENDENT_CODE ON")
cpmaddpackage("gh:charlesneimog/onsetsds_tilde#main")
cpmaddpackage("gh:gabime/spdlog@1.17.0")
cpmaddpackage("gh:charlesneimog/libonnx#master")

target_compile_definitions(spdlog PUBLIC SPDLOG_NO_EXCEPTIONS)

# ╭──────────────────────────────────────╮
# │                OpenScofo             │
# ╰──────────────────────────────────────╯
list(APPEND ScofoSrc Sources/OpenScofo/forward.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/mir.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/onnx.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/score.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/lua.cpp)
list(APPEND ScofoSrc Sources/OpenScofo/OpenScofo.cpp)
add_library(OpenScofo STATIC ${ScofoSrc})

target_include_directories(OpenScofo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Sources/OpenScofo")
set_target_properties(OpenScofo PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(
    OpenScofo
    PUBLIC # score parser
           tree-sitter
           score_parser
           # logging
           spdlog
           # mir
           onsetsds
           onnx
           PFFFT::PFFFT)

if(OPENSCOFO_USE_BESSEL)
    target_compile_definitions(OpenScofo PRIVATE OPENSCOFO_USE_BESSEL=1)
endif()

if(OPENSCOFO_BUILD_WITH_LUA)
    target_link_libraries(OpenScofo PUBLIC OpenScofo_Lua)
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
               -ffast-math
               -DNDEBUG
               -flto=auto
               -ffunction-sections
               -fdata-sections
               >)
    target_compile_options(OpenScofo PRIVATE -fno-exceptions)

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
               -Weffc++
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
    target_compile_options(OpenScofo PRIVATE -fno-exceptions)

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
               -Weffc++
               -Winit-self
               -fstack-protector-strong
               -ftrapv
               -fstrict-aliasing
               $<$<BOOL:${OPENSCOFO_ENABLE_SANITIZERS}>:
               -fsanitize=address
               -fsanitize=undefined
               >
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
    target_link_options(
        OpenScofo PUBLIC
        "$<$<BOOL:${OPENSCOFO_ENABLE_SANITIZERS}>:-fsanitize=address;-fsanitize=undefined>")
    target_compile_options(OpenScofo PRIVATE -fno-exceptions)

else() # Linux (GCC or Clang)
    target_compile_options(
        OpenScofo
        PUBLIC -Wall
               -Wextra
               -Wshadow
               -Wformat=2
               -Wnull-dereference
               $<$<CONFIG:Debug>:
               -ggdb
               -Og
               -fno-omit-frame-pointer
               >
               $<$<CONFIG:RelWithDebInfo>:
               -O2
               -ggdb
               -fno-omit-frame-pointer
               >
               $<$<CONFIG:Release>:
               -O3
               -ffast-math
               -DNDEBUG
               -flto=auto
               -ffunction-sections
               -fdata-sections
               >)

    if(OPENSCOFO_ENABLE_NATIVE_OPTIMIZATIONS)
        target_compile_options(OpenScofo PUBLIC $<$<CONFIG:Release>:-march=native>)
    endif()

    target_link_options(OpenScofo PUBLIC $<$<CONFIG:Release>:-flto=auto;-Wl,--gc-sections;-s >)
    target_compile_options(OpenScofo PRIVATE -fno-exceptions)
endif()

# ╭──────────────────────────────────────╮
# │               Language               │
# ╰──────────────────────────────────────╯
if(NOT EMSCRIPTEN)
    add_library(tree-sitter-openscofo SHARED
                ${CMAKE_CURRENT_BINARY_DIR}/openscofo_grammar/src/parser.c)
endif()

# ╭──────────────────────────────────────╮
# │            Python Module             │
# ╰──────────────────────────────────────╯
if(OPENSCOFO_BUILD_PY_MODULE AND NOT EMSCRIPTEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Sources/Wrappers/Python")
endif()

# ╭──────────────────────────────────────╮
# │           PureData Object            │
# ╰──────────────────────────────────────╯
if(OPENSCOFO_BUILD_PD_OBJECT)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Sources/Wrappers/PureData")
endif()

# ╭──────────────────────────────────────╮
# │             Max Obj                  │
# ╰──────────────────────────────────────╯
if(OPENSCOFO_BUILD_MAX_OBJECT AND NOT LINUX AND NOT EMSCRIPTEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Sources/Wrappers/Max")
endif()

# ╭──────────────────────────────────────╮
# │            Super Collider            │
# ╰──────────────────────────────────────╯
if(OPENSCOFO_BUILD_SC_PLUGIN AND NOT EMSCRIPTEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Sources/Wrappers/SuperCollider")
endif()

# ╭──────────────────────────────────────╮
# │              CSound                  │
# ╰──────────────────────────────────────╯
if(OPENSCOFO_BUILD_CSOUND_PLUGIN AND NOT EMSCRIPTEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Sources/Wrappers/CSound")
endif()

# ╭──────────────────────────────────────╮
# │            Vamp Obj                  │
# ╰──────────────────────────────────────╯
if(OPENSCOFO_BUILD_VAMP_PLUGIN AND NOT EMSCRIPTEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Sources/Wrappers/Vamp")
endif()

# ╭──────────────────────────────────────╮
# │            WASM Obj                  │
# ╰──────────────────────────────────────╯
if(EMSCRIPTEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Sources/Wrappers/Wasm")
endif()

# ╭──────────────────────────────────────╮
# │                Tests                 │
# ╰──────────────────────────────────────╯
if(OPENSCOFO_BUILD_TESTS AND NOT EMSCRIPTEN)
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Tests")
endif()

if(LINUX)
    include(InstallRequiredSystemLibraries)
    set(CPACK_GENERATOR "STGZ")
    set(CPACK_PACKAGE_NAME "OpenScofo")
    set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
    include(CPack)
endif()
