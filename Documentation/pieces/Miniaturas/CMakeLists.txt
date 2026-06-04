cmake_minimum_required(VERSION 3.25)
project("Miniaturas")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(FetchContent)

# ╭──────────────────────────────────────╮
# │               pd.cmake               │
# ╰──────────────────────────────────────╯
set(PDCMAKE_FILE
    "${CMAKE_BINARY_DIR}/pd.cmake"
    CACHE INTERNAL "pd cmake file path")
include("${PDCMAKE_FILE}")

# ╭──────────────────────────────────────╮
# │                Nanovg                │
# ╰──────────────────────────────────────╯
FetchContent_Declare(nanovg SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/nanovg")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/nanovg" "${CMAKE_CURRENT_BINARY_DIR}/nanovg-build")

# ╭──────────────────────────────────────╮
# │              Pd sources              │
# ╰──────────────────────────────────────╯
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -matomics -mbulk-memory")
set(PD4WEB_EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Externals/")

set(PD_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/pure-data/src/")
include("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/libpd.cmake")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/pure-data/src")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/")

add_compile_definitions(PDTHREADS PDINSTANCE)

add_compile_options(-Wno-everything)


# ╭──────────────────────────────────────╮
# │       Debug or Release options       │
# ╰──────────────────────────────────────╯
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "Building in Debug mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -flto -pthread -matomics -mbulk-memory -msimd128")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -flto -pthread -matomics -mbulk-memory -msimd128")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -flto -pthread -matomics -mbulk-memory -msimd128")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -flto -pthread -matomics -mbulk-memory -msimd128")
endif()

# ╭──────────────────────────────────────╮
# │          Pd4Web executable           │
# ╰──────────────────────────────────────╯
add_executable(pd4web "${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/pd4web.cpp"
                      "${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/pd4web_externals.cpp")
target_include_directories(pd4web PUBLIC "${nanovg_SOURCE_DIR}/src")

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Externals/pdlua")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Externals/pdlua")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Externals/pdlua/luas/lua")
endif()

target_include_directories(pd4web PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/pure-data/src")
set_target_properties(pd4web PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/WebPatch")

target_link_libraries(
    pd4web
    PRIVATE embind
            libpd
            nanovg
            pdlua)
target_link_options(
    pd4web
    PRIVATE
    -sMODULARIZE=1
    -sEXPORT_NAME='Pd4WebModule'
    -sEXPORTED_RUNTIME_METHODS=["FS"]
    -sINITIAL_MEMORY=512MB
    -sEXPORT_ES6=${PD4WEB_AS_ES6}
    -sUSE_PTHREADS=1
    -sPTHREAD_POOL_SIZE=4
    -sWASMFS=0
    -sWASM=1
    -sWASM_WORKERS=1
    -sAUDIO_WORKLET=1
    -sUSE_WEBGL2=1
    -sMAX_WEBGL_VERSION=2
    -sMIN_WEBGL_VERSION=2
    -sOFFSCREENCANVAS_SUPPORT
    -sOFFSCREEN_FRAMEBUFFER)

# Externals includes
include("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Externals/pdlua.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Externals/xlab.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Externals/OpenScofo.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Externals/else.cmake")


# Project Externals Libraries
target_link_libraries(pd4web PRIVATE 
                                pdlua
                                openscofo_tilde
                                pv_tilde
                                saf.encoder_tilde
                                infinite.record_tilde
                                dragonfly_tilde
                                lb
                                order
                                args
                                default
                                routeall
                                lfreeze_tilde
                                saf.binaural_tilde
                                eq_tilde
                                mar_tilde)

# Preload Files
get_target_property(EMCC_LINK_FLAGS pd4web LINK_FLAGS)
if(NOT EMCC_LINK_FLAGS)
    set(EMCC_LINK_FLAGS "")
endif()
set_target_properties(pd4web PROPERTIES LINK_FLAGS "${EMCC_LINK_FLAGS}
	--preload-file \"${CMAKE_CURRENT_SOURCE_DIR}/WebPatch/index.pd@/index.pd\"
	--preload-file \"${CMAKE_CURRENT_SOURCE_DIR}/.tmp@/Libs/\"
	--preload-file \"${CMAKE_CURRENT_SOURCE_DIR}/Audios@/Audios/\"
	--preload-file \"${CMAKE_CURRENT_SOURCE_DIR}/Extras@/Extras/\"
	--preload-file \"${CMAKE_CURRENT_SOURCE_DIR}/Pd4Web/Gui@/Gui/\"")

