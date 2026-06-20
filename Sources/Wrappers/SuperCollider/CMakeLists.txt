cmake_minimum_required(VERSION 3.15)
project(SCOpenScofo C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cpmaddpackage(URI "gh:supercollider/supercollider#Version-3.14.1" NAME supercollider DOWNLOAD_ONLY
              YES)

file(GLOB_RECURSE OPENSCOFO_SC_SOURCES "*.cpp" "*.h" "*.hpp")
add_library(SCOpenScofo MODULE ${OPENSCOFO_SC_SOURCES})
target_include_directories(
    SCOpenScofo SYSTEM PRIVATE "${supercollider_SOURCE_DIR}/include/plugin_interface"
                               "${supercollider_SOURCE_DIR}/include/common")
target_link_libraries(SCOpenScofo PRIVATE OpenScofo)

# ╭──────────────────────────────────────╮
# │         SuperCollider Target         │
# ╰──────────────────────────────────────╯
if(APPLE OR WIN32)
    set_target_properties(SCOpenScofo PROPERTIES PREFIX "" SUFFIX ".scx" OUTPUT_NAME "OpenScofo")
else()
    set_target_properties(SCOpenScofo PROPERTIES PREFIX "" SUFFIX ".so" OUTPUT_NAME "OpenScofo")
    target_link_options(SCOpenScofo PRIVATE -Wl,-allow-multiple-definition -flto=auto)
    target_compile_options(SCOpenScofo PRIVATE -flto=auto)
endif()

# ╭──────────────────────────────────────╮
# │         Installation/Data            │
# ╰──────────────────────────────────────╯
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/examples/1-score-follow.scd
               ${CMAKE_CURRENT_BINARY_DIR}/examples/1-score-follow.scd COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/examples/2-descriptors-input.scd
               ${CMAKE_CURRENT_BINARY_DIR}/examples/2-descriptors-input.scd COPYONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo.sc ${CMAKE_CURRENT_BINARY_DIR}/OpenScofo.sc
               COPYONLY)
