project(o.scofo~ C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET "14")

file(GLOB_RECURSE OBJ_SOURCES "*.cpp" "*.h")
pd_add_external(o.scofo~ "${OBJ_SOURCES}" LINK_LIBRARIES OScofo)

# ╭──────────────────────────────────────╮
# │           PureData Objects           │
# ╰──────────────────────────────────────╯
if(NOT APPLE AND NOT EMSCRIPTEN)
    target_compile_options(o.scofo_tilde PRIVATE -static-libstdc++ -static-libgcc
                                                 -Wl,-allow-multiple-definition -static)
    target_link_options(o.scofo_tilde PRIVATE -static-libstdc++ -static-libgcc
                        -Wl,-allow-multiple-definition)
endif()

pd_add_datafile(o.scofo_tilde "${CMAKE_CURRENT_SOURCE_DIR}/../../README.md")
pd_add_datafile(o.scofo_tilde "${CMAKE_CURRENT_SOURCE_DIR}/../../Resources/o.scofo~-help.pd")
pd_add_datafile(o.scofo_tilde "${CMAKE_CURRENT_SOURCE_DIR}/../../Resources/Audios/bwv-1013.mp3")
pd_add_datafile(o.scofo_tilde "${CMAKE_CURRENT_SOURCE_DIR}/../../Resources/tests/bwv-1013.txt")
