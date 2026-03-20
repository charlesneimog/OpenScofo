project(o.scofo~ C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE OPENSCOFO_PDOBJ_SOURCES "*.cpp" "*.h")
pd_add_external(o.scofo~ "${OPENSCOFO_PDOBJ_SOURCES}" LINK_LIBRARIES OpenScofo)

# ╭──────────────────────────────────────╮
# │           PureData Objects           │
# ╰──────────────────────────────────────╯
if(UNIX AND NOT APPLE AND NOT EMSCRIPTEN)
    target_link_options(o.scofo_tilde PRIVATE -Wl,-allow-multiple-definition)
    target_compile_options(o.scofo_tilde PRIVATE -flto=auto)
    target_link_options(o.scofo_tilde PRIVATE -flto=auto)
endif()

pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/README.md")
pd_add_datafile(o.scofo_tilde "${CMAKE_CURRENT_SOURCE_DIR}/o.scofo~-help.pd")
pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.mp3")
pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt")