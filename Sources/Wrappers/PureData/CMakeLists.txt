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
elseif(WIN32)
    set_target_properties(o.scofo_tilde PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/README.md")
pd_add_datafile(o.scofo_tilde "${CMAKE_CURRENT_SOURCE_DIR}/o.scofo~-help.pd")
pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.wav")
pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt")

if(WIN32 AND MINGW)
    find_package(Threads REQUIRED)

    target_compile_definitions(o.scofo_tilde PRIVATE _POSIX_SEM_VALUE_MAX=32767)
    target_compile_options(o.scofo_tilde PRIVATE -pthread)
    target_link_options(o.scofo_tilde PRIVATE
        -static-libgcc
        -static-libstdc++
        -pthread
        -Wl,--enable-auto-import
        -Wl,--enable-runtime-pseudo-reloc
        -Wl,--no-undefined
    )
    target_link_libraries(o.scofo_tilde PRIVATE
        Threads::Threads
        -Wl,--push-state,-Bstatic
        -lwinpthread
        -Wl,--pop-state
    )
endif()
