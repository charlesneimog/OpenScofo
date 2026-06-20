cpmaddpackage("gh:vamp-plugins/vamp-plugin-sdk#master")

add_library(VampOpenScofo SHARED "${CMAKE_CURRENT_SOURCE_DIR}/VampOpenScofo.cpp")
target_include_directories(VampOpenScofo PRIVATE ${vamp-plugin-sdk_SOURCE_DIR})
target_link_libraries(VampOpenScofo PRIVATE vamp-plugin-sdk::vamp-sdk OpenScofo)
if(WIN32)
    set_target_properties(VampOpenScofo PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
if(UNIX AND NOT WIN32)
    target_compile_options(VampOpenScofo PRIVATE -fPIC)
endif()

if(MINGW)
    target_link_options(VampOpenScofo PRIVATE "-Wl,-Bsymbolic"
                        "-Wl,--retain-symbols-file=${CMAKE_CURRENT_SOURCE_DIR}/vamp-plugin.list")
elseif(APPLE)
    target_link_options(
        VampOpenScofo PRIVATE
        "-Wl,-exported_symbols_list,${CMAKE_CURRENT_SOURCE_DIR}/vamp-plugin-macos.list")
elseif(UNIX)
    target_link_options(VampOpenScofo PRIVATE
                        "-Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/vamp-plugin.map")
endif()

set_target_properties(VampOpenScofo PROPERTIES PREFIX "" CXX_VISIBILITY_PRESET default
                                            VISIBILITY_INLINES_HIDDEN OFF OUTPUT_NAME "OpenScofo")

if(MINGW OR WIN32)
    set_target_properties(VampOpenScofo PROPERTIES SUFFIX ".dll")
elseif(APPLE)
    set_target_properties(VampOpenScofo PROPERTIES SUFFIX ".dylib")
else()
    set_target_properties(VampOpenScofo PROPERTIES SUFFIX ".so")
endif()
