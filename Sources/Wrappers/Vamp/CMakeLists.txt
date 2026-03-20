cpmaddpackage("gh:vamp-plugins/vamp-plugin-sdk#master")

add_library(VampOScofo SHARED "${CMAKE_CURRENT_SOURCE_DIR}/VampOpenScofo.cpp")
target_include_directories(VampOScofo PRIVATE ${vamp-plugin-sdk_SOURCE_DIR})
target_link_libraries(VampOScofo PRIVATE vamp-plugin-sdk::vamp-sdk OpenScofo)
if(UNIX AND NOT WIN32)
	target_compile_options(VampOScofo PRIVATE -fPIC)
endif()

if(MINGW)
	target_link_options(VampOScofo PRIVATE
		"-Wl,-Bsymbolic"
		"-Wl,--retain-symbols-file=${CMAKE_CURRENT_SOURCE_DIR}/vamp-plugin.list"
	)
elseif(APPLE)
	target_link_options(VampOScofo PRIVATE
		"-Wl,-exported_symbols_list,${CMAKE_CURRENT_SOURCE_DIR}/vamp-plugin-macos.list"
	)
elseif(UNIX)
    target_link_options(VampOScofo PRIVATE
        "-Wl,--version-script=${CMAKE_CURRENT_SOURCE_DIR}/vamp-plugin.map"
    )
endif()

set_target_properties(VampOScofo PROPERTIES
	PREFIX ""
	CXX_VISIBILITY_PRESET default
	VISIBILITY_INLINES_HIDDEN OFF
)

if(MINGW OR WIN32)
    set_target_properties(VampOScofo PROPERTIES SUFFIX ".dll")
elseif(APPLE)
    set_target_properties(VampOScofo PROPERTIES SUFFIX ".dylib")
else()
    set_target_properties(VampOScofo PROPERTIES SUFFIX ".so")
endif()
