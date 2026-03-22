set(BUILD_UTILITIES OFF CACHE BOOL "Whether to build utility programs (e.g., csound, csound6, etc.)" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "Whether to build tests" FORCE)

cpmaddpackage(URI "gh:csound/csound#6.18.1" NAME csound DOWNLOAD_ONLY YES)

if(csound_ADDED)
    include_directories(${csound_SOURCE_DIR}/include)
    configure_file(${csound_SOURCE_DIR}/include/float-version.h.in ${CMAKE_CURRENT_BINARY_DIR}/float-version.h @ONLY)
    configure_file(${csound_SOURCE_DIR}/include/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h @ONLY)
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
endif()

file(GLOB_RECURSE OPENSCOFO_CSOUND_SOURCES "*.cpp" "*.h" "*.hpp")
add_library(OScofoCSound MODULE ${OPENSCOFO_CSOUND_SOURCES})
target_link_libraries(OScofoCSound PRIVATE OpenScofo)

if(APPLE)
    set_target_properties(OScofoCSound PROPERTIES PREFIX "" SUFFIX ".dylib")
elseif(WIN32)
    set_target_properties(OScofoCSound PROPERTIES PREFIX "" SUFFIX ".dll")
else()
    set_target_properties(OScofoCSound PROPERTIES PREFIX "" SUFFIX ".so")
endif()

if(UNIX AND NOT APPLE AND NOT EMSCRIPTEN)
    target_link_options(OScofoCSound PRIVATE -Wl,-allow-multiple-definition)
    target_compile_options(OScofoCSound PRIVATE -flto=auto)
    target_link_options(OScofoCSound PRIVATE -flto=auto)
endif()

