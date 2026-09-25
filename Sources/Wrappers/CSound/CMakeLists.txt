set(BUILD_UTILITIES OFF
    CACHE BOOL "Whether to build utility programs (e.g., csound, csound6, etc.)" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "Whether to build tests" FORCE)

FetchContent_Declare(csound GIT_REPOSITORY https://github.com/csound/csound.git GIT_TAG 6.18.1
                     SOURCE_SUBDIR __do_not_add__)
FetchContent_MakeAvailable(csound)

configure_file(${csound_SOURCE_DIR}/include/float-version.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/float-version.h @ONLY)
configure_file(${csound_SOURCE_DIR}/include/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h
               @ONLY)

file(GLOB_RECURSE OPENSCOFO_CSOUND_SOURCES "*.cpp" "*.h" "*.hpp")
add_library(OpenScofoCSound MODULE ${OPENSCOFO_CSOUND_SOURCES})
target_link_libraries(OpenScofoCSound PRIVATE OpenScofo)
target_compile_definitions(OpenScofoCSound PUBLIC USE_DOUBLE)
target_include_directories(OpenScofoCSound SYSTEM PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(OpenScofoCSound SYSTEM PRIVATE ${csound_SOURCE_DIR}/include)

if(APPLE)
    set_target_properties(OpenScofoCSound PROPERTIES PREFIX "" SUFFIX ".dylib" OUTPUT_NAME
                                                                               "OpenScofo")
    target_link_options(OpenScofoCSound PRIVATE LINKER:-undefined,dynamic_lookup)
elseif(WIN32)
    set_target_properties(OpenScofoCSound PROPERTIES PREFIX "" SUFFIX ".dll" OUTPUT_NAME
                                                                             "OpenScofo")
else()
    set_target_properties(OpenScofoCSound PROPERTIES PREFIX "" SUFFIX ".so" OUTPUT_NAME "OpenScofo")
endif()

if(UNIX AND NOT APPLE AND NOT EMSCRIPTEN)
    target_link_options(OpenScofoCSound PRIVATE -Wl,-allow-multiple-definition)
    target_compile_options(OpenScofoCSound PRIVATE -flto=auto)
    target_link_options(OpenScofoCSound PRIVATE -flto=auto)
endif()
