project(o.scofo~ C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PDCMAKE_FILE ${CMAKE_BINARY_DIR}/pd.cmake)
set(PDCMAKE_VERSION "v0.2.12")
if(NOT EXISTS "${PDCMAKE_FILE}")
    file(DOWNLOAD
         https://raw.githubusercontent.com/charlesneimog/pd.cmake/refs/heads/flags-test/pd.cmake
         ${PDCMAKE_FILE})
endif()
include(${PDCMAKE_FILE})

# Keep Pd installs prefix-relative so scikit-build wheel installs do not try to write to /usr/local.
set(PD_LIB_DIR "lib/pd-externals" CACHE PATH "Path where lib will be installed" FORCE)

file(GLOB_RECURSE OPENSCOFO_PDOBJ_SOURCES "*.cpp" "*.h")
pd_add_external(o.scofo~ "${OPENSCOFO_PDOBJ_SOURCES}" LINK_LIBRARIES OpenScofo)

pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/README.md")
pd_add_datafile(o.scofo_tilde "${CMAKE_CURRENT_SOURCE_DIR}/o.scofo~-help.pd")
pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.wav")
pd_add_datafile(o.scofo_tilde "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt")

target_compile_definitions(o.scofo_tilde PRIVATE _POSIX_SEM_VALUE_MAX=32767)
