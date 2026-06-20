find_package(Python REQUIRED COMPONENTS Interpreter Development Development.SABIModule)

if(NOT nanobind_DIR)
    execute_process(COMMAND ${Python_EXECUTABLE} -m nanobind --cmake_dir
                    OUTPUT_VARIABLE NANOBIND_CMAKE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(NANOBIND_CMAKE_DIR)
        set(nanobind_DIR "${NANOBIND_CMAKE_DIR}")
    else()
        message(FATAL_ERROR "Nanobind not found")
    endif()
endif()

find_package(nanobind CONFIG REQUIRED)

# Build python module using stable ABI.
nanobind_add_module(_OpenScofo NB_SUPPRESS_WARNINGS STABLE_ABI PyOpenScofo.cpp)
target_compile_options(_OpenScofo PRIVATE -Wno-shadow)
target_link_libraries(_OpenScofo PRIVATE OpenScofo)

target_include_directories(_OpenScofo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../")
set_target_properties(_OpenScofo PROPERTIES SUFFIX ".${Python_SOSABI}${CMAKE_SHARED_MODULE_SUFFIX}")

# 1. Hide all symbols natively compiled into the wrapper
set_target_properties(_OpenScofo PROPERTIES CXX_VISIBILITY_PRESET hidden C_VISIBILITY_PRESET hidden
                                            VISIBILITY_INLINES_HIDDEN YES)

if(UNIX AND NOT APPLE)
    target_link_options(_OpenScofo PRIVATE "-Wl,--exclude-libs,ALL")
endif()

# collect python sources
file(GLOB OPENSCOFO_PY_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo/*.py)
install(TARGETS _OpenScofo LIBRARY DESTINATION OpenScofo)
install(FILES ${OPENSCOFO_PY_MODULE} DESTINATION OpenScofo)
