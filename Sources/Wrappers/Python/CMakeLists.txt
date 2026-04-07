# Use Python stable ABI (abi3) for the extension module.
# Development.Module is enough for extension builds and avoids requiring embed libs.
find_package(Python REQUIRED COMPONENTS Interpreter Development.SABIModule)

# Prefer a pip-installed nanobind when nanobind_DIR isn't set explicitly.
if(NOT nanobind_DIR)
    execute_process(COMMAND ${Python_EXECUTABLE} -m nanobind --cmake_dir
                    OUTPUT_VARIABLE NANOBIND_CMAKE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(NANOBIND_CMAKE_DIR)
        set(nanobind_DIR "${NANOBIND_CMAKE_DIR}" CACHE PATH "nanobind CMake directory")
    endif()
endif()
find_package(nanobind CONFIG REQUIRED)

# Build python module using stable ABI.
nanobind_add_module(_OpenScofo STABLE_ABI PyOpenScofo.cpp)
target_compile_options(_OpenScofo PRIVATE -Wno-shadow)
target_link_libraries(_OpenScofo PRIVATE OpenScofo Python::SABIModule)

target_link_libraries(_OpenScofo PRIVATE OpenScofo)
target_include_directories(_OpenScofo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../")

# collect python sources
file(GLOB OSCOFO_PY_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo/*.py)
install(TARGETS _OpenScofo LIBRARY DESTINATION OpenScofo)
install(FILES ${OSCOFO_PY_MODULE} DESTINATION OpenScofo)
