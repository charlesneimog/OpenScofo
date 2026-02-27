find_package(
    Python
    COMPONENTS Interpreter Development.Module
    REQUIRED)

find_package(pybind11 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${pybind11_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Sources)

pybind11_add_module(_OpenScofo PyOpenScofo.cpp)
target_link_libraries(_OpenScofo PRIVATE OpenScofo ${PYTHON_LIBRARIES})
target_include_directories(_OpenScofo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../")
# set_target_properties(_OpenScofo PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES
# CXX_EXTENSIONS NO)

file(GLOB OSCOFO_PY_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo/*.py)

# if(NOT APPLE AND NOT MSVC) target_compile_options(_OpenScofo PRIVATE -static-libstdc++
# -Wl,--exclude-libs,ALL) target_link_options(_OpenScofo PRIVATE -static-libstdc++
# -Wl,--exclude-libs,ALL) endif()

install(TARGETS _OpenScofo LIBRARY DESTINATION OpenScofo)
install(FILES "${OSCOFO_PY_MODULE}" DESTINATION OpenScofo)
