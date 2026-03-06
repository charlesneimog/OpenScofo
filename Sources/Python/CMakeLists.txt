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

file(GLOB OSCOFO_PY_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo/*.py)

install(TARGETS _OpenScofo LIBRARY DESTINATION OpenScofo)
install(FILES "${OSCOFO_PY_MODULE}" DESTINATION OpenScofo)
