find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 REQUIRED)

# build python module
pybind11_add_module(_OpenScofo PyOpenScofo.cpp)

target_link_libraries(_OpenScofo PRIVATE OpenScofo ${PYTHON_LIBRARIES})
target_include_directories(_OpenScofo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../")

# directories for generated files
set(STUB_DIR ${CMAKE_CURRENT_BINARY_DIR}/stubs)
set(TYPED_FILE ${CMAKE_CURRENT_BINARY_DIR}/py.typed)

# generate stubs and py.typed after build
add_custom_command(
    TARGET _OpenScofo POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${STUB_DIR}
    COMMAND ${Python_EXECUTABLE} -m pybind11_stubgen _OpenScofo -o ${STUB_DIR} --ignore-all-errors
    COMMAND ${CMAKE_COMMAND} -E touch ${TYPED_FILE} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# collect python sources
file(GLOB OSCOFO_PY_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo/*.py)

# install compiled module
install(TARGETS _OpenScofo LIBRARY DESTINATION OpenScofo)

# install python files
install(FILES ${OSCOFO_PY_MODULE} DESTINATION OpenScofo)

# install generated stub
install(FILES ${STUB_DIR}/_OpenScofo.pyi DESTINATION OpenScofo)

# install typing marker (PEP 561)
install(FILES ${TYPED_FILE} DESTINATION OpenScofo)
