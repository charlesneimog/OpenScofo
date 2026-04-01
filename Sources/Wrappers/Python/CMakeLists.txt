# Need Development.Module so CMake defines Python_add_library for extension targets.
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 REQUIRED)

# Build python module
pybind11_add_module(_OpenScofo PyOpenScofo.cpp)

target_link_libraries(_OpenScofo PRIVATE OpenScofo)
target_include_directories(_OpenScofo PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../")

# directories for generated files
message(STATUS "Python executable: ${Python_EXECUTABLE}")

set(STUB_DIR ${CMAKE_CURRENT_BINARY_DIR}/stubs)
set(TYPED_FILE ${CMAKE_CURRENT_BINARY_DIR}/py.typed)
get_filename_component(MINGW_BIN_DIR ${CMAKE_CXX_COMPILER} DIRECTORY)

if(WIN32)
    # Keep ';' as-is for Windows PATH. Quoting at call site keeps this as one argument.
    set(STUBGEN_PATH_ENV "PATH=${CMAKE_CURRENT_BINARY_DIR};${CMAKE_BINARY_DIR};${MINGW_BIN_DIR};$ENV{PATH}")
else()
    set(STUBGEN_PATH_ENV "PATH=$ENV{PATH}")
endif()

if (NOT WIN32)
    # generate stubs and py.typed after build
    add_custom_command(
        TARGET _OpenScofo
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${STUB_DIR}
            COMMAND ${CMAKE_COMMAND} -E env "${STUBGEN_PATH_ENV}" "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}" ${Python_EXECUTABLE} -m
                pybind11_stubgen _OpenScofo -o ${STUB_DIR} --ignore-all-errors
        COMMAND ${CMAKE_COMMAND} -E touch ${TYPED_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

# collect python sources
file(GLOB OSCOFO_PY_MODULE ${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo/*.py)
install(TARGETS _OpenScofo LIBRARY DESTINATION OpenScofo)
install(FILES ${OSCOFO_PY_MODULE} DESTINATION OpenScofo)

if (NOT WIN32)
    install(FILES ${STUB_DIR}/_OpenScofo.pyi DESTINATION OpenScofo)
    install(FILES ${TYPED_FILE} DESTINATION OpenScofo)
endif()

# add_custom_target( run_descriptors_validation ALL COMMAND ${CMAKE_COMMAND} -E env
# PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR} ${Python_EXECUTABLE}
# ${CMAKE_SOURCE_DIR}/Tests/librosa-essentia-validation.py DEPENDS _OpenScofo)
