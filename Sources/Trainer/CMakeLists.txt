set(CATBOOST_VERSION
    "v1.2.10"
    CACHE STRING "CatBoost release used by the native trainer")
set(OPENSCOFO_CATBOOST_SOURCE_DIR
    ""
    CACHE PATH "Use an existing CatBoost source checkout instead of downloading it")
set(OPENSCOFO_CATBOOST_TOOLS_DIR
    ""
    CACHE PATH "Directory containing CatBoost host build tools (ragel, yasm, and swig)")
set(OPENSCOFO_CONAN_EXECUTABLE
    ""
    CACHE FILEPATH "Conan executable used to prepare CatBoost host build tools")

if(NOT WIN32 AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "CatBoost's native C++ training engine must be built with Clang on this platform. "
                        "Configure OpenScofo with -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++.")
endif()

find_package(Threads REQUIRED)

if(NOT TARGET signalsmith-stretch)
    FetchContent_Declare(
        signalsmith-stretch
        GIT_REPOSITORY https://github.com/Signalsmith-Audio/signalsmith-stretch.git
        GIT_TAG 1.4.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(signalsmith-stretch)
endif()

# The public CatBoost C API only applies existing models.  Training, early-stop callbacks, native model ownership, and
# ONNX export are provided by these C++ targets from the CatBoost source distribution.
set(CATBOOST_COMPONENTS
    libs
    CACHE STRING "CatBoost components required by OpenScofo" FORCE)
set(HAVE_CUDA
    OFF
    CACHE BOOL "OpenScofo native training uses CatBoost CPU training" FORCE)

if(OPENSCOFO_CATBOOST_SOURCE_DIR)
    set(openscofo_catboost_SOURCE_DIR "${OPENSCOFO_CATBOOST_SOURCE_DIR}")
    set(openscofo_catboost_BINARY_DIR "${CMAKE_BINARY_DIR}/openscofo-catboost-build")
else()
    FetchContent_Declare(
        openscofo_catboost
        GIT_REPOSITORY https://github.com/catboost/catboost.git
        GIT_TAG "${CATBOOST_VERSION}"
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
        EXCLUDE_FROM_ALL
        # FetchContent obtains the sources but CatBoost is added explicitly below so it can remain excluded from the
        # default non-training build.
        SOURCE_SUBDIR
        __openscofo_source_only__)
    FetchContent_MakeAvailable(openscofo_catboost)
endif()

if(NOT EXISTS "${openscofo_catboost_SOURCE_DIR}/catboost/libs/train_lib/train_model.h")
    message(FATAL_ERROR "OPENSCOFO_CATBOOST_SOURCE_DIR does not contain a CatBoost source tree")
endif()

# CatBoost 1.2.10 is compatible with the OpenSSL range represented by its official 3.0.15 Conan requirement.  Prefer an
# existing compatible install, but do not accidentally compile CatBoost against OpenSSL 4.1's changed API.
find_package(OpenSSL 3.0...<4.1 QUIET)

# CatBoost's generated CMake build expects its code-generation tools in <build>/bin.  Its official standalone build puts
# them there from the tool_requires in CatBoost's conanfile.py.  Reproduce only that build-time preparation here; Conan
# and Python are never linked to or invoked by the OpenScofo training library at runtime.
file(MAKE_DIRECTORY "${openscofo_catboost_BINARY_DIR}/bin")

function(openscofo_copy_catboost_tool tool)
    unset(source)
    if(OPENSCOFO_CATBOOST_TOOLS_DIR)
        find_program(
            source
            NAMES "${tool}"
            HINTS "${OPENSCOFO_CATBOOST_TOOLS_DIR}" "${OPENSCOFO_CATBOOST_TOOLS_DIR}/bin"
            NO_DEFAULT_PATH NO_CACHE)
    else()
        find_program(source NAMES "${tool}" NO_CACHE)
    endif()
    if(source)
        file(COPY "${source}" DESTINATION "${openscofo_catboost_BINARY_DIR}/bin")
    endif()
endfunction()

openscofo_copy_catboost_tool(ragel)
openscofo_copy_catboost_tool(yasm)
openscofo_copy_catboost_tool(swig)
if(OPENSCOFO_CATBOOST_TOOLS_DIR AND EXISTS "${OPENSCOFO_CATBOOST_TOOLS_DIR}/swiglib")
    file(COPY "${OPENSCOFO_CATBOOST_TOOLS_DIR}/swiglib" DESTINATION "${openscofo_catboost_BINARY_DIR}/bin")
elseif(OPENSCOFO_CATBOOST_TOOLS_DIR AND EXISTS "${OPENSCOFO_CATBOOST_TOOLS_DIR}/bin/swiglib")
    file(COPY "${OPENSCOFO_CATBOOST_TOOLS_DIR}/bin/swiglib" DESTINATION "${openscofo_catboost_BINARY_DIR}/bin")
endif()

set(openscofo_catboost_ragel "${openscofo_catboost_BINARY_DIR}/bin/ragel${CMAKE_EXECUTABLE_SUFFIX}")
set(openscofo_catboost_yasm "${openscofo_catboost_BINARY_DIR}/bin/yasm${CMAKE_EXECUTABLE_SUFFIX}")
set(openscofo_catboost_conan_generators "${openscofo_catboost_BINARY_DIR}/conan/build/${CMAKE_BUILD_TYPE}/generators")
if(NOT CMAKE_BUILD_TYPE)
    set(openscofo_catboost_conan_generators "${openscofo_catboost_BINARY_DIR}/conan/build/Release/generators")
endif()
set(openscofo_catboost_need_conan FALSE)
if(NOT EXISTS "${openscofo_catboost_ragel}" OR NOT EXISTS "${openscofo_catboost_yasm}")
    set(openscofo_catboost_need_conan TRUE)
endif()
if(NOT OpenSSL_FOUND AND NOT EXISTS "${openscofo_catboost_conan_generators}/OpenSSLConfig.cmake")
    set(openscofo_catboost_need_conan TRUE)
endif()

if(openscofo_catboost_need_conan)
    if(NOT OPENSCOFO_CONAN_EXECUTABLE)
        find_program(openscofo_conan NAMES conan)
        if(openscofo_conan)
            set(OPENSCOFO_CONAN_EXECUTABLE
                "${openscofo_conan}"
                CACHE FILEPATH "Conan executable used to prepare CatBoost host build tools" FORCE)
        endif()
    endif()
    if(NOT OPENSCOFO_CONAN_EXECUTABLE)
        message(
            FATAL_ERROR
                "CatBoost requires Ragel and Yasm while compiling its native training engine. "
                "Install Conan 2, set OPENSCOFO_CONAN_EXECUTABLE, or point "
                "OPENSCOFO_CATBOOST_TOOLS_DIR at a CatBoost tools/bin directory.")
    endif()

    execute_process(
        COMMAND "${OPENSCOFO_CONAN_EXECUTABLE}" profile path default
        RESULT_VARIABLE openscofo_conan_profile_result
        OUTPUT_QUIET ERROR_QUIET)
    if(NOT openscofo_conan_profile_result EQUAL 0)
        execute_process(COMMAND "${OPENSCOFO_CONAN_EXECUTABLE}" profile detect
                        RESULT_VARIABLE openscofo_conan_profile_result COMMAND_ERROR_IS_FATAL ANY)
    endif()

    set(openscofo_catboost_conan_build_type "${CMAKE_BUILD_TYPE}")
    if(NOT openscofo_catboost_conan_build_type)
        set(openscofo_catboost_conan_build_type Release)
    endif()
    message(STATUS "Preparing CatBoost host build tools with Conan")
    execute_process(
        COMMAND
            "${OPENSCOFO_CONAN_EXECUTABLE}" install "${openscofo_catboost_SOURCE_DIR}"
            "-of=${openscofo_catboost_BINARY_DIR}/conan" "-s" "build_type=${openscofo_catboost_conan_build_type}"
            "--build=missing"
        RESULT_VARIABLE openscofo_conan_install_result COMMAND_ERROR_IS_FATAL ANY)
endif()

if(NOT OpenSSL_FOUND)
    if(NOT EXISTS "${openscofo_catboost_conan_generators}/OpenSSLConfig.cmake")
        message(FATAL_ERROR "Conan did not generate CatBoost's pinned OpenSSL package configuration")
    endif()
    set(openscofo_saved_prefix_path "${CMAKE_PREFIX_PATH}")
    list(PREPEND CMAKE_PREFIX_PATH "${openscofo_catboost_conan_generators}")
    find_package(
        OpenSSL
        3.0.15
        EXACT
        CONFIG
        REQUIRED
        PATHS
        "${openscofo_catboost_conan_generators}"
        NO_DEFAULT_PATH)
    set(CMAKE_PREFIX_PATH "${openscofo_saved_prefix_path}")
endif()

# CatBoost links the lowercase aggregate target produced by its Conan build. Map a compatible system OpenSSL
# installation to that spelling when Conan was not needed.
if(NOT TARGET openssl::openssl)
    add_library(openssl::openssl INTERFACE IMPORTED GLOBAL)
    set_property(TARGET openssl::openssl PROPERTY INTERFACE_LINK_LIBRARIES "OpenSSL::SSL;OpenSSL::Crypto")
endif()

if(NOT EXISTS "${openscofo_catboost_ragel}" OR NOT EXISTS "${openscofo_catboost_yasm}")
    message(FATAL_ERROR "CatBoost host build tools were not installed in ${openscofo_catboost_BINARY_DIR}/bin")
endif()

# CatBoost deliberately links its bundled C++ runtime with -nodefaultlibs. Its official Clang toolchain supplies libc
# and libm explicitly; an embedded build must do the same because the parent toolchain has already initialized.
set(openscofo_saved_c_standard_libraries "${CMAKE_C_STANDARD_LIBRARIES}")
set(openscofo_saved_cxx_standard_libraries "${CMAKE_CXX_STANDARD_LIBRARIES}")
set(openscofo_saved_exe_linker_flags "${CMAKE_EXE_LINKER_FLAGS}")
set(openscofo_saved_shared_linker_flags "${CMAKE_SHARED_LINKER_FLAGS}")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    find_program(openscofo_lld NAMES ld.lld REQUIRED)
    set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES} -lc -lm")
    set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} -lc -lm")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fuse-ld=lld")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=lld")
endif()
add_subdirectory("${openscofo_catboost_SOURCE_DIR}" "${openscofo_catboost_BINARY_DIR}" EXCLUDE_FROM_ALL)
if(NOT TARGET catboost-libs-train_lib)
    message(FATAL_ERROR "CatBoost's native C++ training target is unavailable on this platform")
endif()

# CatBoost's generated build uses its own libc++ implementation. Keep that runtime behind a C ABI so OpenScofo and
# clients of the public Trainer.hpp API continue to use the parent project's standard-library ABI. The bridge itself is
# static so consumers such as the Pure Data external do not need a companion libopenscofo_catboost_training shared file.
add_library(openscofo_catboost_bridge STATIC CatBoostBridge.cpp CatBoostVersion.c)
target_compile_definitions(
    openscofo_catboost_bridge
    PUBLIC OPENSCOFO_CATBOOST_STATIC
    PRIVATE "CATBOOST_VERSION=\"${CATBOOST_VERSION}\"")
target_include_directories(
    openscofo_catboost_bridge
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" SYSTEM
    PRIVATE "${openscofo_catboost_SOURCE_DIR}" "${openscofo_catboost_BINARY_DIR}")
target_link_libraries(
    openscofo_catboost_bridge
    PRIVATE catboost-libs-train_lib
            libs-model-model_export
            catboost-libs-model
            catboost-libs-data
            catboost-libs-helpers
            cpp-threading-local_executor
            library-cpp-json
            yutil
            contrib-libs-cxxsupp
            build-cow-on)
if(COMMAND target_allocator)
    target_allocator(openscofo_catboost_bridge cpp-malloc-tcmalloc)
endif()
set_target_properties(
    openscofo_catboost_bridge
    PROPERTIES OUTPUT_NAME openscofo_catboost_training
               POSITION_INDEPENDENT_CODE ON
               C_VISIBILITY_PRESET hidden
               CXX_VISIBILITY_PRESET hidden
               VISIBILITY_INLINES_HIDDEN ON)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_options(openscofo_catboost_bridge INTERFACE -fuse-ld=lld)
    target_link_libraries(openscofo_catboost_bridge PRIVATE c m dl rt)
endif()

set(CMAKE_C_STANDARD_LIBRARIES "${openscofo_saved_c_standard_libraries}")
set(CMAKE_CXX_STANDARD_LIBRARIES "${openscofo_saved_cxx_standard_libraries}")
set(CMAKE_EXE_LINKER_FLAGS "${openscofo_saved_exe_linker_flags}")
set(CMAKE_SHARED_LINKER_FLAGS "${openscofo_saved_shared_linker_flags}")

add_library(openscofo_training STATIC Audio.cpp CatBoostBackend.cpp Trainer.cpp)
add_library(OpenScofo::Trainer ALIAS openscofo_training)
target_include_directories(openscofo_training PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(
    openscofo_training
    PUBLIC OpenScofo
    PRIVATE openscofo_catboost_bridge AudioFile PFFFT::PFFFT signalsmith-stretch Threads::Threads)
set_target_properties(openscofo_training PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(OPENSCOFO_BUILD_PD_OBJECT)
    pd_add_external(openscofo.train "${CMAKE_CURRENT_SOURCE_DIR}/PureData/openscofo.train.cpp" TARGET openscofo_train
                    LINK_LIBRARIES openscofo_training)
    pd_add_datafile(openscofo_train "${CMAKE_CURRENT_SOURCE_DIR}/PureData/openscofo.train-help.pd")
    if(LINUX)
        target_link_options(openscofo_train PRIVATE -fuse-ld=lld -Wl,--exclude-libs,ALL -Wl,-z,notext)
        target_compile_options(openscofo_train PRIVATE -Wno-cast-function-type -Wno-cast-function-type-mismatch)
    elseif(APPLE)
        target_compile_options(openscofo_train PRIVATE -fvisibility=hidden -Wno-cast-function-type-mismatch)
    endif()
    set_target_properties(openscofo_train PROPERTIES SKIP_BUILD_RPATH ON INSTALL_RPATH "")
endif()

if(OPENSCOFO_BUILD_TESTS)
    add_subdirectory(tests)
endif()
