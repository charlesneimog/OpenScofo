include_guard(GLOBAL)

# Configure CatBoost's native CPU training engine and expose it through OpenScofo::CatBoost.
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
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 20)
    message(FATAL_ERROR "CatBoost's bundled libc++ requires Clang 20 or newer on Linux")
endif()

find_package(OpenSSL 3.0...<4.1 QUIET)

# The public CatBoost C API only applies existing models. Training, callbacks, model ownership, and ONNX export require
# the native C++ targets from the CatBoost source distribution.
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
    set(openscofo_catboost_fetch_options)
    if(WIN32)
        list(APPEND openscofo_catboost_fetch_options SOURCE_DIR "${CMAKE_BINARY_DIR}/cb-src")
        add_compile_definitions(WIN32_LEAN_AND_MEAN _USE_MATH_DEFINES NOMINMAX)
    endif()

    FetchContent_Declare(
        openscofo_catboost
        GIT_REPOSITORY https://github.com/catboost/catboost.git
        GIT_TAG "${CATBOOST_VERSION}"
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
        EXCLUDE_FROM_ALL SOURCE_SUBDIR __openscofo_source_only__ ${openscofo_catboost_fetch_options})
    FetchContent_MakeAvailable(openscofo_catboost)
endif()

if(NOT EXISTS "${openscofo_catboost_SOURCE_DIR}/catboost/libs/train_lib/train_model.h")
    message(FATAL_ERROR "OPENSCOFO_CATBOOST_SOURCE_DIR does not contain a CatBoost source tree")
endif()

file(MAKE_DIRECTORY "${openscofo_catboost_BINARY_DIR}/bin")
cmake_path(GET openscofo_catboost_BINARY_DIR PARENT_PATH openscofo_catboost_parent_binary_dir)
set(openscofo_catboost_conan_tools_dirs "${CMAKE_BINARY_DIR}/bin" "${openscofo_catboost_parent_binary_dir}/bin")

function(openscofo_copy_catboost_tool tool)
    unset(source)
    if(OPENSCOFO_CATBOOST_TOOLS_DIR)
        find_program(
            source
            NAMES "${tool}"
            HINTS "${OPENSCOFO_CATBOOST_TOOLS_DIR}" "${OPENSCOFO_CATBOOST_TOOLS_DIR}/bin"
                  ${openscofo_catboost_conan_tools_dirs}
            NO_DEFAULT_PATH NO_CACHE)
    else()
        find_program(
            source
            NAMES "${tool}"
            HINTS ${openscofo_catboost_conan_tools_dirs} NO_CACHE)
    endif()
    if(source)
        # Homebrew exposes these tools through relative symlinks. Resolve them before copying into CatBoost's tree.
        file(REAL_PATH "${source}" resolved_source)
        file(
            COPY "${resolved_source}"
            DESTINATION "${openscofo_catboost_BINARY_DIR}/bin"
            USE_SOURCE_PERMISSIONS)
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
set(openscofo_catboost_conan_build_type "${CMAKE_BUILD_TYPE}")
if(NOT openscofo_catboost_conan_build_type)
    set(openscofo_catboost_conan_build_type Release)
endif()

function(openscofo_find_catboost_conan_generators output_variable)
    foreach(candidate "${openscofo_catboost_BINARY_DIR}/conan/build/${openscofo_catboost_conan_build_type}/generators"
                      "${openscofo_catboost_BINARY_DIR}/conan/build/generators")
        if(EXISTS "${candidate}/OpenSSLConfig.cmake")
            set(${output_variable}
                "${candidate}"
                PARENT_SCOPE)
            return()
        endif()
    endforeach()
    set(${output_variable}
        ""
        PARENT_SCOPE)
endfunction()

openscofo_find_catboost_conan_generators(openscofo_catboost_conan_generators)
set(openscofo_catboost_need_conan FALSE)
if(NOT EXISTS "${openscofo_catboost_ragel}" OR NOT EXISTS "${openscofo_catboost_yasm}")
    set(openscofo_catboost_need_conan TRUE)
endif()
if(NOT OpenSSL_FOUND AND NOT openscofo_catboost_conan_generators)
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
        execute_process(COMMAND "${OPENSCOFO_CONAN_EXECUTABLE}" profile detect COMMAND_ERROR_IS_FATAL ANY)
    endif()

    message(STATUS "Preparing CatBoost host build tools with Conan")
    execute_process(
        COMMAND
            "${OPENSCOFO_CONAN_EXECUTABLE}" install "${openscofo_catboost_SOURCE_DIR}"
            "-of=${openscofo_catboost_BINARY_DIR}/conan" "-s" "build_type=${openscofo_catboost_conan_build_type}"
            "--build=missing"
        COMMAND_ERROR_IS_FATAL ANY)

    # CatBoost may copy generated host tools beside the dependency build rather than into its own binary directory.
    openscofo_copy_catboost_tool(ragel)
    openscofo_copy_catboost_tool(yasm)
    openscofo_copy_catboost_tool(swig)
    foreach(openscofo_catboost_conan_tools_dir IN LISTS openscofo_catboost_conan_tools_dirs)
        if(EXISTS "${openscofo_catboost_conan_tools_dir}/swiglib")
            file(COPY "${openscofo_catboost_conan_tools_dir}/swiglib"
                 DESTINATION "${openscofo_catboost_BINARY_DIR}/bin")
            break()
        endif()
    endforeach()
endif()

if(NOT OpenSSL_FOUND)
    openscofo_find_catboost_conan_generators(openscofo_catboost_conan_generators)
    if(NOT openscofo_catboost_conan_generators)
        message(FATAL_ERROR "Conan did not generate CatBoost's pinned OpenSSL package configuration")
    endif()
    set(openscofo_saved_prefix_path "${CMAKE_PREFIX_PATH}")
    list(PREPEND CMAKE_PREFIX_PATH "${openscofo_catboost_conan_generators}")
    # Force the configuration selected above so repeated Debug/Release configurations cannot reuse stale targets.
    set(OpenSSL_DIR
        "${openscofo_catboost_conan_generators}"
        CACHE PATH "Directory containing CatBoost's pinned OpenSSL package configuration" FORCE)
    set(ZLIB_DIR
        "${openscofo_catboost_conan_generators}"
        CACHE PATH "Directory containing CatBoost's pinned OpenSSL dependency configuration" FORCE)
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

# CatBoost expects the lowercase aggregate target produced by its Conan build.
if(NOT TARGET openssl::openssl)
    add_library(openssl::openssl INTERFACE IMPORTED GLOBAL)
    set_property(TARGET openssl::openssl PROPERTY INTERFACE_LINK_LIBRARIES "OpenSSL::SSL;OpenSSL::Crypto")
endif()

if(NOT EXISTS "${openscofo_catboost_ragel}" OR NOT EXISTS "${openscofo_catboost_yasm}")
    message(FATAL_ERROR "CatBoost host build tools were not installed in ${openscofo_catboost_BINARY_DIR}/bin")
endif()

# CatBoost links its bundled C++ runtime with -nodefaultlibs, so supply the platform runtime while configuring it.
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
elseif(APPLE)
    set(CMAKE_C_STANDARD_LIBRARIES "${CMAKE_C_STANDARD_LIBRARIES} -lSystem")
    set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_STANDARD_LIBRARIES} -lSystem")
endif()

# CatBoost's generated subdirectories call find_package(OpenSSL) again. Keep them on the selected configuration.
set(openscofo_saved_find_package_prefer_config "${CMAKE_FIND_PACKAGE_PREFER_CONFIG}")
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
add_subdirectory("${openscofo_catboost_SOURCE_DIR}" "${openscofo_catboost_BINARY_DIR}" EXCLUDE_FROM_ALL)
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG "${openscofo_saved_find_package_prefer_config}")
if(NOT TARGET catboost-libs-train_lib)
    message(FATAL_ERROR "CatBoost's native C++ training target is unavailable on this platform")
endif()

if(APPLE AND TARGET contrib-libs-tbb)
    target_compile_definitions(contrib-libs-tbb PUBLIC _LIBCPP_ENABLE_CXX20_REMOVED_UNCAUGHT_EXCEPTION)
endif()

# CatBoost registers excluded upstream tests. Disable only those dependency tests so unfiltered CTest runs stay valid.
function(openscofo_disable_dependency_tests directory)
    get_property(
        dependency_tests
        DIRECTORY "${directory}"
        PROPERTY TESTS)
    if(dependency_tests)
        set_tests_properties(${dependency_tests} DIRECTORY "${directory}" PROPERTIES DISABLED TRUE)
    endif()
    get_property(
        dependency_subdirectories
        DIRECTORY "${directory}"
        PROPERTY SUBDIRECTORIES)
    foreach(dependency_subdirectory IN LISTS dependency_subdirectories)
        openscofo_disable_dependency_tests("${dependency_subdirectory}")
    endforeach()
endfunction()
openscofo_disable_dependency_tests("${openscofo_catboost_BINARY_DIR}")

if(APPLE AND TARGET libs-cxxsupp-builtins)
    # CatBoost includes the generic ELF SME ABI implementation on Darwin arm64, where it is unused and cannot link.
    get_target_property(openscofo_catboost_builtin_sources libs-cxxsupp-builtins SOURCES)
    list(FILTER openscofo_catboost_builtin_sources EXCLUDE REGEX "/aarch64/sme-abi(\\.S|-assert\\.c)$")
    set_property(TARGET libs-cxxsupp-builtins PROPERTY SOURCES "${openscofo_catboost_builtin_sources}")
endif()

# Isolate CatBoost's bundled libc++ behind the bridge's C ABI.
cmake_path(SET openscofo_trainer_source_dir NORMALIZE "${CMAKE_CURRENT_LIST_DIR}/../../Sources/Trainer")
add_library(
    openscofo_catboost_bridge STATIC "${openscofo_trainer_source_dir}/CatBoostBridge.cpp"
                                     "${openscofo_trainer_source_dir}/CatBoostVersion.c")
add_library(OpenScofo::CatBoost ALIAS openscofo_catboost_bridge)
target_compile_definitions(
    openscofo_catboost_bridge
    PUBLIC OPENSCOFO_CATBOOST_STATIC
    PRIVATE "CATBOOST_VERSION=\"${CATBOOST_VERSION}\"")
target_include_directories(
    openscofo_catboost_bridge
    PRIVATE "${openscofo_trainer_source_dir}" SYSTEM
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
if(COMMAND target_allocator AND TARGET cpp-malloc-tcmalloc)
    target_allocator(openscofo_catboost_bridge cpp-malloc-tcmalloc)
elseif(COMMAND target_allocator AND TARGET system_allocator)
    target_allocator(openscofo_catboost_bridge system_allocator)
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
