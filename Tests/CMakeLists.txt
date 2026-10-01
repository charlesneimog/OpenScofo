# Build GoogleTest with the target toolchain for universal macOS, Windows and WebAssembly. A native Linux installation
# can be reused when available.
if(NOT APPLE
   AND NOT WIN32
   AND NOT CMAKE_CROSSCOMPILING)
    find_package(GTest QUIET)

    if(NOT TARGET GTest::gtest_main)
        set(INSTALL_GTEST
            OFF
            CACHE BOOL "" FORCE)
        set(BUILD_GMOCK
            OFF
            CACHE BOOL "" FORCE)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.17.0)
        FetchContent_MakeAvailable(googletest)
    endif()

    FetchContent_Declare(
        AudioFile
        GIT_REPOSITORY https://github.com/adamstark/AudioFile.git
        GIT_TAG master)
    FetchContent_MakeAvailable(AudioFile)

    FetchContent_Declare(
        minimp3
        GIT_REPOSITORY https://github.com/lieff/minimp3.git
        GIT_TAG master)
    FetchContent_MakeAvailable(minimp3)

    function(openscofo_add_test Target Source Name)
        add_executable(${Target} "${CMAKE_CURRENT_SOURCE_DIR}/05-tests/${Source}")
        target_link_libraries(${Target} PRIVATE OpenScofo GTest::gtest_main)
        target_compile_definitions(${Target} PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/02-score")
        add_test(NAME ${Name} COMMAND ${Target})
        set_tests_properties(${Name} PROPERTIES LABELS "openscofo" TIMEOUT 120)
        set_property(GLOBAL APPEND PROPERTY OPENSCOFO_TEST_TARGETS ${Target})
        target_compile_definitions(${Target} PRIVATE BENCHMARK_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
        target_link_libraries(${Target} PRIVATE AudioFile)
        target_include_directories(${Target} SYSTEM PRIVATE "${minimp3_SOURCE_DIR}")
        target_include_directories(${Target} SYSTEM PRIVATE "${audiofile_SOURCE_DIR}")
    endfunction()

    openscofo_add_test(main_class test_openscofo.cpp OpenScofo)
    openscofo_add_test(score_tests test_score.cpp Score)
    openscofo_add_test(configuration_tests test_configuration.cpp Configuration)
    openscofo_add_test(trill_tests test_trill.cpp Trill)
    openscofo_add_test(microstate_tests test_microstates.cpp MicroStates)
    openscofo_add_test(math_tests test_math.cpp Math)
    openscofo_add_test(onnx_tests test_onnx.cpp Onnx)
    openscofo_add_test(inference_tests test_benchmark.cpp Inference)

    get_property(TestTargets GLOBAL PROPERTY OPENSCOFO_TEST_TARGETS)
    if(OPENSCOFO_RUN_TESTS_ON_BUILD)
        set(CheckAll ALL)
    endif()

    add_custom_target(
        check
        ${CheckAll}
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure --no-tests=error -L openscofo -C "$<CONFIG>"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
        DEPENDS ${TestTargets}
        USES_TERMINAL
        COMMENT "Running OpenScofo tests"
        VERBATIM)

    add_executable(openscofo_perf "${CMAKE_CURRENT_SOURCE_DIR}/03-performance/performance.cpp")
    target_link_libraries(openscofo_perf PRIVATE OpenScofo AudioFile)
    target_compile_options(openscofo_perf PRIVATE -O2 -g -fno-omit-frame-pointer -Wno-shadow -Wno-unused-variable)

    #
    add_executable(openscofo_perf_folder "${CMAKE_CURRENT_SOURCE_DIR}/03-performance/performance-folder.cpp")
    target_link_libraries(openscofo_perf_folder PRIVATE OpenScofo AudioFile)
    target_compile_options(openscofo_perf_folder PRIVATE -O2 -g -fno-omit-frame-pointer -Wno-shadow
                                                         -Wno-unused-variable)

endif()
