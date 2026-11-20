find_package(GTest QUIET)
if(NOT TARGET GTest::gtest_main)
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.17.0
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(googletest)
endif()

function(openscofo_add_training_test Target Source TestName)
    add_executable(${Target} "${Source}")
    target_link_libraries(${Target} PRIVATE openscofo_training GTest::gtest_main AudioFile)
    target_include_directories(${Target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
    add_test(NAME ${TestName} COMMAND ${Target})
    set_tests_properties(${TestName} PROPERTIES LABELS "openscofo;training" TIMEOUT 300)
    set_property(GLOBAL APPEND PROPERTY OPENSCOFO_TEST_TARGETS ${Target})
endfunction()

openscofo_add_training_test(training_api_tests test_trainer.cpp TrainingApi)
openscofo_add_training_test(training_e2e_tests test_training_e2e.cpp TrainingEndToEnd)
