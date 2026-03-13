enable_testing()
find_package(GTest REQUIRED)

# Executable for the test
add_executable(main_class "${CMAKE_CURRENT_SOURCE_DIR}/test_openscofo.cpp")
target_link_libraries(main_class PRIVATE OpenScofo GTest::gtest_main)
add_test(NAME OScofo COMMAND main_class)

add_executable(score_tests "${CMAKE_CURRENT_SOURCE_DIR}/test_score.cpp")
target_link_libraries(score_tests PRIVATE OpenScofo GTest::gtest_main)
add_test(NAME Score COMMAND score_tests)

add_executable(math_tests "${CMAKE_CURRENT_SOURCE_DIR}/test_math.cpp")
target_link_libraries(math_tests PRIVATE OpenScofo GTest::gtest_main)
add_test(NAME Math COMMAND math_tests)

# Passa o caminho da pasta 'assets' para o código C++
target_compile_definitions(main_class PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/assets")
target_compile_definitions(score_tests PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/assets")
target_compile_definitions(math_tests PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/assets")

cpmaddpackage("gh:lieff/minimp3#master")
if(minimp3_ADDED)
    include_directories("${minimp3_SOURCE_DIR}")
endif()

add_executable(openscofo_perf "${CMAKE_CURRENT_SOURCE_DIR}/performance.cpp")
target_link_libraries(openscofo_perf PRIVATE OpenScofo)
target_compile_options(openscofo_perf PRIVATE -O3 -g -fno-omit-frame-pointer -Wno-shadow)

add_executable(score "${CMAKE_CURRENT_SOURCE_DIR}/score.cpp")
target_link_libraries(score PRIVATE OpenScofo)
target_compile_options(score PRIVATE -O3 -g -fno-omit-frame-pointer -Wno-shadow)
