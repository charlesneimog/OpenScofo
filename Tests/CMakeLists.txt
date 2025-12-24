project(OScofoTests CXX)

cpmaddpackage("gh:pschatzmann/arduino-minimp3#master")

if(arduino-minimp3_ADDED)
    include_directories("${arduino-minimp3_SOURCE_DIR}")
endif()

add_executable(perf "${CMAKE_CURRENT_SOURCE_DIR}/perf.cpp")
target_link_libraries(perf PRIVATE OScofo)
