project(OScofoTests CXX)

cpmaddpackage("gh:pschatzmann/arduino-minimp3#master")

if(arduino-minimp3_ADDED)
    include_directories("${arduino-minimp3_SOURCE_DIR}")
endif()

add_executable(test1 "${CMAKE_CURRENT_SOURCE_DIR}/test1.cpp")
target_link_libraries(test1 PRIVATE OScofo)
