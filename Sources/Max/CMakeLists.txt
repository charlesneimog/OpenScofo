cmake_minimum_required(VERSION 3.25)
project(max_o.scofo)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

cpmaddpackage("gh:Cycling74/max-sdk#v8.2.0")

set(C74_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/max")
include("${max-sdk_SOURCE_DIR}/source/max-sdk-base/script/max-pretarget.cmake")

set("${PROJECT_NAME}_EXTERN_OUTPUT_NAME" "o.scofo~")

add_library(${PROJECT_NAME} MODULE o.scofo~.cpp maxlua.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenScofo)
target_include_directories(${PROJECT_NAME} PRIVATE ${MAX_SDK_INCLUDES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MAX_SDK_MSP_INCLUDES})

include("${max-sdk_SOURCE_DIR}/source/max-sdk-base/script/max-posttarget.cmake")
