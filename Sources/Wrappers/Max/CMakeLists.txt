cmake_minimum_required(VERSION 3.25)
project(max_openscofo)

cpmaddpackage("gh:Cycling74/max-sdk#v8.2.0")

set(C74_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/max")
set(OPENSCOFO_MAX_MACOS_DEPLOYMENT_TARGET "${CMAKE_OSX_DEPLOYMENT_TARGET}")
include("${max-sdk_SOURCE_DIR}/source/max-sdk-base/script/max-pretarget.cmake")
set("${PROJECT_NAME}_EXTERN_OUTPUT_NAME" "openscofo~")

add_library(${PROJECT_NAME} MODULE openscofo~.cpp maxlua.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenScofo)
target_include_directories(${PROJECT_NAME} PRIVATE ${MAX_SDK_INCLUDES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MAX_SDK_MSP_INCLUDES})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

include("${max-sdk_SOURCE_DIR}/source/max-sdk-base/script/max-posttarget.cmake")