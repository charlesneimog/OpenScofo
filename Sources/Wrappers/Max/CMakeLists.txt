cmake_minimum_required(VERSION 3.25)
project(max_openscofo)

FetchContent_Declare(
    max-sdk
    GIT_REPOSITORY https://github.com/Cycling74/max-sdk.git
    GIT_TAG v8.2.0
    EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(max-sdk)

set(C74_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/max")
include("${max-sdk_SOURCE_DIR}/source/max-sdk-base/script/max-pretarget.cmake")
if(APPLE)
    # Max SDK 8.2 forces this cache entry to 10.11. Restore both forms of the variable before creating the Max target so
    # it is not compiled against a macOS version where std::filesystem is unavailable.
    set(CMAKE_OSX_DEPLOYMENT_TARGET
        "${OPENSCOFO_MACOS_DEPLOYMENT_TARGET}"
        CACHE STRING "Minimum supported macOS version" FORCE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "${OPENSCOFO_MACOS_DEPLOYMENT_TARGET}")
endif()
set("${PROJECT_NAME}_EXTERN_OUTPUT_NAME" "openscofo~")

add_library(${PROJECT_NAME} MODULE openscofo~.cpp maxlua.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenScofo)
target_include_directories(${PROJECT_NAME} PRIVATE ${MAX_SDK_INCLUDES})
target_include_directories(${PROJECT_NAME} PRIVATE ${MAX_SDK_MSP_INCLUDES})
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

if(WIN32)
    target_compile_definitions(${PROJECT_NAME} PRIVATE NOMINMAX)
endif()

include("${max-sdk_SOURCE_DIR}/source/max-sdk-base/script/max-posttarget.cmake")

if(TARGET OpenScofo::Trainer)
    add_subdirectory("${CMAKE_SOURCE_DIR}/Sources/Trainer/Max/max_openscofo_train"
                     "${CMAKE_CURRENT_BINARY_DIR}/openscofo_train")
endif()
