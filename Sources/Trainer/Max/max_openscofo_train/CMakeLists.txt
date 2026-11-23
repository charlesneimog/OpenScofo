include("${max-sdk_SOURCE_DIR}/source/max-sdk-base/script/max-pretarget.cmake")

if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET
        "${OPENSCOFO_MACOS_DEPLOYMENT_TARGET}"
        CACHE STRING "Minimum supported macOS version" FORCE)

    set(CMAKE_OSX_DEPLOYMENT_TARGET "${OPENSCOFO_MACOS_DEPLOYMENT_TARGET}")
endif()

set("${PROJECT_NAME}_EXTERN_OUTPUT_NAME" "openscofo.train")

add_library(${PROJECT_NAME} MODULE ../openscofo.train.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE OpenScofo::Trainer)

target_include_directories(${PROJECT_NAME} PRIVATE
    ${MAX_SDK_INCLUDES}
    ${MAX_SDK_MSP_INCLUDES}
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_20)

if(APPLE)
    find_library(COREFOUNDATION_FRAMEWORK CoreFoundation REQUIRED)
    target_link_libraries(${PROJECT_NAME} PRIVATE ${COREFOUNDATION_FRAMEWORK})
endif()

if(WIN32)
    target_compile_definitions(${PROJECT_NAME} PRIVATE NOMINMAX)
endif()

include("${max-sdk_SOURCE_DIR}/source/max-sdk-base/script/max-posttarget.cmake")