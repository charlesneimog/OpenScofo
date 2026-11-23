find_package(Threads REQUIRED)

if(NOT TARGET nlohmann_json::nlohmann_json)
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.12.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(nlohmann_json)
endif()

if(NOT TARGET signalsmith-stretch)
    FetchContent_Declare(
        signalsmith-stretch
        GIT_REPOSITORY https://github.com/Signalsmith-Audio/signalsmith-stretch.git
        GIT_TAG 1.4.0
        GIT_SHALLOW TRUE
        EXCLUDE_FROM_ALL)
    FetchContent_MakeAvailable(signalsmith-stretch)
endif()

include("${PROJECT_SOURCE_DIR}/Resources/Cmake/CatBoost.cmake")

add_library(openscofo_training STATIC Audio.cpp CatBoostBackend.cpp Trainer.cpp)
add_library(OpenScofo::Trainer ALIAS openscofo_training)
target_include_directories(openscofo_training PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(
    openscofo_training
    PUBLIC OpenScofo
    PRIVATE OpenScofo::CatBoost AudioFile PFFFT::PFFFT signalsmith-stretch Threads::Threads
            nlohmann_json::nlohmann_json)
set_target_properties(openscofo_training PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(OPENSCOFO_BUILD_PD_OBJECT)
    pd_add_external(openscofo.train "${CMAKE_CURRENT_SOURCE_DIR}/PureData/openscofo.train.cpp" TARGET openscofo_train
                    LINK_LIBRARIES openscofo_training)
    pd_add_datafile(openscofo_train "${CMAKE_CURRENT_SOURCE_DIR}/PureData/openscofo.train-help.pd")

    set_target_properties(openscofo_train PROPERTIES SKIP_BUILD_RPATH ON INSTALL_RPATH "")
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_compile_options(
            openscofo_train PRIVATE -Os -ffunction-sections -fdata-sections -g0 -flto -Wno-cast-function-type
                                    -Wno-cast-function-type-mismatch)
        target_link_options(
            openscofo_train
            PRIVATE
            -fuse-ld=lld
            -flto
            -Wl,--gc-sections
            -Wl,--as-needed
            -Wl,--exclude-libs,ALL
            -Wl,-z,notext
            -Wl,-s)
    elseif(APPLE)
        target_compile_options(openscofo_train PRIVATE -fvisibility=hidden -Wno-cast-function-type-mismatch)
    endif()
endif()

if(OPENSCOFO_BUILD_TESTS)
    add_subdirectory(tests)
endif()
