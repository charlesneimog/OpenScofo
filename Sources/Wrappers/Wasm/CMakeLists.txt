if(DEFINED EMSCRIPTEN)
    # Set CXX_FLAGS to include necessary Emscripten flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -sALLOW_MEMORY_GROWTH=1")

    # Create a WebAssembly executable target with a different internal name
    add_executable(OpenScofoWasm "${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo.cpp")

    # Link OpenScofo with Boost and ensure the OpenScofo library is treated as part of it
    target_link_libraries(OpenScofoWasm PRIVATE OpenScofo embind)

    # Link additional Emscripten options
    target_link_options(OpenScofoWasm PRIVATE -sMODULARIZE=1 -sEXPORT_NAME='OpenScofo'
                        -sUSE_PTHREADS=1 -sWASM_WORKERS=1 -sWASM=1)

    target_sources(OpenScofoWasm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/OpenScofo.cpp")
    set_target_properties(OpenScofoWasm PROPERTIES OUTPUT_NAME "OpenScofo" SUFFIX ".js")
    target_compile_options(OpenScofoWasm PUBLIC -matomics -mbulk-memory)
endif()
