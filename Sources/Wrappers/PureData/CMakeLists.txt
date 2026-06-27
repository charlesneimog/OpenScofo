project(openscofo~ C CXX)

file(GLOB_RECURSE OPENSCOFO_PDOBJ_SOURCES "*.cpp" "*.h")
pd_add_external(openscofo~ "${OPENSCOFO_PDOBJ_SOURCES}" LINK_LIBRARIES OpenScofo)

pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../README.md")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/openscofo~-help.pd")
pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/README.md")
pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/LICENSE")
pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/Documentation/assets/icon.png")
pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/Tests/miniaturas/Extras/ai-flute-model.onnx")

pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/Resources/Pieces/miniatura1.scofo")
pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/Resources/Pieces/miniatura1.wav")

pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/Resources/Train/py.o.train.pd_py")
pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/Resources/Train/py.o.train-help.pd")

pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/openscofo~.cpp")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/pdlua.cpp")

if(LINUX)
    target_link_options(openscofo_tilde PRIVATE -Wl,--exclude-libs,ALL)
elseif(APPLE)
    target_compile_options(openscofo_tilde PRIVATE -fvisibility=hidden)
endif()
