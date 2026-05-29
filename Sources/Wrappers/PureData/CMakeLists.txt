project(openscofo~ C CXX)

file(GLOB_RECURSE OPENSCOFO_PDOBJ_SOURCES "*.cpp" "*.h")
pd_add_external(openscofo~ "${OPENSCOFO_PDOBJ_SOURCES}" LINK_LIBRARIES OpenScofo)

pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../README.md")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/openscofo~-help.pd")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../Tests/assets/bwv-1013.wav")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../Tests/assets/bwv-1013.txt")

if(LINUX)
    target_link_options(openscofo_tilde PRIVATE -Wl,--exclude-libs,ALL)
elseif(APPLE)
    target_compile_options(openscofo_tilde PRIVATE -fvisibility=hidden)
elseif(WIN32)
    # non needed?

endif()
