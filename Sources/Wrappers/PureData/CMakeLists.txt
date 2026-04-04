project(openscofo~ C CXX)

file(GLOB_RECURSE OPENSCOFO_PDOBJ_SOURCES "*.cpp" "*.h")
pd_add_external(openscofo~ "${OPENSCOFO_PDOBJ_SOURCES}" LINK_LIBRARIES OpenScofo)

pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/README.md")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/openscofo~-help.pd")
pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.wav")
pd_add_datafile(openscofo~ "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt")

set_target_properties(openscofo_tilde PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
