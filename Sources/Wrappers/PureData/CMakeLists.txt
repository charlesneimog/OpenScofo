project(openscofo~ C CXX)

file(
    GLOB_RECURSE
    OPENSCOFO_PDOBJ_SOURCES
    "*.cpp"
    "*.h")
pd_add_external(
    openscofo~
    "${OPENSCOFO_PDOBJ_SOURCES}"
    LINK_LIBRARIES
    OpenScofo)

# Prepare package-only example assets in the build tree.  The repository keeps
# the original MP3 and test score; installers receive a WAV and a score whose
# model path matches the flat installed package layout.
set(OPENSCOFO_PACKAGE_ASSET_DIR "${CMAKE_BINARY_DIR}/OpenScofoPackageAssets")
set(OPENSCOFO_MINIATURA_SCORE_SOURCE
    "${CMAKE_CURRENT_SOURCE_DIR}/../../../Tests/04-miniaturas/Extras/miniatura1.scofo")
set(OPENSCOFO_MINIATURA_MP3_SOURCE
    "${CMAKE_CURRENT_SOURCE_DIR}/../../../Tests/04-miniaturas/Audios/miniatura1.mp3")
set(OPENSCOFO_MINIATURA_SCORE "${OPENSCOFO_PACKAGE_ASSET_DIR}/miniatura1.scofo")
set(OPENSCOFO_MINIATURA_WAV "${OPENSCOFO_PACKAGE_ASSET_DIR}/miniatura1.wav")

file(MAKE_DIRECTORY "${OPENSCOFO_PACKAGE_ASSET_DIR}")
# Remove the legacy package copy so incremental builds cannot keep shipping it.
file(REMOVE "${CMAKE_BINARY_DIR}/OpenScofo/miniatura1.mp3")
file(READ "${OPENSCOFO_MINIATURA_SCORE_SOURCE}" OPENSCOFO_MINIATURA_SCORE_CONTENT)
string(REPLACE "ONNXMODEL \"../../flute.onnx\"" "ONNXMODEL flute.onnx" OPENSCOFO_MINIATURA_SCORE_CONTENT
               "${OPENSCOFO_MINIATURA_SCORE_CONTENT}")
if(OPENSCOFO_MINIATURA_SCORE_CONTENT MATCHES "\.\./\.\./flute\\.onnx")
    message(FATAL_ERROR "Failed to rewrite miniatura1.scofo model path for packaging")
endif()
file(WRITE "${OPENSCOFO_MINIATURA_SCORE}" "${OPENSCOFO_MINIATURA_SCORE_CONTENT}")

find_program(OPENSCOFO_FFMPEG_EXECUTABLE NAMES ffmpeg REQUIRED)
execute_process(
    COMMAND
        "${OPENSCOFO_FFMPEG_EXECUTABLE}" -hide_banner -loglevel error -nostdin -y -i
        "${OPENSCOFO_MINIATURA_MP3_SOURCE}" -map_metadata -1 -c:a pcm_s16le -ar 48000 -ac 1
        "${OPENSCOFO_MINIATURA_WAV}"
    RESULT_VARIABLE OPENSCOFO_FFMPEG_RESULT
    ERROR_VARIABLE OPENSCOFO_FFMPEG_ERROR)
if(NOT OPENSCOFO_FFMPEG_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to convert miniatura1.mp3 to WAV: ${OPENSCOFO_FFMPEG_ERROR}")
endif()

pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../README.md")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/openscofo~-help.pd")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../LICENSE")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../Documentation/assets/icon.png")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../Tests/flute.onnx")

pd_add_datafile(openscofo~ "${OPENSCOFO_MINIATURA_SCORE}")
pd_add_datafile(openscofo~ "${OPENSCOFO_MINIATURA_WAV}")

pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../Resources/Train/py.o.train.pd_py")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/../../../Resources/Train/py.o.train-help.pd")

pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/openscofo~.cpp")
pd_add_datafile(openscofo~ "${CMAKE_CURRENT_SOURCE_DIR}/pdlua.cpp")

if(LINUX)
    target_link_options(openscofo_tilde PRIVATE -Wl,--exclude-libs,ALL)
    target_compile_options(openscofo_tilde PRIVATE -Wno-cast-function-type)
elseif(APPLE)
    target_compile_options(openscofo_tilde PRIVATE -fvisibility=hidden)
endif()
