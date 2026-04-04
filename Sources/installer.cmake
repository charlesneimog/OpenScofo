include_guard(GLOBAL)

include(CPackComponent)

set(CPACK_PACKAGE_NAME "OpenScofo")
set(CPACK_PACKAGE_VENDOR "Charles K. Neimog")
set(CPACK_PACKAGE_DESCRIPTION "OpenScofo is an open-source library focused on research in score following (automatic tracking of musical scores) in contemporary and electroacoustic music contexts. Developed by Charles K. Neimog, the system uses AI and probabilistic models to synchronize digital scores in real time, and is compatible with the web.
")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://charlesneimog.github.io/OpenScofo/")
set(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/Documentation/assets/logo.svg")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "OpenScofo")
set(CPACK_MONOLITHIC_INSTALL OFF)

if(WIN32)
    set(OSCOFO_PD_DEST "$ENV{USERPROFILE}/Documents/Pd/externals/OpenScofo")
    set(OSCOFO_MAX_DEST "$ENV{USERPROFILE}/Documents/Max 9/Packages/OpenScofo")
    set(OSCOFO_SC_DEST "$ENV{APPDATA}/SuperCollider/Extensions/OpenScofo")
    set(OSCOFO_CSOUND_DEST "$ENV{PROGRAMFILES}/Csound/plugins64/OpenScofo")
    set(OSCOFO_VAMP_DEST "$ENV{PROGRAMFILES}/Vamp Plugins")
elseif(APPLE)
    set(OSCOFO_PD_DEST "$ENV{HOME}/Documents/Pd/externals/OpenScofo")
    set(OSCOFO_MAX_DEST "$ENV{HOME}/Documents/Max 8/Packages/OpenScofo")
    set(OSCOFO_SC_DEST "$ENV{HOME}/Library/Application Support/SuperCollider/Extensions/OpenScofo")
    set(OSCOFO_CSOUND_DEST "lib/csound/plugins64/OpenScofo")
    set(OSCOFO_VAMP_DEST "$ENV{HOME}/Library/Audio/Plug-Ins/Vamp")
else()
    set(OSCOFO_PD_DEST "lib/pd-externals/OpenScofo")
    set(OSCOFO_MAX_DEST "lib/max/OpenScofo")
    set(OSCOFO_SC_DEST "lib/SuperCollider/OpenScofo")
    set(OSCOFO_CSOUND_DEST "lib/csound/plugins64/OpenScofo")
    set(OSCOFO_VAMP_DEST "lib/vamp")
endif()

cpack_add_component_group(Integrations
    DISPLAY_NAME "Integrations"
    DESCRIPTION "OpenScofo plugins and externals for host environments."
)

cpack_add_component(max
    DISPLAY_NAME "Max (Cycling 74)"
    DESCRIPTION "Install the OpenScofo external and help patch for Max/MSP."
    GROUP Integrations
)

cpack_add_component(puredata
    DISPLAY_NAME "Pure Data"
    DESCRIPTION "Install the OpenScofo external and help patch for Pure Data."
    GROUP Integrations
)

cpack_add_component(csound
    DISPLAY_NAME "CSound"
    DESCRIPTION "Install the OpenScofo CSound plugin and example patch."
    GROUP Integrations
)

cpack_add_component(supercollider
    DISPLAY_NAME "Super Collider"
    DESCRIPTION "Install the OpenScofo SuperCollider plugin, class, and examples."
    GROUP Integrations
)

cpack_add_component(vamp
    DISPLAY_NAME "Vamp and Sonic Visualizer"
    DESCRIPTION "Install the OpenScofo Vamp plugin for Vamp hosts and Sonic Visualizer."
    GROUP Integrations
)

if(TARGET max_o.scofo)
    install(TARGETS max_o.scofo
        LIBRARY DESTINATION "${OSCOFO_MAX_DEST}"
        RUNTIME DESTINATION "${OSCOFO_MAX_DEST}"
        BUNDLE DESTINATION "${OSCOFO_MAX_DEST}"
        COMPONENT max
    )
endif()

install(FILES
    "${CMAKE_SOURCE_DIR}/Sources/Wrappers/Max/o.scofo~.maxhelp"
    DESTINATION "${OSCOFO_MAX_DEST}"
    COMPONENT max
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.wav"
    DESTINATION "${OSCOFO_MAX_DEST}"
    RENAME "bwv1013.wav"
    COMPONENT max
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt"
    DESTINATION "${OSCOFO_MAX_DEST}"
    RENAME "bwv1013.txt"
    COMPONENT max
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.wav"
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.txt"
    DESTINATION "${OSCOFO_MAX_DEST}"
    COMPONENT max
)

if(TARGET o.scofo_tilde)
    install(TARGETS o.scofo_tilde
        LIBRARY DESTINATION "${OSCOFO_PD_DEST}"
        RUNTIME DESTINATION "${OSCOFO_PD_DEST}"
        BUNDLE DESTINATION "${OSCOFO_PD_DEST}"
        COMPONENT puredata
    )
endif()

install(FILES
    "${CMAKE_SOURCE_DIR}/Sources/Wrappers/PureData/o.scofo~-help.pd"
    DESTINATION "${OSCOFO_PD_DEST}"
    COMPONENT puredata
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.wav"
    DESTINATION "${OSCOFO_PD_DEST}"
    RENAME "bwv1013.wav"
    COMPONENT puredata
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt"
    DESTINATION "${OSCOFO_PD_DEST}"
    RENAME "bwv1013.txt"
    COMPONENT puredata
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.wav"
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.txt"
    DESTINATION "${OSCOFO_PD_DEST}"
    COMPONENT puredata
)

if(TARGET OScofoCSound)
    install(TARGETS OScofoCSound
        LIBRARY DESTINATION "${OSCOFO_CSOUND_DEST}"
        RUNTIME DESTINATION "${OSCOFO_CSOUND_DEST}"
        BUNDLE DESTINATION "${OSCOFO_CSOUND_DEST}"
        COMPONENT csound
    )
endif()

install(FILES
    "${CMAKE_SOURCE_DIR}/Sources/Wrappers/CSound/examples/1-score-follow.csd"
    DESTINATION "${OSCOFO_CSOUND_DEST}"
    COMPONENT csound
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.wav"
    DESTINATION "${OSCOFO_CSOUND_DEST}"
    RENAME "bwv1013.wav"
    COMPONENT csound
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt"
    DESTINATION "${OSCOFO_CSOUND_DEST}"
    RENAME "bwv1013.txt"
    COMPONENT csound
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.wav"
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.txt"
    DESTINATION "${OSCOFO_CSOUND_DEST}"
    COMPONENT csound
)

if(TARGET SCOpenScofo)
    install(TARGETS SCOpenScofo
        LIBRARY DESTINATION "${OSCOFO_SC_DEST}"
        RUNTIME DESTINATION "${OSCOFO_SC_DEST}"
        BUNDLE DESTINATION "${OSCOFO_SC_DEST}"
        COMPONENT supercollider
    )
endif()

install(FILES
    "${CMAKE_SOURCE_DIR}/Sources/Wrappers/SuperCollider/OpenScofo.sc"
    "${CMAKE_SOURCE_DIR}/Sources/Wrappers/SuperCollider/examples/1-score-follow.scd"
    "${CMAKE_SOURCE_DIR}/Sources/Wrappers/SuperCollider/examples/2-descriptors-input.scd"
    DESTINATION "${OSCOFO_SC_DEST}"
    COMPONENT supercollider
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.wav"
    DESTINATION "${OSCOFO_SC_DEST}"
    RENAME "bwv1013.wav"
    COMPONENT supercollider
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt"
    DESTINATION "${OSCOFO_SC_DEST}"
    RENAME "bwv1013.txt"
    COMPONENT supercollider
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.wav"
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.txt"
    DESTINATION "${OSCOFO_SC_DEST}"
    COMPONENT supercollider
)

if(TARGET VampOScofo)
    install(TARGETS VampOScofo
        LIBRARY DESTINATION "${OSCOFO_VAMP_DEST}"
        RUNTIME DESTINATION "${OSCOFO_VAMP_DEST}"
        BUNDLE DESTINATION "${OSCOFO_VAMP_DEST}"
        COMPONENT vamp
    )
endif()

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.wav"
    DESTINATION "${OSCOFO_VAMP_DEST}"
    RENAME "bwv1013.wav"
    COMPONENT vamp
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/bwv-1013.txt"
    DESTINATION "${OSCOFO_VAMP_DEST}"
    RENAME "bwv1013.txt"
    COMPONENT vamp
)

install(FILES
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.wav"
    "${CMAKE_SOURCE_DIR}/Tests/assets/canticos.txt"
    DESTINATION "${OSCOFO_VAMP_DEST}"
    COMPONENT vamp
)

set(CPACK_COMPONENTS_ALL max puredata csound supercollider vamp)
set(CPACK_NSIS_COMPONENT_INSTALL ON)

if(WIN32)
    set(CPACK_GENERATOR "NSIS")
    set(CPACK_NSIS_DISPLAY_NAME "OpenScofo")
elseif(APPLE)
    set(CPACK_GENERATOR "productbuild")
endif()

include(CPack)
