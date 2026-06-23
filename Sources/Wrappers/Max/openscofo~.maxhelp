{
    "patcher": {
        "fileversion": 1,
        "appversion": {
            "major": 9,
            "minor": 1,
            "revision": 4,
            "architecture": "x64",
            "modernui": 1
        },
        "classnamespace": "box",
        "rect": [ 156.0, 192.0, 1295.0, 683.0 ],
        "showrootpatcherontab": 0,
        "showontab": 0,
        "boxes": [
            {
                "box": {
                    "id": "obj-3",
                    "maxclass": "newobj",
                    "numinlets": 0,
                    "numoutlets": 0,
                    "patcher": {
                        "fileversion": 1,
                        "appversion": {
                            "major": 9,
                            "minor": 1,
                            "revision": 4,
                            "architecture": "x64",
                            "modernui": 1
                        },
                        "classnamespace": "box",
                        "rect": [ 0.0, 26.0, 1295.0, 657.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "fontsize": 16.0,
                                    "id": "obj-6",
                                    "linecount": 5,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 12.0, 198.0, 1264.0, 114.0 ],
                                    "text": "To recognize these sounds, the model uses the audio descriptors provided by OpenScofo, including amplitude, onset, time-domain, pitch, and spectral features. Before using the model in performance, you must first train it with labeled examples of the techniques or sound categories you want it to identify. OpenScofo currently provides training workflows through either Python or Pure Data, allowing you to build custom classifiers tailored to your instrument, repertoire, and performance context. Once trained, the model can listen to incoming audio and recognize the techniques and sound categories you are performing in real time. For complete training instructions, dataset requirements, and implementation details, please refer to the documentation below.\n"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 16.0,
                                    "id": "obj-2",
                                    "linecount": 3,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 12.0, 38.0, 1264.0, 78.0 ],
                                    "text": "The AI model integrated into OpenScofo is designed for musicians who explore sounds beyond traditional instrumental performance. Whether you're working with extended techniques, percussive gestures, breath sounds, vocalizations, body percussion, or other unconventional sound sources, the model can recognize and track them in real time. Its goal is to support contemporary performance practices where musical expression extends well beyond conventional pitched material.\n"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 37.5, 587.0, 256.0, 20.0 ],
                                    "text": "Check complete AI docs"
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.592156862745098, 0.32941176470588235, 0.32941176470588235, 1.0 ],
                                    "id": "obj-31",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 11.5, 587.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-29",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 11.5, 615.0, 420.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/descriptors/ai/"
                                }
                            }
                        ],
                        "lines": [
                            {
                                "patchline": {
                                    "destination": [ "obj-29", 0 ],
                                    "source": [ "obj-31", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 154.0, 271.0, 133.0, 22.0 ],
                    "presentation": 1,
                    "presentation_linecount": 2,
                    "presentation_rect": [ 189.0, 185.0, 100.0, 35.0 ],
                    "text": "p Extended Techniques"
                }
            },
            {
                "box": {
                    "id": "obj-5",
                    "maxclass": "newobj",
                    "numinlets": 0,
                    "numoutlets": 0,
                    "patcher": {
                        "fileversion": 1,
                        "appversion": {
                            "major": 9,
                            "minor": 1,
                            "revision": 4,
                            "architecture": "x64",
                            "modernui": 1
                        },
                        "classnamespace": "box",
                        "rect": [ 0.0, 26.0, 1295.0, 657.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 37.5, 587.0, 256.0, 20.0 ],
                                    "text": "Check complete docs"
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.592156862745098, 0.32941176470588235, 0.32941176470588235, 1.0 ],
                                    "id": "obj-31",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 11.5, 587.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-29",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 11.5, 615.0, 401.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/score/intro"
                                }
                            },
                            {
                                "box": {
                                    "fontface": 0,
                                    "fontsize": 14.0,
                                    "id": "obj-22",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 464.0, 494.0, 429.0, 22.0 ],
                                    "text": "SCORE EXAMPLE",
                                    "textcolor": [ 1.0, 0.5215686274509804, 0.5215686274509804, 1.0 ],
                                    "textjustification": 1
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-21",
                                    "linecount": 5,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 464.0, 530.0, 428.0, 85.0 ],
                                    "text": "NOTE C4 2\n     sendto myreceiver [hello world 1 2 3]  \n\nNOTE D4 2\n     delay 2000 ms sendto myreceiver [hello world 1 2 3]  ",
                                    "textcolor": [ 1.0, 0.5215686274509804, 0.5215686274509804, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-20",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 95.0, 173.0, 339.0, 20.0 ],
                                    "text": "Receive message hello world 1 2 3."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-19",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 166.0, 413.0, 339.0, 20.0 ],
                                    "text": "Receive message hello world 1 2 3 after 2000 ms."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-14",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 23.5, 455.0, 96.0, 22.0 ],
                                    "text": "hello world 1 2 3"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-15",
                                    "maxclass": "newobj",
                                    "numinlets": 0,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 23.5, 413.0, 133.0, 22.0 ],
                                    "text": "r myreceiver-after-2000"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-11",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 577.0, 366.0, 428.0, 22.0 ],
                                    "text": "delay 2000 ms sendto <RECEIVER_NAME> [<LIST OF DATA>]"
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "fontsize": 14.0,
                                    "id": "obj-12",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 19.0, 362.0, 69.0, 22.0 ],
                                    "text": "delay",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-13",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 86.0, 362.0, 419.0, 22.0 ],
                                    "text": "delay 2000 ms sendto myreceiver-after-2000 [hello world 1 2 3]"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-6",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 577.0, 130.0, 389.0, 22.0 ],
                                    "text": "sendto <RECEIVER_NAME> [<LIST OF DATA>]"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-5",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 13.0, 214.0, 96.0, 22.0 ],
                                    "text": "hello world 1 2 3"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-3",
                                    "maxclass": "newobj",
                                    "numinlets": 0,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 13.0, 172.0, 75.0, 22.0 ],
                                    "text": "r myreceiver"
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "fontsize": 14.0,
                                    "id": "obj-2",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 16.0, 130.0, 69.0, 22.0 ],
                                    "text": "sendto",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-1",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 86.0, 130.0, 250.0, 22.0 ],
                                    "text": "sendto myreceiver [hello world 1 2 3]"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 16.0,
                                    "id": "obj-9",
                                    "linecount": 2,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 16.0, 10.0, 616.0, 42.0 ],
                                    "text": "OpenScofo support ACTIONS as consequences of EVENTS. In another words, when event is detect it trigger the actions."
                                }
                            }
                        ],
                        "lines": [
                            {
                                "patchline": {
                                    "destination": [ "obj-14", 1 ],
                                    "source": [ "obj-15", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-5", 1 ],
                                    "source": [ "obj-3", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-29", 0 ],
                                    "source": [ "obj-31", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 122.0, 240.0, 57.0, 22.0 ],
                    "presentation": 1,
                    "presentation_rect": [ 174.0, 170.0, 100.0, 22.0 ],
                    "text": "p Actions"
                }
            },
            {
                "box": {
                    "id": "obj-2",
                    "maxclass": "newobj",
                    "numinlets": 0,
                    "numoutlets": 0,
                    "patcher": {
                        "fileversion": 1,
                        "appversion": {
                            "major": 9,
                            "minor": 1,
                            "revision": 4,
                            "architecture": "x64",
                            "modernui": 1
                        },
                        "classnamespace": "box",
                        "rect": [ 0.0, 26.0, 1295.0, 657.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 37.5, 587.0, 256.0, 20.0 ],
                                    "text": "Check complete docs"
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.592156862745098, 0.32941176470588235, 0.32941176470588235, 1.0 ],
                                    "id": "obj-31",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 11.5, 587.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-29",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 11.5, 615.0, 401.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/score/intro"
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-70",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 308.8000046014786, 296.0, 289.0, 24.0 ],
                                    "text": "Check all descriptors on Documentation!"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-66",
                                    "maxclass": "multislider",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 440.0000065565109, 402.40000599622726, 322.0, 106.0 ],
                                    "setminmax": [ -100.0, 100.0 ],
                                    "size": 40
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 14.0,
                                    "id": "obj-62",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 300.8000044822693, 111.20000165700912, 619.200009226799, 22.0 ],
                                    "text": "Audio Descriptors",
                                    "textjustification": 1
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 12.0,
                                    "id": "obj-61",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 372.00000554323196, 252.80000376701355, 548.0000081658363, 20.0 ],
                                    "text": "Spectral centroid indicates the center of mass of a sound's spectrum.  Same a librosa on Python."
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 12.0,
                                    "id": "obj-60",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 300.8000044822693, 252.80000376701355, 65.0, 20.0 ],
                                    "text": "centroid"
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 12.0,
                                    "id": "obj-59",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 300.8000044822693, 206.8000007867813, 64.80000096559525, 20.0 ],
                                    "text": "logmel"
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 12.0,
                                    "id": "obj-58",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 300.8000044822693, 181.20000040531158, 64.80000096559525, 20.0 ],
                                    "text": "spread"
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 12.0,
                                    "id": "obj-57",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 300.8000044822693, 154.0, 64.80000096559525, 20.0 ],
                                    "text": "flatness"
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 12.0,
                                    "id": "obj-56",
                                    "linecount": 2,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 372.00000554323196, 206.8000007867813, 648.0000096559525, 33.0 ],
                                    "text": "The log-mel spectrum represents how the energy of a sound is distributed across perceptual frequency bands. Same as librosa on Python."
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 12.0,
                                    "id": "obj-55",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 372.00000554323196, 181.20000040531158, 606.4000090360641, 20.0 ],
                                    "text": "Spectral spread quantifies how dispersed the energy is around the spectral centroid. Same a librosa on Python."
                                }
                            },
                            {
                                "box": {
                                    "fontname": "Arial",
                                    "fontsize": 12.0,
                                    "id": "obj-9",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 372.00000554323196, 154.0, 496.00000739097595, 20.0 ],
                                    "text": "Spectral flatness indicates how noisy versus tonal a sound is. Same a librosa on Python."
                                }
                            },
                            {
                                "box": {
                                    "format": 6,
                                    "id": "obj-51",
                                    "maxclass": "flonum",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 521.8, 374.40000557899475, 68.00000101327896, 22.0 ]
                                }
                            },
                            {
                                "box": {
                                    "format": 6,
                                    "id": "obj-49",
                                    "maxclass": "flonum",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 358.6, 374.40000557899475, 68.00000101327896, 22.0 ]
                                }
                            },
                            {
                                "box": {
                                    "format": 6,
                                    "id": "obj-48",
                                    "maxclass": "flonum",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 276.8000041246414, 374.40000557899475, 68.00000101327896, 22.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-6",
                                    "maxclass": "newobj",
                                    "numinlets": 5,
                                    "numoutlets": 5,
                                    "outlettype": [ "", "", "", "", "" ],
                                    "patching_rect": [ 276.8000041246414, 335.0, 345.66666116714475, 22.0 ],
                                    "text": "route flatness spread logmel centroid"
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-27",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 65.0, 71.0, 113.0, 24.0 ],
                                    "text": "Play the music."
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-23",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 193.0, 14.0, 289.0, 24.0 ],
                                    "text": "Load record, on concerts, replace with mic."
                                }
                            },
                            {
                                "box": {
                                    "hiderwff": 1,
                                    "id": "obj-18",
                                    "maxclass": "playbar",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "int" ],
                                    "patching_rect": [ 42.0, 100.0, 196.0, 19.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-1",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 10.0, 15.0, 144.0, 22.0 ],
                                    "text": "open miniatura1.mp3"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-2",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 2,
                                    "outlettype": [ "signal", "bang" ],
                                    "patching_rect": [ 10.0, 153.0, 47.0, 22.0 ],
                                    "text": "sfplay~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-21",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 0,
                                    "patching_rect": [ 31.5, 212.0, 35.0, 22.0 ],
                                    "text": "dac~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-10",
                                    "maxclass": "toggle",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 113.0, 214.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-7",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 113.0, 245.0, 103.0, 22.0 ],
                                    "text": "set description $1"
                                }
                            },
                            {
                                "box": {
                                    "color": [ 0.08235294117647059, 0.4745098039215686, 0.08627450980392157, 1.0 ],
                                    "id": "obj-3",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "float", "float", "" ],
                                    "patching_rect": [ 10.400000154972076, 297.0, 285.40000396966934, 22.0 ],
                                    "text": "openscofo~ flatness spread logmel centroid"
                                }
                            },
                            {
                                "box": {
                                    "background": 1,
                                    "bgcolor": [ 0.792294779904983, 0.561867873693578, 0.668313324973563, 1.0 ],
                                    "bgoncolor": [ 0.55, 0.55, 0.55, 1.0 ],
                                    "fontface": 1,
                                    "hint": "",
                                    "id": "obj-53",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 143.0, 216.0, 20.0, 20.0 ],
                                    "rounded": 60.0,
                                    "saved_attribute_attributes": {
                                        "bgcolor": {
                                            "expression": "themecolor.lesson_downloading"
                                        }
                                    },
                                    "text": "1",
                                    "textcolor": [ 0.34902, 0.34902, 0.34902, 1.0 ],
                                    "textoncolor": [ 1.0, 1.0, 1.0, 1.0 ],
                                    "textovercolor": [ 0.2, 0.2, 0.2, 1.0 ],
                                    "usebgoncolor": 1,
                                    "usetextovercolor": 1
                                }
                            },
                            {
                                "box": {
                                    "background": 1,
                                    "bgcolor": [ 0.792294779904983, 0.561867873693578, 0.668313324973563, 1.0 ],
                                    "bgoncolor": [ 0.55, 0.55, 0.55, 1.0 ],
                                    "fontface": 1,
                                    "hint": "",
                                    "id": "obj-17",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 39.0, 73.0, 20.0, 20.0 ],
                                    "rounded": 60.0,
                                    "saved_attribute_attributes": {
                                        "bgcolor": {
                                            "expression": "themecolor.lesson_downloading"
                                        }
                                    },
                                    "text": "3",
                                    "textcolor": [ 0.34902, 0.34902, 0.34902, 1.0 ],
                                    "textoncolor": [ 1.0, 1.0, 1.0, 1.0 ],
                                    "textovercolor": [ 0.2, 0.2, 0.2, 1.0 ],
                                    "usebgoncolor": 1,
                                    "usetextovercolor": 1
                                }
                            },
                            {
                                "box": {
                                    "background": 1,
                                    "bgcolor": [ 0.792294779904983, 0.561867873693578, 0.668313324973563, 1.0 ],
                                    "bgoncolor": [ 0.55, 0.55, 0.55, 1.0 ],
                                    "fontface": 1,
                                    "hint": "",
                                    "id": "obj-15",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 163.0, 16.0, 20.0, 20.0 ],
                                    "rounded": 60.0,
                                    "saved_attribute_attributes": {
                                        "bgcolor": {
                                            "expression": "themecolor.lesson_downloading"
                                        }
                                    },
                                    "text": "2",
                                    "textcolor": [ 0.34902, 0.34902, 0.34902, 1.0 ],
                                    "textoncolor": [ 1.0, 1.0, 1.0, 1.0 ],
                                    "textovercolor": [ 0.2, 0.2, 0.2, 1.0 ],
                                    "usebgoncolor": 1,
                                    "usetextovercolor": 1
                                }
                            }
                        ],
                        "lines": [
                            {
                                "patchline": {
                                    "destination": [ "obj-2", 0 ],
                                    "source": [ "obj-1", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-7", 0 ],
                                    "source": [ "obj-10", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-2", 0 ],
                                    "source": [ "obj-18", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-21", 1 ],
                                    "order": 0,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-21", 0 ],
                                    "order": 1,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-3", 0 ],
                                    "order": 2,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-6", 0 ],
                                    "source": [ "obj-3", 2 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-29", 0 ],
                                    "source": [ "obj-31", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-48", 0 ],
                                    "source": [ "obj-6", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-49", 0 ],
                                    "source": [ "obj-6", 1 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-51", 0 ],
                                    "source": [ "obj-6", 3 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-66", 0 ],
                                    "source": [ "obj-6", 2 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-3", 0 ],
                                    "source": [ "obj-7", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 45.0, 116.0, 104.0, 22.0 ],
                    "presentation": 1,
                    "presentation_linecount": 2,
                    "presentation_rect": [ 144.0, 140.0, 100.0, 35.0 ],
                    "text": "p Listener Module"
                }
            },
            {
                "box": {
                    "id": "obj-4",
                    "maxclass": "newobj",
                    "numinlets": 0,
                    "numoutlets": 0,
                    "patcher": {
                        "fileversion": 1,
                        "appversion": {
                            "major": 9,
                            "minor": 1,
                            "revision": 4,
                            "architecture": "x64",
                            "modernui": 1
                        },
                        "classnamespace": "box",
                        "rect": [ 0.0, 26.0, 1295.0, 657.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 37.5, 587.0, 256.0, 20.0 ],
                                    "text": "Check complete docs"
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.592156862745098, 0.32941176470588235, 0.32941176470588235, 1.0 ],
                                    "id": "obj-31",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 11.5, 587.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-29",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 11.5, 615.0, 401.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/score/intro"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-18",
                                    "linecount": 10,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 82.5, 352.0, 764.0, 163.0 ],
                                    "text": "UTECH events describes events from extended/percussive techniques **that has NO defined pitches**. It is defined as UTECH <LABEL> <DURATION>. Where:\n\nUTECH <LABEL> <DURATION>\n\n    <LABEL>: Name of the technique (e.g., slap, jet-whistle).\n    <DURATION>: Duration in beats (integer or float).\n\nThese events require an ONNX model to recognize the technique. Training the model involves placing audio examples in labeled folders and running the py.train-onnx object, generating a .onnx file to load using TIMBREMODEL. "
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-16",
                                    "linecount": 9,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 81.5, 167.0, 766.0, 147.0 ],
                                    "text": "PTECH events describes events from extended/percussive techniques that has defined pitches (Pitched Techniques). It is defined as PTECH <LABEL> <PITCH> <DURATION>. Where:\n\n    <LABEL>: Name of the technique (e.g., pizz, tongue-ram).\n    <PITCH>: The pitch of the note (e.g., C4, D#5). Required for pitched techniques.\n    <DURATION>: Duration in beats (integer or float).\n\nThese events require an ONNX model to recognize the technique. Training the model involves placing audio examples in labeled folders and running the py.train-onnx object, generating a .onnx file to load using TIMBREMODEL. "
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-14",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 81.5, 125.0, 904.0, 22.0 ],
                                    "text": "REST events describes rests. It must be defined as REST <DURATION>."
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-13",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 81.5, 79.0, 905.0, 22.0 ],
                                    "text": "TRILL events describes trill and tremolo events. It must be defined as TRILL (<PITCH1> <PITCH2>) <DURATION>."
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-11",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 81.0, 41.0, 906.0, 22.0 ],
                                    "text": "CHORD events describes chords, stable multiphonics, and others events. It must be defined as CHORD (<PITCH1> <PITCH2>) <DURATION>. "
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 14.0,
                                    "id": "obj-9",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 81.0, 7.0, 632.0, 22.0 ],
                                    "text": "NOTE events describes tradicional pitches. It must be defined as NOTE <PITCH> <DURATION>."
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "fontsize": 14.0,
                                    "id": "obj-7",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 11.5, 352.0, 70.0, 22.0 ],
                                    "text": "UTECH",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "fontsize": 14.0,
                                    "id": "obj-6",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.5, 167.0, 70.0, 22.0 ],
                                    "text": "PTECH",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "fontsize": 14.0,
                                    "id": "obj-5",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.5, 125.0, 69.0, 22.0 ],
                                    "text": "REST",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "fontsize": 14.0,
                                    "id": "obj-4",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.5, 79.0, 72.0, 22.0 ],
                                    "text": "TRILL",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "fontsize": 14.0,
                                    "id": "obj-3",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.0, 41.0, 71.0, 22.0 ],
                                    "text": "CHORD",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "fontsize": 14.0,
                                    "id": "obj-2",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.0, 7.0, 69.0, 22.0 ],
                                    "text": "NOTE",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            }
                        ],
                        "lines": [
                            {
                                "patchline": {
                                    "destination": [ "obj-29", 0 ],
                                    "source": [ "obj-31", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 105.0, 205.0, 50.0, 22.0 ],
                    "presentation": 1,
                    "presentation_rect": [ 159.0, 155.0, 100.0, 22.0 ],
                    "text": "p Score"
                }
            },
            {
                "box": {
                    "id": "obj-1",
                    "maxclass": "newobj",
                    "numinlets": 0,
                    "numoutlets": 0,
                    "patcher": {
                        "fileversion": 1,
                        "appversion": {
                            "major": 9,
                            "minor": 1,
                            "revision": 4,
                            "architecture": "x64",
                            "modernui": 1
                        },
                        "classnamespace": "box",
                        "rect": [ 156.0, 218.0, 1295.0, 657.0 ],
                        "toolbarvisible": 0,
                        "enablehscroll": 0,
                        "enablevscroll": 0,
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 37.5, 587.0, 256.0, 20.0 ],
                                    "text": "Check complete docs"
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.592156862745098, 0.32941176470588235, 0.32941176470588235, 1.0 ],
                                    "id": "obj-31",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 11.5, 587.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-29",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 11.5, 615.0, 401.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/score/intro"
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-27",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 65.0, 71.0, 113.0, 24.0 ],
                                    "text": "Play the music."
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-23",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 193.0, 14.0, 289.0, 24.0 ],
                                    "text": "Load record, on concerts, replace with mic."
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-20",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 156.0, 348.0, 114.0, 24.0 ],
                                    "text": "Start to follow"
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-19",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 243.0, 313.0, 114.0, 24.0 ],
                                    "text": "Load the score"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-4",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 83.0, 242.5, 91.0, 22.0 ],
                                    "text": "print @popup 1"
                                }
                            },
                            {
                                "box": {
                                    "hiderwff": 1,
                                    "id": "obj-18",
                                    "maxclass": "playbar",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "int" ],
                                    "patching_rect": [ 42.0, 100.0, 196.0, 19.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-7",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patcher": {
                                        "fileversion": 1,
                                        "appversion": {
                                            "major": 9,
                                            "minor": 1,
                                            "revision": 4,
                                            "architecture": "x64",
                                            "modernui": 1
                                        },
                                        "classnamespace": "box",
                                        "rect": [ 332.0, 219.0, 1000.0, 780.0 ],
                                        "boxes": [
                                            {
                                                "box": {
                                                    "id": "obj-5",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 793.0, 174.0, 91.0, 22.0 ],
                                                    "text": "print @popup 1"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-3",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 665.0, 561.0, 91.0, 22.0 ],
                                                    "text": "print @popup 1"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "comment": "",
                                                    "id": "obj-27",
                                                    "index": 1,
                                                    "maxclass": "outlet",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 17.0, 724.0, 30.0, 30.0 ]
                                                }
                                            },
                                            {
                                                "box": {
                                                    "comment": "",
                                                    "id": "obj-2",
                                                    "index": 1,
                                                    "maxclass": "inlet",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 17.0, 7.0, 30.0, 30.0 ]
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.6313725490196078, 0.12156862745098039, 0.12156862745098039, 1.0 ],
                                                    "id": "obj-58",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 781.0, 502.0, 134.0, 22.0 ],
                                                    "text": "r buffer-pitchshift-space"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.6313725490196078, 0.12156862745098039, 0.12156862745098039, 1.0 ],
                                                    "id": "obj-52",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 754.0, 474.0, 98.0, 22.0 ],
                                                    "text": "r buffer-pitchshift"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.6313725490196078, 0.12156862745098039, 0.12156862745098039, 1.0 ],
                                                    "id": "obj-54",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 728.0, 446.0, 63.0, 22.0 ],
                                                    "text": "r pitchshift"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.6313725490196078, 0.12156862745098039, 0.12156862745098039, 1.0 ],
                                                    "id": "obj-56",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 665.0, 412.0, 99.0, 22.0 ],
                                                    "text": "r pitchshift-space"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.08235294117647059, 0.6666666666666666, 0.13725490196078433, 1.0 ],
                                                    "id": "obj-43",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 838.0, 88.0, 81.0, 22.0 ],
                                                    "text": "r delay-space"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.08235294117647059, 0.6666666666666666, 0.13725490196078433, 1.0 ],
                                                    "id": "obj-45",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 809.0, 59.0, 59.0, 22.0 ],
                                                    "text": "r delay"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.08235294117647059, 0.6666666666666666, 0.13725490196078433, 1.0 ],
                                                    "id": "obj-47",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 781.0, 31.0, 81.0, 22.0 ],
                                                    "text": "r delay-space"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.19607843137254902, 0.1607843137254902, 0.7607843137254902, 1.0 ],
                                                    "id": "obj-37",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 518.0, 146.0, 72.0, 22.0 ],
                                                    "text": "r buffer-play"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.19607843137254902, 0.1607843137254902, 0.7607843137254902, 1.0 ],
                                                    "id": "obj-35",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 491.0, 116.0, 108.0, 22.0 ],
                                                    "text": "r buffer-play-space"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "color": [ 0.19607843137254902, 0.1607843137254902, 0.7607843137254902, 1.0 ],
                                                    "id": "obj-28",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 472.0, 93.0, 84.0, 22.0 ],
                                                    "text": "r buffer-record"
                                                }
                                            }
                                        ],
                                        "lines": [
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-28", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-35", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-37", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-43", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-45", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-47", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-52", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-54", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-56", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-58", 0 ]
                                                }
                                            }
                                        ]
                                    },
                                    "patching_rect": [ 83.0, 200.0, 77.0, 22.0 ],
                                    "text": "p processing"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-63",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "patching_rect": [ 565.0, 231.0, 58.0, 22.0 ],
                                    "text": "loadbang"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-62",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 565.0, 260.0, 145.0, 22.0 ],
                                    "text": "filename miniatura1.scofo"
                                }
                            },
                            {
                                "box": {
                                    "code": "BPM 100\n\n// ╭───────────────────────────────╮\n// │     Check Docs for details    │\n// ╰───────────────────────────────╯\nONNXMODEL ai-flute-model.onnx\nONNXDESCRIPTORS mfcc logmel centroid flatness hfr flux zcr irregularity kurtosis\n\n\n// ╭───────────────────────────────╮\n// │          Measure number 1     │\n// ╰───────────────────────────────╯\nUTECH jet_whistle 1\nREST 0.5\t\nNOTE Bb4 1.5 // tied\n    sendto buffer-record [1 bb4-a4 1]\n    sendto buffer-record [2 bb4 1]\n    delay 1 tempo sendto buffer-record [2 bb4 0]\nNOTE A4 0.5\nREST 0.5\n    sendto buffer-record [1 bb4-a4 0]\n    sendto buffer-play-space [1 init -30 0]\n    sendto buffer-play [1 bb4-a4 0.4]\n\n// ╭───────────────────────────────╮\n// │          Measure number 2     │\n// ╰───────────────────────────────╯\nUTECH jet_whistle 1\nREST 0.5\nNOTE Bb4 1.5 // tied\nNOTE A4 0.5\nREST 0.5\n    sendto buffer-play-space [1 init 30 0]\n    sendto buffer-play [1 bb4 0.2]\n\n// ╭───────────────────────────────╮\n// │          Measure number 3     │\n// ╰───────────────────────────────╯\nPTECH pizzicato D4 0.5\n    sendto delay-space [1 init 0 0]\n    sendto delay-space [1 90 0 5000]\n    sendto delay [1 0.3 250 0.5 1]\nPTECH pizzicato A4 0.5\nREST 0.5\n    sendto delay [1 0.3 250 0.5 0]\n\nPTECH pizzicato D4 0.5\n    sendto delay-space [2 init -90 0]\n    sendto delay-space [2 -40 0 5000]\n    sendto delay [2 0.3 250 0.5 1]\nPTECH pizzicato A4 0.5\nREST 0.5\n    sendto delay [2 0.3 230 0.9 0]\nPTECH pizzicato D4 0.5\n    sendto buffer-play-space [3 init 50 0]\n    sendto buffer-play [3 bb4 0.2]\n\nPTECH pizzicato A4 0.5\n\n// Measure number 4\nREST 0.5\nPTECH pizzicato D4 0.5\nPTECH pizzicato A4 0.5\nREST 0.5\nUTECH jet_whistle 1\n    sendto buffer-play-space [1 init -40 0]\n    sendto buffer-play [1 bb4 0.4]\nREST 1\n    sendto buffer-pitchshift-space [1 init 60 0]\n    sendto buffer-pitchshift [1 bb4 -400 0.05]\n    sendto buffer-pitchshift-space [1 init -50 0]\n    delay 0.5 tempo sendto buffer-pitchshift [1 bb4-a4 -1200 0.05]\n\n// Measure number 5\nPTECH pizzicato D4 0.5\nPTECH pizzicato A4 0.5\n    sendto delay-space [1 init 0 0]\n    sendto delay-space [1 90 0 5000]\n    sendto delay [1 0.2 250 0.5 1]\nREST 0.5\nPTECH pizzicato D4 0.5\n    sendto delay-space [1 init 90 0]\n    sendto delay-space [1 -90 0 1200]\n    sendto delay [1 0.2 250 0.5 0.1]\nPTECH pizzicato Bb4 0.5\nREST 0.5\nPTECH pizzicato D4 0.5\nPTECH pizzicato B4 0.5\n\n// Measure number 6\nREST 0.5\nPTECH pizzicato D4 0.5\nPTECH pizzicato Bb4 0.5\n    sendto delay [1 0.2 250 0.5 0]\nREST 0.5\nUTECH jet_whistle 1\nUTECH jet_whistle 1\n\n// Measure number 7\nUTECH jet_whistle 1\nREST 3\n\n// Measure number 8\nNOTE D4 2\n\tsendto pitchshift-space [1 init 90 0]\n\tsendto pitchshift-space [2 init -90 0]\n\tsendto pitchshift-space [3 init 40 0]\n\tsendto pitchshift-space [4 init -40 0]\n\n\n\tsendto pitchshift [1 381 0.4]\n\tdelay 0.33 tempo sendto pitchshift [2 702 0.4]\n\tdelay 0.66 tempo sendto pitchshift [3 969 0.4]\n\tdelay 1 tempo sendto pitchshift [4 -969 0.4]\n\nREST 2\n\tsendto pitchshift [1 381 0]\n\tsendto pitchshift [2 702 0]\n\tsendto pitchshift [3 969 0]\n\tsendto pitchshift [4 -969 0]\n\n",
                                    "filename": "miniatura1.scofo",
                                    "fontface": 0,
                                    "fontname": "<Monospaced>",
                                    "fontsize": 12.0,
                                    "id": "obj-59",
                                    "maxclass": "text.codebox",
                                    "numinlets": 2,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "" ],
                                    "patching_rect": [ 565.0, 286.0, 639.0, 349.0 ],
                                    "saved_object_attributes": {
                                        "parameter_enable": 0,
                                        "parameter_mappable": 0
                                    },
                                    "textfile": {
                                        "filename": "miniatura1.scofo",
                                        "flags": 1,
                                        "embed": 0,
                                        "autowatch": 1
                                    }
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-48",
                                    "maxclass": "ezdac~",
                                    "numinlets": 2,
                                    "numoutlets": 0,
                                    "patching_rect": [ 27.0, 231.0, 45.0, 45.0 ]
                                }
                            },
                            {
                                "box": {
                                    "autofit": 1,
                                    "data": [ 88664, "png", "IBkSG0fBZn....PCIgDQRA..DPL..D.UHX....vxtAtG....DLmPIQEBHf.B7g.YHB..f.PRDEDU3wI6cdGlTTk919opp6dXHMjAQPIHRPjnDDEPEDEyJBlPbwzOcUDWiXXMf5ZbMfxZb00.XB00.fA.kfh.iRNNCHYPFfggYXlo6tp578Gy9d3T0zCLHwO4495hqo6pq5bdOmS0TU8zuAKkRo.gPHDBgPHDBgPHDxgIXev1.HDBgPHDBgPHDBgPNPBEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHj8C366G38JkB..dddk6iWNFyWWdaCy92yyCJkR+uxCl8guuef2qTpxzFTJUoF6GLorFylaSFeGJY2jCsI74UdddHYxjk6ueEts34dDBgPHG3wR8G4J2DBgPHjcKxC4ZaaqeussMbccQjHQJ2sgsss9AsUJkt81cHBWUd6qx53sssCzmJkBVVV5WumXSGJfo8mJ777fiiyAPKh7mIJOm+XdNn7cbgjIShnQite0FIDBgPHTPLBgPHj8qXJjkkk0tTHFSjGR1TPLKKqxsXMx9mJweJOOvcpNNSA9jwk4C0KiuvOf+ACRk8G9VdDOuyxxhBfQ1iv7bbSu6ZO463.PK37gJeugPHDB4vInfXDBgPH6Gnr7BrcricfJVwJtaOdSOvxba6IOzrH1iqqqVbLwlJOW9WNFQrq8z98PM1c1kondDR4AWWW333n+th3sj6tueYtex4ahX2Gp98GBgPHj+rAEDiPHDBY+HgEzZOgvheIdzj3QI6NRjHAhEKVorm8lG1VrA4A2MEQ5P0GhOr2rIaKU1ap1WBILk04Tk2vGNrGlY98G5oXDBgPHGXfWskPHDBY+DgeP2hJpH.TRHKt6HbNFRBsp8jP6KVrXZQqjjfuuue4NAd655hDIRn2eY7333D30lg7U4snAbfhck.WxbiYQKfBgQJOX98ZkRo8BSILmKOHeWQ9tib9HECiPHDB4.C7JtDBgPH6mPx4WxqSO8zgRoJWILayvnJ7CYWddfaWWW8wC.333f3wiWpDj+thHQhfXwho2eQ3K4g+EaQ5iC0dP9xxKd777fmmmV.ixSHtQHgIQhDZgqMCO5xy4Rld4o3wkRXWRHDBgPNvverxNEgPHDBYWhjOfDQhTJEhGONF4HGIRjHwtU7HI7qZZSaJtnK5hzhnEOdbjVZosa6e4AzkGz1yyC+6+8+FEVXgADzprv00E8nG8.coKcAwhECtttZwij1tvBKD4jSNvxxB0u90GQhD4PdOrx22GtttXkqbknV0pVnF0nFADl3P4P+jbnEx2shGONl5TmJl8rmMRjHARKsz1sdJokkkt3V355h1zl1f9zm9TpPblPHDBgr+ClCwHDBgP1OP3pBoDVUUnBUHfGd433ne34HQhnSF+hGdUm5TGrnEsHTyZVy8nv5SB+JGGG355h0t10hi4XNl.dblY+HgAoYx99kdoWB+0+5eE.AyqQJkBqcsqE23MdiXbiabPoT39u+6G2+8e+kKw5NPfXugqJmyblyDCZPCBYkUVH8zSGu0a8V3htnKBdddHVrXA7ntC1XVXFLycajC9DtnY7LOyyf65ttqR4gWQiFUGhzx26RUtr6tu66FO9i+3k6pHKgPHDBYuGJHFgPHDxAH1wN1AdgW3EPwEWLhGON1zl1DF+3GOxM2bQhDI.vNe.ZSuV5a9luA8oO84OTBwW7rqQMpQgAMnAUpOOZzn5vdrm8rmnksrknF0nFv11F8oO8AcqacqTOftuuOtu669vS7DOQfDC9XG6XwYcVm0ejol8KXVTAjva6DNgS.KXAKPuOQhDAqZUqB0u90ub68cGHHb0Gr7VHEHGXQ9tyzl1zvTm5TgmmGxM2bwhVzhvDm3D06movXNNNnScpS33O9iGMrgMDISlDsu8sGWvEbAbMlPHDB4.HTPLBgPHj8CDNo3Kdck4C7544g7yOeb5m9oiLyLSTwJVQbjG4QhrxJK..cNE5JuxqDuy67N51p7DVegET4bNmyAe8W+033NtiCKbgKLfmoYYYg+8+9eiAO3AuKaOOOOsWwT+5Wer4MuY8C4WgJTALrgML7fO3CtWWIK2WhoHWYmc1nYMqY5OSleG23FG5ae6K.Rck47fAhG9MzgNT344gm5odJTkpTE355pCUOxAWRlLoNrnM8fOQD5O3C9.bUW0UoKbChGY9zO8Si+1e6uoygc.676qlglLgPHDBY+K7psDBgPH6Gv72aRBQRWW2REtjUoJUA8oO8A..UspUEW60dsABEK4AqW25VmV7oxqXHhMjUVYgu8a+V..boW5kpCgS.fzRKMsmSI1nzOwiGGISlTmOyL8TrK7BuPsfZQiFEEWbw33O9iWmr5OXiosITm5TGzzl1T86kp1Yiabi0qKGJHFFPIqwEWbwXjibj3UdkWQulwbb1gNDIRj.EbBfR9Nmr8t10tBWWW820jvQtyctyZQMkPVVVSMaSBgPHDx9W3UbIDBgP1Of4C0JBIIgmH.Bj6tDgYRlLI5Uu5k9y.JQLsjIShoN0otG44HhvIJkBe+2+8v22GMoIMAmwYbFvyySaGwiGOf8FIRDDMZT344gzRKMDMZTsmrHsmuuOtoa5lPKZQKz66UcUWENqy5rBH11ASDOfS7LOeeeT0pVU7BuvKfJVwJp+7+4+7ehVzhVDHInu6RH5GnXMqYM.XmEHAKKqCYleObGSgek0C46HRgaHZzn5sEIRD84i1115BaQjHQzdPFPvp1JgPHDBY+KrJSRHDBgreFIg5aVgFiDIBTJERlLIpTkpjNot2wN1Qzst0ML8oO8.daxa9luItzK8RA.1iR71VVV3C+vOD..WxkbInpUsp5OKVrX5JdoDVggaaOOOsmsX5MLspUsBScpSEKYIKAUoJUAG2wcbGx44RlUxOgy9rOar7kubrjkrDbLGywfFzfFTJue6PEOzQDDyyyKfHojC9XFNzxqk+JeWQDh07ybccQkpTkBHPlmmm9bTV3DHDBgPNvwgF2wGgPHDxexv00U6oG111ZOsBnDANRjHArrrzhRIa222GW8Ue0ZQYj8869tuCKdwKVW4HKOnTJrfEr.L4IOYXYYgy+7OeTTQEoEkSrAQXNQLrjISFPTHyv3RBeR.fZTiZft10thV25Vq8DlCkBYRkRoCmSYb555h5Uu5gS4TNEzfFz.cALHbtb5PA9ke4W.PI1T5omN.PopRgjCNHmuXFBqlddIPIqaUnBUP+Yhf3EWbwZguDuDC.5uaa1FDBgPHj8ePAwHDBgP1OPjHQzObqovXRRQOVrXAdvWQvFaaabZm1ooSf8h.MQiFEe228ckauWRDNQDCqwMtw33NtiSmXuMC2KoukvwyT7NQ7MwCkLEpSde37h1gBgzmDFpJkRmWvLEe..5wpo8dnfsC.TbwEi25sdKsfKEUTQ..ZOJhbvEoHSHgHoo.WhnwwhECEWbwABkRyBYgHZqo2lQHDBgPNvAEDiPHDBY+HhvLBlBxDNGDIz3F2Xzu90uRI7wq+5uN.Pox+WlXFhiwiGGiZTiBJkBW20ccnxUtxZwvj9TdP8vdcloMkLYxRIjjX6Rk1CXmhlIBDHX9f9RB6WrUo8LyiRxwJikTIXU34Ly9VFWhfExbg49ENGPYNejJ61bsHUBWjpwqrsxZ7GN2QYdLSbhSDKe4KWKHloPjgmuKq11zlD6W7FwvywgIUi8vy8gOOz7yCe9S34nT0OxqS0bc3buW38UZ6vyABx3Nb6+GEIetArybGlrc4ex20LWKL+tS3pEp7crCUBYWBgPHj+rCuhKgPHDxgfLvANv.ODcxjIwBVvBvLm4L0Ovrj2uL87DyGxdIKYIX5Se5..nW8pWHYxj5iY2gHnDvNKF.ldIlYtORdfe4A5M8TFY+EAIDOxRD1QDHPZSIDFEOoyT.ASOayzqo.BJlPxjI0h7IgBZZoklNOtI1nHHkLVBKVnjGnDAFEB6kVhMK1nY3zIGq432b+LEDQ5+7yOe7rO6ypCAOwFE61bL644EXcxrHIH4PMSu+KVrX58Wlik8SpBpR+oTJjHQBcAHHrngx9Hyslm6Hm+366GXLKheZJbqrelq6.6TnToujwr47g44CRhrWN+RNWv00EwhESa+LGcQHDBgP.nfXDBgPHGRxIcRmjN4aWgJTAs..e8W+0ZgBLEqQDsvzaTl3DmHrrrvwcbGGZUqZkVrfxSXAFN4dG1KxL+bwaqjbzkovIh3KQhDQKNizVhXWwiGW6YOlEc.Q3CyOSD.xLuZArSOGRBCRgnQip+LybglYE+S9LoOMwzSe.JwSiLEoSDETDgIbUCTF+RXaJsgLVj4VI+RUTQEg65ttKLoIMo.1q4XTNFQTQy0kzRKMjLYRsvXR3uZl6wLGyRNpKZzn59yrJHFKVLsW.JBfIXJDXXg8DgzLsewVDuqxzC8LqFix5c37qkHJqYX5ZJ.pYHFaVYWMqRmldJFgPHDB4vaXUljPHDB4PPpcsqMF7fGLdoW5kzh.333f27MeSb629siJUoJUpiwTbhDIRfW+0ecDIRDbYW1koEOpvBKDQiFMPnkkJL83IQrnved73wQZoklV.pvBXXl6jLEpRDqP5iTUgKUJEl6bmKVvBV.xImbv5V25PKaYKQSZRSPO5QOBDRjldml3kPxmKUPSYbXNWVbwEiJTgJnE7wzitDaPD7RZKy7Q1RVxRvpW8pwl27lwF23Fwl27lQEpPEPMqYMwQcTGkNusEtxcJdqj47fqqK90e8Ww8ce2GlzjlToDSLQhDnhUrh..5hdfYnSJy0lyyxXFnjvZzT.QSArhEKl16vjh7frsvgumzGgyGahPWx3R9LYMR9mbdgXaaXCa.yXFy.YkUVn3hKFMu4MGWvEbAAFW.6LW0IhdYaaihKtXrxUtRrksrErpUsJroMsIrgMrADKVLT6ZWabjG4QhVzhVfV1xVpOVwC4nWhQHDBgPnfXDBgPHGhxEcQWDdoW5kzBXYYYgUspUgYNyYhS8TOU89IgYlD9YJkBKbgKDKYIKA..8oO8QKnSkqbk2shgIHBd455h+6+8+hbxIGjLYRr8sucjat4hhJpHr5UuZje94ie629M7y+7Oi5Tm5.fcJhSrXwv1291wDlvDvV25VQt4lK1111l9uabiaDqe8qGu5q9pnm8rm..Xtyct3QezGEe1m8YA7vJopX15V2Z7vO7CiK7BuP.f.g.n3sVVVVHmbxA+vO7CHu7xCEUTQ32+8eGtttHmbxAacqaE+9u+6HVrXXpScp5wrHThL1My2SJkBEVXg3G9ge.u1q8ZXricrZOTpr75td26digNzghd0qdoEjJYxj5vaLszRCu8a+13ttq6BaZSapTEaAQ7npW8pGXax1A.lwLlA5bm6LTJEJnfBv3F23v1291Q94mO1111ld8ZyadyX4Ke4XDiXDnO8oOAr4HQhfEu3Eie4W9ETTQEgBJn.76+9uijIShst0sh0t10hMu4MiS+zOc73O9iqEgT7JPWWW7ke4WhBKrPrksrEje94iMrgMfjIShbxIGjUVYg669tObkW4UpEHczidz31tsaCaYKaIPHT1oN0I7QezGgi5nNp.gzpr1TbwEie3G9AL5QOZ7du26EX91TnOQDtt28tia+1ucblm4YVpb5GgPHDB4vWnfXDBgPHGBhuuONwS7DQiZTivJW4J0gqlkkE9jO4SP26d2CDdgBhH.icriE..snEs.sqcsSG5fEUTQkpxPlJj7MluuO1111FFv.F.bbbBD1dR9CyrhSZ5oVh2TUPAEf90u9UppYoLNA.RO8zgRov+8+9ewUcUWExO+70UpSOOOsffVVVX9ye9n+8u+3AevGD2y8bOkZdP79nryNabYW1ko8LLQTMYecccQW6ZWCLmKdwTpx0ZYkUV3FuwaDSZRSB111nKcoK3bO2yEsoMsAUpRUBaZSaBKe4KGiabiCSaZSCQhDASXBS.SXBS.25sdqX3Ce3nxUtxZO3xL7F2zl1D.1Yxe2rxCZlurLSt7xbi7dKKKroMsIb4W9kqmekwr49ZtdIu2wwAScpSECYHCQaChWmIgBpqqKNyy7LC3YZh2Ct7kubbIWxknsSQvL47EIzNki4Idhm.OvC7...37NuyCadyaFScpSEVVVXFyXF39u+6Gu9q+55yAk9bYKaY3ZtlqA+3O9iPoT3jNoSB8su8Est0sFYjQFXCaXCXUqZU3K+xuD+zO8S..XpScp3G9ge.W+0e8XDiXDHVrXADNiPHDBgb3ITPLBgPHjCAw11FUnBU.CdvCFOxi7HARp4u7K+x3u+2+6nt0stAxmUh3IEWbw3Mey2DQhDAWy0bMZAnjbD0tSLLfcJ1kXGO1i8X.njvHbYKaY3ce22MPnFJ8qXGR9oxwwAokVZX3Ce3HZznXyadyXDiXDAx0XQiFEEWbwXricr3hu3KF111XPCZP3xu7KGYjQFXyadyXJSYJ3oe5mNPRb+AdfG.mvIbB3LOyyLPBZWBwwFzfFfG+web.ThPSSaZSCe228c.XmBJYNNM8DIy7ykmmGxN6rw.Fv.vBVvB..vvG9vwcbG2ARKszzdBm72a61tM79u+6iAO3Aq29y+7OOZXCaHt0a8VCHrkLV6W+5GRO8zQhDIva+1uMt0a8VCjW3V+5WuVHGybAlqqKpZUqpVLvpW8pqWqxO+7wK8RujtvBXdtkYHuJEef10t1gG5gdH333fBJn.L1wNVLu4MO87kjOw.fVjLYdpt0st3odpmBwiGGEVXg3i+3OFYmc1ZQwL8vqO8S+Tbe228gt0stgQLhQf5V25hN1wNpGuVVV3C9fO.uvK7BHszRSKb2JVwJv4bNmCxN6rgRovi9nOpdMPZeQvs69tua7lu4aha7FuQsM+u+2+aznF0HbO2y8Tt+dHgPHDB4OwnHDBgPHGTv22WoTJ0C8POjxxxRU6ZW6RsOyd1yVA.ksssB.J.nrrrTe1m8YJWW2R0dtttpe5m9IUjHQT.PkYlYF3ym+7mutcj1B.pkrjkDvljW644EnOj2O7gO7.sC.TqYMqoTsQprwO9i+3.8skkk5Ye1mUUu5UOU8qe8USZRSpT8qRoTe4W9kpZVyZp6OaaaUyadyUEVXgobd0zN777TISlT0idzCksssd9rqcsqkZ+877z8exjIUJkRca21sEXNKyLyTEOdbc6Ku1rs9q+0+pB.pnQip.fpZUqZpst0slx4Fyw6a7Fug19hDIhJRjHpMsoMUliO4Xccc0sommmx22W8Ue0Worrrzy0NNNpwMtwkx9UdsztaYKaQkVZoosA.nti63NJ05R30XkRol6bmqdL333n.f5sdq2Rs5UuZUrXwTsnEsP8a+1uoTJk5q+5uVO2J6q44ShcMzgNz.eWHyLyTONM2ujISpW2t9q+5CbLUqZUSkSN4TJ6c+Aqd0qVOua9uoO8oe.o+IDBgPH6ZXUljPHDB4PPT+upAYqacqQW5RWBD1XVVV3ce22sT4BIYe91u8agqqKZe6aOZcqas9yCWYE2cHd9kIx6aYKaYJ2eUnbokYRxW9rl27lGns..ty67NQd4kG9xu7KwodpmZopZg..my4bNA7tGkRgryNaLyYNy.yAlEW.S6wwwQOWJdnkoGWY54UxqiDIB15V2JdgW3EBzuO2y8bAB6Syjsuv0bMWiNGiEIRDrsssMsGWYN2H8ozdldPlmmGbcc0ITey9HbUizrhSJaqEsnEAl6M8zvvqAlglIPIda1wcbGmNGgIm6TVqwlzvF1v.gbokkERKszva7FuARjHAF4HGIZTiZToR1+x4nspUsBUqZUSaO4lat50.o+e1m8Y0y+l+UBcWWWWb8W+0q+rHQhfsu8sirxJqRYuDBgPHjC+fBhQHDBgbHHJiJ42.G3.0aWDK4y+7OGYmc1..5PUTB0sW8UeUXYYg92+9qyEVJib5T4IjI2cDV7D40lBLUdQ8+pliO8S+zn8su851xL+SIg724cdmWf9vyyC+7O+y5j6dYERnxwrmHFnzNETPAAB4uHQhfQMpQgbyM2.hCkHQBsfSVVVnwMtwAD1BnjvXTxQWGpikkEpXEqndNSDRr7r9VoJUo.gEqRovhVzhvS+zOMt5q9pwocZmFJt3hQZokF5Tm5DtvK7BCD5uO4S9jnxUtx59aqacq50O0+KrYG8nGMV+5WefbVlxHzYcbbPSaZSA.zgSoTXJHDBgPHDlCwHDBgPNDDyDdtjHyEumQDPZJSYJ3XNliAUnBU.ISlD111X9ye9XCaXCvxxBmy4bN51.XmUJwT4IW+Qsu8DB6YSx3w22G0pV0BCbfCDVVV5pvnHvgkkEbbbfuuOZVyZFZUqZEV3BWntcl+7me.QZR0XyruKuieYLV8pWcjQFYfBJn..rSOVK7bo3kXhGQI4rMwKo.ftRRdffTIN4dpfkokVZAlqJum2Hh9kd5oi3wiCKKK7nO5ihHQhfa3FtAnTJTgJTAnTJjd5oi268dOLkoLEjLYRz7l2bzrl0L.ry0fFzfFfpV0ph7yOe.TRU.0L+uADr.DHiUAOOOjVZog3wiGnvPPHDBgPN7E5gXDBgPHGBhH7fssMZRSZBNqy5rzdsjD5Z+m+y+A.6Tb.GGG7EewW..fS3DNAzpV0p.IGdoc2WTg81SEVIrmjIHdzSe5SePUqZUCjT7S09B.zt10t.s8N1wN..B3YPkWaYWg3QTUoJUAuxq7J5i011FCe3CGUu5UW2uRhaGnDgGk0HypvoTL.1W3gdGHPrSYdWJfBkGjB.PQEUjdNvxxBssssEcnCcHf2f444gJVwJhd26diy8bOWzrl0r.BuIU8yW9ke4.gE4ce22MpW8pmd+L8FMQ.YyvwLd73..nnhJ5O5TBgPHDB4OQPODiPHDB4PPj7AlH7zfFzfvW9keI.1Y93ZFyXFXVyZVnScpS..Hu7xCu5q9pv11FCX.CP6UUlTdqxjGHP7tJkRgS5jNIsPF.Pa6hG+HdMF.PFYjA.1o.HISlTGZcgGaoR.mxiWNY1e..W5kdoncsqcX4Ke4nd0qdnicri..ZA3hFMpVrQopJt0stU.TRHVJU0w8UBRVd4.YeYRhDIzquldv0fFzfBTEOMyWalq4xZoH.oqqK5e+6O5XG6HVwJVApQMpA5bm6rtcM8nRQXTKKKrhUrBs3jhmgcvZNgPHDBgbnETPLBgPHjCAQDFR7tnt10thJVwJhBKrP89HgMYG5PGfuuOxLyLQd4kG788wYdlmoVPAkRoEXyTjm8VJKgkRknOg2V3vVr5Uu5AD3xzNEuFSvLjCMy0XRXJFNz4BaGkGAwj9SrEeeezhVzBz7l27.d2jHtyZW6ZwRVxRvhVzhvrl0rvDm3Dwl1zl..BDhd+QxwZ6qo71+l4INILc2S5Cy4Ios5PG5..1442x9XdNpbLhnXNNN5vm8XO1iEMsoMUuewiGWe9vZW6ZwhVzhvJVwJvzm9zw27MeC1xV1RoDCaOYbPHDBgP9yKTPLBgPHjCgQ7Vli7HORLfAL.7tu66FvKXdsW60vPFxPPrXwvm+4eN..5RW5BZUqZUJqff6OEC3OhPOxwToJUIsMJhQENw8mLYR333TJAwL8BIffBvjJg3JO4PLQPFwyujb2ljiyhDIBl6bmKl4LmIF0nFEl5Tmp1lEuT6DOwSDYlYlHYxjH8zSGEUTQ61v47PMhGOt9bMfx+ZbZokF777BTHGbbbPqacqCD1i.POeJ43NWWWsW0433fhJpHDKVr.dVlreqd0qFe629s3C9fO.SaZSKv7ussMZe6aOVvBVf19k9gPHDBgPnfXDBgPHGBhH9hHjfRovkdoWJdm24c.PIB1DKVLrrksLjYlYhV0pVg27MeScXoYJhTznQ0h.kHQBsHE6KorDJY2EhfhXKUspUsT6uYX0466qEJSRR6hWEIGmH.SjH682dS3jytXKwiGGSXBS.O8S+zXxSdxZQZN5i9nwUbEWAZUqZEZSaZCZXCaHrssQ0qd0QZokFJpnhRYHct+h8UglojCv.1YtQq7fjS0LK3.MsoMEUtxUtTg2nH5nf3Qdh3il4ALGGGTXgEhINwIhm+4edLoIMIsHmMpQMBCX.C.mvIbB3XO1iEMtwMFNNNnxUtxA7xM5gXDBgPHD.JHFgPHDxAMLySRVVVHQhDAxgR.PK3B.P25V2Psqcswu+6+N.JIT7rsswXFyXPu5Uuzhb06d26.GuY+IhgEVrDQzi8FQT1SpbkhfMxXybLK1iosmJucSDJSlCSkPX6twS3JTnH7k7d4uNNNXxSdx31u8aG+xu7K5i+Vu0aE8qe8CcricLfmq455p8Lp3wiqqvgh3lhsENTAMsKeees2oIg8WjHQzdTkYQVvbdyza1j00XwhoyqWx9EN7MC2VtttXIKYIvwwQetkXS6t40vdkWhDIvIexmbf8I74ll+UVKDww788QjHQvjlzjvcdm2I90e8WAPId70Mdi2Hthq3JPm5TmJU30JUkRo8DgjAPf0gT898Fj7gm42IjuKevNjYIDBgPHk.EDiPHDB4fHVVVnvBKD999nfBJPKlgovXhvAUoJUACZPCBO+y+75vGLYxj3ce22EqYMqAwiGGm5odpnIMoIkYdzJ7C+KOj99BupJber6XWIfVYkL7MSD+BliyTIxQpZKyvlLU6ioXX9993se62FW8Ue05i8nO5iF+m+y+A8nG8Hv7WhDIPrXwzBWI6uTIPMWaELyaVlq+xZiDdfUnBU.tttAD7QFChXYREuLZzn5PzLrfVl8e3vKU1lyc87r...H.jDQAQkXyldHl4ms6vzK+j4yZW6ZWtOOyb7Kdl1a7FuAt4a9l0m6W25VW7AevGft28tqs+vBnUgJTAscHB.JUaRS6L7q2awTr2vmytmHbLgPHDBY+GkO+dmPHDBgrOkhKtX.ryjBu3AKRn0IUXRfc58VJkBW7EewHYxj5vIC.XaaaaXLiYL..3xu7KWKLBPPOfx78gwzijJOXJ7T4Imbs6vTLhxJo7WVDVPr8TgMLOFQPJQTJkRgu+6+dsXXhW+7QezGgd1ydVJwDMmuEgWLCeUWW2.qsoRPFy1PlakvdMr3UqYMqAUu5UGm5odpv22OffYomd55wfYHmF1dEAjDg4DgX25V2ZoDOq7l+sRUngVm5Tmx7yBirNH8227MeCttq65P73wguuO777v3F23vIexmLJpnhBHVJvNEYLu7xS2lg+Lyss+fT8cNYbQHDBgPN3CEDiPHDB4f.UnBUP64OabiaT6QKEUTQ58wLLFEQZZW6ZGZaaaq9grMyGXQiFEcu6cOkOvcY4cJlhssmfoGXYdr+Qam8j8SDgKb39smzdkkcZ11111Hu7xCCaXCKf.V2vMbCnScpS.XmhqHqAQhDQmarjJBpofahGKEFSOBLUgGpHtlD5jR3PtksrEr8sucje94GHzJE6QB+SfcVTBRjHgtMMCISIGyAThGUM+4O+.iw8DwFMCISQzsJUoJUtNVAYMpfBJ.28ce2.njywiDIBF5PGJZcqaMrrrP5omdopnkxwJUgRQLPQfv82hRIyUlgspY38RHDBgPN3CEDiPHDB4fDh2EsrksL.ThnJETPAAd.Z.DPXgXwhgq5ptJc3XIgJG.PO5QOPSaZSKUXQZxtR3nv4wqcEg2OSgj1W5sXk01ME0XWIxvtZ7H1ZpluDgUxJqrPlYlY.AWNkS4TJUaYFhjhWjshUrBs3Lxb61291Sosjp4LyvY0LzEk0aaaa76+9uCeeez3F2XsmFJEV.4bjnQipE0RxiV.kTIHMwL+Z433fu5q9JDMZT83prryTQplSkvYr7jX9EQE888w7m+7wbm6bQZokl1K65ZW6Zo5CIGsIGK.PVYkERjHQfwf44pkUn4t2hYNDKrPhrJWRHDBgbnATPLBgPHjCBHdD15W+5wrl0rzhbLu4Mu.d+E.BTMESjHAN2y8bAPIOXcznQ0dLz.Fv.JyG11TPjvIybIbMMss8F1c4kovIUdysUV6ap1lYXxY5AWkGAMjvR0LImaZqhXLKe4KG.kLmjVZoAeeeT4JWYc3IJGS3b6E.vq9pup90Rt.KbXvZ1eR6HBUYJLp3EWx1Dgsl8rmM..ZVyZl97ESuCKVrXHYxjABsy7xKO344oC8VAybV11111vnF0nBHx2dhfoldokX+lgx6tCyBav5V25.vNqtnNNNH8zSOPwC.nDA9jPRU3Mey2D.PGxo.Hf3Xoh8EdvkLNS04VDBgPHjCMfBhQHDBgbPjwO9wG38iYLiIf.X.HfmrDKVLzzl1TbVm0YoE0HVrXPoTnO8oOAN1T4kNRHwEVHJI79Jujp9XOIGjE1tJqPuLUBIDVLovBbjp1w7uluVBEwv1s3UWYjQFZaPRL9KbgKTKviHLlH1jH9z29seKdu268BLFbbbPQEUDRlLYoxsalBmTqZUKsfPRH2szktT.rSg5jBpvK9huHRO8zQm6bmCD1l.kHLVMqYMQspUszu211FewW7EZwTM8XLwFsrrvS+zOM1wN1QfD9uRovV25VK05QpH77YXQW2SHZzn5PMURV8qXEqPOdk10TzOOOO7Ye1mgO3C9.cgFPoTHRjHXKaYKAFy6uBgRI2lYhssM1wN1w979hPHDBgrmCEDiPHDB4f.VVV3m9oeB2y8bO5p.niiCd629swq8ZulVjGyGb27A+u7K+x0d2UhDIvoe5mNZTiZToxAUlGmRoPwEWbfJtmHjUgEVXfJB3tiTItg3cQ6JJKwGJOhRH1pYeqTJc3DJBjTdxgZJkRWDCLGOlBT444gF0nFoERQBev+0+5eg0t10paaQXLOOODMZTL9wOdbAWvE..fZVyZpCeNOOO7we7GGPjOwlEaxwwA0t10NvZVjHQvW8Uek19DQwd1m8Ywl1zlPKZQKvYe1ms1NDrrrPrXwv0ccWG.1YHKNlwLFsmSIddlzWIRj.O8S+z3e7O9Gne8qenMsoMAD8bYKaYAr4xBw6yLETSDer75AhhGfcrG6whhKt3.hD9O+m+SrwMtQceHBXBTxZ428ceGF3.GHRlLoVPPI+p8Ue0WoymWx+1axCdoBopwBDLOtoTp.I5eBgPHDxAObdnG5gdnC1FAgPHDxe1HUgy2u+6+NxLyLwDm3DwHFwHvPG5PC3sHh..iabiCSaZSC6XG6.4latv22G0pV0Jf2AUyZVS7rO6ypELYXCaXnCcnCADaY4Ke4XEqXEXCaXCXYKaYXJSYJ3odpmBqcsqU+P+hHJYmc1HRjHX6ae6XcqacXiabi..nhUrhA7tGKKKTPAEfQLhQf4Lm4DnZV13F2XzoN0o.IT7vIW7DIRfO4S9D7ce22EH2ZU4JWYzyd1Ss2tENg4KswJW4JwC9fOHxM2b06yF23FQe6aeQ8qe8KUnpIuVN1G8QeT76+9uq+7srksfd26diF0nFUp7ak3gXSaZSCqd0qVKDTAET.xLyLQqZUqP8pW8fssMhGONl+7mOdhm3Ivse62NRlLI9rO6yP26d2wW7EegVPj0u90iksrkgbxIG7Fuwafa+1ucbK2xsnmerssQ0qd0wjlzjvZVyZzi+e7G+QTXgEhnQihUtxUhW+0ec7vO7CC.f+0+5egi+3O9.m2YNGbjG4Qh25sdKcBl2yyCewW7EXaaaaXG6XGHYxjXwKdwXxSdx3ltoaBiZTiBcricDu8a+13q+5uFqXEqPKl1pV0pPFYjApUspEJt3hgqqKhEKVo73sILgIfQO5QWpp+3ocZmFpRUpB.JwaEMCaUQDRy0.WWWTspUMLqYMKjc1YqGSae6aGyd1yFst0sF0rl0DQhDAEVXgXwKdw3Idhm.21scaHd73XLiYLnm8rm3y+7OWK33ZW6Z0eu3Mey2D28ce23FtgaHf3c6NQgMyMXdddXMqYMH6ryFqd0qFKcoKESe5SGO9i+3XkqbkPoT54cKKKL24NWToJUIr4MuYroMsIr90udjHQBT0pVUsnck2PKkPHDBgrWfhPHDBgregDIRD3uSbhSTA.UrXwTNNN5WCf.+y11V+ZKKK0i+3OtRoTJWWWkmmmt8GzfFjJZznJ.nV9xWdo5+W5kdIcaHsq4qM+LSaQ11HG4H0s0jm7jU+k+xeQ06d2aUUpRUTNNN5wf7OGGGUcpScTWzEcQpANvApxJqrTJkRsl0rF00dsWq5htnKRU+5W+R0ex+pRUph5zNsSScsW60pVvBVfxyyS466q9fO3CTWwUbEpN1wNVp9KRjH522111V0kdoWpZjibjJWWWkRoTiXDiPMnAMHU25V2Br+hsKuuksrkp90u9ot1q8Z0GaxjIUJkRsxUtR0QbDGgxxxRUwJVQ8bjiiippUsppt28tqmS.fpYMqYpu9q+ZkRoT4latpi+3O9.1o456K+xubf0Leee84Jx9433nWmk9VVit1q8Z0GmXuddd5ySj16i9nOJv7tkkktMcbbzqEokVZpN24NqV5RWpRoTpd1ydpsc47G48xwmHQBUgEVn5ZtlqQc4W9kqZbiar1FE6T9qiiip6cu6pANvApFxPFhZsqcsAF+xZtLVDVwJVgpgMrgAV+E6nN0oNpS5jNIUUpRUzeVyZVyTSaZSS444o1zl1jpCcnCk53j+JqAx59dBx77y7LOSoVaMOOy70lyIx79vF1v1i6aBgPHDxdGTPLBgPHj8SHOfsHJwDlvDB7vxhvDh3BlBVYJ.x+3e7OBHDVxjIUdddpu5q9JE.T8qe8S2GhfBJkR8hu3KFPvKKKKceHO3djHQz8s49433ndtm64z86a+1usJszRKkOzuoXOR+UwJVQ07l27zBRHet43JUaK8zSWA.0rm8r0yc+8+9eOP+DVjkzRKs.1zMbC2fdtZnCcnkRXB43LEjyTfOYdTo1oXlKdwKVcsW60FX9zb8xwwQ03F2X0S9jOoZG6XGAVqxJqrTmwYbFADi7jO4SV8oe5mp788KkXVBicriUcDGwQTlBp77O+yqxO+70q4lmCjJl8rms5LOyyrLW6pbkqr5du26UkWd4oOlS6zNsTtlatFTPAEn788Crsv6uYaXaaqGK4jSNZgdkuuH+s3hKNv3YoKcopgNzgphDIRoNWHRjHJaaaUSZRSTu3K9hpbyM2.i84Lm4n5ae6afw6odpmp5y+7OOv5rzm6NBOm+LOyynOex11Vail8W34cyykt+6+9UEWbw+gDkiPHDBg7GCKkZ+PVDkPHDB4vbLqlixkZMCKtTsexqk7gUrXwBzllgNn4wK4MLI4tKsijWxLqrgp+WXiYlv1EayztLsS0+KzvRjHQorIkJXkSTB6SIY+CTR9jRp.ilu1brWbwEqSb5l4iK0+KugkpJbXQEUDRO8zQxjIQjHQBDhlxXcWUYDk4Kyjhu41MsUYaKZQKBKe4KGqcsqEqd0qFMnAM.YjQF3nO5iFctycVWoCMyeVtttvyyCqacqC4kWdnZUqZ3HOxiLvboX2R+HuOmbxAyYNyAYkUVXsqcsn5Uu5nwMtwnacqanN0oNAr6vqE.Hv4ERx7eoKcoXtyctX8qe8H2byE0pV0BsnEs.G+we7nd0qd5vXz11FyadyCtttnRUpRHRjH5v5KVrXvxxBUrhUrLyacx1k0uvyylmCH8mY9JK75lb7KXAK.qXEq.qe8qGqd0qF0pV0B0t10FMtwMFsqcsCUtxUV21.Hv4Xqe8qGadyaFMnAM.0t10VW8PMIU8cpPx+XlEw.yuKJmCjp1Jd73kpeIDBgPHGXgBhQHDBgreBUn7HlHLfHfiqqqVvESwlTF4lHGGGc0MTxkWxCb633Tl4ZHSAsDgcjGBurxQRhHF.Pmr9cbbzBnYJtAPIBRIIa9vsioXPhfAx3JYxjZa2be2UigT89TMeCf.y4l8oo.ag6yvhXYhLukJgRJrvBQEqXESoMlJAfjsC.svWhMmJAGMGWhHWl1XXgjBioHkx9H4nJybMlYeaJDX347xRHmsu8sipTkpTJQeME3yT7KeeeDKVrRMlk8yrpcZddtY+adrl8ozOliiv6Sp5WSgTKqyIER0wKykohTcti4Zm4ZRpN+jPHDBgruGJHFgPHDx9ABKjf5+kDtEOexDSQZLev4vOjcXODKbEEzb+MaCwSoLssvd0krulu1ruME.vz6e.BJvhzFgE.HrPTlBrIaWDxvbdRr+vdFWp7tqTMVS0XTDcPZSSuSyT7wvd5l43P9bQfSIoya5wPgGmEVXgnxUtxkRfTSDgNk4YosBOFLWaL6ewtBKDnoXXRaF1KFMm6BKTYpDdKbkvz731UBcFdtQDHS9dQXAgE6VVuJq02xRTSSAHKKwZS0Xd2g45cXwsDQ.MO+wTj1xxqGIDBgPHG3fkvFBgPHj8CXJr.vNCIQQ3hjISp8FE4g3kGfFXmB.X96VYJTf4CTaJTf4wIu1zNLqfcl8qoPClU8OoeiEKVoDfKUgaooXMhPDlBUHiOysKg6novBlB0X5cPBhPFxwXZWliWSgjhGOd.Q6D62LTMSKszBL+XJlgztRXvIsgHbmsssVjNSaWFmUtxUVaClsYwEWrVPJosEQlj1RD5xbrIuVDRRp1khsI1u4ZkqqqtcLE+QVuS0mIu2TTMwVcbbzh+X5YURH6ZNNk+VTQEoE8Tr6vhDKsg7cjHQhn8TLYalywx9aZ2x2AhDIhNjQCOtj4EyvLt7hofklqahPjlm+XZyl6uY+smz2DBgPHj8dnGhQHDBgb..Su1A.ADWHrmqX95T4cSBldhhov.gyCUg85rc0qSkmHY116pP6zzieBK7P3vWSvza0B6wSoZtwreMCkQSuOKUyekU3qlp4fTk6vBayldFloMuqVuDw9LGqox9L+bSuwqrlGB6YZgOeHUdNWY4kZl8Q3virrxWbkk2tEd9xb9LUy+lgqXpN2EnzgpX39vruLOezz6sJq750t572vXNNDaNbNKyb+BOmIymk04JDBgPHj8+PAwHDBgPHDBgPHDBgbXELjIIDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHDBgPHDBgbXETPLBgPHDBgPHDBgPHGVAEDiPHDBgPHDBgPHDxgUPAwHDBgPHG1gRo..fmmG..7880elrs8Fj1vrsL6CBgPHj8FjqiI+E.v00Mv6S00c1UetRoBbcQkRkx9YugTcMVwNbcc2mzG6o1fuuud7J1BuN9gGPAwHDBgPHG1hssM7880+E.vwwYutck1v22WeyzVVV5sQHDBgr2fHPkkkERlLIbccQjHQfkkkVXIaaaTbwECfTesGaa6.smkkEbbbzhCYYYAKKK344Enc2aPt9X73w0aSt9XjHQ1qa+xqM355pmSrssgssMrrrfssMTJEbbbPwEWr1d877BLeQ9yAVp8UR8RHDBgPH++QHBgATxM5J2zaxjIQznQ2qaaOOOc6TbwEiJTgJr2YvDBgPHFXJnCPPu3RDYRP79IQzmvelRozhAIGe31eeAlse73wQZokF.12bs2xa+aN2jLYR..DMZT355pEEzb+DwFMuWAxeNfRbRHDBgPNrD4W8VoTA9Uu2WbC4111HZzn5eA7nQipuoaBgPHj8FRjHA.J4ZMdddAByQQzqDIRDX6h2dArSAvJKeiQDCyzqnJpnh1mD1jVVV51QDCKQhDHZzn6SRYAkm9OQhDZw9hFMpVjqHQh.GGm.dDmuuu9yoXX+4C5gXDBgPHjC6vLLIM+ku2W4IWl+R2l8Q39iPHDBYOEQfKgvdClrO111HYxjvwwIkdEl4wJdDk3ETldH09Zu2x22GtttHVrXk40i2egoWdIBwIheIaWlGjvnLUGK4OGPAwHDBgPHG1g4CBH+J4RHNFNbJ9ihby8gaOIzKHDBgP1SQDmxT.IIY5KWCyyyCQhDoTWuIrHOlWex70IRj.whES6wVhHP6K9QiDgm1Uh3s+lxRrOY6g+grDuqi+nV+4CJHFgPHDB4vNja5MU+p2ETPAnxUtx6SZ+DIRfHQhn+0ukeQbBgPHj8FLy+Vh.VJkBKYIKAKdwKFEVXgnoMsonIMoIn10t1v11FyblyDst0sFomd561e3m824zKy7xUpxqY6uP7xqDIRf0u90irxJKrt0sNznF0HT+5Wezrl0L344gUspUgJUoJg5Uu5UpPNk7mGn7lDBgPHjCKY6ae63Idhm.CYHCACbfCDcu6cG0qd0C8t28dutssrrv67NuCdfG3Av0e8WO5ae6KZQKZARKszvBVvB1GX8DBgPNbDwqsjebEeeeDKVLLqYMKbpm5ohq3JtBr90ud366iwO9wiN24NiQNxQhu8a+VzktzETXgEpyokkUdDKmbxAO7C+v3m+4eVueRes2hzmhmq433fbxIGz0t1UbO2y8rW296Nrrrv5V25vPFxPP25V2vO8S+Drsswrl0rvUdkWIFxPFBVvBV.5ae6q950hXXzWh9yGTPLBgPHDxeJQBajx5FX888QqacqQaZSavu9q+J9oe5mPN4jC5W+5WoRrulODfmmWoJ87l8g755V25h1zl1fhJpH7se62hrxJKzpV0Jzzl1TchONbeHaK7CejpGBQ1GoxgkJLsyTYiluOUySk2sQHDB4.ClghuDBje5m9onycty3XNliASYJSA27MeyXfCbfX3Ce33W+0eE+7O+y3LNiy.VVVnfBJP6MYlEWF45H4kWd3u7W9K3cdm2Q6szl8k40aLuNj72vWez7ykB.fomVoTJr7kubLqYMK7oe5mlxBPizldddo75hguNnYaD954KcoKEm7IexX1yd1XFyXF3AevGDW4Udk3Nuy6DScpSEGwQbDnKcoKHqrxBEWbw5wubM5Tcs6v8kTzCBO1k6KwbLIGu7dy8OU8QpFSoZej1ZWcMay6QJU82t5dcL+rv1SpteiTMOkp6EortumTM2uu.l.KHDBgPH+oCITLj7dhbi+lkV9pUspgK7BuP.TRkk5pu5qF..8nG8PebopjyalPck7Ih4ubrbiymwYbF..3jO4SFidziFQiFEW3EdgH8zSOvwKsuzGl41D4uoJYH666q+U1MGaxMbGKVr.kIdyjkrztR+GNTPRUNVS1GFxHDBgbvC4ZX.k7+QOsoMMzu90Ob9m+4iQNxQpCiRYepV0pFFwHFAxJqrvrl0rPQEUjtcLuFlkkExImbvfG7fQaaaaw69tuKpQMpgtekDsuonXx0HLu9RjHQPhDIfiiSfqkZ9dYaJkBNNNnCcnCX7ie7noMsooLWdJUSSyi2z9E6RNtnQiFHupIG6pW8pw4dtmK1912N94e9mQcpScz1BPIdc28du2K1vF1.F4HGIxO+7APoy8mlgSpHTnzGRdVSRYBxbT3JUozdxwFIRj.EX.45vR3vJ+UNFYrJ1lssMbcc0smTr.REhsFd9zb9122WaiReX1exmIiMyJRprcYbImqjp0SY+j6CJ78+Hio8WgTK8PLBgPHDxe5v7FWMuIOyaHS9EbcccwF23FA.PFYjAZW6Zm9XLuILIGfI+07gRj1NUBXs10tVc9J6TNkSIvupZptIOyDbrzu999A9EXsrrzO3gruxMT533fXwho+ExSUEwRtY6v2.ZpFKRtcQto9T8KSSHDB4.C111nfBJ...acqaE+e+e+evxxBCdvCF.67ZFISlTm.6yHiLv8bO2CrrrvN1wNzsiT8IKpnhv1111vC8PODdjG4Qvi9nOJpV0pVf1B.At1n4ORhHtgqqqNWYZ9CRYdcH45nR+mHQBTgJTAzm9zGzzl1T.ryqCJdnkruBhfalh.IB8IHBiIGuuuOdrG6wvJW4JwUdkWIpUspk99.BKT2sdq2p95nEWbwADFx22W2thvXlh0IhgEKVrc40NM8RJQPIYesss0EEAQDLwq9jwlLuK+vWl4rTypbs44DBgGuh3Wg83O4XL2eY+TJERlLYf8wL+rZYYgHQhnO2QDSS1GGGGctYSVKCe+Zlio8WdmN8PLBgPHDxeJIrmMIH+xjxM6FMZT7ge3GBaaabMWy0TpJCobi6lhWE9WAMU8i7qn9q+5uBfRtIxtzktTpJLoYHBD9WGMbEsxyySayxMSJ2Lc31U9bosja5Lr8KsgY+Z5YXlayrMHDBgbvgJW4JCkRgwMtwgEsnEA.fS3DNAs2gYJTiHZxodpmJpScpiV3njISpEhH8zSGomd5Xjibjor+L8fY.DvakLC4NwKmjPwKUWyPDTQt1iHdS3JkokkU.wkL8THouEuwRoTA79ZSApj4hErfEfW60dM..z+92+T5wZx1N5i9nQ+6e+Qt4lafppYXOlxLOnkLYRsWTI1qLdRkGfG9ZzlW+WFqhfRR6JBHIqIhHjhGwG1apB+izYVHCDaPtmHY+M8jdy6AvrBkJsaXujS5GY7GIRjRIdnLWK26RxjIKkfmx9Z5YboRbt8EPAwHDBgPH+oC4lFCeS7.6z6pjaTcwKdwXtyct..nW8pWAppVgCaCyahStgdIDDk2aVd1A.9vO7CgRovfFzfPkpTkBjGVDaL7M3I2LuX+xMDF9F3C6cYxMQFN7Nj1zLQFa9K8GdtStY3vgQYpD9iPHDxANL+gRF8nGM..5d26Npe8qu95PQiFEIRjHPXyUspUMLnAMHjLYRcH0IjYlYhoO8oirxJKDMZTT25VWb1m8YiVzhVnuliTUK888wLlwLv27MeC1912NZaaaKt3K9hwG+weLN8S+zwQezGstcyN6rwG8QeDVyZVCZUqZE5cu6MV8pWMpQMpAZW6ZGRjHA98e+2wl27lQQEUDbbbPm6bm0BfHhKsxUtRL4IOYrzktTjWd4gFzfFfN1wNh9zm9n6Kwa3LuFOvNEsYBSXB.njq40vF1v.GmHbibc6nQihy+7OerrksL.DrxTNiYLCLwINQrwMtQTiZTCbTG0Qg92+9ipW8pqu1ehDIvJW4JwN1wNPAET.xHiLPKZQKvZVyZv68duGhFMJtrK6xPSaZSwF23FwpV0pfmmGJnfBP25V2PznQwXFyXPlYlI5RW5BN+y+70+PXae6aGSXBS.KbgKDabiaDG4QdjnIMoI3rO6yFUoJUIv4HdddXG6XG3UdkWAqbkqDMrgMDm1ocZnt0st3G9ge.CdvCNfGp466ioLkofYO6YiUrhUfV25Viy5rNKzvF1Ps.dtttHZznX7ie73W+0eE4me9nqcsq3TNkSAYjQFZwyxM2bwzl1zzhjJyeVVVHszRCQiFE8t28NffYScpSEyXFy.qe8qGUspUEMnAM.WzEcQnt0stADPbeJJBgPHDB4Ow344o788S46SjHg50e8WWkd5oqhFMpZsqcsJOOOkRoTtttAZGWWW8wE9yLaeyOO6ryVA.ksssZTiZToz9788UISlTEOdbUxjIU999kY6444o+WxjI0elXyJkRkLYRkmmmx00U2dled31122W455lxwjr8vGCgPHjCNH++wacqaUYaaqrrrT+k+xeIk6ifbMf7yOeUhDIBbsrG3Ad.04bNmiZZSaZpku7kqlwLlg5ltoaREKVL0a8VuUfi2yyS8lu4apthq3JTScpSUszktT0G7AefpIMoIJGGG07m+706alYlo5jNoSR8Ue0Wo9se62f3wnjD...B.IQTPTUe8W+0p9zm9npV0plZzidzJkRol27lm5ltoaR0ktzEE.Tmy4bNJkpjqiIWS5K+xuTU25VW0W9keoZYKaYpEsnEodhm3ITQhDQcK2xsnJpnhBLNMG6Ra366qN+y+7UokVZJaaaUwEWr9ZnlGi40OSjHgJd73phKtXkRoT4kWdp92+9qtwa7FUyadySsl0rF0Tm5TUmzIcRppTkpnVvBVfdrujkrD0PG5PUsnEsP433ntu669TKbgKT0291W08ce2mJszRS0gNzAkuuu5Ue0WUMfAL.UjHQT111pksrkotka4VTOwS7Dp1zl1nbbbTey27MJkRoVwJVgp0st0pQLhQnV3BWnZoKcopO9i+XU0qd0UcqacSssssM83MYxjpbyMWUaaaaUe3G9gpktzkpl0rlk55ttqSkd5oqt0a8VCr1toMsI0McS2jpAMnApu9q+Z0pW8pUOzC8Ppi3HNB0BVvBz2Kw1291UCZPCRcdm24olwLlgZdyadpgNzgpNkS4TTYmc1AV+SO8zUssssU8rO6yp93O9iUevG7ApV25Vq.f5tu66VuNsicrC04cdmm51tsaS8q+5upV1xVlZRSZRpt28tqpQMpgZdyad50l80PAwHDBgPH+oCQLHSBeyTxM.2u90u.2HtPXQlBuc4lCMu4aSwr788Uu268dJKKKU5omtZYKaYkRzov8iocYdi7l+MU6qmmm9l12U6WhDIRYeVdHUB0QHDB4.Odddp0rl0n.fB.pa8VuUU73wUJUo+AaJnfBzGiRsyqC355pl+7muJszRS0u90u.WWXSaZSJ.npe8qutcUpRtVXFYjgZMqYMAtlv2+8euB.pEtvEp210ccWm569tuSe8Seee05V25TUqZUSMlwLl.13BVvBT111p9129p+wXDNhi3HT0st0Ukat4FX7e9m+4q.fJ6ryNv3R5qvz111V87kvN1wN1EyxAayQO5QqrrrTuxq7JZannhJR8y+7Oq.f5u829aA5eeee0W+0esB.p63NtC0EcQWjZiabipO4S9DE.TW5kdoAryN24Nq.fZjibjpm64dNkuuupEsnEJaaa0zl1zTtttpa4VtEE.TyblyLv8k7O+m+SUjHQTidziNv8n7du26oF9vGtVzPeeeUhDITm+4e9pa5lto.BcN7gObE.TSYJSQO1u268dU111pG8QeT89d228cq.fZIKYIJkRoEMr8su8pdzidnERbJSYJpi+3Od0l27lUJkRUbwEqF+3GuB.pNzgNn13F2nRoJ4bp268dOUrXwTu5q9pAlSl7jmbola2WCCYRBgPHDxe5vLbIJpnhv7l27vrl0rv1111PFYjA5Se5CZdyaN1vF1.9jO4SPZokFN6y9rCzFgy6X..aaaaCyadyCyYNyAaaaaC0st0EmxobJn4Mu4.nDW9WB0Caaa78e+2CkRg1111hl0rloaOkQnXJscd4kGlyblClyblCxO+7QqZUqPu5UuPFYjARlLIV7hWLpPEp.pRUpBpW8pWfvxbxSdxXtyctnm8rmnCcnC.Xm4fDIDG+we7GwblybvEewWLpacqa.aHd73329seCyYNyAYmc1vwwAMu4MGsu8sGG8Qez62ptSDBgP1yQEJGOJ4ESybmkDd9UpRUB.6rxDZFV8QiFEwiGGKdwKVmT2srrPMqYMQiabiwu8a+FJrvB0s8l27lQd4kGxLyLQCZPC.PIW2q6cu63ptpqJPEcdYKaYXRSZRnW8pW.njv3q90u93oe5mFEVXg5wgDhil4wRID9bccQcqacwblybPt4lKpbkqrNWUcrG6whnQih0st0gl1zllxbKk404jDwujW0788QEqXECjJDD6TRV+l47JIuks7kub89TgJTATyZVSXYYoCuRYsQ9bfRRcBO5i9nn10t13bO2yEie7iGsoMsIPtBUpnmiXDi.yctyEtttXbiabXMqYMnacqa5vG0xxBqbkqDcpScRulcLGyw.WWWjSN4DXdbkqbkXFyXFkpPFbO2y8f24cdGc5aXwKdw3AdfG.csqcEcsqcE.kD1nacqaE9993HOxiDVVVXQKZQ3odpmBm4Ydln4Mu4564w00EWxkbIXXCaX3S+zOEWxkbIHd73XvCdvnl0rl..Xyady3Ftga.1113EewWD0st0Um24RO8zQhDIvxV1xBTIKOxi7HA.vu8a+VfbP19RnfXDBgPHj+Tgj2Q788wO9i+HF1vFFl9zmNNwS7DwodpmJxLyLwvF1vvXG6XQd4kGhDIBhGON5YO6YoJe7lIy1u3K9BbW20cgLxHCboW5khsrksfG4QdD355hwN1whS+zOc8MLaYYg7xKO7QezGA.fq3Jth.1nbCqxM3M24NW7+8+8+gLyLy.UZpt10thu7K+Rrt0sNzt10N..LwINQTu5UO.Thveu+6+93xu7KGNNNnN0oNHyLyD0u90Wm3bsssQd4kGN2y8bQgEVH777vMey2r9F9yM2bwi8XOFdm24cvcdm2IZZSaJJnfBzs667NuCFv.Ff1lMefKBgPHG3wwwAUrhUTKnwpW8pAvNE8JbdzR9+sMEKwyyCMu4MGKXAK.0t10Vm2wVyZVCV6ZWq9GbYG6XGnZUqZvwwA0nF0.MoIMAW3Edg3xu7KGcqacCsu8sGsnEs.u5q9pADa63O9iGO4S9jX5Se53bNmyAG2wcb33O9iGWy0bM5pSn7isXVQJEjbM029seK17l2LZTiZD..1zl1DxN6rwZVyZfqqafJQo5+UMGkJznYtvroMsoXIKYIvwwAacqaE0t10V2OgyUnlU4QoJGd1m8Yi4N24hi5nNJXaai7yOergMrALsoMMnTJTTQEE3dHLSJ8qacqCcsqcUKlWu6cu04ySaa6.UqwK7BuPjVZoAkRgF23FiF23FqW+dfG3Av.G3.QqacqgqqKxO+7wpW8pQlYlIrssQ73wCbdxwdrGKt+6+9Qm5TmP+6e+Q6ZW6PKaYKQaZSavK9huntOW3BWH..N4S9j0E0..fm3IdB7W+q+UzxV1RXYYg4N24BkRoueDIuqEIRDznF0HXaaiYLiYfALfAf1291iV1xVpO+6tu66FqZUqBO6y9rnyctyADg7bNmyAyctyEMtwMVaSKXAK.yblyD..EVXg62tuCJHFgPHDB4+uCSQYLuoZyeAwW9keYLjgLDTkpTE7QezGgK9huX8ul4nF0nvPFxPvIdhmHbccQqZUqPiZTiBT8mT+uDPukkEti63Nvy8bOGt268dwe+u+20dAVMqYMw+3e7OvkcYWFV+5Wu9AJhEKl1Su..5XG6XJs+nQihoLkofK9huXjSN4fy+7OeL7gObbrG6whhJpH7vO7Cim7IeRcRJNVrXnssss5a3N6ryF25sdqnG8nGXJSYJXSaZS38e+2G+s+1eSON.J4WFeaaaa..nfBJPei51113ltoaB+xu7K5eweozt2912dLlwLFbIWxkfi63NNzpV0JVkIIDB4P.777PMqYMwIcRmD9oe5mvu8a+VfpUr7CqX9Z.fYO6YiMu4Mid0qdout4wcbGGl9zmNdkW4UvW7EeAN4S9jQaZSazdTkY0jLVrX3+7e9O3Ftga.idziFe3G9g5Dz+vF1vvvF1vPEqXEA.vsca2FxJqrvDm3DwO7C+f116e+6OdgW3EP8pW8zEA.wCljq8JdmkqqKpV0pFhGONdkW4Uv6+9uOZTiZDZaaaKJt3h0WGCXmW+2LIta5I08rm8De0W8UvyyC4me9nl0rlAJBAx9IBUkSN4fO4S9DbsW60hjIShzSOcz7l2b78e+2iwMtwgoMsog9129h5Tm5TpJ.o7CRIuuEsnE3XO1iMP0iTVGAP.Oe5DOwSTKrlqqafq6VkpTEzfFz.7we7GiwMtwg0st0g9129hcricDnn9HqWWvEbAXnCcn3EdgW.Oxi7H542N24Nim4YdFsmmshUrBDMZTjQFY.fc9i1kQFYfV0pVo6+ryNaXaaq2Ooe..pQMpA788QlYlIrrrPMpQMzig25sdKLpQMJbFmwYfq4Ztl.yOx4UMu4MGScpSEe1m8YXNyYNnW8pWZuKSD8b+B62BFSBgPHDBY+HgST9ttt57jwHG4HU.PEIRD0m+4edfbsg7ZI+i..0i8XOltcLSX8999pa9luYkkkk5FuwaTEOd7.8qjOLbbbTYkUVAru65ttKkkkk5nO5iNP96P5GeeeU1YmspwMtwJ.7+i8tuCOppx+ef+9duSIID58dIDDhBJUwFpnTUYATQjlrhEPTWXQwxJfJHBr6ph+T5HOfthJxBJhJJnnRKBhqzAMzRnEBFHsobKme+Q1ygyLI3pO52cAy6WOO9jjYt86r6L7d9b9bD2wcbGhPgBoV1HQhHJnfBD0qd0Sz111VgkkkXfCbfwzePFwHFgXQKZQhPgBI5cu6s..hl27lq56Ixlmrqqqnd0qdB.HV6ZWq5XX4Ke4BKKKwTm5TKw02QO5QKLMME986Wr7ku7Rb8iHhn+2ZNyYNp2C5Tm5TBgH19dkNaaawjm7jEKaYKSHDE+d.4lath669tOQ0qd0EKYIKQTTQEoV+l0rlI.f3XG6Xp2yJRjHBGGGQ94mu3S+zOULsoMMQW6ZWUuW5e+u+2EBQwu+k78x1zl1j3UdkWQLnAMH0wZO5QOTGWQhDQbnCcHA.D8pW8RHDw1CMkue9jlzjDm5TmR8dfie7iW..wm7IeRLmi5+TxyySjd5oK74ym..h0u90Gyyo+dyxy+O7C+PwS7DOg5499u+6EctycVz912dwF23FUq+d26dEFFFhN24NKhFMZL8+LYuUqu8suk5DZibe555J5V25lvvvPjd5oWh6cxqkexm7IhF1vFJF3.Gnp2oEJTHw68dumHgDRPL4IO4XVuhJpHgmmmXm6bmhEtvEJFwHFgnF0nFp9CW1YmsvwwQLkoLEA.DiabiKliS89Npqqq3kdoWR..wy+7OeI5UcqZUqR8YYzWmcricHRJojDUspUUr28t2RsWpdfCb.wUe0Ws3xu7KWr4MuY0w8AO3AE.PbcW20Ey0reKwFBAQDQDcAK8uQbY+5Xiabi3AevGD986GidziF8pW8R8sQJKseOOOzpV0J014xu7KWMDGzGpDu268dXVyZVPHDXjibjpovdYOGYm6bmp9vh921cnPgvG8QeD74yG5Uu5kpOhne7ZXXfYO6YiLyLSjPBIfm64dN3ymO0xFHP.DLXPzhVzBr0stU355hd0qdoFdEG5PGBomd53Nti6...2wcbGp9XRVYkk5asGn3u47yblyf1111ppJNSSSrm8rG0vKQtbxiypW8pCOOOXaai10t1odN4xRDQz+6HDBLvANPbEWwU.WWW7Ye1mESeuBn3dnoTQEUDV7hWLZW6ZG.J98g96+8+NlyblCdgW3EPe5SePhIlnZHyIqpX46K533fssssgQNxQhjSNYbi23MhG4QdD7we7Gi8su8gdzidfUspUo1usoMsA4me9ncsqcXjibj30e8WGm3Dm.O7C+v3i+3OFm7jmD.E+dcx2ewwwIlJq5C9fO.ibjiDCe3CG+k+xeAUpRURc9blybF.fXd+K4uK+YznQAPwU7TqacqwHFwH..vxV1x.vYawB111p2Wzue+v11FqXEqP0JEhFMJF7fGL9rO6yvblybP6ae6A.hoJm.JtRujUgk7yTXYYgFzfFDyvxTRVg65Cw0xW9xq11RIjPB3.G3.nO8oOvyyCyd1yFMoIMQ8bETPAv11NliEOOObm24ch8rm8fzRKMLjgLDLiYLCr+8ueLu4MOb7iebr8sucXYYglzjlf.ABfzSOc0qsjmO..m5TmBlllnEsnEvxxR0qQk2yDBAxJqrfooIZYKao53111FO5i9nnnhJBSe5SGMqYMScdJGJnNNNXPCZPX8qe8X1yd1ncsqcHwDST05Ije9M48weqw.wHhHhnK3n+gPk+mggAhFMJdlm4Y.PwePL4G9E.wLzJLLLPAET...RJojvUcUWk5Co444AWWWDMZT7nO5iBGGG7fO3ChK9hu3X5MHG3.G.O+y+7vwwASaZSCImbxpOrVFYjA14N2IrsswsbK2hJ.M49WHDXu6cuX5Se5v00EO4S9jHkTRQ8A+bbb.Pwen9LyLS.T7+nj1111pFJmm4LmAScpSEABD.IjPB3RuzKE.E+gf27l2r57Q1XcyO+7wC8PODBDH.74yGhFMJxKu7..vq7JuB95u9qU+CgLMMwC8POD1111F1291GpV0p1+mzLaIhH5WFYehR1isdlm4YPBIj.dgW3ETewO.E+9LIlXh.n3PKlwLlA5e+6OpcsqsZasjkrD..zjlzDUCmG.XqacqH6ryVst.EGNRgEVHdq25svgO7gQ3vgUCQuF0nFgQLhQfDRHAXaaCKKKr6cuargMrgXBNoF0nFX3Ce3vxxJlfNjAWYYYAKKK0W7yl1zlfggAZdyaN.Na3b4kWdXkqbkpqG5eAPRw+kB42ueL1wNVT6ZWa729a+MUSvG.p2WTdb9ke4WhbxIGz4N2Y..jYlYhsrks.KKKjRJoDS6UXMqYMvxxBQhDQcdHuFHO+RLwDUWGkgMIGRq.E+91x6qxIAg3aOAaYKaAgBEBst0sFkqbkKlf7V6ZWqZxBPeB6wwwAKcoKUMIF.T7m44tu66FUspUEIlXhv00EsqcsC974CqYMqA6bm6T8Yc..99u+6wnG8nQznQQaaaaQZokF97O+yQ3vgi4KTacqacvvv.21scapOGwjm7jwG8QeDF1vFFFv.F.bccgPHvd1ydTCq08u+8iMtwMBCCCT+5We0wommmZxEP94x9+BLPLhHhH5BNxOnndnXdddXSaZSXUqZUPHDXvCdvngMrgknupHWusrks..fAMnAgfAChPgBo9.pVVVXUqZUHiLx...cqacC.P0DX2912N5W+5GJpnhvS7DOAt268dUUSkooI13F2nZVoRNSRo2aOLLLvxV1xPznQgggA5ZW6p5CIK+GJDJTHjWd4g8rm8.e97g1zl1fTRIEUHWspUsB23MdiPHDnnhJBolZppdM1t28tU6GWWWroMsIzpV0Jzu90O00i.ABfa4VtEXXXfCbfCfN1wNhALfAfO7C+PTXgEhjRJIzhVzBjZpohfACpVO4OIhH5+9jMbcfhe+ha3FtALyYNSrwMtQbK2xsnZv9x2uIu7xCuvK7BXSaZSXricrwDTirhedm24cPN4jCrrrvl27lwi9nOppwuul0rF7Mey2fhJpHjPBIfyblyfQLhQfbxIG..UnHqXEq.29se6wzOtt268dwW+0esZYDBAV25VGF5PGppQ8655hvgCCSSS78e+2qBFxue+nQMpQPHD3Mdi2.YkUV..X+6e+3ttq6BWwUbE.n3YY4CbfCficrignQiBOOuX5gXxvYbbbP8pW8vxW9xQCaXCQ25V2vl27lUUmFPwAY89u+6iG9geX7W+q+U.T7mcnV0pVHgDR.tttXIKYIpi424cdGr5UuZTqZUKr0stU7ce22gsssso56Z4kWdv00Em3DmnTqtZ8IDfbxIG3ymObxSdRUElIO1hDIhpgy+9u+6isrks.GGGDJTH7rO6yptlst0sNr6cua7C+vOn1tSXBS.qXEqHlyycsqcg5Tm5n5OXMtwMFie7iGFFFXXCaXX8qe8HZznXO6YO3du26EO5i9nvxxBUspUE+8+9eGG7fGDu9q+5pOK1V25VwBVvBvTlxTPpolJLLLvm+4eNl5TmJZZSaJlxTlh5dxO7C+.F6XGKpacqKbbbPcpScTUE268duGhFMJJpnhvxV1xvG7Ae.pUspE1912tpWm9ateyGDlDQDQD8eA588CoG4QdDUOBYkqbkBg3r8dCcYkUVB.HBFLnXIKYIp9MRjHQDBQw8ujgLjgHLMME.Pb3CeXQznQE6XG6PLgILAUOsXkqbkhvgCGy1NZznhdzidHLMME8t28Nl9Ql7mEUTQhpUspI.f3xtrKSstxdqgrGorqcsK0wvzm9ziY+XaaWhqAiXDiPXZZJtwa7FU8UMWWWwUdkWoXQKZQwrrxy0O9i+XQ0qd0E.PsupQMpgX9ye9p8i9OkGaDQD8+Ntttw79Oae6aWL7gOb06qLlwLFQu5UuDojRJh4O+4KrssEgBERXaaqdulCdvCJt7K+xE.PzfFz.QpolpXDiXDhyblyH1111lnksrkhZVyZJd7G+wEdddhO+y+bw0dsWqXkqbkhdzidHF4HGo34e9mWbS2zMIdtm64DQiFU355JN4IOont0sthUrhUHFxPFhXHCYHhoN0oJF7fGr3dtm6Qb7iebgqqqHu7xSbG2wcHZPCZfpued4W9kKdwW7EEBgPbzidTw8bO2ivvvPXYYIZUqZknacqahctycJJnfBDcu6cW..w.Fv.D4latkZe5pz5mXYmc1hm4YdFQ0pV0DolZphgNzgJF7fGrHkTRQ7LOyyHxM2bUuGo787du268DMpQMR..QKaYKE0st0U729a+MgPHDu4a9lhZVyZJRM0TEu4a9lhrxJKQe6aeEUspUUXYYI.f3Zu1qsD8qTgPHl0rlkncsqcp2+sIMoIha61tMUOBS+bZRSZRBSSSgOe9DssssUzzl1Tw68dumvyyS7XO1iI.fnScpShMu4MKDBgne8qehErfEHdpm5oDcu6cW7rO6yJd7G+wEWwUbEhu8a+1XdsjPHDKXAKPTkpTE08hl0rlIVwJVgZ4jWS9pu5qDsrksTbO2y8HF6XGqnEsnEh29sea0m2JZznha4VtEA.DUrhUTzvF1PQyadyEFFFp9MW26d2UmeqXEqPTm5TGA.DWxkbIhFzfFHdgW3EDQiFU7tu66JpQMpgnQMpQhktzk9y9+8wOWFBg1baJQDQDQWfPnM8lKq9pzRKMrm8rG..b3CeXT25VW02LrdePYYKaY3Nuy6DtttHyLyD0qd0C.msehHzlNvSJojvS9jOIV0pVEN1wNF5ae6K5cu6MZcqaMBFLn5agWNyQd7iebTu5UO355hEtvEhgLjgn11xY2wcricnlsHG23FGF23FWLyZlxs467NuC5e+6OLMMwW8UeEtxq7JiY1CSNkvK2+u669tn+8u+.n39qREpPEvG7Ae.l3DmH9pu5qT8.McNNN3fG7f30dsWCuzK8Rv11VU8.+s+1eCiZTiRMrNzmcOIhH5+9juWh9Lrr9+ey4kWd3HG4HnvBKD0oN0A0pV0pDuegbcjQAru8sOHDBzvF1PjbxIq1ddddnfBJ.UnBUPM6OFIRDjbxIinQihLyLSTXgEhpW8piZW6ZqdeY.fbxIGT0pVU344gie7iiryNaT8pWcTiZTCUETAfXVmR68XDBA1+92OxO+7Q8pW8P0pV0TOtbXBFLXPDIRDDLXvexsa7yDiNNNpsc0pV0PspUsPf.AhYlZTxvv.m9zmFYlYlHPf.p1bfblkNb3vHPf.wzZDjqe7W+k+t9LicoM6YWZyh1Ymc13fG7fnt0stnt0stwbbVXgEFy8uSdxShpW8pCfh6CXYkUVHojRBojRJkX6JudUXgEhicrighJpHzjlzD0vyTdLqecMiLx.111nwMtwp9ul9yqe9o+6.PccSd7GIRDru8sODLXPjRJonpvOaaaHDhXZoC+VhAhQDQDQWvQ+CbAT7GzJRjHpxtGn3gXf7CaoGzkooIt8a+1w69tuK5V25F93O9ii4Cl444gbxIGTm5TG344gN1wNhQMpQgV0pVg5Tm5n9GFnu+0+PqKaYKC8su8E.E26MjePQYy32zzDqYMqAcsqcEdddXQKZQXfCbfpskbncZYYg9zm9fku7kiZVyZhCcnCoFNk56S8881291UCQyzSOczt10NzktzELlwLFzidziXt9433nlJ4kmKm3Dm.qbkqDO4S9j3G+weD1113XG6XpIL.8ObOQDQDQWHieZFhHhH5BNxvkjeC25eyg974CcnCcPEvkrudICN569tuCu669tHXvfnqcsqv00EABD.QiFEie7iGm3Dm.ABDPss6PG5.t0a8VQyadyKQXXxJoxzzT0CO1vF1.LMMwkdoWJZPCZf5a4blybl369tuS0v9kgJU8pWcU0qIeLYCId4Ke4HXvfn+8u+HXvfvmOeXcqacXgKbgkZnTxuwZgPfLxHCr3EuXDHP.zidzCUC8011FaaaaCojRJnScpS3.G3.pi8pUspg+3e7Ohu5q9JU+FQ1OZjmmLLLhHhH52C3mngHhHhtfjrpqjytjACFDsnEs.NNN36+9uWMcyKmsprrrP3vgwbm6bAPwUPVm5TmTM02MrgMfu669NTyZVSToJUIUyFNqrxRst.mMDLY.baZSaRMjMbccwBVvBfmmG5e+6uZ5j+q+5uFKXAK.MrgMDVVVn10t1pgwYgEVnpxqhFMpZVp5se62F986GQiFEW20cc.n3Fv6K+xurZlXRtMjgClTRIopDrO5i9HLgILALgILAUCM1yyC986GKcoKEG4HGAe629sXUqZUpgYhb3ZjZpohJW4Ji.ABfpTkpnN1.XS0mHhHh98AFHFQDQDcAIWWWU0J42ueDHP.zyd1SDLXPjat4hcric.fhC4Q12Jdy27MQ5omtppsRIkTfiiCxLyLwfG7fwy+7OupWYHmhvW8pWMNzgNDRHgD..TAP444gYO6Yi6+9ueb5SeZ..bzidTblybFXYYg1zl1..ficrigALfAfoO8oiJUoJAfhqjqFzfF.OOO7Ye1mo5SFABD.ABD.qcsqEe9m+4pdmwkbIWB..V3BWHpZUqJt5q9pA.TCeRYOOyyyS8bKZQKB8nG8.sqcsKlYkSfhmV2MMMQEqXEQu6cuUCgR4xjYlYh7xKObEWwUfF1vFBGGG0LaFqPLhHhH52C3mngHhHhtfirGaIqVK4PnbLiYLnBUnB..XzidzXiabiHTnP3HG4H34dtmCyYNyAKXAKPMbH2zl1DV25VGFv.F.l6bmKRKszfPHfOe9v841mxWC..f.PRDEDUdu2KRKszvYNyYvC7.O.9lu4aPjHQPQEUD17l2LFxPFBV3BWHVwJVApbkqr53xwwQMrHW6ZWK5V25Fl7jmLthq3JTARU9xWdLwINQ..L6YOaL8oOcbjibD78e+2i4N24hwLlwfW9keYXZZBKKKjd5oi4Mu4g4N24hwMtwUhouc8.uRIkT..PUpRUvi8XOVLAXIuN0yd1SDHP.bYW1ko5gXxlCblYlId7G+wQUpRUvLlwLfe+9UWqiug7SDQDQzEpXS0mHhHhtfi9rwjrmdImAi1yd1CF6XGKVwJVQLqyHFwHvDlvDP0pV0v.Fv.v69tuK777PyZVyvK+xuL5Tm5DBFLnpBvDBAN1wNFF+3GOV3BWnJLIgPfZW6ZiG7AePL7gOb0PJTt+6ZW6JVyZVipg7OoIMIbC2vMDyLLkbXGtnEsHb+2+8Caaa0LNUe6aew+u+e++PUpRUPu6cuwpV0pfggA9C+g+.l4LmIpUspULWGjy3jxIFfcu6ciK4RtD7Juxqf6+9ue0rakmmGBDHfZl4ZSaZSXbiabXW6ZW3lu4aFUu5UGYjQF3i9nOBCcnCEOxi7HngMrg..pIj.8eRDQDQzExXfXDQDQzErzmcHAfZZW211FG7fGDG4HGABg.0u90GMsoMUsb4kWd3a+1uEABD.okVZnhUrhpfdzCaSFx0gO7gwANvAPznQQkpTkPiabiUS86NNNvvvPUEUgCGFaXCa.UrhUDMsoMUEXlj9zstOe9vYNyYvd1ydPt4lKpcsqMZUqZkpptxImbvt10tPsqcsQCaXCUM6e4viTeFoD.pdX1a8VuEV4JWIBFLXLAwEMZT32u+XlTBN7gOLxLyLQAET.pcsqMpe8qOpQMpQLWai+5LQDQDQWniAhQDQDQWPRFziLXH4eKCHRHDwzyqj88q3Gtg5OtL.L8sujLvL49Stt5gRou9wudkVnRxmK9pGSNYAHI2tk19VtcDBAxLyLQG6XGwG8QeDtjK4RhoGi444oBNS+bO9dBl9xpWAZ5AuQDQDQzE5XODiHhHhtfjLrGYfNx+VueZoG1iooYIBCK9GWOLq3WVY.UwuM02mwGFl95UZUXk74z2WxYOS4vbTe6C.bhSbB7Zu1qg+0+5eEy1onhJBiYLiASXBS.W1kcYwrMkUvVoctWZGSxk0xxpDmiDQD8+cJsYxW4v0W+4j8PS8eWNa.quNxmW+2AN6WjR7jSvJmqiE45HG98k19y11tTWd41Se8heeE+youcieaF+4l94j91I9qc5WmJssq90.4equO0mHZh+w+k3mZVaVeaoe7HOuzuOF+xD+4q9xIu2H+x.OW6S4e+ScNI2OmqWWJWe8mW+dU7qSocu9+qwJDiHhHhnyCoWQVxJHywwA2y8bOXgKbgnQMpQXW6ZWHgDR.gCGFO9i+3vxxBSYJSgCuQhH5BX5SbLw2uLk+7vG9vXW6ZWHmbxAUoJUAMoIMAMu4MGBg.6XG6.MnAM.UrhUDG+3GGG5PGBm4LmA4me9nJUoJnScpSw7kFICtH9pS9bcbo6mpuR533.SSyX9xWxLyLwgO7gQt4lKJpnhPpolJtzK8RgooYLaq3qx632WwWwz.PsMbcciohwkrssi4K54nG8nHZznHb3vvxxBACFTUk1974C0t10F.EOaQGJTHUOFMPf.HZznHXvfn90u9p10vOWYjQFpuzJYPVojRJp8srJs877P1YmMJpnhT6aYEaW9xWdTiZTC04Uo8ExIulAfRbOJ91sfrx5ksNB4DLjbaTZUYdosujegd.PMidKet3uepeeSmmmGrssUWSieDA7aIVgXDQDQDcdFgPDyv9TuZ3BGNL..N9wONN8oOM10t1EdfG3APf.Avy8bOGCCiHht.mLLj3Cfv11FG+3GG228cenCcnCX8qe8PHDXO6YOXnCcn3ge3GFae6aG21scaXKaYKPHD3zm9z3S9jOASdxSF21sca3i9nOB.nDsS.YXR.mMDEOOuXpNJ8v4jzqPH4OkOuLLLgPfnQiBOOOr+8ue7Ye1mgG+webLnAMHr0stUUHJx.+jAqD+wnd0Foe7J+c82yTOvG8psVFnjqqK1xV1BF0nFEt4a9lQyZVyPiZTivEcQWDtu669vq+5uNLLLPjHQvnF0nPW5RWPyZVyPpolJZZSaJRM0Tw3G+3gmmmpec9ygPHva8VuE9i+w+HZTiZDZZSaJRKszvW9keYIpHaWWWLkoLEzt10Nz3F2XjRJofF23Fi1291iu7K+R04lLLrRqVmzC4zyyCyctyEsnEs.YjQFwDbkLLNYng58dT8JGWu50kaSY.Zx8ib80qVO48SgP.aaa0q4jmyxsqqqKLMMQvfAUUNlba8acXXxKZDQDQDQmGxwwoDO1V1xVDst0sV..A.DssssUL+4OegsssPHDBOOu+aeXRDQzuQ777Dtttw7XQiFUHDBw2+8eunoMsoh1zl1HxLyLiYYxKu7DSbhST8dCqbkqTschDIh3S9jOQXZZJF0nFkPHJ98Whe+7S83dddhnQiJ777DNNN+rduF46KUZO1BVvBD.PLm4LmRrLtttk56+UZGqk1wg73OTnPpGKb3vBg3rWKkquPHD6YO6Q32ue00M8kW5HG4HB.HLLLDyZVyRbzidT04h917mCGGGgiii3K9huPc+5RtjKQjc1YKDhhuFoedc7iebQ+5W+D.PLwINQQ3vgU2ODhhu2DIRD049O08la+1ucgkkk3sdq2pDmm522k+94ZaIuG8Sc8O96g+T2S022ETPA+rVmeKvJDiHhHhnyC433DyrcIPweSrssssEaXCa.aaaaCYjQFX8qe83tu66NldmFQDQWXRehLA3ryNvm7jmDctycFm3Dm.u268dnd0qdwT8Vku7kGO4S9jXzidz..pgYGPw8vxxUtxEyPNSuZojjUklb34Izp3HCCC0rTrdELoutRxJORVYOQhDQ8bxGSVMyxgnmqqqZeZZZFS0c433Dy9R9dh5CIO8iAYkHkPBIndNYUPIqlJ8Jvq7ku7pJiJwDSD..ACFLlJR6cdm2AMu4MGabiaD2+8e+n10t1wz+QOW8kr3EMZTXYYAKKKzfFz.XZZhpUspg8t28hG9geX0PyTdcA.nl0rlXHCYHvzzDCYHCQMTB862upx5BDHf59a7WSzecxTm5TwJW4JQu5UuTC+SgVk1IOOzGNs.PcOx11Nl6Qxp5Se+Ee+OUuBxhu22Iesh70b..kqbkC.kd0n8aMFHFQDQDQmGR+Cip+A3A.RHgDPKaYKQCZPCPf.ATCGAhHh98AYujRFXwi9nOJN9wONt+6+9QcpScTKibVL111Flll3dtm6AlllH2byMlFotbHwEe3R5u2gdfE5AqHGJa.mMrqRaxZQt+zmTWjCoP8m211F111pdioLfny09zmOewruz6WV5Mve44m9vWLRjHpYRZ8PWh+Zf7XVNiSqe8XYKaYXNyYNX4Ke43xu7KGgCGVsLxgB5O21UP7KmmmGVzhVDrsswRVxRv7m+7UAUpG.TspUsfmmGpYMqo5do73UOHv3CPTNwAICtpwMtwnacqaHwDSTc8VuueIO9jgApGfpkkkJTTgVC2W9ZUYOiSdsTdNHCNS+9q3e2ZHh+3WexBvmOevue+psy+WLjIYfXDQDQDcdF429p9+f.4GnWeYjeHW8Oj3O2ukZhHhN+i9rHn7+u8+0+5egEtvEBaaabS2zMUhdMUf.Afe+9gPHPZokFtq65tPQEUjJnC4yIDBUEPoWoOR5gboWcT5AIE+LvrLvI8YHY.DSk+355p5EUNNNHPf.pfWjUElLDG88odCvWtMiuWc42ueUSXWFrkd+qRFFW7ydkw2iqjueqL3GoktzkhG+web7du26gTRIE.T7WJkb+52u+RTMU+bnWoXcu6cGicriEttt3gdnGBomd5wTkb111nJUoJHojRRErm7ZmdPQ5ARIqHN8yaGGG0mQPu2eoe+RVIX5AkpG5nb8huZsj6W41I9Iv.8PFi+9gdPoVVVppdStLVVVHTnP+ru19KQIm6wIhHhHh9eJ8YRL42dp7CQJ+V0ie4APLCCEhHhtvikkULyZgQiFEqcsqUMqSdwW7EWhYgR82i..nacqaX+6e+wzf4kumRnPgfkkE9lu4avG7Ae.JrvBQG5PGPW5RWPUpRUhoYzKCyZm6bmXsqcsX+6e+vzzDMnAM.ctycFspUsJliyCdvChvgCiryNaXXXfa3FtAjat4h23MdCb7iebbS2zMgq5ptJ.b1vsBDHf585zCT5nG8n3K+xuDe629sH4jSFW20cc3JuxqTMgyHmIIyO+7wG7Ae.RO8zQRIkDZaaaKt1q8ZwLlwLTM8d8ltu7bReHYpOCFJChzyyCKcoKEO6y9rXYKaYH0TSE.wNDB0Ol+kTAStttHPf.p6gddd3oe5mFey27MXMqYMXXCaX3y9rOCUqZUC.mcHYV9xW9XNt0ChasqcsXaaaa3fG7fn8su83ptpqBMrgMTcu4jm7jHmbxAgBEB4kWdHszRC0qd0Cttt3nG8n3zm9zH+7yGEUTQ3xtrKCNNN3Mdi2.YjQFnksrknW8pWnt0stvzzDqd0qFqYMqA..Wy0bMnScpSH4jSVU0WxiY.fst0shu3K9BbzidTz7l2bb8W+0ilzjlnt1ICg0wwA6ZW6BqYMqAYkUVH4jSVcd7JuxqfILgIfgMrgg4Mu48y5Z7OWrBwHhHhH57P5yLT5Cy.4vUnzFhjbXSRDQWXSHDpPl.JNfr0u90qpZmJVwJpBcPFNDPr+++2291W7POzCEyWVRBIjfJfr4Mu4gQO5QiJW4JiCe3CiALfAfgO7giBJnfXBXyyyCScpSEiYLiAsqcsC+4+7eFCYHCA6bm6DW5kdo3Mdi2..mMbnW8UeU7m+y+YzktzE7m9S+Ib5SeZbm24cBWWW77O+yiq4ZtFb7ie7XBNI9groggARO8zQm6bmwt28twvF1vvMbC2.dhm3IvDlvDP3vgiouUMpQMJjQFYfgNzgh92+9iCcnCgK5htHLgILgX5aU5uuodEfICKSd9JCNa4Ke4ne8qevvv.0pV0JlJbReHVJqhI8q0+mnOqbJ2uIlXh3Ue0WE0nF0.6bm6DO5i9nwTQb5ULkdUjWXgEhm3IdBLpQMJjVZogQLhQfSdxShd1ydhMsoMAfhqfvEu3EioO8oit28tit0stgLyLSUvfqd0qFSYJSAW8Ue0nKcoK3q9puB20ccWvwwAImbxXjibjnm8rmH6ryFSYJSAKZQKBUoJUAaaaaC8pW8BSaZSSctn2S3l9zmN5QO5AZdyaNt+6+9QjHQv0bMWC97O+yUWykWul1zlFt5q9pQUpRUvnF0nPu5UuvTlxTP0pV0PqZUqvce22M9jO4SPjHQJQ09oO7M+E6WaW4mHhHhHhHhH5WO8Y2O4u2xV1R0rQnsssZlDTNC+o+6x0QtLxYouUu5UKrrrDokVZh+ze5OIBGNrZYu669tElllhW+0ec01wwwQr28tWA.D8qe8SM6S544IxImbDVVVhF23FKJpnhh43+nG8nB.Ht3K9hESYJSQrgMrAwO7C+f..hDRHAwO9i+nvyyS7lu4aJLLLDyadySMaE555JxN6rEUqZUSL3AOX0wtqqq3a9luQ..wa7Fug5Xbaaaah1zl1nlgEkG2u3K9hB.nlAE0uNUZN5QOpvvvP..wZW6ZE+y+4+Tz7l2bQkqbkE.PL7gOb01tz1FxY6wety9lxiE4LWo79PznQEKcoKUcudNyYNpYyw8su8IpcsqcINel5Tmp..hzSO8X118u+8Wzl1zFwoO8oiYFb79u+6W..wW8Uek5XVdM9ltoaRXYYIZSaZi3XG6XpmaxSdxBSSSQ26d2ESbhSTssNyYNinAMnAB.H18t2cLyvkKe4KWXXXHd228ci4ZynF0nDIlXhhCcnCo19x6sSdxSVHDB0qwSO8zE.PLsoMMQ94mu3Tm5Tk50X8W++KEqPLhHhHhHhHhNOfrhkzGpew+7xpFVVoQh+c01HedYujRNbzjCQOWWWrqcsK7bO2yEyvq+Vu0aEddd3q+5uVU0Nx8ie+9wN1wNhYFErJUoJnN0oN3.G3.Hb3v.3r8gKe97g.ABfibjifCdvChK+xubzjlzDroMsI7Mey2fxW9xGS0XIOGj8CyoO8oiScpSgALfAnN1MLLPqacqgggAl7jmLxO+7A.vINwIvN1wNvwN1wTGC..CYHCAMpQMRUocx9S1Ompn9se62FO9i+3XoKco3ke4WF974CyZVyBKe4KWcMH9d0odqMH98Q72+j+c7U+krej0qd0KLtwMN..be228gzSOcUk.J2ux8wgNzgvy7LOCtga3FP6ae6U8hMfhG1racqaEqbkqLl93UMpQMJwjmfbalXhIBWWW7m+y+YTqZUKUkz0st0M344gO9i+XbW20coNGpPEp.5e+6OLMMUCmVOOOTPAEfwN1wBfhGRkxyaCCC7G9C+ADJTH77O+yq19qe8qGFFF3xtrKSsrxglqggAd629sQxImLpRUpRIlYQEkRER9KACDiHhHhHhHhnySDeHJMsoME.EGdQAETPIBAPe35IG9hRxYeR4xzm9zGUi0WF3VMqYMQf.APQEUTLA6zzl1TrsssMr5UuZ0PULyLyDqacqS03ykyxhx.zbbbfssMJrvBwMdi2nZ+d4W9kiVzhVnN1jmCxvO.JNjnktzkBgPfF0nFoZ77FFF3nG8nH4jSF6ZW6B4jSNHRjHnt0stvyyCWxkbI3O8m9SXIKYIXG6XGvmOeXu6cupPj762upOU8exJW4Jw67NuCRKszv.G3.Q+6e+gggAF9vGNxJqrfsssJLQ8dyl77++D81ef7d.vYGFk974CidziFWy0bMvvv.CaXCSMLSk6W4DhvN1wNPQEUDRKszTgUJar8ETPAvmOeH8zSWcs1ue+3Lm4LwbMWeV8zue+vmOepIO.4rWZRIkD.J90g0u90OlWuU6ZWa344g7xKO04SlYlIxHiLPspUsP0qd0U6GGGGTiZTC00YWWWXaaiDRHA00S40.IgPfJUoJUpCIxesggAv.wHhHhHhHhH57BxvJjUMT7UNSt4lqphmhO.FYOg5Tm5TXFyXF.3r8cRYyzuN0oNwTISx8iL7o3Cn4htnKBG6XGCuvK7Bn0st0XZSaZ36+9uWElld0n444g.AB.gP.aaazgNzgXpjs36WXwWsa4latXO6YOvmOe31u8aGcricDst0sFssssE29se6nYMqYn8su833G+3HXvfnEsnEXdyadvwwAu7K+xXfCbf3xtrKCMu4MGKdwKFFFFHPf.p9+U7SHMkl4Lm4fK8RuT0w3y+7OOpe8qONwINAdjG4QTANI928AMYPT5UmktepvZzqTN4wliiCpTkpDl8rmMRJojvt28twi7HOBDBABGNLDBg5dZFYjA..d+2+8Q6ae6QaaaaQqacqQ6ZW6v7m+7QaZSavIO4Igggg5ZPRIkjpxrjy9m.E2W0rssU8MLY.b986G986GACFDojRJkXVozxxBIjPBHPf.pWu9C+vOnl.HjyrkVVVvmOenpUsp..33G+33Dm3DvxxBW0UcUHojRBqd0qVEvniiCV25VG..5e+6eoVocmq9o5uDbVljHhHhHhHhnyCH+G5KCbA.nG8nGXpScpnnhJBYkUVngMrgwDPfrpvjUA0+5e8uvAO3AgmmGbccge+9QBIjPIpnLIYShWeFN1zzD4lat3u7W9K30e8WGKbgKDaZSaB986GlllXpScppgjo91QdbYZZhxW9xqNOzGFn5UXjdU9nWoXKe4KWUoRxYcSYEDICiIZzn31u8aG8rm8Domd5XSaZS38e+2Gae6aGCcnCE0u90GctyctTGZomKIlXhv11F986GBg.0qd0CyXFy.27MeyXIKYInKcoK3tu66NlFwuLbP8YNweN2eiu4vqWceolZpXVyZVXvCdvXwKdwpgupb3NZYYg7yOeXYYgALfAfIO4IWpyzkxyc4jkfLLU48E86MxWCI2GxikhJpnRD7o9jJfLnN4PlrpUspvmOeHmbxQUIbxikHQh.fhqxsJVwJBCCCjVZogu3K9BbG2wcfJW4Jit0stge3G9A7DOwSfYNyYhgMrgEy0u+SWa+kfUHFQDQDQDQDQmGPuGfAT7+H+1291igLjg.KKK7Vu0aAaa6XBUROHKe97gku7kia7FuQUHJ..gCGFtttH4jSNl.h777T8pqHQhDSEc8RuzKgYMqYgYMqYg9129pB8vwwQMz6he1UTttku7kGUpRURUsOxvyhOLL4+ImAMadyaNLLLPt4lqpplj8BLKKKDNbX0eOiYLCrpUsJTkpTEby27MiIMoIgu669N7oe5mh5Tm5fssssoV+et8PLYEXIOlhDIBtwa7FwnG8ngmmGF4HGIxHiLhILLfXG1pxp9pz115gRlat4Be97oFtg5qiOe9vfFzfvHG4HA.vq+5udL2u..ZUqZE..N3AOn50NxqYx6Cx9slrZ7zCvD.pWKYYYgPgBo52b974SMjFqPEpPLAjoecz00UMClJ0nF0H333fCdvChvgCqNl..NyYNC..t5q9pQ4JW4fggA18t2Ml6bmKti63NPCaXCwd26dQxImL9hu3KvvG9viYlqTt++sZl1lAhQDQDQDQDQz4IjAGIGBZtttX7ie7nJUoJ3UdkWAG5PGBdddpmWFxgPHvm9oeJNzgND5Tm5D.fpJwjUqUnPgfoooJfH8g7mLvCYEO8lu4aBSSSTu5UO344oB9XaaaaH6ryF.msuOIOVjAD0fFzfXpNo3oOjI06cTCdvCFBgHllAu7bywwACe3CG6ae6ScbNyYNSXaaGSHTctycFCZPCJl8WoMAETZ74ymZXDJDBDLXPDLXP7XO1iglzjlfHQhfQLhQfvgCqpHJ.ficrig0t109SFFlj77MqrxBNNN3vG9vpptRt8jK2Dm3DQG5PGTGax8ommGZSaZCpXEqHV7hWrph.kg4YZZhMu4MqpbL8.KkAfJC+Cn3fwjATEMZTUEj42ueTTQEotOEeUXEHP.DNbXUklYZZhZVyZhG3Ad.b5SeZ7C+vO.+98qBObyadyvzzDCe3CWsMdwW7EQiabiwjlzjvccW2Etq65tPu6cuQiabii40lw2P8kWW+4be8bgAhQDQDQDQDQz4IzmwHkAST0pVUrxUtRzvF1Pz0t1Ur4Mu4X5YUddd3se62FO1i8X3EewWDACFLlJXJmbxAFFF3vG9vwrO777PgEVH.fpmTEMZT0L6nqqK9vO7CQQEUDhFMJ1zl1DF23FG5ZW6J..97O+yw28cemJDI4v36HG4HnfBJPsuhOTl7yOeHDBjUVYESHGibjiDokVZ3oe5mFKbgKTcrUPAEfm5odJzt10Nzrl0LU.bqd0qFSXBSPMaWB.jWd4gO4S9Dbi23Mp1ekVkrEJTHTTQEge3G9AUXNe228cnvBKLlgyXnPg.vYmwDW8pWMl3DmHNwINArssQjHQvfG7fw0e8WOdy27MUgCoGVi7XPHDnfBJ.ey27MXpScp..XhSbhHiLx.ETPAwbewxxBUtxUFu5q9pwb7KqRu5Uu5goLko.+98igMrgg8u+8qVlzSOcLpQMJLxQNxXFdm4jSNv00Em5TmJlq6VVV3Dm3Dvvv.m7jmLli6hJpH..r+8u+RzmzjAidricrXFpuiYLiAMpQMByd1yVsrm3Dm.ye9yG8rm8D8qe8SssrsswLlwLvS+zOMl5TmJ9q+0+Jl5TmJdwW7Ew7l27vt10tNmgd8qtOhIHhHhHhHhHh9eNOOOgqqqv00UHDBgsss52CGNrHqrxRLoIMIQsqcsEolZph68duWQ+5W+DMsoMU7rO6yJN0oNkZa433HxO+7E8t28VT8pWcA.D.Pb8W+0KV0pVk3Tm5Tha8VuUQsqcsE.PXYYIt5q9pEqbkqTXaaKxHiLDW0UcUBCCCQEpPEDWxkbIhG7AePQ1YmsX6ae6hzRKMQkpTkDie7iW344Itm64dDMtwMVss5XG6nn+8u+hPgBoNl91u8aEcsqcU32ue0xcK2xsHV+5WuvwwQ333HxM2bEibjiTXYYIBFLn3JuxqTT+5Wewq8ZulvwwQc83kdoWR7TO0SIVvBVfnicrihG6wdLwy7LOinicrihktzkJDBgHRjHBWWWgmmmPHDhnQipNVl6bmq5ZhkkkvvvP..gOe9DCcnCUs78rm8TDHP.0xYYYoVu64dtGgPHD29se6BKKKwG9genv11tD2S0+ugNzgFy9KXvfB.Hl4LmoZcheaL6YOaQUqZUUaC8ykO9i+XQ8qe8ElllhzRKMwkdoWpnO8oOh8su8oVtku7kK5Tm5j57npUsph9zm9HxHiLDuy67Nh1291qNmZbiar3Nuy6T7i+3OJt0a8VE0st0UccoScpSh4O+4K13F2nnyctyp0o7ku7ha4VtEwAO3AUu1KqrxRz291WwMbC2f3oe5mVzwN1Qwe4u7WDgBERcNXaaK9pu5qDsrksTc8H9+ymOeh0rl0ntu+aICg3WQ8kQDQDQDQDQD8aJgVuAy11FABDP8bxgL2t28tQnPgPcqacQ0qd0UyzexgUmdSpW1SnjCewnQiFyxK2mxdEkiiCrrrPjHQvgO7ggkkEpYMqoZFHTVUNEVXgpYiP41W7uGdixgimdSwW1X1kOF.h43RuI0mWd4gLyLS32uezvF1PUSkW974latH4jSFABD.EUTQ3HG4HHb3vHkTRAIjPBk35f72k8nM41QdtHuln2yy.Jt+qIGtnk1rBo7Z6AO3AQCZPCh4Zg94l9uKOdj6a4833WF49zzzD6d26FokVZHZznvzzTM7NkWWN7gOLN4IOIpQMpApe8qu5XS9S4wf7Zn9qUzqzJYkaoOqbFJTHjXhIFy4b7uVU+4zerie7iihJpHT9xWdT0pV0X1WtttXG6XG3Ye1mE8rm8DW8UeU2LOF5...H.jDQAQ0vwwQMA.rksrEL7gOb7vO7CiW7EeweF+ub9kgAhQDQDQDQDQz4Aj+yyE+6gFm9iKGFcx.Zj9oBnPOfI8PHzCcQNqJJCDR975y1hQhDQEJj9PTSt8ieapSO7q3We8siiiCBDHPIVFYuqRODJYfV.kLbO4wh7wi+ZjdvXx91U7yNi5qS7GOxqK5KW7y3l5G65yrlwecK9is3WW41TN6PJOez2dh+8PyL9qG5Gmx6AwGXp9wg99Sttwe7KCOS1n7kGO5qa7u9L9mW9ZbCCCr28tWz7l2bL6YOabe228UhkE.XnCcnnhUrh3u9W+qwDL7uEXODiHhHhHhHhnyCn2b50av5x.WjASAb1fPjUADvYq.H.npjHGGmXBZyxxJlssL3CYnRx0W93FFFHgDRP8b5qqblMT+wkU0liiCrssQvfAQjHQTmGx9Mld+nxzzLlpyRt7xvQjmixJXSFxmPHTgjHqbNYXJxv8jGS5WikGKRx.dz6cXxPpjgQJWW40EfhCZR1v2kUFm77CH1IHA8iC40A44i9LwnbhPPup474ymJTQ45JWG4qWjWOkgJJCyT7umb.jWm0e8f79T7WijMBe8yKaaaU3X974Sc8011FVVVpJ6RuZ5jWmzm8K0eczQO5QgOe9P4JW4TyjjxskiiCxN6rwl1zlPW5RW9MOLL.VgXDQDQDQDQDcdK8gxHPIqJq3q1mRaHrAT7varbkqbwrciFMJBFLXLA0HW+3qjJ8sqblETe3ApSuxiJsgmo9xnG3z+opoRmdEsA.U3MxYlSY3U5WCzOGhuxphe3UJCOTJ9pkRux2JM562RqB0z2tx6ixvujaeY3Vk14uj94gL7JSSyXFxlxqu5Ug24Z6p+5M8prS+7Ruh3h+Zk9vxT+7UFbl9qohFMJFyXFCdsW60v3G+3QqacqQvfAQQEUD10t1EV7hWLdfG3AvfFzfXfXDQDQDQDQD86Y111pPQhenOBT5CIRfyFrx4ZHQdt115KS7aaYPGkVeHC3rgRIC0R1enjAqEeEsE+xDevYxp+RNLIkUX0O0v+Kdmq9ul99N9giH.T+c7CWS8mK9vBie6Fe.RxkO9gYY7z6QW5C4vRKDQ.Dy9V92mqgNYo0WxzC9TNqhJqJr3Gpt5gkpud5Odoc8W+XH99Xl9PxE.HyLyD6bm6D4me9HTnPH4jSF0st0Est0sFVVVk5qq+s.CDiHhHhHhHhnyCnGxP7ULT7AYnWgUmqp84mpGeE+16b0uqJs9FF.J0pbJ9eW5bUsTx+tz5IXmq8eo0Ozhu53hecjgnE+0pRa4+otVE+yet5IZ+mBgJ9pnRtLxgond3fxq4we+8+TuO6b8XkV+CS+wOWud5+z0X8WW7e5Zz4ZhTH9sao0Ox9sDCDiHhHhHhHhHhnxTXS0mHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfX+NmqqK..777..fiiC..DBw+yNlHhHhHhHhHhn+WhAh86bVVVHZznvzr3a0BgP8eDQDQDQDQDQTYQLPremKZznHPf.HZznvwwA986GNNNp.xHhHhHhHhHhnxZLDrTgJyxwwA9746+0GFDQDQDQDQDQz+UwxD524DBgpugI+oruhwvvHhHhHhHhHhJKhAh86bqYMqASaZSC.EG.lssMrrr.vYaz9DQDQDQDQDQTYILPremyvvPE7kPHT8PLgPv9HFQDQDQDQDQTYRLQjemSFBlqqqJXLe97wYYRhHhHhHhHhnxrXSj524jUHlbXRB.VcXDQDQDQDQDQkowTQHhHhHhHhHhHpLEFHFQDQDQDQDQDQkov.wHhHhHhHhHhHpLEFHFQDQDQDQDQDQkov.wHhHhHhHhHhHpLEFHFQDQDQDQDQDQkov.wHhHhHhHhHhHpLEFHF8aFWWWHDB333ThmywwABg..P8yRimmWI9cOOO355Vhsmj91y00Ml0iHhHhHhHhHhhGCDi9USFBkkkELLLfOe9.vYCASHDvmOevvvP825jKmiiCLMMUgeYXX.WWWXZZBKKK0xGNbX3ymO344AOOO01E.piAgP.SS9xahHhHhHhHhnRhIFP+pYZZpBeJZznpvsjgfEIRD0xZXXn9cYHVxkSFjlL7K8fyhFMpphuRHgDP3vgiY+pGzV7gjQDQDQDQDQDQ578+5C.5BexpCCn3.ujga444ASSSjPBIfUspUgu8a+V7XO1ioBwRFvkdkbkUVYgkrjkf8u+8iF1vFhd0qdgK5htHDHP.344AgP.CCCjPBI.WWWXaai0rl0fMrgMfPgBg1zl1fANvApBjiHhHhHhHhHhhGqPL5WM8Jwxue+Hb3v3HG4HXMqYMXlybln28t2n6cu6Hb3vvvvP0+uzCQy00EqYMqAst0sFUspUEibjiDUnBU.srksDqd0qVMDH0q7qhJpHLxQNRrnEsHz291WLnAMHr3EuXz+92+X5YYDQDQDQDQDQjNVBMzuJxJwRVMX.EOjFm1zlFJpnhPgEVHV8pWMrrrfe+9gqqqpxsDBArrrfqqKN3AOH5cu6Ml3DmHtq65t..PKZQKfqqK5ae6K9hu3KPZokFBDHfZa7zO8Siu3K9BrwMtQT0pVUXXXf23MdCz912d7vO7CiYNyYFSuGiHhHhHhHhHh.XEhQ+JIC2xzzTMaPJDB7RuzKgYMqYg+w+3efm64dNHDB32ueXYYgvgCqVeYnXu5q9pHb3v3du26EtttpFqee6aeQ94mOl+7mOBFLnZHYt28tW7xu7KiQNxQhpW8pqp7rJW4Ji69tuaL24NWrm8rm+mbMgHhHhHhHhH57aLPL5WMgPnpPL4r7nggg52qXEqXIpfLYu.yvv.YlYl3UdkWAW20ccHPf.vxxRMbGqUspEtlq4Zvq9puJNxQNBbbbfkkEd+2+8giiCtrK6xTGG986G..ssssE..Ke4K++xWIHhHhHhHhHhtP.CD6B.5yfhxJvRmrZpNyYNC13F2H9rO6yvV25VwO9i+nZFeTe4heFXT975Ki7mw+3.P0b6kLLLhow3GO4LKoooYIN1EBA15V2JbbbTMOe.nFFl..MsoME974CaXCaPUQZyZVyBlllnV0pVpv0ja6ZVyZB.f4Lm4Dywq99l8WLhHhHhHhHhJ6h8PryyIC6wyySU0U.HlJt5HG4HXRSZRXtyctvue+v11F.EWcU2zMcSnoMso..wzOsjqubFgz11VMjFiDIBBFLXLKmrpsj8uKWWW05+qgggANwINALMMQkpTk.Pwgu433nNFJe4KODBAxImbfmmGJnfBvgNzgfggARLwDUyxkxiE4Pq7HG4HHu7xCku7kuDGmx.z90d7SDQDQDQDQDcgGFH144zqtJGGmXBEKZzn3Ue0WEiabiCgBEB..111vmOevwwA4jSN30dsWCO2y8bpJKSug1aXXfbyMWbe228opjLOOOTtxUNTXgEFSHSxe9hu3KhFzfFDyi8qgqqK1291GLMMQvfAUMoeKKKUHcxf+N0oNELMMwoO8oUUrlOe9Jwwg73100Em5TmBImbxpkQVYXxgqIQDQDQDQDQTYOLPryyIGteVVVpvrjl27lGF6XGqZHOJqNLGGGU.ZxGWuubAb1lgukkENvAN.rssQvfAiY3PJmAHANaeBCnjCmxeMrrrPgEVnJLNe97oBqSNDKSHgDfqqKN5QOJ..xKu7TgjIOWkGi5CcR8pkS+4kUaGQDQDQDQDQTYSLPryyoG7ibnMB.r10tV7POzCUpg+32ueU3UVVVnhUrhp0WFfjooIbccQEpPEvV1xV..TUjkssMrrrhoufou+0G5k+VHkTRAFFFHb3vpyGYHfNNNpgmYspUs..P4JW4hYHiJo2n9877fmmGRN4jUAqwgHIQDQDQDQDQD.CD6BBxg5nmmmZlTb1yd1pPuza.8xgZnjqqKZUqZkJ7H8FfukkkZ8LLLTggI2Gtttpvuj8ML+98inQip1FwW0Z+R455hZTiZ.gP.aaaHDhXB5xmOeH+7yGdddH0TSE..0st0Ucr555pNNkgcIO+MLLPkqbkUaq3G9mLfLhHhHhHhHhJahAhcdN8PojxN6rwRW5RUU7U7CYRfy1z3+i+w+H5PG5fJ3G41RF7kqqK14N2o5wjARo2GtjaKWWWbwW7EijRJoeyBRxxxB0t10F986GG8nGElllpP5jU10QNxQfooIpd0qN..BDH.tpq5pv5V25vIO4IQiabiiYadpScJ..z8t2cjXhIptdHE+rrIQDQDQDQDQTYKLPryyICvReleLTnPpgSHPwADEMZTUXXxd+0e3O7Gvzm9zggggZ8k8BLYUfkWd4gNzgNnVW41RRO7HKKKr6cuajZpopFBhwOrE+kxyyCW+0e8n4Mu4XiabiwLTHkGqaYKaAolZpnScpSp06gdnGBqcsqEG6XGSs7xYaR4rV4.G3.UaG8yi3qTLhHhHhHhHhnxVXfXmmSVgXxFduooIpe8qOpe8qOxJqrfiiSLAXIWm9zm9fYO6YixW9xC.nBSSeVpD.HwDSDKcoKEgBERUwXxgEorA2Cb1gsYMpQMTGG5AJUZUxlLHJ4yEMZTUvTxeJ2NiYLiACcnCE6XG6.spUsRUAaaYKaAG6XGCOyy7LHPf..n3fu5V25Fpe8qO93O9iQu6cuU8bL.f+4+7ehZVyZhd26dqNm0+Y7+NQDQDQDQDQTYKLPryyI6IX986WUITlll3se62F8nG8.4kWdvwwQUYWokVZXhSbhnu8supgd3OkjRJIzktzEjPBI..DSEZEeiqOb3vpkSVgXxfnzCCSHDH+7yGIjPB3zm9zX+6e+v00EYkUVH6ryFIkTRvue+HXvfpvx5Se5CV0pVEd1m8Yw7m+7QEqXEQQEUDdpm5oPm6bmwfFzf.vYat+IlXh3e7O9GnG8nGnW8pWnm8rm..XwKdwX4Ke430e8WGIkTR+FbGfHhHhHhHhH52aLDrYJcdMYnTxppRu+Zkc1YiO8S+TjUVYgxW9xiV1xVhV1xVhJUoJAOOOU0d8SQFvD.T8KLSSyXlUJ0aX8QiFUUoVREVXgkXleT1v9k80LYEgUtxUNDNbXbkW4Uhu7K+R01PN7He3G9gwV1xVv0dsWK9m+y+I5RW5BlzjlDpXEqnph0LLLTGee9m+4XDiXD3Zu1qE4kWd3a+1uESaZSC8pW852r6ADQDQDQDQDQ+9BCD6B.5CwPYvTtttwL6SF+y+Ks+dUZ8tq3q5KCCiXlEJ0CSKRjHwLjFsrrJw1wyyKlYFSgP.GGG0P0TF11wO9wwIO4IQ0qd0Qcqacio+oAb1P4bbbfggABEJDN7gOLLMMQyZVyTULGmEIIhHhHhHhHhJMLPrKPn2D6cbbfooILMMgqqaL8yK8d4kdfUmKx.tJsFOeo037iuh0j6iRK7oRKHOc5A4ouOjay3CoSNSWFeUuIeNYEzIqNsesM7ehHhHhHhHhneehIFbdNY.Ux.dBGNbLgbIGNi111wDFlqq6+wvv.fJbIYvT111wzH80CURF9lrxu.NakdoGXkj9P7T9S4rYodC0WNaYJ+o73Vtej+zzzTc7FJTH09IZznpsk9PJU+XgHhHhHhHhHhjXEhcAB8gMndUZAD6Lln9PZrzl4GKMxgdX7CQR8.whuRzjhuRrzCj5+zPVz11F974KlpHS1yxNWA5oW0ak19V925S..DQDQDQDQDQjNFH144hOzK8.uzE+vi7m6PFL99RldElI+c8iA88SoEHmb6ou9xp7R+u0GJjk1xaZZVhepedFe+RS+2kGSrGhQDQDQDQDQDUZXfXDQDQDQDQDQDUlB6gXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCDiHhHhHhHhHhnxTXfXDQDQDQDQDQDUlBCD6BbBgHle2yyS825+NQTYWQiFsDOl9+eGDQDQDQTYABgP84fcccU+alkOl9Oi+24me92eXfXWfyvv.dddp+AullEeK0yyS86DQkMIei6.AB...GGm++r2Ud7QQUx+ucO2gDRBWhJfhfx4OMBxBBhhfbHhb3sKhJhvhGqJqWHhm38w5Ef.JGBnHhJfxBBKHhJH3F4PEAAjaDIHIPtlyte+9iX8R0cFHCjIISBuue9LeRlY5o62QU0qp5UU8jetllVkYSSAETPAETPAETPgJEP5AysYNPf.xuyvv.ZZZx+2zzT9dEpdAkGSphiPgBAccc31sanooIM3UwrpfBJPKlGNbX..3vgCDIRDk7AETPAETPAETPgSJAoGb3vggKWtjYYkOe9jAZBcMBg.Nb3.Nb3..pLvp5HzDp39qJOLLLjLobnhRLETPABzh9DNZxMTPAETPAETPAETn5JHGfw0ClGTIzmS5JyiNLUVVT8CJukTEGBg.SaZSCspUsBabiaDQhDQlayJlUETPA.fUspUg1zl1fYO6YCfhBIbkyvTPAETPAETPAENYCjSuBEJDLMMgggAb5zIV25VGN+y+7wzm9zA.j5JyinLk80U+fxgXUwgllFV4JWI9ke4Wvu7K+hkv4TwvpfBJXXXfMu4Mie9m+Y7ke4WB..ud8pJJnJnfBJnfBJnfBmzA60XWGNb.gPf0t10he4W9E7i+3OJK2H7ZIlR24pmPkxjUCv1291w1291wEdgWH74ymJMIUPAE.PQg+sSmNQgEVHV9xWNZe6aORO8zgCGNPvfAgGOdprahJnfBJnfBJnfBJTgA6o8Houbf.AvW7EeA5XG6HpacqqLxw3kcDUJSV8CJGhUEGTcBixwYhI0zzTVD.UPAEN4EztfoqqKWvWU+vTPAETPAETPAENYE1O0HiDIBb3vgzNZJ.SHaqEBgzAYJT8BpPIJAGjwrD3mrEFFFRibilwsTQ.j+6LLL.PQmNk7OmdFpSNCETn5EnBGJ.fSmNkxKHmlCTDeucYEGs6APwEdTETPAETPAETPAEppAcccKoAI2tYxoXz+S+kB9D52X2NcNNZ1uqPhGTt3LAGTjcXZZhvgCCOd7HixisrksfBJn.XZZh7yOezvF1Pb1m8YiCe3CiksrkAud8JMzsV0pVn8su8XQKZQ.nnZHTAET.b61Mt7K+xwF23FwN1wNjdJ+TO0SEssssUERnJnPUXrqcsKTPAEfBKrPjc1Yi5V25hLxHCje94iksrkIitTMMMT+5WezxV1RrhUrBKNGysa2n28t2XSaZSXyadyxid55V25h1291WI2CUPAETPAETPAET33Cj9uBg.G7fGDYmc1H2byE986G0st0Est0sFETPAXQKZQvkKWxHDq90u9nScpSX9ye9voSmHTnPvsa2vvv.W9ke4XaaaaXqacqvzzDtb4BMnAM.m24cdVBBEUVZjXAUDhkfCJbM000gGOdjNHC.XJSYJXRSZRnqcsqnacqaXW6ZW..H4jSFNc5D24cdmnu8su3q+5uF0oN0Atc6FojRJXJSYJnW8pWXEqXEn10t1vgCGnl0rlXO6YOne8qeXtyctH4jSV4LLETnJNVwJVAF23FGF3.GH5YO6IV9xWtbw+jRJIbW20cgALfAfu9q+ZjRJofjSNYTyZVSLgILAz+92erhUrBjVZog.AB.ud8hCbfCf92+9iEsnEgzRKsJ6tmBJnfBJnfBJnfBGWfbJkttNb3vA9hu3Kvq9puJFzfFD5V25Fl+7mO.Jx9aud8hG7AePLfAL.7+9e+OjTRIA+98iTSMU71u8aiq9puZrvEtPTu5UO3xkK30qWrqcsKLfAL.Lu4MO3zoSo86zySgDKnbHVBNDBgL8jHOYqqqCgPfW7EeQL9wOdbG2wc..HY3..5W+5GF0nFE..ZXCaHZQKZADBA5bm6Ld3G9gA.vodpmJZe6aODBAZTiZD5e+6ONyy7Lwq+5uNZdyatL8JUPAEpZha7FuQLtwMN7BuvK..ffACBgPfTRIEz0t1U7XO1iAccczjlzDzxV1RXZZhK9huX7HOxi.gPf5Uu5gN1wNBud8hy5rNKzyd1Sb5m9oiW3EdAbNmy4TI26TPAETPAETPAET33C7Tczvv.25sdqXxSdx3EewWzRpN5wiGzm9zGbe228A.f5W+5iy+7Oe3ymOzktzELxQNRoczssssE..MtwMFWwUbEnwMtw30e8WGsrksTFjIzIWoBIVP4PrDbnooAmNcBSSS3vgCDNbXYQ.TWWGNc5D0nF0.ZZZkvAV8rm8D..u268dH+7yWl6y+3O9ivkKW30dsWC986WFwHqZUqB24cdmHkTRA.PcZUpfBUwAch3Pmljtb4R5fcGNbfd0qdASSSLgILATXgEJqMB+zO8S..XBSXBH2byUd+xLyLw8bO2CRIkTTQPpBJnfBJnfBJnPUNPkLDNnRRDPQ5KGNbXosv8oO8AZZZXFyXFHPf.RcoW6ZWKb4xEl3DmHBFLH.Jpbir90udL7gOb3xkKoc1BgPpWtBIVP4wiDbXXXffACJYHc4xkkS8Bx4X.EwHSNPKb3v3LOyyDCYHCAqe8qG+7O+yx64TlxTvPFxPPN4jC91u8akEY62+8eebEWwU.MMMoixTPAEp5BpFBRxI3Gaz..m1ocZXXCaX3m9oeBqcsqElllvvv.yblyD20ccWXe6aeXMqYM..vue+X1yd1nu8supnGUAETPAETPAETnJK3m55ABD.Nc5DNc5TVqb4GnTMrgMDCZPCBqcsqEYlYlxZG1G8QeDFxPFB1wN1AV4JWI.JJXVl7jmLtxq7Jk0VL5jpD.JcnS.gxgXI3vgCGviGOVNQ33oPIcJSJDB31saKdf1zzDCX.C.Nb3.yXFy...e+2+8n8su8X3Ce3Hb3vXNyYNvoSmXqacqnfBJ.MqYMCFFFvmOeUZ8YETPg3CnHIkjePNCyvvPdjRe0W8UC.fO7C+Pnqqi0rl0fLxHCbG2wc.CCCL6YOa..r6cua32uebVm0YAGNbnVPWAETPAETPAETnJGHGcQ0yKud8JKMQFFFvue+vsa2PWWGgBEBZZZ35ttqCNc5De3G9g..3G9ge.ssssECcnCE974Ce5m9oHRjHXqacqHb3vnEsnExrufBlkfACppgXIfP4PrDbPFcxOVXIOZSg0IcxPZOhwz00wkbIWBpe8qOF+3GO9i+3Ov7l27v0dsWKNuy67vUbEWAl7jmL1291GV3BWHt8a+1kdEWc7vpfBUOAEQnzhyW3EdgnAMnAXbiab3.G3.XQKZQ3pu5qFspUsB8su8ESdxSF6YO6A+2+6+ECbfCDtc6F.PsftBJnfBJnfBJnPUNP1JGJTHKeF4TLOd7Hq2WT.mz4N2YoM0G7fGDyadyC8su8EsoMsActycFSZRSB6cu6EKcoKECdvCVd.VA.YvqPkvDERrPBuCw3E8NxIMaYKaA29se6XAKXAVtVx4Q7hgWUcPNmhXbMMMkFzxAU78s+4zIIG.vrm8rw7m+7wEcQWDBEJDtsa61..v7l27vzl1zvkdoWpkiDVBTDowgxgYJTQ.djQZ+yiFsZzjWTcF79azFmHGkSNMmpEg.EeHcjTRIggO7gKSa5O4S9DzoN0I..oLh4O+4i28ceWz8t2c48lFeI4C7w6nIyPAEpL.Oxp4Q0X0I8DTPgDQP7a1WOfy6Q0bG63jg0uUPAEp7.oeLsIu.PZuMOaq39Vnl0rl3Nuy6DBg.u268dXAKXA3htnKB..CcnCE..e1m8YXpScpR8kofUgaetR+iDOjv6PLxPNBqd0qFctycVlatKYIKA.Pl9OzuAn5kQY7iqUJxvHFJp+RLaDyG8a5cu6M..tu669vse62tT.vkdoWJpcsqMt268dQG6XGwobJmhbLTWWG4kWdX4Ke4H+7y2hGz4EYPETn7DzIpJwKGJTHIutc9cZALZGdNYfFkjORo.IAgPH+LJsIoZHlc4ENb3.8nG8.FFF39u+6G2wcbGxPGuicriH4jSFOzC8Pnyctynd0qd.nX4sG4HGAqbkqDETPAxz2lhdUkAMJTYCRd.UnbUJjpfBUbfVC9noyJPwG3KjNqz5FmLr9sBJnPkKHYQABD..PVCw.fkCfJt+D5QO5Az00wC8PODF9vGNb61MhDIB5RW5BpYMqIdfG3APFYjARKszj5nqqqi.ABfku7kirxJKUM5NADI7q3XZZJcFSlYlI5ZW6JxJqrj0BmAO3AictycZYwSSSSDNbXIQcUYPN0xty8BFLnEG+40qWoxDQhDQ5HA.fV1xVht10tBWtbgt28tKMJHkTRACYHCAFFFx5HT3vgkQNxa9luItrK6xvHG4HkdPORjHVpCQJnPEAnzDlbDlggg7U3vgkF9RNE6jEko4oTMwWSGtFNb3vhy.LMMk0FP+98KuGQhDAWvEbA3RuzKEZZZnacqaxcypt0st3dtm6A986GW0UcUxwbZNXbiabnKcoKXDiXDRGgQEmzSVlCTHwEzZjbcAnhaqh9TAEJegllFBEJjLsi3ETZ9FcQ5bx2.bkCqUPAEJuAImwqWuRmxS9bfe.Tw2D91zl1fK5htHHDBz0t1Uo910oN0QVetusa61j1sPx4F6XGK5V25Fd7G+wqf6kJDKHgWiPccc3vgCjat4hG7AePIgU3vggGOdv92+9wDlvDjFARJ55xkqpEQHFuXXSLolllVxAYxyyzIPIUP846N2.G3.Q26d2w4dtmqT4DGNbfq4ZtFzrl0Lz4N2Y.Tj..RQk5W+5CSSSjZpoBfhDbPFVaO7OUPgxKv2UYmNcJc9EUu6niFY.HcBD4Hnp6vgCGkHhMoS9FPqoPFK...B.IQTPTghGyn2SNPymOex+mNUZG5PGJ5V25FZdyatkHy6ptpqBsrksDW3EdgxwW5YmVZoAGNbfzSOco7fvgCCud8VgMFnfBGMv2PMtLB5zTUAETn7Etc6VZXIoOKsFMuN3xWmmeZrofBJnP4AnLnfeHSwKsHd73QV2uHcrIeMLrgMLb4W9kiV0pVY4Tb+ZtlqAmy4bNnMsoMR4Xjiw750Kb61MRKszT5ej.BMQB91vDJTH31sa7xu7KiG5gdni50sgMrAbtm64JcVS0gnCi.24VDn5ElPHvnF0nvK9huHxLyLQaaaas7anHkK2byEG9vGFmwYbF.n3TIwue+XW6ZWn4Mu4VhrDhrXe6ae3zO8S2hvBBJkVTnhBbd.RlPnPgrbhnFNbXoyfsmp0U2AIOHb3vxn4jFmLLLvm+4eNtpq5pvjm7jwfG7fkQMF.jG4zG9vGF4jSNnwMtwxeKPQQS1d26dwYbFmg7ynwZ.fsu8siy7LOyRDkt.pzdQgDGvqsf1WKSAETH9iPgBAWtbIM5jWVOHPFeRo4e0Ic2UPAERbAWND4HLccc7we7Giq8ZuVLsoMMLnAMHodrbYX4jSNHu7xC0u902RMHKPf.Xu6cunoMsoxmAWt11291QiabiU5ej.hDdqUb3vA9se62jNCyiGOkvXWud8hkrjkH8hqSmNOpEpyphfNxWIOUmYlYJOZW000we9m+Ib61MZRSZhLpNHFXRYjjSNYoyv.JNT184yGZRSZB.fkSXRZLtAMnAHb3vVR+J.kyvTnhA7cQgRyBZwG2tcCe97Ic9B4fFdcHo5Nre53PiM6ZW6BezG8QHmbxA555H6ryFBg.MqYMC.vRJUSiwokVZngMrg..RGNB.oLB5daXXHc7nPHPiabigttNLMMk0gA98WAEpr.QCxKw..ViHEETPgxG31s6RTyJCEJTIhBLRmS60tGETPAEJu.4Dre+2+cL8oOcjc1YC.f7yOe..Y1R..KQIF.P5omNZTiZDb61MBGNrzlCOd7fl1zlZoTtvKmImwYbFpHTOAEUIbH1zl1zfa2tgttNBFLXITjMPf.XJSYJnvBKTlVDd73oZiBuABD.tc6VF1k2vMbC35u9qG6d26F4lat3S9jOAO7C+vHkTRwhR+TDz..Ko2Hw3RGirzug6gbxq1jCHhDIB762ukPbWAEJuAULKo72mnQEBAV25VGNvANfzgw.EuCNjicqtCxwUjrApdsLwINQLnAMHrrksLnoogu5q9JjQFYf1111ZIBYH9ad5SZW9.PwQ5E+Dukm5KD750qTFrJ5vTnxFadyaF6ae6SVhADBgkSYUETPgxOvqSX.Es9AswUqcsqE4kWdknzaP0SLETPAEJOAYuvK7Bu.tka4VvxV1x..vRW5RQaZSaPKaYKAfUccI8FH6PnRZjc6M3QFKE.JTM1kdoPhER3sXIu7xCSe5SWd5xooogQMpQgcricfBKrPrt0sNz111VroMsIr8sucYMvha3VUYPELedpf1oN0InqqiUtxUhQNxQhK6xtLb228cKcjEYTqa2tipxFTzb.TbMFxNyL2q1.EYbKEMNUWbznBUM.QGZme9hu3KFOvC7...xnW5jMmvPoHIEEmT8Zgh5ycricfm+4ed7y+7OiYLiYHivVfhWj2dJYycvHwqS61E4XR54Qf6H8S1lCTHwEm24cd3QdjGwRJ7RFnqRMKETn7E7n9hz8TSSCKZQKBsqcsCye9y2xgzBPw0M3SF1PKETPgJOP5BmQFY..f8t28hm4YdFrwMtQLsoMMjRJoHuVtNDjdtTcCCn3MS1tcJjdxb8MneiBIVHtVCw3K3ADaoUG2XLxaqjieLMMwhW7hQe5SejDZu8a+131u8aW9bLLLvJW4JwkdoWJ9fO3Cv0dsWqTo2xKEdiV+JZGmzwqzJzdeI2byE+7O+yX26d2nEsnEnIMoIH4jS1xXY4Y+WAEpHAwGwyeefh30pScpC98e+2kNBmexqVci9m3usONPeN4rZGNb.+98iMu4Mie8W+Ublm4Yhl0rlgzSOc40yiTFETn5JzzzPsqcswANvArDYXUmn84JpS8s3UeLZ00oilbnJKDM49pR5PhC30QLRO9YO6Yia7FuQLrgMLLwINQ.TRaA.N4KSD308yRi2phhFmZK10oJd+7sKW4jId3n0mqL5+GK6Xo4mpS5VS8ECCC7i+3Ohsrksf5W+5iy67NOjVZoUY27T3uvwZ8.trxxZcjNtPUSLJ7SuImNcFSMFxPNZgPdp7nqqiu669NoBe29se63VtkawR935zoSjQFY.SSS7q+5uZ4jUr7BQqeEqe1IBn9BYHaMqYMwe6u82P6ae6APwdfl6ayDAEUUPgxJ3mBL7ZPBEQXG4HGAaaaaCMu4M2R3LSNTu5RzJQodBk531+NfhkkBTTc+5bO2yEm64dtk3Xru5xXhBJTZPSSCG5PGBaYKaAsnEsvhifqNnXOOxL4HdY3h855DEgnqZUqBEVXgVpggUVfjukZpohN1wNZ4vARgJWPoSDoSuciqWxRVhbsc95WmrrFE2fbJBuCFLH73wCV7hWbIpYu1WKuhvgIQhDA974C986Gst0sFm0YcVwMYmT+HZ52bxhyvn46EsnEIKKP7ST8xSP7j.nDNPPHD3LOyyDspUsRd3rAfpUa3.0Wb3vAxHiLPKaYKkkWgpS1OTUGGKZMd8njmMQmHzmkYIZ1O4lrqbVowPSEpc2tcKOsYnvPzkKWX0qd0vkKWHb3v3QezGUV+A3nl0rl3BtfK.qe8q2xBpwizlzdw2kufj85uk8OOdT3d4FwxGasuSMzBm1UBTkdiJTUFbGkSfehtDNbX7a+1ugl0rlAMMMKQnZ0gZQBW9B0uoOieROxMLlTd4XYPnc4WJnP0YrsssMzhVzBDJTH3wiGYTnWUe8Q9oekSmNsDgEwi9GoeFMdQ5UMlwLF7EewWToK+fbtB.vke4WNl27lmkC+CkAMUt3nk9PTI.YW6ZW329seCmy4bNVVOi1zqp57mkF30jWZ8XOd7f.ABf9129JKUADpLiXHccc7TO0SgQNxQBmNcJOYqKKfqeCM2y0uop9FVTZfluCEJDFxPFB1+92uzVOxw9k2vtCwHZQSSS7XO1igG+web4bCISknUqpyeRqeP9cfzMvtimUnxAjtLb9.xttnUhVrGfUGunLKsgy.QKnwCY+XMJw.r5LGWtbg7yOe7e+u+W..7POzCgF1vFZo3NSBkc5zIxJqrPZok1Q8Dgnr1Go6o89cztV9mGObHWzdOulhwWLg9tS11kEEp9CdThQJKHDB7QezGg9zm9..qFvwcbVUUX2A+7EuIddJT6ouytQGQSN7I5NnnfBUk.sl8G+weLtxq7Jk69KudEUUFj9ObCqi1QDeYAbiRo6IUCAqrOXBHCy30OUxnMUDhU4C6aNt80bnZgaSZRSrnCq8C7opqfG0k.Ec.Z4wiGYcC1tA4UzFnSAiPzbJW7vYU1qKSj7L61zTcEjie4mp27.bnhb9lGsXDxKu7jyy75lk8.AopJ3QLt8ZVXUccCpN.JSeH4j1OLBhVjxx8ATkRJSRMHt.xiGOrZ26eTGdW6ZWxc4r28t2VtWzNg5zoSbvCdPr6cuaKBPocBNdEgVGq759n8LhGdP2djmQeFeQTdahuqnJuaqPUcvSWRdpASxZz00w6+9uOd1m8YwobJmh7zXUHDxZWRUYXmGmqbB8hb3EoDEOp5rWeFAr5LLkLBENY.yblyDOyy7LnAMnAxz3JZF.TUCb4A10UfbTdYA7nEfGg5G3.GHgHBboME0vv.6bm6rDmFtU0meqpCdjOwoU3EZ+IO4IiAO3AG05UW084OZMZh+h6HrRyViJBdO6E+6csqcEU4Amnvd.LXWekJa4Kk2f1z1HQhfBKrPo7L.DWFeKMvCbDtCgnO6.G3..n3MPx9gcQU84G6Q5F2IXU2k8TU.180AACCC4IVLO8s4xkpzRYRdMCy9N1EKMJdwDj68uCe3CKENb9m+4aovowMNbiabi..3TO0SU5wcxYXw6HzxNJsvyKdIvfqPAWYunU6Ppt38dETfhvIB1c3Ksf9RW5RwfG7fkQtfttd0hHDijgwMTfGxvbvSeR6QH6wJhVUPgpyvkKWvvv.KbgKD+i+w+vhdJU0A2YCTsngaHU7nOx0KiLd+O9i+P97qLA2f8CdvCJ+eUDxm3.dcCinU405yUu5UiLyLSzl1zFK0mnSFl63FhSofXz1Dqi0us7D7Z5L4z43YMXjSaviTmnoeS0QPNWxoSmHTnPx26wiGDLXvJz1.GDs2l27lAPIKEOUWlahVIOhKaRgJWvSeUxtGJ3HHmyxqAi1KgDGuHtjxjNc5D+u+2+CaaaaS134m5Ywx8v9oBfooI1zl1D..ZSaZCV3BWHLMMk43KEVoNb3.iabiCtc6F4latXAKXAxhRHPYeQCJ8.BEJjjQg6LJSSS3ymOow21y+33kC4nzQkGwLbCi4d4mOtexfREJT8EBg.974CgCGVRWSB.4oq8nG8nQnPgPspUsfggQ0lZPhSmNwQNxQvu9q+JRKszPCaXCgOe9foooz4ejLGdM9gTxgbXOORO3g9dU8wGET3XAZ8a.fm9oeZXZZhzSOcXXX.ud8lPTT3KKvzzDd85EETPAX8qe8nF0nFnYMqYviGOxMNrr.t7AdpYGJTHob3JaPYKPnPgvG7Aef7TlJQHB1NYG7S1UZ8ZSSSr5UuZYzvnqqi6+9uebG2wcHK8J0nF0vhd7UWgSmNQAET.b5zobS7IdN6QhSkw50TafL.MqrxByblyD..IkTRkYm1P52rksrEjZpohF0nFAe97ACCCo9MUmAUu3H9BZ9shr9gY2dU9gaQ94mO9nO5ifgggj1jrAOQP1eYEj9AjcE7Mct5PDjWUG70x000gSmNQvfAw9129vd1ydv.G3.Q6ae6ieAdjHNLiGIRD7HOxifoO8oKKN9wZCjGcYjCtHi3JnfBPN4jCRN4jQMqYMk6hBPwmNDABD.G4HGAQhDAojRJnl0rlE0wXNYqr.mNch7yOejWd4IUd1t.jTSMUjbxICfRZ.Z7fgh52TZdvOZQIGiw8dJOkIhl2+UPgpJfbdCs3qKWtj7gG3.GPxmC.jVZoIcNMoncUc5ehOdO6YOnl0rlnV0pVxEv4NCibTNu94XW1CWVnZgdENY.6e+62hLfZW6ZCWtbIqYKU0M3h36000w9129PxImLRO8zkaHPYU9GoWl8TRX+6e+vgCGR8spr.4XyrxJK3xkKjVZoISih3QJipPYCFFFvsa2RZTJUWhDIB9y+7OsTugSIkTPZoklzA1UGV+tz.eyonwExdn8t28hTSMU4FdaeS9qHbXBEsolllHmbxA0nF0.ojRJVpifkEXW+lZW6ZiPgBcRi8K7Hiau6cuHkTRAImbxvvvPpCaEA3YaglVQGpcABDPxWZOZaBGNbbY8kJavKGK7ZUU75P4SgxFH+CQQtGQele94iibjifm8YeV7HOxiXwwtko4LQkLLMMs7WgPHhDIhPHDh4N24J.fXricrBSSSgggg7ZMLLDgBERz0t1UA.DZZZhoMsoIBGNbItOGKPWC+9auskYlYJ.f..Bccc4+C.QKaYKEYkUVQseEKO+p5fOuIDB43u8wxJSDLXPgPTTahZuwy4F6y0lllQktNZ+N93D+Zouqz99p63nM9YZZVBdxZVyZJ98e+2ky2UWfe+9E.PbIWxkH+riEcUhH3zxbYzU1vNeZ4gL6PgBEUdU67yQhDwxmUZ72Q62qfU..gCGNjxJRM0TE6cu6MpiUgBER9+IRznkFH5.GNbHZaaaa49ywzzTbYW1kIFxPFR41y53ECYHCQzst0MgPT0S1X0cDMYZevG7AkPe5d26dKLLLj7dmrKOygCGhYLiYX4yhE8JKOfooonoMsohwLlwD2u2j9McsqcU9rNYD0t10VL1wN1J6lgb7ussssh6+9u+RnmgBwOvGOI4d1Gi4xAOZ1CFqOG5uU0rczd6eLiYLB.H97O+yiqOmDxDAl7JNchP4ymOYDWEJTHYXcNoIMI7ke4WJ+MMtwM1RQlMV18BxyvQqVfI9KuRVPAEH+L9o4jllFdsW60PsqcskWOPwg6ZU8ceNV.saB.E6M23UA2LdA2tcifACZYGGiWmPVztHXu9GDsHTjB8Sdp2ZuHNRu3oBKsqp7zIVHhemhpIx3X4se60IlBKrPY5SCTxBBqBUtfePnjnL2vOLBn0LH9y3wteRQ1azJVq7mMWFB+DPjjGPxNrmFMTAEm98b4KJTDLLLj5Rjat4h28ceWKyGzXJE00.p5GhBJDuPzzCzqWu..xBirllFV3BWHVvBVfkZBiBJnfBUWAOMnI4czmEIRjRTCkAJVdZrXeMoOC2ly.ABbRQz2dhfDRKpICbSJojjo2.QXP0Enkrjkf69tua4IrVCaXCwEbAW..JlHHVcHCmnixcXfhKnjacqaUds7iF2QMpQgdzidThvqjTlNQwnuxaPoxoCGNjgbpPjXTzhICNc4xEV1xVF1yd1ikTJsrBdHka2IU7hiJ4jP6mVh.E6nLxQZ7z.jTLTSSSR6oDlUryDG3.Gnj2yzzDO6y9r3W+0eEBgnJ+ILY0EPg5LsYF.HgZtgRQESSSrvEtPrqcsqR8vR434dCTrbHN3ohF2o2b4Rj7.R1gcYJNc5zRpgwkunbLFPpolJtoa5ljaHhPHvXFyXvO+y+r7Z3adCWmgDkMzQAEpJCtdfj7KpDiLrgMLKod2S8TOExKu7NoXy9TPAETfuAojcc7ZBLPQ9RfzilBTmXQFIo2MWuQZyHT52TRToupSzbZB4vhjRJIDNbXDNbXKDF+vO7Cnu8sux5mgPHve+u+2QRIkjEmRDKNjg6LLNgI2oXKZQKB.EabiKWtvkdoWJF4HGILMMw3G+3w8du2qk6wIKFjK9qZZF+j4HQpX1Rsmu5q9JbYW1kgIMoIIMxLdTPkIiXOZ0MNdzdA.YwBkFi3NJiCJZZ3ii7SMwS118znI7VWWG8nG8.qYMqAolZpvzzDtc6FuvK7BILzeJTrbX2tcK4Cnc+pxFD+X3vgw28ceG5ae6Kl1zlVTcf0IJ3F7QGND.PVToo0fhVQSlb7Euf6BTrCg4Q2FeCXTQXQQviGOnicriXkqbkx5wiggAdkW4U..rn.JeLKQZMLETn5FRJojfllFd3G9gwjlzjjaj5ZW6ZwG+werkLOPAETPgpyfmU..EomXjHQrDPG7S5YdcU9XAgPH06lqSyICGXEmHnR2gXDrWvHEBARM0TA.jGGsgCGFaZSaB8u+8GQhDwRQ.bPCZPVhRqXUYVdjLwiZHJ0S1+92O9jO4SjQ+DPQDjie7iG0nF0.ewW7EXDiXDXaaaaxTkjtmmL3AVMMM3xkKKEcWBIB8epvIejibDnqqKcpTjHQr3foSTv8RO4bVhdhSOvifCtvMdwTk98jwuTDgvoqCDHf79FOMZOQFQyIA.E6Tf10t1gkrjkf5V25hPgBgoMsogu9q+5DBGtbxNHG9Z2gN7c+pxDD+nSmNwt28tkozIO84JqfWPjImeSNhiONPNpm.cJAFshpLWFAce4a.yICEj3XAjQ0su8sGKaYKCokVZ..XVyZVX4Ke4G0T6OQXsKETnpN359vi5chuysa2XHCYHXpScpxu+wdrGSdHInfBJnP0UP9vfeX3QkYC5Tuj6XL56AhsfhfhlLBjslwS8aqNgDhBkg8zqiLLnV0pV..H6ryF5553a9luAW20ccHqrxRdslllXvCdv3rO6y1RMWgbDQoYzEutLwCqPxAOyd1yVZ.CcO+vO7CQyadywZW6Zw0dsWKBEJjrNHP3jkcnmFuowNJMd.pXNEbJMPN8hNdn4mRowCPBU3Q4k85PCmF6PG5PX26d2XO6YO3vG9vHXvfvmOevqWuHszRCMtwMFMpQMBd73QJHieRLwC20SFnuNVfGgnsqcsCSYJSAW4UdkPWWGO+y+73S9jOA974SdsJTwCRt.EMS7TDLQg9kN0dIZEtirKqf3Ssml1TDDyq4fjxN4jSNX26d2X+6e+HmbxQJGvmOenV0pV3zNsSCMrgMTtgQz3p3uN0ko9gxoNEI6jRqzV25ViYNyYhd26divgCiwLlwfN0oNIWif6vQZ7TI2PAENwQz3erW2bLMMwsdq2J98e+2wi9nOJ1291GlvDl.d7G+wqnatJnfBJTggnUivH32ue3ymOK1pRoRIPI8aRzfttNNzgNDVzhVDt5q9pga2tgSmNkmB0JXEIbNDi++0rl0DmxobJXEqXEn4Mu43ptpqB.EGYHFFFHkTRAOzC8PVbZPr5LrnAtSHxImbv67NuibGsz00wC8POD5Uu5E1+92OF5PGpkztiGcYmLUTd000wa+1uMd5m9owDm3DQe6aekQvWhvt7EHP.IyeZok1ITjDdz.2Ip1SORSSSrt0sN7ke4Wh4N24h0u90C+98COd7HcfpllFBFLnLRUxKu7fooIZbiaLtnK5hPO6YOQO5QOPcqacksa5HkOQHBapL.WFgGOdjKpzm9zG7TO0Sgm7IeR7EewWfErfEfq4ZtF4uSYbaEOnPyVWWGyctyE268duXLiYLXHCYHU1MMInh7O0NqQMpQb6Pqv9lyvq8WjRI+7O+y3+7e9OXwKdwXCaXCH6ryFtb4B0nF0PltjTDgQ+0kKWviGOn28t2nacqanW8pWn90u9vmOeQsVjcxJ3oRpKWtPO5QOvy7LOCF8nGM9pu5qvbm6bw0e8WukzIfuAFpwPETnrA67U1SeGZCCt268dw27MeCV7hWLdhm3Iv0e8WOZVyZVkV6VAETPgxSPx935ZPGLcKbgKDCaXCCuxq7JXvCdv.nn.6vsa2kHaBNZHXvfXCaXC3lu4aFKZQKBScpSUZqIeCpUnHjvYQMeWsc61MZaaaKVxRVhzYXzN2RNo5UdkWAMoIMA.E6PJZG9ikbr0NQEYHQjHQPlYlI9ke4Wj6j+ke4WNF8nGMLLLv8bO2CV6ZWqLhi3mBXQK8AqtBZ9JPf.3fG7fH2byE.Pl1QU1vvv.d850xAy.u30WVAeNlLD8fG7f34e9mGMpQMBcpScBiZTiBqYMqQFkZACFDABDPlFvjy5HmgooogcricfYLiYf+w+3ef5W+5iK5htHrzktTYJplHcJdVdiiV+jjCvmGum64dPG6XGgCGNvK8RuzIMiQIpfG4r+4e9m3fG7fkHLtqLA4vJRFO+TOKdP6P8c9IWoCGNvgO7gwq+5uNZPCZ.tvK7Bwi7HOB9pu5qvQNxQ.PQNRL+7y2R8FjRu5PgBgBJn.b3CeXLyYNSbm24chF0nFgN1wNh4O+4eRUJ6WZfFKnHI0gCG3tu66FcoKcABg.+6+8+VtgD757nJBbUPg3ChlbHRlFcPqnooAud8h27MeS3ymO3vgCLsoMsJ9FqBJnfBUfvdoZhxVsCdvChibjiX4jC2iGOwruMnqOiLx.coKcAyd1yFO1i8XRaLUNCqjnLOhvS4C58bEKKMkJseMjGS0zzv29seKV8pWskvqlJ9llll3Jthq.CYHCwRc6Bn3TUKV1ge9tVQ+kbxwy8bOmr3d17l2bLoIMI3wiG7fO3Ch4N24ZYm+o1MOhgRDhfmx6z9fFuBFLHLMMQMpQMNtuG75gC2XNZrmNEKAJYwONVZe.EaPZRIkjznT6Ql3wB7TZxkKWxn6fGMXQhDAOzC8PX7ie7xhmcZokFNxQNBRJojPKaYKQyadyQ8pW8P5omNRO8zQsqcsge+9wgO7gQ1YmMxImbvd26dwO+y+L1111l7YupUsJzqd0KzjlzD7lu4ahd1ydZg9hhZL9IURh.8W7.79AORahVT9kVZogW+0eczt10N7C+vOfUrhUfK9huXoizse8bPJoaebyNMYztOjyTnHI0tbwSTXWlBusVUYAMZbije61saYpxDKznD+N24FTDAS0VKdsJiOdUZi+bGVQxT3okLeble+h01N2I7.EswAuzK8R3Ye1mURiPJ43ymOz5V2Zzrl0LTm5TGjZpoJ+qooIxJqrPN4jCxKu7vN24Nw1111vl27lgggAb4xEVyZVC5e+6OZRSZBd4W9kQ+5W+JA+RzFypNChtgudQJojBdsW60v4e9mO99u+6wW+0eM5d26t764Ne7nchIyqCRTDuveFGs0n32q3YZCSxF3s230Z+1amwqnqNdAR+A5.zId5HynMVZebkSmPHQYrwtrynImr7Fb5EZLhuAkDb3vAN6y9rwa9luIF5PGJdoW5kv8ce2GpW8pm72Ssa97bzVe195lbcIsW+E4f1L8im0mJqHZ51xy1fnIGphLctowM.XoVFEud910A+3QGe6zC1+c1amz2mHkxX1OLWr+YkmHZx231GQiob8mOQ4KrGIT1aC1+eN3xs.JZcNNcYkAJMdfXkGgOVxkoQiWd850xIINIeJVt+Bg.omd53Mey2DcricDuzK8RnoMso31u8aW98QatLQpjl.XU1O8dfh4WNQzMOZnLaQE8f25V2JxImbjon.YTQoU3uc61MxImbjoxVnPgfggA9ge3GvnF0nPf.APxImLJrvBkJgZZZhjRJIbcW20g0rl0D0H9IV2cb2tcibyMW30qWY5noqqiku7kiu5q9J.TjWVe3G9gwV25VwO9i+Hdi23MJQgpK6ryFYlYlx6AIjIQgnpzbFPYAZZZXW6ZWPHDXm6bmX0qd0wLCUvfAQRIkDz00koSHWX.OclHmN31saKm7nGKPBy2wN1A..1wN1A9ge3GjQlEkZmjPFtCz3ofKcZv4wiGDJTH4eImg9ce22g68duWYQPL0TSE8rm8DsnEs.snEs.MtwMF555HXvfVJr17z2jTpmp2RABD.adyaF6cu6E+m+y+AaXCa.aaaaC8su8EsoMsAOwS7DnV0pVx1BY.eRIkDBDHfz.gpxf6PBpdqQy8QhDAadyaFqd0qVpv25KJAA..f.PRDEDUH4vx67NuSLgILA7bO2ygjRJI43M.j0ZojRJIIuJPwE+b606IhFMb3vVhnnBJn.K0eHWtbgfACJSos3wXOIO0oSmnfBJ.e629svqWuxHGJQ2oXDObjHQvl1zl..vd1ydvW+0eMpQMpQodRhwKx8ACFDFFFvsa2xn3wthZBQQGqzETPAvqWuk5b.4fNSSSr0stUnqqie629Mrl0rFHDEUzmon4ztSA3Jncr.QKst0sN7u9W+K32uenoogjRJIbYW1kg1zl1fV25ViZW6ZCud8JWyLojRB986Wp7yYcVmEz00k8Yhdb6ae6XsqcsHyLyDqd0qFae6aGWy0bMHiLx.Oyy7Lnl0rlxZRHOZonPyu5LBFLH9se62v28ceG.fb8iPgBggNzghoLkofm64dNjbxICmNcJqaGzXEsNA2glNc5Tt1Am9A.RmSxibbZNiudAcsk00iomKUDd+y+7OwO7C+fLp2hGf1.H9AKTN4jCb61MV8pWcb4YbhBRVaVYkEN3AOH1vF1.xKu7j79wqzd1t9jD+ejHQfGOdPvfAk7XjtKDsPkInMpKRjHX26d2nUspUx4xJJmpP5vSkWDe97gMrgM.MMMrgMrAr+8ueo9Utb4BsnEs.YjQFXcqac3EdgW.+8+9eWtlcMpQMfe+9k7T.Pp6CoGjSmNgGOdjG.Qd73Q92.ABH0yiJkFz7KUqAo1REwXC23sPgBYgt0vv.acqaEqZUqRtlejHQfOe9j55Td2FowpBKrPTPAEf8u+8aIHEhGO+fACBmNch7yOerxUtRYeMVb3A4bKRWW2tcK0SviGOPWWWReP8EhWkKOqxBTj532ueXXXfcricfUu5UWg5LVmNcJ4sJrvBgCGNfCGNPAET.16d2KxLyLsnSRf.AJQMg8ngHQhfjRJI4b.omNYakCGNPvfAkY8Boqse+9k54vyJL58DsQkk8MQy+BmH7BTTx50qWHDB4+GJTHrsssMo9cj9Kjcv7M.5XARG0vgCi+9e+uiINwIhgMrggBJn.zgNzA3zoSYZXRxSn1Phftg7M2f6Lzcu6cK8Y..rrgGGOaHdzdfkIDIRDgggg3Vu0aU..gKWtDZZZB.Dyu750q7+000E974S9dWtbI.fPSSSnqqK+eMMMgCGNNtdNQ6E0Vo6M87n+m+d9uwoSmBcccgSmNk2Cde2sa2k41V73E0ur29OdmihkW1mOb5z4IT66n8c76ezlWh1K97.u8PettttvgCGQscXu+Pu+n0lc3vgkmAml8nQ6QzP1eV1osb61sk2GM5e6i2GqmcUoWT+RWW2xXpa2tsPGqooUB9NGNbDUZc904vgC4qnMlwoUh1uOVnciGuRM0Ts7dtbySzWGK4.wKYDb5T5+iU9W6uRJojNpsuXQdi8W70Zr+642eheLZ8oRiNvNcUzjqYessncsGMZYd6mVWxd62kKWk34VcQ9vw5EsVscY2b4mjrUNMfc8.3iY7wU554xko4JN8BcOhFuP7nOR+eMqYMiqysQacI5yOQ4gi2ur2WIZ83oNNz5N1uuQ6+o0kJOzw538kKWtjyg100sh5EW9VznY3qIXeLSSSS9a3+Vhul9LOd7bL6a14Oo4oRiVph3kc8anwC6xrn9XkwKtLz3IcDe9nV0pVV9tim9qGOdJw5cz+y0SxoSmIL1lweQsW65kVQLux+K8hlaiFe6wqNbb5Y65gvWOwNMwQS+N6zgIpuhk1GoSG0m35uQeezli35sFqyy1kqPuHZNtrzJBZuXscSzIj7.N+7S7DOgHu7xqD9kJb3vmP9ypL6PLgPHLMMEBgPDJTHQvfAEgCGV32uegPTjCyNVu762uvzzTTXgEJBEJj3IdhmPN.PN.vNi2zm9zE986WDIRDQ3vgkuBEJj7E8Yk1yOTnPVZGgCGV7nO5iJIbF23FmnvBKTHDBQvfAE6cu6MpSbWwUbEhHQhHu1X4YWQ9hONweEOt2gBERLhQLBA.DKaYKSNNVPAETp+VgPHInMMMEFFFRZmBJn.47oc5qXssYXXHBEJjXwKdwB.Hdu268DgCGVje94Kmmn6anPgDFFFBSSSQjHQjOKSSSQf.AJw3Yf.AD+s+1eSnooIZPCZf3W9keQdc76AQiw+L6unqwNcbnPgjzT986W1Vm9zmtjFckqbkBSSSQAETfHb3vBCCCI8ZkMcW73UAETfr+aXXHGic4xkX5Se5h7yOeIsTgEVnj94e+u+2B.HlwLlgbtklO862uE9f.ABH+LCCC4KZtmt2T6fjWwuuT6yvvHt02o4eMMMQm6bmEACFzBeRY89y4krKSHdHefy27jO4SJzzzDyadySTXgEFS2ehVlded4kmjVmjIP7wABDPN+Fqz+zbZAETfX9ye9B.Hl5TmpPHDh.ABH4Ko4XZskXc7QHDhtzktH.fHkTRQ7K+xuHoUH4A76CQSQumjMDM4114yow.hdcVyZVBcccgKWtDye9yWDLXPQf.AjzOUWjObrdkVZoIdq25sj8a958gCGV7JuxqH.fXJSYJVnyH5BZLlniCFLnbME5yoemoooHXvfx2SywACFTTPAEHJrvBki8wq0dI4.gCGV30qWQFYjgTlj80rNQdwW+knYBEJjnacqahAMnAUoO+Rz+CYHCQbIWxkH0mHd0+o4SRW.hGi9L5YPumlO372UluHZi69tuaA.DyZVyRDJTHorxx6me3vgsL1QiUqXEqP3zoSwN24NkeFIukl65e+6uPWWWrxUtR4uMPf.x0WI8K862ukmAWGbxvH5dxu+zbGWNJMGR5ATdO9XXXXQ+FhlyvvPnqqKl5TmprOPeOYuDWdU4Y6izA5LOyyT7zO8SKkADOV+fzEWWWWzktzEQnPgjxHKrvBi41Wt4lqEc4oumZ6zX3QNxQjzAUDiek1Kh1rvBKTTu5UOwq9pupEZ8Jp1Aw+Pz7ACFTzt10NwccW2kE9EhNjzYJVlen4ER1LMePxF3xWowBt9z70s45bWQLtbzralqCncePb7XaM2uCjMFj7nm7IeRA.De5m9oxweZtHV8uBYORAETfvvvPzyd1SoSlZQKZgX6ae6R5vPgBIuuzZXU17Gb6oI91fACJd5m9oE555hErfEHHv4mHa.NdQYNea3oFGUqmnv.MRLDRk7799e9O+m3cdm2A.EWSd.fLsZz00wq+5uNF3.GXTCmTwIPMbfZ6TnH9Ue0Wgm8YeV3zoS7xu7KigNzgJ6W7T0iR8DJ8Kn1Jcj2yu2UlHZ0xHQbLbmi7WonFkBANc5TN+yqcNGq1WxImrk1IE9nIkTRVdNT8YxdNuerPnPgfKWtrTefb5zojtjSexumT37RsQ5zLTvRktAO3Aiu+6+d7u9W+K7pu5qVhvblnYr2GnOiWaHDGkb4F.xTvgpYYlllXPCZP3ptpqBImbxnScpSXO6YOnAMnAVB29DkZjPYETZDP8eJsZI5Npt0EIRDYpOKDBbcW20gG3Ad.74e9mia5ltI..Y5GRgDLM9GM4TzXIE19TpJPoCNARNG+dDKx9hEPxGE+UZ3PyoF+0gEQYEGqPyOdD197Ze.w+3wiGKxIKs1G86c3vARN4jA.jovJud9Qi2Q9qTdIVn+oqIojRRlZsD8lGOdjWGW9C89XY98ttq6Bey27MXXCaXXhSbhVVej5SQ9qzQfKyg5679PzpMfTZYQxknzKvzzD2vMbCn28t2H0TSE8qe8CacqaEMsoME.vxg4Q0Y31sa3wiGYZOR7xDc4McS2DF0nFEVzhVDt4a9lkeGk1jz7CQKvqSfTpyyqKa.PNWH9qz7vsa2VFqo0JiW0.EpFzQo.FcO4oOvIJH8hHP2acccI+RkI3EHXZcARFAm+8DEz7OkZc.PpaKeNkevLvqsnU1iOz5aTp5t10tVbC2vMTgmRVj7RRud9obF0V3yWNc5D2wcbGXdyadXaaaa3BuvKrDkFER+P60RL9I3KApbKPWC8rHdTZdiyeWQUGcrqeC8LoRE.MGxKCMz5mk2sOh2mmhgw60Mn5sKQWPO2XQGARFZJojRIZ2D8MWNXMqYMkWS7RGsxBH5QpTfv0ijulR4EnwIZtkn0b3vgrrhv4M4oAtOe9hI6635MUiZTCKqSxkSR55R7ob9XRdQkQJny8qPzfcaqo4zXYtirUgzaiNkvMLLPvfAgllFRN4jsHajSOWZ7+TamlqZe6aOV7hWL74yG1zl1Dti63Nvrl0rP5omtTdCYuaEQJiGKf5ijNbjLI5uQ9qxdDW2D.q1yGqnLKMvgCGVTriWP9hEgMTcTZDiXDXJSYJxO2kKWVNtzCGNLdpm5ov+7e9OkWi8N7I5DHwzcnCcH7u9W+Knoog5Tm5fCe3CiwLlwHyWeOd7HOI.AfToX.fssssgm4YdFDNbXnoogt0stgtzktbB0dhmHZLlwSBcRgl7yOe.TTAKlnGhkmCwvxIho+Z2wWbl+Xs9IQKdSySjCkHEd32WhYyt.a6sQhdoKcoKnScpSX3Ce3VdlbCxomA4rPhVhWux3yQ1cpKUuC3KpvWbwzzD25sdqx5bV4QgEtxDz7A2QDjiVoZwDPwzCz0KDBb5m9oia7FuQ79u+6i24cdGYsThuHMoPLudjvcvKWwZRgAtPY6Be4G1BwCiN30AGx3a.XwPhxBNV7nwC4D7hDLwCRN7NV3gIdinUyan0b3iEGuNYfpgEjwzjBezyijAvqUYzyLVvkbIWBZcqaMF9vGdIZa1cjG.r7LI5S67.bPziz3B0Wn2mbxICCCCb228cawwj7Mup5LhDIhzfKtSpnMt4TNkSAW20cc3C+vODG4HGA0pV0xhQB7M.wtAaQqvfSxJ35+P53P07D56hGJ2aeNjjODuL1inE4xAIE1CDHPY99WVAWFdf.AjNoDH93PetwJbvGaIdN.qEf7Dk0go0hLLLvBW3Bwy+7OeEFeuFq1tPsEp9PQ7ODuH.rnuSG6XGgGOdvbm6bwfFzf.Pwzij9Nz7CoqD4zHpl3P5C4ymODHP.4ZO76CeS4iWqaGqfzifnS3iUgBExhSJ375wyM09XAROYRNGWe33wylquMoaFQSDK7O7MDilese+ilt9.ktiNpH.oKIMNPzhUDNCC.xZAFe8JtCm4zf1273Xo8w4i3a.MOPJ3yyTMgiB9D5uTMxjPEE8+w6yHVcDFAZ7gSex2PMZ7lGLFzuIVz00tMzz5TTs2ZwKdw3QezGEu9q+5VzcgnKprcXLeyN3xf3zU1CDAdvRb7h3RukbXEYHAQfGKBzb3vAdlm4YvDm3DsTf3B+WERcRXL.P+6e+kLBDC0wi2aOZf18u.AB.gPf5Tm5.CCC71u8aWhEe4Nmg2d24N2Idy27MkKxmZpohK5htnJcBJRfK2oNwSgH7hanI6DwHVAeGnEBAxJqrvm9oeJxImbP6ae6Q26d2sXzgooI9oe5mPqacqiIC5369Hsy+GMuqSsC62StiP3FFby27MKoAsabO4LL6OCdzgD4uJhg74DR.FMdPBv3sI52Qzfu268dkXLMdE8AU1fK7mO1ATztdPEwdZblDhSJA0oN0I7ge3Gh0t10J4G4NOjeO4JmWXgEhUu5UiUspUg5Tm5fALfAf5Uu5YYA7Ce3CiCdvChy9rO6R3zTthikEPJHPJHRzKkGNyHZJMVVAedIu7xC.E6H2XQ1HI+hKG9.G3.3S+zOEG9vGFcnCc.WxkbIRY3DM++6+8+PaaaaKUCa37XzhsbdNxfNtw.j7nXwvoq65tNoxtztuQ2GdAcFvpys3FFPFowkKvud.HKT1jSS4q+FIRD7FuwaHkSmHXjdEEHmYvkuxc9sa2tQ6ZW6vLm4Lwu9q+J5PG5fbNxtSuICUneWnPgvW8UeExLyLQ5omNtlq4ZPcqacsL2kSN4.+98iS6zNMKJ5BDebXCecG98NdYPu8MihyKFOh.qxJ3NxfhFPBwiw.h14m9oeBKYIKA..8pW8Bst0sVdMTwZeqacq3+6+6+S96RD3y3aPjSmNwF23FwZW6ZQ6ZW6pvzQvttKTDNPGfVj7tvgCawYlImbx3pu5qFevG7AHmbxA0pV0RNmxc9AwWRuO7ecfyLu4MOrwMtQzpV0JbsW60ZIpv0zzvV25VQ5omNpScpSTcTREw7mc8a3yGTjkxMLkn2su4skWfnio063NuLd.RNKu.oSxriU9W2tciBJn.rjkrDr10tVzxV1Rzu90OYzMQieacqaEmxobJxnDKQP+Xd+j3SqnbFF8745gw0OHPf.vmOeVhrSttxwx7C24d79oKWtvF1vFvRW5RgggAtxq7JQyadyk5EATDsw5W+5QFYjg79QQhaEkc0QSOatNaQa953U2b9FpP5VFIRDYPlTiZTCK1SP5BFqzub6Z27l2rjNizW3cdm2AWvEbA3Vu0aUtQfwZ.mTQ.Z7jG7A7fUgzIm13ZBmHNzqL2ioFDexgDdFKKnLgILA7TO0SIIfnNMsiV23MdinqcsqXnCcn3fG7fxeG+4UVDdvYraXCaH9oe5mJQ+i+L1yd1CZTiZjEkkMMMQO5QOvBVvBrnLcE0tLEKfVrwdzHUVE7RF4UPAEHUR93MjEIlyCdvChK9huXrsssM.TzX6XG6XwcbG2A.JZtZ1yd1XfCbfXO6YO3zO8SOlt2.vhgh..aYKaAKdwKF+5u9q3BtfK.W7Eew3rNqyR9a3KTwEVQFwFIRDDLXPbnCcHXXXfS8TOU4yjV.g6ca6vtycruqk7eG0tsyzyWfhhJsCcnCgfAChF0nFUsIBPrKegnkoc9k.Ob9ow+F23FCCCCje94aIRh3QSF.rH+IRjHXJSYJVhF027MeSrxUtRjZpoBGNbfrxJKbK2xsfK6xtLLhQLBIOucCGiGfbDDsn.0tiGysbZMtLs308mFG3m7d7TJJV.OxuxN6rQ25V2vu8a+l798tu66JiRRMMM7AevGfANvAhCbfCf5Uu5ESOCJxvHd7vgCicu6ciO+y+brm8rGbNmy4f9zm9fS+zO8RDolGKjat4hbyMWnqqi5W+5aYLg3OowY97JsKX7nUhfcCLc9Wm5g1+dthMzbqggANvAN.LMMQcpSchKocahL36pH.rrgFzeaVyZFb5zINvANfTtbzh5GxnM52M9wOd7fO3CJid725sdK7ke4Wh5W+5CcccjUVYgAMnAgK4RtDLpQMpxEYw1i5XdDOFuhPpnsVlPHJ0SH1JBPsqfACZIUJhWxd000wpV0pPm5TmjxEevG7Aw28ceG5PG5fj+azidzXMqYM3a9luA.kbm4qr.OxZI4Ke1m8YncsqcUXNriqmCO5RzzJNUwHm4Simj9NsrksDZZZH6ryF0pV0R1OnH8hughjy0b4xEF3.GHlyblibL3y9rOCyXFyP1m+ge3Gv.Fv.vzl1zPW5RWjafV71g0wBrqeCYvqwecRmy0mlZ+UTarAm2mVKIduwYzlNQykGOQHEIe51tsaCezG8Qx07t0a8Vw69tuqz4Iqe8qG8u+8GSZRSB8nG8HgH5W.fkMXka6azhH9xqmOW9f84ZJJfI8hOd4GrWBHH99u9q+ZboW5kJelO7C+vX0qd03u829a.nH44O5i9nX0qd03a+1uU5Of3Yz+VVAWm4il9Ywx8fuINDsOoChllVIroguYcwx8m3q9i+3OvLlwLjovNIuj1rZdjDmHs4oQyeE7LmvdPkPe9IB+cYlhhGwB7IINwJwTwMvWHDXdyad3Nuy67nRX8XO1igoN0oht10tBcccrwMtQKFPPSnbvaCwBQC2q0QCzj.sH.8d9hl71cksBPQCZZZ3Nuy6DW4UdkkKBQBEJDxKu7fPHPpolJ.JZN+3crXNyYNXaaaaVX5efG3AP94mOLMMwDlvDv.G3.Q6ae6kFWVZvdMPYRSZR3e+u+2nYMqY39tu6Cu8a+1XHCYHnyctyHyLyD.V20F58jfpcu6ci24cdGz6d2ajRJofS8TOUbVm0YgzSOc7bO2ygbxImRvKPzr7W1A2C37Oi2G3BBo1ygO7gw6+9uOt4a9lQxImLNsS6zPSaZSgCGNvPG5Pwl27lAf0ZwDgXUnZ4IrG0K.vBuFerLb3vRkB4x.hVTQPeV5omN..NvANfkcUf.oH.PwiG4me9XDiXDxqQWWGaZSaBe1m8YPSSCae6aGCZPCBewW7EnUspUVbZtcG6GO.WFC2IJwy4OMMML5QOZbwW7EaIhMKqfTjB.3O+y+DNb3.olZpG2xFnwy4Lm4fsrksHumZZZ3e9O+mnvBKDlll3cdm2ACZPCBcnCc.0t10tTuu7H.i186wO9wiW60dMzzl1TLhQLB7FuwafgO7gi1291iUu5UK+s7wGNc5N1wNv3G+3QO6YOQpolJZXCaHZXCaHpScpCd9m+4wAO3Ai5ZWQyHfRiNhGN8.EqnHe7Iu7xCe7G+w3ltoaBd73AmwYbFngMrgnF0nFX3Ce3Xiabix6G4X8DEXWV5Qilj1Aa9t2FseOesa5ZSO8zQjHQPVYkk7ysWyNruNet4lKdfG3AjxpBFLH1zl1DVvBV...1291Gt9q+5wRVxRrrC2jisOV8kiWPNAfLfNdZnm80Bo9e7N5Jru9zQarw9ZBzukbZb71HRSSS75u9qKeOw2N9wOd.TTpmbe228gW8UeUzidzCKs8DEcAIConwmoLkoHcRO.JAOCPwqAGK2a9uk+dRFD841+qPHrThBnM3ie+pScpCDBQIJUIQ6dRzle+2+8XNyYNVnam0rlEV+5WO..VwJVA5W+5G16d2KZcqaskZTLeSFqnl+rOVaWGG9mehfiEeEmOJRjHx4L65TQ7UTjHGOGa3NDh6LrXo+poogLyLS7QezGYojR7du26gMrgM..fktzkh92+9icu6cKqgl7MIpxFDcK2or.wW5O61qZeMQp1Xy2bEdagGooGOvt8Jz83Mdi2PdMzyYricrxMu9ge3GFu5q9pn6cu6..xMUk5KUTNCy9Z0zFunqqiG9geXbgW3EZ45I8P3x9hk6OAd+JmbxA..okVZRY31CDmX49SzUTIoxdjXdpm5ohq65tNY+jx5oJaGNBXc8DtbR6aPmccDOQkUFWhPLdM4fCpy30qWoG4ocacsqcsX.CX..n3TUgp2G555XZSaZxBgcCZPCfoooLxg3dwj+rncU2tSqJKfHJ3BF3FVSDaABDnDBRRD1gPdDDP0vGdcLId.ud8hibjiHS2HZNOVAMVkat4JaejvyBKrPLm4LGr7kubLyYNSnoog9zm9bb48ZZtxkKWXtyctXtyctnqcsqne8qe3a+1uEyYNyA+we7GXXCaXXwKdwnt0st.vZHWFNbX7hu3KhG6wdLKQEY3vgQnPgP3vgwnG8nwRVxRvm+4eNRIkTha6hg8ZvAgYO6YiG+webrksrEKQQB8bmxTlB9fO3CvF1vFPSZRSrraPIRd+myuRJCR7YbE446ZA8dtC4Iv2I0ZTiZ.GNbf+7O+S4ygmVkz0SeVjHQPgEVnk2SJC7e9O+G32ue7BuvKfcsqcA.fV1xVVhzmiuqtkU9eZW4nTmjtew6c3jRyjn4nlxB3xYJrvBgPHh4BpOPwJhQ8+bxIGYHky4+lyblCVwJVAl4LmILMMQu5UuhI5aNeLg4O+4iErfEfdzidfd1ydhu+6+dL6YOajc1YigO7giUrhUfTRIEKQvEIqZricrXjibjRZGfhc5ZN4jCd7G+wwxV1xvG+weLpUsp0w4nYIAO5Q4QPHsS9yadyC2+8e+R5U.XYcpINwIhYNyYhu+6+dzxV1RKolSzbReEM3N+kyWCXsMZucRxRH4J7hyK+dCTbjhkc1YK+sGM9LZ7l3WHYOTTEthUrBoLhrxJK3xkKzxV1R48kKqHdI+kK6gtuwK4C7wURoSMshJf+w60OnwbZ8NpODMcLowdhtmVafK6Jdr9qttN1yd1iE9YgPfu8a+VLqYMKL1wNVrl0rF3vgCbdm24Y44VYq6G0VH5AR25rxJKLm4LGLjgLDKiQlllHTnPRc2hE8WI5exwVjrG6ELeN8IORt4F+ycTMw2Sav5ANvArjhLz8vdM9zoSmH6ryVJeCnX4DKXAK.e1m8Y3oe5mFtb4Bsu8sGokVZkv4tUz5samO0djNbrzuoz.mm09mQyUjQ71SGNZbf3G4FIGuVWfzOfWmxHGdGKOCtcCbm0EIRD7Ye1mgO+y+b7jO4SBMMMzgNzAbVm0Ykv3HLfhWOgz8xtclkUYrbYn70y4Q7E2YSTj0XXXHsamtd6HVlir6jOp87G+weXYS9MLLvpV0pvG9geHl5TmJ9xu7KgllFtfK3B.PzOv7pn3S47P1Onk3Nphmx2z3aosFLw6A.KQFqtttLHSpUspUI1DuX0gUzXzO9i+HdzG8Qk8CZMAMMM7RuzKgS8TO0RPukHr9Uz1LNRevnIajecmH5.EWjpYOrH4dwj9NtyvxKu7vvF1vj+FCCCKF8unEsHLnAMH4u2kKWnacqaXMqYMPHDVNcwHDsHLIdnvF47H59uyctSY+jTTC.HyLyD4me9RhM5ZprA+DKguSZwKkYICIOxQNBb3nnSANhfMV6+j.wK6xtL4IFHWHzsca2FlwLlg7yunK5hfa2tQgEVXoduo4GtQuSZRSBKYIKA28ce238du2CcpScBlllXcqacXKaYKVDPATjm56ae6KF8nGsbdmnmArdRNQFkSedYEjA6DBDH.762OF0nFEtwa7FwV25VkiKT6gKzLXvfXzidzVZKTd3mn.RoK6KXyCea9IwGAp9gQf6PJfhDHRm9PG7fGTdO3EbQxPKJ7vc5zINsS6zv0e8WukqWHDXtyct3ttq6B6ZW6BZZZn6cu6nQMpQkn+PFCEOVPQWWWVO536fW7ZGNowOZwXud8ZQY3xJ32CxgCTD4Dq+dtiK5QO5gb9kWmIGxPFBl9zmt7ZujK4Rho6u8cKidduxq7J3y+7OGiXDi.yXFy.ssssE986G+xu7KX8qe8kHxByJqrv0dsWKF4HGo79dz18uku70mr76V...H.jDQAQki4N24FyQgwwBTJxwiJLc8hRwz6+9uebsW60hcsqcUh58D01b5zI762Odhm3IrH+HdGE.mnfVavtCFrevjvuVtx+lllxSYR90Q5NHDBTu5UOnqqiCcnCAfhGahlBUzyKkTRA2xsbKREhCEJDDBAl4LmIt+6+9wu+6+NLLLPG6XGwYdlmo72SzWwJ8eoAZcGpuWdEo5bGanoUTZbDOVCga7qgggzwJ.VSuNmNcJKvyz5uTem3E4QqGcuiGsua9luYKyWtc6F6XG6.2xsbKXUqZURi46XG6XB2Zqjr.9gAgggAdjG4QjNYhmZdjyvNd3+I5AZdh6bE65.QiOTcwzda098kmwA7xF.wOyqivzytCcnCnwMtwknc9jO4Shm9oeZ.TzZGWy0bMx1scGFPO+xaXeMb957.nT0uoz.m+f2G000svCYOx7301MZdhW1ChGqcATrw4TJhwmiikTxVSSCspUsRV+oHcE0zzvS+zOMdpm5ojimWy0bMV1rxJ6M6AvZzgAXszcDOjiPadCuuRxq30mT61qoqqC+98WhSTT9gSVrN9Qa.OPwzy21sca..Vr0Z6ae6XvCdv3K+xuT9aae6au7dv4GqncZMQOR5eSE7e95sTvfPNJKVbFCw6wK8MDxKu7fKWtjOGRe2imM5QSSC6ae6C8u+8W9dJB1b5zI5Uu5Etga3FjeGfUd+JavkyDsMminC3NAijabhv+DWRYRtGa4LJTCKRjHVTz4sdq2BqacqSxzy6nqXEqPlhj7cMhRWkbxIGYnwREvdfhixn3YjOQ2WMMM7e+u+Wbe228gq7JuRYeldoqqiryNazidzCLgILgxzDR7FbkHoBzIP7yYczbZf.AjG24.ViRliE36V94e9mOl4LmIZQKZA.Jhd57NuyCu8a+1XRSZRR5KZwO5X29XAdnU5wiGLxQNRbq25sJo674yG5e+6ujlYm6bmxEPMMMwgO7gw0e8WOV5RWpj4iuaz7iic5yV0pV0w2f3w.7cKj10lm8YeV7RuzKYIJHHE637MD+wbm6bQd4kmzoSjf63gACwCDHP.KJrCTrfPxfS6JzSJQQNDIZJNxitfTRIEKxKHkB3xt3Eo2wMtwg68duWKsoa61tMrzktTjQFY.GNbfd0qdI+sbkUimKTqoUzooIIGxtwfwymSgEVnLJthmfteETPAkv4DwR6xzzTl18YjQFX9ye9nMsoM.nnw6y67NOL4IOYL9wOdYjCSQkSoAtQ2jhLibjiD228ceVLV5lu4aF.EslVVYkkEY+4lat3ltoaBKZQKxx8z9bEozottNV1xVVba8AheldVABD.O1i8X3sdq2R94T8nfnYrGAjyadyyRDDSxJprAs1.WwHtxl146iVj4DNbX4Zfz0w2QbZLfVOgjwX+vOgtdfhbJxy8bOGtm64dr3Hp+w+3efu7K+RjQFY.gPfALfAHURluyl1Sk+STv27QtBsj9PkUvGyHDJTHK7jkUPiKNb3P57e5EPwoxEYX.YfEOcJH5D9tzGO3uzzzvPG5Pwq8Zulz4LgCGF8t28FyblyD28ce2PWWG8t28F0oN0QxKx0csxFQKaExN6rwzm9zkeOM+RGrTGOzlTDxCTrdebmYPNofbNMwyQe1QywGbcTHYw.PV+aHdI5y862Oz00QZokFlyblCFv.Ffzn+5V25hm3IdBoLZGNbfNzgNHae79aE4bG2gjz6o1.Querzuoz.OiVrGEdjLdtcZz7G2QMz5BDOX7Hxk3fjmvKf6.wVZYGLXPbZm1og4Lm4fd26dKksWu5UOLlwLF7oe5mJoOtvK7BSHrIiC9AKAMm.XMvRJKfbRMAxgzb8Q3NIEn3C3AReYZiJHaJ.N9RoZ5daZV7oa8.G3.wq9puJpcsqsby75d26N93O9ik5c2291WTm5TmRLNvk6WdC9yg1HcxAVQyQejtIwptSTjlS7XbmTlat4JGy4Agwwa.2LlwLFru8sOKa5K4KfG3AdfRnekcmzVYhiVekqqCcpECTz3mcaIOdPbwyQ9862xDnKWtPf.APRIkDJrvBkGs7gCGF6XG6PF5dDnN8rl0rv4bNmCNzgNjbgBhYgJZ4aZSaRVnroTug1QPfhEh52u+3xofDwD+BuvKfku7kKWPmTfibRge+9w5V25vHG4Hw+O6ckGVTU11+2YVggMUbEw0xLqHWyBSe0vbKEyVM2JMWxVrbK0L60LSM0L2JxJWJWPyz1TyLszLeKq7qESkxbWLPEDksYed99C59g6ygQX.FfA0eWWbwvvLmyy4Y4deoScpSH7vC2msRbYIn4iLyLSXylMjQFY.a1ropEqVZ.QH+7m+7HrvBCojRJRitjYlYVjqA5zoC1rYCgEVXvtc631u8aG6bm6Dm3Dm.gFZnHxHiD50qGqe8qGJJJn28t2.HulafEKVJRujQJ1jSN4.61sia3FtAjYlYhbxIGDTPAIIHQJrDQDQf+4e9G4Z6a8VuE1wN1gjPN2Bz.PEgLRoKWtbgScpS4SiuhBbhUFLX.ae6aGyXFyP99TZCvijQhXAEB4Nc5DYjQFvhEKvoSmpXBVQCRH4LyLSUo2AOknsZ0JRO8zgQiFQngFpptGz4N24j+ubxIGICEpKTQQzoACFPxImrJOJSo4DUTr4dRyiGO34e9mGCe3CGYlYlnN0oNx1PcJojBb4xEhIlXj6UHudSL+oHXozd9hm1G1saWl5mzyXosV9PQFGWnkye9yKET1eEkid73Q1XCRM0TAfumRSjv7W5RWBd73Ast0sFabiaDokVZHrvBCUoJUAJJJXiabiPQQQJXLQKpv.obFO8dttq65PJojhTPUylMibxIG4XUHD3Lm4LRu6Mm4LGUd0jmNODcE99.hQMMOTZ.wehTrzhEKXm6bmX1yd1.H+HikFOj.Vz9d9m47m+7xhqJwuthVodxS0gFZnRi8Q78yN6rQPAEDb61MrZ0JxImbjJVERHgHadCYmc1H8zSWtml38QmWSN4jAPdzZ+m+4ejuFH+8tlLYB1saWJrKMWOkoLE73O9iCa1rg5Tm5Hijo+4e9GnnnfXhIFblybFYai2pUqxtrm+HJKn8ujLVYkUVH0TSs.Q7Qo45S2CsNZwiGOk58vTFCDd3gCqVsJcxAeOKYjLh2FEcQ70nKdwKBqVsBqVsJaNB9CdbDux9129hd26dizRKMDd3giHhHBXznQ7ke4WBc5zgNzgNfyblyHuu9aGVTZF+zbL24Od73ASXBS.su8sGMrgMTUZNdoKcoBjAFWNPx6avfAjat4Jk2hazFdFc3zoSDbvAKit+KcoKIcxnViRA.bzidT4Ykye9yqZdkdNHijZvfAjYlYBGNbfpW8pi27MeS7hu3KBa1rg5Uu5AKVrfe+2+cow0pQMpwkkGAc1or1o.Tp+SxuPMV.f7LNYQIeSQozGM+PxKx6lkz9TRVK5rEeNknsnnjWZRmYlYJcHDW47RJ3Qmjc61w4O+4gCGNjct0hR+I54ppUspXoKcoH0TSE4latn90u9vfAC33G+3RcEqQMpARN4jgGOdPHgDBxM2bqv6zjDuXf7jwK2byEojRJ..pZfAkT3xkKDZnghKcoKAKVrHiNHf7MTNmNEIKrc61gUqVQVYkEtzktjbrQFsLqrxBAGbvE49OWtbI0yhGvLJJJn+8u+nW8pW3rm8rnwMtwR50ewW7E..3Nuy6DojRJpjil1Sxo8WVCR1dRNPJpsRO8zkYG.2oXjby9x5GOaRrZ0JBIjPPVYkEBIjPPlYlIpYMqIN6YOK.TquhuZWAh9BO.UDh7Rw8F23FilzjlfTRIEXznQX2tcXznQnWudoMBBDhRL.HGSj9aji64YPk14kRDsIQoDd73QHDBwi8XOl..BiFMJ.fPud8p9a5GEEEgNc5juV6+W6OJJJp9a92kdOCFLHue76I8YKs+nSmNgYylk+McuzoSmbbXvfA48VQQQnnnHLYxje496O9QQQQNFCN3fKv7p+35SyGz7NeMov9dWtwh2de5Yf9sur1w+75zoS06YvfgB89n84vfACBc5zI2iw2qwes+Zu2kadgOlzN9tbic5050q2uu96O+gulncujd85Uc9i+7YwhkB7rRyEdi9f18XZ2uRWe9Zo14Ns6eoy+d65UR9guuJhHhP08hSSpz7iQiFkySbZV95Yrh5G95UPAETwdtwaeV9ZM8CcsKNicsqsZmCz92z5M+LHmmG+rO88n+uV9T9ie7E5LbZ+Wt+O8a+83yes2waOq7y17OGeMpvtd7mWtbEbdXz+S68les0JSh1wiISlJv22eMOq84LxHizud84yO7eJqVmKr0XsuG8CmlKeN2eQ+U635xIypd85U8+BzNKwelHY0JryJkDYDn0hK2dDZNQ687x840tVqUdauQevjISp1S3M4B3zD81ye4k7QZoyvGiFMZrHkuon9waeNsxz3s4dR1BNMRR1W563Olinwhd85EgFZnpt29xyn1wtQiFElLYp.qodSto.oymDcCRlgBSunRy0m+bSzA35t5suKWuQuc8708fEldyZWO3W+K2Xu7b8gVWzoSm70WtwB+yVTWWs5NPum1yW76Uw0tBZ0AzrYyEpdOAR+PiSSlLUf4Uc5zIl5TmpH6ryt.1kxsa2kH6YoHDkNSrJ9Wqzkd5oK81if4UUxpdJJJ3bm6bHt3hSUHBBjW2jbnCcnxz5fhRDJpM94e9mwJVwJvt10tTkNhFLX.20ccW3dtm6Ast0sF0pV0RFhuT9t6O7RWFYjA5RW5BrYylLBanno..30e8WG8pW8RNtCJnfj02n.gHDyoSm3kdoWBm7jmDevG7ARuz4O7PsPjWcpJlXhA+m+y+Au268dRu.3KdPhliHOWQoRBENpAETP3jm7jHt3hCtb4Be5m9onUspURukUTf1OdfCb.zu90OLqYMKzidzCoGzO8oOM5YO6oLhfdzG8QkemsrksfwLlwnxiuz9c.u2oOz8uoCU8qe88KoOm3e8Huc61wYNyYPm6bmKv8mFW7nNf73K.vRW5Rks3X9yg1Pksh.TzTQdMVmNcXG6XG3UdkWQ1M9nwYvAGLrZ0J.xqF9jat4h4Lm4f3iOdYjNQdakhrf0st0gIMoIg0st0ISQhryNaDd3gKKR8j2SnHtihPMxiPTzW3wiGYGHbvCdvXhSbhROaSysAGbvpROsRqGVb4xErYyFti63NPyadywJVwJPXgElLpbKszWnHHHmbxAKdwKF+zO8SXiabivjIS3RW5Rk5nrknC31saby27Mia9luY7Ye1mIKdyE09OZ+J48GxSQh+MEnI5GojRJxN9yW+0eMpacqqLxRJL3xkKomSO3AOH5Uu5EVzhVjLJyb61M9m+4eP7wGOxImbvjlzjvPFxPj6Y97O+yw3F23jzH3QgA.75YLc5zgcsqcI4WUZ.ctfhb0idzihd1ydJuuz9VZclROHxKvD8127MeSzqd0KXylMX1rYYTHTVGgDEEDrHNmGUy6e+6GyZVyB+7O+yxOq1HwC.H7vCGiabiCCbfCTU5GQzl0oSG1xV1Bdtm64vZVyZP6ae6QVYkEBKrvT0DKnqKQ+mhjC97JQqdNyYN3sdq2BCZPCBuzK8RxnFihdMhll+nNhQzXb4xEZdyaNhN5nwW8Ueke0y4b9LDM1G+webDUTQgYLiYTpt1z4W61siPBIDX0pUYzmLhQLBbK2xsfm4YdFYpVmat4hfBJHnnnH8n8t28twXG6Xwq7JuB5V25lzq6z+uzN937FzoSmrVGsu8sO7nO5ihZW6ZictycpJUSzoSGxM2b8KYoPoc76wiG7Ue0WgwMtwge4W9ELqYMKjXhIJOO74e9milzjlHi.Z54fhHxBCTza3xkKrvEtPbvCdPLyYNSTspUMI82fCNXjUVYIythie7iiwMtwg10t1IkuBPcs+kRUqW9keY7ge3Gh+3O9CUErZtrkz5Mc1jK+n9+sbKP7wezG8Qwd26dwq+5utLaCn6OOx44WuxRPxgAj27cDQDgLRGaQKZQQJeSQUmsn4A.0oiU1YmMxHiLP+6e+wy+7OO5YO6oLpV.xq9EQoroMa1v5V25vZW6Zwrm8rQyadykzzJs7uDBAxImbPKaYKQ6ae6wxW9xkxm5KQBihlr0gVyn0u90u9g8u+8i4Lm4fd26dC850KydHdIynhBTz5kc1YiN0oNgQLhQfG6wdLIMwRq7cD8KEk76N5icriE268duxmcZ9l2wWW1xVF93O9iwBVvBPiabikQqKktkTDCWTyeTzUQeedl1DbvAKi1WZsa26d2X3Ce3nl0rlXG6XGRYQ32O+gb09JnHBmnkS0qUEEEL6YOa7S+zOgO9i+XoNrAETPxHoi9rEFHY30qWOxJqrPngFpLp1u4a9lQKaYKwG9gen7rNsF3sRFwka7SxE1u90OYl3nSmN7hu3KhAMnAIoQPqG.PFwXUzQ4LOyFn4JgPfDRHArfEr.rt0sNYGxjWqLKwnDYFsBATDiIDpsRma2tEqcsqs.Vkb.CX.B61sq5yRWie3G9AwC9fOnJqBx8pC+2QGczhku7kKb61svkKW90mo+u+u+uKajSzzl1TQZoklpwMeNnhFtc6V3wiGwS9jOoHt3hSHD9+w2e+2+s..hAMnAIu1Nb3nXec73wi7GBtb4Rr5UuZA.D2vMbCx01hqEfO7gOr..h0u90KDBgvoSmBgPHF+3GuvfAChd0qdIxImbjedGNbHF5PGpWsVMOxF05oiEtvEVhsN8kCz7wa7FuQArPtVurpMxIGv.FfvgCGxmW55ocdthB7wPt4lqXhSbhEpWKz97sxUtxBLe6wiG49ugMrgI.f3jm7jB2tcK+rzqo+lFGz9KsyMtb4R3vgCQ8pW8TsOxaOKz0veQGxpUqB.H5XG6ne45oEzbv3G+3EsrksTHD9uwNc8O24NmPQQQbu268JDBQwZ+GeMhOt3y2qYMqQnWudQSZRSDtc6V09ce456vgCwQO5QE.PrgMrAIeDOd7HFyXFiPQQQzm9zGQVYkkpwzHFwHtrd2RqGBo+dAKXApF+kVX2tc4qSHgD75XwaQTA85G6wdLgMa1TMdJIzuKq.MVneupUsJuFo3WNZEKXAKP3xkKuddWHDhAO3AKzoSm3Dm3Dx8M74TsiCgPHrYylp+mGOdDNc5TX2tcQ8qe8EJJJhMtwMJDBgJZNd6ZUZAcc0qWunMsoMpFSkVnklFg69tuawi+3Odo95KDEb8UHx6YJiLxPLjgLDQ+6e+EolZpxwhSmNkOyu8a+1hq65tNwO7C+fpqo2V+JMfd94iwQLhQHTTTDicriUHDpoO4ukAnz.Od7HV8pWsvfACBmNcJb4xk3oe5mVdN4AdfGPE8xhybmVY0VvBVf3ltoaRjTRIImCnqsa2tEe629sh5V25JVwJVQAjCQ60xpUqhpW8pKhIlX7JeYsulteZWCnO2AO3AkQavQNxQ7JuDBbZEk0P6bOce8E4aJNPK+SWtbIRJojDcnCcP7xu7KKra2tjNIgKcoKIF8nGsnG8nGhicrio565uPt4lqvnQih65ttqh80kuVqcd5.G3.x83+8e+2p9eARmOoyaUqZUSr3EuXgP3emeI3wiGwwN1wD2wcbGhoO8oKb3vgbeOMejYlYJF0nFkn6cu6hjSN4BLVH8IoqWQA9mm9Nb5B7+1sa2hm5odJA.DSZRSxq7d72xVWTPKsP9yBWdY5Lr1mGeAZmib4xkH0TSU..w8ce2mPHDpVq32OeADshMu4MqR9HZuF+ZweNBDzOTKcNZL8pu5qJ.fXSaZSBgP8bcoYuge27ebKFystnNc5v1111jeF2tciHhHBL8oOcYM0feM93O9iQrwFK13F2H.x2RpjWUIuESeuTRIE73O9iigO7gK6HMDDrBdGA50BMdARqkmoFB.2aHbOK7Vu0agHhHB46IJi8nTwEj0r4Qqj+1qHTcWoN0oNRO1nslPTXfuVo0qCT8CC.XfCbfp7Vsu.Z8frnd1Ymsbuz69tuKl27lGt0a8VwRW5RgEKVT4MsCbfCnZrviV.u04VTTTvfG7f8qVUmue5fG7fxH9flm0NWwedaaaaKVvBVPAJB09hmcJu.26r+2+6+UV6iz5scsiWhV.O5anWSdVzpUqX8qe835u9qGQGczphlN90lONn6q16md85wO7C+frNj05V25BbVW60nzVeMJL3Ooyv2WWVTzbA.N4IOIzqWur8NWb1CJXQzi2FWJJJXsqcs..xNlC+rPQAxKsbO0Sdcd4Ke4X9ye9n4Mu4HgDR.gFZnx8Qtc6FG5PGxqdw0amOoW+XO1iA.+WiMgeO98e+2U0n..xOJLoWS+1sa2n8su8XFyXFAbMZCBh+MBwHdJaXCa.CZPCRE+a97O+0z9rKbgKnpnUS6kb61MxJqrvm9oeJZPCZ.hN5nkyYdqwOv89HOxen0YCFLfe9m+YbpScJnWudz7l2b48j+asWqRK3q+70O+AMdNMM5ZS+1eb84d8kOmnSmNDZngh24cdGbS2zMgG5gdHbzidTYT+XylML9wOdrt0sNricrCY2HiFe965eBI2AEgComd5xF8SbwEm79VVRyujBJpJ3QY6LlwLPe5Sej0dwu5q9J47UwoomPfpqOiZTiByd1yFcqacC6XG6Pd+83wCV6ZWK5ae6KV0pVEF7fGrpnpk2DOrYyFzqWO1+92ORKszPO5QOT0DJDh7aREZeFATK6.+yQ5ebS2zMgq65tNuRal.W1fxZvoYQQgHP9zu3iC58JNQNDE85ZitF850ia3FtA7IexmfCe3Cim4YdFYT5C.jbxIiANvAB61siDSLQTu5UOUWS+E3q+EG91z2k9MMtonl6q+5uF..spUsBM7e6zuDef.I9bz4MtLQ9a5HzdfF0nFgMu4MijRJI7jO4ShrxJK..Yszr+8u+voSmHwDSD0st0s.zzH4m80nqSqr670IZOLEwmW7hWDu268d.HuFff2lKJOjslCNsP5LD8ZAq9qw0ava76ubftFbaGnnnfTRIEnSmNTyZVS.jejpQvWqQnb4i5YO6or9sGTPAgCcnCI4+5wiGUc92.gnmDHe4jz9rp8u4y0kl8Fk4T7IBmomd5HwDSD.4+vLoIMIznF0H4DOkNTaXCa.OvC7.xCe7znLnfBRFNojvs7TuYMqYM3Idhmn.LU.x+vL.JvgTJDsIFFTpkXvfArm8rGUFhfHZOkoLEz4N2YYnHSJcvUx+pAjYlYB.fpW8pqJEC7UAJ71FXJDZOxQNB17l2LTTTvsca2lJlF9BAAtPWjfzd73AKXAK.OwS7DnO8oOXCaXCnV0pVpFKJJJnJUoJpLvD2vRj.2jf7ssssE5zoC+4e9m..p5pYkFPOud73AQFYjpJHlz6SiIRP.c5zgt0stgMtwMhpTkpTgm1FEFHFK6ZW6Bu9q+5xwJY.Zh4KwXl9MMOzm9zmBjBETXee3CeXjUVYgN24NKEhjtVEGgcol1w27MeCDBAZQKZgrHfdMTzPmt75Butb4B0rl0rXc9E.xhBLf5hvK82m9zmFaYKaAtc6FcnCcPUgJ0WfV5zW3BW.JJJ3sdq2BCaXCC8oO8A6XG6.UspUU97PJlSz7n2mLtFsGiuWqMsoMPud83Dm3DpRMoRK3cUmF0nFIKU.bGAQz+3cwp68duWr10tVT25VWoA9o4BJT6qnAm2cxImLF8nGcAde540aF9SHDnW8pWREM47kzqWON4IOItzktjLMS46y708m7O6m8YeFzoSGZUqZEttq657SyBW4Bc+aItfSSfK6kQiFwDm3DwHG4HQm6bmwu9q+JN8oOMF5PGJRO8zwm7IehzPl7FFSwkFegAJ0g4NO4m9oeRdlJlXhQ9rvM7ZfD3EhZmNchHhHBL24NWoBWu7K+xxRBRwAzbM2H78pW8Bexm7IXzidz38du2CNb3.yZVyByadyC6YO6A20ccWRY8o0ch2LPdJqIDBbzidTnnnfa8VuUUmaId49x4SxwY1rYS5zjG7AeP.D3sFUVAtQIn8wz4LgPfHiLRr7kubT0pVU7HOxifSdxShe629Mzst0MzwN1QrvEtPTkpTEUEIb95UEI3OWz3gJSLKcoKElLYRlxo.d24DWICx3wD+QROhksrkgZW6ZiG3Ad.jbxIi8su8g64dtGz912d71u8aKk0gKCAm1l+d9ynQiXO6YOR5O24cdmATFsrrBb9dbYZRKszjqUZ+bZ+t958voSmXZSaZnCcnCvlMaX0qd0H8zSWdOIi+QNj3pQTlSUfN3PciHtxJ268du.HeENCN3fwYNyYvS9jOI.TGEWT8PgpeCbOFQs+a5+ul0rFrvEtP4XvaGrHBnbkY3eNJWoA.9se62Tcsb5zI5RW5Bl3DmHDBAdm24cvnG8nKfRXWsPzk5BF0t10t.dQtn.mP.Gzb4+6+8+j+exi6j288EvUFzsa2Xtyct39u+6GiabiCO+y+7HwDSD0u90WZ3ERgJOd7fF1vFhnhJJUBBP+e5225sdqXCaXC3oe5mFtc6FW3BW..Ezh9kTvYBEVXgI61pbCtQFJF.nV0pVX9ye9HwDSD0t10t.QLAGAJFzQQQASaZSC.4SKfXdS6C3JyROG8t28FUoJUA.46sIdTiEUTQgl27ligLjgT.idC36sNZylMibyMW7lu4aB.fG8QeT.b0iAu8GfhhTxvy9J3qQjRt7nFTmNc3q9puRddqEsnExumur9PFRkqTWBIj.hO93wnF0nvXFyXvG7Ae.hHhHPPAEjrtQQBLbC2vMfvCObUFakDpi1mdK2xsfssssI4Qb7iebep9u3KfLnLsuNnfBR1gi4F3g9LlMaFgDRHXgKbg3C+vODQGczpLd.2PRUzcfKBzdfUrhUH6jeDeCO+ac4hLRpNc5TIafEKVPSaZSUEUDbkmpacqKtoa5lvi7HORIh9HIWhhhBxImbvJW4JA.vfFzfpvqOiUF.smD.RCiw2yBj2dx92y+FtLK...B.IQTPT89i0st0g90u9gd1ydhF0nFg24cdGT0pVUUQXF2gl9KPFngSOYsqcsPHD3ttq6RJ+.PfCOUN3zG.frl2c8W+0iEu3EC.fe9m+Y7Mey2nxPS9J8I9mijKpUspU3K+xuDexm7Inssss3G+weDaaaaSUGhmhTSNedBJJJnAMnAn4Mu4nacqaRmfwOi5K0fO57YRIkD94e9mghhB5Tm5jzfbWoCs7.46SI9Ct+2tL4rm8rw8e+2uzPXu1q8ZXricrRYtxM2burQPeEETTTJfSPxM2bwe9m+IRJojfCGNPO5QOT84ub5bbkH3NBlmgDAETPXpScpXfCbfnyctynqcsq3Ue0WEiabiSEsNtQZzNm4Oj+ku9sgMrAoA0iHhHtpYMBnfmmnteaTQEE.J4QzH2PZTVPrvEtPXvfAjSN4f2+8eeU7Ko.LxiGOADNDs7FkKQHla2tQFYjgpV56PFxPPyZVyjE3WZx+8du2Comd5pLlDWYWEEEzxV1RzktzE.j+AdhgIQvd5Se5xMUZOfyYJv8pGoXDWofScpSg0rl0n5yT0pVU7FuwafPCMTrksrE7TO0SgCe3CKu2TAN9pAHD4UX5AxyfXdyieE022afVWn49G3Ad.Tm5TGUL+70P1EH+Pe8jm7jXyadyHgDR.yYNyQp.EOZqHAyUTTPCaXCwBW3BkFhhtdCX.C.e7G+wX26d2n28t23Vu0aE.PFgX9y0e9bjISlvm+4eNti63NTEp0wFarXIKYI3m9oeBO6y9rH7vCW022aBcFH.EEEbnCcHr6cua.ju2mH5DjQF3QDDsdMzgNTUyyjvuz2qF0nFXe6ae31u8aWFQoTz6.T7Bs1e7G+Qjd5oCEEEbK2xsHGqWCENn4Zp8YSQj.fuITMmV7kCqcsqENc5D228ceHxHiTFEg9x0mRGOhOhQiFQRIkD15V2JV4JWojNOwChafY61sCa1rgF1vFhDRHAT+5We.jeDs8HOxifO9i+X78e+2iN0oNgq65tN3xkKI8R+gAmH5Vzdd850iPBIDroMsIzpV0JUoJeW6ZWw7m+7we+2+Md1m8YkM8Fdziw8lXf.nmI61siDRHAYAUGHemoX0pUU6S3J+OlwLFDVXgcYSGzvCObbfCb.bm24cV.mY4K7u3Fn+G+weDm8rmUF0PAZzZCDAeci6gZdjQyipQ2tciPBIDoQPI4xTTTjEDXhOf+BD8AZsNkTRAIlXhvjIS3ge3GV0dKNOg.AGlPicxndZ4Wd+2+8igO7gCylMiEu3Eqxgl9B8IpYCv4oRFthNCERHgHiPSsoYjV4soqI.P6ZW6vO9i+Hpd0qtbMktF9ZATltGT5yonnfa5ltIDbvA6qSgUpAmFDQuTKsPty+AxKvDnFbA8+UTTfEKV.Pd7FzVxPpn.OCIn8MVrXAacqaEtb4BlLYBMsoM0uz.SpLBsoKKoCMc9ln6FTPAAmNcV.G4y2mnk1l+h+lACFvYO6YwZW6ZgNc5P26d2qvaFcUTfluSM0TgNc5PcqacU89E2rqfvN1wNvjm7jwHFwHvZW6ZkN1cBSXBXO6YOEHUTCjbHZ4IJy20QGhnH+fXX1ktzEYZm.jmxAm+7mGSe5SWkQp..dtm64vfFzfPt4lKZW6ZGhO93wK+xuLN6YOK9ge3GvK9huHRJoj.P9DEyN6rwZVyZv3G+3UkVE.pUpFP8lM50TGwXSaZSppcYd73Au669tnYMqY3.G3.x5YBwn1tc6xn04pAnnjW2KA.REBI3KGZ8lBuz9hie7iKEjghlvRhGpb4xEhHhH.PdQz0xV1xPqZUq..jcUPRgXR.OZLYznQ7TO0SgQLhQfKbgKfryNaUcvNhvMQ7H0TSU5wM+A3FsgpyOwGe7nm8rmHszRS1wPoPbl6AL52dyPXARdI6Tm5TpRORf7iBP5rG+roSmNQ+5W+jcIJR.ONCa56SyejQ2Ig.3Q4SgARfhst0sB.fpV0ph1111J+eABBEFHCZO1YNyYfa2tkoQlut+iSilSqf96icrigu669N.jOMBse2hBz9AKVr.mNchV1xVhUrhUfa7FuQU6g3uFHu8TgEVXvjISX3Ce33wdrGC4latHiLx.Uu5UGgEVXp1iDRHg.CFLfScpS408rkTPoKF0Yk83wC5Uu5E5V25Fb3vAN0oNEZTiZjzyejhpd73QEuJtPn7yKUjfFqomd53rm8rEPQNheuGOdT0Iqb61Mtwa7Fwy8bOmpnNBP85n14eZdQqRhWNPd31fACxZTTPAEDZSaZy0nM3Cf6nCx3zzdRZd0nQiXMqYMXFyXFXCaXCnAMnA3Ye1mEiXDi.KbgKDgFZnvnQip3E6Oc.DI+Hor3t28tkFiq0st0vkKWpRaYtblUz6A3QLGWt.94gW7EeQrpUsJ7Ue0Wge9m+Yb629s6Sc3O.0orFctRmNc3G+weDO1i8XXZSaZ39tu6Cu1q8Znm8rmXkqbk3Ftgan.ee.HMjMuToPxhQF1hLRpuRWh1esjkrDXznQ7POzCgvBKL.DXICTYE3AD.W1YZ+.MulSN4fYNyYh+5u9KrksrEjQFYfgMrggibjifILgIH0aiNOFn.N+.5YzkKWxxyy8e+2OpRUpRApSwz28J80eN8G9Y5KdwKhW+0ec7m+4ehu3K9BjSN4fAMnAgicrigQO5QCiFMBqVsJMjt2Lbs+JCXzqWO16d2qbMJ1Xi8ph0FB7mU5rIIurVcpITblal27lGF+3GupuqEKVj0KvQNxQhu8a+VTspUM40kbBYfb41or.kKbq83wCxHiLTYnqa3FtAoWsoTdL0TSUx71jISH5niF6d26FyctyEMu4MGsoMsA..G6XGCJJJnl0rl3du26E6bm6TlFSjPKNc5DabiaT50Pt..bhn7TkAH+MZjPUIjPBRioYznQ77O+yi669tOjd5oiANvAhKdwKJIBSdhKPw65kW3HG4HvfACpJ5lZUdzW.+vH.v29seqLBBaVyZlpPrGv27.KuVjL8oOc7EewWfV1xVpJ8HonJxaeW.H8bRMqYMQiZTijsrYZeiMa1jFjJ0TSUUMtxe.mNcJEfgGgbQFYjn10t1RBYNc5r.6s0Bh4VfDyFdceQa6DVa8PRud8nCcnCHgDRPV7yIl0bA8Hkf4dEmTTltN9pAac4xEdu268fhhB5ae6KrXwR.gxNUlvQNxQ..PSZRS.fuktKD3BiQquzeu8sucYJCGSLwH2G3qQmAW4vvBKL7Zu1qgO8S+Tz7l2bX1rYUmk0ZDEmNchrxJKUQzXUpRUPCZPCTYLL5LWXgEFb61sLRC8W6e30zCd3taznQXxjIzrl0LYSlPac3wsa2R5KDMA54MPQwGRAXsmW0xig306zoSTkpTErhUrBYzk3MOcyKt9TzPvSARekOtACFvktzkv69tuKLXv.dzG8QQvAGbINMGtZBTzTpc8kVGxM2bwbm6bwxW9xwV1xVPqZUqPjQFIVwJVAZRSZBhO93woN0o.P9EmY509CPmcoTLA.XSaZS.HuzsswMtwpZfPbdUAJQ4.Muv2ayqEiMnAMPVxB9vO7CAfuO1IZFbYk9jO4SvC8PODV1xVFdnG5gfISlv+8+9ewy8bOGhKt3v1291AP97347m4Nwid+byMW48xtc6EHkZKpm88t28hicriA2tcid26d6y0erqT.mGGWNJZuP5omNFwHFAxM2bwG7Ae.hJpnvMey2L13F2H9i+3OvHG4HQlYlor1rFHkJ3b43H461+92O98e+2ghhB5Se5C.7994.IYfKKAutfJDBjQFYfm3IdBjc1Yiksrkg5Uu5gl1zlhst0sh8t28hQMpQgye9yifCNX4YQsNq1eM2oWud3vgCrgMrA..b8W+0il1zlVrjgqxNn4VtgFI4kIChwCHfhCsqjRJIL9wOdYz0Re+byMWYcvNojRBibjiTkCIH8tuZCk4ZzQJuT6ZWaUBMDczQKKb1jP7m5TmRR3xgCGX0qd0nCcnC.HeiRz8t2cr4MuYUFFnl0rlXIKYIn6cu6p5pL6cu6Em+7mG.ELUIInk3N2fAewW7ExtHlGOdPO6YOwzl1zfSmNwvG9vw92+9Uc8HC6EnnHQ4ADBAN9wONhJpnfEKVJ1gz4ky.M4jSN3C+vODJJ4UKIZcqacALZUwQnWgHulf.0IL0oSmbuFE8h.4G0AjPvTzdQ0lJsonDPdQYFQ7f1C3Og1Z8..j0yAdp8XxjIUBrvUFWKCs.E3zoSz5V2ZU0nMsFdfqveW5RWvG8QeDrXwBLa1rp+GOxO4OqZMDF2CzEE73wC18t2MtzktDDBA5Tm5D.Br5RQUFvgO7gQ0qd0QngFJ.Jdmc0JXD8csZ0J9nO5ifd85kBwq8yTTfR0JZOwDm3DQ8qe8ga2tUQGf32vMthQiFkN.gR8Z5dSuGOpLHARyImb7aFLmLBDcurZ0pzA.7vd2fACxzDhaXcxXD7nFiFmAB6wIuZW0pVUoCwHGOAjO+adJPESLwfMu4Mi63NtCogwIZBZMLForXvAGrLReJN7wIk12291GxJqrfKWtP6ae684T5+pc3zoSYJ1Q7koFNw4O+4wS9jOI9q+5uPhIlHZTiZD.x2IlSXBS.ibjiDcricDe228cx0Oxfv9qHXftNgDRH3Lm4LX8qe8PHDXfCbfxZXoVitFHb1AH+yGZyZAZtVmNcvlMaXXCaXH7vCGevG7Axljju.EEEXylM40ZAKXAXVyZVXaaaa3Nuy6TkRX8qe8Cqe8qGO4S9j3cdm2QkrhdatingYwhE3vgCHDBX1rYUoIaQAWtbgu+6+d4e27l2boNBWMd9j6zTGNbf8u+8id1ydhXiMVL24NWYzy4vgCDUTQg2+8eeDRHgfG9geXogmIYeCDLXg1yY50qGe0W8U.HO56srksr.7ZCTF6kGfz+1kKWvgCG3O+y+Dcu6cGwFarXNyYNR5Wd73AUqZUCqcsqEgDRHX.CX.xl+CP9oMq1qs+.W3BW.qd0qFFLX.csqcEgDRH90qefNzJqpNc5ve8W+EhLxHUk4OkD7S+zOIkek6Depz.Qx37Ye1mg23MdiBHO4UanL2fXzBoEKVTY8wfBJHYZmAj2BPN4jibQn+8u+nicriRARHActsa61vEu3EQZoklpzBKnfBBu1q8Z.Pc9MmVZooJcGo6EAs05JhQoMa1v7m+7kDSZZSaJdy27Mgd85wK+xuL17l2r7d6xkKX1rYog8.Jn.RWoBGNbfzRKMDczQqx3h9pW33BlvKreW3BW.ae6aWJ3ImIlNc5TYDqBCTTIxK50.4KDF2R37zkg9LlMaV1NZo5rCOJKHFNm+7mG5zkWqhme+Jsft9z4.50VrXQ06we1zFwizOZURKPvSeFMZDUu5UGKZQKRtNPQHJM+Qmudq25svl1zlPMpQM7ZqgmSDWa8NfVSIB99pmN0oSmr6RpSmNzt10tqJNW6u.cFKkTRAMpQMRZ.AeMB8tbvsa23rm8r3q+5uFtc6FCYHCQVn5o0Ue85yOOS64HiPqMEi3FBmLrDO5VnmWdZaQQ4XN4jCzqWOBKrv7aQOh11BNUOeHCBvoAwU7jFezObZKzbQffBiDcT850KK.3.PpDNIafUqVghhBF8nGM1111FhM1XKPiPgG8I.4+LxiTXmNcJ4i6K6en09su8sK8vZm5TmBHnsVY.TpNRmyHi3lc1Yi6+9uezv+s97UqZUKUxVoSmNX1rYLvANPjXhIhANvAhcu6cK+e7n1pz.scaUJ0dzqWOZe6aO.xu7ZnMBJBD3SvcvG24Y7TMNnfBBQDQD3kdoWBW3BWPFcB95dXpKuOlwLFr28tWrksrEY8Al.85XiMVryctSrgMrALyYNSIudhVEM2wowB.U0oLtCvJJnSmN79u+6CgPfV0pVgl1zlFPYv+xCvmmnyZ1saGG+3GG8rm8DSXBS.Oyy7Lpb..I+kISlvBW3Bw8e+2OhKt3voO8ok0i0.A9CZidIWtbgUtxUBiFMhXhIFb8W+0W.Y+4ok6UKPQIul4Sm6bmwDlvDvy8bOmpH.hWh.l27lGdjG4QPbwEGN24NmTWKsNZ1er96wiGr28tW.jG8FpAH3qkzjJ6fa+C9emRJof5Tm5nh9MWtE.ei9Gk5k7LmvfACX1yd13Lm4LH0TSEIkTRXyadyH4jSVFEz9KGJUYCJB+.UA2tcim8YeVrjkrDIiKtmwIPdmfHPQodFGzgOxiVbg04sUaJkODBgpZW.kFUDyeCFL.ylMKutTMihTxwgCGRuCy8RNQ.ktVzXl1TQFif+7QQrhc61QvAGb.wlJtwAowLYwXeQnQJWhonjhXZxe9IqOyqQMda8+xAtvvbkT301KuEl79ZT93K2adJP5vgCY8ng5HRz3iHP4zoSYM8BPcmRM2byUV7mKsBkqUPVGNbTfyFU1AOszHCVAjufKdqv8RmAKqYZRFGmtuz9esqKkkPmNcvpUqHxHiDW5RWREsuRK4aRvChNJWYS+gR8jxsjAlb4xErXwRA5XXkTX2tc47Awzu7TPWx3SDMJNOCh9DmlAP96e4F8sjBsFCmFOTGEpxtP+bdCZM5HOUqo4RxIFJJ4Ur8KqiVaZsk2cboTKhFqkG2egP.GNbfHhHBX0pUogkJsmgIddT2XkhRYa1ro54rrDbmZQQzIMVJqSactRHTjmQ667lSYBz.+LBQWfV+3F9kmkCTTu5KzmHZcD8NsocC2.WjisnWWdXTE5bAP9FeUKM5JRPNOgzMh6DFpF9nMJi40+K+g7kDMSJcTo6U444Ka1rgZTiZfLxHC..e1f1Z2mwanYjgVCjga240kcc4xEra2tp5Wb4QI4f1SQmSny9bZ9kkPqN67FmQffSEn4Gt8An0EpIAxitbJaC3ApPgA99VxgtzyLoiYoAb8W3APDoSafdjRxMhGM2R1swsa2XZSaZ3EdgWPJqCu7kTRfewM050qGO0S8TnW8pWpV.IlsBQds40oLkoHCs8O5i9HDd3gqRn8Se5Sim3IdBHDBLtwMNbW20cIuNjRZ+3O9iX7ie73EdgW.wEWbRCw3zoS3vgC76+9uiW3EdAoBGKZQKB0qd0CgDRHvsa2RAooM0lMaFYmc1xnUyjIS3bm6bXbiabH4jSF50qGuvK7BncsqcvhEKH0TSEOxi7HpRiFEEEb629siYLiYHWLHipTQ6oXhQ27m+7QJojBdsW60jLUowtu.sdhgtF6ZW6ByblyDCaXCCCX.CPlBDVrXA4jSNEoPibASn8ATz6PJZpM+m4iG+Anw.2fgNc5DKYIKAm5TmByblyTxffDJm9Lzg0ctychYLiYfQNxQh3iOdovykVFJjgf0qWORLwDwW+0eMV1xVlJCYTYFjxrNb3P1f.xLyLkFgkhjThFhc61gSmNk0wqxitGDQ7kLx.MtKOhhFR32dzidfniNZrl0rFYTovSuzRJn88d7jWyB4nG8n30e8WWJndok9ka2twAO3Aw3F23P+6e+wPFxPjFpveHzGk5vzODJuDDNgDR.G9vGFKZQKRRCyjISvpUqREVnz.66+9uGyXFy.idziFwGe7RF7kFPoTHInvG9geH13F2HRLwD8K6OBT.WwZGNb.61sKib2pV0pJOmXxjIjat4BgPfvCObeNRhKMfL1fVZDz3trDzZrNc40ctpQMpAd629skJZTZMHHs+kaXWOd7foLkofvCObLwINQ+zSh2A8bDTPAIS6mPCMTox6k07+zd8oRmfVCeFnBxQte4W9k3Mdi2.adyaVxuUQQQkRWYlYlnu8sunUspU3Ue0W0mn8yqetZkiianBiFMJOKxkYurd9ijCiLPLsmo7h+suL9.fJCEBjmxv8oO8Aie7iGwFarxLB.H+Zmn+vgvjNONb3PZXFtCdJqUXlbr48bO2CpW8pmL04.fj9dgAsqkjdCzdtxCiJUZfACFfMa1fACFPe6aeQ+6e+Qe5SeJWLFIPdmeyImbTUWb454Udb9fLrDOp7ITQe9TqNojAvrYyFVwJVA9i+3Ova7FugzI.7LxgbbSQgbyMWTiZTCoyt2291GlxTlBFzfFjr1nWRwEtvEvXG6Xw4N24jN43UdkWAsrksLf+rAf5fdg3aPkSoDSLQzfFz.oLNz7NMmysyiuB+hAwb5zIt4a9lwMey2rWSKQZvMfAL.rm8rGnnnfnhJJzt10NUeF61sisrksfMtwMhHiLRz4N2Y4+iHPT25VWLoIMInWudb228cCf7K715zkWgFmDbygCGn6cu6nt0stphnBtwUzF0S50qGSbhSDImbxPQQAKdwKFibjiTtnbjibDIyef78PasqcsQ6ZW6jQuC4gkJ5BqJovz5W+5Qt4lK5ZW6ZwZih1OKYjS549a+1uE..8su8EcricT0l2R6ANdceprDWNAKV+5WON24NG5V25lW+7zynKWtvV25VgACFPm6bmQG5PG7ad+jOO9y+7Oist0sJ226Oliqng1yjD3Q+Eu1sQmmHiDVd87WQMWyi9mvBKL4dQesKfUTfu2+K+xuDm3Dmvuu+Ze6ae..3AevGT5jC+gv7.U7mA1912NNwINAhKt3.P97o71y1V1xVfhhB5QO5AZe6auWqKGkTPmW1+92Ob4xE5bm6bE9bi+.bAR0B54iWC236q7WmQ7UTQMeSOmtc6FgFZnxyu9Kvi3df7js6UdkWAMrgMzueuJrw.4Ec5uqnLnQE48t3BZuQJojBb61M5ZW6J.JH8W5LTyadyQZokFhM1X84tktVY94mC3724ue4cGr0ewuor.ZkmjfSmNQSZRSPO6YOU84z987WfGQZkW5tP2KEEEDVXgI2eVZncWX7LBz.WNWWtbgq+5udoduk26Y0FY6kW3xoCPfBtb70+xu7Kwe8W+E5d26tJdBkj4O99fe+2+c..zyd1S+B+UCFLfG3Ad.oQ5ZSaZC5bm6b.QWD2WAY3bZ79K+xu.mNchpW8pC.3Um+VR1C62hPL9ffSbm1bnSmN7fO3ChW5kdIb9yedrsssMU0iGJEJehm3Ivl1zlvpV0pvS+zOsT.HRIinhJJXvfAYAsWKgaxRtFLX.MoIMAQGczp7JDeLQQ5CoDid85wN24Nwa7FuA..d629swPG5PUMwx8dAYTHZbP0RAZrFHzManCxbOWqMZqJLvS2RdAImLJVRIkD..t0a8VUst6qoLhVCTxumdSnSseF+UJiweFIlwlMaVUp3v8HMWHOgPf0rl0.WtbgXhIF40veZHOxnuTaPF3JiBen1zJkREGdM+hhPL5rp1TWrrd7A.uJnP4QJqpSmNjYlYB850Ki3Ttg.KsO+bCXa1rYYpQqsaeVZt9G3.G..4UPi4B+6KQ3QQM+pUPEdDJ3Ke+RK3oWJu4WPdWWQQQFgKe1m8YPHDnwMtwpRMkRK3Oiz9DNutJyfu+jKKAf5ltAWoB9d3xiHHhFGZoMTdL+ymOnz3SqAIJsWeNOc554uRI2h57IecjjEjRsHe4rs+b9mKWRYw0ur.TDh4wiGXwhEoRPz5JslR7SaTiZD13F2H.7M9a7Hh.PcSHf9+ZkYi6D5xi4Othp7wQfvZGervUllpkgb5Z7yz9546hyYD5dQ5t3Kveb9OmbxQJuMOEb80qu10y.o02hB7wZPAEjzA675pn+35e4.mWAOsm8Um0UVL9JKxBnRCzNFIcRBN3fUEQc7yLEG4OzZLve629Mnnnf11115Wd9u268dw7l27vy+7OOTTTvoO8okz7CDleKJPx3n0dRz9Csm6KN6e0hRsAwzlxaZObw2TW0pVUL1wNV7JuxqfUu5UiwLlwH6xEDw3tzktfQNxQhEu3Ei268dO7zO8SKYxRgCWW6ZWw27MeCrZ0pL+cIkkO8oOsjnxfG7fUYTBNQVsB2nnnfye9yiwN1wBWtbg5Tm5fie7iiW8UeUYjjYxjIjd5oKIdYylMIyiCcnCgoN0oJu1sqcsCcqacqB2KE7MI7vKt35MasFRjVS+ke4WPUpRUjVpUaHlVbLJjVBfda7oc9rzN+5Mga0VelHhcdySF5zoC+xu7K3bm6bHt3hC0st0U9Y7Gq8zbBkhlDiRtfsUlAktEzYXtA+1yd1Cl7jmLpW8pGl6bmKpd0qtrF.xEdq7BZuWkWdzhRcQtPwdSI7RJn8qz7J2vikV3wSds04PCMTTqZUKuJ7s+Dk2QtgUqVAfZmBQ0cQJUjzqWO1291GNyYNC5d26NZbiareKJSzdc3MMhqDhPLBdKhyEBA16d2Kl7jmLpUspEVvBV.hJpnTYP2JJ5Cdy3IkEfyOhROLs7pKsfOWR6q3MxfxRP78LZzHNyYNCF0nFERO8zw7l27jcczBCkEiOusNGnCEEEYpYQvaNknJUoJviGOHqrxBUu5U2md9zZPdBtb4BqcsqEyYNyAwGe73Ue0WU0ms7nFfd4Fmd6uqn.2orzdcJxW4x+Qe1xBZ6jdL6d26FSdxSFMnAM.u1q8Znd0qdE42qzBpFpoU9lh6yn2FKAJqwWNnnnHMvh2x5gx5w+u9q+JF4HGIZRSZBl6bmqLap.7sn3zeoiykilZE85GeOHErKjSmHds7fxgno4qxevoAR6C9se62fISlPTQEke44Wud8XricrHmbxAu7K+xXoKco3gdnGRZPu.cvsWC8Zs0fX+UzrVpuJbEyHifPDv4D5o++PG5PwpV0pPRIkDVxRVBlzjljJCZonnf4N24BylMim4YdFnnnfQLhQnJ5ltoa5lvV1xVPpolJZXCanJux+Ue0WAOd7fniNZL7gOboRddKLF0RvMyLyD1rYCUu5UGFLX.Ke4KW98Ma1rrtvPBWPLqb5zIN5QOJdm24cjBKFd3gitzktTgeflaDLsFSwWsNLIzM2vjBg.m9zmFm3Dm.OvC7.x+OuFKURh.jByiA7+Vqg4Jofq7.smkS.jLBE82biEPLv2912NzqWO5W+5mp8S9KquyGeTc3f2E6pLCsBnSBBQ0olu669N..zpV0JLlwLlBDQnkWQH1kCkGQfhIS9z14X...f.PRDEDUljzcz5ED+UDRRzI30MF+wy14O+4wu+6+NhKt3j0AHx3xkGd.srd8gSGky6iW6DMZzH17l2LLZzHt+6+9KVd3unfV5MjmznntrxfG.KLb4LrC8b+BuvKfe3G9A3xkKDarwJipbscuoxJTQu+it+TDdxiLQ+w8mnOvihWx347FcRYEHOY61sarl0rF7oe5mBgPfwO9wiu9q+5q3WeKsfjekKGN2wG.4W1OTTTjxSmVZogZTiZ3SxuQeFNcMmNchKcoKgG8QeTDTPAgYMqYgN0oNgtzktnx4hWoy+tn.EsCjA945MQyOZinN+oLlbYtb4xEdoW5kvd26dw+6+8+PKaYKw3F23JUW+hBj7MZctquJeSf95aQAdjuR6An0ixCY6G0nFE9i+3Ov9129PKaYKwXG6XgNc5j0Ttq1OexMLKO6gHdf.Pk9gbm66Kqg75cEU+x+y+7OQG5PG7Ky+b88eoW5kPSaZSwfG7fQO5QOvm+4eNhHhHJUW+xZnUGA50bCmqMnVzVmmKNnTaPLdGOghfKR3I9FHBgGd3X4Ke43+7e9O3EdgW.21scappUXJJ4E1+Se5SG228ceXVyZVXcqacX0qd0n90u9PQQA27Mey.HuhveCZPCjaVSJojv67NuC..l7jmLpRUphJAC4gam2Piabiko.HYsVBzghicrigl1zlJsPIUb86YO6I17l2L.J+qcIEFHuIQ4fK2aSZYtVTf2YmTTTvN24NA.P6ae6UMew2.6KaH0FxiZeuBynY9ClF7PTlSPi7TGeuC+dZvfAjQFYf27MeS3zoSzwN1Qovl7qSoEDQAZsjST9JAnU3GR.9ZVyZBf7d9qRUphpm2bxIGDRHgTtHzf18yE1dS+MTTTP1YmM.fLJw3q8k16u1TxfnU3u7b+d1yd..PG6XGUYLS5d3um+3JmAT9XPDsQ8FuIcnSmNjVZogUtxUBmNch1291Kel8GzG3OmZqGlUVpODEFzt+P6YwnhJJowuhHhHjOykWQdhVC6Ru1eFgVEF3JzRcIVtRU9q8+dKB8.J68fOmOZXgElTVyFzfFTt2vH7lbMA5Njhbblc61kFhR67F2.LYlYlvsa2H7vC2mifKZOu1Hm2nQinN0oNHkTRA..0pV0R02q7pFUQiQu8dABqebm6RAU.m+tVEtITRT3SKzqWuTgYSlLgZVyZpJidJONemSN4HM7WIQ9FsqwEW8OpHAMVIGL3vgihsdYkFTu5UOYMdsZUqZR5.7taZYI7lAeBjhPL.nx4Z7yolMaV1QS4xLWbiPad1DrqcsKHDBzktzE+17OISjGOdvC+vOLhIlXv27MeCxLyLkYnWfL3mE3F6BPcozPq7Nkj4N+RDhoU.edzBwiF.pEyd629si28ceWLrgMLz6d2aroMsIDWbwox6llMaFwFar3y+7OGYjQFnpUspxqeiZTifa2twwN1wv+4+7efPHvEtvEvfG7fA.vHG4Hw.Fv.7pxqTTPPF.g6ATdTPwMtCuETSchE5YmTjjGk.ARc1K95.uyRp0feEF7lBxd73QVP8iM1XkWK59UblG7Vjhn0yiz64s+tz.thidqnlRD8nmOsgQ75V25v+7O+CF3.GHZXCanT3E+kwvHCTPoNhYylU48g.AFFkVP6S3FmwfAC3Mdi2.soMsAQFYjXfCbfpXPDRHgTtM93LAKuA4AQZeJOTr8WJ7xoEvqeG9CZY6ZW6BJJJ3tu661qL17GvaJ3TdtVoMp5H9Ajvkexm7I33G+33QdjGAMsoMU9c7WJDR6EH5BjgwBz3EUR.2KqdidGk5bUspUECbfCTJmQg43qxJvGmkWd1l7lOEsVjyD7WJUQJLySkGGNbffCNXox.k0f1CO3AOXoWze7G+wKWt2dCU134pSmNu1Yw4cnTRtlKcoKA.00b1hBbdizq0qWOhHhHvV25Vwm+4eNZYKaIZdyaN.T6vfxZ5SUFVq3NjUKeApjR.jutV7ZZq+.bdEKZQKB21scaHxHiDCX.Cvub8KLvi.TR9lKmBtEWTYXsWmNcpLHh2zmtrDye9yGspUsB0t10FO7C+v.n7sbKTVnSm+Dz9SB7f7g2gd0V1n7U46HahHD4UGx9we7Ggd85wcbG2geYtf6vNZ8rYMqYnIMoIALAsSQAsNXldONsB.0FNqjt+sTOiPGd.PA79u1tOC2PICX.C.0st0ECX.C.coKcAScpSESbhST5QKtUVov5iX.DczQC.fCe3CC.f+5u9K7bO2yg8su8ggMrggYO6YeY8vE2iKbAr4iQsedJ+tIB2DnmORID5fBm.WEM3FY..dcySgAsFCiG58acqaEgEVXn0st0phRhRRDl3sC+WNka8mDM8lmCo8d5zoCVsZUUwLmG55okVZXpScpvfACXDiXDpHbp8vZIE70Ip9Yo88uRCz4yniNZLgILgJ3QSdnhxqUbE74Jf5O8hHoLsVlH9xYXtP.786Dyqu9q+ZnWudYABsrn6UUQJDEmwL+YiFSW3BW.ScpSE5zoCO4S9jpT9ze.szH8FOsJynndFpacqKF+3Gup2qhntDocOX44bOc9k2LGJt2+KmR1bGFQWeJRFJudFo41fCNXLzgNzxk6YgMNJuAWwXeg1g2j6kajA5+woaQo9cpolJBMzPQ0pV07Y4WJLdiMu4MWZHLBWNYsKuP445HOX.7VFy3s5OD+0bcpnwdYgiPo6acpScJ2k4hGUsz7UwgGYYo9AEEJpHPDv6mG4AN.uydx00n7vfEQEUTXhSbhE38qnblTfF3NCmbNGuQ7Q1CfmZyEGmgw4WqWud7oe5mBc5zg67NuS+h9iWN4SprXLL.0osJe8fzElf+n1wVp20y6JH7PKzlMaph1.p1GQQUkACFvccW2E10t1E5Se5Cl5TmJtka4VPBIj.N6YOqpB1JsYi96niNZXvfA7Ye1mgW+0ecz5V2ZricrCLsoMMLu4MODd3g62R6GZyMsQ+PG5PplzMXv.zqWO9ge3GfMa1fKWtBXLFFP9EYYtmGz1QLJLvquATzJY0pU7G+wefTSMUbO2y8T.lzUlTDiFq7B5IP9BdvqyFzbIQXb1yd13hW7hHt3hCwFarxqC+ZTZgc61KvYHgPHe80vU1farofBJH.j29.uwPnjBZ+rd85kEId58KJPmSH5y7Tt7.G3.HojRB8nG8P5Abs0LtJ6fqnBuitAj25zrl0rP5omNhKt3vcdm2I.xmWB24JkTPoHGQqwgCGEHRPtFtxE78e75ZI.748WTwAlTPinqvMvFWdJGNbHih+qgxVvqgo.pSARBjgDHPQHHs9SQBl1NeMsNSq61rYCe228cHt3hSRu9Zv+.5bDwi.H+0ExQW7TtmuFRJeyq6jWor1vijWRNGpNL6ujuorDb8dnwNAskRAtCM4YEAuQInUOjqgJVPqM7RED43HRtW57HEvLZkysnt9jN1+xu7K3Dm3D3dtm6AAGbvWwbFuz.d.fvMZNkIhZK0NDJoyc9kTljagSJ5nBJnfT0U7HgmHA2nGfl0rlgMtwMhctych4Lm4fQMpQgm8YeVYsEqYMqYnpUspvhEKviGO3bm6b3Tm5TvkKW3PG5PXxSdxXnCcnXbiabxTViThzeDxoDyost0shu3K9Br90udoxejg5.xqaiEarwhm64dNLrgMr.l5HF4wAJR2JIiKxChbO0tjkrDnWudL3AOXUDv0VOnBzOTqk.F2yTTphP6u4oU1JW4JwBVvBfa2twjlzjjOy75yl+HJdLa1Lra2NLZzHrXwhjALkJmWITmftFJbPdQNmbxA.PRK0eSewsa2HzPCE.9dDnQzS30hD58W25VmrYS.nl4eYQ8Cqh.j.8ZqCg..aaaaCyadyC..u3K9hpDLlWnqKMPaTGDbvAKoEckv760PgCRPQ2tcKk4hLDRwc+EOBw3M5HRAU.HihAOd7.a1r4eeXtFJ.zVNGnea1rYYT5wkcg9NDMVGNbHkcSud8vlMaHnfBRxSgqLwe8W+ERO8zQaZSaJVkTiqgKO3QGF2vVbmRyijHtLj7FPCf5x2B++WYF1saWVPw4MsqJiOaZqaTz5FcVhbDFY3LsMGJhtJOh4tlwwpXA2Xs.ELZcICWwsABAeIkv4F4dyadyPmNcH93iWd8uZW+NdzyRzFoyGjig3+Me9tjHCrhvOTrKb3vAlwLlA13F2HLXvfzK0T81hrbJs.SoUHQDvkKWp5feYlYlH2byEtc6VVvE0NAQOnQFYjnF0nFxMkDSG+UNPSJ7jbxIiKbgKHeOZLQoO..jd0nwMtwH7vCWFkbUjfhXsSe5SibxIGzrl0Ljat4hfCNX3vgCeRnYiFMJK1izm+.G3.vnQiH5niFgDRHxMijQPqrzAD4dixkKWH7vCGW7hWDlLYBokVZvpUqn10t1RgHoClm3Dm.Nb3.QDQDH5niVU8miDBRaJ1VRAMmmbxIiryNaT+5WeXwhE+x09ZHvFTJ18m+4ehZUqZgZUqZA61sKMnRoMJqn8+Nb3.W7hWDm6bmCMsoMs.EuReAjwYLZzHra2NN5QOJb4xEtwa7Fgc61QXgElToMpYjTYfFQggzRKMjQFYfXhIFX2tcY8UJmbxAG4HGA..UoJUQV72UTTPvAGrjdSok+fCGNPUqZUQFYjALa1L9m+4efMa1PCZPC.Pkqn08Zn3Cdsp6.G3.nZUqZnd0qdxn9pnN+Rm+o8tTw0lT9lL9NwijjY3vG9vHhHhn.EJ8qA+Kn0PWtbgfBJH3vgC45RVYkErXwhz4Btc6FVrXAYjQFHjPBA1rYSRqIiLx.ojRJ3VtkaQpfNYXd56lQFYfzSOcT6ZWaDZng52je4pYPFPVHDH6ryFFLXPtNFbvAK0WwpUqHjPBAtb4RFTAG7fGTxyO6ryV1EYuRfuIASlLA61siCcnCgZVyZhpW8pqxAzkW0ovRJ3oxEoyCoKfCGNfCGNPjQFIrYylTWQd.CPcyQ.fCdvChZW6ZiZTiZHiJ2qTVmqrBR9agPffCNXHDBjUVYgvBKLb5SeZjc1YiF23FKclbt4laApY2EF3MfvCcnCAOd7fa4VtEow1Bz2+WVCxdRTPqPNA5bm6bH8zSGyadyCiZTip.YdRIsoRUpMAOYTqN1wNJMb.23HbijPCV5PtISlP1YmsL5uzlqtjf.TZ7PE3tPCMTrm8rGr7kubLfAL.oBWzFRxR79h.gEEHOx3vgCjYlYhPBID4lTtQ3HKWFbvAC.nJ7IqHAY7wDSLQjVZogm7IeRoBaVrXwmOvQLgMXv.99u+6wAO3AQrwFKd7G+wUEpuZK16kWEW3RJLXvfz.gzbAE8UqbkqDojRJ3EewWT1kl96+9uwa+1usLx.G+3GOhJpnjL5HE8ohtXoUgTNC2u7K+RrqcsKLkoLEYWV7pcBlWoCxSh6bm6DgGd3nksrk.v6cToRJHlGqZUqBUqZUCie7iG1rYClMa1mL3lCGNfEKVjmAHAbo1G8i8XOFBJnffUqVkOOVsZEgFZnU52+t5UuZjRJofgLjgHUz4e9m+AIjPB.HOApF6XGKpacqqbsh64vR65GuNcZxjIr8sucrgMrAL9wO9BjJUWCW4AZ80iGO3m9oeBgDRHHlXhwmKYDjWWSN4jwBVvBvC9fOHt8a+1AP9Q2h15i2N1wNPFYjA5V25FZW6ZWY2C20fjVAIKqMa1PDQDAxJqrvl27lwYNyYv.Fv.jcgYdGxlTd3Tm5TXoKconssssXnCcnxHxg7vNYnkksrkgTSMU73O9iiF0nFAqVsVoKJcBzfhhBrYyFrXwB9se62vG9geHhO93QKZQKPVYkkL.AH4mc5zIb4xERLwDQ8pW8vfFzfP8pW8j5WPJmSQAXkcP5WsqcsKTspUMzhVzB4yo+pr2TVBdPBP5vBj2ykc61wrm8rwMcS2D5XG6nrAnDRHgHirEylMibxIG7Mey2fyblyf9129hl0rloJRBuFpXAIyJUNAHCYut0sNbjibDL4IOY3zoSU1JPmNc9TZ+RNj5m9oeBG3.G.wFarXnCcnpBthqlAw+SaMb6jm7j3vG9vnksrkEHRK4oWdwFhRIb3vQAdO2tcKb4xk7u4ul9+tc6V3wiG46Y2tcgSmNEd73Q99tc6V02woSmx2669tuSXxjIwBW3BKvm0a+coEzyfc61KvXm9+7mSs++JJPiim4YdFQ26d2kiQWtb4yiQ9ykSmNEMsoMUXvfAwwO9wkyy70s.kmceA7wsPj2yJs+8odpmRzoN0IgGOdDYmc1hEsnEIzqWuPmNcB.Hl8rmsbtgtF1saWdc7Gfecl9zmtnd0qdBgv+u+9ZnxEb61sWo8VR.s2cTiZTh1111Vf2unFGdCspUsR..we7G+g7LkGOdDNc5rH+tUlvjlzjDwFarBOd7HrYylXoKcoB.HTTTD50qWLyYNSU7F3zc4yEkFvmGW3BWnnF0nF9kq60PkCv2ao887E3wiGga2tEG6XGSz111VwTm5TEYkUVBgP8dqKdwKJdlm4YD268duhibji3mF8WCEFHZENb3n.zKb4xkXwKdwhXhIFwu8a+l78I9Bd73Q7se62JZTiZj38e+2u.zyo+1kKWB61sKVwJVgXnCcn9M9JWC4CZs6K+xuTzfFz.wm9oep7+QxL51saQxImr3tu66VL4IOY4YP92WHJdxtWYCb8IpL.93Tqtq1saW72+8eK5PG5fXlyblBqVsJDB0qkokVZhwMtwI5d26t3Tm5TxuK+2WCUbfuejudX2tcwDlvDDst0sVd9kuWvlMa978voSmhV0pVILXvf3W9keo.2qq1Ae9fjkVK8AGNbH73wSArUQwEk57ovnQiRq4y8VM2xlTp8b4JPjjkUo7kWv51e7qI8ic61QiZTifCGNPRIkj7+SVvk71hvO4cbp3WRdVid9Hu0HD42tyIHBP7NOEoc1saGVsZUNuvCa2BCZ8RSBIj.NxQNB5W+5mpndfaMadJkFnCgH+1RK8rR+lpeJaYKaA8nG8.iYLiQ9d21scaX3Ce3xHhgG0i.9u4.9YGJb6826uuFBrA2KgNb3Pt96upu.zd+vCOb.jew31W7NEuMHSWqUtxUhe629Mz6d2abK2xsHOKQd2gRQ9qDRmO61sCmNchssssgN24NigMrgAf7lGZQKZAd1m8YUUaIHuc4uh9BuMORkpfqgqN.UqnHd5TJUWbht.c5zgF0nFgu3K9BbxSdRLlwLFbgKbAIOlyblyfANvABEEErpUsJbcW20cM9OkCfnASQRDf55S5y7LOCdsW60Pu5UuvV25VkeVOd7fDSLQ7vO7Ciku7kiAMnAI8xNP9YpAO8zFzfFDV5RWZ.irqWIAR2lt0stgMsoMgoLkofEsnEISeUgPf8su8g10t1gG3Ad.L0oNUDZngJadB754GwC4Jg0Hw+lEPzyBWehJCOe7H0WaiEyjIS35ttqCaZSaB+9u+63odpmBW3BWPdN9jm7jXHCYHvlMa3i9nOB0st0UUYT3JA4iprCdzJxsGgISlj0PSJ0doZemGOdfYyl84zMOwDSD+5u9qn6cu6nksrkxztsxv9+xZvs+C.j1YgG85TM1ll2H6xTRju1uzkIMZznpNMIs4Qv5zMjRAjfZbBH7B1p1PEkByT.nxnC0pV0BUoJUA+vO7CpxCWCFLHYd3uB2PtfHDn5gF2fJ74DsFErhF50qGlMaVNm6q0XMZsvsa233G+3XZSaZH3fCFyblyDFMZzqs2a.euCWUQCZ+Eua2P0.sidzih+u+u+ODe7wiu+6+dUOSIjPBHhHhPVi1HEP4Bq5uV+4JTekT8i3Zv2.uC9xcbf+Jb5IZ0VsZUkg88UClSWCOd7fye9yiwLlw.CFLf27MeSIsdd8ngJB+WIfibjifCbfCfdzidfu+6+dUzQV3BWnjlKu8bC.YX2WZAesh34Q7jtFt5.tb4R5HFgPHqIM9p.0bkQiLxHwRVxRPspUsv.Fv.PpolJ90e8WQW6ZWQG6XGwq+5uNBKrv.PfeCy4JEvoSvSkFZ9uqcsqXKaYKXbiab38du2C4lat3ke4WFKdwKF6d26FcpScRlNk75qKImB2AHjyoqrXPhJCPqtAwDSLXG6XG3a9luASXBS.1rYCqcsqE8oO8Au669tXDiXDRijQks.xI17NH6UBm+z1Q93Fxux3yGIWF2vdVrXAqZUqB0qd0COxi7H3Tm5TX+6e+nW8pWn8su8X9ye9vhEKxFjg+T1tqgRGH6SPk9IZOIYzKZcljuk59g9RA0GHu5B8Dm3DgISlvhVzh.P9cVwJi6+82f2YO0ZaAxvXj7NTfoPAoUIA9ktLIGZa8kdaSg1uC4kDtv7ZivL55QLyMXv.hKt3vG+werrH2QL44Jk3OA+YgLhGce3iW+88la7Ju08U71gO92glWIKZy6DjZEJhtNZulYmc1nO8oOHyLyDIjPBH5ni1qFUilGpLU6I3JSvixsyctygbyMW.nNpY9u+2+KZSaZC.fJinQ+e+MHgeHO+xY1dMb0Azd9lSiTac6SK36q81mkmm9jvXbEmtbWW94eJReiO93QFYjAl+7mOhJpnTQe3JEifwQpolppBlKMuMsoMMDarwpxyV.vq70Jsf2Mx3Et2qgqN.mWaw87l1yvDeloO8oiUrhUfN24NiKbgKf24cdGDe7w60tw70PYK3c5ZuslpWudDSLwfu9q+ZL3AOXLu4MOzpV0J7oe5mhZUqZIoOPN5SqLqzZp1Zux0junzCsmunWGYjQhDSLQLkoLEzl1zF3wiG7EewWfVzhVT.CBoMJj.txiWZYo9SkGPqdq7mGx3ISaZSCKaYKC8nG8.Imbx3C9fO.wGe7EHZdo8I9ZPKbMTxA+Lk2zsluFP0DLxnWz2m.u6dqkmJ8atwsyN6rQe6aeQpolJl+7mOpe8qO.xuPxest7a9vWjalLLF.JwycADm1nAOsIhmJh7zXSaW4poMso..H8zSu.WyqjXlyObo8.Ko7C2hnG8nGEG7fGT0mQmNcxh9LW4LJDCoHCgGQeDNxQNBtq65tv92+9Q+5W+vvG9vupgXMmfmACFPqacqwK7BuP4Zw.m7LXt4lqTvV.bMuHcMnJZDov0lCdprmc1YK+rG9vGF6e+6WlN3.4GoX7Htktt.4aTX5u4g08YNyYPKZQKvu9q+J5e+6Od5m9ouhSncuAxX0TTl5wiGDarwhQO5QWtPej3IP+1pUqEH8StFtFtbf6LLtCEICiw6rc7VfNUz1uFJaAQ+9+u8tWiNppN6Cf++LyYRlo.gXLj.g.FTBxMKPj6TIRMTzJJdgKBDkpVPAVAwZcwBnh1EdopnnBZaVq10pBpfWZsEofbq.pw.EMt.rBZTHPTBrHk.kPxLy4Lm86GR26rOSBV5KPtLy+eegjYNyjgYly9rOOmm8ySzW3A.2WjNYwa2mOeH0TSU0bozy.G8.12XGqftvSt+krAKo2Pt.p6Xrssss0U1lD8xcO5LBjZ4nwJuDxeVde5KeN47ERN4jcctXxRAjdPSiGN+pla5Kib4JKSedSxOCJszRwt28tUWTgHQh3pT4n+8.8kit9eG4ymiiC91u8aQt4lK93O9iwTlxTvbm6bUATSlktb74ldM66wIOf+27MeCl27lGFzfFD5bm6Lt+6+9QIkThqAGjS7Gntu.18t2c.T2Ii0XQlOVfdMevwwQErPIYvDkKeyu7K+RbK2xsf7xKOrqcsKUaeVtitb6kctS8GKP86XGNbXr90udLqYMKzu90O7Ye1mgBJn.7pu5qF2DLL.35jCrsswxW9xgOe9Nm5fHWHEIRDDHP.DJTHUckfH8k1hrSTATeGORlUABg.ssssENNNpwHF8nGM1yd1i54RVWqjWgJI4DEjALW1AUcbbv1111vbm6bwUdkWI95u9qwjlzjve3O7GbMVcrL8kdjLXBKcoKEIkTRMYSnQtLpkc.Ie974ZxZDctPFTWKKKLu4MO75u9qiMsoMgMtwMhkrjkfm64dNW0TTdEraZnWmF0qCi5Wvgcricfa9luYL9wOdr0stULzgNTjWd4gu3K9B0bGkedEcVNPWbIm+nLPxQhDAUUUUXxSdxvvv.adyaFO9i+33Ftga.+0+5eUcxv7Dha4qw1GRuFXK+4ZqsVrfEr.r90udr0stUrt0sNrjkrD7LOyynlKkr7JvO6aZoOOU446ouL028t2Ml3DmHFyXFCJpnhTyoV+hEI+rV+BJqub+jcczssssg4O+4iryNa7Ye1mg69tua75u9q6JAKjmaIGetYfnEfRJoDQ5omtp6bA.A.Dd85UrksrEWcbAYm7vwwQ7tu66JLLLDu8a+1td9ZszgRNWzXcKJ8NegrCWHDBQokVpn28t2p26RM0TE6bm6THDBwzm9zEibjirAOOxmqfACJF7fGrXHCYHhrxJKQRIkjvqWupmmUrhUHDh56bCWn5hhsjDcGi4G9C+gpNJI.DKdwKtY40j785m8YeVQlYlYCdsRwuz2OTuqToOlYs0Vq56K6e+6Wzyd1SgOe9D.PjVZoIJt3hE111hG9geXwPFxPD0TSMt9aH65XQhDQjat4Jt5q9pEcu6cWjPBIH750q..hjSNYQgEVn5wnOtTrrAMnAILMMUiQ7q+0+5l7WCxNriPHDKe4KWjZpoJDhXywnoK7zOFSEUTg3Nti6PL6YOaQM0Ti56UG8nGUba21sItq65tDUUUUrSD1DRe+6F68827MeSQO5QODacqasA2dZoklXyadyBgv8bmEBN9PSM4bm28t2snW8pWtNdYjHQD6YO6Qz291WwRVxRT2F05f94jEc2s6XG6XhILgIHl8rmsn5pqVsO7QNxQD29se6hoO8oK9W+q+kPHb2MY4b7u3S+83n6jq111hRKsTQ+6e+E.P3wiGQpolpnjRJQHDBwu3W7KDCZPCpAeNo2Q0G4HGoXfCbfhryNaQ6ZW6TySL0TSU7a+s+1FzQIkiIyO6adzrGBxvgCi4Lm4ficriA.2coOt4pT5..PFWlDQAQkHQhfEsnEoxHFw+4JvKKtcIkTRPHD3zm9ztdNiktx3xqBnLiMDZoTKP8ues28tWba21sgu3K9B.TWFEcxSdRLtwMNTTQEgDSLQUMnRuaeJ9OQG2qWunjRJA6bm6D0TSMnKcoK39tu6CqcsqEkUVY3Nuy6TckE0qiUw5j++racqapkJJPSWSCPVG9hDIBps1ZQaZSaTclGhzq4essss00UmTxue+vvv.kVZoX7ie7nzRKUsDNprxJwMey2rJSw73wCBDHfZbB.2MJhO9i+X7oe5mhyblyfd0qdgYO6Yi0rl0fJqrRbO2y83Ja0hGXYYoJp4cricDyYNyoIcojIWt6xi4ENbX068wCKYU57mLqu1+92OF6XGKxImbvK7Bu.BDH.LLLPnPgP5omNdi23MPFYjAl3DmHpnhJZteYG2Pe+a8lSxYNyYvy8bOGdxm7IwZVyZv0dsWK.peI2MwINQrt0sNL8oOcTXgE5p9QJ9O0wF5hKg1JuviGOXCaXCXzidz3EdgW.yXFyPctM111n28t2XKaYKXqacqnfBJ.0TSMwUy0NVf9RTNXvfnzRKE+ze5OECaXCCKcoKEsoMsQU6E6Tm5DVwJVAtzK8RwjlzjP4kWtpax1TVRVhmE8Jr.n9xPzgO7gwMcS2D16d2K.p63jUVYk3m7S9InjRJQkU95Kgc85tqWudwt10tvt28twoN0oPm6bmwLm4LwZW6ZwgNzgv8e+2uq5JF.b0Q1oldM6m0xG8QeDJpnhPf.APs0VqZodH+BRQEUD9lu4aPe5SeTSbSFDH8ZOl7j1hUWq85czxfACpRu1vgCiu669Nbm24ch8t28pVixxNsvQO5Qwsdq2JF4HGoZBP5clL80vd0UWspvzc15RFxSL9+kNoQqYx+uNu4MOjd5o651aJHeOV989ybly3JE7i0e+m99o2wv.pOPs5AjRtLImvDl.1291GbbbTK6F.fie7iiwLlwfQMpQoZ7Fx5Kndw8TV2DjKaxfACB+98650i76jwKiO30qW32ueDLXP7nO5ihjSNY08IGu9hI4jqkKyB+982poC+RM+jeG8K9hu.iYLiAO+y+7XBSXBttOYMvIwDSDKZQKBqbkqD4kWdX8qe83JthqnY9+AwOzqsiBg.Oxi7H3PG5P3u+2+6HkTRQsc5M9fq9puZroMsILyYNSbpScJ7POzCodNnK9zqUXqd0qFO9i+3Xyadyn28t2.n95KkbYRlVZogW60dM7q9U+JLkoLE7tu66554Ano43Jz4F4mE5kLAf512qhJp.iXDi.EVXg3Vu0a0UGnT9X9A+fe.drG6wvJW4JwHFwHv1291Q25V2TIu.+b9hO44QKKKSIjPB3fG7f3lu4aFe4W9kt1V4EQdzidzXricrvmOep4GKClkdbHjIqSzMqN8K3rdMgVuwlDOL+4VZZ1CH1wO9wgPHPs0VqJHNh+SABVtNZqpppZvfN1113Tm5T..nssssMXfiXk5bkLqrj++q28t23fG7fHb3vpS.JXvftpcP.PkEcBg.G+3GGqYMqQUr9hDIB762OBEJDJojRPe6aec09RkCDKeuVHDHTnPvue+tNw63gcViDIB5XG6Hl5Tmp520KxvWrOfkrgIHG3syctyeucbJJ9hdWNMTnPvzzDW20ccX6ae6psIwDSTUmJzK.+50llpppJ7m9S+IWAf0xxBlllXKaYKXjibjpwBjG.2ue+tp0EwqemLXvfnqcsq3m8y9Yp2iap5RPxOijSjxmOeHiLx..wNGCjt3xiGOncsqc3sdq2BCYHCA.vU1fGc1lNiYLCzu90Oz9129lqWxwUjmjk9pEvmOe3Ftga.Ce3CGsoMsQ84jrdB50qWUvxu7K+xwpW8pwm+4etZ7a8tcFcwi99OcqacCqe8qGcoKcQc+x8yz6TvImbx3kdoWBabiaD.wVq3kXMQWL80+71qWuXCaXCX.CX.pyoRerzPgBoZbI2y8bOX.CX.ttfV7X2W7433fQLhQf8rm8fvgC6pF4pSuCr62uebhSbBrhUrB08I6t2986GqacqSk.JQGHL8ycTO.pQ2Y24XyMOZ1CH1YNyYTQCWt7OhDIhqCTzgNzA01qmt2G6XGSMYN48IEqbx.xS9TtySe6aeQaZSaT2lkkEJqrxTQhVtSk9AYkYBxUcUWkpa1HG71ue+p+FmshspgggJXXxSlVuPaGKyvv.Se5SWkILx++1T1x4k+sdvG7Awbm6bA.uJgTczO.ZhIlHbbbPW5RWv.Fv.TW04yblyfxKubb5SeZWE7S816tssMRLwDwUdkWI762u564111HkTR4r93zyHW8wChGFa.n9t73Tm5TUiQHuh+MEzyTTGGGbu268he9O+m2nsPbhhlrzTzktzE0IpKu.L.v0I2I2dKKKL3AOXd7mlHxOKjK6N4XKidzitAiKquOudYwHkTRAWy0bMpOu0yDA5hG4mKFFFMHXyx8yzuHS5t9q+5OqEWcN+uVNzWJx5Ks4tzktft10t1f4BI2dYBHH+rb.CX.t1mjYHVSi9zm9fZqsV0JaHTnPn7xKGm7jmz0pUCntiCJWEEssssEcqacSEDSYyLpcsqctB7k9wSkiQqmrI5zSDkXkXXzZRy9Lluhq3JT0yFYvvjhDIBxKu7PVYkkZPE85YyF1vFfiiC5e+6uZ6kaSr1IiIGX7sdq2pA2WYkUFlxTlBJt3hcckBAp6fuYlYl38e+2G8oO8468Dk9ucBTxWCwR0nln6HH.0+8nq9puZr28tWLhQLBWKyT410bLoDd.RRm7.l5s04UtxU5ZarssQEUTAti63NPwEWrqNgi7wlYlYh0st0g9129p9d8YKnVmsCRquswBiMbtX.CX.Xu6cuXTiZTpaqodxL5K6cIFLL5bkbR4x86i9hco2gYkYdDqqQM8NaiEKGOWd+QmA4Q2A1iWFatkhnONfbrY8xfh9E0P9umsfVx4.1xh9mQQmwX.0uT5jmiQzk0mn2+rwddoKN73wCJrvBcc7M.fCe3CioLkofhJpH.T+E9Tdd0W4Udk3Mey2D8qe86r9bquet9eOfF9Ydi8c.94eSul82wyN6rA.bslZko7c25V2vS8TOkJiEj0WLOd7f8u+8i25sdKLsoMsFj59xnrFOLoMGGGjUVYg24cdGLjgLDU11AT2.w8nG8.aZSaRUC1zquDDbMoC42YjQ2OkTRAYlYlX3Ce3pIsHylFd0anVCjEY8N24NiW60dM0UoFntCT60qWzktzEroMsIz291W0xyQlstz2u10t1gNzgNfgMrg4p9QnmIcD0Rl9w+zWxcQ2.YjGWTdrOVq5ZYv111UvKkKY63g4+1Zm947HOlq9RphZcSVn0kimpO1ImePKKxjoITnPnqcsq3Ue0WECaXCC.0uJJ750K5ZW6Jdm24c9dCFF05jgnYdT2vgCiErfEfW7EeQ0ABt7K+xQAET.F23Fmq0au7f8VVVXxSdxXMqYMn3hKFCbfCD.t6FZ5owZrL8zq7PG5P3Nti6.+i+w+.lllHiLx.u669tn+8u+p2KhWVJS+uPVOkzOw.CCCTTQEgvgCiQMpQ4J6CkGfiosN0ZfLkscbbvwN1wv3F23vt10t..PFYjA1vF1.5d26NRHgDbcUo32u+uaG6XGv11F+nezORcaxINw26nV5ZrLFWO.XRxZhGWFGsboms6.b76VSzWBU5m7M+7K1PnPgTkzBN9YKKx80Zryqq7xKGSbhSD6XG6...YkUV38du2y0JofhczrGPL.fidzihe7O9Gi8su8o98TSMUWe4TNwsCdvChEtvEh27MeS77O+yi4Lm4.f5S0P8ZpS73WVkAEq7xKGqYMqA4jSNp6SF3GYVPEu9dTzhdYPpu7C.ZXgOswdbD0Rjd.ysssQBIj.Ju7xwDlvDPEUTAV6ZWqpoZD8igAO+bm96UrFPPsVnWLeAfq4NUd4ki4Mu4gCbfCfW4UdE0bIh9j2olWxr1yiGOXQKZQ3u829aXdyadXBSXBMo05T5+cezG8Q3QdjGAolZp3kdoWBcpScRkIQbYuGaPe9.6ZW6BO3C9fHiLx.uxq7JHszRq49kWbO4Xm5WHHYMMOwDSDG3.G.SYJSAG4HGAqcsqE8t28lAqNFUyd.wjC9e7iebbu268h2+8eeLzgNTL6YOajQFYf.ABfvgCiJqrRrsssM76+8+dXZZhBKrPL4IOY.zvqxY7zIxI2YNTnPvqWuvzzDe0W8UvxxB8oO8QsMxNfgNV3kcGTK8SJP9dirSvDcGNMd+8Mp0inCLissMN7gOLpolZPe6aec88dYlfvuietq1ZqEABDvUcBAH1owtPw1jkBfnuPhOyy7LXdyadviGOXjibjXKaYKp6y11lAaoE.8wX13F2HFyXFC.par7JpnBdB2svkat4hO3C9...rjkrDL24NW0wPZp5TwzEO5A2LRjHH2byEEWbw..3oe5mF+xe4urY9UX7M8UUFvYO69KqrxP0UWsqLCiyQN1Sy9mlxkpVm5Tmve4u7Wv5V25ve9O+mQ94mup.tFNbX333fbxIGr3EuXLoIMIjQFYnB7Uzeor4rnm2TStirrikDIRDzidzC.T2xbviGOvzzT0ILjc.Mtybczu5M.tWtA.0+9p9jSh0+NEE6Qul2YZZhK6xtrFsXrKGafN2333f.AB.f565l.7hMPsdHCFlj7JkK6l0x5lh91vLLukAYP3iDIBRJojTEG5zSO8FbAPoVdRO8zAPc6ClRJonlKpWudYvvhAHuHC.0MlYlYloZtXIkTRMyu5H.2I.goooqZ8lbr0rxJKHDBDLXP03pb9cwdZ1yPLf5uJW111PHDvmOenlZpAUTQE3Tm5TnssssnicriMX.D4K8yV.vhWlvVzKwonuxRQmoBxZJDOos5EclhIWVo5sLW8TpkW8Np0fFaeb8h1q7JQaZZ5JCw.huxz1++RVbb84yGrrrfooIWNYTqF5yMH582qolZvpV0pvQO5Qw8bO2C5Tm5jqwS32waYPFTL.f268dOTRIkfa61tMbUW0UwOeZgqhJp.u9q+5H4jSF4me9tpim5y8jZcReIoGIRDTYkUhUspUgjSNYLkoLEFz5lY5G+Sdrsnq+3MVIBH5LKihMzhJfX5s2a8kdhjrKHIKPzQW77iNKehWo+9o78Q8BGOOQ25E86E5umo2JritXiCvLEiZ8H5NWk7.7xeVea.3I69+J8263EZfZs3rset9UB+bY6oVV3x0t0iFa93b+rXCM1EbfMLgVVzCD1YqoRoegh4XqwtZQDPLhHhHhnVlzOQfFqwQH+Y4IOvSpmHhHhZMfg4jHhHhHxE4RCA.txzAYcFL5NTokkkJXXxamHhHhnVxXFhQDQDQD0njA3J5RYgboUpu7q06ZkDQDQD0RGmwBQDQDQjK5A5RVmij2lssspNiIqkqQmMYDQDQD0RGm0BQDQDQTiJ5L.CnttsbnPgZv1ICblb6HhHhHpkLFPLhHhHhHWzqnFNNNvwwQ0AUqnhJvMdi2HV3BWnZajcwZ.vNYMQDQD0p.CHFQDQDQjKxBiugggq5BVs0VKNxQNB1xV1Bdm24cTA+JTnPLPXDQDQTqJrn5SDQDQD0nrrrfOe9.PcKERud8hfAChhJpHjUVYgq3JtBUQ2ORjHMn36SDQDQTKULfXDQDQDEGRFHKgP.CCC0uC.0so24HOWB7k7wI+2yF8sKRjHpkiodMKS96BgP820qWup+kHhHhnyGbISRDQDQTbFaaa3wiGXYYoB9j7esssggggZaz6bj111vqWuvvv.VVVM34L5.YIuc4sI+YcxfgIeMHCRl72k+883wi5uOQDQDQmuL+uuIDQDQDQwRLMMwd1ydfiiCppppPf.Av.G3.wl27lwG9geHbbbv0e8WOxM2bUOlpqtZTYkUhu8a+VUvtzu+ZqsV79u+6CSSS30qWWEZ+HQhfjRJIjWd4AKKKb5SeZrsssMXXXfDRHAXaaCGGGDIRDDHP.32uebcW204ZIaJec+eK6yHhHhH5bAyPLhHhHhhyTUUUg+3e7Oh7yOejWd4gUrhUf4N24hsu8siAO3AipqtZbsW60he2u62oB9UYkUFV9xWNVvBV.t1q8ZwK7BufJPVNNN3.G3.XpScp3odpmBG9vGFACFDm9zmFKbgKDSdxSFaZSaBQhDA974CG5PGBSXBS.OwS7DnrxJCACFDVVV3QezGEicriEqe8qG..974CgCGF.tWFkDQDQDc9hYHFQDQDQwYtjK4Rvy9rOKtka4VPt4lKdi23MvS7DOAl8rmM..toa5lfWudwLm4LQu5UuPt4lKtpq5pvS+zOM95u9qQe5SebsrFMLLPUUUE5Uu5EV6ZWKRKszfkkE9vO7Cw+7e9OQN4jCdnG5gfWudgssMBGNL5YO6I1vF1.RN4jgiiC15V2J1291GxImbvi9nOJrrrfooIRHgDX8CiHhHhtfiYHFQDQDQwYbbbfWudQf.AfWudgkkEt8a+1A.TA5RFTrBKrP.TWfu74ymZYK50qWXZZ5p36OsoMMjVZoA.fJpnBbW20cAOd7fktzkptcSSSTYkUhYMqYgK8RuT.TWFqce228AgPfksrkg.AB.e97opkYxffIesQDQDQz4KFPLhHhHhhyHKT8ACFDQhDAiZTiBcricTkAVd73AYmc1viGOXUqZU3Dm3D.vcWerlZpQUX8EBAxN6rQ94mODBABGNLVvBV.NxQNBVxRVBF4HGoJXVVVVX3Ce3XbiabpmuG3Ad.TVYkgm9oeZL7gObWE3eYQ22wwAIjPBbISRDQDQWPvkLIQDQDQwYjA9xue+..HiLxPk0X.0E7ojSNYU.uNyYNCZe6aupqS533.+98Ce97gHQh.SSSjQFYn5zjqbkqDqZUqB23MdiXFyXFtdt84yGZe6aOtjK4RfiiCVwJVAV8pWMtwa7FwC7.OfJ6yzyLLYw0mACiHhHhtPgYHFQDQDQwY750qpX36wiGb5SeZUViIusScpSo1VYvnLLLT2ukkkqfVAT2xp7y+7OGETPAHojRBO6y9rnMsoMpGm9yuPHvW8UeEJnfBP6ae6wRW5RUKSRY.6jOFYmljADiHhHhtPgADiHhHhn3PFFFPHDvwwAG5PGxUQxu1ZqEm3Dm.lllnicriHkTRQszEkamPHfooIrssQvfAAPcKwwBJn.DNbXrrksLjc1YCf5BBV26d2UON.fpqtZLiYLCDLXPr7kubb4W9kCOd7.aaajc1YCCCC30qWHDB0eSOd7nVNkDQDQDc9fADiHhHhn3LNNNtx1phKtXrm8rG333.CCCDHP.ryctSXaai4O+4iDRHA.TW.uRLwDUY3krFe42ue333feyu42fssssgoMsog7yOeUFds+8ue3333paQ9xu7KiO7C+Pb228ci7yOe0qku7K+RUv5rrrfggA95u9qwt10tPjHQTOmDQDQDc9fADiHhHhn3L50BLCCCLzgNTr3EuXbxSdRDLXP7oe5mhEsnEgbxIGL9wOd..DNbXXZZhScpSAaaa7ce224ZYVt0stU7XO1igd1ydhm7IeRXaaCe97gRKsT7vO7CiLyLS01+AevGf4O+4id1ydhm4YdF0Rube6ae3QdjGAcnCcP0UKO7gOL5YO6Itlq4Zvt28tUOGDQDQDc9fEUehHhHhhyHyrKYWdb3Ce3XDiXDn6cu63xtrKC6cu6EicriEKaYKComd5prIa7ie73S9jOA974Ce9m+4XDiXDXRSZRXtyct3EewWD974CkWd4XXCaXvmOe3q9puB..lllXTiZTpfjsrksL3wiGTYkUh92+9iDRHATd4kCaaa333f7xKOUgzOojRBcpScBm3Dm.ImbxMmusQDQDQwPLDr5jRDQDQTbEY.wJt3hwvG9vwrl0rvK+xuLBEJD9lu4avkbIWBRKszTKuQYgyW93zK58xLMSOysz+cYP2jcOR4TOMLLTOeM1iS+0pssM92+6+MRKszTE2ehHhHhNevLDiHhHhn3LlllvxxREnKYWjzzzD8t28tAAzR9uxrJS93jaiL.UgCGFIjPBviGOHTnPHwDSDFFFpZBl9ymgggqfoIW9l.P8XsssgggA762O762uqfyQDQDQz4CVDFHhHhHJNib4HZYYAf5BjkrqNBTe.thDIRC5tj5+rkkkJisbbbTEeeGGGjXhI5pCQJCJlr9kI+YI4yikkk5w50qWUlnI2dt3FHhHhnKDX.wHhHhHJNiGOdvm7IeB1zl1D..13F2H1912NNxQNhqL0ROarjAFSOvY974Sss5A2ROywzeNzWpiQ+3zeNk2u91H+ctbIIhHhnKDXMDiHhHhn3L0TSM3se62F.0EnqZpoFzl1zFjd5oiq65ttl4WcDQDQDcwGCHFQDQDQwgps1ZQf.A.P8KWQVitHhHhn3EbISRDQDQTbFaaaUvvBEJjZYHxqSJQDQDEufYHFQDQDQwgz61i.vUmkjHhHhnXcbVODQDQDEGxxxxUApWFLLaa6lqWRDQDQD0jgYHFQDQDQwYzqUXVVVPHDvzzjcwQhHhHJtACHFQDQDQwozCLVzKgRhHhHhhkw.hQDQDQDQDQDQTbEVCwHhHhHhHhHhHJtBCHFQDQDQDQDQDQwUX.wHhHhHhHhHhHJtBCHFQDQDQDQDQDQwUX.wHhHhHhHhHhHJtBCHFQDQDQDQDQDQwUX.wHhHhHhHhHhHJtBCHFQDQDQDQDQDQwUX.wHhHhHhHhHhHJtBCHFQDQDQDQDQDQwUX.wHhHhHhHhHhHJtx+GbmRCyQYG0EJ.....IUjSD4pPfIH" ],
                                    "embed": 1,
                                    "forceaspect": 1,
                                    "id": "obj-22",
                                    "maxclass": "fpic",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "jit_matrix" ],
                                    "patching_rect": [ 565.0, 14.0, 639.0, 178.08196721311475 ],
                                    "pic": "miniatura1.jpeg"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 12.0,
                                    "id": "obj-12",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 205.0, 439.0, 41.0, 20.0 ],
                                    "text": "BPM"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 12.0,
                                    "id": "obj-5",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 109.0, 441.0, 74.0, 20.0 ],
                                    "text": "EVENT"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 20.0,
                                    "format": 6,
                                    "id": "obj-11",
                                    "maxclass": "flonum",
                                    "numdecimalplaces": 2,
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 251.0, 433.0, 77.0, 31.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 24.0,
                                    "id": "obj-13",
                                    "maxclass": "number",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 10.0, 433.0, 95.0, 35.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-1",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 10.0, 15.0, 144.0, 22.0 ],
                                    "text": "open miniatura1.mp3"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-2",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 2,
                                    "outlettype": [ "signal", "bang" ],
                                    "patching_rect": [ 10.0, 153.0, 47.0, 22.0 ],
                                    "text": "sfplay~"
                                }
                            },
                            {
                                "box": {
                                    "color": [ 0.08235294117647059, 0.4745098039215686, 0.08627450980392157, 1.0 ],
                                    "id": "obj-3",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "float", "float" ],
                                    "patching_rect": [ 10.0, 389.0, 260.0, 22.0 ],
                                    "text": "openscofo~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-9",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 72.0, 314.0, 136.0, 22.0 ],
                                    "text": "score miniatura1.scofo"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-10",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 92.0, 348.0, 32.0, 22.0 ],
                                    "text": "start"
                                }
                            },
                            {
                                "box": {
                                    "background": 1,
                                    "bgcolor": [ 0.792294779904983, 0.561867873693578, 0.668313324973563, 1.0 ],
                                    "bgoncolor": [ 0.55, 0.55, 0.55, 1.0 ],
                                    "fontface": 1,
                                    "hint": "",
                                    "id": "obj-17",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 39.0, 73.0, 20.0, 20.0 ],
                                    "rounded": 60.0,
                                    "saved_attribute_attributes": {
                                        "bgcolor": {
                                            "expression": "themecolor.lesson_downloading"
                                        }
                                    },
                                    "text": "4",
                                    "textcolor": [ 0.34902, 0.34902, 0.34902, 1.0 ],
                                    "textoncolor": [ 1.0, 1.0, 1.0, 1.0 ],
                                    "textovercolor": [ 0.2, 0.2, 0.2, 1.0 ],
                                    "usebgoncolor": 1,
                                    "usetextovercolor": 1
                                }
                            },
                            {
                                "box": {
                                    "background": 1,
                                    "bgcolor": [ 0.792294779904983, 0.561867873693578, 0.668313324973563, 1.0 ],
                                    "bgoncolor": [ 0.55, 0.55, 0.55, 1.0 ],
                                    "fontface": 1,
                                    "hint": "",
                                    "id": "obj-15",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 163.0, 16.0, 20.0, 20.0 ],
                                    "rounded": 60.0,
                                    "saved_attribute_attributes": {
                                        "bgcolor": {
                                            "expression": "themecolor.lesson_downloading"
                                        }
                                    },
                                    "text": "3",
                                    "textcolor": [ 0.34902, 0.34902, 0.34902, 1.0 ],
                                    "textoncolor": [ 1.0, 1.0, 1.0, 1.0 ],
                                    "textovercolor": [ 0.2, 0.2, 0.2, 1.0 ],
                                    "usebgoncolor": 1,
                                    "usetextovercolor": 1
                                }
                            },
                            {
                                "box": {
                                    "background": 1,
                                    "bgcolor": [ 0.792294779904983, 0.561867873693578, 0.668313324973563, 1.0 ],
                                    "bgoncolor": [ 0.55, 0.55, 0.55, 1.0 ],
                                    "fontface": 1,
                                    "hint": "",
                                    "id": "obj-6",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 130.0, 350.0, 20.0, 20.0 ],
                                    "rounded": 60.0,
                                    "saved_attribute_attributes": {
                                        "bgcolor": {
                                            "expression": "themecolor.lesson_downloading"
                                        }
                                    },
                                    "text": "2",
                                    "textcolor": [ 0.34902, 0.34902, 0.34902, 1.0 ],
                                    "textoncolor": [ 1.0, 1.0, 1.0, 1.0 ],
                                    "textovercolor": [ 0.2, 0.2, 0.2, 1.0 ],
                                    "usebgoncolor": 1,
                                    "usetextovercolor": 1
                                }
                            },
                            {
                                "box": {
                                    "background": 1,
                                    "bgcolor": [ 0.792294779904983, 0.561867873693578, 0.668313324973563, 1.0 ],
                                    "bgoncolor": [ 0.55, 0.55, 0.55, 1.0 ],
                                    "fontface": 1,
                                    "hint": "",
                                    "id": "obj-93",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 212.0, 315.0, 20.0, 20.0 ],
                                    "rounded": 60.0,
                                    "saved_attribute_attributes": {
                                        "bgcolor": {
                                            "expression": "themecolor.lesson_downloading"
                                        }
                                    },
                                    "text": "1",
                                    "textcolor": [ 0.34902, 0.34902, 0.34902, 1.0 ],
                                    "textoncolor": [ 1.0, 1.0, 1.0, 1.0 ],
                                    "textovercolor": [ 0.2, 0.2, 0.2, 1.0 ],
                                    "usebgoncolor": 1,
                                    "usetextovercolor": 1
                                }
                            }
                        ],
                        "lines": [
                            {
                                "patchline": {
                                    "destination": [ "obj-2", 0 ],
                                    "source": [ "obj-1", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-3", 0 ],
                                    "source": [ "obj-10", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-2", 0 ],
                                    "source": [ "obj-18", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-3", 0 ],
                                    "order": 3,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-48", 1 ],
                                    "order": 1,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-48", 0 ],
                                    "order": 2,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-7", 0 ],
                                    "order": 0,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-11", 0 ],
                                    "source": [ "obj-3", 1 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-13", 0 ],
                                    "source": [ "obj-3", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-29", 0 ],
                                    "source": [ "obj-31", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-59", 0 ],
                                    "source": [ "obj-62", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-62", 0 ],
                                    "source": [ "obj-63", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-4", 0 ],
                                    "source": [ "obj-7", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-3", 0 ],
                                    "source": [ "obj-9", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 11.0, 11.0, 79.0, 22.0 ],
                    "presentation": 1,
                    "presentation_rect": [ 114.0, 110.0, 100.0, 22.0 ],
                    "text": "p Basic"
                }
            }
        ],
        "lines": [],
        "autosave": 0,
        "oscreceiveudpport": 0
    }
}