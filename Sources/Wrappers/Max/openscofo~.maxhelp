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
        "rect": [ 118.0, 154.0, 1266.0, 852.0 ],
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
                        "rect": [ 0.0, 26.0, 1266.0, 826.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-23",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 42.159996032714844, 466.0, 188.0, 24.0 ],
                                    "text": "Check how to train AI models",
                                    "textcolor": [ 0.5176470588235295, 0.8901960784313725, 0.7647058823529411, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 16.0,
                                    "id": "obj-3",
                                    "linecount": 6,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.0, 291.0, 843.0, 114.0 ],
                                    "text": "So far, OpenScofo has been successfully used to recognize flute extended techniques as well as speech phonemes. Since the training process is generic, it can also be adapted to other instruments and sound sources by providing an appropriate labeled dataset. In general, the more diverse the training data for each label (e.g., different performers, microphones, instruments, or flute brands, etc...), the better the model will generalize and the more reliable the recognition will be in real-world performance. High-quality, balanced datasets are usually far more important than increasing model complexity."
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.592156862745098, 0.32941176470588235, 0.32941176470588235, 1.0 ],
                                    "id": "obj-8",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 10.0, 466.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-11",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 10.0, 494.0, 433.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/descriptors/ai/"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 16.0,
                                    "id": "obj-12",
                                    "linecount": 10,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.0, 84.0, 843.0, 185.0 ],
                                    "text": "The model distributed with OpenScofo is trained to recognize only a small set of techniques (currently tongue-ram, pizzicato, jet whistle, and key click). If you want to recognize additional techniques or use your own dataset, you'll need to train a new model.\n\nThis process is straightforward: organize your audio recordings into labeled folders (one folder per technique), then train the model using either the OpenScofo Python library or Pure Data with py4pd. OpenScofo automatically extracts the selected audio descriptors and builds the model for you.\n\nFor complete setup instructions, software installation, and a step-by-step training guide, refer to the online documentation."
                                }
                            },
                            {
                                "box": {
                                    "border": 0,
                                    "filename": "helpname.js",
                                    "id": "obj-4",
                                    "ignoreclick": 1,
                                    "jsarguments": [ "openscofo~" ],
                                    "maxclass": "jsui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 10.0, 10.0, 252.3199920654297, 57.599853515625 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 40.0, 859.0, 256.0, 20.0 ],
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
                                    "patching_rect": [ 14.0, 859.0, 24.0, 24.0 ]
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
                                    "patching_rect": [ 14.0, 887.0, 401.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/score/intro"
                                }
                            }
                        ],
                        "lines": [
                            {
                                "patchline": {
                                    "destination": [ "obj-29", 0 ],
                                    "source": [ "obj-31", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-11", 0 ],
                                    "source": [ "obj-8", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 84.0, 145.0, 212.0, 22.0 ],
                    "presentation": 1,
                    "presentation_linecount": 4,
                    "presentation_rect": [ 159.0, 155.0, 100.0, 62.0 ],
                    "text": "p Training Extended Techniquel Model"
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
                        "rect": [ 0.0, 26.0, 1266.0, 826.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-7",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 40.0, 739.0, 256.0, 20.0 ],
                                    "text": "Check complete docs"
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.592156862745098, 0.32941176470588235, 0.32941176470588235, 1.0 ],
                                    "id": "obj-8",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 14.0, 739.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-10",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 14.0, 767.0, 401.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/score/intro"
                                }
                            },
                            {
                                "box": {
                                    "border": 0,
                                    "filename": "helpname.js",
                                    "id": "obj-4",
                                    "ignoreclick": 1,
                                    "jsarguments": [ "openscofo~" ],
                                    "maxclass": "jsui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 10.0, 10.0, 252.3199920654297, 57.599853515625 ]
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
                                    "patching_rect": [ 597.5, 681.0, 429.0, 22.0 ],
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
                                    "patching_rect": [ 597.5, 717.0, 428.0, 85.0 ],
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
                                    "patching_rect": [ 94.0, 249.0, 339.0, 20.0 ],
                                    "text": "Receive message hello world 1 2 3."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-19",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 160.0, 490.0, 339.0, 20.0 ],
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
                                    "patching_rect": [ 25.0, 531.0, 96.0, 22.0 ],
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
                                    "patching_rect": [ 25.0, 489.0, 133.0, 22.0 ],
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
                                    "patching_rect": [ 597.5, 442.0, 428.0, 22.0 ],
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
                                    "patching_rect": [ 17.0, 438.0, 69.0, 22.0 ],
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
                                    "patching_rect": [ 80.0, 438.0, 419.0, 22.0 ],
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
                                    "patching_rect": [ 597.5, 206.0, 389.0, 22.0 ],
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
                                    "patching_rect": [ 14.0, 290.0, 96.0, 22.0 ],
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
                                    "patching_rect": [ 14.0, 248.0, 75.0, 22.0 ],
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
                                    "patching_rect": [ 10.0, 206.0, 69.0, 22.0 ],
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
                                    "patching_rect": [ 80.0, 206.0, 250.0, 22.0 ],
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
                                    "patching_rect": [ 10.0, 91.0, 616.0, 42.0 ],
                                    "text": "OpenScofo supports ACTIONS that are triggered by EVENTS. In other words, when an event is detected, it triggers one or more actions."
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
                                    "destination": [ "obj-10", 0 ],
                                    "source": [ "obj-8", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 42.5, 79.0, 57.0, 22.0 ],
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
                        "rect": [ 0.0, 26.0, 1266.0, 826.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-13",
                                    "linecount": 3,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 301.0, 486.5, 465.0, 51.0 ],
                                    "text": "To use OpenScofo's descriptors, you need to know the ID assigned to each one. A complete list of descriptor IDs is available in the OpenScofo Descriptors documentation.",
                                    "textcolor": [ 1.0, 1.0, 1.0, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 16.0,
                                    "id": "obj-12",
                                    "linecount": 3,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.0, 84.0, 839.0, 60.0 ],
                                    "text": "OpenScofo also provides a set of audio descriptors, primarily intended for training and running the AI model for extended technique recognition. However, they can also be used independently in your own applications. Most of these descriptors are compatible with those available in the Python library librosa."
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.592156862745098, 0.32941176470588235, 0.32941176470588235, 1.0 ],
                                    "id": "obj-8",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 28.0, 711.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-11",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 28.0, 739.0, 407.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/descriptors/"
                                }
                            },
                            {
                                "box": {
                                    "border": 0,
                                    "filename": "helpname.js",
                                    "id": "obj-4",
                                    "ignoreclick": 1,
                                    "jsarguments": [ "openscofo~" ],
                                    "maxclass": "jsui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 10.0, 10.0, 252.3199920654297, 57.599853515625 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 40.0, 859.0, 256.0, 20.0 ],
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
                                    "patching_rect": [ 14.0, 859.0, 24.0, 24.0 ]
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
                                    "patching_rect": [ 14.0, 887.0, 401.0, 35.0 ],
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
                                    "patching_rect": [ 77.0, 711.0, 289.0, 24.0 ],
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
                                    "patching_rect": [ 440.40000396966934, 618.0, 322.0, 106.0 ],
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
                                    "patching_rect": [ 301.0, 269.0, 619.200009226799, 22.0 ],
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
                                    "patching_rect": [ 372.0, 410.0, 548.0000081658363, 20.0 ],
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
                                    "patching_rect": [ 301.0, 410.0, 65.0, 20.0 ],
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
                                    "patching_rect": [ 301.0, 364.0, 64.80000096559525, 20.0 ],
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
                                    "patching_rect": [ 301.0, 339.0, 64.80000096559525, 20.0 ],
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
                                    "patching_rect": [ 301.0, 312.0, 64.80000096559525, 20.0 ],
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
                                    "patching_rect": [ 372.0, 364.0, 648.0000096559525, 33.0 ],
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
                                    "patching_rect": [ 372.0, 339.0, 606.4000090360641, 20.0 ],
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
                                    "patching_rect": [ 372.0, 312.0, 496.00000739097595, 20.0 ],
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
                                    "patching_rect": [ 522.4000039696693, 590.0, 68.00000101327896, 22.0 ]
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
                                    "patching_rect": [ 359.40000396966934, 590.0, 68.00000101327896, 22.0 ]
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
                                    "patching_rect": [ 277.40000396966934, 590.0, 68.00000101327896, 22.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-6",
                                    "maxclass": "newobj",
                                    "numinlets": 5,
                                    "numoutlets": 5,
                                    "outlettype": [ "", "", "", "", "" ],
                                    "patching_rect": [ 277.40000396966934, 551.0, 345.66666116714475, 22.0 ],
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
                                    "patching_rect": [ 65.0, 252.0, 113.0, 24.0 ],
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
                                    "patching_rect": [ 193.0, 195.0, 342.0, 24.0 ],
                                    "text": "Load recording. For live performances, replace with [adc~]"
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
                                    "patching_rect": [ 42.0, 281.0, 196.0, 19.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-1",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 10.0, 196.0, 144.0, 22.0 ],
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
                                    "patching_rect": [ 10.0, 334.0, 47.0, 22.0 ],
                                    "text": "sfplay~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-21",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 0,
                                    "patching_rect": [ 32.0, 393.0, 35.0, 22.0 ],
                                    "text": "dac~"
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
                                    "patching_rect": [ 11.0, 501.0, 285.40000396966934, 22.0 ],
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
                                    "id": "obj-17",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 39.0, 254.0, 20.0, 20.0 ],
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
                                    "id": "obj-15",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 163.0, 197.0, 20.0, 20.0 ],
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
                                    "destination": [ "obj-11", 0 ],
                                    "source": [ "obj-8", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 56.0, 114.0, 104.0, 22.0 ],
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
                        "rect": [ 0.0, 26.0, 1266.0, 826.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-23",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 480.0, 694.0, 526.0, 24.0 ],
                                    "text": "Import Tradicional scores using MusicXML from apps like Sibelius, MuseScore and others!",
                                    "textcolor": [ 0.5176470588235295, 0.8901960784313725, 0.7647058823529411, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "bgcolor": [ 0.5176470588235295, 0.8901960784313725, 0.7647058823529411, 1.0 ],
                                    "id": "obj-15",
                                    "maxclass": "button",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 452.0, 694.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-17",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 452.0, 722.0, 401.0, 35.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/Editor/"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 16.0,
                                    "id": "obj-8",
                                    "linecount": 2,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 10.0, 84.0, 839.0, 42.0 ],
                                    "text": "OpenScofo uses a text file to define the score. In this file, you specify the musical events and the actions associated with each event. This section explains how to define those events."
                                }
                            },
                            {
                                "box": {
                                    "border": 0,
                                    "filename": "helpname.js",
                                    "id": "obj-1",
                                    "ignoreclick": 1,
                                    "jsarguments": [ "openscofo~" ],
                                    "maxclass": "jsui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 10.0, 10.0, 252.3199920654297, 57.599853515625 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 39.0, 696.0, 256.0, 20.0 ],
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
                                    "patching_rect": [ 13.0, 696.0, 24.0, 24.0 ]
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
                                    "patching_rect": [ 13.0, 724.0, 401.0, 35.0 ],
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
                                    "patching_rect": [ 82.0, 520.0, 764.0, 163.0 ],
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
                                    "patching_rect": [ 81.0, 335.0, 766.0, 147.0 ],
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
                                    "patching_rect": [ 81.0, 293.0, 904.0, 22.0 ],
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
                                    "patching_rect": [ 81.0, 247.0, 905.0, 22.0 ],
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
                                    "patching_rect": [ 81.0, 209.0, 906.0, 22.0 ],
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
                                    "patching_rect": [ 81.0, 175.0, 632.0, 22.0 ],
                                    "text": "NOTE events describes tradicional pitches. It must be defined as NOTE <PITCH> <DURATION>. "
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
                                    "patching_rect": [ 11.0, 520.0, 70.0, 22.0 ],
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
                                    "patching_rect": [ 10.0, 335.0, 70.0, 22.0 ],
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
                                    "patching_rect": [ 10.0, 293.0, 69.0, 22.0 ],
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
                                    "patching_rect": [ 10.0, 247.0, 72.0, 22.0 ],
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
                                    "patching_rect": [ 10.0, 209.0, 71.0, 22.0 ],
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
                                    "patching_rect": [ 10.0, 175.0, 69.0, 22.0 ],
                                    "text": "NOTE",
                                    "textcolor": [ 0.8156862745098039, 0.12941176470588237, 0.12941176470588237, 1.0 ]
                                }
                            }
                        ],
                        "lines": [
                            {
                                "patchline": {
                                    "destination": [ "obj-17", 0 ],
                                    "source": [ "obj-15", 0 ]
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
                    "patching_rect": [ 25.5, 44.0, 50.0, 22.0 ],
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
                        "rect": [ 118.0, 180.0, 1266.0, 826.0 ],
                        "toolbarvisible": 0,
                        "enablehscroll": 0,
                        "enablevscroll": 0,
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-25",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 106.0, 414.5, 113.0, 24.0 ],
                                    "text": "Turn Audio On"
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-22",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 290.0, 361.0, 203.0, 24.0 ],
                                    "text": "Actions trigger from the score"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-21",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 266.0, 397.0, 141.0, 22.0 ],
                                    "text": "print pitch @popup 1"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-4",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 237.0, 427.0, 123.0, 22.0 ],
                                    "text": "print delay @popup 1"
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-14",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 772.0, 43.5, 452.0, 24.0 ],
                                    "text": "You can convert musicxml to this text format using OpenScofo converter tool."
                                }
                            },
                            {
                                "box": {
                                    "filename": "score_renderer.js",
                                    "id": "obj-16",
                                    "maxclass": "jsui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 620.0, 442.0, 639.0, 314.0 ]
                                }
                            },
                            {
                                "box": {
                                    "border": 0,
                                    "filename": "helpdetails.js",
                                    "id": "obj-8",
                                    "ignoreclick": 1,
                                    "jsarguments": [ "openscofo~" ],
                                    "maxclass": "jsui",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 10.0, 10.0, 597.0, 112.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 45.0, 749.0, 256.0, 20.0 ],
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
                                    "patching_rect": [ 19.0, 749.0, 24.0, 24.0 ]
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
                                    "patching_rect": [ 19.0, 777.0, 401.0, 35.0 ],
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
                                    "patching_rect": [ 72.0, 233.0, 113.0, 24.0 ],
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
                                    "patching_rect": [ 200.0, 176.0, 349.0, 24.0 ],
                                    "text": "Load recording. For live performances, replace with [adc~]"
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-20",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 163.0, 510.0, 313.0, 24.0 ],
                                    "text": "Start to follow (press again if you want to play again)"
                                }
                            },
                            {
                                "box": {
                                    "bubble": 1,
                                    "id": "obj-19",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 250.0, 475.0, 114.0, 24.0 ],
                                    "text": "Load the score"
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
                                    "patching_rect": [ 49.0, 262.0, 196.0, 19.0 ]
                                }
                            },
                            {
                                "box": {
                                    "color": [ 0.11764705882352941, 0.17647058823529413, 0.7764705882352941, 1.0 ],
                                    "id": "obj-7",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "signal", "", "" ],
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
                                        "rect": [ 191.0, 203.0, 1149.0, 780.0 ],
                                        "boxes": [
                                            {
                                                "box": {
                                                    "bubble": 1,
                                                    "id": "obj-14",
                                                    "linecount": 2,
                                                    "maxclass": "comment",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 280.0, 3.5, 334.0, 37.0 ],
                                                    "text": "Check miniatura1.scofo, you will see \"sendto delay\" with the parameters received by this receiver."
                                                }
                                            },
                                            {
                                                "box": {
                                                    "bubble": 1,
                                                    "id": "obj-23",
                                                    "linecount": 2,
                                                    "maxclass": "comment",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 768.0, 3.5, 334.0, 37.0 ],
                                                    "text": "Check miniatura1.scofo, you will see \"sendto pitchshift\" with the parameters received by this receiver."
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-12",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
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
                                                        "rect": [ 84.0, 149.0, 1000.0, 780.0 ],
                                                        "boxes": [
                                                            {
                                                                "box": {
                                                                    "id": "obj-22",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 533.0, 240.0, 70.0, 22.0 ],
                                                                    "text": "loadmess 1"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-21",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 382.0, 365.0, 42.0, 22.0 ],
                                                                    "text": "mute~"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-12",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 382.0, 240.0, 134.0, 22.0 ],
                                                                    "text": "if $f1 == 0 then 1 else 0"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-8",
                                                                    "maxclass": "toggle",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "int" ],
                                                                    "parameter_enable": 0,
                                                                    "patching_rect": [ 382.0, 318.0, 24.0, 24.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-6",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 2,
                                                                    "outlettype": [ "float", "float" ],
                                                                    "patching_rect": [ 341.0, 175.0, 60.0, 22.0 ],
                                                                    "text": "t f f"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-3",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 196.0, 175.0, 59.0, 22.0 ],
                                                                    "text": "transratio"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-1",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 2,
                                                                    "outlettype": [ "float", "float" ],
                                                                    "patching_rect": [ 196.0, 103.0, 164.0, 22.0 ],
                                                                    "text": "unpack f f"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-4",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "signal" ],
                                                                    "patching_rect": [ 50.0, 365.0, 29.5, 22.0 ],
                                                                    "text": "*~"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "fontname": "Arial",
                                                                    "fontsize": 13.0,
                                                                    "id": "obj-14",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "signal" ],
                                                                    "patching_rect": [ 50.0, 290.0, 165.0, 23.0 ],
                                                                    "text": "pfft~ gizmo_loadme 2048 2"
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
                                                                    "patching_rect": [ 50.0, 40.0, 30.0, 30.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "comment": "",
                                                                    "id": "obj-9",
                                                                    "index": 2,
                                                                    "maxclass": "inlet",
                                                                    "numinlets": 0,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 196.0, 40.0, 30.0, 30.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "comment": "",
                                                                    "id": "obj-18",
                                                                    "index": 1,
                                                                    "maxclass": "outlet",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 0,
                                                                    "patching_rect": [ 50.0, 421.0, 30.0, 30.0 ]
                                                                }
                                                            }
                                                        ],
                                                        "lines": [
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-3", 0 ],
                                                                    "source": [ "obj-1", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-6", 0 ],
                                                                    "source": [ "obj-1", 1 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-8", 0 ],
                                                                    "source": [ "obj-12", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-4", 0 ],
                                                                    "source": [ "obj-14", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 0 ],
                                                                    "source": [ "obj-2", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 0 ],
                                                                    "midpoints": [ 391.5, 397.0, 235.45703125, 397.0, 235.45703125, 261.984375, 59.5, 261.984375 ],
                                                                    "source": [ "obj-21", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-8", 0 ],
                                                                    "source": [ "obj-22", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 1 ],
                                                                    "source": [ "obj-3", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-18", 0 ],
                                                                    "source": [ "obj-4", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-12", 0 ],
                                                                    "source": [ "obj-6", 1 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-4", 1 ],
                                                                    "midpoints": [ 350.5, 339.94921875, 70.0, 339.94921875 ],
                                                                    "source": [ "obj-6", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-21", 0 ],
                                                                    "source": [ "obj-8", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-1", 0 ],
                                                                    "source": [ "obj-9", 0 ]
                                                                }
                                                            }
                                                        ]
                                                    },
                                                    "patching_rect": [ 886.0, 226.0, 53.0, 22.0 ],
                                                    "text": "p voice4"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "comment": "",
                                                    "id": "obj-10",
                                                    "index": 3,
                                                    "maxclass": "outlet",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 644.0, 48.0, 30.0, 30.0 ]
                                                }
                                            },
                                            {
                                                "box": {
                                                    "comment": "",
                                                    "id": "obj-8",
                                                    "index": 2,
                                                    "maxclass": "outlet",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 298.0, 58.0, 30.0, 30.0 ]
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-7",
                                                    "maxclass": "message",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 171.0, 232.0, 45.0, 22.0 ],
                                                    "text": "$1 200"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-5",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 2,
                                                    "outlettype": [ "signal", "bang" ],
                                                    "patching_rect": [ 171.0, 274.0, 34.0, 22.0 ],
                                                    "text": "line~"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-21",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
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
                                                        "rect": [ 84.0, 149.0, 1000.0, 780.0 ],
                                                        "boxes": [
                                                            {
                                                                "box": {
                                                                    "id": "obj-22",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 533.0, 240.0, 70.0, 22.0 ],
                                                                    "text": "loadmess 1"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-21",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 382.0, 365.0, 42.0, 22.0 ],
                                                                    "text": "mute~"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-12",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 382.0, 240.0, 134.0, 22.0 ],
                                                                    "text": "if $f1 == 0 then 1 else 0"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-8",
                                                                    "maxclass": "toggle",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "int" ],
                                                                    "parameter_enable": 0,
                                                                    "patching_rect": [ 382.0, 318.0, 24.0, 24.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-6",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 2,
                                                                    "outlettype": [ "float", "float" ],
                                                                    "patching_rect": [ 341.0, 175.0, 60.0, 22.0 ],
                                                                    "text": "t f f"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-3",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 196.0, 175.0, 59.0, 22.0 ],
                                                                    "text": "transratio"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-1",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 2,
                                                                    "outlettype": [ "float", "float" ],
                                                                    "patching_rect": [ 196.0, 103.0, 164.0, 22.0 ],
                                                                    "text": "unpack f f"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-4",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "signal" ],
                                                                    "patching_rect": [ 50.0, 365.0, 29.5, 22.0 ],
                                                                    "text": "*~"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "fontname": "Arial",
                                                                    "fontsize": 13.0,
                                                                    "id": "obj-14",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "signal" ],
                                                                    "patching_rect": [ 50.0, 290.0, 165.0, 23.0 ],
                                                                    "text": "pfft~ gizmo_loadme 2048 2"
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
                                                                    "patching_rect": [ 50.0, 40.0, 30.0, 30.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "comment": "",
                                                                    "id": "obj-9",
                                                                    "index": 2,
                                                                    "maxclass": "inlet",
                                                                    "numinlets": 0,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 196.0, 40.0, 30.0, 30.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "comment": "",
                                                                    "id": "obj-18",
                                                                    "index": 1,
                                                                    "maxclass": "outlet",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 0,
                                                                    "patching_rect": [ 50.0, 421.0, 30.0, 30.0 ]
                                                                }
                                                            }
                                                        ],
                                                        "lines": [
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-3", 0 ],
                                                                    "source": [ "obj-1", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-6", 0 ],
                                                                    "source": [ "obj-1", 1 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-8", 0 ],
                                                                    "source": [ "obj-12", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-4", 0 ],
                                                                    "source": [ "obj-14", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 0 ],
                                                                    "source": [ "obj-2", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 0 ],
                                                                    "midpoints": [ 391.5, 397.0, 235.45703125, 397.0, 235.45703125, 261.984375, 59.5, 261.984375 ],
                                                                    "source": [ "obj-21", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-8", 0 ],
                                                                    "source": [ "obj-22", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 1 ],
                                                                    "source": [ "obj-3", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-18", 0 ],
                                                                    "source": [ "obj-4", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-12", 0 ],
                                                                    "source": [ "obj-6", 1 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-4", 1 ],
                                                                    "midpoints": [ 350.5, 339.94921875, 70.0, 339.94921875 ],
                                                                    "source": [ "obj-6", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-21", 0 ],
                                                                    "source": [ "obj-8", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-1", 0 ],
                                                                    "source": [ "obj-9", 0 ]
                                                                }
                                                            }
                                                        ]
                                                    },
                                                    "patching_rect": [ 813.0, 226.0, 53.0, 22.0 ],
                                                    "text": "p voice3"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-20",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
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
                                                        "rect": [ 84.0, 149.0, 1000.0, 780.0 ],
                                                        "boxes": [
                                                            {
                                                                "box": {
                                                                    "id": "obj-22",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 533.0, 240.0, 70.0, 22.0 ],
                                                                    "text": "loadmess 1"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-21",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 382.0, 365.0, 42.0, 22.0 ],
                                                                    "text": "mute~"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-12",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 382.0, 240.0, 134.0, 22.0 ],
                                                                    "text": "if $f1 == 0 then 1 else 0"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-8",
                                                                    "maxclass": "toggle",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "int" ],
                                                                    "parameter_enable": 0,
                                                                    "patching_rect": [ 382.0, 318.0, 24.0, 24.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-6",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 2,
                                                                    "outlettype": [ "float", "float" ],
                                                                    "patching_rect": [ 341.0, 175.0, 60.0, 22.0 ],
                                                                    "text": "t f f"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-3",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 196.0, 175.0, 59.0, 22.0 ],
                                                                    "text": "transratio"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-1",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 2,
                                                                    "outlettype": [ "float", "float" ],
                                                                    "patching_rect": [ 196.0, 103.0, 164.0, 22.0 ],
                                                                    "text": "unpack f f"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-4",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "signal" ],
                                                                    "patching_rect": [ 50.0, 365.0, 29.5, 22.0 ],
                                                                    "text": "*~"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "fontname": "Arial",
                                                                    "fontsize": 13.0,
                                                                    "id": "obj-14",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "signal" ],
                                                                    "patching_rect": [ 50.0, 290.0, 165.0, 23.0 ],
                                                                    "text": "pfft~ gizmo_loadme 2048 2"
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
                                                                    "patching_rect": [ 50.0, 40.0, 30.0, 30.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "comment": "",
                                                                    "id": "obj-9",
                                                                    "index": 2,
                                                                    "maxclass": "inlet",
                                                                    "numinlets": 0,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 196.0, 40.0, 30.0, 30.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "comment": "",
                                                                    "id": "obj-18",
                                                                    "index": 1,
                                                                    "maxclass": "outlet",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 0,
                                                                    "patching_rect": [ 50.0, 421.0, 30.0, 30.0 ]
                                                                }
                                                            }
                                                        ],
                                                        "lines": [
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-3", 0 ],
                                                                    "source": [ "obj-1", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-6", 0 ],
                                                                    "source": [ "obj-1", 1 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-8", 0 ],
                                                                    "source": [ "obj-12", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-4", 0 ],
                                                                    "source": [ "obj-14", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 0 ],
                                                                    "source": [ "obj-2", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 0 ],
                                                                    "midpoints": [ 391.5, 397.0, 235.45703125, 397.0, 235.45703125, 261.984375, 59.5, 261.984375 ],
                                                                    "source": [ "obj-21", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-8", 0 ],
                                                                    "source": [ "obj-22", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 1 ],
                                                                    "source": [ "obj-3", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-18", 0 ],
                                                                    "source": [ "obj-4", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-12", 0 ],
                                                                    "source": [ "obj-6", 1 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-4", 1 ],
                                                                    "midpoints": [ 350.5, 339.94921875, 70.0, 339.94921875 ],
                                                                    "source": [ "obj-6", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-21", 0 ],
                                                                    "source": [ "obj-8", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-1", 0 ],
                                                                    "source": [ "obj-9", 0 ]
                                                                }
                                                            }
                                                        ]
                                                    },
                                                    "patching_rect": [ 740.0, 226.0, 53.0, 22.0 ],
                                                    "text": "p voice2"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-19",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
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
                                                        "rect": [ 84.0, 149.0, 1000.0, 780.0 ],
                                                        "boxes": [
                                                            {
                                                                "box": {
                                                                    "id": "obj-22",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 533.0, 240.0, 70.0, 22.0 ],
                                                                    "text": "loadmess 1"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-21",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 382.0, 365.0, 42.0, 22.0 ],
                                                                    "text": "mute~"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-12",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 382.0, 240.0, 134.0, 22.0 ],
                                                                    "text": "if $f1 == 0 then 1 else 0"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-8",
                                                                    "maxclass": "toggle",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "int" ],
                                                                    "parameter_enable": 0,
                                                                    "patching_rect": [ 382.0, 318.0, 24.0, 24.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-6",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 2,
                                                                    "outlettype": [ "float", "float" ],
                                                                    "patching_rect": [ 341.0, 175.0, 60.0, 22.0 ],
                                                                    "text": "t f f"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-3",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 196.0, 175.0, 59.0, 22.0 ],
                                                                    "text": "transratio"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-1",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 2,
                                                                    "outlettype": [ "float", "float" ],
                                                                    "patching_rect": [ 196.0, 103.0, 164.0, 22.0 ],
                                                                    "text": "unpack f f"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "id": "obj-4",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "signal" ],
                                                                    "patching_rect": [ 50.0, 365.0, 29.5, 22.0 ],
                                                                    "text": "*~"
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "fontname": "Arial",
                                                                    "fontsize": 13.0,
                                                                    "id": "obj-14",
                                                                    "maxclass": "newobj",
                                                                    "numinlets": 2,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "signal" ],
                                                                    "patching_rect": [ 50.0, 290.0, 165.0, 23.0 ],
                                                                    "text": "pfft~ gizmo_loadme 2048 2"
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
                                                                    "patching_rect": [ 50.0, 40.0, 30.0, 30.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "comment": "",
                                                                    "id": "obj-9",
                                                                    "index": 2,
                                                                    "maxclass": "inlet",
                                                                    "numinlets": 0,
                                                                    "numoutlets": 1,
                                                                    "outlettype": [ "" ],
                                                                    "patching_rect": [ 196.0, 40.0, 30.0, 30.0 ]
                                                                }
                                                            },
                                                            {
                                                                "box": {
                                                                    "comment": "",
                                                                    "id": "obj-18",
                                                                    "index": 1,
                                                                    "maxclass": "outlet",
                                                                    "numinlets": 1,
                                                                    "numoutlets": 0,
                                                                    "patching_rect": [ 50.0, 421.0, 30.0, 30.0 ]
                                                                }
                                                            }
                                                        ],
                                                        "lines": [
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-3", 0 ],
                                                                    "source": [ "obj-1", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-6", 0 ],
                                                                    "source": [ "obj-1", 1 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-8", 0 ],
                                                                    "source": [ "obj-12", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-4", 0 ],
                                                                    "source": [ "obj-14", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 0 ],
                                                                    "source": [ "obj-2", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 0 ],
                                                                    "midpoints": [ 391.5, 397.0, 235.45703125, 397.0, 235.45703125, 261.984375, 59.5, 261.984375 ],
                                                                    "source": [ "obj-21", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-8", 0 ],
                                                                    "source": [ "obj-22", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-14", 1 ],
                                                                    "source": [ "obj-3", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-18", 0 ],
                                                                    "source": [ "obj-4", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-12", 0 ],
                                                                    "source": [ "obj-6", 1 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-4", 1 ],
                                                                    "midpoints": [ 350.5, 339.94921875, 70.0, 339.94921875 ],
                                                                    "source": [ "obj-6", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-21", 0 ],
                                                                    "source": [ "obj-8", 0 ]
                                                                }
                                                            },
                                                            {
                                                                "patchline": {
                                                                    "destination": [ "obj-1", 0 ],
                                                                    "source": [ "obj-9", 0 ]
                                                                }
                                                            }
                                                        ]
                                                    },
                                                    "patching_rect": [ 667.0, 226.0, 53.0, 22.0 ],
                                                    "text": "p voice1"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-17",
                                                    "maxclass": "newobj",
                                                    "numinlets": 5,
                                                    "numoutlets": 5,
                                                    "outlettype": [ "", "", "", "", "" ],
                                                    "patching_rect": [ 701.0, 66.0, 311.0, 22.0 ],
                                                    "text": "route 1 2 3 4"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-16",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 701.0, 11.0, 63.0, 22.0 ],
                                                    "text": "r pitchshift"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-4",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 17.0, 384.0, 29.5, 22.0 ],
                                                    "text": "*~"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-3",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 3,
                                                    "outlettype": [ "float", "float", "float" ],
                                                    "patching_rect": [ 222.5, 62.0, 67.0, 22.0 ],
                                                    "text": "unpack f f f"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-1",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 222.5, 11.0, 45.0, 22.0 ],
                                                    "text": "r delay"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-11",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 249.5, 341.0, 40.0, 22.0 ],
                                                    "text": "*~ 0.6"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "comment": "",
                                                    "id": "obj-9",
                                                    "index": 1,
                                                    "maxclass": "outlet",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 16.75, 619.0, 30.0, 30.0 ]
                                                }
                                            },
                                            {
                                                "box": {
                                                    "fontname": "Arial",
                                                    "fontsize": 13.0,
                                                    "id": "obj-22",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 17.0, 267.0, 80.0, 23.0 ],
                                                    "text": "tapout~ 200"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "fontname": "Arial",
                                                    "fontsize": 13.0,
                                                    "id": "obj-30",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "tapconnect" ],
                                                    "patching_rect": [ 17.0, 176.0, 79.0, 23.0 ],
                                                    "text": "tapin~ 1000"
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
                                            }
                                        ],
                                        "lines": [
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-3", 0 ],
                                                    "order": 1,
                                                    "source": [ "obj-1", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-8", 0 ],
                                                    "order": 0,
                                                    "source": [ "obj-1", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-30", 0 ],
                                                    "midpoints": [ 259.0, 373.0, 130.5, 373.0, 130.5, 113.0, 26.5, 113.0 ],
                                                    "source": [ "obj-11", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-9", 0 ],
                                                    "midpoints": [ 895.5, 559.7578125, 26.25, 559.7578125 ],
                                                    "source": [ "obj-12", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-10", 0 ],
                                                    "order": 1,
                                                    "source": [ "obj-16", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-17", 0 ],
                                                    "order": 0,
                                                    "source": [ "obj-16", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-12", 1 ],
                                                    "source": [ "obj-17", 3 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-19", 1 ],
                                                    "source": [ "obj-17", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-20", 1 ],
                                                    "source": [ "obj-17", 1 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-21", 1 ],
                                                    "source": [ "obj-17", 2 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-9", 0 ],
                                                    "midpoints": [ 676.5, 559.03515625, 26.25, 559.03515625 ],
                                                    "source": [ "obj-19", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-12", 0 ],
                                                    "midpoints": [ 26.5, 98.25390625, 895.5, 98.25390625 ],
                                                    "order": 0,
                                                    "source": [ "obj-2", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-19", 0 ],
                                                    "midpoints": [ 26.5, 96.9140625, 676.5, 96.9140625 ],
                                                    "order": 3,
                                                    "source": [ "obj-2", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-20", 0 ],
                                                    "midpoints": [ 26.5, 96.53515625, 749.5, 96.53515625 ],
                                                    "order": 2,
                                                    "source": [ "obj-2", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-21", 0 ],
                                                    "midpoints": [ 26.5, 97.734375, 822.5, 97.734375 ],
                                                    "order": 1,
                                                    "source": [ "obj-2", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-30", 0 ],
                                                    "order": 4,
                                                    "source": [ "obj-2", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-9", 0 ],
                                                    "midpoints": [ 749.5, 559.38671875, 26.25, 559.38671875 ],
                                                    "source": [ "obj-20", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-9", 0 ],
                                                    "midpoints": [ 822.5, 559.77734375, 26.25, 559.77734375 ],
                                                    "source": [ "obj-21", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-11", 0 ],
                                                    "midpoints": [ 26.5, 315.5, 259.0, 315.5 ],
                                                    "order": 0,
                                                    "source": [ "obj-22", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-4", 0 ],
                                                    "order": 1,
                                                    "source": [ "obj-22", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-11", 1 ],
                                                    "midpoints": [ 280.0, 204.0, 280.0, 204.0 ],
                                                    "source": [ "obj-3", 2 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-30", 0 ],
                                                    "midpoints": [ 232.0, 158.90625, 26.5, 158.90625 ],
                                                    "source": [ "obj-3", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-7", 0 ],
                                                    "midpoints": [ 256.0, 200.921875, 180.5, 200.921875 ],
                                                    "source": [ "obj-3", 1 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-22", 0 ],
                                                    "source": [ "obj-30", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-9", 0 ],
                                                    "midpoints": [ 26.5, 447.0, 26.25, 447.0 ],
                                                    "source": [ "obj-4", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-4", 1 ],
                                                    "midpoints": [ 180.5, 340.0, 37.0, 340.0 ],
                                                    "source": [ "obj-5", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-5", 0 ],
                                                    "source": [ "obj-7", 0 ]
                                                }
                                            }
                                        ]
                                    },
                                    "patching_rect": [ 208.0, 362.0, 77.0, 22.0 ],
                                    "text": "p ACTIONS"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-63",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "bang" ],
                                    "patching_rect": [ 620.0, 15.5, 58.0, 22.0 ],
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
                                    "patching_rect": [ 620.0, 44.5, 145.0, 22.0 ],
                                    "text": "filename miniatura1.scofo"
                                }
                            },
                            {
                                "box": {
                                    "code": "BPM 100\n\nONNXMODEL ai-flute-model.onnx\nONNXDESCRIPTORS mfcc logmel centroid flatness hfr flux zcr irregularity kurtosis\n\n// Measure number 1\nUTECH jet_whistle 1\nREST 0.5\t\nNOTE Bb4 1.5 // tied\n    sendto pitchshift [1 3.81 0.4]\n\tdelay 0.33 tempo sendto pitchshift [2 7.02 0.4]\n\tdelay 0.66 tempo sendto pitchshift [3 9.69 0.4]\n\tdelay 1 tempo sendto pitchshift [4 -9.69 0.4]\nNOTE A4 0.5\nREST 0.5\n    sendto pitchshift [1 3.81 0]\n\tsendto pitchshift [2 7.02 0]\n\tsendto pitchshift [3 9.69 0]\n\tsendto pitchshift [4 -9.69 0]\n\n// Measure number 2\nUTECH jet_whistle 1\nREST 0.5\nNOTE Bb4 1.5 // tied\n\tsendto pitchshift [1 0.81 0.4]\n\tdelay 0.33 tempo sendto pitchshift [2 2.02 0.4]\n\tdelay 0.66 tempo sendto pitchshift [3 3.69 0.4]\n\tdelay 1 tempo sendto pitchshift [4 -3.69 0.4]\nNOTE A4 0.5\nREST 0.5\n    sendto pitchshift [1 3.81 0]\n\tsendto pitchshift [2 7.02 0]\n\tsendto pitchshift [3 9.69 0]\n\tsendto pitchshift [4 -9.69 0]\n\n// Measure number 3\nPTECH pizzicato D4 0.5\n\tsendto delay [150 1 0.5]\nPTECH pizzicato A4 0.5\nREST 0.5\n\tsendto delay [150 0 0]\n\nPTECH pizzicato D4 0.5\n\tsendto delay [123 1 0.5]\nPTECH pizzicato A4 0.5\nREST 0.5\n\tsendto delay [113 0 0]\nPTECH pizzicato D4 0.5\n\tsendto delay [133 1 0.5]\nPTECH pizzicato A4 0.5\n\n// Measure number 4\nREST 0.5\n\tsendto delay [133 0 0]\n\nPTECH pizzicato D4 0.5\n\tsendto delay [50 1 0.5]\nPTECH pizzicato A4 0.5\nREST 0.5\n\tsendto delay [50 0 0]\nUTECH jet_whistle 1\nREST 1\n\n// Measure number 5\nPTECH pizzicato D4 0.5\n\tsendto delay [66 1 0.5]\nPTECH pizzicato A4 0.5\nREST 0.5\n\tsendto delay [66 1 0]\nPTECH pizzicato D4 0.5\n\tsendto delay [88 1 0.5]\nPTECH pizzicato Bb4 0.5\nREST 0.5\n\tsendto delay [88 1 0]\nPTECH pizzicato D4 0.5\n\tsendto delay [100 1 0.5]\nPTECH pizzicato B4 0.5\n\n// Measure number 6\nREST 0.5\n\tsendto delay [100 1 0]\nPTECH pizzicato D4 0.5\n\tsendto delay [133 1 0.5]\nPTECH pizzicato Bb4 0.5\nREST 0.5\n\tsendto delay [100 0 0]\nUTECH jet_whistle 1\nUTECH jet_whistle 1\n\n// Measure number 7\nUTECH jet_whistle 1\nREST 3\n\n// Measure number 8\nNOTE D4 2\n\tsendto pitchshift [1 3.81 0.4]\n\tdelay 0.33 tempo sendto pitchshift [2 7.02 0.4]\n\tdelay 0.66 tempo sendto pitchshift [3 9.69 0.4]\n\tdelay 1 tempo sendto pitchshift [4 -9.69 0.4]\n\nREST 2\n\tsendto pitchshift [1 3.81 0]\n\tsendto pitchshift [2 7.02 0]\n\tsendto pitchshift [3 9.69 0]\n\tsendto pitchshift [4 -9.69 0]\n\n",
                                    "filename": "miniatura1-max.scofo",
                                    "fontface": 0,
                                    "fontname": "<Monospaced>",
                                    "fontsize": 12.0,
                                    "id": "obj-59",
                                    "maxclass": "text.codebox",
                                    "numinlets": 2,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "" ],
                                    "patching_rect": [ 620.0, 70.5, 639.0, 349.0 ],
                                    "saved_object_attributes": {
                                        "parameter_enable": 0,
                                        "parameter_mappable": 0
                                    },
                                    "textfile": {
                                        "filename": "miniatura1-max.scofo",
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
                                    "patching_rect": [ 33.5, 404.0, 45.0, 45.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 12.0,
                                    "id": "obj-12",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 212.0, 601.0, 41.0, 20.0 ],
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
                                    "patching_rect": [ 116.0, 603.0, 74.0, 20.0 ],
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
                                    "patching_rect": [ 258.0, 595.0, 77.0, 31.0 ]
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
                                    "patching_rect": [ 17.0, 595.0, 95.0, 35.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-1",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 17.0, 177.0, 144.0, 22.0 ],
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
                                    "patching_rect": [ 17.0, 315.0, 47.0, 22.0 ],
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
                                    "patching_rect": [ 17.0, 551.0, 260.0, 22.0 ],
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
                                    "patching_rect": [ 50.5, 476.0, 156.0, 22.0 ],
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
                                    "patching_rect": [ 99.0, 510.0, 32.0, 22.0 ],
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
                                    "id": "obj-24",
                                    "ignoreclick": 1,
                                    "legacytextcolor": 1,
                                    "maxclass": "textbutton",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "", "", "int" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 84.0, 416.5, 20.0, 20.0 ],
                                    "rounded": 60.0,
                                    "saved_attribute_attributes": {
                                        "bgcolor": {
                                            "expression": "themecolor.lesson_downloading"
                                        }
                                    },
                                    "text": "5",
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
                                    "patching_rect": [ 46.0, 235.0, 20.0, 20.0 ],
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
                                    "patching_rect": [ 170.0, 178.0, 20.0, 20.0 ],
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
                                    "patching_rect": [ 137.0, 512.0, 20.0, 20.0 ],
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
                                    "patching_rect": [ 219.0, 477.0, 20.0, 20.0 ],
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
                                    "destination": [ "obj-16", 0 ],
                                    "midpoints": [ 26.5, 654.8046875, 492.23828125, 654.8046875, 492.23828125, 427.796875, 629.5, 427.796875 ],
                                    "source": [ "obj-13", 0 ]
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
                                    "destination": [ "obj-21", 0 ],
                                    "source": [ "obj-7", 2 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-4", 0 ],
                                    "source": [ "obj-7", 1 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-48", 1 ],
                                    "order": 0,
                                    "source": [ "obj-7", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-48", 0 ],
                                    "order": 1,
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