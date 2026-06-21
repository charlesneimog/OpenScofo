{
    "patcher": {
        "fileversion": 1,
        "appversion": {
            "major": 9,
            "minor": 1,
            "revision": 3,
            "architecture": "x64",
            "modernui": 1
        },
        "classnamespace": "box",
        "rect": [ 331.0, 258.0, 666.0, 452.0 ],
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
                            "revision": 3,
                            "architecture": "x64",
                            "modernui": 1
                        },
                        "classnamespace": "box",
                        "rect": [ 331.0, 284.0, 666.0, 426.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 390.0, 341.0, 256.0, 20.0 ],
                                    "text": "OPEN OpenScofo Documentation Website"
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
                                    "patching_rect": [ 364.0, 341.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-29",
                                    "linecount": 3,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 364.0, 369.0, 300.0, 49.0 ],
                                    "presentation_linecount": 3,
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/descriptors/"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-17",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 0,
                                    "patching_rect": [ 32.0, 129.0, 35.0, 22.0 ],
                                    "text": "dac~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-13",
                                    "maxclass": "multislider",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 182.0, 305.0, 144.0, 113.0 ],
                                    "size": 4
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "id": "obj-26",
                                    "linecount": 8,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 352.0, 12.5, 303.0, 114.0 ],
                                    "text": "You can also use ONNX models. These are AI models trained with the OpenScofo Python module, and in PureData, you can use the object o.train. Unfortunately, this is not possible in Max due to the absence of the py4pd object.\n\nONNX models can be used to identify complex timbres and extended techniques."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-12",
                                    "linecount": 2,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 183.0, 197.0, 372.0, 35.0 ],
                                    "text": "set onnxmodel flute.onnx mfcc logmelspectrum zcr centroid flatness hfr"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-4",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "" ],
                                    "patching_rect": [ 182.0, 277.0, 65.0, 22.0 ],
                                    "text": "route onnx"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-14",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 82.0, 158.0, 45.0, 27.0 ],
                                    "text": "1)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-10",
                                    "maxclass": "toggle",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "int" ],
                                    "parameter_enable": 1,
                                    "patching_rect": [ 53.0, 159.0, 24.0, 24.0 ],
                                    "saved_attribute_attributes": {
                                        "valueof": {
                                            "parameter_enum": [ "off", "on" ],
                                            "parameter_longname": "toggle",
                                            "parameter_mmax": 1,
                                            "parameter_modmode": 0,
                                            "parameter_shortname": "toggle",
                                            "parameter_type": 2
                                        }
                                    },
                                    "valuepopuplabel": 1,
                                    "varname": "toggle"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-7",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 53.0, 197.0, 121.0, 22.0 ],
                                    "text": "set justdescription $1"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-20",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 53.0, 49.5, 45.0, 27.0 ],
                                    "text": "4)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-19",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 156.0, 12.5, 45.0, 27.0 ],
                                    "text": "3)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
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
                                    "text": "open bwv-1013.wav"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-2",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 2,
                                    "outlettype": [ "signal", "bang" ],
                                    "patching_rect": [ 10.0, 96.0, 47.0, 22.0 ],
                                    "text": "sfplay~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-3",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "float", "float", "" ],
                                    "patching_rect": [ 10.0, 251.0, 191.0, 22.0 ],
                                    "text": "openscofo~ onnx"
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
                                    "patching_rect": [ 27.0, 51.0, 24.0, 24.0 ]
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
                                    "destination": [ "obj-3", 0 ],
                                    "source": [ "obj-12", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-17", 1 ],
                                    "order": 0,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-17", 0 ],
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
                                    "destination": [ "obj-4", 0 ],
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
                                    "destination": [ "obj-13", 0 ],
                                    "source": [ "obj-4", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-3", 0 ],
                                    "source": [ "obj-7", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-2", 0 ],
                                    "source": [ "obj-8", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 12.0, 75.0, 74.0, 22.0 ],
                    "presentation": 1,
                    "presentation_rect": [ 144.0, 140.0, 100.0, 22.0 ],
                    "text": "p AI | ONNX"
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
                            "revision": 3,
                            "architecture": "x64",
                            "modernui": 1
                        },
                        "classnamespace": "box",
                        "rect": [ 0.0, 26.0, 666.0, 426.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 390.0, 341.0, 256.0, 20.0 ],
                                    "text": "OPEN OpenScofo Documentation Website"
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
                                    "patching_rect": [ 364.0, 341.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-29",
                                    "linecount": 3,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 364.0, 369.0, 300.0, 49.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/descriptors/"
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "id": "obj-26",
                                    "linecount": 3,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 362.0, 12.5, 302.0, 47.0 ],
                                    "text": "To check complet list of descriptors on OpenScofo Documentation. Some of the descriptors are: mfcc, crest, centroid, flatness, flux, etc..."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-25",
                                    "linecount": 3,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 236.0, 158.0, 210.0, 47.0 ],
                                    "text": "To use just the descriptors of OpenScofo, you need to set just description to 1."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-21",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 0,
                                    "patching_rect": [ 27.0, 158.0, 35.0, 22.0 ],
                                    "text": "dac~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-18",
                                    "maxclass": "multislider",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "" ],
                                    "parameter_enable": 0,
                                    "patching_rect": [ 104.0, 267.0, 178.0, 88.0 ],
                                    "size": 13
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-15",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 2,
                                    "outlettype": [ "", "" ],
                                    "patching_rect": [ 104.0, 238.0, 65.0, 22.0 ],
                                    "text": "route mfcc"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-14",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 142.0, 125.0, 45.0, 27.0 ],
                                    "text": "1)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
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
                                    "patching_rect": [ 113.0, 127.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-7",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 113.0, 158.0, 121.0, 22.0 ],
                                    "text": "set justdescription $1"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-20",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 53.0, 49.5, 45.0, 27.0 ],
                                    "text": "3)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-19",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 156.0, 12.5, 45.0, 27.0 ],
                                    "text": "2)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
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
                                    "text": "open bwv-1013.wav"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-2",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 2,
                                    "outlettype": [ "signal", "bang" ],
                                    "patching_rect": [ 10.0, 96.0, 47.0, 22.0 ],
                                    "text": "sfplay~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-3",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 3,
                                    "outlettype": [ "float", "float", "" ],
                                    "patching_rect": [ 10.0, 210.0, 113.0, 22.0 ],
                                    "text": "openscofo~ mfcc"
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
                                    "patching_rect": [ 27.0, 51.0, 24.0, 24.0 ]
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
                                    "destination": [ "obj-18", 0 ],
                                    "source": [ "obj-15", 0 ]
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
                                    "destination": [ "obj-15", 0 ],
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
                                    "destination": [ "obj-3", 0 ],
                                    "source": [ "obj-7", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-2", 0 ],
                                    "source": [ "obj-8", 0 ]
                                }
                            }
                        ]
                    },
                    "patching_rect": [ 12.0, 42.0, 79.0, 22.0 ],
                    "presentation": 1,
                    "presentation_rect": [ 129.0, 125.0, 100.0, 22.0 ],
                    "text": "p Descriptors"
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
                            "revision": 3,
                            "architecture": "x64",
                            "modernui": 1
                        },
                        "classnamespace": "box",
                        "rect": [ 0.0, 26.0, 666.0, 426.0 ],
                        "showontab": 1,
                        "boxes": [
                            {
                                "box": {
                                    "id": "obj-33",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 390.0, 341.0, 256.0, 20.0 ],
                                    "text": "OPEN OpenScofo Documentation Website"
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
                                    "patching_rect": [ 364.0, 341.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-29",
                                    "linecount": 3,
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 364.0, 369.0, 300.0, 49.0 ],
                                    "text": ";\rmax launchbrowser https://charlesneimog.github.io/OpenScofo/descriptors/"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-16",
                                    "linecount": 2,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 83.0, 49.5, 156.0, 33.0 ],
                                    "text": "Play the music (or use a record)."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-15",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 160.0, 152.0, 156.0, 20.0 ],
                                    "text": "Check a minimal example."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-14",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 156.0, 274.5, 120.0, 20.0 ],
                                    "text": "Start the inference."
                                }
                            },
                            {
                                "box": {
                                    "fontface": 1,
                                    "id": "obj-7",
                                    "linecount": 11,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 355.5, 12.5, 295.0, 154.0 ],
                                    "text": "The score file is a plain text file that must follow the structure of the musical score. To create or edit a score file, you can refer to the OpenScofo documentation for detailed syntax and guidelines.\n\n\nFor convenience, OpenScofo provides an online editor where you can create score files and import existing MusicXML files directly. This allows for seamless conversion from standard music notation to the OpenScofo text format."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-25",
                                    "linecount": 2,
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 241.0, 234.0, 143.0, 33.0 ],
                                    "text": "Load the score and start the inference process."
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-6",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 1,
                                    "outlettype": [ "signal" ],
                                    "patcher": {
                                        "fileversion": 1,
                                        "appversion": {
                                            "major": 9,
                                            "minor": 1,
                                            "revision": 3,
                                            "architecture": "x64",
                                            "modernui": 1
                                        },
                                        "classnamespace": "box",
                                        "rect": [ 48.0, 105.0, 962.0, 687.0 ],
                                        "boxes": [
                                            {
                                                "box": {
                                                    "fontname": "Lato",
                                                    "fontsize": 12.0,
                                                    "id": "obj-1",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 339.0, 437.5, 124.0, 23.0 ],
                                                    "text": "gen~ @gen freeverb"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-26",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 339.0, 374.0, 173.0, 22.0 ],
                                                    "text": "*~"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-13",
                                                    "maxclass": "message",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 493.0, 173.0, 45.0, 22.0 ],
                                                    "text": "$1 250"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-22",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 2,
                                                    "outlettype": [ "signal", "bang" ],
                                                    "patching_rect": [ 493.0, 205.0, 34.0, 22.0 ],
                                                    "text": "line~"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-9",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 2,
                                                    "outlettype": [ "bang", "float" ],
                                                    "patching_rect": [ 339.0, 107.0, 206.0, 22.0 ],
                                                    "text": "t b f"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-27",
                                                    "maxclass": "toggle",
                                                    "numinlets": 1,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "int" ],
                                                    "parameter_enable": 0,
                                                    "patching_rect": [ 195.5, 112.0, 24.0, 24.0 ]
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-24",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 306.0, 576.0, 29.5, 22.0 ],
                                                    "text": "+~"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-23",
                                                    "maxclass": "message",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 339.0, 148.0, 42.0, 22.0 ],
                                                    "text": "start 0"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-21",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 2,
                                                    "outlettype": [ "signal", "bang" ],
                                                    "patching_rect": [ 339.0, 175.0, 103.0, 22.0 ],
                                                    "text": "play~ somesound"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "fontname": "Arial",
                                                    "fontsize": 13.0,
                                                    "id": "obj-25",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 2,
                                                    "outlettype": [ "float", "bang" ],
                                                    "patching_rect": [ 179.0, 39.5, 167.0, 23.0 ],
                                                    "text": "buffer~ somesound 100000"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-20",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 338.91, 70.0, 38.0, 22.0 ],
                                                    "text": "r play"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-19",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 195.5, 70.0, 50.0, 22.0 ],
                                                    "text": "r record"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-12",
                                                    "maxclass": "newobj",
                                                    "numinlets": 3,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 179.0, 173.0, 115.0, 22.0 ],
                                                    "text": "record~ somesound"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "autowrite": 1,
                                                    "code": "/* Generated by OpenScofo online editor */\n\nBPM 80\n\n// Measure number 1\n\nNOTE E5 0.25\nNOTE A5 0.25\nNOTE G#5 0.25\nNOTE A5 0.25\nNOTE C6 0.25\nNOTE A5 0.25\nNOTE E5 0.25\n    sendto trans [1 7.02]\n\nNOTE A4 0.25\nNOTE E5 0.25\nNOTE A5 0.25\nNOTE G#5 0.25\nNOTE A5 0.25\nNOTE C6 0.25\nNOTE A5 0.25\nNOTE E5 0.25\n\n// Measure number 2\nNOTE A4 0.25\nNOTE C5 0.25\nNOTE E5 0.25\nNOTE F5 0.25\nNOTE G#4 0.25\nNOTE F5 0.25\nNOTE E5 0.25\nNOTE D5 0.25\nNOTE C5 0.25\nNOTE E5 0.25\nNOTE G#5 0.25\nNOTE A5 0.25\nNOTE E4 0.25\nNOTE D5 0.25\nNOTE C5 0.25\nNOTE B4 0.25\n\tsendto trans [0 0]\n\n// Measure number 3\nNOTE A4 0.25\nNOTE C5 0.25\nNOTE E5 0.25\nNOTE F5 0.25\nNOTE G#4 0.25\nNOTE F5 0.25\nNOTE E5 0.25\nNOTE D5 0.25\nNOTE C5 0.25\nNOTE E5 0.25\nNOTE G#5 0.25\nNOTE A5 0.25\nNOTE E4 0.25\nNOTE D5 0.25\nNOTE C5 0.25\nNOTE B4 0.25\n\n// Measure number 4\nNOTE A4 0.25\nNOTE C5 0.25\nNOTE G#4 0.25\nNOTE A4 0.25\nNOTE F4 0.25\nNOTE D5 0.25\nNOTE G#4 0.25\n    sendto record [1]\nNOTE A4 0.25\nNOTE E4 0.25\nNOTE C5 0.25\nNOTE G#4 0.25\nNOTE A4 0.25\nNOTE D4 0.25\nNOTE C5 0.25\nNOTE B4 0.25\n    sendto record [0]\nNOTE A4 0.25\n\n// Measure number 5\nNOTE G#4 0.25\nNOTE E4 0.25\nNOTE G#4 0.25\nNOTE B4 0.25\nNOTE E5 0.25\nNOTE G#5 0.25\n    sendto play [1]\nNOTE B5 0.25\nNOTE D6 0.25\nNOTE F5 0.25\nNOTE D6 0.25\nNOTE B5 0.25\nNOTE G#5 0.25\nNOTE E5 0.25\nNOTE D5 0.25\nNOTE C5 0.25\nNOTE B4 0.25\n\n// Measure number 6\nNOTE C5 0.25\nNOTE A4 0.25\nNOTE C5 0.25\nNOTE E5 0.25\nNOTE A5 0.25\nNOTE C6 0.25\nNOTE E6 0.25\n\tsendto play [0]\nNOTE F5 0.25\nNOTE G5 0.25\nNOTE B5 0.25\nNOTE G5 0.25\nNOTE E5 0.25\n\tsendto play [1]\nNOTE C5 0.25\nNOTE Bb4 0.25\nNOTE A4 0.25\nNOTE G4 0.25\n\n// Measure number 7\nNOTE A4 0.25\nNOTE F4 0.25\nNOTE A4 0.25\nNOTE C5 0.25\nNOTE D5 0.25\nNOTE F5 0.25\nNOTE A5 0.25\nNOTE C5 0.25\nNOTE B4 0.25\n\tsendto play [1]\nNOTE G4 0.25\nNOTE B4 0.25\nNOTE D5 0.25\nNOTE E5 0.25\nNOTE G5 0.25\nNOTE B5 0.25\nNOTE D5 0.25\n\n// Measure number 8\nNOTE C5 0.25\nNOTE A4 0.25\nNOTE C5 0.25\nNOTE E5 0.25\nNOTE F5 0.25\nNOTE A5 0.25\nNOTE C6 0.25\nNOTE E5 0.25\nNOTE D5 0.25\nNOTE B4 0.25\nNOTE D5 0.25\nNOTE F#5 0.25\nNOTE G5 0.25\nNOTE B5 0.25\nNOTE D6 0.25\nNOTE F5 0.25\n\n// Measure number 9\nNOTE E5 0.25\nNOTE G5 0.25\nNOTE C6 0.25\nNOTE B5 0.25\nNOTE C6 0.25\nNOTE E6 0.25\nNOTE C6 0.25\nNOTE G5 0.25\nNOTE C5 0.25\nNOTE G5 0.25\nNOTE C6 0.25\nNOTE B5 0.25\nNOTE C6 0.25\nNOTE E6 0.25\nNOTE C6 0.25\nNOTE G5 0.25\n\n// Measure number 10\nNOTE C5 0.25\nNOTE G5 0.25\nNOTE Bb5 0.25\nNOTE A5 0.25\nNOTE Bb5 0.25\nNOTE E6 0.25\nNOTE Bb5 0.25\nNOTE G5 0.25\nNOTE C5 0.25\nNOTE G5 0.25\nNOTE Bb5 0.25\nNOTE A5 0.25\nNOTE Bb5 0.25\nNOTE E6 0.25\nNOTE Bb5 0.25\nNOTE G5 0.25\n\n// Measure number 11\nNOTE C#5 0.25\nNOTE G5 0.25\nNOTE Bb5 0.25\nNOTE A5 0.25\nNOTE Bb5 0.25\nNOTE E6 0.25\nNOTE Bb5 0.25\nNOTE G5 0.25\nNOTE C#5 0.25\nNOTE G5 0.25\nNOTE E5 0.25\nNOTE C#5 0.25\nNOTE A4 0.25\nNOTE G4 0.25\nNOTE F4 0.25\nNOTE E4 0.25\n\n// Measure number 12\nNOTE D4 0.25\nNOTE G5 0.25\nNOTE F5 0.25\nNOTE E5 0.25\nNOTE F5 0.25\nNOTE A5 0.25\nNOTE F5 0.25\nNOTE D5 0.25\nNOTE B4 0.25\nNOTE G#4 0.25\nNOTE B4 0.25\nNOTE D5 0.25\nNOTE E5 0.25\nNOTE D5 0.25\nNOTE C5 0.25\nNOTE B4 0.25\n\n// Measure number 13\nNOTE C5 0.25\nNOTE E5 0.25\nNOTE A5 0.25\nNOTE G#5 0.25\nNOTE A5 0.25\nNOTE C6 0.25\nNOTE A5 0.25\nNOTE F#5 0.25\nNOTE D#5 0.25\nNOTE B4 0.25\nNOTE D#5 0.25\nNOTE F#5 0.25\nNOTE B5 0.25\nNOTE A5 0.25\nNOTE G5 0.25\nNOTE F#5 0.25\n\n// Measure number 14\nNOTE G5 0.25\nNOTE F#5 0.25\nNOTE E5 0.25\nNOTE G5 0.25\nNOTE A4 0.25\nNOTE E5 0.25\nNOTE G5 0.25\nNOTE C6 0.25\nNOTE F#5 0.25\nNOTE E5 0.25\nNOTE D5 0.25\nNOTE F#5 0.25\nNOTE G4 0.25\nNOTE D5 0.25\nNOTE F#5 0.25\nNOTE B5 0.25\n\n// Measure number 15\nNOTE E5 0.25\nNOTE D5 0.25\nNOTE C5 0.25\nNOTE E5 0.25\nNOTE F#4 0.25\nNOTE E5 0.25\nNOTE F#5 0.25\nNOTE A5 0.25\nNOTE D#5 0.25\nNOTE B4 0.25\nNOTE C5 0.25\nNOTE A4 0.25\nNOTE D#4 0.25\nNOTE C5 0.25\nNOTE B4 0.25\nNOTE A4 0.25\n\n// Measure number 16\nNOTE G4 0.25\nNOTE E5 0.25\nNOTE F5 0.25\nNOTE D5 0.25\nNOTE G#4 0.25\nNOTE F5 0.25\nNOTE E5 0.25\nNOTE D5 0.25\nNOTE A4 0.25\nNOTE F#5 0.25\nNOTE G5 0.25\nNOTE E5 0.25\nNOTE A#4 0.25\nNOTE G5 0.25\nNOTE F#5 0.25\nNOTE E5 0.25\n\n// Measure number 17\nNOTE D#5 0.25\nNOTE B5 0.25\nNOTE G#5 0.25\nNOTE D5 0.25\nNOTE C#5 0.25\nNOTE A5 0.25\nNOTE F#5 0.25\nNOTE C5 0.25\nNOTE B4 0.25\nNOTE G5 0.25\nNOTE E5 0.25\nNOTE Bb4 0.25\nNOTE A4 0.25\nNOTE F5 0.25\nNOTE D#5 0.25\nNOTE A4 0.25\n\n// Measure number 18\nNOTE G4 0.25\nNOTE E5 0.25\nNOTE C5 0.25\nNOTE A4 0.25\nNOTE F#4 0.25\nNOTE C5 0.25\nNOTE A4 0.25\nNOTE F#4 0.25\nNOTE D#4 0.25\nNOTE F#4 0.25\nNOTE A4 0.25\nNOTE C5 0.25\nNOTE B4 0.25\nNOTE A5 0.25\nNOTE G5 0.25\nNOTE F#5 0.25\n",
                                                    "filename": "bwv-1013.txt",
                                                    "fontface": 0,
                                                    "fontname": "<Monospaced>",
                                                    "fontsize": 12.0,
                                                    "id": "obj-18",
                                                    "maxclass": "text.codebox",
                                                    "numinlets": 2,
                                                    "numoutlets": 3,
                                                    "outlettype": [ "", "", "" ],
                                                    "patching_rect": [ 677.0, 118.0, 621.0, 551.0 ],
                                                    "saved_object_attributes": {
                                                        "parameter_enable": 0,
                                                        "parameter_mappable": 0
                                                    },
                                                    "textfile": {
                                                        "filename": "bwv-1013.txt",
                                                        "flags": 1,
                                                        "embed": 0,
                                                        "autowatch": 1
                                                    }
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-17",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 31.0, 521.0, 29.5, 22.0 ],
                                                    "text": "+~"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-16",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 53.0, 438.0, 29.5, 22.0 ],
                                                    "text": "*~"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-15",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 199.0, 308.0, 59.0, 22.0 ],
                                                    "text": "transratio"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-14",
                                                    "maxclass": "message",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 677.0, 88.0, 125.0, 22.0 ],
                                                    "text": "filename bwv-1013.txt"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-10",
                                                    "maxclass": "message",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 266.0, 302.0, 52.0, 22.0 ],
                                                    "text": "$1 2000"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-8",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 2,
                                                    "outlettype": [ "signal", "bang" ],
                                                    "patching_rect": [ 266.0, 334.0, 34.0, 22.0 ],
                                                    "text": "line~"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-7",
                                                    "maxclass": "newobj",
                                                    "numinlets": 1,
                                                    "numoutlets": 2,
                                                    "outlettype": [ "float", "float" ],
                                                    "patching_rect": [ 199.0, 266.0, 61.0, 22.0 ],
                                                    "text": "unpack f f"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "fontname": "Arial",
                                                    "fontsize": 13.0,
                                                    "id": "obj-6",
                                                    "maxclass": "newobj",
                                                    "numinlets": 2,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 53.0, 357.0, 165.0, 23.0 ],
                                                    "text": "pfft~ gizmo_loadme 2048 4"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "id": "obj-4",
                                                    "maxclass": "newobj",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "" ],
                                                    "patching_rect": [ 199.0, 237.0, 43.0, 22.0 ],
                                                    "text": "r trans"
                                                }
                                            },
                                            {
                                                "box": {
                                                    "comment": "",
                                                    "id": "obj-2",
                                                    "index": 1,
                                                    "maxclass": "outlet",
                                                    "numinlets": 1,
                                                    "numoutlets": 0,
                                                    "patching_rect": [ 306.0, 619.0, 30.0, 30.0 ]
                                                }
                                            },
                                            {
                                                "box": {
                                                    "comment": "",
                                                    "id": "obj-3",
                                                    "index": 1,
                                                    "maxclass": "inlet",
                                                    "numinlets": 0,
                                                    "numoutlets": 1,
                                                    "outlettype": [ "signal" ],
                                                    "patching_rect": [ 31.0, 40.0, 31.0, 31.0 ]
                                                }
                                            }
                                        ],
                                        "lines": [
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-24", 1 ],
                                                    "source": [ "obj-1", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-8", 0 ],
                                                    "source": [ "obj-10", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-22", 0 ],
                                                    "source": [ "obj-13", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-18", 0 ],
                                                    "source": [ "obj-14", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-6", 1 ],
                                                    "source": [ "obj-15", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-17", 1 ],
                                                    "midpoints": [ 62.5, 490.5, 51.0, 490.5 ],
                                                    "source": [ "obj-16", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-24", 0 ],
                                                    "source": [ "obj-17", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-27", 0 ],
                                                    "source": [ "obj-19", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-9", 0 ],
                                                    "source": [ "obj-20", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-26", 0 ],
                                                    "source": [ "obj-21", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-26", 1 ],
                                                    "source": [ "obj-22", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-21", 0 ],
                                                    "source": [ "obj-23", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-2", 0 ],
                                                    "source": [ "obj-24", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-1", 0 ],
                                                    "source": [ "obj-26", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-12", 0 ],
                                                    "source": [ "obj-27", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-12", 0 ],
                                                    "midpoints": [ 40.5, 97.5, 188.5, 97.5 ],
                                                    "order": 0,
                                                    "source": [ "obj-3", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-17", 0 ],
                                                    "order": 2,
                                                    "source": [ "obj-3", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-6", 0 ],
                                                    "midpoints": [ 40.5, 163.5, 62.5, 163.5 ],
                                                    "order": 1,
                                                    "source": [ "obj-3", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-7", 0 ],
                                                    "source": [ "obj-4", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-16", 0 ],
                                                    "source": [ "obj-6", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-10", 0 ],
                                                    "source": [ "obj-7", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-15", 0 ],
                                                    "source": [ "obj-7", 1 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-16", 1 ],
                                                    "midpoints": [ 275.5, 397.0, 73.0, 397.0 ],
                                                    "source": [ "obj-8", 0 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-13", 0 ],
                                                    "source": [ "obj-9", 1 ]
                                                }
                                            },
                                            {
                                                "patchline": {
                                                    "destination": [ "obj-23", 0 ],
                                                    "source": [ "obj-9", 0 ]
                                                }
                                            }
                                        ]
                                    },
                                    "patching_rect": [ 58.0, 151.0, 100.0, 22.0 ],
                                    "text": "p simple-iteration"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 12.0,
                                    "id": "obj-12",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 205.0, 364.0, 41.0, 20.0 ],
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
                                    "patching_rect": [ 109.0, 366.0, 74.0, 20.0 ],
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
                                    "patching_rect": [ 251.0, 358.0, 77.0, 31.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-4",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 0,
                                    "patching_rect": [ 58.0, 191.0, 35.0, 22.0 ],
                                    "text": "dac~"
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-20",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 53.0, 49.5, 28.0, 27.0 ],
                                    "text": "4)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-19",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 156.0, 12.5, 45.0, 27.0 ],
                                    "text": "3)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-18",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 126.0, 271.0, 28.0, 27.0 ],
                                    "text": "2)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
                                }
                            },
                            {
                                "box": {
                                    "fontsize": 18.0,
                                    "id": "obj-17",
                                    "maxclass": "comment",
                                    "numinlets": 1,
                                    "numoutlets": 0,
                                    "patching_rect": [ 212.0, 237.0, 27.0, 27.0 ],
                                    "text": "1)",
                                    "textcolor": [ 1.0, 0.0, 0.0, 1.0 ]
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
                                    "patching_rect": [ 10.0, 358.0, 95.0, 35.0 ]
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
                                    "text": "open bwv-1013.wav"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-2",
                                    "maxclass": "newobj",
                                    "numinlets": 2,
                                    "numoutlets": 2,
                                    "outlettype": [ "signal", "bang" ],
                                    "patching_rect": [ 10.0, 96.0, 47.0, 22.0 ],
                                    "text": "sfplay~"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-3",
                                    "maxclass": "newobj",
                                    "numinlets": 1,
                                    "numoutlets": 2,
                                    "outlettype": [ "float", "float" ],
                                    "patching_rect": [ 10.0, 314.0, 260.0, 22.0 ],
                                    "text": "openscofo~"
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
                                    "patching_rect": [ 27.0, 51.0, 24.0, 24.0 ]
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-9",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 72.0, 239.0, 136.0, 22.0 ],
                                    "text": "score bwv-1013.txt"
                                }
                            },
                            {
                                "box": {
                                    "id": "obj-10",
                                    "maxclass": "message",
                                    "numinlets": 2,
                                    "numoutlets": 1,
                                    "outlettype": [ "" ],
                                    "patching_rect": [ 92.0, 273.0, 32.0, 22.0 ],
                                    "text": "start"
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
                                    "destination": [ "obj-3", 0 ],
                                    "order": 1,
                                    "source": [ "obj-2", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-6", 0 ],
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
                                    "destination": [ "obj-4", 1 ],
                                    "order": 0,
                                    "source": [ "obj-6", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-4", 0 ],
                                    "order": 1,
                                    "source": [ "obj-6", 0 ]
                                }
                            },
                            {
                                "patchline": {
                                    "destination": [ "obj-2", 0 ],
                                    "source": [ "obj-8", 0 ]
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
                    "patching_rect": [ 12.0, 11.0, 79.0, 22.0 ],
                    "presentation": 1,
                    "presentation_rect": [ 114.0, 110.0, 100.0, 22.0 ],
                    "text": "p Basic"
                }
            }
        ],
        "lines": [],
        "parameters": {
            "obj-3::obj-10": [ "toggle", "toggle", 0 ],
            "parameterbanks": {
                "0": {
                    "index": 0,
                    "name": "",
                    "parameters": [ "-", "-", "-", "-", "-", "-", "-", "-" ],
                    "buttons": [ "-", "-", "-", "-", "-", "-", "-", "-" ]
                }
            },
            "inherited_shortname": 1
        },
        "autosave": 0,
        "oscreceiveudpport": 0
    }
}
