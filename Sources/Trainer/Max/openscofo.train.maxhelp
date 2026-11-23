{
    "patcher": {
        "fileversion": 1,
        "appversion": {
            "major": 9,
            "minor": 0,
            "revision": 0,
            "architecture": "x64",
            "modernui": 1
        },
        "classnamespace": "box",
        "rect": [ 120.0, 120.0, 900.0, 600.0 ],
        "boxes": [
            {
                "box": {
                    "id": "obj-title",
                    "maxclass": "comment",
                    "fontsize": 22.0,
                    "patching_rect": [ 30.0, 20.0, 500.0, 32.0 ],
                    "text": "openscofo.train"
                }
            },
            {
                "box": {
                    "id": "obj-description",
                    "maxclass": "comment",
                    "linecount": 2,
                    "patching_rect": [ 30.0, 55.0, 760.0, 38.0 ],
                    "text": "Build an OpenScofo classification model and cut a Max buffer~ into labeled training files. Relative paths are resolved from the current patcher folder."
                }
            },
            {
                "box": {
                    "id": "obj-replace",
                    "maxclass": "message",
                    "numinlets": 2,
                    "numoutlets": 1,
                    "outlettype": [ "" ],
                    "patching_rect": [ 30.0, 125.0, 54.0, 22.0 ],
                    "text": "replace"
                }
            },
            {
                "box": {
                    "id": "obj-buffer",
                    "maxclass": "newobj",
                    "numinlets": 1,
                    "numoutlets": 2,
                    "outlettype": [ "float", "bang" ],
                    "patching_rect": [ 105.0, 125.0, 164.0, 22.0 ],
                    "text": "buffer~ training_audio"
                }
            },
            {
                "box": {
                    "id": "obj-cut",
                    "maxclass": "message",
                    "numinlets": 2,
                    "numoutlets": 1,
                    "outlettype": [ "" ],
                    "patching_rect": [ 30.0, 180.0, 168.0, 22.0 ],
                    "text": "cut training_audio violin"
                }
            },
            {
                "box": {
                    "id": "obj-cut-comment",
                    "maxclass": "comment",
                    "patching_rect": [ 215.0, 180.0, 390.0, 22.0 ],
                    "text": "cut <buffer> <label> writes label/label_0000.wav ..."
                }
            },
            {
                "box": {
                    "id": "obj-folder",
                    "maxclass": "message",
                    "numinlets": 2,
                    "numoutlets": 1,
                    "outlettype": [ "" ],
                    "patching_rect": [ 30.0, 230.0, 153.0, 22.0 ],
                    "text": "folder training_dataset"
                }
            },
            {
                "box": {
                    "id": "obj-descriptors",
                    "maxclass": "message",
                    "numinlets": 2,
                    "numoutlets": 1,
                    "outlettype": [ "" ],
                    "patching_rect": [ 200.0, 230.0, 250.0, 22.0 ],
                    "text": "descriptors mfcc logmel centroid flux"
                }
            },
            {
                "box": {
                    "id": "obj-analyze",
                    "maxclass": "message",
                    "numinlets": 2,
                    "numoutlets": 1,
                    "outlettype": [ "" ],
                    "patching_rect": [ 30.0, 270.0, 55.0, 22.0 ],
                    "text": "analyze"
                }
            },
            {
                "box": {
                    "id": "obj-train",
                    "maxclass": "message",
                    "numinlets": 2,
                    "numoutlets": 1,
                    "outlettype": [ "" ],
                    "patching_rect": [ 100.0, 270.0, 41.0, 22.0 ],
                    "text": "train"
                }
            },
            {
                "box": {
                    "id": "obj-export",
                    "maxclass": "message",
                    "numinlets": 2,
                    "numoutlets": 1,
                    "outlettype": [ "" ],
                    "patching_rect": [ 155.0, 270.0, 116.0, 22.0 ],
                    "text": "export model.onnx"
                }
            },
            {
                "box": {
                    "id": "obj-cancel",
                    "maxclass": "message",
                    "numinlets": 2,
                    "numoutlets": 1,
                    "outlettype": [ "" ],
                    "patching_rect": [ 285.0, 270.0, 48.0, 22.0 ],
                    "text": "cancel"
                }
            },
            {
                "box": {
                    "id": "obj-trainer",
                    "maxclass": "newobj",
                    "numinlets": 1,
                    "numoutlets": 3,
                    "outlettype": [ "", "", "" ],
                    "patching_rect": [ 30.0, 330.0, 112.0, 22.0 ],
                    "text": "openscofo.train"
                }
            },
            {
                "box": {
                    "id": "obj-progress",
                    "maxclass": "newobj",
                    "numinlets": 1,
                    "numoutlets": 0,
                    "patching_rect": [ 30.0, 385.0, 119.0, 22.0 ],
                    "text": "print TRAIN-PROGRESS"
                }
            },
            {
                "box": {
                    "id": "obj-status",
                    "maxclass": "newobj",
                    "numinlets": 1,
                    "numoutlets": 0,
                    "patching_rect": [ 180.0, 385.0, 112.0, 22.0 ],
                    "text": "print TRAIN-STATUS"
                }
            },
            {
                "box": {
                    "id": "obj-result",
                    "maxclass": "newobj",
                    "numinlets": 1,
                    "numoutlets": 0,
                    "patching_rect": [ 325.0, 385.0, 119.0, 22.0 ],
                    "text": "print TRAIN-RESULT"
                }
            },
            {
                "box": {
                    "id": "obj-outlets",
                    "maxclass": "comment",
                    "patching_rect": [ 30.0, 435.0, 630.0, 22.0 ],
                    "text": "Outlets: progress | status and errors | evaluation metrics and cut segment count"
                }
            }
        ],
        "lines": [
            { "patchline": { "source": [ "obj-replace", 0 ], "destination": [ "obj-buffer", 0 ] } },
            { "patchline": { "source": [ "obj-cut", 0 ], "destination": [ "obj-trainer", 0 ] } },
            { "patchline": { "source": [ "obj-folder", 0 ], "destination": [ "obj-trainer", 0 ] } },
            { "patchline": { "source": [ "obj-descriptors", 0 ], "destination": [ "obj-trainer", 0 ] } },
            { "patchline": { "source": [ "obj-analyze", 0 ], "destination": [ "obj-trainer", 0 ] } },
            { "patchline": { "source": [ "obj-train", 0 ], "destination": [ "obj-trainer", 0 ] } },
            { "patchline": { "source": [ "obj-export", 0 ], "destination": [ "obj-trainer", 0 ] } },
            { "patchline": { "source": [ "obj-cancel", 0 ], "destination": [ "obj-trainer", 0 ] } },
            { "patchline": { "source": [ "obj-trainer", 0 ], "destination": [ "obj-progress", 0 ] } },
            { "patchline": { "source": [ "obj-trainer", 1 ], "destination": [ "obj-status", 0 ] } },
            { "patchline": { "source": [ "obj-trainer", 2 ], "destination": [ "obj-result", 0 ] } }
        ],
        "autosave": 0
    }
}
