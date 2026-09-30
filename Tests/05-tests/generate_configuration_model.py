"""Regenerate the tiny, deterministic ONNX test fixture (requires onnx).

No Python packages or model training are needed to build/run the C++ tests.
The checked-in model classifies RMS <= 0.1 as quiet, otherwise loud. Its second
input is centroid, so a descriptor-order regression changes the prediction.
"""

from pathlib import Path

import onnx
from onnx import TensorProto, helper

node = helper.make_node(
    "TreeEnsembleClassifier",
    ["features"],
    ["label", "probabilities"],
    domain="ai.onnx.ml",
    classlabels_strings=["quiet", "loud"],
    nodes_treeids=[0, 0, 0],
    nodes_nodeids=[0, 1, 2],
    nodes_featureids=[0, 0, 0],
    nodes_modes=["BRANCH_LEQ", "LEAF", "LEAF"],
    nodes_values=[0.1, 0.0, 0.0],
    nodes_truenodeids=[1, 0, 0],
    nodes_falsenodeids=[2, 0, 0],
    nodes_missing_value_tracks_true=[0, 0, 0],
    class_treeids=[0, 0, 0, 0],
    class_nodeids=[1, 1, 2, 2],
    class_ids=[0, 1, 0, 1],
    class_weights=[1.0, 0.0, 0.0, 1.0],
    post_transform="NONE",
)
graph = helper.make_graph(
    [node],
    "configuration-test",
    [helper.make_tensor_value_info("features", TensorProto.FLOAT, [1, 2])],
    [
        helper.make_tensor_value_info("label", TensorProto.STRING, [1]),
        helper.make_tensor_value_info("probabilities", TensorProto.FLOAT, [1, 2]),
    ],
)
model = helper.make_model(graph, opset_imports=[helper.make_opsetid("ai.onnx.ml", 1)], ir_version=8)
# Deliberately omit openscofo.descriptors: this fixture must use the score list.
helper.set_model_props(model, {
    "openscofo.sample_rate": "44100",
    "openscofo.fft_size": "4096",
    "openscofo.hop_size": "256",
    "openscofo.labels": '["quiet", "loud"]',
})
onnx.checker.check_model(model)
onnx.save(model, Path(__file__).resolve().parents[1] / "02-score" / "configuration.onnx")
