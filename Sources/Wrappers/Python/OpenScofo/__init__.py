"""
Copyright (c) 2024-2026 Charles K. Neimog
Website: charlesneimog.github.io

This file is part of a project licensed under the
GNU General Public License v3.0 or later (GPL-3.0-or-later).
See the LICENSE file for details.
"""

from ._OpenScofo import (
    OpenScofo,
    EventType,
    HMMType,
    AudioState,
    MarkovState,
    Description,
    Descriptors,
    Configuration,
)

__all__ = [
    "OpenScofo",
    "EventType",
    "HMMType",
    "AudioState",
    "MarkovState",
    "Description",
    "Descriptors",
    "Configuration",
    "ExtendedTechniqueClassifier",
]


def __getattr__(name):
    if name == "ExtendedTechniqueClassifier":
        from .train import ExtendedTechniqueClassifier

        return ExtendedTechniqueClassifier
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")