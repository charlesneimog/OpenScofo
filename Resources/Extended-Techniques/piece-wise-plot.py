#!/usr/bin/env python3
"""
Interactive Streamlit dashboard for OpenScofo benchmark results.

Run:
    streamlit run benchmark_dashboard.py
"""

from __future__ import annotations

import json
import os
from OpenScofo import OpenScofo

os.chdir(os.path.dirname(__file__))

scofo = OpenScofo(48000, 2048, 512)

ok = scofo.load_score("./audios/score-1.txt")
if not ok:
    raise Exception("Score load failed")



print(dir(scofo.get_configuration()))
