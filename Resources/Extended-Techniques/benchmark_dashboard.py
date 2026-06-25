#!/usr/bin/env python3
"""
Interactive Streamlit dashboard for OpenScofo benchmark results.

Run:
    streamlit run benchmark_dashboard.py
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import altair as alt
import pandas as pd
import streamlit as st

DEFAULT_RESULTS_PATH = Path(__file__).with_name("follower_validation.json")


def _as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _as_int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


@st.cache_data(show_spinner=False)
def load_results(path_text: str, _mtime: float) -> Dict[str, Any]:
    """Load and validate benchmark JSON data."""
    results_path = Path(path_text).expanduser()
    with open(results_path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError("results JSON must contain a top-level object")
    return payload


def build_run_dataframe(storage: Dict[str, Any]) -> Tuple[Dict[int, Dict[str, Any]], pd.DataFrame]:
    """Create a DataFrame view over run_history with normalized fields."""
    run_history = storage.get("run_history", [])
    if not isinstance(run_history, list):
        return {}, pd.DataFrame()

    run_map: Dict[int, Dict[str, Any]] = {}
    rows: List[Dict[str, Any]] = []

    for run_id, run in enumerate(run_history):
        if not isinstance(run, dict):
            continue

        metrics = run.get("global_metrics", {})
        if not isinstance(metrics, dict):
            metrics = {}

        timestamp = str(run.get("timestamp", ""))
        implementation = str(run.get("implementation", "unknown"))
        protocol = str(run.get("protocol", "unknown"))
        tolerance_ms = _as_float(run.get("tolerance_ms"), 0.0)
        files_processed = _as_int(
            run.get("files_processed"),
            _as_int(metrics.get("files_processed"), 0),
        )

        overall_precision_pct = _as_float(metrics.get("overall_precision_pct"), 0.0)
        piecewise_precision_pct = _as_float(metrics.get("piecewise_precision_pct"), 0.0)
        mean_abs_offset_ms = _as_float(metrics.get("global_mean_abs_offset_ms"), 0.0)
        total_missed_notes = _as_int(metrics.get("total_missed_notes"), 0)
        total_false_positives = _as_int(metrics.get("total_false_positives"), 0)

        run_map[run_id] = run
        rows.append(
            {
                "run_id": run_id,
                "timestamp": timestamp,
                "timestamp_dt": pd.to_datetime(timestamp, errors="coerce"),
                "implementation": implementation,
                "protocol": protocol,
                "tolerance_ms": tolerance_ms,
                "files_processed": files_processed,
                "overall_precision_pct": overall_precision_pct,
                "piecewise_precision_pct": piecewise_precision_pct,
                "global_mean_abs_offset_ms": mean_abs_offset_ms,
                "total_missed_notes": total_missed_notes,
                "total_false_positives": total_false_positives,
                "label": (
                    f"{timestamp} | {implementation} | {protocol} | "
                    f"P={overall_precision_pct:.2f}% | files={files_processed}"
                ),
            }
        )

    if not rows:
        return run_map, pd.DataFrame()

    runs_df = pd.DataFrame(rows).sort_values("timestamp", ascending=True).reset_index(
        drop=True
    )
    return run_map, runs_df


def _entry_match_for_run(entry: Dict[str, Any], run: Dict[str, Any]) -> bool:
    """Match one history entry to one run record."""
    return (
        str(entry.get("timestamp", "")) == str(run.get("timestamp", ""))
        and str(entry.get("implementation", ""))
        == str(run.get("implementation", ""))
        and str(entry.get("protocol", "")) == str(run.get("protocol", ""))
        and abs(
            _as_float(entry.get("tolerance_ms"), 0.0)
            - _as_float(run.get("tolerance_ms"), 0.0)
        )
        < 1e-9
    )


def _extract_piece_order(audio_file: str, fallback_order: int) -> int:
    stem = Path(audio_file).stem
    match = re.search(r"(\d+)$", stem)
    if match:
        return int(match.group(1))
    return fallback_order


def build_piece_dataframe(storage: Dict[str, Any], run: Dict[str, Any]) -> pd.DataFrame:
    """Create per-piece DataFrame for a selected run."""
    history = storage.get("history", [])
    if not isinstance(history, list):
        return pd.DataFrame()

    entries = [
        item
        for item in history
        if isinstance(item, dict) and _entry_match_for_run(item, run)
    ]

    if not entries:
        files_processed = _as_int(run.get("files_processed"), 0)
        entries = [item for item in history if isinstance(item, dict)]
        if files_processed > 0:
            entries = entries[-files_processed:]

    ordered_entries: List[Dict[str, Any]] = []
    pieces = run.get("pieces", [])

    if isinstance(pieces, list) and pieces:
        by_key: Dict[Tuple[str, str], Dict[str, Any]] = {}
        for entry in entries:
            key = (str(entry.get("audio_file", "")), str(entry.get("score_file", "")))
            by_key[key] = entry

        used_keys = set()
        for piece in pieces:
            if not isinstance(piece, dict):
                continue
            key = (str(piece.get("audio_file", "")), str(piece.get("score_file", "")))
            if key in by_key:
                ordered_entries.append(by_key[key])
                used_keys.add(key)

        for entry in entries:
            key = (str(entry.get("audio_file", "")), str(entry.get("score_file", "")))
            if key not in used_keys:
                ordered_entries.append(entry)
    else:
        ordered_entries = entries

    rows: List[Dict[str, Any]] = []
    for index, entry in enumerate(ordered_entries, start=1):
        metrics = entry.get("metrics", {})
        if not isinstance(metrics, dict):
            metrics = {}

        audio_file = str(entry.get("audio_file", ""))
        score_file = str(entry.get("score_file", ""))
        piece_name = Path(audio_file).stem if audio_file else f"piece_{index:02d}"

        missing_positions = metrics.get("missing_positions", [])
        misaligned_positions = metrics.get("misaligned_positions", [])
        unexpected_positions = metrics.get("unexpected_positions", [])

        if not isinstance(missing_positions, list):
            missing_positions = []
        if not isinstance(misaligned_positions, list):
            misaligned_positions = []
        if not isinstance(unexpected_positions, list):
            unexpected_positions = []

        rows.append(
            {
                "piece": piece_name,
                "piece_order": _extract_piece_order(audio_file, 10000 + index),
                "audio_file": audio_file,
                "score_file": score_file,
                "precision_pct": _as_float(metrics.get("precision_pct"), 0.0),
                "mean_abs_offset_ms": _as_float(metrics.get("mean_abs_offset_ms"), 0.0),
                "mean_offset_ms": _as_float(metrics.get("mean_offset_ms"), 0.0),
                "std_offset_ms": _as_float(metrics.get("std_offset_ms"), 0.0),
                "missed_notes": _as_int(metrics.get("missed_notes"), 0),
                "false_positives": _as_int(metrics.get("false_positives"), 0),
                "reference_events": _as_int(metrics.get("reference_events"), 0),
                "matched_events": _as_int(metrics.get("matched_events"), 0),
                "missed_notes_pct": _as_float(metrics.get("missed_notes_pct"), 0.0),
                "unexpected_event_tags": _as_int(metrics.get("unexpected_event_tags"), 0),
                "missing_positions": missing_positions,
                "misaligned_positions": misaligned_positions,
                "unexpected_positions": unexpected_positions,
            }
        )

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows).sort_values(["piece_order", "piece"], ascending=True)


def select_best_run(runs_df: pd.DataFrame, selected: pd.Series) -> Tuple[pd.Series, str]:
    """Pick best run in comparable scope using benchmark ranking priorities."""
    same_scope = runs_df[
        (runs_df["protocol"] == selected["protocol"])
        & ((runs_df["tolerance_ms"] - selected["tolerance_ms"]).abs() < 1e-9)
    ]

    if same_scope.empty:
        candidates = runs_df
        scope_label = "all runs"
    else:
        candidates = same_scope
        scope_label = "same protocol and tolerance"

    ranked = candidates.sort_values(
        by=[
            "overall_precision_pct",
            "global_mean_abs_offset_ms",
            "total_missed_notes",
            "files_processed",
            "timestamp",
        ],
        ascending=[False, True, True, False, False],
    )
    return ranked.iloc[0], scope_label


def render_metric_card(title: str, value: str, delta: str | None = None) -> None:
    if delta is None:
        st.metric(title, value)
    else:
        st.metric(title, value, delta=delta)


def main() -> None:
    st.set_page_config(page_title="OpenScofo Benchmark Dashboard", layout="wide")
    st.title("OpenScofo Benchmark Dashboard")
    st.caption("Interactive explorer for run history and per-piece alignment quality.")

    with st.sidebar:
        st.header("Data Source")
        results_path_text = st.text_input("Results JSON", str(DEFAULT_RESULTS_PATH))
        if st.button("Reload JSON"):
            st.cache_data.clear()
            st.rerun()

    results_path = Path(results_path_text).expanduser()
    if not results_path.exists():
        st.error(f"Results file not found: {results_path}")
        st.stop()

    try:
        storage = load_results(str(results_path), results_path.stat().st_mtime)
    except Exception as exc:
        st.error(f"Failed to read results JSON: {exc}")
        st.stop()

    run_map, runs_df = build_run_dataframe(storage)
    if runs_df.empty:
        st.warning("No run history found in this JSON file.")
        st.stop()

    with st.sidebar:
        st.header("Run Filters")
        implementation_options = sorted(runs_df["implementation"].dropna().unique().tolist())
        protocol_options = sorted(runs_df["protocol"].dropna().unique().tolist())

        selected_implementations = st.multiselect(
            "Implementation",
            options=implementation_options,
            default=implementation_options,
        )
        selected_protocols = st.multiselect(
            "Protocol",
            options=protocol_options,
            default=protocol_options,
        )

        files_min = int(runs_df["files_processed"].min())
        files_max = int(runs_df["files_processed"].max())
        if files_min == files_max:
            selected_files_range = (files_min, files_max)
            st.caption(f"Files processed is fixed at {files_min} in this dataset.")
        else:
            selected_files_range = st.slider(
                "Files processed",
                min_value=files_min,
                max_value=files_max,
                value=(files_min, files_max),
            )

    filtered_runs_df = runs_df[
        runs_df["implementation"].isin(selected_implementations)
        & runs_df["protocol"].isin(selected_protocols)
        & runs_df["files_processed"].between(
            selected_files_range[0],
            selected_files_range[1],
        )
    ]

    if filtered_runs_df.empty:
        st.warning("No runs match the active filters.")
        st.stop()

    run_labels = {
        int(row.run_id): str(row.label)
        for row in filtered_runs_df[["run_id", "label"]].itertuples(index=False)
    }
    run_id_options = filtered_runs_df["run_id"].tolist()

    with st.sidebar:
        st.header("Run Selection")
        selected_run_id = st.selectbox(
            "Choose run",
            options=run_id_options,
            index=len(run_id_options) - 1,
            format_func=lambda run_id: run_labels.get(run_id, str(run_id)),
        )

    selected_run = run_map[int(selected_run_id)]
    selected_row = runs_df.loc[runs_df["run_id"] == int(selected_run_id)].iloc[0]
    global_metrics = selected_run.get("global_metrics", {})
    if not isinstance(global_metrics, dict):
        global_metrics = {}

    best_row, best_scope = select_best_run(runs_df, selected_row)

    st.subheader(
        f"Run {selected_row['timestamp']} | {selected_row['implementation']} | {selected_row['protocol']}"
    )
    st.caption(f"Comparison scope for best run: {best_scope}")

    overall_delta = selected_row["overall_precision_pct"] - best_row["overall_precision_pct"]
    offset_delta = (
        selected_row["global_mean_abs_offset_ms"] - best_row["global_mean_abs_offset_ms"]
    )
    missed_delta = selected_row["total_missed_notes"] - best_row["total_missed_notes"]

    metric_cols = st.columns(6)
    with metric_cols[0]:
        render_metric_card(
            "Overall Precision",
            f"{_as_float(global_metrics.get('overall_precision_pct')):.2f}%",
            f"{overall_delta:+.2f}% vs best",
        )
    with metric_cols[1]:
        render_metric_card(
            "Piecewise Precision",
            f"{_as_float(global_metrics.get('piecewise_precision_pct')):.2f}%",
        )
    with metric_cols[2]:
        render_metric_card(
            "Global Mean |Offset|",
            f"{_as_float(global_metrics.get('global_mean_abs_offset_ms')):.2f} ms",
            f"{offset_delta:+.2f} ms vs best",
        )
    with metric_cols[3]:
        render_metric_card(
            "Missed Notes",
            str(_as_int(global_metrics.get("total_missed_notes"))),
            f"{missed_delta:+d} vs best",
        )
    with metric_cols[4]:
        render_metric_card(
            "False Positives",
            str(_as_int(global_metrics.get("total_false_positives"))),
        )
    with metric_cols[5]:
        render_metric_card(
            "Files Processed",
            str(_as_int(global_metrics.get("files_processed"))),
        )

    piece_df = build_piece_dataframe(storage, selected_run)
    if piece_df.empty:
        st.warning("No piece-level history entries found for the selected run.")
        st.stop()

    st.markdown("### Piece Quality Controls")
    quality_cols = st.columns(2)
    with quality_cols[0]:
        precision_threshold = st.slider(
            "Precision alert threshold (%)",
            min_value=0.0,
            max_value=100.0,
            value=70.0,
            step=0.5,
        )
    with quality_cols[1]:
        offset_threshold = st.slider(
            "Mean |offset| alert threshold (ms)",
            min_value=0.0,
            max_value=500.0,
            value=100.0,
            step=1.0,
        )

    piece_df = piece_df.copy()
    piece_df["quality_flag"] = piece_df.apply(
        lambda row: "attention"
        if row["precision_pct"] < precision_threshold
        or row["mean_abs_offset_ms"] > offset_threshold
        else "ok",
        axis=1,
    )

    attention_df = piece_df[piece_df["quality_flag"] == "attention"].copy()

    chart_cols = st.columns(2)

    with chart_cols[0]:
        st.markdown("### Precision by Piece")
        precision_chart = (
            alt.Chart(piece_df)
            .mark_bar()
            .encode(
                x=alt.X("piece:N", sort=piece_df["piece"].tolist(), title="Piece"),
                y=alt.Y("precision_pct:Q", title="Precision (%)", scale=alt.Scale(domain=[0, 100])),
                color=alt.condition(
                    alt.datum.quality_flag == "attention",
                    alt.value("#d94841"),
                    alt.value("#2b8cbe"),
                ),
                tooltip=[
                    alt.Tooltip("piece:N", title="Piece"),
                    alt.Tooltip("precision_pct:Q", title="Precision (%)", format=".2f"),
                    alt.Tooltip("mean_abs_offset_ms:Q", title="Mean |Offset| (ms)", format=".2f"),
                    alt.Tooltip("missed_notes:Q", title="Missed Notes"),
                    alt.Tooltip("false_positives:Q", title="False Positives"),
                ],
            )
            .properties(height=320)
        )
        st.altair_chart(precision_chart, width="stretch")

    with chart_cols[1]:
        st.markdown("### Offset vs Precision")
        scatter_chart = (
            alt.Chart(piece_df)
            .mark_circle(opacity=0.85)
            .encode(
                x=alt.X("mean_abs_offset_ms:Q", title="Mean |Offset| (ms)"),
                y=alt.Y("precision_pct:Q", title="Precision (%)", scale=alt.Scale(domain=[0, 100])),
                size=alt.Size("missed_notes:Q", title="Missed Notes", scale=alt.Scale(range=[60, 800])),
                color=alt.condition(
                    alt.datum.quality_flag == "attention",
                    alt.value("#d94841"),
                    alt.value("#2b8cbe"),
                ),
                tooltip=[
                    alt.Tooltip("piece:N", title="Piece"),
                    alt.Tooltip("mean_abs_offset_ms:Q", title="Mean |Offset| (ms)", format=".2f"),
                    alt.Tooltip("precision_pct:Q", title="Precision (%)", format=".2f"),
                    alt.Tooltip("missed_notes:Q", title="Missed Notes"),
                ],
            )
            .properties(height=320)
        )
        st.altair_chart(scatter_chart, width="stretch")

    st.markdown("### Attention Pieces")
    if attention_df.empty:
        st.success("All pieces are above the selected thresholds.")
    else:
        st.dataframe(
            attention_df[
                [
                    "piece",
                    "precision_pct",
                    "mean_abs_offset_ms",
                    "missed_notes",
                    "false_positives",
                    "missed_notes_pct",
                ]
            ]
            .rename(
                columns={
                    "piece": "Piece",
                    "precision_pct": "Precision (%)",
                    "mean_abs_offset_ms": "Mean |Offset| (ms)",
                    "missed_notes": "Missed Notes",
                    "false_positives": "False Positives",
                    "missed_notes_pct": "Missed (%)",
                }
            )
            .sort_values(["Precision (%)", "Mean |Offset| (ms)"], ascending=[True, False]),
            width="stretch",
            hide_index=True,
        )

    st.markdown("### Piece Metrics Table")
    sort_key = st.selectbox(
        "Sort table by",
        options=[
            "piece_order",
            "precision_pct",
            "mean_abs_offset_ms",
            "missed_notes",
            "false_positives",
        ],
        format_func=lambda value: {
            "piece_order": "Piece",
            "precision_pct": "Precision (%)",
            "mean_abs_offset_ms": "Mean |Offset| (ms)",
            "missed_notes": "Missed Notes",
            "false_positives": "False Positives",
        }.get(value, value),
    )

    ascending = sort_key == "piece_order"
    display_df = piece_df.sort_values(sort_key, ascending=ascending).copy()

    display_df = display_df[
        [
            "piece",
            "audio_file",
            "score_file",
            "precision_pct",
            "mean_abs_offset_ms",
            "mean_offset_ms",
            "std_offset_ms",
            "missed_notes",
            "false_positives",
            "reference_events",
            "matched_events",
            "missed_notes_pct",
            "unexpected_event_tags",
        ]
    ].rename(
        columns={
            "piece": "Piece",
            "audio_file": "Audio",
            "score_file": "Score",
            "precision_pct": "Precision (%)",
            "mean_abs_offset_ms": "Mean |Offset| (ms)",
            "mean_offset_ms": "Mean Offset (ms)",
            "std_offset_ms": "Offset Std (ms)",
            "missed_notes": "Missed Notes",
            "false_positives": "False Positives",
            "reference_events": "Reference Events",
            "matched_events": "Matched Events",
            "missed_notes_pct": "Missed (%)",
            "unexpected_event_tags": "Unexpected Tags",
        }
    )

    st.dataframe(display_df, width="stretch", hide_index=True)

    st.download_button(
        label="Download piece table as CSV",
        data=display_df.to_csv(index=False),
        file_name="openscofo_piece_metrics.csv",
        mime="text/csv",
    )

    st.markdown("### Run History Trends")
    trend_df = runs_df.copy().sort_values("timestamp")
    trend_df["run_label"] = trend_df["timestamp"].fillna("") + " | " + trend_df[
        "implementation"
    ].fillna("")

    trend_chart = (
        alt.Chart(trend_df)
        .transform_fold(
            ["overall_precision_pct", "global_mean_abs_offset_ms", "total_missed_notes"],
            as_=["metric", "value"],
        )
        .mark_line(point=True)
        .encode(
            x=alt.X("timestamp:T", title="Run timestamp"),
            y=alt.Y("value:Q", title="Metric value"),
            color=alt.Color("metric:N", title="Metric"),
            tooltip=[
                alt.Tooltip("timestamp:N", title="Timestamp"),
                alt.Tooltip("implementation:N", title="Implementation"),
                alt.Tooltip("protocol:N", title="Protocol"),
                alt.Tooltip("metric:N", title="Metric"),
                alt.Tooltip("value:Q", title="Value", format=".3f"),
            ],
        )
        .properties(height=320)
    )

    st.altair_chart(trend_chart, width="stretch")


if __name__ == "__main__":
    main()
