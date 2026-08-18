#!/usr/bin/env python3
"""Explain which score contexts are associated with OpenScofo failures.

The analysis deliberately distinguishes a failure *trigger* from a failure that
occurs after the follower is already lost.  A trigger is a missed event whose
preceding scored event was matched.  This avoids making a long alignment
cascade look like dozens of independent failures caused by note properties.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Callable


DEFAULT_RESULTS = Path(__file__).with_name("follower_validation.json")
Event = dict[str, Any]
Predicate = Callable[[Event], bool]


def as_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def pitch_to_midi(pitch: str | None) -> int | None:
    if not pitch:
        return None
    match = re.fullmatch(r"([A-G])([#b]?)(-?\d+)", pitch)
    if not match:
        return None
    pitch_classes = {"C": 0, "D": 2, "E": 4, "F": 5, "G": 7, "A": 9, "B": 11}
    accidental = {"#": 1, "b": -1}.get(match.group(2), 0)
    return 12 * (int(match.group(3)) + 1) + pitch_classes[match.group(1)] + accidental


def parse_score(score_path: Path) -> list[Event]:
    """Parse the NOTE/PTECH/REST subset used by the benchmark generator."""
    tokens: list[Event] = []
    measure: int | None = None

    for line_number, raw_line in enumerate(
        score_path.read_text(encoding="utf-8").splitlines(), 1
    ):
        line = raw_line.strip()
        measure_match = re.match(r"//\s*Measure\s+(\d+)", line, re.IGNORECASE)
        if measure_match:
            measure = int(measure_match.group(1))
        if not line or line.startswith("//"):
            continue

        fields = line.split()
        try:
            if fields[0] == "REST" and len(fields) >= 2:
                tokens.append(
                    {
                        "kind": "REST",
                        "duration": float(fields[1]),
                        "line": line_number,
                        "measure": measure,
                    }
                )
            elif fields[0] == "NOTE" and len(fields) >= 3:
                tokens.append(
                    {
                        "kind": "NOTE",
                        "technique": None,
                        "pitch": fields[1],
                        "duration": float(fields[2]),
                        "line": line_number,
                        "measure": measure,
                    }
                )
            elif fields[0] == "PTECH" and len(fields) >= 4:
                tokens.append(
                    {
                        "kind": "PTECH",
                        "technique": fields[1],
                        "pitch": fields[2],
                        "duration": float(fields[3]),
                        "line": line_number,
                        "measure": measure,
                    }
                )
        except ValueError:
            continue

    events: list[Event] = []
    previous_event: Event | None = None
    rest_before = 0.0

    for token in tokens:
        if token["kind"] == "REST":
            rest_before += token["duration"]
            continue

        event = dict(token)
        event["position"] = len(events) + 1
        event["previous_kind"] = previous_event["kind"] if previous_event else "START"
        event["after_rest"] = rest_before > 0.0
        event["rest_before"] = rest_before
        event["ioi_beats"] = (
            (previous_event["duration"] if previous_event else 0.0) + rest_before
        )

        previous_pitch = previous_event.get("pitch") if previous_event else None
        event["repeated_pitch"] = bool(
            previous_pitch and event.get("pitch") == previous_pitch
        )
        current_midi = pitch_to_midi(event.get("pitch"))
        previous_midi = pitch_to_midi(previous_pitch)
        event["absolute_interval"] = (
            abs(current_midi - previous_midi)
            if current_midi is not None and previous_midi is not None
            else None
        )

        events.append(event)
        previous_event = event
        rest_before = 0.0

    return events


def select_run(storage: dict[str, Any], selector: str) -> tuple[int, dict[str, Any]]:
    runs = storage.get("run_history", [])
    if not isinstance(runs, list) or not runs:
        raise ValueError("No run_history entries were found")

    if selector == "latest":
        return len(runs) - 1, runs[-1]

    try:
        index = int(selector)
    except ValueError:
        matches = [
            (index, run)
            for index, run in enumerate(runs)
            if str(run.get("timestamp", "")) == selector
        ]
        if not matches:
            raise ValueError(f"No run has timestamp {selector!r}")
        return matches[-1]

    if index < 0:
        index += len(runs)
    if index < 0 or index >= len(runs):
        raise ValueError(f"Run index {selector} is outside 0..{len(runs) - 1}")
    return index, runs[index]


def history_for_run(storage: dict[str, Any], run: dict[str, Any]) -> list[dict[str, Any]]:
    def matches(entry: dict[str, Any]) -> bool:
        return (
            str(entry.get("timestamp", "")) == str(run.get("timestamp", ""))
            and str(entry.get("implementation", ""))
            == str(run.get("implementation", ""))
            and str(entry.get("protocol", "")) == str(run.get("protocol", ""))
            and abs(
                as_float(entry.get("tolerance_ms"))
                - as_float(run.get("tolerance_ms"))
            )
            < 1e-9
        )

    entries = [
        entry
        for entry in storage.get("history", [])
        if isinstance(entry, dict) and matches(entry)
    ]
    by_files = {
        (str(entry.get("audio_file", "")), str(entry.get("score_file", ""))): entry
        for entry in entries
    }
    ordered = []
    for piece in run.get("pieces", []):
        key = (str(piece.get("audio_file", "")), str(piece.get("score_file", "")))
        if key in by_files:
            ordered.append(by_files[key])
    return ordered or entries


def load_events(
    entries: list[dict[str, Any]], results_path: Path
) -> tuple[list[Event], list[str]]:
    rows: list[Event] = []
    warnings: list[str] = []

    for entry in entries:
        score_name = str(entry.get("score_file", ""))
        score_path = Path(score_name).expanduser()
        if not score_path.is_absolute():
            score_path = results_path.parent / score_path
        if not score_path.exists():
            warnings.append(f"Score not found: {score_path}")
            continue

        metrics = entry.get("metrics", {})
        events = parse_score(score_path)
        reference_events = int(metrics.get("reference_events", 0))
        if len(events) != reference_events:
            warnings.append(
                f"{score_path.name}: parsed {len(events)} events, but the run records "
                f"{reference_events}; skipped (the score may have changed since this run)"
            )
            continue

        missing = {int(position) for position in metrics.get("missing_positions", [])}
        misaligned = {
            int(position) for position in metrics.get("misaligned_positions", [])
        }
        previous_miss = False
        for event in events:
            event["piece"] = score_path.name
            event["miss"] = event["position"] in missing
            event["misaligned"] = event["position"] in misaligned
            event["undetected"] = event["miss"] and not event["misaligned"]
            event["previous_miss"] = previous_miss
            rows.append(event)
            previous_miss = bool(event["miss"])

    return rows, warnings


def wilson_interval(successes: int, total: int) -> tuple[float, float]:
    """Return an approximate 95% Wilson binomial confidence interval."""
    if total == 0:
        return 0.0, 0.0
    z = 1.959963984540054
    proportion = successes / total
    denominator = 1.0 + z * z / total
    center = (proportion + z * z / (2.0 * total)) / denominator
    margin = (
        z
        * math.sqrt(
            proportion * (1.0 - proportion) / total
            + z * z / (4.0 * total * total)
        )
        / denominator
    )
    return 100.0 * (center - margin), 100.0 * (center + margin)


def pattern_definitions() -> list[tuple[str, Predicate]]:
    return [
        ("PTECH event", lambda event: event["kind"] == "PTECH"),
        ("key_click event", lambda event: event.get("technique") == "key_click"),
        ("tongue_ram event", lambda event: event.get("technique") == "tongue_ram"),
        ("event immediately after REST", lambda event: event["after_rest"]),
        (
            "PTECH immediately after REST",
            lambda event: event["kind"] == "PTECH" and event["after_rest"],
        ),
        (
            "key_click after NOTE + REST",
            lambda event: event.get("technique") == "key_click"
            and event["previous_kind"] == "NOTE"
            and event["after_rest"],
        ),
        (
            "tongue_ram after NOTE + REST",
            lambda event: event.get("technique") == "tongue_ram"
            and event["previous_kind"] == "NOTE"
            and event["after_rest"],
        ),
        ("after NOTE", lambda event: event["previous_kind"] == "NOTE"),
        ("after PTECH", lambda event: event["previous_kind"] == "PTECH"),
        ("repeated pitch", lambda event: event["repeated_pitch"]),
        (
            "small pitch interval (<= 2 semitones)",
            lambda event: event["absolute_interval"] is not None
            and event["absolute_interval"] <= 2,
        ),
        (
            "large pitch leap (>= octave)",
            lambda event: event["absolute_interval"] is not None
            and event["absolute_interval"] >= 12,
        ),
        ("short IOI (<= 0.5 beat)", lambda event: event["ioi_beats"] <= 0.5),
        ("long IOI (>= 1.5 beats)", lambda event: event["ioi_beats"] >= 1.5),
    ]


def analyze_patterns(events: list[Event], minimum_support: int) -> list[dict[str, Any]]:
    eligible = [
        event
        for event in events
        if event["previous_kind"] != "START" and not event["previous_miss"]
    ]
    patterns = []

    for name, predicate in pattern_definitions():
        selected = [event for event in eligible if predicate(event)]
        complement = [event for event in eligible if not predicate(event)]
        if len(selected) < minimum_support:
            continue

        misses = sum(bool(event["miss"]) for event in selected)
        complement_misses = sum(bool(event["miss"]) for event in complement)
        rate = misses / len(selected)
        complement_rate = complement_misses / len(complement) if complement else 0.0
        ci_low, ci_high = wilson_interval(misses, len(selected))
        patterns.append(
            {
                "pattern": name,
                "support": len(selected),
                "trigger_misses": misses,
                "trigger_miss_pct": 100.0 * rate,
                "ci95_low_pct": ci_low,
                "ci95_high_pct": ci_high,
                "complement_pct": 100.0 * complement_rate,
                "risk_ratio": rate / complement_rate if complement_rate else None,
            }
        )

    return sorted(
        patterns,
        key=lambda pattern: pattern["risk_ratio"]
        if pattern["risk_ratio"] is not None
        else math.inf,
        reverse=True,
    )


def analyze_cascades(events: list[Event]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    run_lengths: list[int] = []
    examples: list[dict[str, Any]] = []

    for piece in sorted({str(event["piece"]) for event in events}):
        piece_events = [event for event in events if event["piece"] == piece]
        index = 0
        while index < len(piece_events):
            if not piece_events[index]["miss"]:
                index += 1
                continue
            end = index + 1
            while end < len(piece_events) and piece_events[end]["miss"]:
                end += 1
            length = end - index
            run_lengths.append(length)
            if length >= 2:
                start = piece_events[index]
                examples.append(
                    {
                        "piece": piece,
                        "start_position": start["position"],
                        "end_position": piece_events[end - 1]["position"],
                        "length": length,
                        "starting_event": " ".join(
                            value
                            for value in (
                                start["kind"],
                                start.get("technique"),
                                start.get("pitch"),
                            )
                            if value
                        ),
                        "after_rest": bool(start["after_rest"]),
                        "measure": start.get("measure"),
                        "score_line": start["line"],
                    }
                )
            index = end

    misses = [event for event in events if event["miss"]]
    after_miss = [event for event in events if event["previous_miss"]]
    after_match = [
        event
        for event in events
        if event["previous_kind"] != "START" and not event["previous_miss"]
    ]
    misses_in_cascades = sum(length for length in run_lengths if length >= 2)
    summary = {
        "total_events": len(events),
        "total_misses": len(misses),
        "isolated_misses": sum(length == 1 for length in run_lengths),
        "misses_in_cascades": misses_in_cascades,
        "misses_in_cascades_pct": (
            100.0 * misses_in_cascades / len(misses) if misses else 0.0
        ),
        "miss_after_prior_miss_pct": (
            100.0 * sum(bool(event["miss"]) for event in after_miss) / len(after_miss)
            if after_miss
            else 0.0
        ),
        "miss_after_prior_match_pct": (
            100.0
            * sum(bool(event["miss"]) for event in after_match)
            / len(after_match)
            if after_match
            else 0.0
        ),
        "longest_cascade": max(run_lengths, default=0),
    }
    return summary, sorted(examples, key=lambda item: item["length"], reverse=True)


def analyze_mechanisms(events: list[Event]) -> list[dict[str, Any]]:
    rows = []
    for kind in ("NOTE", "PTECH"):
        selected = [event for event in events if event["kind"] == kind]
        misses = [event for event in selected if event["miss"]]
        rows.append(
            {
                "event_type": kind,
                "events": len(selected),
                "misses": len(misses),
                "miss_pct": 100.0 * len(misses) / len(selected) if selected else 0.0,
                "never_detected": sum(bool(event["undetected"]) for event in misses),
                "detected_outside_tolerance": sum(
                    bool(event["misaligned"]) for event in misses
                ),
            }
        )
    return rows


def build_analysis(
    storage: dict[str, Any], results_path: Path, selector: str, minimum_support: int
) -> dict[str, Any]:
    run_index, run = select_run(storage, selector)
    entries = history_for_run(storage, run)
    events, warnings = load_events(entries, results_path)
    if not events:
        raise ValueError("No score events could be associated with the selected run")
    cascade_summary, cascade_examples = analyze_cascades(events)
    return {
        "run": {
            "index": run_index,
            "timestamp": run.get("timestamp"),
            "implementation": run.get("implementation"),
            "protocol": run.get("protocol"),
            "tolerance_ms": run.get("tolerance_ms"),
            "files_analyzed": len({event["piece"] for event in events}),
        },
        "method": {
            "minimum_support": minimum_support,
            "trigger_definition": (
                "A missed event whose preceding scored event was matched"
            ),
            "score_context_source": "Current score files on disk",
        },
        "cascade_summary": cascade_summary,
        "patterns": analyze_patterns(events, minimum_support),
        "mechanisms": analyze_mechanisms(events),
        "longest_cascades": cascade_examples,
        "warnings": warnings,
    }


def format_text(analysis: dict[str, Any], top: int) -> str:
    run = analysis["run"]
    cascade = analysis["cascade_summary"]
    lines = [
        "OpenScofo contextual failure analysis",
        "=======================================",
        (
            f"Run {run['index']}: {run['timestamp']} | {run['implementation']} | "
            f"{run['protocol']} ({as_float(run['tolerance_ms']):.0f} ms)"
        ),
        f"Files analyzed: {run['files_analyzed']}",
        "",
        "Failure cascades",
        "----------------",
        f"Total events: {cascade['total_events']}",
        f"Total misses: {cascade['total_misses']}",
        (
            f"Misses in multi-event cascades: {cascade['misses_in_cascades']} "
            f"({cascade['misses_in_cascades_pct']:.1f}%)"
        ),
        f"Miss after a prior match: {cascade['miss_after_prior_match_pct']:.1f}%",
        f"Miss after a prior miss:  {cascade['miss_after_prior_miss_pct']:.1f}%",
        f"Longest cascade: {cascade['longest_cascade']} events",
        "",
        "Strongest failure triggers",
        "--------------------------",
    ]

    for pattern in analysis["patterns"][:top]:
        ratio = pattern["risk_ratio"]
        ratio_text = f"{ratio:.2f}x" if ratio is not None else "n/a"
        lines.append(
            f"- {pattern['pattern']}: {pattern['trigger_miss_pct']:.1f}% "
            f"({pattern['trigger_misses']}/{pattern['support']}), "
            f"baseline {pattern['complement_pct']:.1f}%, risk {ratio_text}, "
            f"95% CI {pattern['ci95_low_pct']:.1f}-{pattern['ci95_high_pct']:.1f}%"
        )

    neutral = [
        pattern
        for pattern in analysis["patterns"]
        if pattern["risk_ratio"] is not None and 0.8 <= pattern["risk_ratio"] <= 1.2
    ]
    if neutral:
        lines.extend(
            [
                "",
                "Patterns not meaningfully worse in this run",
                "-----------------------------------------",
            ]
        )
        for pattern in neutral:
            lines.append(
                f"- {pattern['pattern']}: {pattern['trigger_miss_pct']:.1f}% vs "
                f"{pattern['complement_pct']:.1f}% baseline "
                f"({pattern['risk_ratio']:.2f}x)"
            )

    lines.extend(["", "Miss mechanism by event type", "----------------------------"])
    for mechanism in analysis["mechanisms"]:
        lines.append(
            f"- {mechanism['event_type']}: {mechanism['misses']}/{mechanism['events']} "
            f"missed ({mechanism['miss_pct']:.1f}%); "
            f"never detected={mechanism['never_detected']}, "
            f"detected outside tolerance={mechanism['detected_outside_tolerance']}"
        )

    examples = analysis["longest_cascades"][: min(top, 10)]
    if examples:
        lines.extend(["", "Longest cascade examples", "------------------------"])
        for example in examples:
            location = f"measure {example['measure']}, line {example['score_line']}"
            lines.append(
                f"- {example['piece']} positions {example['start_position']}-"
                f"{example['end_position']}: {example['length']} misses; starts with "
                f"{example['starting_event']} ({location}, after REST="
                f"{'yes' if example['after_rest'] else 'no'})"
            )

    if analysis["warnings"]:
        lines.extend(["", "Warnings", "--------"])
        lines.extend(f"- {warning}" for warning in analysis["warnings"])

    lines.extend(
        [
            "",
            "Interpretation note: associations identify where to investigate; "
            "they do not by themselves prove causation.",
            "Historical runs are reliable only when their score files have not "
            "changed since the run.",
        ]
    )
    return "\n".join(lines) + "\n"


def format_csv(patterns: list[dict[str, Any]]) -> str:
    if not patterns:
        return ""
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=list(patterns[0]))
    writer.writeheader()
    writer.writerows(patterns)
    return output.getvalue()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results",
        type=Path,
        default=DEFAULT_RESULTS,
        help=f"Benchmark JSON (default: {DEFAULT_RESULTS})",
    )
    parser.add_argument(
        "--run",
        default="latest",
        help="Run index, exact timestamp, or 'latest' (default: latest)",
    )
    parser.add_argument(
        "--min-support",
        type=int,
        default=20,
        help="Minimum eligible events required for a pattern (default: 20)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Maximum ranked patterns and examples in text output (default: 10)",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json", "csv"),
        default="text",
        help="Output format; CSV contains the pattern table only",
    )
    parser.add_argument(
        "--output", type=Path, help="Write output to a file instead of stdout"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.min_support < 1:
        print("error: --min-support must be at least 1", file=sys.stderr)
        return 2
    if args.top < 1:
        print("error: --top must be at least 1", file=sys.stderr)
        return 2

    results_path = args.results.expanduser().resolve()
    try:
        storage = json.loads(results_path.read_text(encoding="utf-8"))
        analysis = build_analysis(
            storage, results_path, args.run, args.min_support
        )
    except (OSError, json.JSONDecodeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    if args.format == "json":
        rendered = json.dumps(analysis, indent=2) + "\n"
    elif args.format == "csv":
        rendered = format_csv(analysis["patterns"])
    else:
        rendered = format_text(analysis, args.top)

    if args.output:
        args.output.expanduser().write_text(rendered, encoding="utf-8")
    else:
        sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
