#!/usr/bin/env python3

import argparse
import csv
import json
import re
from pathlib import Path

import librosa
import matplotlib.pyplot as plt
import numpy as np

from OpenScofo import OpenScofo, Descriptors


def get_score_sample_rate(score_path: str) -> int:
    text = Path(score_path).read_text(encoding="utf-8")
    match = re.search(r"(?mi)^\s*SR\s+([0-9]+(?:\.[0-9]+)?)\b", text)
    if match:
        return int(float(match.group(1)))
    return 48000


def load_ground_truth(json_path: str):
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)
    events = sorted(data["events"], key=lambda x: x["timestamp_seconds"])
    return data, events


def expected_events_at_times(events, times):
    event_times = np.asarray([x["timestamp_seconds"] for x in events], dtype=np.float64)
    event_numbers = np.asarray([x["event"] for x in events], dtype=np.int32)
    indices = np.searchsorted(event_times, times, side="right") - 1
    result = np.zeros(len(times), dtype=np.int32)
    valid = indices >= 0
    result[valid] = event_numbers[indices[valid]]
    return result


def sparse_state_ticks(ax, state_info, maximum=20):
    count = len(state_info)
    if count == 0:
        return
    step = max(1, int(np.ceil(count / maximum)))
    positions = np.arange(0, count, step)
    labels = [str(state_info[i]["score_pos"]) for i in positions]
    ax.set_yticks(positions)
    ax.set_yticklabels(labels)


def add_ground_truth_lines(ax, events, end_time):
    for event in events:
        t = float(event["timestamp_seconds"])
        if 0.0 <= t <= end_time:
            ax.axvline(t, linestyle=":", linewidth=0.7, alpha=0.22)


def plot_heatmap(values, times, state_info, events, title, output, log=False, inspection_time=None):
    if log:
        image = np.log10(np.maximum(values, 1e-12)).T
        label = "log10(value)"
        vmin, vmax = -8.0, 0.0
    else:
        image = values.T
        label = "value"
        vmin, vmax = 0.0, 1.0

    fig, ax = plt.subplots(figsize=(16, 9))
    im = ax.imshow(
        image,
        origin="lower",
        aspect="auto",
        extent=[times[0], times[-1], -0.5, len(state_info) - 0.5],
        interpolation="nearest",
        vmin=vmin,
        vmax=vmax,
    )

    add_ground_truth_lines(ax, events, times[-1])

    if inspection_time is not None:
        ax.axvline(inspection_time, linestyle="--", linewidth=1.5)

    sparse_state_ticks(ax, state_info)
    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Score event")
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def plot_path(times, predicted, expected, inspection_time, output):
    fig, ax = plt.subplots(figsize=(16, 7))
    ax.step(times, expected, where="post", linewidth=2.0, label="Ground truth")
    ax.step(times, predicted, where="post", linewidth=1.2, label="OpenScofo")

    if inspection_time is not None:
        ax.axvline(inspection_time, linestyle="--", linewidth=1.5, label="Inspection time")

    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Score event")
    ax.set_title("OpenScofo path vs ground truth")
    ax.grid(True, alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def plot_onset(times, onset, events, inspection_time, output):
    fig, ax = plt.subplots(figsize=(16, 6))
    ax.plot(times, onset, linewidth=1.0, label="Onset")

    for event in events:
        t = float(event["timestamp_seconds"])
        if times[0] <= t <= times[-1]:
            ax.axvline(t, linestyle=":", linewidth=0.8, alpha=0.45)

    if inspection_time is not None:
        ax.axvline(inspection_time, linestyle="--", linewidth=1.5, label="Inspection time")

    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Onset value")
    ax.set_title("Onset descriptor")
    ax.grid(True, alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def plot_local_states(times, forward, state_info, events, inspection_time, previous_states, next_states, output):
    if inspection_time is None:
        return

    row = int(np.argmin(np.abs(times - inspection_time)))

    expected_event = 0
    for event in events:
        if event["timestamp_seconds"] <= inspection_time:
            expected_event = int(event["event"])
        else:
            break

    center = None
    if expected_event > 0:
        for col, info in enumerate(state_info):
            if info["score_pos"] == expected_event:
                center = col
                break

    if center is None:
        center = int(np.argmax(forward[row]))

    first = max(0, center - previous_states + 1)
    last = min(len(state_info), center + next_states + 1)

    fig, ax = plt.subplots(figsize=(16, 9))

    for col in range(first, last):
        info = state_info[col]
        ax.plot(
            times,
            forward[:, col],
            linewidth=1.2,
            label=f"Event {info['score_pos']} [state {info['internal_index']}]",
        )

    ax.axvline(inspection_time, linestyle="--", linewidth=1.5, label="Inspection time")

    for event in events:
        t = float(event["timestamp_seconds"])
        if times[0] <= t <= times[-1]:
            ax.axvline(t, linestyle=":", linewidth=0.6, alpha=0.2)

    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Forward probability")
    ax.set_ylim(0.0, 1.0)
    ax.set_title(f"Local states around {inspection_time:.3f} s")
    ax.grid(True, alpha=0.2)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    plt.close(fig)


def print_inspection(times, forward, exit_prob, best_obs, onset, silence, state_info, events, inspection_time, top_n):
    if inspection_time is None:
        return

    row = int(np.argmin(np.abs(times - inspection_time)))

    expected = 0
    for event in events:
        if event["timestamp_seconds"] <= inspection_time:
            expected = int(event["event"])
        else:
            break

    ranking = np.argsort(forward[row])[::-1]

    print()
    print("=" * 88)
    print("INSPECTION")
    print("=" * 88)
    print(f"Requested time : {inspection_time:.6f} s")
    print(f"Frame time     : {times[row]:.6f} s")
    print(f"Ground truth   : event {expected}")
    print(f"Onset          : {onset[row]:.8f}")
    print(f"SilenceProb    : {silence[row]:.8f}")
    print()
    print("Top states:")

    for col in ranking[:top_n]:
        info = state_info[int(col)]
        print(
            f"  event={info['score_pos']:4d} "
            f"state={info['internal_index']:4d} "
            f"F={forward[row, col]:.10f} "
            f"Exit={exit_prob[row, col]:.10f} "
            f"Obs={best_obs[row, col]:.10f}"
        )


def write_timeline_csv(output, times, predicted, expected, onset, silence):
    with open(output, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "time_seconds",
            "predicted_event",
            "expected_event",
            "onset",
            "silence_probability",
        ])

        for i in range(len(times)):
            writer.writerow([
                f"{times[i]:.9f}",
                int(predicted[i]),
                int(expected[i]),
                f"{onset[i]:.10f}",
                f"{silence[i]:.10f}",
            ])


def main():
    parser = argparse.ArgumentParser(
        description="Record the complete OpenScofo state probability history."
    )

    parser.add_argument("audio", help="Audio file")
    parser.add_argument("score", help="OpenScofo score")
    parser.add_argument("ground_truth", help="Ground-truth JSON")
    parser.add_argument("--time", type=float, default=None, help="Optional time in seconds to inspect")
    parser.add_argument("--states", type=int, default=10, help="Previous states in local inspection plot")
    parser.add_argument("--next", type=int, default=3, help="Following states in local inspection plot")
    parser.add_argument("--top", type=int, default=15, help="Number of states printed at inspection time")
    parser.add_argument(
        "--include-silence-states",
        action="store_true",
        help="Include internal inter-event silence states in matrices",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("openscofo-history"),
    )

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    _, events = load_ground_truth(args.ground_truth)
    score_sr = get_score_sample_rate(args.score)

    print(f"Loading audio at {score_sr} Hz...")

    audio, _ = librosa.load(args.audio, sr=score_sr, mono=True)
    audio = np.ascontiguousarray(audio, dtype=np.float32)

    follower = OpenScofo(float(score_sr), 2048.0, 512.0)

    if not follower.load_score(args.score):
        raise RuntimeError("Could not load score")

    # Make onset available even if the score itself has no explicit ONSET observation.
    follower.request_descriptor(Descriptors.ODSONSET)

    sr = int(follower.get_sr())
    hop = int(follower.get_hop_size())
    fft_size = int(follower.get_fft_size())

    all_states = follower.get_states()

    if args.include_silence_states:
        state_indices = list(range(len(all_states)))
    else:
        state_indices = [
            i
            for i, state in enumerate(all_states)
            if state.score_pos > 0 and not state.inter_event_silence
        ]

    if not state_indices:
        raise RuntimeError("No score states found")

    state_info = []
    for i in state_indices:
        state = all_states[i]
        state_info.append({
            "column": len(state_info),
            "internal_index": i,
            "score_pos": int(state.score_pos),
            "line": int(state.line),
            "inter_event_silence": bool(state.inter_event_silence),
        })

    frame_count = len(audio) // hop
    state_count = len(state_indices)

    print(f"SR             : {sr}")
    print(f"FFT            : {fft_size}")
    print(f"Hop            : {hop}")
    print(f"Audio duration : {len(audio) / sr:.3f} s")
    print(f"Frames         : {frame_count}")
    print(f"Recorded states: {state_count}")

    times = np.zeros(frame_count, dtype=np.float64)
    forward_history = np.zeros((frame_count, state_count), dtype=np.float32)
    exit_history = np.zeros((frame_count, state_count), dtype=np.float32)
    obs_history = np.zeros((frame_count, state_count), dtype=np.float32)
    predicted_event = np.zeros(frame_count, dtype=np.int32)
    onset_history = np.zeros(frame_count, dtype=np.float32)
    silence_history = np.zeros(frame_count, dtype=np.float32)

    # State.forward / exit_prob / best_obs are circular buffers.
    # Snapshot the current slot on every processed frame to preserve full history.
    # The write index is obtained BEFORE process_block(), because processing advances m_Tau.
    for frame_index in range(frame_count):
        start = frame_index * hop
        block = np.ascontiguousarray(audio[start:start + hop], dtype=np.float32)

        write_index = follower.get_current_buffer_index()

        ok = follower.process_block(block)
        if not ok:
            raise RuntimeError(f"OpenScofo failed at frame {frame_index}")

        times[frame_index] = (start + hop) / sr
        predicted_event[frame_index] = int(follower.get_current_score_position())

        desc = follower.get_description()
        onset_history[frame_index] = float(desc.onset)
        silence_history[frame_index] = float(desc.silence)

        states = follower.get_states()

        for col, state_index in enumerate(state_indices):
            state = states[state_index]

            if write_index < len(state.forward):
                forward_history[frame_index, col] = state.forward[write_index]

            if write_index < len(state.exit_prob):
                exit_history[frame_index, col] = state.exit_prob[write_index]

            if write_index < len(state.best_obs):
                obs_history[frame_index, col] = state.best_obs[write_index]

        if frame_index % 1000 == 0:
            print(
                f"\rProcessing {frame_index:6d}/{frame_count:6d}",
                end="",
                flush=True,
            )

    print(f"\rProcessing {frame_count:6d}/{frame_count:6d}")

    expected_event = expected_events_at_times(events, times)

    np.savez_compressed(
        args.output_dir / "state-history.npz",
        times=times,
        forward=forward_history,
        exit_prob=exit_history,
        best_obs=obs_history,
        predicted_event=predicted_event,
        expected_event=expected_event,
        onset=onset_history,
        silence=silence_history,
    )

    with open(args.output_dir / "state-map.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                "audio": args.audio,
                "score": args.score,
                "sample_rate": sr,
                "fft_size": fft_size,
                "hop_size": hop,
                "states": state_info,
            },
            f,
            indent=2,
        )

    write_timeline_csv(
        args.output_dir / "timeline.csv",
        times,
        predicted_event,
        expected_event,
        onset_history,
        silence_history,
    )

    plot_path(
        times,
        predicted_event,
        expected_event,
        args.time,
        args.output_dir / "state-path.png",
    )

    plot_onset(
        times,
        onset_history,
        events,
        args.time,
        args.output_dir / "onset-history.png",
    )

    plot_heatmap(
        forward_history,
        times,
        state_info,
        events,
        "Forward probability — complete history",
        args.output_dir / "forward-history.png",
        log=False,
        inspection_time=args.time,
    )

    plot_heatmap(
        forward_history,
        times,
        state_info,
        events,
        "Forward probability — complete history (log)",
        args.output_dir / "forward-history-log.png",
        log=True,
        inspection_time=args.time,
    )

    plot_heatmap(
        exit_history,
        times,
        state_info,
        events,
        "Exit probability — complete history (log)",
        args.output_dir / "exit-history-log.png",
        log=True,
        inspection_time=args.time,
    )

    plot_heatmap(
        obs_history,
        times,
        state_info,
        events,
        "Observation likelihood — complete history",
        args.output_dir / "observation-history.png",
        log=False,
        inspection_time=args.time,
    )

    plot_local_states(
        times,
        forward_history,
        state_info,
        events,
        args.time,
        args.states,
        args.next,
        args.output_dir / "local-forward.png",
    )

    print_inspection(
        times,
        forward_history,
        exit_history,
        obs_history,
        onset_history,
        silence_history,
        state_info,
        events,
        args.time,
        args.top,
    )

    print()
    print(f"Complete history written to: {args.output_dir}")
    print("Files:")
    print("  state-history.npz")
    print("  state-map.json")
    print("  timeline.csv")
    print("  state-path.png")
    print("  forward-history.png")
    print("  forward-history-log.png")
    print("  exit-history-log.png")
    print("  observation-history.png")
    print("  onset-history.png")
    if args.time is not None:
        print("  local-forward.png")


if __name__ == "__main__":
    main()
