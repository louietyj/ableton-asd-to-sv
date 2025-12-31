#!/usr/bin/env python3
"""
Ableton ASD to Sonic Visualizer Beat CSV Converter

Extracts warp marker data from Ableton .asd files and generates a CSV file
with beat positions that can be imported into Sonic Visualizer as time instants.
"""

import argparse
import csv
import sys
from pathlib import Path

from abletonparsing import Clip  # type: ignore[import-untyped]


def parse_time_signature(ts_string: str) -> tuple[int, int]:
    """Parse time signature string like '4/4' into (numerator, denominator)."""
    try:
        parts = ts_string.split("/")
        return int(parts[0]), int(parts[1])
    except (ValueError, IndexError):
        raise ValueError(f"Invalid time signature: {ts_string}. Expected format: 4/4")


def beats_to_bar_beat(beat: int, beats_per_bar: int) -> str:
    """
    Convert beat number to bar.beat notation (Ableton style).

    Beat 0 = 1 (bar 1, beat 1 - just the bar number)
    Beat 1 = 1.2
    Beat 4 = 2 (in 4/4 time)
    Beat -1 = -1.4
    Beat -4 = -1
    Beat -5 = -2.4

    No bar 0 - goes from -1 to 1.
    First beat of each bar omits the .1 suffix.
    """
    if beat >= 0:
        bar = beat // beats_per_bar + 1
    else:
        # For negative beats, no +1 offset (skips bar 0)
        bar = beat // beats_per_bar
    beat_in_bar = beat % beats_per_bar + 1

    # Ableton style: first beat of bar is just the bar number
    if beat_in_bar == 1:
        return str(bar)
    return f"{bar}.{beat_in_bar}"


def beat_to_seconds(target_beat: float, warp_markers: list) -> float:
    """
    Convert a beat position to seconds using warp markers.

    Warp markers define a piecewise linear mapping from beats to seconds.
    We interpolate/extrapolate based on surrounding markers.
    """
    if len(warp_markers) < 2:
        raise ValueError("Need at least 2 warp markers for interpolation")

    # Sort markers by beats
    markers = sorted(warp_markers, key=lambda m: m.beats)

    # Find the two markers surrounding the target beat
    for i in range(len(markers) - 1):
        m1, m2 = markers[i], markers[i + 1]

        if m1.beats <= target_beat <= m2.beats:
            # Interpolate between these markers
            if m2.beats == m1.beats:
                return m1.seconds
            seconds_per_beat = (m2.seconds - m1.seconds) / (m2.beats - m1.beats)
            return m1.seconds + (target_beat - m1.beats) * seconds_per_beat

    # Extrapolate if target_beat is outside marker range
    if target_beat < markers[0].beats:
        # Extrapolate backwards using first two markers
        m1, m2 = markers[0], markers[1]
        seconds_per_beat = (m2.seconds - m1.seconds) / (m2.beats - m1.beats)
        return m1.seconds + (target_beat - m1.beats) * seconds_per_beat
    else:
        # Extrapolate forwards using last two markers
        m1, m2 = markers[-2], markers[-1]
        seconds_per_beat = (m2.seconds - m1.seconds) / (m2.beats - m1.beats)
        return m2.seconds + (target_beat - m2.beats) * seconds_per_beat


def seconds_to_beat(target_seconds: float, warp_markers: list) -> float:
    """
    Convert a time position to beats using warp markers.
    Extrapolates using first two markers if before the first marker.
    """
    if len(warp_markers) < 2:
        raise ValueError("Need at least 2 warp markers for interpolation")

    markers = sorted(warp_markers, key=lambda m: m.beats)
    m1, m2 = markers[0], markers[1]

    # Calculate beats per second from first two markers
    beats_per_second = (m2.beats - m1.beats) / (m2.seconds - m1.seconds)

    # Extrapolate: beat = m1.beats + (target_seconds - m1.seconds) * beats_per_second
    return m1.beats + (target_seconds - m1.seconds) * beats_per_second


def get_beat_range(warp_markers: list, end_marker: float) -> tuple[int, int]:
    """
    Determine the range of beats to generate.

    Returns (start_beat, end_beat) as integers.
    Extrapolates to time 0 for start, uses end_marker for end.
    """
    if not warp_markers:
        raise ValueError("No warp markers found")

    # Calculate the beat at time 0 (start of audio)
    beat_at_zero = seconds_to_beat(0.0, warp_markers)

    # Round down to get the first complete beat at or after time 0
    # (we'll filter out negative times later)
    min_beat = int(beat_at_zero)
    if beat_at_zero < 0:
        min_beat = int(beat_at_zero)  # Include partial negative beat

    # Always use end_marker for the end
    return min_beat, int(end_marker) + 1


def convert_asd_to_csv(asd_path: str, output_path: str, time_signature: str = "4/4"):
    """
    Convert an Ableton .asd file to a Sonic Visualizer compatible CSV.

    Args:
        asd_path: Path to the .asd file
        output_path: Path for the output CSV file
        time_signature: Time signature string (e.g., "4/4", "3/4", "6/8")
    """
    beats_per_bar, _ = parse_time_signature(time_signature)

    # Parse the .asd file (dummy sr and num_samples since we don't need get_time_map)
    clip = Clip(asd_path, sr=44100, num_samples=0)

    if not clip.warp_markers:
        raise ValueError("No warp markers found in the .asd file")

    print(f"Found {len(clip.warp_markers)} warp markers")
    for i, m in enumerate(clip.warp_markers):
        print(f"  [{i}] seconds={m.seconds:.3f}, beats={m.beats:.3f}")

    # Determine beat range (uses end_marker when warp markers only define BPM)
    start_beat, end_beat = get_beat_range(clip.warp_markers, clip.end_marker)
    print(f"end_marker: {clip.end_marker:.3f} beats")
    print(f"\nGenerating beats from {start_beat} to {end_beat - 1}")
    print(f"Time signature: {time_signature} ({beats_per_bar} beats per bar)")

    # Generate beat events
    events = []
    for beat in range(start_beat, end_beat):
        time_seconds = beat_to_seconds(beat, clip.warp_markers)
        label = beats_to_bar_beat(beat, beats_per_bar)
        events.append((time_seconds, label))

    # Filter out negative times (before audio starts)
    events = [(t, l) for t, l in events if t >= 0]

    # Sort by time
    events.sort(key=lambda x: x[0])

    # Write CSV
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Time", "Label"])
        for time_sec, label in events:
            writer.writerow([f"{time_sec:.6f}", label])

    print(f"\nWrote {len(events)} beat events to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert Ableton .asd warp data to Sonic Visualizer CSV"
    )
    parser.add_argument("asd_file", help="Path to the Ableton .asd file")
    parser.add_argument(
        "-o",
        "--output",
        help="Output CSV file path (default: same name as input with .csv extension)",
    )
    parser.add_argument(
        "-t", "--time-signature", default="4/4", help="Time signature (default: 4/4)"
    )

    args = parser.parse_args()

    # Validate input file
    asd_path = Path(args.asd_file)
    if not asd_path.exists():
        print(f"Error: File not found: {asd_path}")
        sys.exit(1)

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        name = asd_path.name
        output_name = name + ".csv"
        output_path = asd_path.parent / output_name

    try:
        convert_asd_to_csv(str(asd_path), str(output_path), args.time_signature)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
