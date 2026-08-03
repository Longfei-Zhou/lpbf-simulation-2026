from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


FLOAT_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


@dataclass
class PathEntity:
    kind: str
    segments: list[tuple[float, float, float, float]]


@dataclass
class Layer:
    raw_z: float
    entities: list[PathEntity] = field(default_factory=list)

    @property
    def hatch_count(self) -> int:
        return sum(len(e.segments) for e in self.entities if e.kind == "hatch")

    @property
    def contour_segment_count(self) -> int:
        return sum(len(e.segments) for e in self.entities if e.kind == "contour")


def numbers(text: str) -> list[float]:
    return [float(x) for x in FLOAT_RE.findall(text)]


def parse_cli(path: Path) -> tuple[float, list[Layer]]:
    units: float | None = None
    layers: list[Layer] = []
    current: Layer | None = None

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line_no, raw_line in enumerate(f, 1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith("$$UNITS/"):
                vals = numbers(line.split("/", 1)[1])
                if not vals:
                    raise ValueError(f"Cannot parse $$UNITS on line {line_no}")
                units = vals[0]
                continue

            if line.startswith("$$LAYER/"):
                vals = numbers(line.split("/", 1)[1])
                if not vals:
                    raise ValueError(f"Cannot parse $$LAYER on line {line_no}")
                current = Layer(raw_z=vals[0])
                layers.append(current)
                continue

            if line.startswith("$$HATCHES/"):
                if current is None:
                    raise ValueError(f"$$HATCHES before $$LAYER on line {line_no}")

                vals = numbers(line.split("/", 1)[1])
                if len(vals) < 2:
                    raise ValueError(f"Malformed $$HATCHES on line {line_no}")

                # CLI hatch records contain an ID, a count, then endpoint quartets.
                declared = int(round(vals[1]))
                coords = vals[2:]
                available = len(coords) // 4
                count = min(declared, available)

                if count != declared:
                    print(
                        f"Warning: line {line_no}: declared {declared} hatches "
                        f"but found {available}; using {count}.",
                        file=sys.stderr,
                    )

                segments = [
                    tuple(coords[i:i + 4])  # type: ignore[arg-type]
                    for i in range(0, count * 4, 4)
                ]
                current.entities.append(PathEntity("hatch", segments))
                continue

            if line.startswith("$$POLYLINE/"):
                if current is None:
                    raise ValueError(f"$$POLYLINE before $$LAYER on line {line_no}")

                vals = numbers(line.split("/", 1)[1])
                if len(vals) < 3:
                    raise ValueError(f"Malformed $$POLYLINE on line {line_no}")

                # CLI polyline records contain an ID, direction, count, and points.
                direction = int(round(vals[1]))
                declared_points = int(round(vals[2]))
                coords = vals[3:]
                available_points = len(coords) // 2
                point_count = min(declared_points, available_points)

                if point_count < 2:
                    continue

                pts = [
                    (coords[i], coords[i + 1])
                    for i in range(0, point_count * 2, 2)
                ]

                segments = [
                    (pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1])
                    for i in range(len(pts) - 1)
                ]

                # Closed contour indicators are normally 0 or 1. If Netfabb
                # omitted a repeated final point, close it explicitly.
                if direction in (0, 1) and pts[-1] != pts[0]:
                    segments.append((pts[-1][0], pts[-1][1], pts[0][0], pts[0][1]))

                current.entities.append(PathEntity("contour", segments))

    if units is None:
        raise ValueError("The CLI file does not contain $$UNITS.")
    if not layers:
        raise ValueError("The CLI file does not contain any $$LAYER records.")

    return units, layers


def nonempty_layers(layers: Iterable[Layer]) -> list[Layer]:
    return [layer for layer in layers if layer.entities]


def choose_layer(
    layers: list[Layer],
    units: float,
    layer_number: int | None,
    layer_z_mm: float | None,
) -> list[Layer]:
    """Select exactly one non-empty layer. layer_number is 1-based."""
    populated = nonempty_layers(layers)

    if layer_number is not None and layer_z_mm is not None:
        raise ValueError("Use either --layer-number or --layer-z-mm, not both.")

    if layer_number is None and layer_z_mm is None:
        raise ValueError("Choose one layer with --layer-number N or --layer-z-mm Z.")

    if layer_number is not None:
        if not 1 <= layer_number <= len(populated):
            raise ValueError(f"--layer-number must be between 1 and {len(populated)}.")
        return [populated[layer_number - 1]]

    assert layer_z_mm is not None
    return [min(populated, key=lambda layer: abs(layer.raw_z * units - layer_z_mm))]


def collect_segments(
    selected_layers: list[Layer],
    units: float,
    include_contours: bool,
    max_hatches: int | None,
) -> list[tuple[str, float, float, float, float, float]]:
    """Return ``(kind, x1_mm, y1_mm, x2_mm, y2_mm, z_mm)`` records."""
    output: list[tuple[str, float, float, float, float, float]] = []
    hatch_seen = 0

    for layer in selected_layers:
        z_mm = layer.raw_z * units

        for entity in layer.entities:
            if entity.kind == "contour" and not include_contours:
                continue

            for x1, y1, x2, y2 in entity.segments:
                if entity.kind == "hatch":
                    if max_hatches is not None and hatch_seen >= max_hatches:
                        return output
                    hatch_seen += 1

                output.append(
                    (
                        entity.kind,
                        x1 * units,
                        y1 * units,
                        x2 * units,
                        y2 * units,
                        z_mm,
                    )
                )

    return output

CONTIGUOUS_EPS_MM = 1.0e-3
MIN_POSITIONING_TIME_S = 1.0e-9


def write_path(
    path: Path,
    segments: list[tuple[str, float, float, float, float, float]],
    scan_speed: float,
    contour_speed: float | None,
    positioning_time: float | None,
    jump_speed: float,
    jump_delay: float,
    z_to_zero: bool,
    recenter_xy: bool,
    no_header: bool,
) -> dict[str, float]:
    if not segments:
        raise ValueError("No path segments were selected.")

    x_shift = segments[0][1] if recenter_xy else 0.0
    y_shift = segments[0][2] if recenter_xy else 0.0
    z_shift = segments[0][5] if z_to_zero else 0.0
    prev_x = prev_y = 0.0

    n_jumps = 0
    jump_dist_mm = 0.0
    jump_time_s = 0.0
    melt_time_s = 0.0
    melt_dist_mm = 0.0
    max_jump_mm = 0.0

    with path.open("w", encoding="utf-8", newline="\n") as f:
        if not no_header:
            f.write("Mode X(mm) Y(mm) Z(mm) Pmod Vel(m/s)/Time(s)\n")

        for kind, x1, y1, x2, y2, z in segments:
            speed = (
                contour_speed
                if kind == "contour" and contour_speed is not None
                else scan_speed
            )

            x1 -= x_shift
            y1 -= y_shift
            x2 -= x_shift
            y2 -= y_shift
            z -= z_shift

            distance_mm = math.hypot(x1 - prev_x, y1 - prev_y)

            if positioning_time is not None:
                t_jump = positioning_time
            elif distance_mm <= CONTIGUOUS_EPS_MM:
                # Contiguous segments only require a laser-state transition.
                t_jump = MIN_POSITIONING_TIME_S
            else:
                t_jump = distance_mm * 1.0e-3 / jump_speed + jump_delay

            if distance_mm > CONTIGUOUS_EPS_MM:
                n_jumps += 1
                jump_dist_mm += distance_mm
                max_jump_mm = max(max_jump_mm, distance_mm)
            jump_time_s += t_jump

            seg_mm = math.hypot(x2 - x1, y2 - y1)
            melt_dist_mm += seg_mm
            melt_time_s += seg_mm * 1.0e-3 / speed

            f.write(
                f"1 {x1:.9f} {y1:.9f} {z:.9f} "
                f"0 {t_jump:.9e}\n"
            )
            f.write(
                f"0 {x2:.9f} {y2:.9f} {z:.9f} "
                f"1 {speed:.9f}\n"
            )

            prev_x, prev_y = x2, y2

    return {
        "n_jumps": n_jumps,
        "jump_dist_mm": jump_dist_mm,
        "max_jump_mm": max_jump_mm,
        "jump_time_s": jump_time_s,
        "melt_dist_mm": melt_dist_mm,
        "melt_time_s": melt_time_s,
        "total_time_s": jump_time_s + melt_time_s,
    }


def print_summary(units: float, layers: list[Layer]) -> None:
    populated = nonempty_layers(layers)
    first = populated[0]
    last = populated[-1]

    total_hatches = sum(layer.hatch_count for layer in populated)
    total_contour_segments = sum(
        layer.contour_segment_count for layer in populated
    )

    print(f"CLI unit scale:          {units:g} mm per CLI unit")
    print(f"All $$LAYER records:     {len(layers)}")
    print(f"Non-empty layers:        {len(populated)}")
    print(f"First non-empty Z:       {first.raw_z * units:.6f} mm")
    print(f"Last non-empty Z:        {last.raw_z * units:.6f} mm")
    print(f"Total hatch segments:    {total_hatches}")
    print(f"Total contour segments:  {total_contour_segments}")



def print_layer_table(units: float, layers: list[Layer]) -> None:
    populated = nonempty_layers(layers)
    print("Layer number | Z (mm)     | Hatches | Contour segments")
    print("-------------+------------+---------+-----------------")
    for i, layer in enumerate(populated, 1):
        print(f"{i:12d} | {layer.raw_z * units:10.6f} | {layer.hatch_count:7d} | {layer.contour_segment_count:16d}")

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert one Netfabb CLI layer to the official 3DThesis Path.txt format."
    )
    parser.add_argument("input_cli", type=Path)
    parser.add_argument("output_path", type=Path)

    layer = parser.add_mutually_exclusive_group()
    layer.add_argument(
        "--layer-number",
        type=int,
        help="One-based number among non-empty layers, e.g. 500 means the 500th non-empty layer.",
    )
    layer.add_argument(
        "--layer-z-mm",
        type=float,
        help="Select the non-empty layer nearest this physical Z value in mm.",
    )

    parser.add_argument(
        "--max-hatches",
        type=int,
        default=None,
        help="Keep only the first N hatch segments. Omit it to export the complete selected layer.",
    )
    parser.add_argument(
        "--whole-layer",
        action="store_true",
        help="Export every hatch in the selected layer. This is also the default when --max-hatches is omitted.",
    )
    parser.add_argument(
        "--include-contours",
        action="store_true",
        help="Also convert $$POLYLINE contour segments.",
    )
    parser.add_argument(
        "--scan-speed",
        type=float,
        default=0.8,
        help="Hatch scan speed in m/s. Default: 0.8",
    )
    parser.add_argument(
        "--contour-speed",
        type=float,
        default=None,
        help="Contour speed in m/s; defaults to --scan-speed.",
    )
    parser.add_argument(
        "--jump-speed",
        type=float,
        default=5.0,
        help=(
            "Galvanometer jump speed in m/s, used to estimate Mode-1 travel "
            "time from distance. Use the measured machine value. Default: 5.0"
        ),
    )
    parser.add_argument(
        "--jump-delay",
        type=float,
        default=2.0e-4,
        help=(
            "Fixed per-jump overhead in seconds, including skywriting, laser "
            "switching, and settling. Use the measured machine value. "
            "Default: 2.0e-4"
        ),
    )
    parser.add_argument(
        "--positioning-time",
        "--jump-time",
        dest="positioning_time",
        type=float,
        default=None,
        help=(
            "Override jump speed and delay with one fixed Mode-1 duration in "
            "seconds. The legacy value 1e-7 approximates instantaneous motion "
            "and should only be used to reproduce historical results."
        ),
    )
    parser.add_argument(
        "--keep-z",
        action="store_true",
        help="Keep original physical Z instead of shifting the first selected layer to Z=0.",
    )
    parser.add_argument(
        "--recenter-xy",
        action="store_true",
        help="Shift the first selected path start point to X=0, Y=0.",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Do not write the column-name header line.",
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Inspect the CLI and print a summary without writing Path.txt.",
    )
    parser.add_argument(
        "--list-layers",
        action="store_true",
        help="Print all non-empty layer numbers, physical Z values and path counts.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    if args.scan_speed <= 0:
        raise ValueError("--scan-speed must be positive.")
    if args.contour_speed is not None and args.contour_speed <= 0:
        raise ValueError("--contour-speed must be positive.")
    if args.positioning_time is not None and args.positioning_time <= 0:
        raise ValueError("--positioning-time must be positive.")
    if args.jump_speed <= 0:
        raise ValueError("--jump-speed must be positive.")
    if args.jump_delay < 0:
        raise ValueError("--jump-delay must not be negative.")
    if args.max_hatches is not None and args.max_hatches <= 0:
        raise ValueError("--max-hatches must be positive.")
    if args.whole_layer and args.max_hatches is not None:
        raise ValueError("--whole-layer and --max-hatches cannot be used together.")

    units, layers = parse_cli(args.input_cli)
    print_summary(units, layers)

    if args.list_layers:
        print()
        print_layer_table(units, layers)
        return 0

    if args.summary_only:
        return 0

    selected = choose_layer(
        layers,
        units,
        args.layer_number,
        args.layer_z_mm,
    )
    segments = collect_segments(
        selected,
        units,
        args.include_contours,
        args.max_hatches,
    )

    hatch_output = sum(1 for s in segments if s[0] == "hatch")
    contour_output = sum(1 for s in segments if s[0] == "contour")

    stats = write_path(
        args.output_path,
        segments,
        args.scan_speed,
        args.contour_speed,
        args.positioning_time,
        args.jump_speed,
        args.jump_delay,
        z_to_zero=not args.keep_z,
        recenter_xy=args.recenter_xy,
        no_header=args.no_header,
    )

    selected_z = selected[0].raw_z * units
    print(f"Selected layer Z:           {selected_z:.6f} mm")
    print(f"Written hatch segments:     {hatch_output}")
    print(f"Written contour segments:   {contour_output}")
    print(f"Output rows:                 {2 * len(segments)}")
    print(f"Created:                     {args.output_path}")

    total = stats["total_time_s"]
    print()
    print("Time budget for this layer")
    if args.positioning_time is not None:
        print(f"  Jump model:               fixed {args.positioning_time:.3e} s per Mode-1 row")
    else:
        print(f"  Jump model:               dist / {args.jump_speed:g} m/s + {args.jump_delay:.3e} s")
    print(f"  Melting  {stats['melt_time_s'] * 1e3:9.3f} ms"
          f"   ({stats['melt_dist_mm']:8.3f} mm at {args.scan_speed:g} m/s)")
    print(f"  Jumping  {stats['jump_time_s'] * 1e3:9.3f} ms"
          f"   ({stats['jump_dist_mm']:8.3f} mm over {stats['n_jumps']} jumps,"
          f" longest {stats['max_jump_mm']:.3f} mm)")
    print(f"  Total    {total * 1e3:9.3f} ms"
          f"   (jumping is {100.0 * stats['jump_time_s'] / total:.1f}% of layer time)")
    print()
    print("Official path semantics:")
    print("  Mode 1 -> final column is positioning time in seconds")
    print("  Mode 0 -> final column is scan velocity in m/s")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(2)
