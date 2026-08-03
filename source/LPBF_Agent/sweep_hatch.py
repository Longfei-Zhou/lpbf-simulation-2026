#!/usr/bin/env python3
"""Run a controlled hatch-spacing sweep at 200 W, 0.8 m/s, and 30 µm.

Each hatch value follows the complete CLI-to-score pipeline. All other inputs
remain fixed, while the snapshot window and sampling times are recomputed for
each scan path.

Usage::
    python3 sweep_hatch.py \\
        --cli 70=/path/to/70um_3dbenchy.cli \\
        --cli 80=/path/to/80um_3dbenchy.cli \\
        --cli 100=/path/to/100um_3dbenchy.cli \\
        --cli 120=/path/to/120um_3dbenchy.cli \\
        --thesis-bin /path/to/3DThesis \\
        --workdir source/Test/316L/sweep_hatch \\
        --output-dir source/LPBF_Agent/results/sweep_hatch

    python3 sweep_hatch.py ... --dry-run
    python3 sweep_hatch.py ... --resume
    python3 sweep_hatch.py ... --report-only
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent.parent
SOURCE_DIR = HERE.parent

# Experimental anchor with matching power, speed, and layer thickness.
LITERATURE_HATCH_UM = 120.0
LITERATURE_LABEL = "Metals 2021 11(5) 832 · 200 W / 800 mm/s / 30 µm → 99.8% dense"
LITERATURE_DOI = "10.3390/met11050832"

SNAPSHOT_WINDOW_M = (1.9e-3, 1.8e-3)
# Keep the approximately 334 µm melt-pool tail inside the snapshot window.
SNAPSHOT_MARGIN_M = 0.30e-3
SNAPSHOT_COUNT = 20


def load_cli_trans(path: Path):
    spec = importlib.util.spec_from_file_location("cli_trans_for_sweep", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Can't load {path}")
    module = importlib.util.module_from_spec(spec)
    # dataclass resolves module metadata through sys.modules during import.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@dataclass
class LayerGeometry:
    z_mm: float
    hatch_count: int
    hatch_angle_deg: float
    hatch_spacing_um: float
    contour_entities: int
    contour_segments: int
    bbox_mm: tuple[float, float, float, float]
    layer_count: int = 0

    def describe(self) -> str:
        return (
            f"z={self.z_mm:.3f} mm | {self.hatch_count} infill @ "
            f"{self.hatch_angle_deg:.2f}° / {self.hatch_spacing_um:.2f} µm · "
            f"{self.contour_entities} contours ({self.contour_segments} segs)"
        )


def measure_layer(cli_trans, cli_path: Path, layer_z_mm: float | None,
                  layer_number: int | None = None) -> LayerGeometry:
    units, layers = cli_trans.parse_cli(cli_path)
    nonempty = cli_trans.nonempty_layers(layers)
    if not nonempty:
        raise ValueError(f"{cli_path.name}: no non-empty layers were found")

    if layer_number is not None:
        index = layer_number - 1 if layer_number > 0 else len(nonempty) + layer_number
        if not 0 <= index < len(nonempty):
            raise ValueError(
                f"{cli_path.name}: only {len(nonempty)} non-empty layers are available; "
                f"layer {layer_number} cannot be selected"
            )
        target = nonempty[index]
    else:
        target = min(nonempty,
                     key=lambda layer: abs(layer.raw_z * units - layer_z_mm))
        if abs(target.raw_z * units - layer_z_mm) > 1e-6:
            raise ValueError(
                f"{cli_path.name}: the nearest layer to z={layer_z_mm} mm is at "
                f"z={target.raw_z * units:.6f} mm. If this is a uniform offset caused by "
                "a different Netfabb part-start height, select the equivalent section with "
                "--layer-number."
            )
    z_mm = target.raw_z * units

    hatches = np.array(
        [seg for e in target.entities if e.kind == "hatch" for seg in e.segments],
        dtype=float,
    ) * units
    if len(hatches) < 3:
        raise ValueError(
            f"{cli_path.name}: z={z_mm} mm contains only {len(hatches)} hatch segments"
        )

    angles = np.degrees(
        np.arctan2(hatches[:, 3] - hatches[:, 1], hatches[:, 2] - hatches[:, 0])
    ) % 180.0
    angle = float(np.median(angles))
    spread = float(np.max(np.abs(((angles - angle + 90.0) % 180.0) - 90.0)))
    if spread > 1.0:
        raise ValueError(
            f"{cli_path.name}: hatch directions are inconsistent (spread {spread:.2f}°). "
            "This script assumes one hatch direction; island or checkerboard strategies "
            "require separate handling."
        )

    # Cluster collinear fragments before measuring spacing because holes and
    # islands can split one physical track into nearly coincident pieces.
    theta = math.radians(angle)
    normal = np.array([-math.sin(theta), math.cos(theta)])
    projected = np.sort(hatches[:, 0:2] @ normal)
    cluster_tol_mm = 5e-3
    centers = [projected[0]]
    current = [projected[0]]
    for value in projected[1:]:
        if value - current[-1] <= cluster_tol_mm:
            current.append(value)
        else:
            centers[-1] = float(np.mean(current))
            centers.append(value)
            current = [value]
    centers[-1] = float(np.mean(current))
    gaps = np.diff(np.asarray(centers, dtype=float))
    spacing_um = float(np.median(gaps)) * 1e3 if len(gaps) else float("nan")

    contours = [e for e in target.entities if e.kind == "contour"]
    all_segments = np.array(
        [seg for e in target.entities for seg in e.segments], dtype=float
    ) * units

    return LayerGeometry(
        layer_count=len(nonempty),
        z_mm=z_mm,
        hatch_count=len(hatches),
        hatch_angle_deg=angle,
        hatch_spacing_um=spacing_um,
        contour_entities=len(contours),
        contour_segments=sum(len(e.segments) for e in contours),
        bbox_mm=(
            float(all_segments[:, [0, 2]].min()),
            float(all_segments[:, [0, 2]].max()),
            float(all_segments[:, [1, 3]].min()),
            float(all_segments[:, [1, 3]].max()),
        ),
    )


def validate_inputs(
    cli_trans,
    cli_map: dict[float, Path],
    layer_z_mm: float | None,
    spacing_tol_um: float,
    layer_number: int | None = None,
) -> dict[float, LayerGeometry]:
    """Validate spacing and cross-section equivalence before simulation."""
    print("=" * 72)
    print("Stage 1/5  Validating input CLI files")
    print("=" * 72)

    geometries: dict[float, LayerGeometry] = {}
    problems: list[str] = []
    for hatch_um in sorted(cli_map):
        cli_path = cli_map[hatch_um]
        if not cli_path.is_file():
            problems.append(f"h={hatch_um:g} um: file not found: {cli_path}")
            continue
        geometry = measure_layer(cli_trans, cli_path, layer_z_mm, layer_number)
        geometries[hatch_um] = geometry
        deviation = abs(geometry.hatch_spacing_um - hatch_um)
        status = "OK " if deviation <= spacing_tol_um else "BAD "
        print(f"  [{status}] h={hatch_um:6.1f} µm  {geometry.describe()}")
        print(f"          {cli_path}")
        if deviation > spacing_tol_um:
            problems.append(
                f"h={hatch_um:g} um: measured hatch spacing {geometry.hatch_spacing_um:.2f} µm, "
                f"deviation {deviation:.2f} um exceeds tolerance {spacing_tol_um:g} um"
            )

    if len(geometries) >= 2:
        reference_h = sorted(geometries)[0]
        reference = geometries[reference_h]
        for hatch_um, geometry in geometries.items():
            if geometry.contour_segments != reference.contour_segments:
                problems.append(
                    f"h={hatch_um:g} um: contour segment count {geometry.contour_segments} != "
                    f"h={reference_h:g} um: {reference.contour_segments}. The section itself "
                    "changed, so this is not a controlled hatch-only comparison."
                )
            # Compare dimensions rather than position because paths are recentered.
            def span(box):
                return (box[1] - box[0], box[3] - box[2])

            deltas = [abs(a - b) for a, b in zip(span(geometry.bbox_mm),
                                                 span(reference.bbox_mm))]
            if max(deltas) > 1e-3:
                problems.append(
                    f"h={hatch_um:g} um: section size "
                    f"{span(geometry.bbox_mm)[0]:.3f}×{span(geometry.bbox_mm)[1]:.3f} mm "
                    f"versus h={reference_h:g} um: "
                    f"{span(reference.bbox_mm)[0]:.3f}×{span(reference.bbox_mm)[1]:.3f} mm "
                    f" differs by {max(deltas) * 1e3:.1f} um. Not the same cross-section."
                )
            if geometry.layer_count != reference.layer_count:
                problems.append(
                    f"h={hatch_um:g} um: non-empty layer count {geometry.layer_count} != "
                    f"h={reference_h:g} um: {reference.layer_count}. Different layer counts can "
                    "make the same layer number refer to different sections."
                )
        angles = {h: g.hatch_angle_deg for h, g in geometries.items()}
        if max(angles.values()) - min(angles.values()) > 1.0:
            problems.append(
                f"Hatch angle differs between inputs: {angles}. With a 113° interlayer "
                "rotation, different starting angles produce different selected-layer angles."
            )

        z_values = {h: round(g.z_mm, 6) for h, g in geometries.items()}
        if len(set(z_values.values())) > 1:
            print(f"\n  Note: absolute z differs between inputs: {z_values}")
            print("  Normally a uniform shift from a different part start height "
                  "setting in Netfabb.")
            print("  Harmless: CLI_Trans moves the selected layer to z=0; absolute height never")
            print("  enters the simulation. The section itself was checked item by item above.")

    if problems:
        print("\nValidation failed:")
        for item in problems:
            print(f"  ✗ {item}")
        raise SystemExit(
            "\nInputs do not satisfy the controlled-comparison precondition; "
            "stopped before running any simulation."
            "\nIf this is acceptable, bypass with --skip-validation "
            "(not recommended)."
        )
    print(
        "\n  Same layer, same section, same contours across all inputs; "
        "only hatch differs"
        "\n  -- controlled-comparison precondition holds."
    )
    return geometries


@dataclass
class PathTiming:
    """Expanded path timing and trajectory used for snapshot selection."""

    xy: np.ndarray
    duration: np.ndarray
    time: np.ndarray
    powered: np.ndarray
    total_time: float

    # Dense beam trajectory columns: time, x, y, interval, powered.
    samples: np.ndarray = field(repr=False, default=None)


def read_path_timing(path_txt: Path, sample_dt: float = 2.0e-5) -> PathTiming:
    frame = pd.read_csv(path_txt, sep=r"\s+")
    xy = frame[["X(mm)", "Y(mm)"]].to_numpy(float) * 1e-3
    mode = frame["Mode"].to_numpy(int)
    pmod = frame["Pmod"].to_numpy(float)
    speed_or_time = frame.iloc[:, 5].to_numpy(float)

    step = np.r_[0.0, np.linalg.norm(np.diff(xy, axis=0), axis=1)]
    # Mode 0 stores velocity; Mode 1 stores duration directly.
    duration = np.where(
        mode == 0,
        np.divide(step, np.where(speed_or_time > 0, speed_or_time, 1.0),
                  out=np.zeros_like(step), where=speed_or_time > 0),
        speed_or_time,
    )
    duration[0] = 0.0
    time = np.cumsum(duration)
    powered = (mode == 0) & (pmod > 0)

    chunks = []
    for index in range(1, len(xy)):
        count = max(2, int(duration[index] / sample_dt))
        fraction = (np.arange(count) + 0.5) / count
        points = xy[index - 1] + (xy[index] - xy[index - 1]) * fraction[:, None]
        chunks.append(
            np.c_[
                time[index - 1] + fraction * duration[index],
                points,
                np.full(count, duration[index] / count),
                np.full(count, 1.0 if powered[index] else 0.0),
            ]
        )

    return PathTiming(
        xy=xy,
        duration=duration,
        time=time,
        powered=powered,
        total_time=float(time[-1]),
        samples=np.vstack(chunks),
    )


def powered_bbox(timing: PathTiming) -> tuple[float, float, float, float]:
    index = np.where(timing.powered)[0]
    points = np.vstack([timing.xy[index - 1], timing.xy[index]])
    return (
        float(points[:, 0].min()), float(points[:, 0].max()),
        float(points[:, 1].min()), float(points[:, 1].max()),
    )


def choose_snapshot_window(
    timing: PathTiming,
    window_m: tuple[float, float] = SNAPSHOT_WINDOW_M,
    margin_m: float = SNAPSHOT_MARGIN_M,
    step_m: float = 50e-6,
    dwell_tolerance: float = 0.95,
) -> tuple[float, float, float, float, dict]:
    """Select a snapshot window using identical criteria for every hatch.

    Maximize powered dwell inside the margin first. Among candidates within
    ``dwell_tolerance`` of that maximum, choose the largest dwell-weighted
    10th-to-90th-percentile time span so snapshots cover heat accumulation
    instead of clustering in a long-tail visit.
    """
    width, height = window_m
    samples = timing.samples
    x_lo, x_hi = samples[:, 1].min(), samples[:, 1].max()
    y_lo, y_hi = samples[:, 2].min(), samples[:, 2].max()

    candidates = []
    for x_start in np.arange(x_lo - 0.2e-3, x_hi - width + 0.2e-3, step_m):
        for y_start in np.arange(y_lo - 0.2e-3, y_hi - height + 0.2e-3, step_m):
            inside = (
                (samples[:, 1] >= x_start + margin_m)
                & (samples[:, 1] <= x_start + width - margin_m)
                & (samples[:, 2] >= y_start + margin_m)
                & (samples[:, 2] <= y_start + height - margin_m)
                & (samples[:, 4] > 0)
            )
            if not inside.any():
                continue
            picked = samples[inside]
            order = np.argsort(picked[:, 0])
            times = picked[order, 0]
            cumulative = np.cumsum(picked[order, 3])
            total = cumulative[-1]
            low, high = np.searchsorted(cumulative, [0.10 * total, 0.90 * total])
            spread = float(
                times[min(high, len(times) - 1)] - times[min(low, len(times) - 1)]
            )
            candidates.append((float(total), spread, float(x_start), float(y_start)))
    if not candidates:
        raise ValueError(
            "No window position leaves sufficient beam margin; check the path or reduce margin"
        )

    table = np.array(candidates)
    best_dwell = table[:, 0].max()
    eligible = table[table[:, 0] >= dwell_tolerance * best_dwell]
    dwell, spread, x_start, y_start = eligible[np.argmax(eligible[:, 1])]

    diagnostics = {
        "candidate_count": len(table),
        "best_dwell_s": float(best_dwell),
        "selected_dwell_s": float(dwell),
        "selected_dwell_p10_p90_span_s": float(spread),
        "selected_dwell_p10_p90_fraction": float(spread / timing.total_time),
        "center_mm": [
            float((x_start + width / 2) * 1e3),
            float((y_start + height / 2) * 1e3),
        ],
    }
    return x_start, x_start + width, y_start, y_start + height, diagnostics


def choose_scan_fracs(
    timing: PathTiming,
    window: tuple[float, float, float, float],
    count: int = SNAPSHOT_COUNT,
    margin_m: float = SNAPSHOT_MARGIN_M,
    min_separation_frac: float = 0.002,
) -> tuple[list[float], dict]:
    """Choose margin-safe snapshots by cumulative powered dwell.

    Dwell quantiles distribute samples across intermittent window visits; a
    minimum time separation rejects near-duplicate melt-pool states. ScanFracs
    is used because the scorer reads this field when validating evidence.
    """
    x_min, x_max, y_min, y_max = window
    samples = timing.samples
    eligible = samples[
        (samples[:, 1] >= x_min + margin_m)
        & (samples[:, 1] <= x_max - margin_m)
        & (samples[:, 2] >= y_min + margin_m)
        & (samples[:, 2] <= y_max - margin_m)
        & (samples[:, 4] > 0)
    ]
    if len(eligible) < count:
        raise ValueError(
            f"The window contains only {len(eligible)} eligible samples; cannot select "
            f"{count} snapshot times"
        )

    order = np.argsort(eligible[:, 0])
    times = eligible[order, 0]
    weights = eligible[order, 3]
    cumulative = np.cumsum(weights)
    total_dwell = float(cumulative[-1])

    targets = (np.arange(count) + 0.5) / count * total_dwell
    picked = times[np.searchsorted(cumulative, targets).clip(0, len(times) - 1)]

    minimum_gap = min_separation_frac * timing.total_time
    chosen: list[float] = []
    for value in np.sort(picked):
        if not chosen or value - chosen[-1] >= minimum_gap:
            chosen.append(float(value))
    # Refill rejected slots from the largest available time gap.
    pool = np.sort(times)
    while len(chosen) < count:
        gaps = np.diff(chosen) if len(chosen) > 1 else np.array([np.inf])
        insert_at = int(np.argmax(gaps)) if len(chosen) > 1 else 0
        low = chosen[insert_at] if len(chosen) else float(pool[0])
        high = chosen[insert_at + 1] if len(chosen) > 1 else float(pool[-1])
        window_pool = pool[(pool > low + minimum_gap) & (pool < high - minimum_gap)]
        if not len(window_pool):
            break
        chosen.append(float(window_pool[len(window_pool) // 2]))
        chosen.sort()

    fractions = sorted(float(t / timing.total_time * 100.0) for t in chosen[:count])
    separations = np.diff(fractions) if len(fractions) > 1 else np.array([0.0])
    diagnostics = {
        "eligible_sample_count": int(len(eligible)),
        "eligible_dwell_s": total_dwell,
        "eligible_time_span_s": float(times[-1] - times[0]),
        "requested_count": count,
        "selected_count": len(fractions),
        "min_separation_percent": float(separations.min()),
        "first_frac_percent": fractions[0],
        "last_frac_percent": fractions[-1],
        "total_scan_time_s": timing.total_time,
    }
    return fractions, diagnostics


CASE_TEMPLATE_FILES = (
    "Material.txt", "Beam.txt", "Settings.txt", "Output.txt", "ParamInput.txt",
)


def copy_template_files(template: Path, case_dir: Path,
                        max_threads: int | None = None) -> None:
    """Copy a case template and optionally override its OpenMP thread cap."""
    for name in CASE_TEMPLATE_FILES:
        shutil.copy2(template / name, case_dir / name)
    if max_threads is None:
        return
    settings = case_dir / "Settings.txt"
    text = settings.read_text(encoding="utf-8")
    patched, count = re.subn(r"(MaxThreads\s+)\d+", rf"\g<1>{max_threads}", text)
    if count == 0:
        patched = text.rstrip() + f"\n\nCompute\n{{\n    MaxThreads {max_threads}\n}}\n"
    settings.write_text(patched, encoding="utf-8")


def write_domain_regular(
    path: Path, x: tuple[float, float], y: tuple[float, float],
    z: tuple[float, float], res: tuple[float, float, float], header: str,
) -> None:
    path.write_text(
        f"{header}\n"
        f"X\n{{\n    Min {x[0]:.9g}\n    Max {x[1]:.9g}\n    Res {res[0]:.9g}\n}}\n\n"
        f"Y\n{{\n    Min {y[0]:.9g}\n    Max {y[1]:.9g}\n    Res {res[1]:.9g}\n}}\n\n"
        f"Z\n{{\n    Min {z[0]:.9g}\n    Max {z[1]:.9g}\n    Res {res[2]:.9g}\n}}\n",
        encoding="utf-8",
    )


def build_snapshot_case(
    case_dir: Path, template: Path, path_txt: Path,
    window: tuple[float, float, float, float], fractions: list[float],
    z_depth_m: float, res_xy_m: float, res_z_m: float, hatch_um: float,
    max_threads: int | None = None,
) -> None:
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "Data").mkdir(exist_ok=True)
    copy_template_files(template, case_dir, max_threads)
    shutil.copy2(path_txt, case_dir / "Path.txt")

    x_min, x_max, y_min, y_max = window
    write_domain_regular(
        case_dir / "Domain.txt",
        (x_min, x_max), (y_min, y_max), (-z_depth_m, 0.0),
        (res_xy_m, res_xy_m, res_z_m),
        f"# Generated by sweep_hatch.py for h = {hatch_um:g} µm.\n"
        f"# The window maximizes powered dwell beyond a "
        f"{SNAPSHOT_MARGIN_M * 1e6:.0f} µm boundary margin, then time coverage.",
    )

    (case_dir / "Mode.txt").write_text(
        f"# Generated by sweep_hatch.py for h = {hatch_um:g} µm.\n"
        f"# {len(fractions)} dwell-weighted snapshots keep the beam at least "
        f"{SNAPSHOT_MARGIN_M * 1e6:.0f} µm from the window boundary.\n"
        f"Snapshots\n{{\n    ScanFracs "
        + ",".join(f"{value:.6f}" for value in fractions)
        + "\n    Tracking None\n}\n",
        encoding="utf-8",
    )


def evaluation_box(
    bbox: tuple[float, float, float, float], margin_m: float,
    region_m: float | None, center_m: tuple[float, float] | None,
) -> tuple[tuple[float, float, float, float], bool]:
    """Return the full powered bounds or a resolution-preserving subregion.

    The complete scan path still supplies heat sources; only the area where
    fusion is evaluated is reduced, so a restricted result is local evidence.
    """
    x_min, x_max, y_min, y_max = bbox
    if region_m is None:
        return (x_min - margin_m, x_max + margin_m,
                y_min - margin_m, y_max + margin_m), False

    cx, cy = center_m if center_m else ((x_min + x_max) / 2, (y_min + y_max) / 2)
    half = region_m / 2
    # Keep the requested subregion inside the powered layer bounds.
    cx = min(max(cx, x_min + half), x_max - half) if x_max - x_min > region_m else (x_min + x_max) / 2
    cy = min(max(cy, y_min + half), y_max - half) if y_max - y_min > region_m else (y_min + y_max) / 2
    return (cx - half, cx + half, cy - half, cy + half), True


def build_solidification_case(
    case_dir: Path, template: Path, path_txt: Path,
    bbox: tuple[float, float, float, float], margin_m: float,
    grid_xy_m: float, layer_thickness_m: float, z_levels: int, hatch_um: float,
    region_m: float | None = None,
    region_center_m: tuple[float, float] | None = None,
    max_threads: int | None = None,
) -> None:
    """Build the solidification case before generating coverage points."""
    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "Data").mkdir(exist_ok=True)
    copy_template_files(template, case_dir, max_threads)
    shutil.copy2(template / "Mode.txt", case_dir / "Mode.txt")
    shutil.copy2(path_txt, case_dir / "Path.txt")

    (x_min, x_max, y_min, y_max), restricted = evaluation_box(
        bbox, margin_m, region_m, region_center_m)
    res_z = layer_thickness_m / max(z_levels - 1, 1)
    scope = (
        f"# Evaluation subregion: {(x_max - x_min) * 1e3:.2f} x "
        f"{(y_max - y_min) * 1e3:.2f} mm, centered at "
        f"({(x_min + x_max) / 2 * 1e3:.3f}, "
        f"{(y_min + y_max) / 2 * 1e3:.3f}) mm.\n"
        "# The full scan path remains active; conclusions are local to this region.\n"
        if restricted else
        "# Evaluation region: full powered-layer bounds plus the margin.\n"
    )
    header = (scope + (
        f"# Generated by sweep_hatch.py for h = {hatch_um:g} µm.\n"
        "# 3DThesis evaluates Custom points, while the scorer uses X/Y/Z to "
        "reconstruct the target grid.\n"
        "# X/Y bounds and resolution must match coverage_target_points.txt.\n"
    ))
    (case_dir / "Domain.txt").write_text(
        header
        + f"X\n{{\n    Min {x_min:.9g}\n    Max {x_max:.9g}\n"
          f"    Res {grid_xy_m:.9g}\n}}\n\n"
        + f"Y\n{{\n    Min {y_min:.9g}\n    Max {y_max:.9g}\n"
          f"    Res {grid_xy_m:.9g}\n}}\n\n"
        + f"Z\n{{\n    Min {-layer_thickness_m:.9g}\n    Max 0.0\n"
          f"    Res {res_z:.9g}\n}}\n\n"
        + "Custom\n{\n    File coverage_target_points.txt\n}\n",
        encoding="utf-8",
    )


def run(command: list[str], cwd: Path | None = None, dry_run: bool = False) -> None:
    printable = " ".join(str(part) for part in command)
    print(f"    $ {printable}" + (f"   (cwd={cwd})" if cwd else ""))
    if dry_run:
        return
    started = time.time()
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if result.returncode != 0:
        sys.stdout.write(result.stdout[-4000:])
        sys.stderr.write(result.stderr[-4000:])
        raise SystemExit(f"command failed (returncode={result.returncode}): {printable}")
    print(f"      done in {time.time() - started:.1f} s")


def has_thesis_output(case_dir: Path, kind: str) -> bool:
    data = case_dir / "Data"
    if not data.is_dir():
        return False
    if kind == "snapshots":
        return len(list(data.glob("*.Snapshot.*.csv"))) >= SNAPSHOT_COUNT
    return any(data.glob("*.Solidification.Final*.csv"))


def run_one_hatch(args, cli_trans, hatch_um: float, cli_path: Path,
                  geometry: LayerGeometry) -> dict:
    tag = f"h{hatch_um:g}"
    workdir = args.workdir / tag
    workdir.mkdir(parents=True, exist_ok=True)
    path_txt = workdir / "Path.txt"
    coverage_points = workdir / "coverage_target_points.txt"
    case_a = workdir / "A_snapshots"
    case_b = workdir / "B_solidification"
    result_dir = args.output_dir / tag

    print(f"\n{'-' * 72}\nh = {hatch_um:g} µm   ({cli_path.name})\n{'-' * 72}")

    if path_txt.is_file() and args.resume:
        print("  [skip] Path.txt already exists")
    else:
        command = [
            sys.executable, str(SOURCE_DIR / "CLI_Trans.py"),
            str(cli_path), str(path_txt),
            # Equivalent slices can differ by a global CLI start-height offset.
            "--layer-z-mm", str(geometry.z_mm),
            "--whole-layer", "--include-contours", "--recenter-xy",
            "--scan-speed", str(args.scan_speed),
            "--contour-speed", str(args.contour_speed),
            "--jump-speed", str(args.jump_speed),
            "--jump-delay", str(args.jump_delay),
        ]
        run(command, dry_run=args.dry_run)
    if args.dry_run and not path_txt.is_file():
        print("  [dry-run] Path.txt not generated yet; later stages only print commands")
        return {"hatch_um": hatch_um, "dry_run": True}

    timing = read_path_timing(path_txt)
    bbox = powered_bbox(timing)
    print(f"  Path: {int(timing.powered.sum())} powered segments, total "
          f"{timing.total_time * 1e3:.1f} ms")

    x_min, x_max, y_min, y_max, window_info = choose_snapshot_window(timing)
    fractions, frac_info = choose_scan_fracs(timing, (x_min, x_max, y_min, y_max))
    print(f"  Snapshot window centre ({window_info['center_mm'][0]:.3f}, "
          f"{window_info['center_mm'][1]:.3f}) mm | dwell "
          f"{window_info['selected_dwell_s'] * 1e3:.1f} ms | dwell spread "
          f"{window_info['selected_dwell_p10_p90_fraction'] * 100:.1f}% of scan")
    print(f"  Snapshot times {frac_info['first_frac_percent']:.2f}% .. "
          f"{frac_info['last_frac_percent']:.2f}%")

    # Build Domain.txt first because coverage generation reads its grid.
    build_snapshot_case(
        case_a, args.template_a, path_txt, (x_min, x_max, y_min, y_max),
        fractions, args.snapshot_depth_um * 1e-6,
        args.snapshot_res_xy_um * 1e-6, args.snapshot_res_z_um * 1e-6, hatch_um,
        max_threads=args.max_threads,
    )
    build_solidification_case(
        case_b, args.template_b, path_txt, bbox,
        args.margin_mm * 1e-3, args.grid_um * 1e-6,
        args.layer_thickness_um * 1e-6, args.z_levels, hatch_um,
        region_m=(args.region_mm * 1e-3) if args.region_mm else None,
        # Center both evidence sources on the same high-dwell material region.
        region_center_m=((x_min + x_max) / 2, (y_min + y_max) / 2),
        max_threads=args.max_threads,
    )
    print(f"  Cases assembled: {case_a.name} / {case_b.name}")

    case_points = case_b / "coverage_target_points.txt"
    if case_points.is_file() and args.resume:
        print("  [skip] coverage_target_points.txt already exists")
    else:
        run([
            sys.executable, str(HERE / "generate_coverage_points.py"),
            "--case-dir", str(case_b),
            "--domain-file", str(case_b / "Domain.txt"),
            "--path-file", str(path_txt),
            "--source-cli", str(cli_path),
            "--layer-thickness-um", str(args.layer_thickness_um),
            "--hatch-spacing-um", str(hatch_um),
            "--output", str(case_points),
        ], dry_run=args.dry_run)
    if case_points.is_file():
        shutil.copy2(case_points, coverage_points)

    for case_dir, kind in ((case_a, "snapshots"), (case_b, "solidification")):
        if args.resume and has_thesis_output(case_dir, kind):
            print(f"  [skip] {case_dir.name} already has output")
            continue
        print(f"  Running 3DThesis: {case_dir.name}  (~20 min expected)")
        run([str(args.thesis_bin), "ParamInput.txt"], cwd=case_dir,
            dry_run=args.dry_run)

    assessment_json = result_dir / "assessment.json"
    if args.resume and assessment_json.is_file():
        print("  [skip] assessment.json already exists")
    else:
        run([
            sys.executable, str(HERE / "run_all.py"),
            "--solidification-dir", str(case_b),
            "--snapshots-dir", str(case_a),
            "--source-cli", str(cli_path),
            "--layer-thickness-um", str(args.layer_thickness_um),
            "--hatch-spacing-um", str(hatch_um),
            "--layer-id", f"L1600-v{args.scan_speed}-h{hatch_um:g}",
            "--output-dir", str(result_dir),
        ], dry_run=args.dry_run)

    record = {
        "hatch_um": hatch_um,
        "cli": str(cli_path),
        "layer_z_mm": geometry.z_mm,
        "cli_measured_hatch_um": geometry.hatch_spacing_um,
        "hatch_count": geometry.hatch_count,
        "hatch_angle_deg": geometry.hatch_angle_deg,
        "powered_segments": int(timing.powered.sum()),
        "scan_time_s": timing.total_time,
        "snapshot_window": window_info,
        "scan_fracs": frac_info,
        "assessment_json": str(assessment_json),
    }
    (workdir / "sweep_stage.json").write_text(
        json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return record


def collect_results(output_dir: Path, hatches: list[float]) -> pd.DataFrame:
    rows = []
    for hatch_um in sorted(hatches):
        assessment_path = output_dir / f"h{hatch_um:g}" / "assessment.json"
        if not assessment_path.is_file():
            continue
        data = json.loads(assessment_path.read_text(encoding="utf-8"))
        metrics = data.get("metrics", {})
        coverage = data.get("coverage", {})
        rows.append({
            "hatch_um": hatch_um,
            "decision": data.get("decision"),
            "quality_score": data.get("quality_score"),
            "grade": data.get("grade"),
            "width_median_um": _um(metrics.get("melt_pool_width_median_m")),
            "depth_median_um": _um(metrics.get("melt_pool_depth_median_m")),
            "width_p10_um": _um(_lof(metrics, "nominal", "width_m")),
            "depth_p10_um": _um(_lof(metrics, "nominal", "depth_m")),
            "lof_nominal": metrics.get("lof_index_nominal"),
            "lof_midpoint": metrics.get("lof_index_midpoint"),
            "lof_conservative": metrics.get("lof_index_conservative"),
            "lof_resolution_sensitive": metrics.get("lof_resolution_sensitive"),
            "cells_across_pool": metrics.get("minimum_cells_across_p10_pool_dimension"),
            "top_coverage": coverage.get("coverage_fraction"),
            "interface_coverage": metrics.get("interface_melt_coverage_fraction"),
            "depth_over_width": metrics.get("keyhole_aspect_ratio"),
            "remelt_ge_2": metrics.get("remelt_fraction_ge_2"),
            "flags": "; ".join(data.get("flags", [])),
        })
    return pd.DataFrame(rows)


def _um(value):
    return None if value is None else float(value) * 1e6


def _lof(metrics: dict, label: str, key: str):
    return (metrics.get("lof_sensitivity") or {}).get(label, {}).get(key)


def critical_hatch(frame: pd.DataFrame, column: str) -> float | None:
    """Interpolate the first measured crossing of LOF(h) = 1."""
    if column not in frame.columns or "hatch_um" not in frame.columns:
        return None
    usable = frame.dropna(subset=[column]).sort_values("hatch_um")
    if len(usable) < 2:
        return None
    h = usable["hatch_um"].to_numpy(float)
    y = usable[column].to_numpy(float) - 1.0
    for i in range(len(h) - 1):
        if y[i] == 0.0:
            return float(h[i])
        if y[i] * y[i + 1] < 0.0:
            return float(h[i] + (h[i + 1] - h[i]) * (-y[i]) / (y[i + 1] - y[i]))
    return None


def _make_plot_english(frame: pd.DataFrame, output: Path, plt) -> bool:
    """Render the hatch-sweep plot with portable English labels."""
    figure, axis = plt.subplots(figsize=(8.5, 5.5))
    h = frame["hatch_um"]
    band = frame.dropna(subset=["lof_nominal", "lof_conservative"])
    if len(band):
        axis.fill_between(band["hatch_um"], band["lof_nominal"],
                          band["lof_conservative"], alpha=0.18, color="tab:blue",
                          label="nominal-conservative band (+-1 grid cell)")
    for column, style, label in (("lof_nominal", "o--", "nominal"),
                                 ("lof_midpoint", "s-", "midpoint (scored)"),
                                 ("lof_conservative", "^:", "conservative")):
        if frame[column].notna().any():
            axis.plot(h, frame[column], style, label=f"LOF {label}")
    axis.axhline(1.0, color="crimson", lw=1.6, label="Tang criterion = 1")
    axis.axvline(LITERATURE_HATCH_UM, color="seagreen", ls="--", lw=1.6,
                 label=f"literature {LITERATURE_HATCH_UM:.0f} um (99.8% dense)")
    axis.set_xlabel("hatch spacing (um)")
    axis.set_ylabel("LOF index  $(h/W)^2+(t/D)^2$")
    axis.set_title("316L / 200 W / 0.8 m/s / 30 um layer - benchy layer 1600")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(output, dpi=160)
    plt.close(figure)
    return True


def make_plot(frame: pd.DataFrame, output: Path) -> bool:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  [skip] matplotlib not installed; no plot")
        return False
    return _make_plot_english(frame, output, plt)


def write_report(frame: pd.DataFrame, output: Path, plot_name: str | None,
                 stages: list[dict], args) -> None:
    critical = {
        label: critical_hatch(frame, column)
        for label, column in (
            ("Nominal", "lof_nominal"),
            ("Half-cell", "lof_midpoint"),
            ("Conservative", "lof_conservative"),
        )
    }

    lines: list[str] = []
    w = lines.append
    w("# Hatch Sweep: Model-Predicted Window vs Literature Experiment")
    w("")
    w(f"> 316L · 200 W · 0.8 m/s · 30 µm layer · benchy layer 1600 ("
      f"{_report_z(stages, args)})")
    w(">")
    w("> Full-chain rerun: each hatch spacing is sliced, converted to a path, simulated "
      "with 3DThesis, and scored independently. Melt-pool geometry therefore changes with "
      "hatch spacing and includes adjacent-track preheating.")
    w("")
    w("---")
    w("")
    w("## 1. Question Addressed")
    w("")
    w("The model previously estimated a critical hatch spacing near 81 µm and consequently "
      "suggested 70 µm. This sweep tests that conclusion directly.")
    w(f"The external anchor, {LITERATURE_LABEL} (doi:{LITERATURE_DOI}), reports "
      f"99.8% density at h={LITERATURE_HATCH_UM:.0f} µm with matched power, speed, and "
      "layer thickness.")
    w("")
    w("The approximately 1.5× difference can reflect conservative model behavior or hidden "
      "experimental differences. The sweep plots the model response against the literature point.")
    w("")
    w("**Why this matters:** it requires no machine experiment and informs later interpretation "
      "of single-track metallography. A systematically conservative result implies that "
      "`Depth_Z` and/or `Efficiency` will likely need calibration toward a larger melt pool.")
    w("")
    w("## 2. Results")
    w("")
    if frame.empty:
        w("_No results are available. Rerun with `--report-only` after all cases finish._")
    else:
        w("| h (µm) | Decision | Quality | Width P10 | Depth P10 | Nominal LOF | "
          "Half-cell LOF | Conservative LOF | Top coverage | Interface fusion |")
        w("|---|---|---|---|---|---|---|---|---|---|")
        for _, row in frame.iterrows():
            w(
                f"| {row['hatch_um']:.0f} | {_text(row['decision'])} | "
                f"{_num(row['quality_score'], '{:.1f}')} | "
                f"{_num(row['width_p10_um'], '{:.1f}')} | "
                f"{_num(row['depth_p10_um'], '{:.1f}')} | "
                f"{_num(row['lof_nominal'], '{:.3f}')} | "
                f"{_num(row['lof_midpoint'], '{:.3f}')} | "
                f"{_num(row['lof_conservative'], '{:.3f}')} | "
                f"{_pct(row['top_coverage'])} | {_pct(row['interface_coverage'])} |"
            )
        w("")
        w("The nominal LOF uses the P10 melt-pool dimensions directly. The conservative LOF "
          "subtracts one grid cell from both width and depth. The half-cell result lies between "
          "them and is used for scoring. Reporting all three makes discretization uncertainty explicit.")
    w("")

    if plot_name:
        w(f"![LOF vs hatch]({plot_name})")
        w("")

    w("## 3. Critical Hatch Spacing")
    w("")
    w("Linear interpolation between measured sweep points is used to solve LOF(h) = 1:")
    w("")
    w("| Estimate | Critical h (µm) | Ratio to literature 120 µm |")
    w("|---|---|---|")
    for label, value in critical.items():
        ratio = f"{value / LITERATURE_HATCH_UM:.2f}×" if value else "—"
        w(f"| {label} | {f'{value:.1f}' if value else 'No crossing of 1'} | {ratio} |")
    w("")
    w("> `W` and `D` vary with hatch spacing because of adjacent-track preheating. The "
      "critical value therefore cannot be solved analytically from `(h/W)²+(t/D)²=1`; it must "
      "be interpolated between simulated points. Add hatch values for a more precise crossing.")
    w("")
    w("## 4. Interpretation")
    w("")
    reference = critical.get("Half-cell")
    if reference is None:
        w("The LOF curve does not cross 1 within the sweep. If every value exceeds 1, the "
          "model predicts inadequate fusion throughout the range and an undersized `Depth_Z` "
          "should be investigated. If every value is below 1, extend the sweep to larger h.")
    elif reference < 0.85 * LITERATURE_HATCH_UM:
        w(f"Critical h ≈ {reference:.0f} µm, substantially below the literature value of "
          f"{LITERATURE_HATCH_UM:.0f} µm.")
        w("")
        w("**This is direct evidence of conservative model behavior.** The most likely source "
          "is an undersized predicted melt pool, controlled by two uncalibrated parameters:")
        w("")
        w("- `Depth_Z = 3.5e-5` (volumetric heat-source penetration depth): changing "
          "20→60 µm moves predicted depth by 41%")
        w("- `Efficiency = 0.35` (absorptivity): Trapp 2017 measured approximately twofold "
          "variation with power and speed")
        w("")
        w("Single-track metallography is expected to move one or both parameters toward a "
          "larger melt pool. Until then, FAIL results and the 70 µm hatch suggestion are "
          "**lower-bound outputs of an uncalibrated model**, not direct process recommendations.")
    elif reference > 1.15 * LITERATURE_HATCH_UM:
        w(f"Critical h ≈ {reference:.0f} µm, substantially above the literature value of "
          f"{LITERATURE_HATCH_UM:.0f} µm. The model is optimistic and likely overpredicts melt-pool "
          "size. This is the more dangerous direction because it can pass parameters with physical "
          "LOF. Do not treat PASS as release evidence before calibration.")
    else:
        w(f"Critical h ≈ {reference:.0f} µm, within 15% of the literature value of "
          f"{LITERATURE_HATCH_UM:.0f} µm.")
        w("")
        w("**This is the strongest current external validation of the pipeline.** It supports "
          "the relative trend and overlap criterion even while `Depth_Z` and `Efficiency` retain "
          "literature defaults. It does not validate absolute melt depth or peak temperature.")
    w("")

    w("## 5. Controlled-Comparison Checks")
    w("")
    w("The following quantities remain fixed across cases:")
    w("")
    w(f"- Slice: {_report_z(stages, args)}; same STL and placement. The script verifies "
      "matching contour-segment count and bounding-box dimensions.")
    w(f"- Timing: scan {args.scan_speed} m/s, contour {args.contour_speed} m/s, "
      f"jump {args.jump_speed} m/s plus {args.jump_delay:g} s delay")
    w("- Material / Beam / Settings: copied unchanged from template cases")
    w(f"- Coverage points: {args.grid_um:g} µm grid and {args.z_levels} z levels; each "
      "case passes its hatch spacing explicitly. Automatic inference is unsafe because 119 "
      "contours can make `infer_hatch_spacing()` return 21.33 µm.")
    w("")
    w("Two quantities **must be recomputed for each h**:")
    w("")
    if stages:
        w("| h (µm) | Snapshot-window center (mm) | Window dwell (ms) | Dwell span | "
          "Snapshot-time range | Powered segments | Layer time (ms) |")
        w("|---|---|---|---|---|---|---|")
        for record in sorted(stages, key=lambda item: item.get("hatch_um", 0)):
            if record.get("dry_run"):
                continue
            window = record.get("snapshot_window", {})
            fracs = record.get("scan_fracs", {})
            center = window.get("center_mm", [float("nan")] * 2)
            w(
                f"| {record['hatch_um']:.0f} | "
                f"({center[0]:.2f}, {center[1]:.2f}) | "
                f"{window.get('selected_dwell_s', 0) * 1e3:.1f} | "
                f"{window.get('selected_dwell_p10_p90_fraction', 0) * 100:.0f}% | "
                f"{fracs.get('first_frac_percent', 0):.1f}–"
                f"{fracs.get('last_frac_percent', 0):.1f}% | "
                f"{record.get('powered_segments', 0)} | "
                f"{record.get('scan_time_s', 0) * 1e3:.1f} |"
            )
        w("")
    w("`choose_snapshot_window()` and `choose_scan_fracs()` apply the **same rule** to every h: "
      "maximize powered dwell while keeping the beam more than 300 µm from window boundaries, "
      "then prefer the widest visit-time span among near-ties. Reusing a manually tuned h=100 "
      "window can leave other cases with no scanned track in the window.")
    w("")

    w("## 6. Questions This Sweep Cannot Answer")
    w("")
    w("- **Whether absolute melt depth is correct:** it still depends on two uncalibrated "
      "parameters and requires single-track metallography")
    w("- **Physical density:** the LOF criterion excludes keyhole porosity, spatter, and balling")
    w("- **Multilayer accumulation:** one simulated layer cannot capture the long-range effect "
      "of 113° interlayer rotation")
    w("- **Cross-section representativeness:** the benchy top layer is an isolated 6.4×6.4 mm "
      "section with unusually weak heat sinking")
    w("")
    w("---")
    w("")
    w(f"_Generated by `sweep_hatch.py` · {time.strftime('%Y-%m-%d %H:%M')}_")

    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _report_z(stages: list[dict], args) -> str:
    """Describe equivalent slices that may have different global Z offsets."""
    values = sorted({round(float(s["layer_z_mm"]), 6)
                     for s in stages if s.get("layer_z_mm") is not None})
    if len(values) == 1:
        return f"z = {values[0]:.3f} mm"
    if not values:
        if args.layer_number is not None:
            return f"non-empty layer {args.layer_number}"
        return f"z = {args.layer_z_mm:.3f} mm"
    listed = " / ".join(f"{v:.3f}" for v in values)
    label = (f"non-empty layer {args.layer_number}" if args.layer_number is not None
             else f"z ≈ {args.layer_z_mm:.3f} mm")
    return (f"{label} (absolute CLI z values: {listed} mm; uniform translation only; "
            "sections verified equivalent)")


def _text(value):
    return "—" if value is None or (isinstance(value, float) and math.isnan(value)) else str(value)


def _num(value, fmt="{:.3f}"):
    try:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return "—"
        return fmt.format(float(value))
    except (TypeError, ValueError):
        return "—"


def _pct(value):
    try:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return "—"
        return f"{float(value) * 100:.2f}%"
    except (TypeError, ValueError):
        return "—"


def parse_cli_map(entries: list[str]) -> dict[float, Path]:
    mapping: dict[float, Path] = {}
    for entry in entries:
        if "=" not in entry:
            raise SystemExit(
                f"--cli requires the format '<hatch_um>=<path>'; received: {entry}"
            )
        key, _, value = entry.partition("=")
        try:
            hatch_um = float(key)
        except ValueError as error:
            raise SystemExit(f"--cli hatch value is not numeric: {key}") from error
        mapping[hatch_um] = Path(value).expanduser()
    return mapping


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run a full-chain hatch-spacing sweep at 200 W / 0.8 m/s / 30 µm.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--cli", action="append", required=True, metavar="H=PATH",
        help="Repeatable '<hatch_um>=<CLI path>'; example: --cli 120=/path/120um.cli",
    )
    parser.add_argument("--thesis-bin", type=Path, default=None,
                        help="3DThesis executable; optional with --dry-run.")
    parser.add_argument("--workdir", type=Path,
                        default=SOURCE_DIR / "Test" / "316L" / "sweep_hatch",
                        help="Directory for Path.txt and generated case directories.")
    parser.add_argument("--output-dir", type=Path,
                        default=HERE / "results" / "sweep_hatch",
                        help="Directory for assessment results and the aggregate report.")
    parser.add_argument("--template-a", type=Path,
                        default=SOURCE_DIR / "Test" / "316L" / "A_200W_316L",
                        help="Snapshot-case template (Material/Beam/Settings/Output).")
    parser.add_argument("--template-b", type=Path,
                        default=SOURCE_DIR / "Test" / "316L" / "B_200W_316L",
                        help="Solidification-case template.")

    layer_group = parser.add_mutually_exclusive_group()
    layer_group.add_argument("--layer-z-mm", type=float, default=51.0,
                             help="Select by absolute z in mm (default: 51.0).")
    layer_group.add_argument("--layer-number", type=int, default=None,
                             help="Select the Nth non-empty layer (1-based; negative values "
                                  "count from the end; -1 is the top layer). Use when CLI files "
                                  "have equivalent sections at uniformly shifted absolute z values.")
    parser.add_argument("--scan-speed", type=float, default=0.8)
    parser.add_argument("--contour-speed", type=float, default=0.6)
    parser.add_argument("--jump-speed", type=float, default=5.0)
    parser.add_argument("--jump-delay", type=float, default=2.0e-4)

    parser.add_argument("--layer-thickness-um", type=float, default=30.0)
    parser.add_argument("--grid-um", type=float, default=25.0)
    parser.add_argument("--z-levels", type=int, default=3)
    parser.add_argument("--margin-mm", type=float, default=1.5)
    parser.add_argument(
        "--region-mm", type=float, default=None,
        help="Assess only an N x N mm region centered where the beam dwells longest. "
             "Use this for large sections that would require millions of coverage points. "
             "The full layer path is still simulated, but the conclusion applies only to the region.",
    )

    parser.add_argument("--snapshot-res-xy-um", type=float, default=10.0)
    parser.add_argument("--snapshot-res-z-um", type=float, default=5.0)
    parser.add_argument("--snapshot-depth-um", type=float, default=80.0)

    parser.add_argument(
        "--max-threads", type=int, default=None,
        help="Override MaxThreads in generated Settings files. Templates may contain local "
             "machine values (A=8, B=13), which otherwise limit server CPU use.")
    parser.add_argument("--spacing-tol-um", type=float, default=2.0,
                        help="Allowed difference between measured and declared CLI hatch spacing.")
    parser.add_argument("--skip-validation", action="store_true")
    parser.add_argument("--resume", action="store_true",
                        help="Skip stages with existing output.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate inputs and print commands without running 3DThesis.")
    parser.add_argument("--report-only", action="store_true",
                        help="Regenerate only the aggregate outputs from completed cases.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.workdir = args.workdir.expanduser().resolve()
    args.output_dir = args.output_dir.expanduser().resolve()
    cli_map = parse_cli_map(args.cli)

    # Fail before assembling cases if the simulator is unavailable.
    if not args.report_only and not args.dry_run:
        if args.thesis_bin is None:
            raise SystemExit(
                "--thesis-bin must point at the 3DThesis executable (or use --dry-run).\n"
                "In this repository it is normally 3DThesis/build/bin/3DThesis"
            )
        thesis = Path(args.thesis_bin).expanduser()
        if not thesis.is_file():
            hint = ""
            for guess in (PROJECT_ROOT / "3DThesis" / "build" / "bin" / "3DThesis",
                          PROJECT_ROOT / "3DThesis" / "build" / "install" / "bin" / "3DThesis",
                          PROJECT_ROOT / "build" / "bin" / "3DThesis"):
                if guess.is_file():
                    hint = f"\nFound one in this repository, try:\n  --thesis-bin {guess}"
                    break
            raise SystemExit(f"--thesis-bin does not exist: {thesis}{hint}")
        import os
        if not os.access(thesis, os.X_OK):
            raise SystemExit(
                f"{thesis} is not executable; run:\n  chmod +x {thesis}")
        args.thesis_bin = thesis.resolve()

    cli_trans = load_cli_trans(SOURCE_DIR / "CLI_Trans.py")

    stages: list[dict] = []
    if not args.report_only:
        geometries = (
            {h: measure_layer(cli_trans, p, args.layer_z_mm, args.layer_number)
             for h, p in cli_map.items()}
            if args.skip_validation
            else validate_inputs(cli_trans, cli_map, args.layer_z_mm,
                                 args.spacing_tol_um, args.layer_number)
        )
        print(f"\nStage 2-4/5  Full chain per hatch"
              f"({len(cli_map)} points x 2 cases, ~{len(cli_map) * 40} min expected)")
        for hatch_um in sorted(cli_map):
            stages.append(
                run_one_hatch(args, cli_trans, hatch_um, cli_map[hatch_um],
                              geometries[hatch_um])
            )
    else:
        for hatch_um in sorted(cli_map):
            stage_file = args.workdir / f"h{hatch_um:g}" / "sweep_stage.json"
            if stage_file.is_file():
                stages.append(json.loads(stage_file.read_text(encoding="utf-8")))

    print(f"\n{'=' * 72}\nStage 5/5  Aggregating\n{'=' * 72}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = collect_results(args.output_dir, list(cli_map))

    if frame.empty and not args.dry_run:
        print("  No assessment.json available to aggregate yet")
    csv_path = args.output_dir / "hatch_sweep.csv"
    frame.to_csv(csv_path, index=False)
    print(f"  wrote {csv_path}")

    plot_name = None
    if not frame.empty:
        plot_path = args.output_dir / "hatch_sweep_lof.png"
        if make_plot(frame, plot_path):
            plot_name = plot_path.name
            print(f"  wrote {plot_path}")

    report_path = args.output_dir / "Hatch_Sweep_Report.md"
    write_report(frame, report_path, plot_name, stages, args)
    print(f"  wrote {report_path}")

    if not frame.empty:
        print("\n" + frame[[
            "hatch_um", "decision", "lof_nominal", "lof_midpoint", "lof_conservative",
        ]].to_string(index=False))
        reference = critical_hatch(frame, "lof_midpoint")
        if reference:
            print(f"\n  Critical hatch (half-cell) ≈ {reference:.1f} µm  ·  "
                  f"literature {LITERATURE_HATCH_UM:.0f} µm  ·  "
                  f"ratio {reference / LITERATURE_HATCH_UM:.2f}×")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
