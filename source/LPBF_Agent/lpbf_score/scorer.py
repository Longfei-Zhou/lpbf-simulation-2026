"""Physics-informed, field-adaptive scoring for one 3DThesis LPBF layer.

The output is a screening assessment, not an experimental qualification.
One 3DThesis run is treated as one physical build layer; z is depth within
the simulated domain and is never interpreted as a build-layer identifier.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml

from .diagnostics import build_process_diagnostics
from .reporting import write_intuitive_outputs


PROJECT_DIR = Path(__file__).resolve().parent.parent
ENGINE_VERSION = "2.1.0"

COMPONENT_LABELS = {
    "coverage": "Top-surface scan coverage",
    "fusion": "Interlayer fusion completeness",
    "fusion_margin": "LOF geometric margin",
    "keyhole_margin": "Keyhole-risk margin",
    "pool_consistency": "Melt-pool consistency",
    "thermal_uniformity": "Thermal-field consistency",
    "remelt": "Remelt control",
}


def finite_positive(values: Iterable[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return array[np.isfinite(array) & (array > 0)]


def clipped(value: float, low: float = 0.0, high: float = 100.0) -> float:
    return float(np.clip(value, low, high))


def logistic_score(margin: float, steepness: float) -> float:
    """Map a signed safety margin to 0–100; the boundary maps to 50."""
    exponent = float(np.clip(-steepness * margin, -60.0, 60.0))
    return 100.0 / (1.0 + math.exp(exponent))


def robust_log_mad(values: Iterable[float]) -> float | None:
    array = finite_positive(values)
    if len(array) < 3:
        return None
    logs = np.log10(array)
    median = np.median(logs)
    return float(np.median(np.abs(logs - median)))


def coefficient_of_variation(values: Iterable[float]) -> float | None:
    array = finite_positive(values)
    if len(array) < 2:
        return None
    mean = float(np.mean(array))
    if mean == 0:
        return None
    return float(np.std(array, ddof=1) / mean)


def numeric_summary(series: pd.Series) -> dict[str, float | int | None]:
    values = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan)
    valid = values.dropna()
    if valid.empty:
        return {"count": 0, "missing": int(len(series))}
    return {
        "count": int(valid.count()),
        "missing": int(values.isna().sum()),
        "min": float(valid.min()),
        "p05": float(valid.quantile(0.05)),
        "p10": float(valid.quantile(0.10)),
        "median": float(valid.median()),
        "mean": float(valid.mean()),
        "p90": float(valid.quantile(0.90)),
        "p95": float(valid.quantile(0.95)),
        "max": float(valid.max()),
        "std": float(valid.std(ddof=1)) if len(valid) > 1 else 0.0,
    }


def parse_grouped_text(path: Path) -> dict[str, float | str]:
    """Parse the simple nested-brace format used by 3DThesis input files."""
    if not path.is_file():
        return {}
    values: dict[str, float | str] = {}
    stack: list[str] = []
    pending_group: str | None = None
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line == "{":
            if pending_group:
                stack.append(pending_group)
                pending_group = None
            continue
        if line.endswith("{"):
            group = line[:-1].strip()
            if group:
                stack.append(group)
            continue
        if line == "}":
            if stack:
                stack.pop()
            continue
        tokens = line.split()
        if len(tokens) == 1:
            pending_group = tokens[0]
            continue
        key, raw_value = tokens[0], tokens[1]
        full_key = ".".join(stack + [key])
        try:
            values[full_key] = float(raw_value.rstrip(","))
        except ValueError:
            values[full_key] = raw_value
    return values


def find_key(values: dict[str, Any], suffix: str) -> Any | None:
    suffix_lower = suffix.lower()
    for key, value in values.items():
        if key.lower() == suffix_lower or key.lower().endswith("." + suffix_lower):
            return value
    return None


def data_directory(path: Path) -> tuple[Path, Path]:
    resolved = path.expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Input path does not exist: {resolved}")
    if resolved.is_file():
        return resolved.parent, resolved.parent
    child = resolved / "Data"
    return (child, resolved) if child.is_dir() else (resolved, resolved)


def find_solidification_files(path: Path) -> tuple[list[Path], Path]:
    data_dir, case_root = data_directory(path)
    if path.is_file():
        return [path.resolve()], case_root
    files = sorted(
        file
        for file in data_dir.glob("*.csv")
        if ".Solidification.Final" in file.name
    )
    if not files:
        files = sorted(
            file
            for file in data_dir.glob("*.csv")
            if "Solidification" in file.name and "Final" in file.name
        )
    if not files:
        raise FileNotFoundError(
            f"No *Solidification.Final*.csv file found under {data_dir}"
        )
    return files, case_root


def find_snapshot_files(path: Path | None) -> tuple[list[Path], Path | None]:
    if path is None:
        return [], None
    data_dir, case_root = data_directory(path)
    if path.is_file():
        return [path.resolve()], case_root
    files = sorted(file for file in data_dir.glob("*.csv") if ".Snapshot." in file.name)
    return files, case_root


def find_rdf_files(case_root: Path) -> list[Path]:
    data_dir = case_root / "Data" if (case_root / "Data").is_dir() else case_root
    return sorted(
        file
        for file in data_dir.glob("*.csv")
        if ".RDF.Final" in file.name
    )


def load_csvs(files: list[Path]) -> pd.DataFrame:
    frames = [pd.read_csv(file) for file in files]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def analyse_rdf(files: list[Path]) -> dict[str, Any] | None:
    """Summarise the ExaCA-compatible RDF output without mixing it into scoring."""
    if not files:
        return None
    frame = load_csvs(files)
    required = {"x", "y", "z", "tm", "tl", "cr"}
    missing = sorted(required - set(frame.columns))
    result: dict[str, Any] = {
        "files": [str(file) for file in files],
        "row_count": int(len(frame)),
        "available_columns": list(frame.columns),
        "missing_required_columns": missing,
        "used_in_quality_score": False,
    }
    if missing:
        return result
    melt_time = pd.to_numeric(frame["tm"], errors="coerce")
    solidification_time = pd.to_numeric(frame["tl"], errors="coerce")
    duration = solidification_time - melt_time
    finite_duration = pd.to_numeric(duration, errors="coerce").dropna()
    result.update(
        {
            "melt_time_s": numeric_summary(melt_time),
            "solidification_time_s": numeric_summary(solidification_time),
            "liquid_duration_s": numeric_summary(duration),
            "cooling_rate_k_s": numeric_summary(frame["cr"]),
            "negative_liquid_duration_fraction": (
                float(finite_duration.lt(0).mean())
                if len(finite_duration)
                else None
            ),
        }
    )
    return result


def audit_requested_outputs(
    case_root: Path,
    actual_columns: Iterable[str],
    rdf_files: list[Path],
) -> dict[str, Any]:
    """Compare Output.txt requests with fields that actually arrived."""
    requested = parse_grouped_text(case_root / "Output.txt")
    actual = set(actual_columns)
    field_map: dict[str, tuple[str, ...]] = {
        "Grid.x": ("x",),
        "Grid.y": ("y",),
        "Grid.z": ("z",),
        "Temperature.T": ("T",),
        "Solidification.tSol": ("tSol",),
        "Solidification.G": ("G",),
        "Solidification.Gx": ("Gx",),
        "Solidification.Gy": ("Gy",),
        "Solidification.Gz": ("Gz",),
        "Solidification.V": ("V",),
        "Solidification.dTdt": ("dTdt",),
        "Solidification.eqFrac": ("eqFrac",),
        "Solidification.depth": ("depth",),
        "Solidification.numMelt": ("numMelt",),
        "Solidification.RDF": ("__RDF_FILE__",),
        "Solidification.MP_Stats": (
            "MP_width",
            "MP_length",
        ),
        "Solidification+.H": ("H",),
        "Solidification+.Hx": ("Hx",),
        "Solidification+.Hy": ("Hy",),
        "Solidification+.Hz": ("Hz",),
    }
    enabled = {
        key
        for key, value in requested.items()
        if isinstance(value, (int, float)) and float(value) == 1.0
    }
    missing: list[str] = []
    satisfied: list[str] = []
    for key in sorted(enabled):
        expected = field_map.get(key)
        if expected is None:
            continue
        present = (
            bool(rdf_files)
            if expected == ("__RDF_FILE__",)
            else all(name in actual for name in expected)
        )
        (satisfied if present else missing).append(key)
    mode = parse_grouped_text(case_root / "Mode.txt")
    secondary = find_key(mode, "Solidification.Secondary")
    if secondary is None:
        secondary = find_key(mode, "Secondary")
    return {
        "output_file": str(case_root / "Output.txt"),
        "enabled_requests": sorted(enabled),
        "satisfied_requests": satisfied,
        "missing_requested_outputs": missing,
        "secondary_mode": secondary,
        "secondary_requested_but_disabled": bool(
            any(key.startswith("Solidification+.") for key in enabled)
            and secondary is not None
            and float(secondary) == 0.0
        ),
    }


def parse_path_file(path: Path) -> pd.DataFrame | None:
    if not path.is_file():
        return None
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 2:
        return None
    columns = lines[0].split()
    required = {"Mode", "X(mm)", "Y(mm)", "Pmod"}
    if len(columns) < 6 or not required.issubset(columns):
        return None
    records: list[list[float]] = []
    short_lines = 0
    trailing_token_lines = 0
    discarded_tokens = 0
    for line in lines[1:]:
        tokens = line.split()
        if not tokens:
            continue
        if len(tokens) < len(columns):
            short_lines += 1
            continue
        if len(tokens) > len(columns):
            trailing_token_lines += 1
            discarded_tokens += len(tokens) - len(columns)
        try:
            # 3DThesis reads one path record from each physical line. Taking
            # the first six tokens mirrors that behaviour when an extreme test
            # file contains concatenated trailing tokens on the same line.
            records.append([float(token) for token in tokens[: len(columns)]])
        except ValueError:
            short_lines += 1
    if not records:
        return None
    frame = pd.DataFrame(records, columns=columns)
    frame.attrs["path_parse_diagnostics"] = {
        "physical_record_count": int(len(records)),
        "short_or_non_numeric_line_count": int(short_lines),
        "trailing_token_line_count": int(trailing_token_lines),
        "discarded_trailing_token_count": int(discarded_tokens),
        "parsing_policy": "first record per physical line, matching 3DThesis",
        "extreme_concatenated_line_format_detected": bool(trailing_token_lines),
    }
    return frame


INFILL_DIRECTION_TOLERANCE_DEG = 10.0


def dominant_direction_mask(
    segments: np.ndarray,
    tolerance_deg: float = INFILL_DIRECTION_TOLERANCE_DEG,
) -> np.ndarray:
    """Select the parallel-infill subset of an (N, 4) x0,y0,x1,y1 segment array.

    Contour/border segments follow the part outline, so their directions are
    spread over the full circle while infill hatches are mutually parallel.
    Routines that assume "all powered segments are parallel hatches" (hatch
    spacing inference, CLI layer matching) break once contours are present:
    with 119 contour segments mixed into 94 hatches the inferred spacing
    collapsed from 100 um to 21 um.

    The dominant direction is a length-weighted circular mean on 2*theta (the
    doubling removes the 180-degree ambiguity of an undirected line), so long
    parallel hatches dominate short contour chords. Segments within
    `tolerance_deg` of that direction are kept. If fewer than three survive the
    input probably is not a hatch pattern at all and everything is kept, which
    reproduces the previous behaviour.
    """
    if len(segments) < 3:
        return np.ones(len(segments), dtype=bool)
    vectors = segments[:, 2:4] - segments[:, 0:2]
    lengths = np.linalg.norm(vectors, axis=1)
    valid = lengths > 0
    if valid.sum() < 3:
        return np.ones(len(segments), dtype=bool)
    angles = np.arctan2(vectors[:, 1], vectors[:, 0])
    weights = np.where(valid, lengths, 0.0)
    theta = 0.5 * math.atan2(
        float((weights * np.sin(2 * angles)).sum()),
        float((weights * np.cos(2 * angles)).sum()),
    )
    # Fold the deviation into [0, pi/2]: a line and its reverse are the same.
    deviation = np.abs(((angles - theta + math.pi / 2) % math.pi) - math.pi / 2)
    mask = valid & (deviation <= math.radians(tolerance_deg))
    if mask.sum() < 3:
        return np.ones(len(segments), dtype=bool)
    return mask


def infer_hatch_spacing(path_frame: pd.DataFrame | None) -> float | None:
    """Infer spacing between parallel powered line segments in metres."""
    if path_frame is None or len(path_frame) < 3:
        return None
    segments: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    for index in range(1, len(path_frame)):
        row = path_frame.iloc[index]
        if int(row["Mode"]) != 0 or float(row["Pmod"]) <= 0:
            continue
        start = path_frame.iloc[index - 1][["X(mm)", "Y(mm)"]].to_numpy(float) / 1000
        end = row[["X(mm)", "Y(mm)"]].to_numpy(float) / 1000
        vector = end - start
        if np.linalg.norm(vector) > 0:
            segments.append((start, end, vector))
    if len(segments) < 3:
        return None

    # Drop contour segments before estimating spacing; see
    # dominant_direction_mask() for why mixing them in destroys the estimate.
    stacked = np.array([[a[0], a[1], b[0], b[1]] for a, b, _ in segments], dtype=float)
    infill = dominant_direction_mask(stacked)
    segments = [seg for seg, keep in zip(segments, infill) if keep]

    angles = np.array([math.atan2(v[1], v[0]) for _, _, v in segments])
    lengths = np.array([np.linalg.norm(v) for _, _, v in segments])
    theta = 0.5 * math.atan2(
        float((lengths * np.sin(2 * angles)).sum()),
        float((lengths * np.cos(2 * angles)).sum()),
    )
    normal = np.array([-math.sin(theta), math.cos(theta)])
    offsets = np.sort(np.array([((a + b) / 2) @ normal for a, b, _ in segments]))
    differences = np.diff(offsets)
    differences = differences[differences > 1e-7]
    if len(differences) < 2:
        return None
    q25, q75 = np.quantile(differences, [0.25, 0.75])
    central = differences[(differences >= q25 * 0.75) & (differences <= q75 * 1.25)]
    return float(np.median(central if len(central) else differences))


def scan_speed(path_frame: pd.DataFrame | None) -> float | None:
    if path_frame is None:
        return None
    speed_column = next(
        (column for column in path_frame.columns if "Vel" in column or "Time" in column),
        None,
    )
    if speed_column is None:
        return None
    mask = (path_frame["Mode"] == 0) & (path_frame["Pmod"] > 0)
    speeds = finite_positive(path_frame.loc[mask, speed_column])
    return float(np.median(speeds)) if len(speeds) else None


def parse_cli_layer_thickness(path: Path | None) -> dict[str, Any] | None:
    """Read the median physical layer spacing from an ASCII CLI slice file."""
    if path is None or not path.is_file():
        return None
    text = path.read_text(encoding="utf-8", errors="replace")
    units_match = re.search(
        r"^\s*\$\$UNITS\s*/\s*([-+0-9.eE]+)",
        text,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    raw_layers = re.findall(
        r"^\s*\$\$LAYER\s*/\s*([-+0-9.eE]+)",
        text,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    if units_match is None or len(raw_layers) < 2:
        return None
    units_mm = float(units_match.group(1))
    layer_positions_m = np.sort(
        np.unique(np.asarray(raw_layers, dtype=float) * units_mm * 1e-3)
    )
    differences = np.diff(layer_positions_m)
    differences = differences[differences > 1e-12]
    if not len(differences):
        return None
    return {
        "path": str(path.resolve()),
        "units_mm_per_cli_unit": units_mm,
        "layer_count": int(len(layer_positions_m)),
        "layer_thickness_m": float(np.median(differences)),
        "minimum_spacing_m": float(np.min(differences)),
        "maximum_spacing_m": float(np.max(differences)),
    }


def parse_cli_layer_geometry(path: Path | None) -> dict[str, Any] | None:
    """Read CLI hatch segments and closed contours for solid/void masking."""
    if path is None or not path.is_file():
        return None
    number_pattern = re.compile(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
    )
    units_mm: float | None = None
    layers: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if line.startswith("$$UNITS/"):
                values = number_pattern.findall(line.split("/", 1)[1])
                units_mm = float(values[0]) if values else None
            elif line.startswith("$$LAYER/"):
                values = number_pattern.findall(line.split("/", 1)[1])
                if values:
                    current = {
                        "raw_z": float(values[0]),
                        "hatches": [],
                        "contours": [],
                    }
                    layers.append(current)
            elif current is not None and line.startswith("$$HATCHES/"):
                values = [
                    float(value)
                    for value in number_pattern.findall(line.split("/", 1)[1])
                ]
                if len(values) < 2:
                    continue
                count = min(int(round(values[1])), (len(values) - 2) // 4)
                current["hatches"].extend(
                    tuple(values[2 + 4 * index : 6 + 4 * index])
                    for index in range(count)
                )
            elif current is not None and line.startswith("$$POLYLINE/"):
                values = [
                    float(value)
                    for value in number_pattern.findall(line.split("/", 1)[1])
                ]
                if len(values) < 7:
                    continue
                count = min(int(round(values[2])), (len(values) - 3) // 2)
                points = np.asarray(values[3 : 3 + 2 * count], dtype=float).reshape(
                    (-1, 2)
                )
                if len(points) >= 3:
                    if not np.allclose(points[0], points[-1], rtol=0.0, atol=0.0):
                        points = np.vstack((points, points[0]))
                    current["contours"].append(points)
    if units_mm is None:
        return None
    scale = units_mm * 1e-3
    populated = []
    for layer in layers:
        if not layer["hatches"]:
            continue
        populated.append(
            {
                "raw_z": layer["raw_z"],
                "z_m": layer["raw_z"] * scale,
                "hatches": np.asarray(layer["hatches"], dtype=float) * scale,
                "contours": [
                    np.asarray(contour, dtype=float) * scale
                    for contour in layer["contours"]
                ],
            }
        )
    return {
        "path": str(path.resolve()),
        "units_mm_per_cli_unit": units_mm,
        "layers": populated,
    }


def match_cli_solid_geometry(
    path: Path | None,
    path_segments: np.ndarray,
    tolerance_m: float = 5e-8,
) -> tuple[dict[str, Any] | None, dict[str, Any]]:
    """Match a recentered Path to its CLI layer and return aligned contour rings."""
    parsed = parse_cli_layer_geometry(path)
    if parsed is None or not len(path_segments):
        return None, {
            "available": False,
            "applied": False,
            "reason": "source CLI geometry or powered path is unavailable",
        }
    # The matcher compares the Path against the CLI's $$HATCHES sequence, so
    # contour segments must be dropped first. With them included the 30 um
    # benchy layer (94 hatches + 119 contours) failed to match any layer and
    # the solid/void mask was silently unavailable; hatch-only matched fine.
    path_segments = np.asarray(path_segments, dtype=float)
    infill = dominant_direction_mask(path_segments)
    if infill.sum() >= 3:
        path_segments = path_segments[infill]
    best: tuple[float, int, dict[str, Any], np.ndarray] | None = None
    for layer_number, layer in enumerate(parsed["layers"], 1):
        hatches = layer["hatches"]
        if len(hatches) < len(path_segments):
            continue
        comparison_count = min(len(path_segments), len(hatches), 20)
        translation = path_segments[0, :2] - hatches[0, :2]
        aligned = hatches[:comparison_count].copy()
        aligned[:, :2] += translation
        aligned[:, 2:4] += translation
        rmse = float(
            np.sqrt(
                np.mean(
                    (aligned - path_segments[:comparison_count]) ** 2
                )
            )
        )
        candidate = (rmse, layer_number, layer, translation)
        if best is None or rmse < best[0]:
            best = candidate
    if best is None or best[0] > tolerance_m:
        return None, {
            "available": parsed is not None,
            "applied": False,
            "reason": "no CLI layer matches the powered Path segment sequence",
            "best_path_match_rmse_m": best[0] if best is not None else None,
            "match_tolerance_m": tolerance_m,
        }
    rmse, layer_number, layer, translation = best
    rings = [
        contour + translation
        for contour in layer["contours"]
    ]
    if not rings:
        return None, {
            "available": True,
            "applied": False,
            "reason": "matched CLI layer has no closed contour rings",
            "matched_layer_number": layer_number,
            "matched_layer_z_m": layer["z_m"],
            "best_path_match_rmse_m": rmse,
        }
    geometry = {"rings": rings}
    metadata = {
        "available": True,
        "applied": True,
        "method": "path corridor intersected with matched CLI even-odd contours",
        "source_cli": parsed["path"],
        "matched_layer_number": layer_number,
        "matched_layer_z_m": layer["z_m"],
        "matched_layer_hatch_count": int(len(layer["hatches"])),
        "matched_layer_contour_count": int(len(rings)),
        "path_segment_count": int(len(path_segments)),
        "path_match_rmse_m": rmse,
        "translation_m": [float(translation[0]), float(translation[1])],
    }
    return geometry, metadata


def cli_solid_region_mask(
    xy: np.ndarray,
    geometry: dict[str, Any] | None,
) -> np.ndarray:
    """Apply the CLI even-odd fill rule so holes and exterior voids are excluded."""
    if geometry is None or not len(xy):
        return np.ones(len(xy), dtype=bool)
    solid = np.zeros(len(xy), dtype=bool)
    x, y = xy[:, 0], xy[:, 1]
    for ring in geometry.get("rings", []):
        inside = np.zeros(len(xy), dtype=bool)
        x1, y1 = ring[-1]
        for x2, y2 in ring:
            crosses = (y1 > y) != (y2 > y)
            intersection_x = (
                (x2 - x1) * (y - y1) / (y2 - y1 + 1e-300) + x1
            )
            inside ^= crosses & (x < intersection_x)
            x1, y1 = x2, y2
        solid ^= inside
    return solid


def discover_source_cli(roots: Iterable[Path]) -> Path | None:
    """Find one unambiguous source CLI next to the selected 3DThesis cases."""
    candidates: set[Path] = set()
    for root in roots:
        for parent in (root, root.parent):
            candidates.update(path.resolve() for path in parent.glob("*.cli"))
            candidates.update(path.resolve() for path in parent.glob("*.CLI"))
    return next(iter(candidates)) if len(candidates) == 1 else None


def powered_segments(path_frame: pd.DataFrame | None) -> np.ndarray:
    """Return powered scan segments as (x0, y0, x1, y1) in metres."""
    if path_frame is None or len(path_frame) < 2:
        return np.empty((0, 4), dtype=float)
    segments: list[list[float]] = []
    for index in range(1, len(path_frame)):
        row = path_frame.iloc[index]
        try:
            powered = int(row["Mode"]) == 0 and float(row["Pmod"]) > 0
        except (KeyError, TypeError, ValueError):
            powered = False
        if not powered:
            continue
        start = path_frame.iloc[index - 1][["X(mm)", "Y(mm)"]].to_numpy(float) * 1e-3
        end = row[["X(mm)", "Y(mm)"]].to_numpy(float) * 1e-3
        if np.linalg.norm(end - start) > 1e-12:
            segments.append([start[0], start[1], end[0], end[1]])
    return np.asarray(segments, dtype=float).reshape((-1, 4))


def infill_segments(path_frame: pd.DataFrame | None) -> np.ndarray:
    """Powered segments with contour/border passes removed.

    Use this wherever the code assumes "one parallel hatch pattern":
    hatch-spacing inference and CLI layer matching. Do NOT use it for the
    coverage corridor — contour passes melt material too, so target points
    along the outline must stay in the target region.

    Contours are separated by scan speed first, because build strategies
    almost always expose the border at a different speed than the infill
    (here 0.6 vs 0.8 m/s). When several speeds are present the group with the
    largest total scanned length is taken as the infill. Direction filtering
    then removes any stragglers. If speed carries no information the direction
    filter alone is used.
    """
    if path_frame is None or len(path_frame) < 2:
        return np.empty((0, 4), dtype=float)
    rows: list[list[float]] = []
    speeds: list[float] = []
    for index in range(1, len(path_frame)):
        row = path_frame.iloc[index]
        try:
            powered = int(row["Mode"]) == 0 and float(row["Pmod"]) > 0
        except (KeyError, TypeError, ValueError):
            continue
        if not powered:
            continue
        start = path_frame.iloc[index - 1][["X(mm)", "Y(mm)"]].to_numpy(float) * 1e-3
        end = row[["X(mm)", "Y(mm)"]].to_numpy(float) * 1e-3
        if np.linalg.norm(end - start) <= 1e-12:
            continue
        rows.append([start[0], start[1], end[0], end[1]])
        try:
            speeds.append(float(row.iloc[5]))
        except (IndexError, TypeError, ValueError):
            speeds.append(float("nan"))
    segments = np.asarray(rows, dtype=float).reshape((-1, 4))
    if len(segments) < 3:
        return segments

    speed_array = np.round(np.asarray(speeds, dtype=float), 6)
    lengths = np.linalg.norm(segments[:, 2:4] - segments[:, 0:2], axis=1)
    unique = np.unique(speed_array[np.isfinite(speed_array)])
    if len(unique) > 1:
        totals = {
            float(value): float(lengths[speed_array == value].sum())
            for value in unique
        }
        infill_speed = max(totals, key=totals.get)
        by_speed = speed_array == infill_speed
        if by_speed.sum() >= 3:
            segments = segments[by_speed]

    direction = dominant_direction_mask(segments)
    if direction.sum() >= 3:
        segments = segments[direction]
    return segments


def target_region_mask(
    xy: np.ndarray,
    segments: np.ndarray,
    corridor_half_width_m: float,
) -> np.ndarray:
    """Rasterise the intended scan area as the union of path corridors."""
    if len(xy) == 0 or len(segments) == 0:
        return np.zeros(len(xy), dtype=bool)
    mask = np.zeros(len(xy), dtype=bool)
    radius = float(corridor_half_width_m)
    if not math.isfinite(radius) or radius <= 0:
        return mask

    # Restrict each distance calculation to the segment's padded bounding
    # box. This preserves the exact capsule geometry while avoiding an
    # O(domain_points * segments) full-domain temporary array per segment.
    for x0, y0, x1, y1 in segments:
        candidates = np.flatnonzero(
            (~mask)
            & (xy[:, 0] >= min(x0, x1) - radius)
            & (xy[:, 0] <= max(x0, x1) + radius)
            & (xy[:, 1] >= min(y0, y1) - radius)
            & (xy[:, 1] <= max(y0, y1) + radius)
        )
        if not len(candidates):
            continue
        start = np.array([x0, y0])
        vector = np.array([x1 - x0, y1 - y0])
        length_squared = float(vector @ vector)
        if length_squared <= 0:
            continue
        candidate_xy = xy[candidates]
        relative = candidate_xy - start
        # Explicit 2-D dot product avoids platform BLAS edge cases seen for
        # very small LPBF coordinates.
        projection = relative[:, 0] * vector[0] + relative[:, 1] * vector[1]
        fraction = np.clip(projection / length_squared, 0.0, 1.0)
        closest = start + fraction[:, None] * vector
        distance_squared = np.sum((candidate_xy - closest) ** 2, axis=1)
        mask[candidates[distance_squared <= radius**2]] = True
    return mask


def scan_fraction_metadata(case_root: Path | None) -> dict[str, Any]:
    if case_root is None:
        return {"values": [], "count": None, "source": None}
    mode_path = case_root / "Mode.txt"
    if not mode_path.is_file():
        return {"values": [], "count": None, "source": None}
    text = mode_path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"^\s*ScanFracs\s+([^\n#}]+)",
        text,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    if match is None:
        return {"values": [], "count": None, "source": str(mode_path)}
    values = [
        float(value)
        for value in re.findall(
            r"[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?", match.group(1)
        )
    ]
    return {"values": values, "count": len(values), "source": str(mode_path)}


def case_provenance(
    case_root: Path | None,
    data_files: list[Path],
) -> dict[str, Any]:
    """Detect whether present input text files post-date the supplied CSV results."""
    if case_root is None or not data_files:
        return {"consistent": None, "newer_configuration_files": []}
    oldest_data_time = min(path.stat().st_mtime for path in data_files)
    config_files = [
        case_root / name
        for name in (
            "Mode.txt",
            "Output.txt",
            "Material.txt",
            "Domain.txt",
            "Beam.txt",
            "Path.txt",
        )
        if (case_root / name).is_file()
    ]
    newer = [
        str(path)
        for path in config_files
        if path.stat().st_mtime > oldest_data_time + 60.0
    ]
    return {
        "consistent": not newer,
        "newer_configuration_files": newer,
        "oldest_data_mtime": oldest_data_time,
    }


def cross_case_compatibility(
    solid_root: Path,
    snapshot_root: Path | None,
) -> dict[str, Any]:
    """Check that solidification and snapshot cases use the same physics/path inputs."""
    if snapshot_root is None or solid_root.resolve() == snapshot_root.resolve():
        return {"consistent": True, "mismatched_files": []}
    mismatched: list[str] = []
    missing: list[str] = []
    parameter_differences: dict[str, dict[str, Any]] = {}
    for filename in ("Material.txt", "Beam.txt", "Path.txt"):
        solid_file = solid_root / filename
        snapshot_file = snapshot_root / filename
        if not solid_file.is_file() or not snapshot_file.is_file():
            missing.append(filename)
        elif solid_file.read_bytes() != snapshot_file.read_bytes():
            mismatched.append(filename)
            if filename in ("Material.txt", "Beam.txt"):
                solid_values = parse_grouped_text(solid_file)
                snapshot_values = parse_grouped_text(snapshot_file)
                for key in sorted(set(solid_values) | set(snapshot_values)):
                    solid_value = solid_values.get(key)
                    snapshot_value = snapshot_values.get(key)
                    if solid_value != snapshot_value:
                        parameter_differences[f"{filename}:{key}"] = {
                            "solidification": solid_value,
                            "snapshots": snapshot_value,
                        }
    return {
        "consistent": not mismatched and not missing,
        "mismatched_files": mismatched,
        "missing_files": missing,
        "parameter_differences": parameter_differences,
    }


def infer_grid_step(values: Iterable[float], fallback: float) -> float:
    unique = np.sort(np.unique(np.asarray(values, dtype=float)))
    if len(unique) < 2:
        return fallback
    differences = np.diff(unique)
    differences = differences[differences > max(fallback * 1e-4, 1e-12)]
    if not len(differences):
        return fallback
    return float(np.min(differences))


def interpolated_pool_depth(
    pool: pd.DataFrame,
    shoulder: pd.DataFrame,
    liquidus_k: float,
    dx: float,
    dy: float,
    dz: float,
    domain_z_min: float | None = None,
) -> tuple[float, bool, float]:
    """Interpolate maximum depth to the liquidus isotherm in each XY column.

    This avoids grid-quantized ``z.max() - z.min() + dz`` depths. The return
    values are depth, domain-bottom contact, and the interpolated-column ratio;
    a bottom-limited depth is only a lower bound.
    """
    if pool.empty:
        return 0.0, False, 0.0

    def index_of(values: pd.Series, origin: float, step: float) -> np.ndarray:
        return np.rint((values.to_numpy(float) - origin) / step).astype(np.int64)

    x0 = float(pool["x"].min())
    y0 = float(pool["y"].min())
    z0 = float(pool["z"].min())

    pool_ix = index_of(pool["x"], x0, dx)
    pool_iy = index_of(pool["y"], y0, dy)
    pool_iz = index_of(pool["z"], z0, dz)
    pool_t = pool["T"].to_numpy(float)

    # Use common integer indices to find the cold cell below each molten cell.
    below: dict[tuple[int, int, int], float] = {}
    if len(shoulder):
        s_ix = index_of(shoulder["x"], x0, dx)
        s_iy = index_of(shoulder["y"], y0, dy)
        s_iz = index_of(shoulder["z"], z0, dz)
        for key, value in zip(zip(s_ix, s_iy, s_iz), shoulder["T"].to_numpy(float)):
            below[key] = float(value)

    deepest: dict[tuple[int, int], tuple[int, float]] = {}
    for ix, iy, iz, temperature in zip(pool_ix, pool_iy, pool_iz, pool_t):
        key = (int(ix), int(iy))
        current = deepest.get(key)
        if current is None or iz < current[0]:
            deepest[key] = (int(iz), float(temperature))

    surface_z = float(pool["z"].max())

    crossings: list[float] = []
    bottom_limited = False
    interpolated = 0
    for (ix, iy), (iz, hot_t) in deepest.items():
        hot_z = z0 + iz * dz
        if domain_z_min is not None and hot_z <= domain_z_min + 0.5 * dz:
            crossings.append(hot_z)
            bottom_limited = True
            continue
        cold_t = below.get((ix, iy, iz - 1))
        if cold_t is None:
            # A missing shoulder cell implies a steep gradient; use a half-cell.
            crossings.append(hot_z - 0.5 * dz)
            continue
        span = hot_t - cold_t
        if abs(span) < 1e-12:
            crossings.append(hot_z - 0.5 * dz)
        else:
            ratio = float(np.clip((liquidus_k - cold_t) / span, 0.0, 1.0))
            crossings.append(hot_z - (1.0 - ratio) * dz)
            interpolated += 1

    if not crossings:
        return 0.0, False, 0.0
    depth = surface_z - min(crossings)
    return max(depth, 0.0), bottom_limited, interpolated / len(deepest)


def component_for_geometry(
    molten: pd.DataFrame,
    dx: float,
    dy: float,
    dz: float,
) -> tuple[pd.DataFrame, int]:
    """Return the current/largest 6-connected molten component."""
    if molten.empty:
        return molten, 0
    origin = molten[["x", "y", "z"]].min().to_numpy(float)
    steps = np.array([dx, dy, dz], dtype=float)
    integer = np.rint(
        (molten[["x", "y", "z"]].to_numpy(float) - origin) / steps
    ).astype(int)
    coordinates = [tuple(row) for row in integer]
    coordinate_rows: dict[tuple[int, int, int], list[int]] = {}
    for row_index, coordinate in enumerate(coordinates):
        coordinate_rows.setdefault(coordinate, []).append(row_index)
    unseen = set(coordinate_rows)
    components: list[list[tuple[int, int, int]]] = []
    while unseen:
        seed = unseen.pop()
        stack = [seed]
        component = [seed]
        while stack:
            current = stack.pop()
            for axis in range(3):
                for direction in (-1, 1):
                    neighbour = list(current)
                    neighbour[axis] += direction
                    neighbour_tuple = tuple(neighbour)
                    if neighbour_tuple in unseen:
                        unseen.remove(neighbour_tuple)
                        stack.append(neighbour_tuple)
                        component.append(neighbour_tuple)
        components.append(component)

    hottest_row = int(molten["T"].astype(float).idxmax())
    hottest_position = molten.index.get_loc(hottest_row)
    hottest_coordinate = coordinates[hottest_position]
    selected = next((part for part in components if hottest_coordinate in part), [])
    largest = max(components, key=len)
    if len(selected) < max(2, int(0.1 * len(largest))):
        selected = largest
    row_positions = [
        row_position
        for coordinate in selected
        for row_position in coordinate_rows[coordinate]
    ]
    return molten.iloc[row_positions].copy(), len(components)


def snapshot_geometry(
    file: Path,
    liquidus_k: float,
    resolution: tuple[float, float, float],
) -> dict[str, Any]:
    molten_parts: list[pd.DataFrame] = []
    shoulder_parts: list[pd.DataFrame] = []
    # Retain a temperature shoulder for liquidus interpolation without loading
    # the near-ambient bulk of the domain.
    shoulder_floor = 0.5 * liquidus_k
    row_count = 0
    finite_temperature = 0
    domain_z_min = float("inf")
    for chunk in pd.read_csv(file, chunksize=250_000):
        required = {"x", "y", "z", "T"}
        if not required.issubset(chunk.columns):
            raise ValueError(f"{file.name} is missing snapshot columns: {required}")
        row_count += len(chunk)
        z_values = pd.to_numeric(chunk["z"], errors="coerce")
        if np.isfinite(z_values).any():
            domain_z_min = min(domain_z_min, float(np.nanmin(z_values)))
        temperatures = pd.to_numeric(chunk["T"], errors="coerce")
        finite_temperature += int(np.isfinite(temperatures).sum())
        mask = temperatures >= liquidus_k
        if mask.any():
            molten_parts.append(chunk.loc[mask, ["x", "y", "z", "T"]])
        band = (temperatures >= shoulder_floor) & (~mask)
        if band.any():
            shoulder_parts.append(chunk.loc[band, ["x", "y", "z", "T"]])

    if not molten_parts:
        return {
            "file": file.name,
            "rows": row_count,
            "finite_fraction": finite_temperature / max(row_count, 1),
            "molten_voxels": 0,
            "components": 0,
            "width_m": None,
            "length_m": None,
            "depth_m": None,
            "molten_volume_m3": 0.0,
        }

    molten = pd.concat(molten_parts, ignore_index=True)
    # Snapshot and solidification cases may use different grids. Infer the
    # geometry grid from snapshot coordinates and use Domain values as fallback.
    dx = infer_grid_step(molten["x"], resolution[0] or 50e-6)
    dy = infer_grid_step(molten["y"], resolution[1] or 50e-6)
    dz = infer_grid_step(molten["z"], resolution[2] or 25e-6)
    current_pool, component_count = component_for_geometry(molten, dx, dy, dz)
    xy = current_pool[["x", "y"]].to_numpy(float)
    if len(xy) >= 3:
        centered = xy - np.mean(xy, axis=0)
        covariance = np.cov(centered, rowvar=False)
        _, eigenvectors = np.linalg.eigh(covariance)
        projected = centered @ eigenvectors
        extents = np.ptp(projected, axis=0) + max(dx, dy)
    else:
        extents = np.ptp(xy, axis=0) + np.array([dx, dy])
    width, length = sorted(float(value) for value in extents)
    shoulder = (
        pd.concat(shoulder_parts, ignore_index=True)
        if shoulder_parts
        else pd.DataFrame(columns=["x", "y", "z", "T"])
    )
    depth, bottom_limited, interpolated_columns = interpolated_pool_depth(
        current_pool, shoulder, liquidus_k, dx, dy, dz,
        domain_z_min if math.isfinite(domain_z_min) else None,
    )
    return {
        "file": file.name,
        "rows": row_count,
        "finite_fraction": finite_temperature / max(row_count, 1),
        "molten_voxels": int(len(molten)),
        "current_pool_voxels": int(len(current_pool)),
        "components": int(component_count),
        "width_m": width,
        "length_m": length,
        "depth_m": depth,
        "depth_bottom_limited": bool(bottom_limited),
        "depth_interpolated_column_fraction": float(interpolated_columns),
        "grid_x_m": float(dx),
        "grid_z_m": float(dz),
        "molten_volume_m3": float(len(current_pool) * dx * dy * dz),
        "max_temperature_k_diagnostic_only": float(molten["T"].max()),
    }


def analyse_thermal_histories(
    directory: Path | None,
    liquidus_k: float,
    maximum_files: int,
) -> dict[str, Any] | None:
    if directory is None:
        return None
    data_dir, _ = data_directory(directory)
    files = sorted(data_dir.glob("*T_hist*.csv"))
    if not files:
        return None
    selected = files
    sampled = False
    if len(files) > maximum_files:
        indices = np.linspace(0, len(files) - 1, maximum_files, dtype=int)
        selected = [files[index] for index in indices]
        sampled = True
    peak_temperatures: list[float] = []
    cooling_rates: list[float] = []
    cycle_counts: list[int] = []
    failures = 0
    for file in selected:
        try:
            frame = pd.read_csv(file)
            if not {"t", "T"}.issubset(frame.columns):
                failures += 1
                continue
            time = pd.to_numeric(frame["t"], errors="coerce").to_numpy(float)
            temperature = pd.to_numeric(frame["T"], errors="coerce").to_numpy(float)
            valid = np.isfinite(time) & np.isfinite(temperature)
            time, temperature = time[valid], temperature[valid]
            if len(time) < 2:
                failures += 1
                continue
            peak_temperatures.append(float(np.max(temperature)))
            delta_time = np.diff(time)
            rates = np.diff(temperature)[delta_time != 0] / delta_time[delta_time != 0]
            cooling_rates.append(float(abs(np.min(rates))) if len(rates) else 0.0)
            above = temperature >= liquidus_k
            cycles = int(np.sum(above & np.concatenate(([True], ~above[:-1]))))
            cycle_counts.append(cycles)
        except Exception:
            failures += 1
    return {
        "available_files": len(files),
        "analysed_files": len(selected) - failures,
        "sampled": sampled,
        "failed_files": failures,
        "PeakT": numeric_summary(pd.Series(peak_temperatures, dtype=float)),
        "CoolingRate": numeric_summary(pd.Series(cooling_rates, dtype=float)),
        "ThermalCycle": numeric_summary(pd.Series(cycle_counts, dtype=float)),
    }


def generate_domain_xy(
    parsed_values: dict[str, Any],
    grid_x_m: float,
    grid_y_m: float,
) -> np.ndarray:
    def value(suffix: str) -> float | None:
        found = find_key(parsed_values, suffix)
        return float(found) if found is not None else None

    x_min, x_max = value("X.Min"), value("X.Max")
    y_min, y_max = value("Y.Min"), value("Y.Max")
    if None in (x_min, x_max, y_min, y_max):
        return np.empty((0, 2), dtype=float)
    x_count = int(math.floor((x_max - x_min) / grid_x_m + 0.5)) + 1
    y_count = int(math.floor((y_max - y_min) / grid_y_m + 0.5)) + 1
    if x_count <= 0 or y_count <= 0 or x_count * y_count > 20_000_000:
        return np.empty((0, 2), dtype=float)
    x_values = x_min + np.arange(x_count) * grid_x_m
    y_values = y_min + np.arange(y_count) * grid_y_m
    xx, yy = np.meshgrid(x_values, y_values, indexing="ij")
    return np.column_stack((xx.ravel(), yy.ravel()))


def build_coverage_target_points(
    parsed_values: dict[str, Any],
    path_frame: pd.DataFrame | None,
    parameters: "ProcessParameters",
    source_cli: Path | None = None,
) -> pd.DataFrame:
    """Build a compact 3DThesis Custom point file at z=0 and z=-layer thickness."""
    xy = generate_domain_xy(
        parsed_values, parameters.grid_x_m, parameters.grid_y_m
    )
    segments = powered_segments(path_frame)
    corridor_half_width = 0.5 * max(
        parameters.hatch_spacing_m,
        parameters.grid_x_m,
        parameters.grid_y_m,
    )
    mask = target_region_mask(xy, segments, corridor_half_width)
    # Coverage includes powered contours; CLI hatch matching uses infill only.
    cli_geometry, cli_metadata = match_cli_solid_geometry(
        source_cli, infill_segments(path_frame)
    )
    if cli_geometry is not None:
        corridor_mask = mask.copy()
        solid_mask = cli_solid_region_mask(xy, cli_geometry)
        intersection_count = int((corridor_mask & solid_mask).sum())
        solid_count = int(solid_mask.sum())
        # When the source CLI is available, its complete solid cross-section is
        # the required denominator.  Using corridor∩solid can hide a solid island
        # whose scan path is missing — exactly the defect coverage is meant to find.
        mask = solid_mask
        cli_metadata.update(
            {
                "method": "complete matched CLI solid cross-section; path corridor is diagnostic only",
                "path_corridor_point_count": int(corridor_mask.sum()),
                "solid_target_point_count": solid_count,
                "path_corridor_and_solid_point_count": intersection_count,
                "solid_points_outside_path_corridor": solid_count - intersection_count,
                "path_corridor_coverage_of_solid_fraction": (
                    intersection_count / solid_count if solid_count else None
                ),
                "excluded_void_point_count": int(len(xy) - solid_count),
            }
        )
    target_xy = xy[mask]
    if not len(target_xy):
        empty = pd.DataFrame(columns=["x_mm", "y_mm", "z_mm"])
        empty.attrs["cli_solid_mask"] = cli_metadata
        return empty
    top = np.column_stack((target_xy, np.zeros(len(target_xy))))
    interface = np.column_stack(
        (
            target_xy,
            np.full(len(target_xy), -parameters.layer_thickness_m),
        )
    )
    points_mm = np.vstack((top, interface)) * 1e3
    result = pd.DataFrame(
        points_mm, columns=["x_mm", "y_mm", "z_mm"]
    ).drop_duplicates()
    result.attrs["cli_solid_mask"] = cli_metadata
    return result


def target_plane_num_melt(
    solid: pd.DataFrame,
    target_xy: np.ndarray,
    target_mask: np.ndarray,
    target_z_m: float,
    parameters: "ProcessParameters",
    segments: np.ndarray | None = None,
    corridor_half_width_m: float | None = None,
) -> dict[str, Any]:
    target = target_xy[target_mask]
    result: dict[str, Any] = {
        "requested_z_m": target_z_m,
        "target_point_count": int(len(target)),
        "available": False,
        "coordinate_completeness_fraction": 0.0,
        "melted_fraction_with_missing_as_unmelted": None,
    }
    if not len(target) or "numMelt" not in solid.columns:
        return result
    z_values = pd.to_numeric(solid["z"], errors="coerce").dropna().unique()
    if not len(z_values):
        return result
    selected_z = float(z_values[np.argmin(np.abs(z_values - target_z_m))])
    result["selected_z_m"] = selected_z
    result["z_offset_m"] = abs(selected_z - target_z_m)
    if abs(selected_z - target_z_m) > 0.51 * parameters.grid_z_m:
        return result

    plane = solid[
        np.isclose(
            pd.to_numeric(solid["z"], errors="coerce"),
            selected_z,
            atol=max(parameters.grid_z_m * 0.1, 1e-12),
        )
    ].copy()
    if plane.empty:
        return result
    target_keys = pd.DataFrame(
        {
            "ix": np.rint(target[:, 0] / parameters.grid_x_m).astype(np.int64),
            "iy": np.rint(target[:, 1] / parameters.grid_y_m).astype(np.int64),
            "x": target[:, 0],
            "y": target[:, 1],
        }
    ).drop_duplicates(subset=["ix", "iy"])
    plane["ix"] = np.rint(
        pd.to_numeric(plane["x"], errors="coerce") / parameters.grid_x_m
    ).astype("Int64")
    plane["iy"] = np.rint(
        pd.to_numeric(plane["y"], errors="coerce") / parameters.grid_y_m
    ).astype("Int64")
    plane["numMelt_numeric"] = pd.to_numeric(plane["numMelt"], errors="coerce")
    reported = (
        plane.dropna(subset=["ix", "iy", "numMelt_numeric"])
        .groupby(["ix", "iy"], as_index=False)["numMelt_numeric"]
        .max()
    )
    joined = target_keys.merge(reported, on=["ix", "iy"], how="left")
    present = joined["numMelt_numeric"].notna()
    result.update(
        {
            "available": True,
            "target_point_count": int(len(joined)),
            "reported_target_point_count": int(present.sum()),
            "coordinate_completeness_fraction": float(present.mean()),
            "melted_fraction_with_missing_as_unmelted": float(
                joined["numMelt_numeric"].fillna(0).ge(1).mean()
            ),
        }
    )
    result["spatial_diagnostics"] = coverage_spatial_diagnostics(
        joined,
        segments if segments is not None else np.empty((0, 4), dtype=float),
        corridor_half_width_m,
    )
    return result


def coverage_spatial_diagnostics(
    joined: pd.DataFrame,
    segments: np.ndarray,
    corridor_half_width_m: float | None,
) -> dict[str, Any]:
    """Explain whether unmelted target points occur at edges, ends, or track interiors."""
    if joined.empty:
        return {}
    melted = joined["numMelt_numeric"].fillna(0).ge(1).to_numpy(bool)
    keys = list(zip(joined["ix"].astype(int), joined["iy"].astype(int)))
    key_set = set(keys)
    boundary = np.asarray(
        [
            any(
                neighbour not in key_set
                for neighbour in (
                    (ix + 1, iy),
                    (ix - 1, iy),
                    (ix, iy + 1),
                    (ix, iy - 1),
                )
            )
            for ix, iy in keys
        ],
        dtype=bool,
    )
    unseen = {
        key for key, is_melted in zip(keys, melted) if not is_melted
    }
    component_sizes: list[int] = []
    while unseen:
        seed = unseen.pop()
        stack = [seed]
        size = 1
        while stack:
            ix, iy = stack.pop()
            for neighbour in (
                (ix + 1, iy),
                (ix - 1, iy),
                (ix, iy + 1),
                (ix, iy - 1),
            ):
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    stack.append(neighbour)
                    size += 1
        component_sizes.append(size)

    result: dict[str, Any] = {
        "target_count": int(len(joined)),
        "melted_count": int(melted.sum()),
        "unmelted_count": int((~melted).sum()),
        "boundary_target_count": int(boundary.sum()),
        "boundary_unmelted_count": int((boundary & ~melted).sum()),
        "interior_unmelted_count": int((~boundary & ~melted).sum()),
        "unmelted_cluster_count": len(component_sizes),
        "largest_unmelted_cluster_count": (
            max(component_sizes) if component_sizes else 0
        ),
    }
    if not len(segments) or corridor_half_width_m is None:
        return result

    xy = joined[["x", "y"]].to_numpy(float)
    nearest_distance = np.full(len(xy), np.inf, dtype=float)
    nearest_is_endpoint = np.zeros(len(xy), dtype=bool)
    for x0, y0, x1, y1 in segments:
        vector = np.array([x1 - x0, y1 - y0], dtype=float)
        length_squared = float(vector @ vector)
        if length_squared <= 0:
            continue
        relative = xy - np.array([x0, y0])
        raw_fraction = (
            relative[:, 0] * vector[0] + relative[:, 1] * vector[1]
        ) / length_squared
        fraction = np.clip(raw_fraction, 0.0, 1.0)
        closest = np.array([x0, y0]) + fraction[:, None] * vector
        distance = np.linalg.norm(xy - closest, axis=1)
        better = distance < nearest_distance
        nearest_distance[better] = distance[better]
        nearest_is_endpoint[better] = (
            (raw_fraction[better] <= 0.0) | (raw_fraction[better] >= 1.0)
        )

    track_interior = ~nearest_is_endpoint
    result.update(
        {
            "endpoint_target_count": int(nearest_is_endpoint.sum()),
            "endpoint_melted_count": int((nearest_is_endpoint & melted).sum()),
            "endpoint_melted_fraction": (
                float(melted[nearest_is_endpoint].mean())
                if nearest_is_endpoint.any()
                else None
            ),
            "track_interior_target_count": int(track_interior.sum()),
            "track_interior_melted_count": int((track_interior & melted).sum()),
            "track_interior_melted_fraction": (
                float(melted[track_interior].mean())
                if track_interior.any()
                else None
            ),
        }
    )
    radius = float(corridor_half_width_m)
    band_edges = (0.0, 0.25, 0.50, 0.75, 1.000001, np.inf)
    bands: list[dict[str, Any]] = []
    normalized_distance = nearest_distance / max(radius, 1e-15)
    for low, high in zip(band_edges, band_edges[1:]):
        selected = (normalized_distance >= low) & (normalized_distance < high)
        if not selected.any():
            continue
        bands.append(
            {
                "distance_fraction_low": low,
                "distance_fraction_high": None if np.isinf(high) else high,
                "distance_um_low": low * radius * 1e6,
                "distance_um_high": (
                    None if np.isinf(high) else high * radius * 1e6
                ),
                "target_count": int(selected.sum()),
                "melted_count": int((selected & melted).sum()),
                "melted_fraction": float(melted[selected].mean()),
            }
        )
    result["distance_bands"] = bands
    return result


def build_problem_diagnosis(
    coverage: dict[str, Any],
    metrics: dict[str, Any],
    component_scores: dict[str, float | None],
    parameters: "ProcessParameters",
    snapshot_count: int,
    minimum_geometry_snapshots: int,
    hard_fail_fraction: float,
    pass_fraction: float,
) -> list[dict[str, Any]]:
    """Turn numerical findings into bounded, evidence-linked engineering hypotheses."""
    issues: list[dict[str, Any]] = []

    def add(
        severity: str,
        category: str,
        where: str,
        finding: str,
        evidence: str,
        likely_causes: list[str],
        actions: list[str],
        confidence: str,
    ) -> None:
        issues.append(
            {
                "severity": severity,
                "category": category,
                "where": where,
                "finding": finding,
                "evidence": evidence,
                "likely_causes": likely_causes,
                "recommended_actions": actions,
                "confidence": confidence,
            }
        )

    planes = coverage.get("num_melt_target_grid", {}).get("planes", {})
    for key, label, category in (
        ("previous_layer_interface", "Previous-layer interface", "Interlayer fusion"),
        ("top_surface", "Current-layer top surface", "Scan coverage"),
    ):
        plane = planes.get(key, {})
        fraction = plane.get("melted_fraction_with_missing_as_unmelted")
        if fraction is None or fraction >= pass_fraction:
            continue
        spatial = plane.get("spatial_diagnostics", {})
        unmelted = spatial.get("unmelted_count", 0)
        interior_unmelted = spatial.get("interior_unmelted_count", 0)
        boundary_unmelted = spatial.get("boundary_unmelted_count", 0)
        track_fraction = spatial.get("track_interior_melted_fraction")
        endpoint_fraction = spatial.get("endpoint_melted_fraction")
        bands = spatial.get("distance_bands", [])
        centre_fraction = bands[0].get("melted_fraction") if bands else None
        outer_fraction = bands[-1].get("melted_fraction") if bands else None
        causes: list[str] = []
        actions: list[str] = []
        if key == "previous_layer_interface":
            if (
                centre_fraction is not None
                and outer_fraction is not None
                and centre_fraction - outer_fraction > 0.30
            ):
                causes.append(
                    "The melt pool narrows substantially at layer-depth, leaving no "
                    "continuous overlap between adjacent tracks"
                )
            if interior_unmelted > boundary_unmelted:
                causes.append(
                    "Lack of fusion is concentrated inside the target region, not only "
                    "at a potentially imperfect mask boundary"
                )
            causes.extend(
                [
                    "Layer thickness is large relative to effective penetration depth",
                    "Hatch spacing is large relative to effective melt width at the interface",
                ]
            )
            actions.extend(
                [
                    "Test a smaller layer thickness and hatch spacing together before "
                    "increasing power further",
                    "Refine the Z and XY grids and export complete interface numMelt data",
                ]
            )
        else:
            if unmelted and boundary_unmelted / unmelted >= 0.80:
                causes.append(
                    "Unmelted points are concentrated at the buffered path boundary; "
                    "the target mask may be conservative"
                )
                actions.append(
                    "Intersect the scan-path target region with the actual CLI layer section"
                )
            if (
                track_fraction is not None
                and endpoint_fraction is not None
                and track_fraction - endpoint_fraction > 0.20
            ):
                causes.append("Coverage at scan endpoints is substantially below track interiors")
                actions.append(
                    "Check lead-in/lead-out segments, endpoint power, and target end-cap definitions"
                )
            if not causes:
                causes.append("A local scan interruption or target-region mismatch is present")
        evidence = (
            f"Coverage {100 * float(fraction):.2f}%; {unmelted} unmelted points "
            f"({interior_unmelted} interior, {boundary_unmelted} boundary)"
        )
        add(
            "CRITICAL" if float(fraction) < hard_fail_fraction else "WARNING",
            category,
            label,
            (
                "Authoritative numMelt data show a large area that never melted"
                if float(fraction) < hard_fail_fraction
                else "Coverage is close to, but below, the release target"
            ),
            evidence,
            causes,
            actions,
            "HIGH",
        )

    if metrics.get("lof_resolution_sensitive"):
        add(
            "WARNING",
            "Numerical resolution",
            "Melt-pool width/depth and LOF",
            "The LOF result crosses the acceptance boundary under grid uncertainty",
            (
                f"Nominal/half-cell/conservative indices are "
                f"{metrics.get('lof_index_nominal', float('nan')):.3f}/"
                f"{metrics.get('lof_index_midpoint', float('nan')):.3f}/"
                f"{metrics.get('lof_index_conservative', float('nan')):.3f}"
            ),
            [
                "The current grid is coarse relative to the melt-pool dimensions",
                "The liquidus-isotherm dimensions have about one cell of uncertainty",
            ],
            [
                "Rerun with a finer grid",
                "Treat LOF only as an interval diagnosis until the refined rerun",
            ],
            "HIGH",
        )

    if snapshot_count < minimum_geometry_snapshots:
        add(
            "WARNING",
            "Evidence adequacy",
            "Melt-pool temporal sampling",
            "The snapshots support a dimension estimate but do not robustly represent the full path",
            f"{snapshot_count} snapshots; the geometry-statistics target is "
            f"{minimum_geometry_snapshots}",
            ["A small sample can miss dimensional changes between tracks, at corners, and at endpoints"],
            ["Retain 20-50 snapshots distributed uniformly over scan progress"],
            "HIGH",
        )

    peak_temperature = metrics.get(
        "snapshot_peak_temperature_max_k_diagnostic_only"
    )
    if (
        peak_temperature is not None
        and peak_temperature > 3.0 * parameters.liquidus_k
    ):
        add(
            "WARNING",
            "Model applicability",
            "Melt-pool hotspot",
            "Peak temperature is far above liquidus; a conduction-only model cannot "
            "reliably infer physical keyholing or evaporation",
            (
                f"Snapshot peak {peak_temperature:.0f} K; "
                f"liquidus {parameters.liquidus_k:.0f} K"
            ),
            [
                "Local energy input is very high",
                "The model does not resolve flow, recoil pressure, or energy loss by evaporation",
            ],
            [
                "Verify beam power, efficiency, and units",
                "Reassess keyhole risk experimentally or with a fluid-flow model",
            ],
            "MEDIUM",
        )

    thermal_score = component_scores.get("thermal_uniformity")
    if thermal_score is not None and thermal_score < 70.0:
        add(
            "ADVISORY",
            "Thermal-field consistency",
            "Melted region",
            "Solidification thermal parameters show moderate variation across the layer",
            f"Thermal-field consistency score {thermal_score:.1f}/100",
            ["Path endpoints, corners, or heat accumulation alter G, V, and cooling rate"],
            [
                "Use Coverage_Diagnostics.csv to locate low-coverage regions before "
                "inspecting their local thermal histories"
            ],
            "MEDIUM",
        )
    return issues


def analyse_layer_coverage(
    solid: pd.DataFrame,
    solid_files: list[Path],
    snapshot_files: list[Path],
    snapshot_features: pd.DataFrame,
    parameters: "ProcessParameters",
    path_frame: pd.DataFrame | None,
    parsed_values: dict[str, Any],
    solid_root: Path,
    snapshot_root: Path | None,
    config: dict[str, Any],
    source_cli: Path | None = None,
) -> dict[str, Any]:
    """Assess cumulative target coverage without mistaking snapshots for Tmax."""
    coverage_config = config.get("coverage", {})
    segments = powered_segments(path_frame)
    segment_lengths = (
        np.linalg.norm(segments[:, 2:4] - segments[:, 0:2], axis=1)
        if len(segments)
        else np.array([], dtype=float)
    )
    total_scan_length = float(np.sum(segment_lengths))
    corridor_half_width = 0.5 * max(
        parameters.hatch_spacing_m,
        parameters.grid_x_m,
        parameters.grid_y_m,
    )
    snapshot_provenance = case_provenance(snapshot_root, snapshot_files)
    solid_provenance = case_provenance(solid_root, solid_files)
    cross_case = cross_case_compatibility(solid_root, snapshot_root)
    snapshot_comparable = cross_case.get("consistent") is True
    coverage_snapshot_files = snapshot_files if snapshot_comparable else []
    result: dict[str, Any] = {
        "status": "INSUFFICIENT_EVIDENCE",
        "authoritative": False,
        "method": None,
        "coverage_fraction": None,
        "minimum_required_fraction": float(
            coverage_config.get("minimum_coverage_fraction", 0.99)
        ),
        "target_definition": {
            "method": "powered-path corridor fallback; complete CLI solid cross-section when matched",
            "powered_segment_count": int(len(segments)),
            "total_scan_length_m": total_scan_length,
            "corridor_half_width_m": corridor_half_width,
        },
        "snapshot_cumulative_temperature": {
            "available": False,
            "authoritative": False,
            "interpretation": "lower bound unless sampling is proven dense",
            "file_count": len(snapshot_files),
            "comparable_to_solidification": snapshot_comparable,
        },
        "num_melt_target_grid": {
            "available": "numMelt" in solid.columns,
            "authoritative": False,
        },
        "provenance": {
            "solidification": solid_provenance,
            "snapshots": snapshot_provenance,
            "cross_case_inputs": cross_case,
        },
    }

    reference_xy = generate_domain_xy(
        parsed_values, parameters.grid_x_m, parameters.grid_y_m
    )
    target_mask = target_region_mask(reference_xy, segments, corridor_half_width)
    # Coverage includes powered contours; CLI hatch matching uses infill only.
    cli_geometry, cli_metadata = match_cli_solid_geometry(
        source_cli, infill_segments(path_frame)
    )
    if cli_geometry is not None:
        corridor_mask = target_mask.copy()
        solid_mask = cli_solid_region_mask(reference_xy, cli_geometry)
        intersection_count = int((corridor_mask & solid_mask).sum())
        solid_count = int(solid_mask.sum())
        target_mask = solid_mask
        cli_metadata.update(
            {
                "method": "complete matched CLI solid cross-section; path corridor is diagnostic only",
                "path_corridor_point_count": int(corridor_mask.sum()),
                "solid_target_point_count": solid_count,
                "path_corridor_and_solid_point_count": intersection_count,
                "solid_points_outside_path_corridor": solid_count - intersection_count,
                "path_corridor_coverage_of_solid_fraction": (
                    intersection_count / solid_count if solid_count else None
                ),
                "excluded_void_point_count": int(len(reference_xy) - solid_count),
            }
        )
    result["target_definition"]["cli_solid_mask"] = cli_metadata
    snapshot_section = result["snapshot_cumulative_temperature"]
    coordinate_grid: np.ndarray | None = None
    maximum_temperature: np.ndarray | None = None
    grid_consistent = True
    read_error: str | None = None
    for file in coverage_snapshot_files:
        try:
            frame = pd.read_csv(file, usecols=["x", "y", "z", "T"])
            coordinates = frame[["x", "y", "z"]].apply(
                pd.to_numeric, errors="coerce"
            ).to_numpy(float)
            temperature = pd.to_numeric(frame["T"], errors="coerce").to_numpy(float)
            if coordinate_grid is None:
                coordinate_grid = coordinates
                maximum_temperature = temperature
            else:
                same_grid = len(coordinates) == len(coordinate_grid) and np.allclose(
                    coordinates, coordinate_grid, rtol=0.0, atol=1e-12, equal_nan=True
                )
                if not same_grid:
                    grid_consistent = False
                    break
                maximum_temperature = np.fmax(maximum_temperature, temperature)
        except Exception as error:
            read_error = f"{file.name}: {error}"
            grid_consistent = False
            break

    snapshot_planes: dict[str, Any] = {}
    if (
        coordinate_grid is not None
        and maximum_temperature is not None
        and grid_consistent
        and len(segments)
    ):
        z_values = np.unique(coordinate_grid[:, 2][np.isfinite(coordinate_grid[:, 2])])
        for label, requested_z in (
            ("top_surface", 0.0),
            ("previous_layer_interface", -parameters.layer_thickness_m),
        ):
            if not len(z_values):
                continue
            selected_z = float(z_values[np.argmin(np.abs(z_values - requested_z))])
            plane_mask = np.isclose(
                coordinate_grid[:, 2],
                selected_z,
                atol=max(parameters.grid_z_m * 0.1, 1e-12),
            )
            xy = coordinate_grid[plane_mask, :2]
            temperatures = maximum_temperature[plane_mask]
            plane_target = target_region_mask(xy, segments, corridor_half_width)
            if cli_geometry is not None:
                plane_target = cli_solid_region_mask(xy, cli_geometry)
            molten = np.isfinite(temperatures) & (
                temperatures >= parameters.liquidus_k
            )
            snapshot_planes[label] = {
                "requested_z_m": requested_z,
                "selected_z_m": selected_z,
                "z_offset_m": abs(selected_z - requested_z),
                "target_point_count": int(plane_target.sum()),
                "observed_molten_fraction_lower_bound": (
                    float(molten[plane_target].mean())
                    if plane_target.any()
                    else None
                ),
            }
            if label == "top_surface":
                reference_xy = xy
                target_mask = plane_target

        pool_lengths = (
            finite_positive(snapshot_features["length_m"])
            if "length_m" in snapshot_features.columns
            else np.array([], dtype=float)
        )
        median_pool_length = float(np.median(pool_lengths)) if len(pool_lengths) else None
        estimated_minimum = (
            int(math.ceil(total_scan_length / median_pool_length))
            if median_pool_length and total_scan_length > 0
            else None
        )
        scan_fractions = scan_fraction_metadata(snapshot_root)
        metadata_matches = (
            scan_fractions["count"] is not None
            and scan_fractions["count"] == len(snapshot_files)
        )
        enough_by_estimate = (
            estimated_minimum is not None and len(snapshot_files) >= estimated_minimum
        )
        allow_dense = bool(
            coverage_config.get("allow_dense_snapshots_as_authoritative", False)
        )
        snapshot_authoritative = bool(
            allow_dense
            and enough_by_estimate
            and metadata_matches
            and snapshot_provenance.get("consistent") is True
            and all(
                plane.get("z_offset_m", np.inf) <= 0.51 * parameters.grid_z_m
                for plane in snapshot_planes.values()
            )
            and len(snapshot_planes) == 2
        )
        lower_bounds = [
            plane["observed_molten_fraction_lower_bound"]
            for plane in snapshot_planes.values()
            if plane.get("observed_molten_fraction_lower_bound") is not None
        ]
        snapshot_section.update(
            {
                "available": True,
                "authoritative": snapshot_authoritative,
                "grid_consistent": True,
                "planes": snapshot_planes,
                "median_observed_pool_length_m": median_pool_length,
                "estimated_minimum_snapshot_count": estimated_minimum,
                "temporal_sampling_ratio": (
                    min(len(snapshot_files) / estimated_minimum, 1.0)
                    if estimated_minimum
                    else None
                ),
                "scan_fraction_metadata": scan_fractions,
                "metadata_file_count_matches": metadata_matches,
                "coverage_fraction_if_dense": min(lower_bounds)
                if lower_bounds
                else None,
            }
        )
    else:
        snapshot_section.update(
            {
                "grid_consistent": grid_consistent,
                "read_error": read_error,
                "reason": (
                    "snapshot physics/path inputs differ from the solidification case"
                    if not snapshot_comparable and snapshot_files
                    else "no powered path/grid/snapshot data"
                    if read_error is None
                    else "snapshot accumulation failed"
                ),
                "scan_fraction_metadata": scan_fraction_metadata(snapshot_root),
                "metadata_file_count_matches": (
                    scan_fraction_metadata(snapshot_root).get("count")
                    == len(snapshot_files)
                    if snapshot_files
                    else None
                ),
            }
        )

    tracking = find_key(parse_grouped_text(solid_root / "Mode.txt"), "Solidification.Tracking")
    if tracking is None:
        tracking = find_key(parse_grouped_text(solid_root / "Mode.txt"), "Tracking")
    num_melt_section = result["num_melt_target_grid"]
    num_melt_section["tracking_mode"] = str(tracking) if tracking is not None else None
    if len(reference_xy) and len(target_mask):
        num_melt_planes = {
            "top_surface": target_plane_num_melt(
                solid,
                reference_xy,
                target_mask,
                0.0,
                parameters,
                segments,
                corridor_half_width,
            ),
            "previous_layer_interface": target_plane_num_melt(
                solid,
                reference_xy,
                target_mask,
                -parameters.layer_thickness_m,
                parameters,
                segments,
                corridor_half_width,
            ),
        }
        completeness = [
            plane["coordinate_completeness_fraction"]
            for plane in num_melt_planes.values()
            if plane.get("available")
        ]
        melted = [
            plane["melted_fraction_with_missing_as_unmelted"]
            for plane in num_melt_planes.values()
            if plane.get("melted_fraction_with_missing_as_unmelted") is not None
        ]
        accepted_modes = {
            str(value).strip().lower()
            for value in coverage_config.get(
                "authoritative_num_melt_tracking_modes", ["None"]
            )
        }
        custom_domain_file = find_key(parsed_values, "Custom.File")
        complete_grid = (
            len(completeness) == 2
            and min(completeness)
            >= float(coverage_config.get("minimum_coordinate_completeness", 0.99))
        )
        num_melt_authoritative = bool(
            "numMelt" in solid.columns
            and tracking is not None
            and (
                str(tracking).strip().lower() in accepted_modes
                or custom_domain_file is not None
            )
            and complete_grid
            and solid_provenance.get("consistent") is True
        )
        num_melt_section.update(
            {
                "planes": num_melt_planes,
                "coordinate_grid_complete": complete_grid,
                "custom_domain_file": custom_domain_file,
                "authoritative": num_melt_authoritative,
                "coverage_fraction": min(melted) if melted else None,
                "reported_num_melt_zero_fraction": float(
                    pd.to_numeric(solid["numMelt"], errors="coerce").fillna(0).eq(0).mean()
                )
                if "numMelt" in solid.columns
                else None,
            }
        )

    if num_melt_section.get("authoritative"):
        result.update(
            {
                "status": "COMPLETE_EVIDENCE",
                "authoritative": True,
                "method": "complete target-grid numMelt",
                "coverage_fraction": num_melt_section.get("coverage_fraction"),
            }
        )
    elif snapshot_section.get("authoritative"):
        result.update(
            {
                "status": "COMPLETE_EVIDENCE",
                "authoritative": True,
                "method": "dense cumulative liquidus snapshots",
                "coverage_fraction": snapshot_section.get(
                    "coverage_fraction_if_dense"
                ),
            }
        )
    elif snapshot_section.get("available"):
        result.update(
            {
                "status": "LOWER_BOUND_ONLY",
                "method": "sparse snapshot Tmax lower bound",
            }
        )
    return result


@dataclass
class ProcessParameters:
    liquidus_k: float
    liquidus_source: str
    layer_thickness_m: float
    layer_thickness_source: str
    hatch_spacing_m: float
    hatch_spacing_source: str
    grid_x_m: float
    grid_y_m: float
    grid_z_m: float
    power_w: float | None
    efficiency: float | None
    scan_speed_m_s: float | None


class LayerAssessment:
    def __init__(
        self,
        solidification_path: Path,
        snapshots_path: Path | None,
        output_dir: Path,
        config_path: Path,
        layer_id: str,
        thermal_history_path: Path | None = None,
        source_cli: Path | None = None,
        layer_thickness_m: float | None = None,
        hatch_spacing_m: float | None = None,
        stress_test: bool = False,
    ) -> None:
        self.solidification_path = solidification_path
        self.snapshots_path = snapshots_path
        self.output_dir = output_dir
        self.config_path = config_path
        self.layer_id = layer_id
        self.thermal_history_path = thermal_history_path
        self.source_cli = source_cli
        self.explicit_layer_thickness_m = layer_thickness_m
        self.explicit_hatch_spacing_m = hatch_spacing_m
        self.stress_test = stress_test
        self.config = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}

    def process_parameters(
        self,
        case_root: Path,
        snapshot_root: Path | None,
    ) -> tuple[ProcessParameters, dict[str, Any]]:
        roots = [case_root]
        if snapshot_root is not None and snapshot_root != case_root:
            roots.append(snapshot_root)
        parsed: dict[str, Any] = {}
        sources: dict[str, str] = {}
        for root in roots:
            for filename in ("Material.txt", "Domain.txt", "Beam.txt"):
                values = parse_grouped_text(root / filename)
                for key, value in values.items():
                    parsed.setdefault(key, value)
                    sources.setdefault(key, str(root / filename))

        override = self.config.get("process", {})
        fallback = self.config.get("fallback", {})

        def choose(
            override_key: str,
            parsed_suffix: str,
            fallback_key: str,
        ) -> tuple[float, str]:
            if override.get(override_key) is not None:
                return float(override[override_key]), f"config.process.{override_key}"
            value = find_key(parsed, parsed_suffix)
            if value is not None:
                matching = next(
                    key
                    for key in parsed
                    if key.lower() == parsed_suffix.lower()
                    or key.lower().endswith("." + parsed_suffix.lower())
                )
                return float(value), sources.get(matching, parsed_suffix)
            return float(fallback[fallback_key]), f"fallback.{fallback_key}"

        liquidus, liquidus_source = choose("liquidus_k", "Constants.T_L", "liquidus_k")

        source_cli = self.source_cli or discover_source_cli(roots)
        cli_layer = parse_cli_layer_thickness(source_cli)
        if self.explicit_layer_thickness_m is not None:
            layer_thickness = float(self.explicit_layer_thickness_m)
            layer_source = "--layer-thickness-um"
        elif override.get("layer_thickness_m") is not None:
            layer_thickness = float(override["layer_thickness_m"])
            layer_source = "config.process.layer_thickness_m"
        elif cli_layer is not None:
            layer_thickness = float(cli_layer["layer_thickness_m"])
            layer_source = str(cli_layer["path"])
        else:
            layer_thickness = float(fallback["layer_thickness_m"])
            layer_source = "fallback.layer_thickness_m"

        path_file = next((root / "Path.txt" for root in roots if (root / "Path.txt").is_file()), None)
        path_frame = parse_path_file(path_file) if path_file else None
        path_diagnostics = (
            path_frame.attrs.get("path_parse_diagnostics", {})
            if path_frame is not None
            else {}
        )
        inferred_hatch = (
            None
            if path_diagnostics.get("extreme_concatenated_line_format_detected")
            else infer_hatch_spacing(path_frame)
        )
        if self.explicit_hatch_spacing_m is not None:
            hatch = float(self.explicit_hatch_spacing_m)
            hatch_source = "--hatch-spacing-um"
        elif override.get("hatch_spacing_m") is not None:
            hatch = float(override["hatch_spacing_m"])
            hatch_source = "config.process.hatch_spacing_m"
        elif inferred_hatch is not None:
            hatch = inferred_hatch
            hatch_source = str(path_file)
        else:
            hatch = float(fallback["hatch_spacing_m"])
            hatch_source = "fallback.hatch_spacing_m"

        grid_x, _ = choose("grid_x_m", "X.Res", "grid_x_m")
        grid_y, _ = choose("grid_y_m", "Y.Res", "grid_y_m")
        grid_z, _ = choose("grid_z_m", "Z.Res", "grid_z_m")
        power = find_key(parsed, "Intensity.Power")
        efficiency = find_key(parsed, "Intensity.Efficiency")
        parameters = ProcessParameters(
            liquidus_k=liquidus,
            liquidus_source=liquidus_source,
            layer_thickness_m=layer_thickness,
            layer_thickness_source=layer_source,
            hatch_spacing_m=hatch,
            hatch_spacing_source=hatch_source,
            grid_x_m=grid_x,
            grid_y_m=grid_y,
            grid_z_m=grid_z,
            power_w=float(power) if power is not None else None,
            efficiency=float(efficiency) if efficiency is not None else None,
            scan_speed_m_s=scan_speed(path_frame),
        )
        return parameters, {
            "parsed_values": parsed,
            "path_file": str(path_file) if path_file else None,
            "path_parse_diagnostics": (
                path_diagnostics
                if path_frame is not None
                else {"error": "Path.txt could not be parsed"}
            ),
            "source_cli": str(source_cli) if source_cli else None,
            "cli_layer_analysis": cli_layer,
        }

    def run(self) -> dict[str, Any]:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        solid_files, case_root = find_solidification_files(self.solidification_path)
        snapshot_files, snapshot_root = find_snapshot_files(self.snapshots_path)
        rdf_files = find_rdf_files(case_root)
        parameters, raw_config = self.process_parameters(case_root, snapshot_root)
        solid = load_csvs(solid_files)
        if solid.empty:
            raise ValueError("Solidification data is empty")

        required_coordinates = {"x", "y", "z"}
        missing_coordinates = sorted(required_coordinates - set(solid.columns))
        if missing_coordinates:
            raise ValueError(f"Missing coordinate columns: {missing_coordinates}")

        summaries = {
            column: numeric_summary(solid[column])
            for column in solid.columns
            if pd.api.types.is_numeric_dtype(solid[column])
        }
        snapshot_rows = [
            snapshot_geometry(
                file,
                parameters.liquidus_k,
                (parameters.grid_x_m, parameters.grid_y_m, parameters.grid_z_m),
            )
            for file in snapshot_files
        ]
        snapshot_frame = pd.DataFrame(snapshot_rows)

        histories = analyse_thermal_histories(
            self.thermal_history_path,
            parameters.liquidus_k,
            int(self.config.get("thermal_history", {}).get("maximum_files", 5000)),
        )
        path_frame = (
            parse_path_file(Path(raw_config["path_file"]))
            if raw_config.get("path_file")
            else None
        )
        coverage = analyse_layer_coverage(
            solid=solid,
            solid_files=solid_files,
            snapshot_files=snapshot_files,
            snapshot_features=snapshot_frame,
            parameters=parameters,
            path_frame=path_frame,
            parsed_values=raw_config["parsed_values"],
            solid_root=case_root,
            snapshot_root=snapshot_root,
            config=self.config,
            source_cli=(
                Path(raw_config["source_cli"])
                if raw_config.get("source_cli")
                else None
            ),
        )
        target_points = build_coverage_target_points(
            raw_config["parsed_values"],
            path_frame,
            parameters,
            source_cli=(
                Path(raw_config["source_cli"])
                if raw_config.get("source_cli")
                else None
            ),
        )
        if not target_points.empty:
            target_path = self.output_dir / "coverage_target_points.txt"
            target_points.to_csv(
                target_path,
                sep=" ",
                header=False,
                index=False,
                float_format="%.9f",
            )
            coverage["recommended_custom_point_file"] = str(target_path.resolve())
            coverage["recommended_custom_point_count"] = int(len(target_points))

        result = self.score(
            solid,
            snapshot_frame,
            parameters,
            solid_files,
            snapshot_files,
            summaries,
            histories,
            coverage,
        )
        result["raw_input_config"] = raw_config
        path_diagnostics = raw_config.get("path_parse_diagnostics", {})
        result["input"]["path_parse_diagnostics"] = path_diagnostics
        if path_diagnostics.get("extreme_concatenated_line_format_detected"):
            result["flags"].append(
                (
                    "STRESS-TEST: Path.txt contains trailing concatenated tokens; "
                    "the first record on each physical line was used to mirror 3DThesis"
                    if self.stress_test
                    else "INFO: Path.txt contains trailing concatenated tokens; "
                    "the first record on each physical line was used to mirror 3DThesis"
                )
            )
        result["rdf"] = analyse_rdf(rdf_files)
        result["output_audit"] = audit_requested_outputs(
            case_root, solid.columns, rdf_files
        )
        result["input"]["rdf_files"] = [str(file) for file in rdf_files]
        solid_data_dir = (
            case_root / "Data" if (case_root / "Data").is_dir() else case_root
        )
        result["input"]["ignored_intermediate_solidification_files"] = [
            str(file)
            for file in sorted(solid_data_dir.glob("*.Solidification.*.csv"))
            if ".Solidification.Final" not in file.name
        ]
        if result["output_audit"].get("secondary_requested_but_disabled"):
            result["data_gaps"].append(
                {
                    "severity": "optional",
                    "item": "secondary solidification H fields",
                    "why": (
                        "Output.txt requests H/Hx/Hy/Hz but Mode.txt has Secondary 0, "
                        "so the fields are correctly absent"
                    ),
                    "required_3dthesis_output": (
                        "leave them disabled for normal scoring; only set Secondary 1 "
                        "for an explicit experimental study"
                    ),
                }
            )
        self.write_outputs(result, snapshot_frame, summaries)
        return result

    def geometry_values(
        self,
        solid: pd.DataFrame,
        snapshots: pd.DataFrame,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
        def column_values(names: tuple[str, ...]) -> np.ndarray:
            for name in names:
                if name in solid.columns:
                    return finite_positive(solid[name])
            return np.array([], dtype=float)

        direct_width = column_values(
            ("MP_width", "MP_Width", "MPWidth", "melt_pool_width")
        )
        direct_depth = column_values(
            ("MP_depth", "MP_Depth", "MPDepth", "melt_pool_depth", "depth")
        )
        direct_length = column_values(
            ("MP_length", "MP_Length", "MPLength", "melt_pool_length")
        )
        snap_width = (
            finite_positive(snapshots["width_m"])
            if "width_m" in snapshots.columns
            else np.array([], dtype=float)
        )
        snap_depth = (
            finite_positive(snapshots["depth_m"])
            if "depth_m" in snapshots.columns
            else np.array([], dtype=float)
        )
        snap_length = (
            finite_positive(snapshots["length_m"])
            if "length_m" in snapshots.columns
            else np.array([], dtype=float)
        )
        width = direct_width if len(direct_width) else snap_width
        depth = direct_depth if len(direct_depth) else snap_depth
        length = direct_length if len(direct_length) else snap_length
        source_parts = []
        if len(direct_width) or len(direct_depth):
            source_parts.append("Solidification MP statistics")
        if (not len(direct_width) and len(snap_width)) or (
            not len(direct_depth) and len(snap_depth)
        ):
            source_parts.append("liquidus isosurface from snapshots")
        return width, depth, length, " + ".join(source_parts) or "unavailable"

    def score(
        self,
        solid: pd.DataFrame,
        snapshots: pd.DataFrame,
        parameters: ProcessParameters,
        solid_files: list[Path],
        snapshot_files: list[Path],
        summaries: dict[str, Any],
        histories: dict[str, Any] | None,
        coverage: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        scoring = self.config.get("scoring", {})
        decision_config = self.config.get("decision", {})
        weights = scoring.get(
            "weights",
            {
                "coverage": 0.20,
                "fusion": 0.30,
                "fusion_margin": 0.10,
                "keyhole_margin": 0.10,
                "pool_consistency": 0.10,
                "thermal_uniformity": 0.12,
                "remelt": 0.08,
            },
        )
        keyhole_limit = float(scoring.get("keyhole_aspect_limit", 0.5))
        cross_case_consistent = (
            (coverage or {})
            .get("provenance", {})
            .get("cross_case_inputs", {})
            .get("consistent")
        )
        geometry_snapshots = (
            pd.DataFrame() if cross_case_consistent is False else snapshots
        )
        width, depth, length, geometry_source = self.geometry_values(
            solid, geometry_snapshots
        )

        # Resolution diagnostics must use the grid that produced the geometry,
        # not the potentially different solidification Domain grid.
        geometry_grid_x = parameters.grid_x_m
        geometry_grid_z = parameters.grid_z_m
        geometry_from_snapshots = "snapshots" in (geometry_source or "")
        if geometry_from_snapshots:
            if "grid_x_m" in geometry_snapshots.columns:
                sampled = finite_positive(geometry_snapshots["grid_x_m"])
                if len(sampled):
                    geometry_grid_x = float(np.median(sampled))
            if "grid_z_m" in geometry_snapshots.columns:
                sampled = finite_positive(geometry_snapshots["grid_z_m"])
                if len(sampled):
                    geometry_grid_z = float(np.median(sampled))

        component_scores: dict[str, float | None] = {
            "coverage": None,
            "fusion": None,
            "fusion_margin": None,
            "keyhole_margin": None,
            "pool_consistency": None,
            "thermal_uniformity": None,
            "remelt": None,
        }
        metrics: dict[str, Any] = {
            "melt_pool_geometry_source": geometry_source,
            "melt_pool_width_median_m": float(np.median(width)) if len(width) else None,
            "melt_pool_depth_median_m": float(np.median(depth)) if len(depth) else None,
            "melt_pool_length_median_m": float(np.median(length)) if len(length) else None,
        }
        if "max_temperature_k_diagnostic_only" in geometry_snapshots.columns:
            peak_values = finite_positive(
                geometry_snapshots["max_temperature_k_diagnostic_only"]
            )
            if len(peak_values):
                metrics["snapshot_peak_temperature_max_k_diagnostic_only"] = float(
                    np.max(peak_values)
                )
        model_config = self.config.get("model", {})
        calibration_status = str(
            model_config.get("calibration_status", "uncalibrated")
        ).strip().upper()
        experimentally_calibrated = calibration_status in {
            "EXPERIMENTALLY_CALIBRATED",
            "EXPERIMENTALLY_VALIDATED",
            "VALIDATED",
        }
        model_calibration = {
            "status": calibration_status,
            "experimentally_calibrated": experimentally_calibrated,
            "basis": str(
                model_config.get("calibration_basis", "engineering prior")
            ),
            "parameters_requiring_calibration": list(
                model_config.get(
                    "parameters_requiring_calibration",
                    ["Beam.Intensity.Efficiency", "Beam.Shape.Depth_Z"],
                )
            ),
            "meaning": (
                "Experimentally constrained; interpret within the validated range"
                if experimentally_calibrated
                else "Conditional prediction; not a physical-build qualification"
            ),
        }
        flags: list[str] = []
        limitations: list[str] = [
            "3DThesis is a conduction-only semi-analytical thermal model; it does not resolve "
            "melt-pool fluid flow, recoil pressure, powder stochasticity, gas pores, or spatter.",
            "PASS means model-screening pass, not experimental proof or part qualification.",
        ]
        if self.stress_test:
            limitations.append(
                "Stress-test mode accepts large spatial/parameter offsets as intentional "
                "inputs; it still preserves physical FAIL criteria and never merges "
                "incompatible cases."
            )
        coverage = coverage or {
            "status": "INSUFFICIENT_EVIDENCE",
            "authoritative": False,
            "method": None,
            "coverage_fraction": None,
            "minimum_required_fraction": float(
                self.config.get("coverage", {}).get(
                    "minimum_coverage_fraction", 0.99
                )
            ),
            "snapshot_cumulative_temperature": {"available": False},
            "num_melt_target_grid": {"available": "numMelt" in solid.columns},
            "provenance": {},
        }
        data_gaps: list[dict[str, str]] = []
        if not experimentally_calibrated:
            data_gaps.append(
                {
                    "severity": "advisory",
                    "item": "heat-source experimental calibration",
                    "why": (
                        "effective absorptivity and volumetric source depth are engineering "
                        "priors, so all thermal dimensions are conditional on those inputs"
                    ),
                    "required_3dthesis_output": (
                        "run the C single-track case over an Efficiency×Depth_Z sweep and "
                        "fit width/depth to same-condition metallography; then set "
                        "model.calibration_status to experimentally_calibrated"
                    ),
                }
            )
        pass_fraction = float(
            coverage.get("minimum_required_fraction", 0.99)
        )
        hard_fail_fraction = float(
            self.config.get("coverage", {}).get("hard_fail_fraction", 0.95)
        )
        coverage_fraction = coverage.get("coverage_fraction")
        if coverage.get("authoritative") and coverage_fraction is not None:
            num_melt_planes = coverage.get("num_melt_target_grid", {}).get(
                "planes", {}
            )
            snapshot_planes_authoritative = coverage.get(
                "snapshot_cumulative_temperature", {}
            ).get("planes", {})

            def plane_fraction(label: str) -> float | None:
                value = num_melt_planes.get(label, {}).get(
                    "melted_fraction_with_missing_as_unmelted"
                )
                if value is None:
                    value = snapshot_planes_authoritative.get(label, {}).get(
                        "observed_molten_fraction_lower_bound"
                    )
                return float(value) if value is not None else None

            top_fraction = plane_fraction("top_surface")
            interface_fraction = plane_fraction("previous_layer_interface")
            # Backward compatibility for callers that provide only one authoritative
            # aggregate coverage fraction.
            if top_fraction is None:
                top_fraction = float(coverage_fraction)
            if interface_fraction is None:
                interface_fraction = float(coverage_fraction)
            component_scores["coverage"] = clipped(100.0 * top_fraction)
            component_scores["fusion"] = clipped(100.0 * interface_fraction)
            metrics.update(
                {
                    "target_region_cumulative_coverage_fraction": float(
                        coverage_fraction
                    ),
                    "top_surface_cumulative_coverage_fraction": top_fraction,
                    "interface_cumulative_fusion_fraction": interface_fraction,
                }
            )
            for label, value in (
                ("top surface", top_fraction),
                ("previous-layer interface", interface_fraction),
            ):
                if value < hard_fail_fraction:
                    flags.append(
                        f"FAIL: authoritative {label} melt coverage is below "
                        f"{100 * hard_fail_fraction:.1f}%"
                    )
                elif value < pass_fraction:
                    flags.append(
                        f"REVIEW: authoritative {label} melt coverage is below "
                        f"the {100 * pass_fraction:.1f}% release target"
                    )
        else:
            flags.append(
                "REVIEW: cumulative melt coverage of the complete target region is not proven"
            )
            tracking_mode = coverage.get("num_melt_target_grid", {}).get(
                "tracking_mode"
            )
            data_gaps.append(
                {
                    "severity": "critical",
                    "item": "authoritative cumulative target coverage",
                    "why": (
                        f"numMelt tracking mode is {tracking_mode!r} and/or does not "
                        "contain the complete target grid; sparse snapshots are only a lower bound"
                    ),
                    "required_3dthesis_output": (
                        "numMelt on the complete top/interface target grid including zero-melt "
                        "points (recommended), or sufficiently dense temperature histories/"
                        "snapshots that cover every powered path interval"
                    ),
                }
            )
        snapshot_coverage = coverage.get("snapshot_cumulative_temperature", {})
        snapshot_planes = snapshot_coverage.get("planes", {})
        for label, metric_name in (
            ("top_surface", "snapshot_top_coverage_lower_bound"),
            (
                "previous_layer_interface",
                "snapshot_interface_coverage_lower_bound",
            ),
        ):
            value = snapshot_planes.get(label, {}).get(
                "observed_molten_fraction_lower_bound"
            )
            if value is not None:
                metrics[metric_name] = float(value)
        metrics["snapshot_estimated_minimum_count_for_path"] = snapshot_coverage.get(
            "estimated_minimum_snapshot_count"
        )
        metrics["snapshot_temporal_sampling_ratio"] = snapshot_coverage.get(
            "temporal_sampling_ratio"
        )
        scan_fraction_info = snapshot_coverage.get("scan_fraction_metadata", {})
        raw_scan_fractions = np.asarray(
            scan_fraction_info.get("values", []), dtype=float
        )
        raw_scan_fractions = raw_scan_fractions[np.isfinite(raw_scan_fractions)]
        if len(raw_scan_fractions):
            # Accept either fractions (0..1) or percentages (0..100).
            normalized_scan_fractions = raw_scan_fractions.copy()
            if np.max(np.abs(normalized_scan_fractions)) > 1.5:
                normalized_scan_fractions = normalized_scan_fractions / 100.0
            scan_fraction_span = float(
                np.clip(
                    np.max(normalized_scan_fractions)
                    - np.min(normalized_scan_fractions),
                    0.0,
                    1.0,
                )
            )
            metrics["snapshot_scan_fraction_min"] = float(
                np.min(normalized_scan_fractions)
            )
            metrics["snapshot_scan_fraction_max"] = float(
                np.max(normalized_scan_fractions)
            )
            metrics["snapshot_scan_fraction_span"] = scan_fraction_span
            minimum_path_span = float(
                self.config.get("evidence", {}).get(
                    "minimum_snapshot_path_span_fraction", 0.80
                )
            )
            if scan_fraction_span < minimum_path_span:
                data_gaps.append(
                    {
                        "severity": "advisory",
                        "item": "snapshot path representativeness",
                        "why": (
                            f"snapshot geometry spans only {100 * scan_fraction_span:.1f}% "
                            f"of the executed path progress; it is a local ROI rather than "
                            "a whole-layer geometry sample"
                        ),
                        "required_3dthesis_output": (
                            "distribute geometry snapshots across hatch interior, contour, "
                            "turns and endpoints over at least 80% of path progress"
                        ),
                    }
                )
        num_melt_coverage = coverage.get("num_melt_target_grid", {})
        for label, metric_name in (
            ("top_surface", "num_melt_top_reported_target_fraction_provisional"),
            (
                "previous_layer_interface",
                "num_melt_interface_reported_target_fraction_provisional",
            ),
        ):
            value = num_melt_coverage.get("planes", {}).get(label, {}).get(
                "melted_fraction_with_missing_as_unmelted"
            )
            if value is not None:
                metrics[metric_name] = float(value)
        if (
            snapshot_coverage.get("available")
            and snapshot_coverage.get("metadata_file_count_matches") is False
        ):
            flags.append(
                "REVIEW: snapshot file count does not match Mode.txt ScanFracs"
            )
            data_gaps.append(
                {
                    "severity": "critical",
                    "item": "snapshot schedule metadata",
                    "why": (
                        "the number of snapshot CSVs differs from the current ScanFracs list"
                    ),
                    "required_3dthesis_output": (
                        "retain the exact snapshot Mode.txt from the same run so every "
                        "snapshot time/fraction is known"
                    ),
                }
            )

        provenance = coverage.get("provenance", {})
        configuration_provenance_failures = [
            name
            for name, details in provenance.items()
            if name in ("solidification", "snapshots")
            and details.get("consistent") is False
        ]
        if configuration_provenance_failures:
            flags.append(
                "REVIEW: present case configuration files are newer than the scored CSV data"
            )
            data_gaps.append(
                {
                    "severity": "critical",
                    "item": "same-run input provenance",
                    "why": (
                        "The current Mode/Output/other case text files do not have a "
                        "reliably matching timestamp with the CSV outputs"
                    ),
                    "required_3dthesis_output": (
                        "rerun into a clean case directory and retain the exact Mode.txt, "
                        "Output.txt, Domain.txt, Material.txt, Beam.txt and Path.txt with its CSVs"
                    ),
                }
            )
        cross_case = provenance.get("cross_case_inputs", {})
        if cross_case.get("consistent") is False:
            flags.append(
                (
                    "STRESS-TEST: solidification and snapshot physics differ; "
                    "snapshot-derived evidence was kept separate"
                    if self.stress_test
                    else "INFO: solidification and snapshot physics differ; "
                    "snapshot-derived evidence was kept separate"
                )
            )
            data_gaps.append(
                {
                    "severity": "expected-stress" if self.stress_test else "advisory",
                    "item": "cross-case physics mismatch",
                    "why": (
                        "Material/Beam/Path inputs differ between Solidification and "
                        "Snapshots, so their fields cannot be fused into one physical score"
                    ),
                    "required_3dthesis_output": (
                        "for a formal combined assessment, rerun both modes with identical "
                        "Material.txt, Beam.txt and Path.txt; no change is required for this "
                        "intentional stress test"
                    ),
                }
            )
        if parameters.layer_thickness_source.startswith("fallback"):
            flags.append(
                "REVIEW: physical powder-layer thickness was not found in CLI/config"
            )
            data_gaps.append(
                {
                    "severity": "critical",
                    "item": "physical layer thickness",
                    "why": "Domain Z.Res is numerical resolution and is not the powder-layer thickness",
                    "required_3dthesis_output": (
                        "provide the source slicer CLI with --source-cli, or pass "
                        "--layer-thickness-um explicitly"
                    ),
                }
            )
        if parameters.hatch_spacing_source.startswith("fallback"):
            flags.append(
                "REVIEW: hatch spacing could not be reliably inferred from the executed path"
            )
            data_gaps.append(
                {
                    "severity": "critical",
                    "item": "physical hatch spacing",
                    "why": (
                        "the extreme/concatenated path has no reliable parallel-hatch "
                        "structure, so the configured fallback is being used"
                    ),
                    "required_3dthesis_output": (
                        "set config.process.hatch_spacing_m to the intended slicer hatch "
                        "spacing for a formal LOF assessment"
                    ),
                }
            )
        cli_solid_mask = coverage.get("target_definition", {}).get(
            "cli_solid_mask", {}
        )
        if cli_solid_mask.get("applied") is not True:
            data_gaps.append(
                {
                    "severity": "advisory",
                    "item": "physical solid/void target mask",
                    "why": (
                        cli_solid_mask.get("reason")
                        or "the executed Path could not be matched to one CLI layer"
                    ),
                    "required_3dthesis_output": (
                        "provide the exact source CLI and preserve the layer/Path "
                        "mapping so intentional holes and exterior voids are excluded"
                    ),
                }
            )

        lof_index = None
        if len(width) and len(depth):
            width_p10 = float(np.quantile(width, 0.10))
            depth_p10 = float(np.quantile(depth, 0.10))
            lof_sensitivity: dict[str, dict[str, float]] = {}
            for label, grid_fraction in (
                ("nominal", 0.0),
                ("midpoint", 0.5),
                ("conservative", 1.0),
            ):
                adjusted_width = max(
                    width_p10 - grid_fraction * geometry_grid_x,
                    0.25 * geometry_grid_x,
                )
                adjusted_depth = max(
                    depth_p10 - grid_fraction * geometry_grid_z,
                    0.25 * geometry_grid_z,
                )
                index = (
                    (parameters.hatch_spacing_m / adjusted_width) ** 2
                    + (parameters.layer_thickness_m / adjusted_depth) ** 2
                )
                lof_sensitivity[label] = {
                    "grid_cell_subtraction_fraction": grid_fraction,
                    "width_m": adjusted_width,
                    "depth_m": adjusted_depth,
                    "index": index,
                    "score": logistic_score(1.0 - index, 8.0),
                }
            nominal = lof_sensitivity["nominal"]
            midpoint = lof_sensitivity["midpoint"]
            conservative = lof_sensitivity["conservative"]
            lof_index = midpoint["index"]
            component_scores["fusion_margin"] = midpoint["score"]
            resolution_sensitive = (
                nominal["index"] <= 1.0 < conservative["index"]
            )
            metrics.update(
                {
                    "lof_index_nominal": float(nominal["index"]),
                    "lof_index_midpoint": float(midpoint["index"]),
                    "lof_index_conservative": float(conservative["index"]),
                    "lof_acceptance_limit": 1.0,
                    "lof_resolution_sensitive": resolution_sensitive,
                    "lof_sensitivity": lof_sensitivity,
                    "conservative_width_m": conservative["width_m"],
                    "conservative_depth_m": conservative["depth_m"],
                    "penetration_ratio_conservative": conservative["depth_m"]
                    / parameters.layer_thickness_m,
                    "track_overlap_ratio_conservative": conservative["width_m"]
                    / parameters.hatch_spacing_m,
                }
            )
            if nominal["index"] > 1.0:
                flags.append(
                    "FAIL: even nominal melt-pool geometry predicts lack of fusion"
                )
            elif resolution_sensitive:
                flags.append(
                    "REVIEW: LOF result crosses the acceptance boundary within "
                    "the grid-resolution uncertainty band"
                )
        else:
            flags.append("REVIEW: width/depth evidence is incomplete; LOF cannot be evaluated")

        keyhole_aspect = None
        if len(width) and len(depth):
            keyhole_aspect = float(np.quantile(depth, 0.90) / np.quantile(width, 0.10))
            peak_temperature = metrics.get(
                "snapshot_peak_temperature_max_k_diagnostic_only"
            )
            peak_review_multiple = float(
                scoring.get("peak_temperature_review_multiple_of_liquidus", 3.0)
            )
            thermal_model_out_of_scope = bool(
                peak_temperature is not None
                and peak_temperature
                > peak_review_multiple * parameters.liquidus_k
            )
            metrics.update(
                {
                    "keyhole_aspect_ratio_p90_over_p10": keyhole_aspect,
                    "keyhole_review_limit": keyhole_limit,
                    "keyhole_severe_limit": float(
                        scoring.get("keyhole_severe_aspect_limit", 0.8)
                    ),
                    "thermal_model_out_of_scope_for_keyhole": thermal_model_out_of_scope,
                    "keyhole_score_status": (
                        "NOT_SCORED_MODEL_OUT_OF_SCOPE"
                        if thermal_model_out_of_scope
                        else "SCORED_GEOMETRY_PROXY"
                    ),
                }
            )
            if thermal_model_out_of_scope:
                # A conduction-only solution far beyond the liquidus cannot
                # resolve recoil pressure, evaporation or a vapour cavity.  A
                # numeric "keyhole safety score" would create false precision.
                component_scores["keyhole_margin"] = None
                flags.append(
                    "REVIEW: snapshot peak temperature is far above liquidus; "
                    "keyhole proxy is reported but deliberately not scored"
                )
            else:
                component_scores["keyhole_margin"] = logistic_score(
                    keyhole_limit - keyhole_aspect, 15.0
                )
                if keyhole_aspect > float(
                    scoring.get("keyhole_severe_aspect_limit", 0.8)
                ):
                    flags.append(
                        "FAIL: melt-pool aspect ratio is in the severe keyhole-risk range"
                    )
                elif keyhole_aspect > keyhole_limit:
                    flags.append(
                        "REVIEW: melt-pool aspect ratio crosses the keyhole proxy boundary"
                    )

        pool_cvs = [
            value
            for value in (
                coefficient_of_variation(width),
                coefficient_of_variation(depth),
            )
            if value is not None
        ]
        if pool_cvs:
            mean_pool_cv = float(np.mean(pool_cvs))
            scale = float(scoring.get("pool_cv_scale", 0.35))
            component_scores["pool_consistency"] = clipped(
                100.0 * math.exp(-((mean_pool_cv / scale) ** 2))
            )
            metrics["melt_pool_dimension_mean_cv"] = mean_pool_cv

        thermal_columns = [
            column for column in ("dTdt", "G", "V") if column in solid.columns
        ]
        thermal_mads = [
            robust_log_mad(solid[column])
            for column in thermal_columns
        ]
        thermal_mads = [value for value in thermal_mads if value is not None]
        if not thermal_mads and histories:
            cooling = histories.get("CoolingRate", {})
            if cooling.get("median", 0) > 0:
                # A summary cannot reproduce point values; use its p10/p90 spread.
                ratio = max(cooling["p90"] / max(cooling["p10"], 1e-12), 1.0)
                thermal_mads = [0.25 * math.log10(ratio)]
        if thermal_mads:
            mean_log_mad = float(np.mean(thermal_mads))
            scale = float(scoring.get("thermal_log_mad_scale", 0.35))
            component_scores["thermal_uniformity"] = clipped(
                100.0 * math.exp(-((mean_log_mad / scale) ** 2))
            )
            metrics["thermal_log10_mad_mean"] = mean_log_mad

        if "numMelt" in solid.columns:
            remelt = pd.to_numeric(solid["numMelt"], errors="coerce").dropna().to_numpy(float)
            if len(remelt):
                excess = float(np.mean(np.clip(remelt - 2.0, 0.0, None)))
                component_scores["remelt"] = clipped(100.0 * math.exp(-excess / 1.5))
                metrics.update(
                    {
                        "remelt_fraction_ge_2": float(np.mean(remelt >= 2)),
                        "excess_remelt_mean_above_2": excess,
                        "num_melt_max": float(np.max(remelt)),
                    }
                )
        def has_positive_column(names: tuple[str, ...]) -> bool:
            return any(
                name in solid.columns and len(finite_positive(solid[name])) > 0
                for name in names
            )

        direct_width_usable = has_positive_column(
            ("MP_width", "MP_Width", "MPWidth", "melt_pool_width")
        )
        direct_depth_usable = has_positive_column(
            (
                "MP_depth",
                "MP_Depth",
                "MPDepth",
                "melt_pool_depth",
                "depth",
            )
        )
        if not (direct_width_usable and direct_depth_usable):
            data_gaps.append(
                {
                    "severity": "advisory",
                    "item": "usable direct melt-pool statistics",
                    "why": (
                        "MP width/depth columns are absent or contain no positive values; "
                        "geometry is reconstructed from liquidus snapshots"
                    ),
                    "required_3dthesis_output": (
                        "enable MP_Stats for direct width/length and MP_depth when the "
                        "installed 3DThesis version provides it; otherwise retain "
                        "liquidus snapshots for depth"
                    ),
                }
            )

        evidence_config = self.config.get("evidence", {})
        minimum_geometry_snapshots = int(
            evidence_config.get("minimum_geometry_snapshots", 20)
        )
        if (
            "liquidus isosurface from snapshots" in geometry_source
            and len(snapshot_files) < minimum_geometry_snapshots
        ):
            data_gaps.append(
                {
                    "severity": "advisory",
                    "item": "melt-pool geometry sampling",
                    "why": (
                        f"only {len(snapshot_files)} snapshots support width/depth "
                        f"statistics; the configured geometry target is "
                        f"{minimum_geometry_snapshots}"
                    ),
                    "required_3dthesis_output": (
                        "retain snapshots at evenly distributed scan fractions; "
                        "20–50 frames are recommended for geometry variability"
                    ),
                }
            )
        recommended_cells = float(
            evidence_config.get("recommended_cells_across_pool_dimension", 5.0)
        )
        geometry_cells = None
        if len(width) and len(depth):
            geometry_cells = min(
                float(np.quantile(width, 0.10)) / geometry_grid_x,
                float(np.quantile(depth, 0.10)) / geometry_grid_z,
            )
            metrics["minimum_cells_across_p10_pool_dimension"] = geometry_cells
            metrics["geometry_grid_x_m"] = float(geometry_grid_x)
            metrics["geometry_grid_z_m"] = float(geometry_grid_z)
            if geometry_cells < recommended_cells:
                data_gaps.append(
                    {
                        "severity": "advisory",
                        "item": "melt-pool geometry resolution",
                        "why": (
                            f"the limiting P10 pool dimension spans only "
                            f"{geometry_cells:.2f} grid cells; "
                            f"{recommended_cells:.1f} are recommended"
                        ),
                        "required_3dthesis_output": (
                            "refine X/Y/Z resolution, especially Z, and rerun "
                            "snapshots before treating the LOF margin as precise"
                        ),
                    }
                )
        missing_gradient_directions = [
            name for name in ("Gx", "Gy", "Gz") if name not in solid.columns
        ]
        if missing_gradient_directions:
            data_gaps.append(
                {
                    "severity": "optional",
                    "item": "solidification-gradient direction",
                    "why": (
                        "grain-orientation consistency cannot be assessed without "
                        + ", ".join(missing_gradient_directions)
                    ),
                    "required_3dthesis_output": "enable Gx, Gy and Gz if microstructure direction matters",
                }
            )

        available_weight = sum(
            float(weights.get(name, 0.0))
            for name, value in component_scores.items()
            if value is not None
        )
        quality_score = (
            sum(
                float(weights.get(name, 0.0)) * float(value)
                for name, value in component_scores.items()
                if value is not None
            )
            / available_weight
            if available_weight
            else 0.0
        )
        score_breakdown = [
            {
                "component": name,
                "label": COMPONENT_LABELS.get(name, name),
                "score": float(value),
                "configured_weight": float(weights.get(name, 0.0)),
                "effective_weight": (
                    float(weights.get(name, 0.0)) / available_weight
                    if available_weight
                    else 0.0
                ),
                "points_contributed": (
                    float(weights.get(name, 0.0))
                    / available_weight
                    * float(value)
                    if available_weight
                    else 0.0
                ),
            }
            for name, value in component_scores.items()
            if value is not None
        ]

        coordinate_values = solid[["x", "y", "z"]].apply(pd.to_numeric, errors="coerce")
        finite_coordinate_fraction = float(
            np.isfinite(coordinate_values.to_numpy(float)).all(axis=1).mean()
        )
        evidence_parts: dict[str, float] = {}
        evidence_parts["valid_rows"] = 5.0 * finite_coordinate_fraction
        evidence_parts["liquidus"] = (
            5.0 if not parameters.liquidus_source.startswith("fallback") else 2.5
        )
        evidence_parts["layer_thickness"] = (
            5.0 if not parameters.layer_thickness_source.startswith("fallback") else 2.0
        )
        evidence_parts["hatch_spacing"] = (
            5.0 if not parameters.hatch_spacing_source.startswith("fallback") else 2.0
        )
        evidence_parts["solidification_fields"] = (
            20.0
            * sum(column in solid.columns for column in ("G", "V", "dTdt"))
            / 3.0
        )
        evidence_parts["melt_pool_width"] = 10.0 if len(width) else 0.0
        evidence_parts["melt_pool_depth"] = 10.0 if len(depth) else 0.0
        evidence_parts["snapshot_geometry"] = (
            5.0 if len(snapshot_files) and cross_case_consistent is not False else 0.0
        )
        evidence_parts["remelt"] = 5.0 if "numMelt" in solid.columns else 0.0
        evidence_parts["coverage_lower_bound"] = (
            5.0 if snapshot_coverage.get("available") else 0.0
        )
        evidence_parts["coverage_authoritative"] = (
            25.0 if coverage.get("authoritative") else 0.0
        )
        evidence_score = float(sum(evidence_parts.values()))

        authoritative_parameters = sum(
            not source.startswith("fallback")
            for source in (
                parameters.liquidus_source,
                parameters.layer_thickness_source,
                parameters.hatch_spacing_source,
            )
        )
        field_fraction = (
            sum(column in solid.columns for column in ("G", "V", "dTdt")) / 3.0
        )
        num_melt_complete = coverage.get("num_melt_target_grid", {}).get(
            "coordinate_grid_complete"
        )
        if num_melt_complete is None:
            num_melt_complete = coverage.get("authoritative", False)
        geometry_sampling_fraction = (
            1.0
            if direct_width_usable and direct_depth_usable
            else min(
                len(snapshot_files) / max(minimum_geometry_snapshots, 1),
                1.0,
            )
            if len(width) and len(depth)
            else 0.0
        )
        geometry_resolution_fraction = (
            min(geometry_cells / max(recommended_cells, 1e-12), 1.0)
            if geometry_cells is not None
            else 0.0
        )
        geometry_path_representativeness_fraction = (
            1.0
            if direct_width_usable and direct_depth_usable
            else float(metrics.get("snapshot_scan_fraction_span", 0.0))
            if len(width) and len(depth)
            else 0.0
        )
        metrics["geometry_path_representativeness_fraction"] = (
            geometry_path_representativeness_fraction
        )
        evidence_adequacy_parts = {
            "authoritative_complete_coverage": (
                30.0
                if coverage.get("authoritative") and num_melt_complete
                else 0.0
            ),
            "same_run_provenance": (
                15.0
                if not configuration_provenance_failures
                and cross_case_consistent is not False
                else 0.0
            ),
            "physical_parameter_authority": 10.0
            * authoritative_parameters
            / 3.0,
            "solidification_field_support": 10.0 * field_fraction,
            "geometry_sampling_support": 10.0 * geometry_sampling_fraction,
            "geometry_resolution_support": 10.0
            * geometry_resolution_fraction,
            "geometry_path_representativeness": 15.0
            * geometry_path_representativeness_fraction,
        }
        evidence_adequacy_score = float(sum(evidence_adequacy_parts.values()))

        validation: dict[str, Any] = {
            "finite_coordinate_fraction": finite_coordinate_fraction,
            "row_count": int(len(solid)),
            "z_plane_count": int(solid["z"].nunique()),
            "z_is_depth_not_layer": True,
            "spatial_offset_policy": (
                "translation-invariant; large absolute coordinates are not corruption"
            ),
            "coordinate_bounds_m": {
                axis: {
                    "min": float(pd.to_numeric(solid[axis], errors="coerce").min()),
                    "max": float(pd.to_numeric(solid[axis], errors="coerce").max()),
                }
                for axis in ("x", "y", "z")
            },
        }
        if {"G", "V", "dTdt"}.issubset(solid.columns):
            g = pd.to_numeric(solid["G"], errors="coerce").to_numpy(float)
            v = pd.to_numeric(solid["V"], errors="coerce").to_numpy(float)
            cooling = pd.to_numeric(solid["dTdt"], errors="coerce").to_numpy(float)
            valid = np.isfinite(g) & np.isfinite(v) & np.isfinite(cooling)
            relative_error = np.abs(cooling[valid] - g[valid] * v[valid]) / np.maximum(
                np.abs(cooling[valid]), 1e-12
            )
            validation["median_relative_error_dTdt_vs_G_times_V"] = (
                float(np.median(relative_error)) if len(relative_error) else None
            )

        if {"Gx", "Gy", "Gz"}.issubset(solid.columns):
            directions = solid[["Gx", "Gy", "Gz"]].apply(
                pd.to_numeric, errors="coerce"
            ).to_numpy(float)
            valid = np.isfinite(directions).all(axis=1)
            directions = directions[valid]
            if len(directions):
                norms = np.linalg.norm(directions, axis=1)
                validation["gradient_direction_norm_median"] = float(np.median(norms))
                validation["gradient_orientation_resultant"] = float(
                    np.linalg.norm(np.mean(directions, axis=0))
                )

        hard_fail = any(flag.startswith("FAIL:") for flag in flags)
        review_flag = any(flag.startswith("REVIEW:") for flag in flags)
        pass_score = float(decision_config.get("pass_score", 75.0))
        minimum_evidence = float(decision_config.get("minimum_pass_evidence", 85.0))
        minimum_adequacy = float(
            decision_config.get("minimum_pass_evidence_adequacy", 70.0)
        )
        if hard_fail:
            decision = "FAIL"
            predicted_success = False
        elif (
            quality_score >= pass_score
            and evidence_score >= minimum_evidence
            and evidence_adequacy_score >= minimum_adequacy
            and not review_flag
        ):
            decision = "PASS"
            predicted_success = True
        else:
            decision = "REVIEW"
            predicted_success = None

        if quality_score >= 90:
            continuous_grade = "A / Excellent continuous score"
        elif quality_score >= 80:
            continuous_grade = "B / Good continuous score"
        elif quality_score >= 70:
            continuous_grade = "C / Moderate continuous score"
        elif quality_score >= 60:
            continuous_grade = "D / Marginal continuous score"
        else:
            continuous_grade = "F / Poor continuous score"
        grade = continuous_grade + " (not a release decision)"

        problem_diagnosis = build_problem_diagnosis(
            coverage=coverage,
            metrics=metrics,
            component_scores=component_scores,
            parameters=parameters,
            snapshot_count=len(snapshot_files),
            minimum_geometry_snapshots=minimum_geometry_snapshots,
            hard_fail_fraction=hard_fail_fraction,
            pass_fraction=pass_fraction,
        )
        path_span = metrics.get("snapshot_scan_fraction_span")
        minimum_path_span = float(
            evidence_config.get("minimum_snapshot_path_span_fraction", 0.80)
        )
        if path_span is not None and float(path_span) < minimum_path_span:
            problem_diagnosis.append(
                {
                    "severity": "WARNING",
                    "category": "Evidence representativeness",
                    "where": "Snapshot path coverage",
                    "finding": "Melt-pool geometry represents only a local scan window and "
                    "cannot be extrapolated directly to the full layer",
                    "evidence": (
                        f"ScanFracs span only {100 * float(path_span):.1f}% of path progress"
                    ),
                    "likely_causes": ["Snapshots are concentrated in one local ROI"],
                    "recommended_actions": [
                        "Cover at least 80% of the path and stratify snapshots across hatch, "
                        "contour, corner, and endpoint segments"
                    ],
                    "confidence": "HIGH",
                }
            )
        if not experimentally_calibrated:
            problem_diagnosis.append(
                {
                    "severity": "ADVISORY",
                    "category": "Model calibration",
                    "where": "Heat-source Efficiency / Depth_Z",
                    "finding": "The heat source has not been constrained by a matched experiment; "
                    "the result is conditional",
                    "evidence": f"calibration_status={calibration_status}",
                    "likely_causes": ["No single-track metallographic width and depth are available"],
                    "recommended_actions": [
                        "Run a two-dimensional sensitivity sweep first, then use case C to "
                        "jointly calibrate width and depth after experiments"
                    ],
                    "confidence": "HIGH",
                }
            )

        line_energy = None
        volumetric_energy = None
        if (
            parameters.power_w is not None
            and parameters.efficiency is not None
            and parameters.scan_speed_m_s
        ):
            line_energy = (
                parameters.efficiency * parameters.power_w / parameters.scan_speed_m_s
            )
            volumetric_energy = line_energy / (
                parameters.hatch_spacing_m * parameters.layer_thickness_m
            )
            metrics["absorbed_line_energy_j_m"] = line_energy
            metrics["absorbed_volumetric_energy_j_m3_diagnostic_only"] = volumetric_energy

        available_columns = list(solid.columns)
        process_diagnostics = build_process_diagnostics(solid, snapshots)
        known_columns = {
            "x", "y", "z", "T", "tSol", "G", "Gx", "Gy", "Gz", "V", "dTdt",
            "eqFrac", "depth", "numMelt", "MP_width", "MP_Width", "MPWidth",
            "MP_length", "MP_Length", "MPLength", "MP_depth", "MP_Depth", "MPDepth",
            "H", "Hx", "Hy", "Hz",
        }
        if self.stress_test:
            data_readiness = "STRESS_TEST"
        elif (
            not coverage.get("authoritative")
            or configuration_provenance_failures
            or evidence_adequacy_score < minimum_adequacy
        ):
            data_readiness = "PROVISIONAL"
        elif not experimentally_calibrated:
            data_readiness = "MODEL_SCREENING"
        else:
            data_readiness = "QUALIFIABLE"

        return {
            "schema_version": "2.1",
            "engine": {
                "name": "LPBF_Agent physics-informed scorer",
                "version": ENGINE_VERSION,
                "scorer_file": str(Path(__file__).resolve()),
                "config_file": str(self.config_path.resolve()),
            },
            "assessment_mode": "STRESS_TEST" if self.stress_test else "STANDARD",
            "layer_id": self.layer_id,
            "decision": decision,
            "predicted_success": predicted_success,
            "quality_score": round(quality_score, 3),
            "grade": grade,
            "continuous_grade": continuous_grade,
            "quality_score_semantics": (
                "continuous process-margin index; hard release gates take precedence"
            ),
            "evidence_completeness_score": round(evidence_score, 3),
            "evidence_adequacy_score": round(evidence_adequacy_score, 3),
            "component_scores": {
                key: round(value, 3) if value is not None else None
                for key, value in component_scores.items()
            },
            "component_weights": weights,
            "score_breakdown": score_breakdown,
            "metrics": metrics,
            "flags": flags,
            "data_readiness": data_readiness,
            "model_calibration": model_calibration,
            "data_gaps": data_gaps,
            "problem_diagnosis": problem_diagnosis,
            "limitations": limitations,
            "parameters": parameters.__dict__,
            "evidence_parts": evidence_parts,
            "evidence_adequacy_parts": evidence_adequacy_parts,
            "validation": validation,
            "input": {
                "solidification_files": [str(file) for file in solid_files],
                "snapshot_files": [str(file) for file in snapshot_files],
                "available_columns": available_columns,
                "recognized_columns": sorted(set(available_columns) & known_columns),
                "additional_numeric_columns_are_summarized": sorted(
                    set(available_columns) - known_columns
                ),
            },
            "thermal_history": histories,
            "process_diagnostics": process_diagnostics,
            "coverage": coverage,
            "feature_summary": summaries,
        }

    def write_outputs(
        self,
        result: dict[str, Any],
        snapshots: pd.DataFrame,
        summaries: dict[str, Any],
    ) -> None:
        json_path = self.output_dir / "assessment.json"
        json_path.write_text(
            json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False),
            encoding="utf-8",
        )
        if not snapshots.empty:
            snapshots.to_csv(self.output_dir / "snapshot_features.csv", index=False)
        pd.DataFrame(result.get("score_breakdown", [])).to_csv(
            self.output_dir / "Score_Breakdown.csv", index=False
        )
        problem_rows = []
        for problem in result.get("problem_diagnosis", []):
            row = dict(problem)
            row["likely_causes"] = "; ".join(row.get("likely_causes", []))
            row["recommended_actions"] = "; ".join(
                row.get("recommended_actions", [])
            )
            problem_rows.append(row)
        pd.DataFrame(
            problem_rows,
            columns=[
                "severity",
                "category",
                "where",
                "finding",
                "evidence",
                "likely_causes",
                "recommended_actions",
                "confidence",
            ],
        ).to_csv(self.output_dir / "Problem_Diagnosis.csv", index=False)
        coverage_rows: list[dict[str, Any]] = []
        planes = (
            result.get("coverage", {})
            .get("num_melt_target_grid", {})
            .get("planes", {})
        )
        for plane_name, plane in planes.items():
            spatial = plane.get("spatial_diagnostics", {})
            coverage_rows.extend(
                [
                    {
                        "plane": plane_name,
                        "scope": "whole_plane",
                        "distance_um_low": None,
                        "distance_um_high": None,
                        "target_count": spatial.get("target_count"),
                        "melted_count": spatial.get("melted_count"),
                        "melted_fraction": plane.get(
                            "melted_fraction_with_missing_as_unmelted"
                        ),
                        "unmelted_count": spatial.get("unmelted_count"),
                    },
                    {
                        "plane": plane_name,
                        "scope": "track_interior",
                        "distance_um_low": None,
                        "distance_um_high": None,
                        "target_count": spatial.get(
                            "track_interior_target_count"
                        ),
                        "melted_count": spatial.get(
                            "track_interior_melted_count"
                        ),
                        "melted_fraction": spatial.get(
                            "track_interior_melted_fraction"
                        ),
                        "unmelted_count": (
                            (spatial.get("track_interior_target_count") or 0)
                            - (spatial.get("track_interior_melted_count") or 0)
                        ),
                    },
                    {
                        "plane": plane_name,
                        "scope": "scan_endpoints",
                        "distance_um_low": None,
                        "distance_um_high": None,
                        "target_count": spatial.get("endpoint_target_count"),
                        "melted_count": spatial.get("endpoint_melted_count"),
                        "melted_fraction": spatial.get(
                            "endpoint_melted_fraction"
                        ),
                        "unmelted_count": (
                            (spatial.get("endpoint_target_count") or 0)
                            - (spatial.get("endpoint_melted_count") or 0)
                        ),
                    },
                ]
            )
            for band in spatial.get("distance_bands", []):
                coverage_rows.append(
                    {
                        "plane": plane_name,
                        "scope": "distance_from_track_centerline",
                        "distance_um_low": band.get("distance_um_low"),
                        "distance_um_high": band.get("distance_um_high"),
                        "target_count": band.get("target_count"),
                        "melted_count": band.get("melted_count"),
                        "melted_fraction": band.get("melted_fraction"),
                        "unmelted_count": (
                            band.get("target_count", 0)
                            - band.get("melted_count", 0)
                        ),
                    }
                )
        pd.DataFrame(coverage_rows).to_csv(
            self.output_dir / "Coverage_Diagnostics.csv", index=False
        )
        (self.output_dir / "feature_summary.json").write_text(
            json.dumps(summaries, ensure_ascii=False, indent=2, allow_nan=False),
            encoding="utf-8",
        )
        score_row = {
            "Layer": self.layer_id,
            "Decision": result["decision"],
            "PredictedSuccess": result["predicted_success"],
            "QualityScore": result["quality_score"],
            "Grade": result["grade"],
            "EvidenceCompletenessScore": result["evidence_completeness_score"],
            "EvidenceAdequacyScore": result.get("evidence_adequacy_score"),
            **{
                f"{name}Score": value
                for name, value in result["component_scores"].items()
            },
            **{
                key: value
                for key, value in result["metrics"].items()
                if isinstance(value, (int, float, bool)) or value is None
            },
        }
        pd.DataFrame([score_row]).to_csv(
            self.output_dir / "Layer_Assessment.csv", index=False
        )
        (self.output_dir / "Assessment_Report.md").write_text(
            self.report_markdown(result), encoding="utf-8"
        )
        write_intuitive_outputs(self.output_dir, result)

    def report_markdown(self, result: dict[str, Any]) -> str:
        metrics = result["metrics"]
        parameters = result["parameters"]
        components = result["component_scores"]
        coverage = result.get("coverage", {})
        success_text = {
            True: "Predicted success",
            False: "Predicted failure",
            None: "Review required",
        }[result["predicted_success"]]

        def number(value: Any, digits: int = 3) -> str:
            return "N/A" if value is None else f"{float(value):.{digits}f}"

        problems = result.get("problem_diagnosis", [])
        primary_problem = next(
            (
                problem
                for problem in problems
                if problem.get("severity") == "CRITICAL"
            ),
            problems[0] if problems else None,
        )
        lines = [
            f"# LPBF Single-Layer Assessment Technical Appendix: {result['layer_id']}",
            "",
            "> Open `Dashboard.html` or `00_READ_ME_FIRST.md` first. This appendix "
            "retains the equations, data provenance, and full technical detail.",
            "",
            f"- Decision: **{result['decision']} / {success_text}**",
            f"- Assessment mode: **{result.get('assessment_mode', 'STANDARD')}**",
            f"- Quality score: **{result['quality_score']:.1f}/100** ({result['grade']})",
            f"- Evidence completeness: **{result['evidence_completeness_score']:.1f}/100** "
            "(whether the required files and fields are present)",
            f"- Evidence adequacy: **{result.get('evidence_adequacy_score', 0):.1f}/100** "
            "(whether sampling count, path span, and grid resolution support a detailed conclusion)",
            f"- Data readiness: **{result.get('data_readiness', 'PROVISIONAL')}**",
            f"- Model calibration: **{result.get('model_calibration', {}).get('status', 'UNKNOWN')}**",
            "",
            "> This is a process-screening result based on 3DThesis thermal-conduction "
            "outputs. It is not CT, metallography, or density inspection and cannot replace "
            "physical-part qualification.",
            "",
            "## Conclusion Summary",
            "",
            (
                f"**Primary issue: {primary_problem['where']} — "
                f"{primary_problem['finding']}.**"
                if primary_problem
                else "**No explicit anomaly reached the reporting threshold.**"
            ),
            "",
            "High evidence completeness means only that the inputs are present. The quality "
            "score and FAIL/PASS decision depend on actual melt coverage, fusion margin, and "
            "thermal-field behavior.",
            "",
            "### Problem Localization and Cause Inference",
            "",
            "| Severity | Location | Finding | Key evidence | Likely causes | Recommended actions | Confidence |",
            "|---|---|---|---|---|---|---|",
        ]
        for problem in problems:
            lines.append(
                f"| {problem['severity']} | {problem['where']} | "
                f"{problem['finding']} | {problem['evidence']} | "
                f"{'; '.join(problem.get('likely_causes', []))} | "
                f"{'; '.join(problem.get('recommended_actions', []))} | "
                f"{problem['confidence']} |"
            )
        if not problems:
            lines.append(
                "| INFO | Full layer | No explicit anomaly found | All current thresholds met | - | - | - |"
            )
        lines.extend(
            [
                "",
                "## How the Score Is Calculated",
                "",
                "| Component | Raw score | Effective weight | Contribution |",
                "|---|---:|---:|---:|",
            ]
        )
        for item in result.get("score_breakdown", []):
            lines.append(
                f"| {item['label']} | {item['score']:.1f} | "
                f"{100 * item['effective_weight']:.1f}% | "
                f"{item['points_contributed']:.2f} |"
            )
        lines.extend(
            [
                "",
                "Top-surface scan coverage and previous-layer interface fusion are scored "
                "separately. LOF carries only 10% as a geometric margin and no longer duplicates "
                "the 30% penalty already represented by authoritative interface `numMelt`.",
                "",
                "## Core Criteria and Layer Values",
                "",
                "### 1. Top-Surface Coverage and Interlayer Fusion (Authoritative numMelt)",
                "",
                r"$$C_{top}=\frac{N_{top}(numMelt\ge1)}{N_{top,target}},\qquad"
                r"C_{interface}=\frac{N_{interface}(numMelt\ge1)}"
                r"{N_{interface,target}}$$",
            "",
            f"- Evidence status: `{coverage.get('status', 'INSUFFICIENT_EVIDENCE')}`",
            f"- Decision method: `{coverage.get('method') or 'N/A'}`",
            f"- Authoritative top-surface coverage: "
            f"`{number(metrics.get('top_surface_cumulative_coverage_fraction'))}`",
            f"- Authoritative previous-layer interface fusion: "
            f"`{number(metrics.get('interface_cumulative_fusion_fraction'))}`",
            f"- CLI solid/void mask applied: "
            f"`{coverage.get('target_definition', {}).get('cli_solid_mask', {}).get('applied', False)}`",
            f"- Target-region method: "
            f"`{coverage.get('target_definition', {}).get('cli_solid_mask', {}).get('method', coverage.get('target_definition', {}).get('method', 'N/A'))}`",
            f"- Matched CLI layer: "
            f"`{coverage.get('target_definition', {}).get('cli_solid_mask', {}).get('matched_layer_number', 'N/A')}`; "
            f"excluded void points: "
            f"`{coverage.get('target_definition', {}).get('cli_solid_mask', {}).get('excluded_void_point_count', 0)}` per plane",
            f"- Snapshot top-coverage lower bound: `{number(metrics.get('snapshot_top_coverage_lower_bound'))}`",
            f"- Snapshot interface-coverage lower bound: `{number(metrics.get('snapshot_interface_coverage_lower_bound'))}`",
            f"- `{coverage.get('num_melt_target_grid', {}).get('tracking_mode') or 'unknown'} "
            "numMelt` top-surface recorded fraction (reference only when non-authoritative): "
            f"`{number(metrics.get('num_melt_top_reported_target_fraction_provisional'))}`",
            f"- `{coverage.get('num_melt_target_grid', {}).get('tracking_mode') or 'unknown'} "
            "numMelt` interface recorded fraction (reference only when non-authoritative): "
            f"`{number(metrics.get('num_melt_interface_reported_target_fraction_provisional'))}`",
            f"- Snapshot count / estimated minimum: `{len(result['input']['snapshot_files'])}` / "
            f"`{metrics.get('snapshot_estimated_minimum_count_for_path') or 'N/A'}`",
            f"- Path-progress span covered by snapshots: "
            f"`{number(metrics.get('snapshot_scan_fraction_span'))}`",
            f"- Recommended Custom coverage-point file: "
            f"`{coverage.get('recommended_custom_point_file') or 'N/A'}`"
            f" ({coverage.get('recommended_custom_point_count') or 0} points)",
            "",
            "Only `numMelt` over the complete target grid, including zero-valued points, "
            "or demonstrably dense cumulative temperature evidence can decide this component. "
            "Absence of melting in a sparse snapshot may simply mean the laser-passage instant "
            "was not captured, so snapshots provide a lower bound rather than a false LOF verdict.",
            "",
            "See `Coverage_Diagnostics.csv` for track-interior, endpoint, boundary, and "
            "centerline-distance groups.",
            "",
            "### 2. Lack-of-Fusion (LOF) Geometric Criterion",
            "",
            r"$$I_{LOF}=\left(\frac{h}{W_c}\right)^2+\left(\frac{t}{D_c}\right)^2$$",
            "",
            r"$I_{LOF}\leq1$ indicates geometric overlap between adjacent tracks and the "
            "previous layer under an elliptical melt-pool approximation. The system reports "
            "nominal, half-cell, and one-cell-corrected results. The score uses the half-cell "
            "value; an interval crossing 1 is marked resolution-sensitive instead of being "
            "forced to fail by the most conservative estimate.",
            "",
            f"- Nominal LOF index (no grid correction): `{number(metrics.get('lof_index_nominal'))}`",
            f"- Half-cell LOF index (used for scoring): `{number(metrics.get('lof_index_midpoint'))}`",
            f"- Conservative LOF index: `{number(metrics.get('lof_index_conservative'))}` (limit 1.0)",
            f"- Resolution sensitive: `{metrics.get('lof_resolution_sensitive', 'N/A')}`",
            f"- Conservative width/hatch-spacing ratio: `{number(metrics.get('track_overlap_ratio_conservative'))}`",
            f"- Conservative depth/layer-thickness ratio: `{number(metrics.get('penetration_ratio_conservative'))}`",
            "",
            "The LOF-margin score uses the half-cell index and the interpretable logistic "
            r"function $S=100/[1+\exp(-8(1-I_{LOF}))]$; the physical boundary maps to 50 points.",
            "",
            "### 3. Keyhole-Risk Proxy",
            "",
            r"$$A_{KH}=\frac{P_{90}(D)}{P_{10}(W)}$$",
            "",
            f"- Layer depth-to-width ratio: `{number(metrics.get('keyhole_aspect_ratio_p90_over_p10'))}`",
            f"- Review boundary: `{number(metrics.get('keyhole_review_limit'))}`",
            f"- Scoring status: `{metrics.get('keyhole_score_status', 'N/A')}`",
            "",
            "A depth-to-width ratio near 0.5 is a common conduction/transition-keyhole "
            "morphology boundary. Because 3DThesis excludes flow and recoil pressure, this "
            "criterion is only a geometric risk proxy.",
            "",
            "### 4. Thermal Field and Solidification",
            "",
            r"3DThesis outputs satisfy $\dot T\approx G\,V$. The ratio $G/V$ is primarily "
            r"associated with solidification morphology, while $G\,V$ (cooling rate) is "
            "primarily associated with microstructural scale. The system does not assume one "
            "microstructure is universally superior; it evaluates within-layer dispersion and "
            "retains absolute values for microstructure-specific objectives.",
            "",
            f"- Median relative error between `dTdt` and `G×V`: "
            f"`{number(result['validation'].get('median_relative_error_dTdt_vs_G_times_V'), 6)}`",
            f"- Mean thermal-parameter log10-MAD: `{number(metrics.get('thermal_log10_mad_mean'))}`",
            "",
            "### 5. Process-Parameter Diagnostics",
            "",
            f"- Liquidus temperature: `{parameters['liquidus_k']:.1f} K`",
            f"- Layer thickness: `{parameters['layer_thickness_m'] * 1e6:.1f} µm`",
            f"  - Source: `{parameters['layer_thickness_source']}`. `Domain.txt/Z.Res` is used "
            "only as the numerical grid spacing, not as powder-layer thickness.",
            f"- Hatch spacing: `{parameters['hatch_spacing_m'] * 1e6:.1f} µm`",
            f"- Melt-pool dimension source: `{metrics.get('melt_pool_geometry_source')}`",
            f"- Median melt-pool width: `{number(metrics.get('melt_pool_width_median_m', None) * 1e6 if metrics.get('melt_pool_width_median_m') is not None else None, 1)} µm`",
            f"- Median melt-pool depth: `{number(metrics.get('melt_pool_depth_median_m', None) * 1e6 if metrics.get('melt_pool_depth_median_m') is not None else None, 1)} µm`",
            "",
            ]
        )
        path_diagnostics = result.get("input", {}).get(
            "path_parse_diagnostics", {}
        )
        cross_case = coverage.get("provenance", {}).get(
            "cross_case_inputs", {}
        )
        lines.extend(
            [
                "## Data-Format and Extreme-Test Diagnostics",
                "",
                "- Coordinate handling: only relative distances and the actual path are used. "
                "The model need not be centered at the origin, and a large coordinate offset "
                "alone is not treated as file corruption.",
                f"- Physical Path records: `{path_diagnostics.get('physical_record_count', 'N/A')}`",
                f"- Path lines with concatenated trailing fields: "
                f"`{path_diagnostics.get('trailing_token_line_count', 0)}`. "
                "Only the first record on each physical line is read, matching 3DThesis line-by-line parsing.",
                f"- Solidification/Snapshots inputs can be merged: "
                f"`{cross_case.get('consistent', 'unknown')}`",
            ]
        )
        for key, difference in cross_case.get(
            "parameter_differences", {}
        ).items():
            lines.append(
                f"  - `{key}`: Solidification=`{difference.get('solidification')}`, "
                f"Snapshots=`{difference.get('snapshots')}`"
            )
        rdf = result.get("rdf")
        if rdf:
            lines.extend(
                [
                    f"- RDF: `{rdf.get('row_count', 0)}` rows, used for ExaCA/solidification-history "
                    "diagnostics and not counted again in the quality score.",
                    f"- Median RDF liquid duration: "
                    f"`{number(rdf.get('liquid_duration_s', {}).get('median'), 6)} s`",
                    f"- Median RDF cooling rate: "
                    f"`{number(rdf.get('cooling_rate_k_s', {}).get('median'), 1)} K/s`",
                ]
            )
        output_audit = result.get("output_audit", {})
        if output_audit:
            lines.extend(
                [
                    f"- Fields requested in Output but absent: "
                    f"`{', '.join(output_audit.get('missing_requested_outputs', [])) or 'None'}`",
                    f"- Intermediate Solidification CSV files found and intentionally ignored: "
                    f"`{len(result.get('input', {}).get('ignored_intermediate_solidification_files', []))}`; "
                    "scoring reads only `.Solidification.Final`.",
                ]
            )
        lines.append("")
        if result["flags"]:
            lines.extend(["## Items Requiring Attention", ""])
            lines.extend(f"- {flag}" for flag in result["flags"])
            lines.append("")
        if result.get("data_gaps"):
            lines.extend(
                [
                    "## Missing Data and Recommended 3DThesis Outputs",
                    "",
                    "| Severity | Missing item | Why it is needed | Recommended output/action |",
                    "|---|---|---|---|",
                ]
            )
            for gap in result["data_gaps"]:
                lines.append(
                    f"| {gap['severity']} | {gap['item']} | {gap['why']} | "
                    f"{gap['required_3dthesis_output']} |"
                )
            lines.append("")
        lines.extend(
            [
                "## How to Read the Outputs",
                "",
                "- `Dashboard.html`: open in a browser for the clearest decision, components, and actions.",
                "- `00_READ_ME_FIRST.md`: one-page English summary; read this first.",
                "- `Action_Plan.csv`: prioritized issues and next steps.",
                "- `Assessment_Report.md`: equations, provenance, and full technical detail.",
                "- `Problem_Diagnosis.csv`: location, evidence, inferred causes, and actions for each issue.",
                "- `Coverage_Diagnostics.csv`: top/interface groups by track interior, endpoint, and centerline distance.",
                "- `Score_Breakdown.csv`: raw component scores, effective weights, and score contributions.",
                "- `assessment.json`: complete machine-readable result, including evidence adequacy and LOF interval.",
                "",
                "## Model Boundaries",
                "",
                "- `z` is simulated depth below the same build layer, not a layer number.",
                "- Spatter, balling, gas-entrapped pores, and physical cracks require a fluid/mechanical "
                "model or experimental labels. This system will not produce falsely certain defect "
                "probabilities from a temperature field alone.",
                "- Volumetric energy density is diagnostic only, not a sole defect criterion.",
                "- Recalibrate thresholds and weights using CT/metallographic relative density, "
                "melt-pool cross sections, and tensile/fatigue data.",
                "",
                "See `feature_summary.json` for detailed field statistics and "
                "`snapshot_features.csv` for per-snapshot melt-pool features.",
                "",
            ]
        )
        return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Physics-informed screening score for one 3DThesis LPBF layer."
    )
    parser.add_argument(
        "--solidification-dir",
        type=Path,
        required=True,
        help="3DThesis case directory, Data directory, or Final CSV.",
    )
    parser.add_argument(
        "--snapshots-dir",
        type=Path,
        default=None,
        help="Optional 3DThesis snapshots case/Data directory.",
    )
    parser.add_argument(
        "--thermal-history-dir",
        type=Path,
        default=None,
        help="Optional directory containing T_hist CSV files.",
    )
    parser.add_argument(
        "--source-cli",
        type=Path,
        default=None,
        help=(
            "Original slicer CLI file used to read the physical layer spacing. "
            "If omitted, one unambiguous *.cli next to the cases is auto-detected."
        ),
    )
    parser.add_argument(
        "--layer-thickness-um",
        type=float,
        default=None,
        help="Explicit physical powder-layer thickness in micrometres.",
    )
    parser.add_argument(
        "--hatch-spacing-um",
        type=float,
        default=None,
        help="Explicit intended hatch spacing in micrometres.",
    )
    parser.add_argument(
        "--stress-test",
        action="store_true",
        help=(
            "Accept intentional large offsets/parameter differences as stress-test "
            "inputs while preserving physical FAIL criteria and evidence isolation."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "results" / "physics_assessment",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=PROJECT_DIR / "config" / "scoring.yaml",
    )
    parser.add_argument("--layer-id", default="Layer-1")
    return parser


def main(argv: list[str] | None = None) -> dict[str, Any]:
    args = build_parser().parse_args(argv)
    assessment = LayerAssessment(
        solidification_path=args.solidification_dir,
        snapshots_path=args.snapshots_dir,
        output_dir=args.output_dir,
        config_path=args.config,
        layer_id=args.layer_id,
        thermal_history_path=args.thermal_history_dir,
        source_cli=args.source_cli,
        layer_thickness_m=(
            args.layer_thickness_um * 1e-6
            if args.layer_thickness_um is not None
            else None
        ),
        hatch_spacing_m=(
            args.hatch_spacing_um * 1e-6
            if args.hatch_spacing_um is not None
            else None
        ),
        stress_test=args.stress_test,
    )
    print("=" * 68)
    print("3DThesis LPBF Physics-Informed Layer Assessment")
    print("=" * 68)
    print(f"Solidification: {args.solidification_dir}")
    print(f"Snapshots:      {args.snapshots_dir or 'not supplied'}")
    print(f"Engine:         v{ENGINE_VERSION} ({Path(__file__).resolve()})")
    result = assessment.run()
    print("-" * 68)
    print(f"Decision:       {result['decision']}")
    print(f"Process score:  {result['quality_score']:.1f}/100 ({result['continuous_grade']})")
    print(
        f"Data complete:  {result['evidence_completeness_score']:.1f}/100"
    )
    print(
        f"Evidence rep.:  {result.get('evidence_adequacy_score', 0):.1f}/100"
    )
    print(f"Mode:           {result['assessment_mode']}")
    print(
        f"Calibration:    {result.get('model_calibration', {}).get('status', 'UNKNOWN')}"
    )
    top = result.get("metrics", {}).get(
        "top_surface_cumulative_coverage_fraction"
    )
    interface = result.get("metrics", {}).get(
        "interface_cumulative_fusion_fraction"
    )
    print(
        "Top/interface:  "
        + (
            f"{100 * top:.1f}% / {100 * interface:.1f}%"
            if top is not None and interface is not None
            else f"not proven ({result['coverage']['status']})"
        )
    )
    print(f"Output:         {args.output_dir.resolve()}")
    primary = next(
        iter(result.get("problem_diagnosis", [])), None
    )
    if primary:
        print(
            f"Primary issue:  {primary.get('where', 'unknown')} - "
            f"{primary.get('finding', '')}"
        )
    print(f"Read first:     {(args.output_dir / 'Dashboard.html').resolve()}")
    print("=" * 68)
    return result


if __name__ == "__main__":
    main()
