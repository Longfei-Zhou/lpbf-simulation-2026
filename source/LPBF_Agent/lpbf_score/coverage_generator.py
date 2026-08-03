"""Generate a compact 3DThesis Custom-domain coverage point file."""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd

from .scorer import (
    ProcessParameters,
    build_coverage_target_points,
    find_key,
    infer_hatch_spacing,
    parse_cli_layer_thickness,
    parse_grouped_text,
    parse_path_file,
    powered_segments,
)


def _positive(value: float, label: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{label} must be a finite positive number.")
    return number


def _domain_value(parsed: dict[str, Any], suffix: str) -> float | None:
    value = find_key(parsed, suffix)
    return float(value) if value is not None else None


def _grid_resolution(parsed: dict[str, Any], axis: str) -> float:
    resolution = _domain_value(parsed, f"{axis}.Res")
    if resolution is not None:
        return _positive(resolution, f"{axis}.Res")
    minimum = _domain_value(parsed, f"{axis}.Min")
    maximum = _domain_value(parsed, f"{axis}.Max")
    count = _domain_value(parsed, f"{axis}.Num")
    if None not in (minimum, maximum, count) and int(count) > 1:
        return _positive(
            (float(maximum) - float(minimum)) / (int(count) - 1),
            f"{axis} derived resolution",
        )
    raise ValueError(
        f"The reference Domain must define {axis}.Res, or {axis}.Min/Max/Num."
    )


def _validate_regular_xy_domain(parsed: dict[str, Any], path: Path) -> None:
    missing = [
        key
        for key in ("X.Min", "X.Max", "Y.Min", "Y.Max")
        if find_key(parsed, key) is None
    ]
    if missing:
        raise ValueError(
            f"{path} is not a regular XY Domain; missing {', '.join(missing)}. "
            "If the Solidification Domain is already Custom, pass the Snapshot "
            "case's original regular Domain with --domain-file."
        )


def _atomic_write_points(points: pd.DataFrame, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_name = temporary.name
            points.to_csv(
                temporary,
                sep=" ",
                header=False,
                index=False,
                float_format="%.9f",
                lineterminator="\n",
            )
            temporary.flush()
            os.fsync(temporary.fileno())
        os.replace(temporary_name, output)
    except Exception:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)
        raise


def generate_coverage_points(
    *,
    case_dir: Path,
    domain_file: Path | None = None,
    path_file: Path | None = None,
    source_cli: Path | None = None,
    layer_thickness_um: float | None = None,
    hatch_spacing_um: float | None = None,
    stress_test: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Generate top/interface target points without requiring simulation CSVs."""
    case_dir = case_dir.expanduser().resolve()
    domain_file = (
        domain_file.expanduser().resolve()
        if domain_file is not None
        else case_dir / "Domain.txt"
    )
    path_file = (
        path_file.expanduser().resolve()
        if path_file is not None
        else case_dir / "Path.txt"
    )
    if not case_dir.is_dir():
        raise FileNotFoundError(f"Case directory does not exist: {case_dir}")
    if not domain_file.is_file():
        raise FileNotFoundError(f"Reference Domain file does not exist: {domain_file}")
    if not path_file.is_file():
        raise FileNotFoundError(f"Path file does not exist: {path_file}")

    parsed_domain = parse_grouped_text(domain_file)
    _validate_regular_xy_domain(parsed_domain, domain_file)
    grid_x_m = _grid_resolution(parsed_domain, "X")
    grid_y_m = _grid_resolution(parsed_domain, "Y")
    grid_z_m = _domain_value(parsed_domain, "Z.Res")
    if grid_z_m is not None:
        grid_z_m = _positive(grid_z_m, "Z.Res")

    path_frame = parse_path_file(path_file)
    if path_frame is None:
        raise ValueError(f"Could not parse a valid 3DThesis Path file: {path_file}")
    diagnostics = dict(path_frame.attrs.get("path_parse_diagnostics", {}))
    if diagnostics.get("extreme_concatenated_line_format_detected") and not stress_test:
        raise ValueError(
            "Path.txt contains concatenated trailing records. Regenerate it with "
            "CLI_Trans.py, or pass --stress-test to deliberately mirror 3DThesis "
            "by using the first record on each physical line."
        )
    segments = powered_segments(path_frame)
    if not len(segments):
        raise ValueError("Path.txt contains no non-zero powered scan segments.")

    if layer_thickness_um is not None:
        layer_thickness_m = _positive(
            layer_thickness_um, "--layer-thickness-um"
        ) * 1e-6
        layer_source = "--layer-thickness-um"
    else:
        cli_result = parse_cli_layer_thickness(
            source_cli.expanduser().resolve() if source_cli is not None else None
        )
        if cli_result is None:
            raise ValueError(
                "Provide --source-cli or --layer-thickness-um. Domain Z.Res is "
                "numerical resolution and is not accepted as physical layer thickness."
            )
        layer_thickness_m = _positive(
            cli_result["layer_thickness_m"], "CLI layer thickness"
        )
        layer_source = str(cli_result["path"])

    if hatch_spacing_um is not None:
        hatch_spacing_m = _positive(
            hatch_spacing_um, "--hatch-spacing-um"
        ) * 1e-6
        hatch_source = "--hatch-spacing-um"
    else:
        inferred_hatch = infer_hatch_spacing(path_frame)
        if inferred_hatch is None:
            raise ValueError(
                "Hatch spacing could not be inferred reliably; pass "
                "--hatch-spacing-um explicitly."
            )
        hatch_spacing_m = _positive(inferred_hatch, "inferred hatch spacing")
        hatch_source = str(path_file)

    parameters = ProcessParameters(
        liquidus_k=0.0,
        liquidus_source="not required for target generation",
        layer_thickness_m=layer_thickness_m,
        layer_thickness_source=layer_source,
        hatch_spacing_m=hatch_spacing_m,
        hatch_spacing_source=hatch_source,
        grid_x_m=grid_x_m,
        grid_y_m=grid_y_m,
        grid_z_m=grid_z_m or layer_thickness_m,
        power_w=None,
        efficiency=None,
        scan_speed_m_s=None,
    )
    points = build_coverage_target_points(
        parsed_domain,
        path_frame,
        parameters,
        source_cli=(
            source_cli.expanduser().resolve()
            if source_cli is not None
            else None
        ),
    )
    if points.empty:
        raise ValueError(
            "No coverage points were generated. Check that the Path lies inside "
            "the reference Domain and that all coordinate units are correct."
        )
    z_counts = {
        f"{float(z):.9f}": int(count)
        for z, count in points.groupby("z_mm", sort=True).size().items()
    }
    summary = {
        "case_dir": str(case_dir),
        "path_file": str(path_file),
        "domain_file": str(domain_file),
        "point_count": int(len(points)),
        "powered_segment_count": int(len(segments)),
        "layer_thickness_um": layer_thickness_m * 1e6,
        "layer_thickness_source": layer_source,
        "hatch_spacing_um": hatch_spacing_m * 1e6,
        "hatch_spacing_source": hatch_source,
        "grid_x_um": grid_x_m * 1e6,
        "grid_y_um": grid_y_m * 1e6,
        "corridor_half_width_um": (
            0.5 * max(hatch_spacing_m, grid_x_m, grid_y_m) * 1e6
        ),
        "z_plane_counts_mm": z_counts,
        "x_bounds_mm": [
            float(points["x_mm"].min()),
            float(points["x_mm"].max()),
        ],
        "y_bounds_mm": [
            float(points["y_mm"].min()),
            float(points["y_mm"].max()),
        ],
        "path_parse_diagnostics": diagnostics,
        "cli_solid_mask": points.attrs.get("cli_solid_mask", {}),
    }
    return points, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a 3DThesis Custom coverage_target_points.txt before "
            "running Solidification; no Final CSV is required."
        )
    )
    parser.add_argument(
        "--case-dir",
        type=Path,
        required=True,
        help="Current layer case directory; Path.txt is read here by default.",
    )
    parser.add_argument(
        "--domain-file",
        type=Path,
        default=None,
        help=(
            "Regular XYZ reference Domain.txt. Required when case-dir/Domain.txt "
            "has already been replaced by a Custom Domain."
        ),
    )
    parser.add_argument(
        "--path-file",
        type=Path,
        default=None,
        help="Override the current layer Path.txt.",
    )
    parser.add_argument(
        "--source-cli",
        type=Path,
        default=None,
        help="Original ASCII CLI used to derive physical layer thickness.",
    )
    parser.add_argument(
        "--layer-thickness-um",
        type=float,
        default=None,
        help="Explicit physical layer thickness in micrometres; overrides CLI.",
    )
    parser.add_argument(
        "--hatch-spacing-um",
        type=float,
        default=None,
        help="Explicit intended hatch spacing in micrometres.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Output point file. Defaults to case-dir/coverage_target_points.txt."
        ),
    )
    parser.add_argument(
        "--stress-test",
        action="store_true",
        help="Allow concatenated Path lines and mirror 3DThesis first-record parsing.",
    )
    return parser


def main(argv: list[str] | None = None) -> dict[str, Any]:
    args = build_parser().parse_args(argv)
    output = (
        args.output.expanduser().resolve()
        if args.output is not None
        else args.case_dir.expanduser().resolve() / "coverage_target_points.txt"
    )
    points, summary = generate_coverage_points(
        case_dir=args.case_dir,
        domain_file=args.domain_file,
        path_file=args.path_file,
        source_cli=args.source_cli,
        layer_thickness_um=args.layer_thickness_um,
        hatch_spacing_um=args.hatch_spacing_um,
        stress_test=args.stress_test,
    )
    _atomic_write_points(points, output)
    summary["output"] = str(output)

    print("=" * 68)
    print("3DThesis Custom Coverage Point Generator")
    print("=" * 68)
    print(f"Path:             {summary['path_file']}")
    print(f"Reference Domain: {summary['domain_file']}")
    print(f"Powered segments: {summary['powered_segment_count']}")
    print(f"Layer thickness:  {summary['layer_thickness_um']:.3f} um")
    print(f"Hatch spacing:    {summary['hatch_spacing_um']:.3f} um")
    print(f"Grid:             {summary['grid_x_um']:.3f} x {summary['grid_y_um']:.3f} um")
    print(f"Points:           {summary['point_count']}")
    print(f"Z planes (mm):    {summary['z_plane_counts_mm']}")
    print(f"Created:          {output}")
    print("-" * 68)
    print("Use in Domain.txt:")
    print("Custom")
    print("{")
    print(f"    File {output.name}")
    print("}")
    print("=" * 68)
    return summary
