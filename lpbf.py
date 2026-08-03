#!/usr/bin/env python3
"""Single entry point for the LPBF screening pipeline.

    python3 lpbf.py <subcommand> [options]

Subcommands
-----------
    layers      list candidate layers, to choose which one to run
    probe       estimate coverage points and wall time (seconds; do this first)
    layer       full chain for one layer: CLI -> Path -> 3DThesis -> scoring
    sweep       compare several slices; emits a report and an LOF curve
    score       re-score without re-simulating
    report      re-render the report only
    calibrate   fit Depth_Z / Efficiency to measured melt-pool sizes
    uq          uncertainty envelope for the uncalibrated parameters
    test        run the unit tests

Every subcommand supports --help, e.g.

    python3 lpbf.py layer --help

The underlying scripts are unchanged; this file only wires them together and
normalises paths. See docs/PIPELINE.md for the mapping.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "source"
AGENT = SOURCE / "LPBF_Agent"
CALIB = SOURCE / "Test" / "316L" / "C_calib_200W_316L"

# Runtime estimate from the archived two-plane, solid-masked layer-1600 case.
REFERENCE_POINTS = 80_498
REFERENCE_MINUTES = 19.0


def run(command: list[str], cwd: Path) -> int:
    printable = " ".join(str(part) for part in command)
    print(f"$ (cd {cwd.relative_to(ROOT)} && {printable})\n", flush=True)
    return subprocess.run([str(part) for part in command], cwd=cwd).returncode


def cmd_probe(args, extra: list[str]) -> int:
    """Estimate point count before a large section (mid-hull is ~30x the top)."""
    import importlib.util
    import tempfile

    cli_path = Path(args.cli).expanduser().resolve()
    if not cli_path.is_file():
        print(f"CLI file not found: {cli_path}")
        return 2

    # Reuse production case assembly so the estimate uses an identical Domain.
    sys.path.insert(0, str(AGENT))
    spec = importlib.util.spec_from_file_location(
        "sweep_hatch_for_probe", AGENT / "sweep_hatch.py")
    sweep = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = sweep
    spec.loader.exec_module(sweep)

    with tempfile.TemporaryDirectory() as temp:
        path_txt = Path(temp) / "probe_path.txt"
        code = run([
            sys.executable, SOURCE / "CLI_Trans.py", cli_path, path_txt,
            "--layer-z-mm", str(args.layer_z_mm),
            "--whole-layer", "--include-contours", "--recenter-xy",
            "--scan-speed", str(args.scan_speed),
            "--contour-speed", str(args.contour_speed),
        ], cwd=SOURCE)
        if code != 0:
            return code

        timing = sweep.read_path_timing(path_txt)
        case = Path(temp) / "probe_case"
        sweep.build_solidification_case(
            case, SOURCE / "Test" / "316L" / "B_200W_316L", path_txt,
            sweep.powered_bbox(timing), args.margin_mm * 1e-3,
            args.grid_um * 1e-6, args.layer_thickness_um * 1e-6,
            args.z_levels, args.hatch_spacing_um,
        )
        points_file = case / "coverage_target_points.txt"
        result = subprocess.run([
            sys.executable, str(AGENT / "generate_coverage_points.py"),
            "--case-dir", str(case),
            "--domain-file", str(case / "Domain.txt"),
            "--path-file", str(path_txt),
            "--source-cli", str(cli_path),
            "--layer-thickness-um", str(args.layer_thickness_um),
            "--hatch-spacing-um", str(args.hatch_spacing_um),
            "--output", str(points_file),
        ], cwd=AGENT, text=True, capture_output=True)
        if result.returncode != 0:
            sys.stdout.write(result.stdout)
            sys.stderr.write(result.stderr)
            return result.returncode
        for line in result.stdout.splitlines():
            if any(key in line for key in
                   ("Powered segments", "Layer thickness", "Hatch spacing",
                    "Grid:", "Points:", "Z planes")):
                print("  " + line.strip())
        points = sum(1 for _ in points_file.open()) if points_file.is_file() else 0

    print("\n" + "=" * 62)
    if not points:
        print("Could not parse a point count; see the statistics above.")
        return 0
    minutes = REFERENCE_MINUTES * points / REFERENCE_POINTS
    print(f"Coverage points   {points:,}")
    print(f"Case B estimate   ~{minutes:.0f} min"
          f" (linear from layer 1600: {REFERENCE_POINTS:,} pts / {REFERENCE_MINUTES:.0f} min)")
    print("=" * 62)
    if points > 500_000:
        print("\n! Over 500k points. Pick one:")
        print("  1. Restrict the evaluation region (recommended): --region-mm 6."
              " No loss of in-plane resolution; conclusion scoped to that box.")
        print("  2. --grid-um 50: a quarter of the points, but the unmelted"
              " ridge is only ~12 um wide -- rough look only.")
        print("  3. Accept the runtime and run overnight.")
    else:
        print("\nOK -- size is normal; go ahead with lpbf.py layer")
    return 0


def cmd_layers(args, extra: list[str]) -> int:
    """List candidate layers. Infill count tracks section area and runtime."""

    cli_path = Path(args.cli).expanduser().resolve()
    if not cli_path.is_file():
        print(f"CLI file not found: {cli_path}")
        return 2

    command = [sys.executable, str(SOURCE / "CLI_Trans.py"),
               str(cli_path), "/dev/null", "--list-layers"]
    result = subprocess.run(command, cwd=SOURCE, text=True, capture_output=True)
    if result.returncode != 0:
        sys.stderr.write(result.stderr)
        return result.returncode

    rows = []
    for line in result.stdout.splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) == 4 and parts[0].isdigit():
            rows.append((int(parts[0]), float(parts[1]), int(parts[2]), int(parts[3])))
    if not rows:
        sys.stdout.write(result.stdout)
        return 0

    print(f"{len(rows)} non-empty layers | z from {rows[0][1]:.3f} to {rows[-1][1]:.3f} mm\n")
    reference = max(r[2] for r in rows)
    print(f"{'layer':>6} {'z (mm)':>10} {'infill':>8} {'contour':>8}  relative size")
    print("-" * 58)
    chosen = sorted(rows, key=lambda r: -r[2])[:args.top] if args.by_size else \
        rows[:: max(len(rows) // args.top, 1)][:args.top]
    for number, z, hatches, contours in sorted(chosen):
        bar = "█" * max(1, round(20 * hatches / reference))
        print(f"{number:>6} {z:>10.3f} {hatches:>8} {contours:>8}  {bar}")
    print("\nThen estimate cost:  python3 lpbf.py probe --cli <cli> --layer-z-mm <z>")
    return 0


def cmd_layer(args, extra: list[str]) -> int:
    """Full chain for one layer -- sweep with a single --cli."""
    return run([
        sys.executable, AGENT / "sweep_hatch.py",
        "--cli", f"{args.hatch_spacing_um:g}={Path(args.cli).resolve()}",
        "--layer-z-mm", str(args.layer_z_mm),
        "--workdir", Path(args.workdir).resolve(),
        "--output-dir", Path(args.output_dir).resolve(),
        *resolve_paths(extra),
    ], cwd=AGENT)


PATH_OPTIONS = (
    "--workdir", "--output-dir", "--thesis-bin", "--binary",
    "--template-a", "--template-b", "--targets", "--output", "--config",
    "--solidification-dir", "--snapshots-dir", "--source-cli",
    "--thermal-history-dir",
)


def resolve_paths(extra: list[str]) -> list[str]:
    """Resolve relative paths in pass-through arguments against the caller's cwd.

    Subcommands cd into the underlying script's directory before running, so
    without this a relative --cli path would resolve against the wrong root.
    """
    resolved: list[str] = []
    index = 0
    while index < len(extra):
        token = extra[index]
        if token == "--cli" and index + 1 < len(extra) and "=" in extra[index + 1]:
            key, _, value = extra[index + 1].partition("=")
            resolved += [token, f"{key}={Path(value).expanduser().resolve()}"]
            index += 2
        elif token in PATH_OPTIONS and index + 1 < len(extra):
            resolved += [token, str(Path(extra[index + 1]).expanduser().resolve())]
            index += 2
        elif token.startswith("--cli=") and "=" in token[6:]:
            key, _, value = token[6:].partition("=")
            resolved.append(f"--cli={key}={Path(value).expanduser().resolve()}")
            index += 1
        else:
            resolved.append(token)
            index += 1
    return resolved


def cmd_sweep(args, extra: list[str]) -> int:
    return run([sys.executable, AGENT / "sweep_hatch.py", *resolve_paths(extra)],
               cwd=AGENT)


def cmd_score(args, extra: list[str]) -> int:
    return run([sys.executable, AGENT / "run_all.py", *resolve_paths(extra)], cwd=AGENT)


def cmd_report(args, extra: list[str]) -> int:
    forwarded = resolve_paths(extra)
    if forwarded and not forwarded[0].startswith("-"):
        forwarded[0] = str(Path(forwarded[0]).expanduser().resolve())
    return run([sys.executable, AGENT / "report_summary.py", *forwarded], cwd=AGENT)


def cmd_calibrate(args, extra: list[str]) -> int:
    return run([sys.executable, CALIB / "calibrate_beam.py", *resolve_paths(extra)], cwd=CALIB)


def cmd_uq(args, extra: list[str]) -> int:
    forwarded = resolve_paths(extra)
    if forwarded and not forwarded[0].startswith("-"):
        forwarded[0] = str(Path(forwarded[0]).expanduser().resolve())
    return run([sys.executable, AGENT / "uq_envelope.py", *forwarded], cwd=AGENT)


def cmd_test(args, extra: list[str]) -> int:
    import os
    env_command = [
        sys.executable, "-m", "unittest",
        "tests.test_physics_scorer", "tests.test_coverage_generator", "-v",
    ]
    print(f"$ (cd source/LPBF_Agent && PYTHONPATH=. python3 "
          f"{' '.join(env_command[1:])})\n", flush=True)
    environment = dict(os.environ, PYTHONPATH=str(AGENT))
    return subprocess.run(env_command, cwd=AGENT, env=environment).returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lpbf.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    probe = subparsers.add_parser(
        "probe", help="estimate coverage points and wall time")
    probe.add_argument("--cli", required=True)
    probe.add_argument("--layer-z-mm", type=float, required=True)
    probe.add_argument("--hatch-spacing-um", type=float, default=100.0)
    probe.add_argument("--grid-um", type=float, default=25.0)
    probe.add_argument("--layer-thickness-um", type=float, default=30.0)
    probe.add_argument("--z-levels", type=int, default=3)
    probe.add_argument("--margin-mm", type=float, default=1.5)
    probe.add_argument("--scan-speed", type=float, default=0.8)
    probe.add_argument("--contour-speed", type=float, default=0.6)
    probe.set_defaults(handler=cmd_probe)

    layers = subparsers.add_parser("layers", help="list candidate layers")
    layers.add_argument("--cli", required=True)
    layers.add_argument("--top", type=int, default=20, help="how many layers to show")
    layers.add_argument("--by-size", action="store_true",
                        help="sort by section size (default: even sampling by height)")
    layers.set_defaults(handler=cmd_layers)

    layer = subparsers.add_parser("layer", help="full chain for one layer")
    layer.add_argument("--cli", required=True)
    layer.add_argument("--layer-z-mm", type=float, required=True)
    layer.add_argument("--hatch-spacing-um", type=float, default=100.0)
    layer.add_argument("--workdir", default=str(SOURCE / "Test" / "316L" / "runs"))
    layer.add_argument("--output-dir", default=str(AGENT / "results" / "layer"))
    layer.set_defaults(handler=cmd_layer)

    for name, handler, description in (
        ("sweep", cmd_sweep, "compare several slices"),
        ("score", cmd_score, "re-score without re-simulating"),
        ("report", cmd_report, "re-render the report only"),
        ("calibrate", cmd_calibrate, "fit Depth_Z / Efficiency"),
        ("uq", cmd_uq, "uncertainty envelope for uncalibrated parameters"),
    ):
        sub = subparsers.add_parser(
            name, help=description, add_help=False,
            description=f"{description}. Arguments are passed through verbatim to the "
                        f"underlying script, including -h.")
        sub.set_defaults(handler=handler)

    tests = subparsers.add_parser("test", help="run the unit tests")
    tests.set_defaults(handler=cmd_test)

    return parser


def main() -> int:
    parser = build_parser()
    known, extra = parser.parse_known_args()
    return known.handler(known, extra)


if __name__ == "__main__":
    raise SystemExit(main())
