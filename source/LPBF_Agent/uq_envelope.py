#!/usr/bin/env python3
"""Propagate uncalibrated beam parameters and apply literature constraints.

Efficiency and effective source depth require calibration against single-track
metallography. Until then, the analysis reports dimensionless h/W, t/D, and
Tang LOF envelopes over a parameter grid and excludes combinations that cannot
reproduce a published dense-process anchor.

Temperature rise is linear in Efficiency in 3DThesis, so each Depth_Z value is
simulated once and the Efficiency axis is recovered by changing the melt
threshold: ``dT_ref >= (T_L - T_0) * eff_ref / eff``.

Usage:
    python3 uq_envelope.py ../Test/316L/A_200W_316L --output-dir results/uq
    python3 uq_envelope.py <case> --depth-um 20,27.5,35,45,60 \\
        --efficiency 0.25,0.35,0.45,0.55,0.70 --snapshots 5
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lpbf_score.scorer import component_for_geometry, infer_grid_step  # noqa: E402

DEFAULT_BIN = Path("/Users/joe/Projects/MIT_Project/3DThesis/build/bin/3DThesis")


def patch_case(src: Path, dst: Path, depth_z_m: float, eff_ref: float,
               n_snapshots: int | None) -> None:
    dst.mkdir(parents=True, exist_ok=True)
    (dst / "Data").mkdir(exist_ok=True)
    for f in src.glob("*.txt"):
        if f.name.endswith(".bak"):
            continue
        shutil.copy(f, dst / f.name)

    beam = (dst / "Beam.txt").read_text(encoding="utf-8")
    beam = re.sub(r"(?m)^(\s*Depth_Z\s+)\S+", rf"\g<1>{depth_z_m:.6e}", beam)
    beam = re.sub(r"(?m)^(\s*Efficiency\s+)\S+", rf"\g<1>{eff_ref:.6f}", beam)
    (dst / "Beam.txt").write_text(beam, encoding="utf-8")

    if n_snapshots:
        mode = (dst / "Mode.txt").read_text(encoding="utf-8")
        match = re.search(r"(?m)^(\s*ScanFracs\s+)(\S+)", mode)
        if match:
            values = match.group(2).split(",")
            step = max(1, len(values) // n_snapshots)
            kept = ",".join(values[::step][:n_snapshots])
            mode = mode.replace(match.group(2), kept)
            (dst / "Mode.txt").write_text(mode, encoding="utf-8")


def run_case(case: Path, binary: Path) -> None:
    result = subprocess.run([str(binary), "ParamInput.txt"], cwd=case,
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"3DThesis failed ({case}):\n{result.stdout[-800:]}")


def load_fields(case: Path, t0: float, keep_above: float):
    """Load snapshot points above the least restrictive melt threshold."""
    fields = []
    for f in sorted((case / "Data").glob("snapshots.Snapshot.*.csv")):
        df = pd.read_csv(f, usecols=["x", "y", "z", "T"])
        df["dT"] = df["T"] - t0
        df = df[df["dT"] >= keep_above]
        if len(df):
            fields.append(df.reset_index(drop=True))
    if not fields:
        raise RuntimeError(f"No points in {case} exceed the threshold; check power and efficiency.")
    return fields


def pool_geometry(molten: pd.DataFrame, dx: float, dy: float, dz: float):
    """Measure the largest connected pool with PCA-aligned XY extents."""
    if len(molten) < 3:
        return None, None
    pool, _ = component_for_geometry(molten, dx, dy, dz)
    xy = pool[["x", "y"]].to_numpy(float)
    if len(xy) < 3:
        return None, None
    centred = xy - xy.mean(axis=0)
    _, vectors = np.linalg.eigh(np.cov(centred, rowvar=False))
    extents = np.ptp(centred @ vectors, axis=0) + max(dx, dy)
    width = float(min(extents))
    depth = float(pool["z"].max() - pool["z"].min() + dz)
    return width, depth


def sweep_efficiency(fields, t0: float, t_liq: float, eff_ref: float,
                     efficiencies) -> dict[float, tuple[float, float]]:
    """Recover an Efficiency axis as ``{efficiency: (W_p10, D_p10)}``."""
    grids = [
        (infer_grid_step(f["x"], 1e-5),
         infer_grid_step(f["y"], 1e-5),
         infer_grid_step(f["z"], 5e-6))
        for f in fields
    ]
    out: dict[float, tuple[float, float]] = {}
    for eff in efficiencies:
        threshold = (t_liq - t0) * eff_ref / eff
        widths, depths = [], []
        for field, (dx, dy, dz) in zip(fields, grids):
            molten = field[field["dT"] >= threshold]
            w, d = pool_geometry(molten, dx, dy, dz)
            if w is not None:
                widths.append(w)
                depths.append(d)
        if len(widths) >= 2:
            out[eff] = (float(np.quantile(widths, 0.10)),
                        float(np.quantile(depths, 0.10)))
        else:
            out[eff] = (float("nan"), float("nan"))
    return out


def lof_index(hatch_m: float, layer_m: float, w: float, d: float) -> float:
    """Tang, Pistorius & Beuth, Addit. Manuf. 2017"""
    if not (w > 0 and d > 0):
        return float("inf")
    return (hatch_m / w) ** 2 + (layer_m / d) ** 2


def main() -> int:
    p = argparse.ArgumentParser(description="Uncalibrated-parameter envelope with literature constraints.")
    p.add_argument("case", help="Snapshot case directory, for example A_200W_316L.")
    p.add_argument("--binary", default=str(DEFAULT_BIN))
    p.add_argument("--output-dir", default="results/uq")
    p.add_argument("--depth-um", default="20,27.5,35,45,60",
                   help="Comma-separated Depth_Z values in µm.")
    p.add_argument("--efficiency", default="0.25,0.35,0.45,0.55,0.70",
                   help="Comma-separated Efficiency values recovered by linear rescaling.")
    p.add_argument("--eff-ref", type=float, default=0.35, help="Efficiency used in each simulation.")
    p.add_argument("--snapshots", type=int, default=5, help="Snapshots retained per simulation.")
    p.add_argument("--hatch-um", type=float, default=100.0, help="Study hatch spacing in µm.")
    p.add_argument("--layer-um", type=float, default=30.0, help="Layer thickness in µm.")
    p.add_argument("--anchor-hatch-um", type=float, default=120.0,
                   help="Hatch spacing of the published dense-build anchor.")
    p.add_argument("--anchor-tol", type=float, default=1.0,
                   help="Maximum I_LOF allowed at the literature anchor.")
    p.add_argument("--t0", type=float, default=298.15)
    p.add_argument("--t-liq", type=float, default=1708.0)
    p.add_argument("--keep-runs", action="store_true", help="Retain intermediate cases.")
    args = p.parse_args()

    case = Path(args.case).expanduser().resolve()
    out = Path(args.output_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    depths = [float(v) * 1e-6 for v in args.depth_um.split(",")]
    effs = [float(v) for v in args.efficiency.split(",")]
    hatch, layer = args.hatch_um * 1e-6, args.layer_um * 1e-6
    anchor = args.anchor_hatch_um * 1e-6

    # Retain every point that can melt at any requested efficiency.
    keep_above = (args.t_liq - args.t0) * args.eff_ref / max(effs)

    rows = []
    workdir = out / "_runs"
    for depth in depths:
        run = workdir / f"Dz{depth*1e6:g}"
        print(f"[run] Depth_Z = {depth*1e6:g} um ...", flush=True)
        patch_case(case, run, depth, args.eff_ref, args.snapshots)
        run_case(run, Path(args.binary))
        fields = load_fields(run, args.t0, keep_above)
        table = sweep_efficiency(fields, args.t0, args.t_liq, args.eff_ref, effs)
        for eff, (w, d) in table.items():
            if not (w > 0 and d > 0):
                continue
            rows.append({
                "depth_z_um": depth * 1e6,
                "efficiency": eff,
                "W_p10_um": w * 1e6,
                "D_p10_um": d * 1e6,
                "t_over_D": layer / d,
                "h_over_W": hatch / w,
                "D_over_W": d / w,
                "I_LOF_study": lof_index(hatch, layer, w, d),
                "I_LOF_anchor": lof_index(anchor, layer, w, d),
                "h_critical_um": w * math.sqrt(max(1 - (layer / d) ** 2, 0.0)) * 1e6,
            })

    if not args.keep_runs and workdir.exists():
        shutil.rmtree(workdir, ignore_errors=True)

    frame = pd.DataFrame(rows)

    # Feasible cases reproduce the dense LOF anchor and remain wider than deep.
    frame["feasible_lof"] = frame["I_LOF_anchor"] <= args.anchor_tol
    frame["feasible_mode"] = frame["D_over_W"] <= 1.0
    frame["feasible"] = frame["feasible_lof"] & frame["feasible_mode"]
    feasible = frame[frame["feasible"]]
    frame.to_csv(out / "uq_grid.csv", index=False)

    lines: list[str] = []
    w = lines.append
    w("# Uncalibrated-Parameter Envelope and Literature Constraints")
    w("")
    w(f"- Case: `{case.name}`")
    w(f"- Grid: Depth_Z in {{{args.depth_um}}} µm × Efficiency in "
      f"{{{args.efficiency}}} ({len(frame)} combinations)")
    w(f"- Study condition: hatch {args.hatch_um:g} µm; layer {args.layer_um:g} µm")
    w(f"- Literature anchor: dense at hatch {args.anchor_hatch_um:g} µm; "
      f"required I_LOF <= {args.anchor_tol:g}")
    w("")
    w("> The Efficiency axis requires no additional simulations: temperature "
      "rise is linear in Efficiency. Each Depth_Z value is simulated once.")
    w("")

    w("## Full dimensionless envelope")
    w("")
    w("| Depth_Z (µm) | Efficiency | W_p10 (µm) | D_p10 (µm) | t/D | h/W | D/W | **I_LOF** | critical h (µm) |")
    w("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for _, r in frame.iterrows():
        w(f"| {r.depth_z_um:g} | {r.efficiency:.2f} | {r.W_p10_um:.1f} | {r.D_p10_um:.1f} "
          f"| {r.t_over_D:.3f} | {r.h_over_W:.3f} | {r.D_over_W:.3f} "
          f"| **{r.I_LOF_study:.3f}** | {r.h_critical_um:.0f} |")
    w("")

    w("## Conclusions valid across the entire envelope")
    w("")
    lo, hi = frame["I_LOF_study"].min(), frame["I_LOF_study"].max()
    w(f"- I_LOF ∈ **[{lo:.3f}, {hi:.3f}]**")
    if lo > 1.0:
        w(f"  - **Every case exceeds 1.0: hatch {args.hatch_um:g} µm has "
          "insufficient geometric overlap independent of calibration.**")
    elif hi < 1.0:
        w("  - Every case is below 1.0: overlap is sufficient independent of calibration.")
    else:
        w("  - The envelope crosses 1.0: **calibration is required before a verdict.**")
    w(f"- h/W ∈ [{frame['h_over_W'].min():.3f}, {frame['h_over_W'].max():.3f}]; "
      "required h/W <= sqrt(1-(t/D)^2)")
    w(f"- Critical hatch spacing ∈ **[{frame['h_critical_um'].min():.0f}, "
      f"{frame['h_critical_um'].max():.0f}] µm**")
    w(f"- D/W ∈ [{frame['D_over_W'].min():.3f}, {frame['D_over_W'].max():.3f}]"
      " (<0.5 conduction-like; >0.8 keyhole risk)")
    w("")

    w("## Feasible region after applying literature constraints")
    w("")
    if len(feasible):
        w(f"{len(feasible)}/{len(frame)} combinations satisfy both constraints:")
        w("")
        w("| Depth_Z (µm) | Efficiency | W_p10 | D_p10 | study I_LOF | anchor I_LOF |")
        w("|---:|---:|---:|---:|---:|---:|")
        for _, r in feasible.iterrows():
            w(f"| {r.depth_z_um:g} | {r.efficiency:.2f} | {r.W_p10_um:.1f} "
              f"| {r.D_p10_um:.1f} | {r.I_LOF_study:.3f} | {r.I_LOF_anchor:.3f} |")
        w("")
        w(f"Within the feasible region, I_LOF ∈ [{feasible['I_LOF_study'].min():.3f}, "
          f"{feasible['I_LOF_study'].max():.3f}]; "
          f"critical hatch ∈ [{feasible['h_critical_um'].min():.0f}, "
          f"{feasible['h_critical_um'].max():.0f}] µm")
    else:
        w("**The feasible region is empty.** No grid combination reproduces the "
          "literature anchor. At least one explanation must hold:")
        w("")
        w("1. The true parameters lie outside the scanned range.")
        w("2. The conduction model systematically under-predicts the pool.")
        w("3. The Tang criterion is conservative here, or the anchor tolerates some LOF.")
        w("")
        need = frame.assign(
            need_W_um=lambda x: anchor / np.sqrt(
                np.maximum(1 - (layer / (x.D_p10_um * 1e-6)) ** 2, 1e-9)) * 1e6)
        w("Width increase required to reproduce the anchor:")
        w("")
        w("| Depth_Z (µm) | Efficiency | current W (µm) | anchor W minimum (µm) | gap |")
        w("|---:|---:|---:|---:|---:|")
        for _, r in need.sort_values("D_p10_um").iterrows():
            w(f"| {r.depth_z_um:g} | {r.efficiency:.2f} | {r.W_p10_um:.1f} "
              f"| {r.need_W_um:.0f} | {100*(r.need_W_um/r.W_p10_um-1):+.0f}% |")
    w("")

    w("## How to use the envelope")
    w("")
    w("- Report dimensionless h/W, t/D, and I_LOF so later calibration can be substituted directly.")
    w("- Report intervals rather than an unsupported single critical-hatch value.")
    w("- An empty feasible region quantifies a testable model-to-literature discrepancy.")
    w("- After calibration, narrow `--depth-um` and `--efficiency` around the fitted values.")
    w("")
    w("---")
    w("")
    w("**Sources:** Tang, Pistorius, and Beuth, *Addit. Manuf.* 2017 for "
      "the LOF criterion; Trapp et al., *Appl. Mater. Today* 9 (2017) 341 "
      "for absorptivity variation; Coleman et al., *Addit. Manuf.* 95 "
      "(2024) 104531 for single-track calibration of source size and absorptivity.")

    (out / "UQ_Report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\n[wrote] {out/'uq_grid.csv'}\n[wrote] {out/'UQ_Report.md'}")
    print(f"\nI_LOF envelope [{lo:.3f}, {hi:.3f}]; feasible {len(feasible)}/{len(frame)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        raise SystemExit(1)
