#!/usr/bin/env python3
"""Read melt-pool depth and width for comparison with metallography.

Usage:
    python3 read_meltpool.py [Data/calib.Solidification.Final.csv]

Why ``MP_depth`` is not used directly:
    Melt.cpp converts a Z cell count using ``xres`` rather than ``zres``:
        const double depth = sim.domain.xres * (*std::max_element(...));
    On an anisotropic grid this makes MP_depth wrong by xres/zres. Here the
    ratio is 5.0e-6/2.5e-6, so the reported value is twice the true depth.

    The ``depth`` column uses the correct Z spacing and sub-cell interpolation.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT = Path(__file__).parent / "Data" / "calib.Solidification.Final.csv"
ZRES = 2.5e-6


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    if not path.is_file():
        sys.exit(f"File not found: {path}")

    df = pd.read_csv(path)
    melted = df[df["numMelt"] > 0] if "numMelt" in df else df
    if not len(melted):
        sys.exit("No points melted; check power, speed, and T_L.")

    print(f"File: {path.name}    Melted points: {len(melted)}\n")

    # Cross-check depth using all available representations.
    depth_col = melted["depth"].max() if "depth" in melted else np.nan
    depth_geo = -melted["z"].min()

    print("Melt depth")
    print(f"  depth column (recommended, sub-cell) {depth_col * 1e6:6.1f} um")
    print(f"  geometric minimum molten z           {depth_geo * 1e6:6.1f} um"
          f"   [quantized by zres={ZRES * 1e6:.1f} um]")
    if "MP_depth" in melted:
        mp = melted["MP_depth"].max()
        print(f"  MP_depth column (known bug; do not use) {mp * 1e6:6.1f} um"
              f"   -> corrected {mp * 1e6 / 2:.1f} um")
    if np.isfinite(depth_col) and abs(depth_col - depth_geo) > 2 * ZRES:
        print("  !! depth differs from geometry by more than two cells; check Tracking Surface.")

    top = melted[np.isclose(melted["z"], 0.0)]
    print("\nMelt width")
    if "MP_width" in melted:
        print(f"  MP_width in scan coordinates          {melted['MP_width'].max() * 1e6:6.1f} um")
    if len(top):
        print(f"  top-surface Y span across all tracks  {np.ptp(top['y']) * 1e6:6.1f} um")

    # Exclude track start and stop transients from the steady-state summary.
    steady = melted[(np.abs(melted["y"]) < 1e-9)
                    & (melted["x"] > 0.4e-3) & (melted["x"] < 1.0e-3)]
    if len(steady):
        print(f"\nSteady region (center track y=0, x=0.4–1.0 mm, {len(steady)} points)")
        print(f"  depth       {-steady['z'].min() * 1e6:6.1f} um")
        print(f"  median G    {steady['G'].median():.2e} K/m")
        print(f"  median V    {steady['V'].median():6.3f} m/s")
        print(f"  median dT/dt {steady['dTdt'].median():.2e} K/s")

    print("\nFor metallography, use the same datum as the simulation: z=0 is the solid surface")
    print("and depth is measured from that original surface to the lowest pool point.")


if __name__ == "__main__":
    main()
