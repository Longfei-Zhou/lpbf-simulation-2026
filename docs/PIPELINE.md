# Pipeline

Everything runs through `lpbf.py`. The underlying scripts are unchanged and can
still be called directly — the mapping is in the last section.

---

## Data flow

```
                    3dbenchy.stl
                         │  Netfabb slicing (only manual step — see SLICING.md)
                         ▼
              <h>µm_3dbenchy.cli
                         │  CLI_Trans.py    select layer, add timing, add contours
                         ▼
                     Path.txt
                         │
        ┌────────────────┼────────────────┐
        ▼                ▼                ▼
   case A            case B          (case C, independent)
   snapshots       solidification     calibration
   Tracking None   Tracking None      Tracking Surface
   10/10/5 µm      custom points      5/5/2.5 µm
   20 fields         ▲                single track vs. micrographs
        │            │ coverage_target_points.txt
        │            │ ← generate_coverage_points.py
        │            │   (needs case B's Domain.txt to exist first)
        └──────┬─────┘
               ▼  lpbf_score
        assessment.json + report + plot
```

**Ordering constraint.** The coverage-point generator reads the sampling grid and
z-plane positions from case B's regular `X/Y/Z` blocks via `--domain-file`, so
`Domain.txt` must be written **before** the points are generated. `sweep_hatch.py`
assembles in that order; do not reverse it by hand.

---

## Commands

```bash
python3 lpbf.py <subcommand> [options]
```

| Subcommand | Purpose | When | Runtime |
|---|---|---|---|
| `layers` | List candidate layers | Choosing which layer to run | seconds |
| `probe` | Estimate coverage points and wall time | **First thing after slicing a new layer** | seconds |
| `layer` | Full chain for one layer | Single-layer assessment | ~40 min |
| `sweep` | Compare several slices | Hatch or parameter sweep | N × 40 min |
| `score` | Re-score only | After changing scoring code | minutes |
| `report` | Re-render report only | After changing the template | seconds |
| `calibrate` | Fit `Depth_Z` / `Efficiency` | Once measured pool sizes exist | hours |
| `uq` | Uncertainty envelope | **Before** calibration | minutes |
| `test` | 27 unit tests | After touching any scoring code | seconds |

`layer` is `sweep` with a single `--cli`. They share the same case-assembly code,
so a single-layer run and a sweep can never diverge through configuration drift.

---

## Parameter files: three classes

Remembering this classification avoids every common mistake.

### 1. Never edited — copied from templates automatically

| File | Contains | Why it never changes |
|---|---|---|
| `Material.txt` | T_L, k, c, ρ | Material is 316L regardless of layer |
| `Beam.txt` | Power, spot, `Depth_Z`, `Efficiency` | Machine and laser are fixed |
| `Settings.txt` | Buffer, thread count | Numerical settings |
| `Output.txt` | Which columns to write | Same |
| `ParamInput.txt` | File manifest | Same |

Copied from `--template-a` / `--template-b` (default `A_200W_316L` /
`B_200W_316L`). Nothing to do.

> The one time you touch this class: writing calibrated `Depth_Z` /
> `Efficiency` back into `Beam.txt` after `calibrate`. Change it once — every
> case follows.

### 2. Regenerated per layer and per hatch — **never edit by hand**

| File | Why it must be recomputed |
|---|---|
| `Path.txt` | Each layer's scan path is entirely different |
| case A `Domain.txt` | The snapshot window must sit where the spot dwells longest |
| case A `Mode.txt` | The 20 snapshot times must keep the spot inside that window |
| case B `Domain.txt` | X/Y = powered-path bounding box ± 1.5 mm |
| `coverage_target_points.txt` | Target points follow the cross-section |

All five are generated automatically. The window in the shipped
`A_200W_316L/Domain.txt` was hand-picked for layer 1600 at h = 100 µm; applying
it to another layer or hatch yields a window with no scan tracks in it.

**Window selection rule** — identical for every hatch and every layer, which is
what makes runs comparable:

> Maximise the powered dwell time with the spot inside the window and more than
> 300 µm from its edge. Break ties by maximising the 10–90 percentile span of the
> dwell-weighted time distribution.

The tie-break deliberately optimises dwell-weighted spread rather than the
first-to-last visit gap. One measured window spanned 98.5% of the build from
first to last visit but had nearly all its dwell in the final 20%, so the 20
snapshots only covered 80.6–92.1% of the scan history.

**Snapshot time selection** samples quantiles of *cumulative eligible dwell*, not
wall-clock time. The spot enters the window in bursts (once per track), so
uniform sampling in time lands mostly in gaps and then snaps to the start of the
next burst — measured, that put 13 of 20 snapshots inside the same 0.06%
interval, leaving only 7 distinct frames.

### 3. Passed on the command line

`--layer-z-mm` or `--layer-number` · `--cli <hatch>=<path>`

**Always pass hatch spacing explicitly.** With 119 contour segments in the path,
automatic inference collapses to 21.33 µm.

---

## Running other layers

```bash
python3 lpbf.py layer \
  --cli source/Test/CLI/100µm_3dbenchy.cli \
  --layer-z-mm 11.4 \
  --thesis-bin 3DThesis/build/bin/3DThesis \
  --workdir source/Test/316L/L_z11.4 \
  --output-dir source/LPBF_Agent/results/L_z11.4
```

Window, snapshot times, coverage points and both domains are recomputed for the
new layer; template parameters are reused unchanged.

### Estimate the cost first

Layer 1600 is an isolated 6.4 × 6.4 mm section → 80,498 coverage points → case B
takes about 19 minutes.

The largest sections in this part are near z ≈ 12.6–13.6 mm at
**53.7 × 28.1 mm** — roughly 30× the area → about **2.04 million points** →
case B would take **8 hours or more**.

```bash
python3 lpbf.py probe --cli <cli> --layer-z-mm 12.750
```

Above roughly 500,000 points, choose one:

1. **Restrict the evaluation region** (recommended).
   `--region-mm 6` evaluates a 6 × 6 mm box centred on the busiest area.
   Point count returns to ~116,000, about 27 minutes, **with no loss of in-plane
   resolution**. The scan path fed to 3DThesis is still the whole layer — the
   `Buffer` setting accounts for heat from tracks outside the box. Only the region
   where fusion is checked shrinks, and the conclusion is scoped to it.
2. **Coarsen the grid.** `--grid-um 50` cuts points to a quarter, but the unmelted
   inter-track ridge is only ~12 µm wide and a 50 µm grid takes 2 points per
   100 µm period. Suitable for a rough look, not for conclusions.
3. **Accept the runtime** and run overnight.

> Larger sections are generally *easier* to print. Layer 1600 is the top of a
> chimney — an isolated small section with very poor heat sinking, among the
> hardest positions in the part. The value of running another layer is showing
> the method holds beyond a small section.

---

## Subcommand → underlying script

| Subcommand | Underlying |
|---|---|
| `layers` | `CLI_Trans.py --list-layers` |
| `probe` | `CLI_Trans.py` + `generate_coverage_points.py` on a temporary case |
| `layer` / `sweep` | `LPBF_Agent/sweep_hatch.py` |
| `score` | `LPBF_Agent/run_all.py` → `lpbf_score/scorer.py` |
| `report` | `LPBF_Agent/report_summary.py` |
| `calibrate` | `Test/316L/C_calib_200W_316L/calibrate_beam.py` |
| `uq` | `LPBF_Agent/uq_envelope.py` |

`generate_coverage_points.py` is an 8-line wrapper; the implementation lives in
`lpbf_score/coverage_generator.py`.

---

## `uq` vs `calibrate`

Both sweep a `Depth_Z × Efficiency` grid, but they belong to different stages.

| | `uq_envelope.py` | `calibrate_beam.py` |
|---|---|---|
| Precondition | **No** measured data | **Has** measured pool sizes |
| What it does | Sweeps the whole grid and reports only conclusions that hold **everywhere** on it; then narrows the feasible region using published process windows | Finds the single combination with the smallest error against measurements |
| Output | An uncertainty band | Concrete values |
| Status here | ✅ usable now | ⏳ waiting on `calib_targets.csv` |

**The order is `uq` → obtain data → `calibrate`.** Do not skip `uq`: before
calibration it is the only route to a defensible conclusion.

> `uq_envelope.py` contains a useful finding: **3DThesis is exactly linear in
> `Efficiency`** (`Init.cpp:868` only does `q = q*eff*2`, and the kernel is linear
> in q; measured ratio 2.000000 when doubling 0.35 → 0.70). Sweeping `Efficiency`
> therefore does not require re-running — it can be rescaled in post-processing.
> `calibrate_beam.py` currently re-runs honestly; there is a 4× speedup available
> here if the grid ever grows too large.
