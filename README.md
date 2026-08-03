# LPBF Single-Layer Printability Screening

A scan-path-resolved screening pipeline for laser powder bed fusion. It takes the
**actual scan path of one layer** from a sliced part, runs it through a
semi-analytical heat-conduction model ([3DThesis](https://github.com/ORNL-MDF/3DThesis)),
and scores whether that layer can be built successfully.

The difference from conventional printability maps is the input. Maps sweep
**single tracks** through (power, velocity) space. This pipeline runs a **whole
layer** with its real hatch ordering, real jump timing, and real contour/infill
speed difference — so it captures inter-track heat accumulation, jump cooling,
and contour–infill interaction, none of which exist in a single-track map.

> ### ⚠️ The model is not calibrated
>
> Two parameters — `Depth_Z` (volumetric heat-source penetration depth) and
> `Efficiency` (absorptivity) — are still literature defaults. **There is no
> experimental data point anywhere in this chain.**
>
> **Trustworthy:** relative results — trends, rankings, coverage fractions.
> For example: reducing layer thickness 50 → 30 µm raises interface fusion from
> 43% to 89%.
>
> **Not trustworthy:** absolute verdicts — the PASS/FAIL decision itself,
> "hatch should be 70 µm", any absolute temperature.
>
> Details in [`docs/PARAMETER_PROVENANCE.md`](docs/PARAMETER_PROVENANCE.md) and
> [`docs/CALIBRATION.md`](docs/CALIBRATION.md).

---

## Requirements

| | Version | Required |
|---|---|---|
| Python | 3.10+ | yes (uses `X \| None` syntax) |
| numpy, pandas, pyyaml | recent | yes |
| matplotlib | recent | no — plots only |
| C++ compiler with OpenMP | — | yes, to build 3DThesis |
| CMake | 3.9+ | no — a plain makefile also ships |

**Do not enable MPI.** `Run.cpp:690` disables melt-pool statistics under MPI
domain decomposition, and the calibration case needs them. Single-node OpenMP
is sufficient.

## Install and build

```bash
git clone <this-repo> lpbf-project
cd lpbf-project

python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install numpy pandas pyyaml matplotlib

JOBS=8 bash 3DThesis/build_example.sh

python3 lpbf.py test
```

Replace `8` with the number of physical CPU cores you want to use. The expected
executable is `3DThesis/build/bin/3DThesis`. Re-run `source .venv/bin/activate`
in each new shell. On macOS, install the compiler dependencies first with
`brew install cmake libomp`; the helper supplies Apple Clang's required OpenMP
flags automatically. On Debian/Ubuntu, use
`sudo apt install build-essential cmake libomp-dev`.

Manual CMake equivalent on Linux:

```bash
cmake -S 3DThesis -B 3DThesis/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE
cmake --build 3DThesis/build --parallel 8
cmake --install 3DThesis/build --prefix 3DThesis/install
```

## Input data

Slice files (`*.cli`) are **not** in this repository — four of them total 200 MB.
Regenerate them with Netfabb from `source/3dbenchy.stl`. Exact settings and the
measured validation values are in [`docs/SLICING.md`](docs/SLICING.md).

Expected in `source/Test/CLI/`:

```
70µm_3dbenchy.cli     hatch 0.070 mm
80µm_3dbenchy.cli     hatch 0.080 mm
100µm_3dbenchy.cli    hatch 0.100 mm
120µm_3dbenchy.cli    hatch 0.120 mm
```

## Run

```bash
# 1. List candidate layers
python3 lpbf.py layers --cli source/Test/CLI/100µm_3dbenchy.cli

# 2. Estimate cost before committing — seconds, do not skip for large sections
python3 lpbf.py probe --cli source/Test/CLI/100µm_3dbenchy.cli --layer-z-mm 51.0

# 3. Validate inputs without simulating — seconds
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 100=source/Test/CLI/100µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 --dry-run

# 4. Run the full chain
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 100=source/Test/CLI/100µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 \
  --thesis-bin 3DThesis/build/bin/3DThesis \
  --max-threads 8 --resume
```

Step 3 must print `controlled-comparison precondition holds` before you spend
hours on step 4. It has already caught two real export errors: a hatch angle off
by one layer rotation, and a hatch distance entered as 1.0 mm instead of 0.1 mm.

Output lands in `source/LPBF_Agent/results/sweep_hatch/` —
`Hatch_Sweep_Report.md`, `hatch_sweep_lof.png`, `hatch_sweep.csv`.

Every command in order: [`docs/RUNBOOK.md`](docs/RUNBOOK.md).
Cluster deployment: [`docs/SERVER.md`](docs/SERVER.md).

---

## Pipeline

```
3dbenchy.stl
   │  Netfabb slicing   (the only manual step — see docs/SLICING.md)
   ▼
<h>µm_3dbenchy.cli
   │  source/CLI_Trans.py    select layer, add timing, add contours
   ▼
Path.txt
   ├─► case A  snapshots      ──3DThesis──► 20 temperature fields → melt-pool W/D
   ├─► case B  solidification ──3DThesis──► tSol/G/V/numMelt      → fusion coverage
   │      ▲ coverage_target_points.txt  (needs case B's Domain.txt first)
   └─► case C  calibration    ──3DThesis──► single-track section vs. micrographs
                                   │
                                   ▼  source/LPBF_Agent/lpbf_score
                         assessment.json + report + plot
```

**Ordering matters.** The coverage-point generator reads the sampling grid and
z-planes from case B's regular `X/Y/Z` blocks, so `Domain.txt` must exist *before*
the points are generated. `sweep_hatch.py` already does this — do not reorder it.

## Commands

| Command | Purpose | Runtime |
|---|---|---|
| `layers` | List candidate layers with cross-section size | seconds |
| `probe` | Estimate coverage-point count and wall time | seconds |
| `layer` | Full chain for one layer | ~40 min |
| `sweep` | Compare several slices (e.g. a hatch sweep) | N × 40 min |
| `score` | Re-score without re-simulating | minutes |
| `report` | Re-render the report only | seconds |
| `calibrate` | Fit `Depth_Z` / `Efficiency` to measured pool sizes | hours |
| `uq` | Uncertainty envelope for the uncalibrated parameters | minutes |
| `test` | 27 unit tests | seconds |

`layer` is `sweep` with a single `--cli`. Both share the same case-assembly code,
so "run one layer" and "sweep a set" cannot diverge through configuration drift.

`uq` and `calibrate` are two stages of one job, and the order matters:
**`uq` before you have measurements, `calibrate` after.** Until calibration data
exists, `uq` is the only route to a defensible conclusion.

### Complete command reference

Run commands from the repository root unless a command starts with `cd`. Use
`python3 lpbf.py <command> --help` for every available option. Relative paths
passed to the wrapper are resolved before it changes into a tool directory.

Inspect and select a layer:

```bash
python3 lpbf.py layers \
  --cli source/Test/CLI/100µm_3dbenchy.cli \
  --top 20

# Sort by cross-section size instead of sampling evenly through the height.
python3 lpbf.py layers \
  --cli source/Test/CLI/100µm_3dbenchy.cli \
  --top 20 --by-size
```

Estimate the point count and case-B runtime before a large run:

```bash
python3 lpbf.py probe \
  --cli source/Test/CLI/100µm_3dbenchy.cli \
  --layer-z-mm 51.0 \
  --hatch-spacing-um 100 \
  --layer-thickness-um 30 \
  --grid-um 25
```

Run one layer end to end. Add `--region-mm 6` after the executable path when a
large cross-section would generate too many coverage points:

```bash
python3 lpbf.py layer \
  --cli source/Test/CLI/100µm_3dbenchy.cli \
  --layer-z-mm 51.0 \
  --hatch-spacing-um 100 \
  --thesis-bin 3DThesis/build/bin/3DThesis \
  --max-threads 8 \
  --resume
```

Validate and run a controlled hatch sweep:

```bash
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 100=source/Test/CLI/100µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 --dry-run

python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 100=source/Test/CLI/100µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 \
  --thesis-bin 3DThesis/build/bin/3DThesis \
  --max-threads 8 --resume
```

Re-score existing 3DThesis outputs without simulation, then rebuild or print the
one-page report:

```bash
python3 lpbf.py score \
  --solidification-dir source/Test/316L/runs/h100/B_solidification \
  --snapshots-dir source/Test/316L/runs/h100/A_snapshots \
  --source-cli source/Test/CLI/100µm_3dbenchy.cli \
  --layer-thickness-um 30 \
  --hatch-spacing-um 100 \
  --layer-id benchy-L1600-h100 \
  --output-dir source/LPBF_Agent/results/layer/h100

python3 lpbf.py report source/LPBF_Agent/results/layer/h100
python3 lpbf.py report source/LPBF_Agent/results/layer/h100 --print
```

Regenerate a completed sweep summary without rerunning its cases. The CLI list
is still required so the four result directories can be identified:

```bash
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 100=source/Test/CLI/100µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --report-only
```

Run the uncalibrated-parameter envelope. Each `Depth_Z` value runs once;
Efficiency values are recovered from the model's linear temperature scaling:

```bash
python3 lpbf.py uq \
  source/Test/316L/A_200W_316L \
  --binary 3DThesis/build/bin/3DThesis \
  --output-dir source/LPBF_Agent/results/uq \
  --depth-um 20,27.5,35,45,60 \
  --efficiency 0.25,0.35,0.45,0.55,0.70 \
  --snapshots 5
```

Calibrate after entering measured `width_um` and `depth_um` values in
`source/Test/316L/C_calib_200W_316L/calib_targets.csv`:

```bash
python3 lpbf.py calibrate \
  --thesis-bin 3DThesis/build/bin/3DThesis \
  --targets source/Test/316L/C_calib_200W_316L/calib_targets.csv \
  --depth-z-um 20,25,30,35,40,50,60 \
  --efficiency 0.30,0.35,0.40,0.45 \
  --resume

# Rebuild calibration output from completed runs only.
python3 lpbf.py calibrate \
  --targets source/Test/316L/C_calib_200W_316L/calib_targets.csv \
  --report-only
```

Run the tests:

```bash
python3 lpbf.py test
```

### Direct low-level commands

The `lpbf.py` wrapper is preferred. These commands expose individual stages for
debugging and custom workflows.

Inspect a CLI, list layers, or convert one layer to `Path.txt`:

```bash
python3 source/CLI_Trans.py \
  source/Test/CLI/100µm_3dbenchy.cli /tmp/unused \
  --summary-only

python3 source/CLI_Trans.py \
  source/Test/CLI/100µm_3dbenchy.cli /tmp/unused \
  --list-layers

python3 source/CLI_Trans.py \
  source/Test/CLI/100µm_3dbenchy.cli Path.txt \
  --layer-z-mm 51.0 \
  --whole-layer --include-contours --recenter-xy \
  --scan-speed 0.8 --contour-speed 0.6 \
  --jump-speed 5.0 --jump-delay 0.0002
```

Generate authoritative top/interface coverage points after the case's regular
`X/Y/Z` Domain blocks and `Path.txt` exist:

```bash
python3 source/LPBF_Agent/generate_coverage_points.py \
  --case-dir source/Test/316L/runs/h100/B_solidification \
  --domain-file source/Test/316L/runs/h100/B_solidification/Domain.txt \
  --path-file source/Test/316L/runs/h100/Path.txt \
  --source-cli source/Test/CLI/100µm_3dbenchy.cli \
  --layer-thickness-um 30 \
  --hatch-spacing-um 100 \
  --output source/Test/316L/runs/h100/B_solidification/coverage_target_points.txt
```

Run 3DThesis from inside a case directory because `ParamInput.txt` contains
relative filenames:

```bash
cd source/Test/316L/runs/h100/A_snapshots
../../../../../../3DThesis/build/bin/3DThesis ParamInput.txt

cd ../B_solidification
../../../../../../3DThesis/build/bin/3DThesis ParamInput.txt
```

The absolute-path form is less error-prone on a server:

```bash
cd /absolute/path/to/case
/absolute/path/to/3DThesis/build/bin/3DThesis ParamInput.txt
```

Call the scorer and report generator directly:

```bash
python3 source/LPBF_Agent/run_all.py \
  --solidification-dir source/Test/316L/runs/h100/B_solidification \
  --snapshots-dir source/Test/316L/runs/h100/A_snapshots \
  --source-cli source/Test/CLI/100µm_3dbenchy.cli \
  --layer-thickness-um 30 --hatch-spacing-um 100 \
  --layer-id benchy-L1600-h100 \
  --output-dir source/LPBF_Agent/results/layer/h100

python3 source/LPBF_Agent/report_summary.py \
  source/LPBF_Agent/results/layer/h100
```

Read the calibration case's melt-pool geometry:

```bash
python3 source/Test/316L/C_calib_200W_316L/read_meltpool.py \
  source/Test/316L/C_calib_200W_316L/Data/calib.Solidification.Final.csv
```

Build and install 3DThesis with the helper, then run its bundled examples
(the example runner uses GNU `find`, so use it on Linux):

```bash
JOBS=8 bash 3DThesis/build_example.sh
cd 3DThesis/examples
bash run_examples.sh
```

---

## Scoring

Seven components; weights live in `source/LPBF_Agent/config/scoring.yaml`.

| Component | Weight | Core formula | Basis |
|---|---|---|---|
| Top-surface coverage | 20% | `N(numMelt≥1) / N_target` | geometric self-consistency |
| Interface fusion | 30% | same, on the z = −t plane | thermal proxy for interlayer bonding |
| LOF geometric margin | 10% | `I = (h/W)² + (t/D)²` | Tang, Pistorius & Beuth (2017) |
| Keyhole margin | 10% | `P90(D) / P10(W)`, boundary 0.5 | Cunningham et al. (2019) |
| Pool consistency | 10% | CV of width and depth | engineering default |
| Thermal uniformity | 12% | robust MAD of log10(G/V/dTdt) | engineering default |
| Remelt control | 8% | excess above `numMelt > 2` | one overlap remelt is not penalised |

Verdict: both coverages ≥ 99% with no hard failure → `PASS`; below 95% → `FAIL`;
otherwise `REVIEW`.

The LOF index is reported at three grid-uncertainty levels — nominal, half-cell
and conservative. When nominal passes but conservative fails, the result is
flagged `RESOLUTION_SENSITIVE` rather than silently passed.

**The seven weights are expert judgement with no data behind them.** State this in
any write-up. Change behaviour in `scoring.yaml`, never in the code.

Derivations, literature basis and model limits:
[`source/LPBF_Agent/docs/SCORING_SYSTEM.md`](source/LPBF_Agent/docs/SCORING_SYSTEM.md).

---

## Layout

```
lpbf.py                          single entry point, 9 subcommands
docs/
  PARAMETER_PROVENANCE.md        where every number comes from, graded A/B/C/D
  SLICING.md                     Netfabb settings + validation values
  RUNBOOK.md                     every command, compile → score
  SERVER.md                      cluster deployment
  PIPELINE.md                    data flow and parameter-reuse rules
  CALIBRATION.md                 literature anchors for the two free parameters
source/
  CLI_Trans.py                   CLI → Path.txt
  3dbenchy.stl                   geometry source
  Test/CLI/                      slice files — not in repo, regenerate
  Test/316L/A_200W_316L/         snapshot template        10/10/5 µm window
  Test/316L/B_200W_316L/         solidification template  custom point set
  Test/316L/C_calib_200W_316L/   calibration case         5/5/2.5 µm single track
  LPBF_Agent/
    sweep_hatch.py               full-chain driver: assembly, batching, reporting
    run_all.py                   scoring entry point
    generate_coverage_points.py  coverage-point entry point
    uq_envelope.py               uncertainty envelope + literature inversion
    lpbf_score/scorer.py         scoring core
    lpbf_score/coverage_generator.py   corridor ∩ CLI solid mask
    config/scoring.yaml          weights and thresholds
    tests/                       27 unit tests
3DThesis/                        upstream source — see 3DThesis/UPSTREAM.md
```

### Parameter files come in three classes

**Never edited.** `Material.txt`, `Beam.txt`, `Settings.txt`, `Output.txt`,
`ParamInput.txt` — copied from the templates automatically. The one time you
touch them is writing calibrated `Depth_Z` / `Efficiency` back into `Beam.txt`;
change it once and every case follows.

**Regenerated per layer and per hatch — never edit by hand.** `Path.txt`, case A's
`Domain.txt` and `Mode.txt`, case B's `Domain.txt`, `coverage_target_points.txt`.
The snapshot window in the shipped `A_200W_316L/Domain.txt` was hand-picked for
one specific path; reusing it for a different layer or hatch gives a window with
no scan tracks in it at all.

**Passed on the command line.** `--layer-z-mm` or `--layer-number`, and
`--cli <hatch>=<path>`. **Always pass hatch spacing explicitly** — with 119
contour segments in the path, automatic inference collapses to 21.33 µm.

---

## Three source-level constraints

Hard facts located in the 3DThesis source, not design choices:

1. The `depth` column requires `Tracking Surface` (`Grid.h:226`); it never appears
   with a custom point set.
2. `MP_Stats` is incompatible with a `Custom` domain (`Melt.cpp` relies on the
   structured grid's `ijk_to_p`) and silently outputs zeros.
3. `MP_depth` converts z-direction cell counts with `xres` instead of `zres`
   (`Melt.cpp:595`), so it is wrong by `xres/zres` on anisotropic grids.
   **Use the `depth` column, never `MP_depth`.**

## Known limitations

1. **`Depth_Z` and `Efficiency` are uncalibrated** — the dominant uncertainty.
   Sweeping `Depth_Z` from 20 to 60 µm swings melt depth by 41% and halves peak
   temperature. Two literature anchors are documented in
   [`docs/CALIBRATION.md`](docs/CALIBRATION.md), but single-track metallography
   is still required.
2. **The seven scoring weights are expert judgement.** A ±50% perturbation study
   is needed to show the conclusions do not depend on them.
3. **Peak temperature sits on the conduction model's validity boundary.**
   Published work on the same model class places that boundary near 6000 °C;
   9 of 20 snapshots on the reference layer exceed it. Beyond it the model
   **under-predicts** melt depth — the real pool is deeper than computed, which
   makes the pipeline conservative rather than optimistic.
4. **The keyhole aspect ratio uses `P90(D)/P10(W)`** with numerator and
   denominator drawn from *different* snapshots. This is a deliberate worst-case
   pairing; report it as an upper bound, not as any real pool's aspect ratio.
5. **Snapshots sample a local window**, so the pool-consistency CV is biased low
   relative to the whole layer.
6. **One layer, one geometry, one parameter point.** The reference layer is an
   isolated 6.4 × 6.4 mm section at the top of a chimney with very poor heat
   sinking — among the hardest positions in the part, and not representative.

## References

| | |
|---|---|
| Model | Stump & Plotkowski, *Comput. Mater. Sci.* (2020) · [3DThesis](https://github.com/ORNL-MDF/3DThesis) |
| Heat source | Nguyen et al., *Welding Journal* (1999) |
| LOF criterion | Tang, Pistorius & Beuth, *Addit. Manuf.* (2017) · [doi:10.1016/j.addma.2016.12.001](https://doi.org/10.1016/j.addma.2016.12.001) |
| Keyhole transition | Cunningham et al., *Science* (2019) · [doi:10.1126/science.aav4687](https://doi.org/10.1126/science.aav4687) |
| 316L properties | Pichler et al., *J. Mater. Sci.* 55 (2020) 4081, NIST SRM 1155a · [doi:10.1007/s10853-019-04261-6](https://doi.org/10.1007/s10853-019-04261-6) |
| Absorptivity anchor | *Integr. Mater. Manuf. Innov.* (2024) · [doi:10.1007/s40192-024-00382-2](https://doi.org/10.1007/s40192-024-00382-2) |
| Experimental comparison | *Metals* 11 (2021) 832 · [doi:10.3390/met11050832](https://doi.org/10.3390/met11050832) |

## License

3DThesis is distributed under its own license — see `3DThesis/LICENSE`.
