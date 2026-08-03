# 3DThesis to LPBF 316L Single-Layer Scoring System

## Scope

This system is a physics-informed process screener, not nondestructive inspection.
It reports:

1. a release state: `PASS`, `REVIEW`, or `FAIL`;
2. a 0–100 continuous quality score;
3. data completeness and evidence adequacy as separate quantities;
4. calibration and model-scope warnings.

A high score is not a calibrated probability of build success. Even `PASS` does
not replace CT, density, metallography, tensile, or fatigue validation.

One 3DThesis case represents one scanned layer. Multiple CSV `z` planes are
depths below that layer, not separate physical build layers.

## Why relative normalization is not used

Older approaches normalized each dataset to itself and treated high temperature
or cooling rate as a direct defect score. That is unsuitable because:

- a one-sample standard deviation can collapse to zero and produce an automatic 100;
- adding another case changes the first case's score;
- self-normalized values have no fixed physical meaning;
- a conduction-only model does not solve recoil pressure, free-surface flow,
  gas entrainment, spatter, or fracture.

This implementation uses fixed physical criteria, geometric overlap, robust
statistics, and explicit missing-evidence handling.

## Inputs

The scorer accepts a final Solidification case and optional Snapshot and thermal
history cases.

| Category | Fields or source |
|---|---|
| Coordinates | `x`, `y`, `z` |
| Temperature | `T`, optional temperature history |
| Solidification | `tSol`, `G`, `Gx`, `Gy`, `Gz`, `V`, `dTdt`, `eqFrac`, `numMelt` |
| Pool statistics | `MP_width`, `MP_length`, optional `MP_depth` or `depth` |
| Secondary solidification | `H`, `Hx`, `Hy`, `Hz`; diagnostic only |
| Process definition | Material, Domain, Beam, Path, Settings, Mode, and Output files |
| Physical slice | source ASCII CLI for layer thickness and solid/void geometry |
| RDF | `tm`, `tl`, `cr`; downstream and cross-check diagnostics only |

Additional numeric fields are summarized in `feature_summary.json`. Final-file
discovery is restricted to `.Solidification.Final` so intermediate schemas are
not silently concatenated into the assessment.

## Provenance and stress tests

`--stress-test` accepts deliberate coordinate shifts, parameter differences, and
concatenated trailing Path tokens as test data. It does not disable physical
failure criteria. Solidification and Snapshot cases with inconsistent Material,
Beam, or Path inputs remain isolated and cannot be combined into a high score.

Path parsing mirrors 3DThesis: only the first six fields on each physical line
are executed. Trailing concatenated tokens are reported as a diagnostic.

## Component scores

Weights are configured in `config/scoring.yaml`.

### 1. Top-surface coverage: 20%

The executed powered Path is matched to one CLI layer. When matching succeeds,
the full even-odd CLI solid cross-section becomes the required target. Holes,
cavities, and exterior voids are excluded; solid regions missing scan paths
remain in the denominator.

At the top plane and previous-layer interface:

$$
C_{top}=\frac{N_{top}(numMelt\ge1)}{N_{top,target}},\qquad
C_{interface}=\frac{N_{interface}(numMelt\ge1)}{N_{interface,target}}
$$

$$
S_{coverage}=100C_{top}
$$

Coverage at or above 0.99 reaches the release target, 0.95–0.99 requires review,
and below 0.95 is a hard failure.

Only a complete target grid containing both `numMelt == 0` and melted points is
authoritative. Surface tracking omits points for multiple possible reasons and
cannot prove complete coverage. `Tracking None` or a Custom top/interface target
domain provides authoritative evidence.

Snapshot temperature fields can estimate sampled coverage with

$$
T_{max,sampled}(x)=\max_j T(x,t_j)
$$

but missing a melt instant does not prove that a location never melted. Sparse
snapshot coverage is therefore a lower bound, not a release gate.

### 2. Interface fusion: 30%

$$
S_{fusion}=100C_{interface}
$$

The plane at `z=-t` is the thermal proxy for fusion into the previous layer. It
uses the same 0.99 release target and 0.95 hard-failure boundary as top coverage.

### 3. Lack-of-fusion geometry: 10%

The Tang, Pistorius, and Beuth elliptical overlap criterion is

$$
I_{LOF}=\left(\frac{h}{W}\right)^2+\left(\frac{t}{D}\right)^2
$$

where `h` is hatch spacing, `t` is physical layer thickness, `W` is melt-pool
width, and `D` is melt depth. `I_LOF <= 1` indicates geometric inter-track and
inter-layer overlap.

Grid uncertainty is reported at three levels:

$$
I_{nominal}: (W,D)=(P_{10}(W),P_{10}(D))
$$

$$
I_{mid}: (W,D)=(P_{10}(W)-0.5\Delta x,P_{10}(D)-0.5\Delta z)
$$

$$
I_{conservative}: (W,D)=(P_{10}(W)-\Delta x,P_{10}(D)-\Delta z)
$$

The continuous score uses the half-cell value:

$$
S_{LOF}=\frac{100}{1+\exp[-8(1-I_{mid})]}
$$

The physical boundary maps to 50. A nominal failure is a hard failure; nominal
pass plus conservative failure becomes `RESOLUTION_SENSITIVE` and `REVIEW`.

### 4. Keyhole geometry proxy: 10%

$$
A_{KH}=\frac{P_{90}(D)}{P_{10}(W)}
$$

$$
S_{KH}=\frac{100}{1+\exp[-15(0.5-A_{KH})]}
$$

An aspect ratio above 0.5 triggers review and above 0.8 is severe. This is only
a geometry proxy because 3DThesis has no vapor cavity or recoil-pressure model.

If peak temperature exceeds the configured conduction-model range, the aspect
ratio remains visible but the safety score is withheld as
`NOT_SCORED_MODEL_OUT_OF_SCOPE`. This prevents false keyhole precision from an
out-of-scope temperature solution.

### 5. Melt-pool consistency: 10%

For width and depth coefficients of variation:

$$
CV_x=\frac{\sigma_x}{\bar{x}},\qquad
S_{pool}=100\exp\left[-\left(\frac{\operatorname{mean}(CV_W,CV_D)}{0.35}\right)^2\right]
$$

The scale 0.35 is an engineering default that should be recalibrated using
measured cross-sections.

### 6. Thermal-field consistency: 12%

3DThesis provides approximately `dT/dt = G V`. `G/V` is related to solidification
morphology while `G·V` is related to microstructural scale. The scorer does not
assume that either direction is universally better. It evaluates dispersion in
`log10(G)`, `log10(V)`, and `log10(|dTdt|)` using

$$
MAD_x=\operatorname{median}|x-\operatorname{median}(x)|
$$

$$
S_{thermal}=100\exp\left[-\left(\frac{\operatorname{mean}(MAD)}{0.35}\right)^2\right]
$$

Absolute G/V, G·V, and `eqFrac` values remain available for later material- and
machine-specific calibration.

### 7. Remelt control: 8%

One overlap remelt is normally beneficial, so `numMelt <= 2` is not penalized.
Repeated remelting is treated as a heat-accumulation proxy:

$$
E=\operatorname{mean}[\max(numMelt-2,0)]
$$

$$
S_{remelt}=100\exp(-E/1.5)
$$

### Process and microstructure diagnostics

G/V, G×V, |dT/dt|, V, `tSol`, and remelt distributions are written under
`assessment.json/process_diagnostics`. `tSol` is an absolute solidification event
time, not liquid residence duration. These quantities are not included in the
release score until same-material, same-machine experimental windows exist.

## Total score and missing components

$$
Q=\frac{\sum_{i\in available}w_iS_i}{\sum_{i\in available}w_i}
$$

Unavailable components are omitted rather than fabricated as zeros. Remaining
weights are renormalized, while completeness and adequacy scores are reduced.
Consequently, a high quality score with missing authoritative coverage still
produces `REVIEW`, never a false `PASS`.

## Verdict rules

- `FAIL`: authoritative top or interface coverage below 0.95, nominal LOF above
  1, or another enabled hard physical failure.
- `PASS`: both authoritative coverages at least 0.99, no hard/review flags,
  quality at least 75, completeness at least 85, and adequacy at least 70.
- `REVIEW`: coverage between 0.95 and 0.99, a grid-sensitive LOF interval,
  missing or weak evidence, aspect ratio beyond the review threshold, or any
  release precondition not satisfied.

Continuous grades are A >= 90, B >= 80, C >= 70, D >= 60, and F < 60. They do
not override hard coverage or LOF gates.

## Energy-density diagnostics

The report includes

$$
E_l=\frac{\eta P}{v},\qquad E_v=\frac{\eta P}{vht}
$$

only as diagnostics. Energy density collapses process combinations with
different pool shapes into one number and cannot distinguish LOF from keyhole
behavior.

## Model boundaries

3DThesis uses a moving three-dimensional Gaussian source in a semi-analytical
conduction solution. It can screen fusion geometry, temperature and
solidification consistency, remelting, and G/V or G·V trends. It can only flag a
keyhole proxy and cannot independently establish spatter, balling, gas porosity,
real cracking, residual stress, distortion, or fatigue life.

## Running and outputs

```bash
python3 run_all.py \
  --solidification-dir /path/to/solidification \
  --snapshots-dir /path/to/snapshots \
  --source-cli /path/to/source.cli \
  --layer-id Layer-1
```

Use `--layer-thickness-um` when no source CLI is available and
`--hatch-spacing-um` when the path does not contain a reliably inferable parallel
hatch structure.

| Output | Purpose |
|---|---|
| `Dashboard.html` | Most direct browser view |
| `00_READ_ME_FIRST.md` | Concise conclusion and next actions |
| `Layer_Assessment.csv` | Verdict and total score |
| `Action_Plan.csv` | Prioritized problems and actions |
| `Assessment_Report.md` | Formulas, evidence, and full technical detail |
| `Problem_Diagnosis.csv` | Location, evidence, inferred cause, action, confidence |
| `Coverage_Diagnostics.csv` | Interior, endpoint, boundary, and centerline groups |
| `Score_Breakdown.csv` | Component weights and contributions |
| `assessment.json` | Complete machine-readable result |
| `snapshot_features.csv` | Per-snapshot pool geometry |
| `feature_summary.json` | Numeric-field statistics |
| `coverage_target_points.txt` | Custom points for an authoritative coverage run |

## References

1. [3DThesis](https://github.com/ORNL-MDF/3DThesis).
2. Stump and Plotkowski, semi-analytical conduction modeling,
   [doi:10.1016/j.commatsci.2020.109861](https://doi.org/10.1016/j.commatsci.2020.109861).
3. Tang, Pistorius, and Beuth, LOF geometry,
   [doi:10.1016/j.addma.2016.12.001](https://doi.org/10.1016/j.addma.2016.12.001).
4. Experimental evaluation of LOF and energy-density limits,
   [doi:10.1007/s00170-023-11163-0](https://doi.org/10.1007/s00170-023-11163-0).
5. Cunningham et al., keyhole transition,
   [doi:10.1126/science.aav4687](https://doi.org/10.1126/science.aav4687).
6. LPBF 316L G/R and cellular solidification,
   [doi:10.3390/met8080643](https://doi.org/10.3390/met8080643).

Production probability calibration requires labeled CT/metallography porosity,
measured pool geometry, roughness, crack/spatter observations, and mechanical
properties where relevant. Preserve these physics features and fit
`P(success | features)` with held-out process conditions or held-out parts.
