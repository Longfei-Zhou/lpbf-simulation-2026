# Calibrating `Depth_Z` and `Efficiency`

These are the only two free parameters in the model, and together they determine
melt depth — which owns 40% of the LOF budget via the `(t/D)²` term.

```
Depth_Z     Effective penetration depth of the volumetric heat source.
            Real absorption depth in steel is tens of nanometres; energy is
            deposited at the surface. What carries it 40 µm down is powder
            multiple reflection, surface depression and the keyhole.
            3DThesis resolves none of them — all of it is lumped into this one
            number.

Efficiency  Effective absorptivity. Trapp et al. (2017) measured it varying by a
            factor of two with power and velocity. Not a material constant.
```

## Why you cannot copy them from a paper

**They are implementation-bound.** Every code defines its heat-source shape
differently (double ellipsoid, cylindrical, exponential decay). A penetration
depth fitted for one kernel means something different inside 3DThesis's
`exp(-3·dx²/φx)`. This project already hit that trap once: the original
`Width` = 4.5e-5 came from not applying the `√6·σ` convention.

**They are only constrained as a pair.** Melt-pool size constrains the
*combination*. Paper A with absorptivity 0.3 and depth X, and paper B with 0.5
and depth Y, can both reproduce the same pool. Copying one without the other
guarantees a mismatch.

**They are not constants.** `Depth_Z` absorbs keyhole depth, which depends on
power, velocity and spot size. "What is 316L's `Depth_Z`" has no answer.

## What you *can* take from a paper

**Measured melt-pool width and depth.** Those are the *targets*, not the
parameters. Take someone's cross-section measurements, then fit your own
`Depth_Z` and `Efficiency` so the model reproduces them. `calibrate_beam.py`
does not care whether the numbers came from your metallography or a publication.

---

## Two literature anchors (verified 2026-07-31)

Both come from the same paper, which validates a **heat conduction model** — the
same model class as 3DThesis — on 316L bare plate with 200+ single tracks:

> *Validity of Thermal Simulation Models for Different Laser Beam Shapes in
> Bead-on-Plate Melting*, **Integr. Mater. Manuf. Innov.** (2024)
> [doi:10.1007/s40192-024-00382-2](https://doi.org/10.1007/s40192-024-00382-2)

### Anchor 1 — absorptivity ≈ 0.40 (grade B)

> "The laser absorptivity **is calibrated to 0.4** based on the result
> corresponding to a Gaussian-shaped laser beam, and the same value is used for
> all other thermal simulations in this paper."

This project's current 0.35 is the 3DThesis README's generic default with no
provenance. 0.40 is at least same-material, same-model-class, experimentally
calibrated — use it as the **centre of the search grid**, not as a final value.
Their model applies a surface flux, not a volumetric source, so it says nothing
about `Depth_Z`.

### Anchor 2 — validity indicator: peak temperature ≈ 6000 °C (grade A)

> "up to a simulated maximum temperature of approximately **6000 °C**, the
> difference in melt pool depth remains **below 20%**. Beyond this temperature,
> the model **under-predicts depth**, with the error increasing linearly."

This one is usable immediately, with no simulation required.

---

## Applying anchor 2 to this project

### Where this layer sits

Peak temperature across the 20 snapshots of layer 1600:

| | K | °C |
|---|---|---|
| min | 6125 | 5852 |
| median | 6248 | 5975 |
| max | 6442 | 6169 |

**Nine of twenty exceed the 6000 °C threshold; the median sits right on it.**
The case is on the validity boundary — neither safely inside nor clearly outside.

### The error direction is known

Past the threshold the model **under-predicts** depth. The real pool is deeper
than computed. That is consistent with what we observe: the model demands hatch
be shrunk to 70 µm while the literature achieves 99.8% density at 120 µm.

### Does the FAIL verdict survive its own error bar?

With P10 width 116.29 µm, P10 depth 46.22 µm, h = 100 µm, t = 30 µm:

```
LOF = (100/116.29)² + (30/46.22)² = 0.7395 + 0.4214 = 1.161   →  FAIL
```

The width term 0.7395 is fixed, so LOF ≤ 1 requires depth ≥ **58.8 µm** —
**27.2% larger** than computed.

| Depth under-prediction | Real depth | LOF | Verdict |
|---|---|---|---|
| 0% (as computed) | 46.2 µm | 1.161 | FAIL |
| 10% | 50.8 µm | 1.088 | FAIL |
| 20% (published in-range bound) | 55.5 µm | 1.032 | FAIL, barely |
| **27%** | **58.8 µm** | **1.000** | **critical** |
| 30% | 60.1 µm | 0.989 | PASS |

**Flipping the verdict needs 27%. The published in-range error bound is 20%, and
half of this layer's snapshots are already past the threshold where the error
grows.** The FAIL verdict lies at the edge of the model's own documented error
band and cannot be reported as a physical conclusion.

### A separate observation: at h = 120 µm, width is the binding constraint

```
h = 120 µm  →  (120/116.29)² = 1.065
```

The width term alone already exceeds 1, independent of depth. So making the
model pass at the literature's 120 µm cannot be done by increasing `Depth_Z` —
melt **width** has to exceed 120 µm, and width is driven mainly by `Efficiency`
and the spot.

This is why both anchors matter together, and why `calibrate_beam.py` fits width
and depth **simultaneously**. Fitting depth alone yields a parameter set that
still fails at 120 µm.

---

## Running the calibration

`calib_targets.csv` lists five conditions from *Micromachines* 15 (2024) 170
(316L, 100 µm spot — matching this project — bare plate and powder plate,
cross-sections measured). The `width_um` and `depth_um` columns are **empty**:
that paper reports melt-pool dimensions only in bar charts (Figures 5–7), with
no numbers in the text or tables. They have to be read off the figures by hand.

When reading them:

1. **Use bare plate (BP) only.** In the simulation z = 0 is the solid top surface
   — there is no powder layer. Powder-plate depths use a different datum and mixing
   them introduces a systematic bias.
2. **Consider excluding N01** (260 W / 0.52 m/s). The paper reports its
   width-to-depth ratio as 0.633, i.e. already in keyhole mode, where a
   conduction model does not apply. Including it drags the fit.
3. **N02** (300 J/m) is closest to this project's 250 J/m — weight it highest.

Then:

```bash
python3 lpbf.py calibrate --thesis-bin 3DThesis/build/bin/3DThesis
```

The script sweeps a `Depth_Z × Efficiency` grid, runs each condition, and reports
the combination with the smallest weighted error. Depth is weighted 2× by default
because it dominates the LOF criterion and is the least certain quantity.

### The peak-temperature check cannot run from this case

`read_pool()` deliberately refuses to report a peak temperature from the
solidification output. The `T` column there maxes out at exactly 1708 K — the
liquidus — because it records *the temperature at which each point solidified*,
not the peak of its thermal history. Using it for the 6000 °C check would report
"within validity range" unconditionally, which is worse than not checking.

To verify anchor 2 for a calibrated parameter set, run the same parameters once
more in Snapshots mode and take `max_temperature`.

---

## After calibration

Write the fitted values back into `Beam.txt` — once, in the template. Every case
picks them up automatically.

Then update the grading in
[`PARAMETER_PROVENANCE.md`](PARAMETER_PROVENANCE.md):

| Parameter | Before | After (using published pool sizes) |
|---|---|---|
| `Efficiency` | D — uncalibrated | C — same material and model class, but not this machine |
| `Depth_Z` | D — uncalibrated | C — same caveat |
| Model validity boundary | not recorded | **A** — 6000 °C, quotable |

Reaching grade **B** requires metallography from your own machine and powder.
A borrowed calibration is a real improvement over a literature default, but it
must be labelled as borrowed.
