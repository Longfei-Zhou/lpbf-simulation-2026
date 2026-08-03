# Calibration Evidence and Literature Constraints

Updated 2026-07-31 for the two uncalibrated parameters, `Depth_Z` and
`Efficiency`.

## What the literature provides

The target melt-pool widths and depths in Micromachines 2024 are shown only in
bar charts (Figures 5–7), not in machine-readable tables. They must be digitized
manually before `calibrate_beam.py` can perform an experimental fit.

A more directly relevant thermal-model study is:

Kollmannsberger et al., *Validity of Thermal Simulation Models for Different
Laser Beam Shapes in Bead-on-Plate Melting*, *Integr. Mater. Manuf. Innov.*
(2024), [doi:10.1007/s40192-024-00382-2](https://doi.org/10.1007/s40192-024-00382-2).

It reports two useful constraints for 316L conduction-model simulations:

1. Laser absorptivity was calibrated to 0.40 for a Gaussian beam. This supports
   0.40 as the center of an Efficiency search grid, not as a machine-specific
   final value, because the cited model uses a surface source while 3DThesis
   uses a volumetric source coupled to `Depth_Z`.
2. Up to a simulated peak temperature of about 6000 °C, melt-depth error remains
   below approximately 20%. Above that point the conduction model increasingly
   under-predicts depth.

## Relevance to the reference layer

The 20 reference snapshots span approximately 6125–6442 K, or 5852–6169 °C,
with a median near 5975 °C. Nine snapshots exceed the cited 6000 °C boundary.
The case therefore sits on the edge of conduction-model validity rather than
comfortably inside it.

For the current P10 width of 116.29 µm, P10 depth of 46.22 µm, 100 µm hatch,
and 30 µm layer:

```text
LOF = (100/116.29)^2 + (30/46.22)^2 = 1.161
```

Passing requires a depth of at least 58.8 µm, 27.2% above the simulated value.
The published error is already close to this magnitude at the validity boundary
and grows beyond it. The simulated FAIL is therefore near the model's documented
uncertainty range and must not be presented as an experimentally established
failure.

At a 120 µm hatch, the width term alone is `(120/116.29)^2 = 1.065`. Increasing
only `Depth_Z` cannot make that point pass; Efficiency and beam width control
must also reproduce the measured pool width. Width and depth must be fitted
together.

## Measurement rules

- Use bare-plate data because the model places z=0 at a solid surface and has no
  powder layer.
- Exclude N01 (260 W, 0.52 m/s) from conduction-model fitting because its reported
  width/depth ratio of 0.633 indicates keyhole behavior.
- Give N02 (300 J/m) additional attention because it is closest to the study's
  250 J/m condition.
- Treat any good fit above the 6000 °C model boundary as extrapolation.

After digitizing the bar-chart values, run:

```bash
python3 calibrate_beam.py --thesis-bin /absolute/path/to/3DThesis --resume
```

## Evidence grades

| Parameter | Current status | Basis |
|---|---|---|
| `Efficiency` | C-level prior centered at 0.40 | Experimentally calibrated for the same material and model class, but not this machine/source formulation |
| `Depth_Z` | D-level uncalibrated | No independent external anchor; requires melt-pool fitting |
| 6000 °C validity boundary | A-level quoted constraint | Directly reported by the cited thermal-model study |

## Sources still requiring manual access

- *Micromachines* 2024, 15(2), 170: width and depth appear only in Figures 5–7;
  [PMC10890519](https://pmc.ncbi.nlm.nih.gov/articles/PMC10890519/).
- *Materials & Design* 2026, "Melt pool geometry and process windows in PBF of
  316L: Comprehensive single-source dataset": check the full text through an
  institutional account for tabulated measurements.
