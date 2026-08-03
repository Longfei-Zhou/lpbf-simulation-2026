# C_calib_200W_316L Single-Track Calibration Case

## Purpose

3DThesis is a semi-analytical conduction model without explicit Marangoni flow,
evaporation, recoil pressure, or a keyhole cavity. Two effective parameters must
therefore be calibrated:

- `Beam.txt / Efficiency`, which primarily controls absolute pool size and peak temperature;
- `Beam.txt / Depth_Z`, which primarily controls depth and aspect ratio.

Use measured single-track cross-section metallography to fit both parameters,
then copy the calibrated values into the A snapshot and B solidification case
templates. Absolute melt-pool dimensions should not be treated as validated
before this step.

## Why the B case cannot be used for calibration

The B case uses a Custom point domain. In 3DThesis:

- `depth` requires `Tracking Surface`;
- melt-pool statistics require structured-grid resolution and `ijk_to_p` indexing;
- Custom points bypass structured-domain parameter construction.

Consequently, `MP_width`, `MP_length`, and `MP_depth` are not valid calibration
outputs in the B case. This C case uses a structured grid with `Tracking Surface`.

## Case configuration

| Item | Value |
|---|---|
| Scan | Three parallel tracks at y = -0.1, 0, +0.1 mm; length 1.4 mm; speed 0.8 m/s |
| Inter-track move | 5.0e-4 s including jump and settling time |
| Grid | 5 µm X × 5 µm Y × 2.5 µm Z; 321 × 101 × 49 points |
| Time step | 2.5e-6 s; beam travel 2 µm per step |
| Tracking | `Tracking Surface` |

The expected width spans roughly 13 cells and the expected 30–50 µm depth spans
roughly 12–20 Z cells.

## Calibration procedure

1. Print three to five single-track conditions at different powers or speeds.
2. Measure width W and depth D from transverse metallographic sections.
3. Enter the measurements in `calib_targets.csv`.
4. Run the grid search:

   ```bash
   python3 calibrate_beam.py \
     --thesis-bin /absolute/path/to/3DThesis \
     --resume
   ```

5. Fit depth and width together. Depth-only calibration can still fail the
   inter-track overlap criterion.
6. Copy the selected `Depth_Z` and `Efficiency` to both A and B templates.
7. Cross-validate against a power/speed condition that was not used for fitting.

Read a finished case with:

```bash
python3 read_meltpool.py Data/calib.Solidification.Final.csv
```

Existing `Data/` files, if present, are verification output generated from this
source tree. They are not measurements from the target machine and may be
regenerated.
