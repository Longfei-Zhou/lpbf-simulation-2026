# Slice files go here

`*.cli` files are not tracked (30–50 MB each, ~200 MB for four). Regenerate them
with Netfabb from `source/3dbenchy.stl`.

**Settings and validation values: [`docs/SLICING.md`](../../../docs/SLICING.md).**
It lists every setting that must stay identical and the measured values the
validator checks against — layer thickness 0.030 mm, 1601 layers, ASCII,
units 0.005 mm, 113°/layer rotation, top-layer hatch angle 80.00°,
2 contour entities with 119 segments.

Expected filenames:

```
70µm_3dbenchy.cli      hatch 0.070 mm
80µm_3dbenchy.cli      hatch 0.080 mm
100µm_3dbenchy.cli     hatch 0.100 mm
120µm_3dbenchy.cli     hatch 0.120 mm
```

Validate before simulating:

```bash
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 100=source/Test/CLI/100µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 --dry-run
```

It must print `controlled-comparison precondition holds`.
