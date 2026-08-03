# Slicing in Netfabb

Slice files are not in this repository (four of them total 200 MB). Regenerate
them from `source/3dbenchy.stl`.

**Change hatch distance and nothing else.** Everything below is measured from a
known-good export and is what the validator checks against.

---

## What to export

| Hatch distance | Filename | Why this point |
|---|---|---|
| 0.070 mm | `70µm_3dbenchy.cli` | The model's own recommendation |
| 0.080 mm | `80µm_3dbenchy.cli` | Brackets the model's critical value of 81 µm |
| 0.100 mm | `100µm_3dbenchy.cli` | Baseline |
| 0.120 mm | `120µm_3dbenchy.cli` | Literature point — *Metals* 11 (2021) 832 measured 99.8% density here under otherwise identical parameters |

70 and 120 bracket the whole disputed range; four points are enough to
interpolate where LOF crosses 1.

## Export settings

`Slices > Export > as CLI`:

- ✅ Export Contours
- ✅ Export Hatches
- ❌ **Convert all Slices to Hatches**

The third one matters. Tutorials often enable it, but it turns contours into
hatch segments. The scorer separates contours from infill by scan speed and
direction (`infill_segments()`); once merged they cannot be told apart, and
`infer_hatch_spacing()` collapses to 21 µm. A correct export keeps `$$HATCHES`
and `$$POLYLINE` separate.

---

## What must stay identical

Every value below is measured from a known-good export. After exporting, check
them — the validator does this automatically, but knowing the numbers helps you
diagnose a mismatch.

| Item | Required value | Where to check |
|---|---|---|
| Layer thickness | **0.030 mm** | header `$$LAYERS` + first layer z |
| Layer count | **1601** | `$$LAYERS/001601` |
| Format | **ASCII** (not binary) | header line 2, `$$ASCII` |
| Units | **0.005 mm** | `$$UNITS/00000000.005000` |
| First layer z | **0.030 mm** | first `$$LAYER` record |
| Bounding box | must match across all four | `$$DIMENSION`, character for character |
| Layer rotation | **113° per layer** | angle increment between consecutive layers |
| Top-layer hatch angle | **80.00°** | determined by start angle + 1600 × 113° |
| Contours on top layer | **2 entities, 119 segments** | `$$POLYLINE` |
| Beam compensation | unchanged | changing it shifts the contours |

### Three traps

**1. Do not touch Orient part.** A previous export had the `$$DIMENSION` Z span
come out as 49.26 instead of 48.000 because the part had been rotated by ~3.31°.
Z span must be exactly 48.000.

**2. Keep start angle and layer rotation unchanged.** The top layer's hatch angle
is `start angle + 1600 × 113°`. If it changes, the scan direction rotates
relative to the cross-section, which changes track lengths, jump count and
inter-track heat accumulation — no longer a controlled comparison. This has
already happened once: one export came out at 147° while the others were 80°,
a difference of exactly one layer rotation (180 − 113 = 67°).

**3. Watch the units on hatch distance.** One export had 1.0 mm instead of
0.1 mm — 10 hatch segments instead of 94, and a 9.7 MB file instead of ~38 MB.

### Absolute z may differ — that is harmless

If "part start height above platform" differs between exports, the whole part
shifts in z (observed: 51.000 vs 48.000, exactly the 3.000 mm start height).
This does **not** affect physics — `CLI_Trans.py` moves the selected layer to
z = 0 and absolute height never enters the simulation. Select layers by index
instead:

```bash
python3 lpbf.py sweep --cli ... --layer-number -1    # -1 = topmost layer
```

The validator reports the z difference as a note, not an error, and still checks
that the cross-sections match.

---

## Scan speed and power do not matter here

The CLI carries no speed or power information. Speeds are added afterwards by
`CLI_Trans.py` (`--scan-speed`, `--contour-speed`) and power lives in `Beam.txt`.
Whatever Netfabb has configured for them is irrelevant to this pipeline.

---

## Validate before simulating

```bash
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 100=source/Test/CLI/100µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 --dry-run
```

Runs in seconds and checks, before any simulation starts:

- each file really has a layer at the requested position
- measured hatch spacing matches the declared value (2 µm tolerance)
- contour segment count, cross-section size, layer count and hatch angle agree
  across all four files

It prints `controlled-comparison precondition holds` when everything passes.

**Do not bypass it with `--skip-validation`.** That discards the premise that
this is a controlled comparison, and the hours of simulation that follow produce
data that cannot be compared.
