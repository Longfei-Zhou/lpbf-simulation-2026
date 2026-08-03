# Parameter Provenance

Where every number in this pipeline comes from, what is solid, and what is a
compromise.

Reference configuration: XM200G / 316L / 200 W / 100 µm spot / 30 µm layer /
layer 1600 (z = 51.000 mm, top of the chimney).

**Grading**

| | Meaning |
|---|---|
| **A — exact** | Vendor data, measured back from the actual files, or same-source experimental literature. Safe to put in a methods section. |
| **B — derived** | Computed from exact data by an explicit rule. The rule is defensible; some approximation remains. |
| **C — compromise** | Generic empirical value or engineering default. Findable, just not yet looked up — or forced by model structure. |
| **D — uncalibrated** | Pure placeholder. Must be fixed by experiment; the current value is not trustworthy. |

---

## 1. Slicing and geometry

| Parameter | Value | Grade | Source |
|---|---|---|---|
| Material | 316L | **A** | Xact official 316L datasheet (PraxAir/Linde 316, 2021-05) |
| Laser power | 200 W | **A** | Same, "Printer Process Specifications" |
| Spot diameter D4σ | 100 µm | **A** | Same (machine also offers 50 µm optics — confirm which is installed) |
| Layer thickness | 30 µm | **A** | Same |
| Hatch spacing | 100.0 µm | **A** | Measured back from the CLI; constant across 1600 sampled layers |
| Layer rotation | 113°/layer | **A** | Same; 12 consecutive layers all increment by 113 |
| Part geometry | 60.003 × 31.006 × 48.000 mm | **A** | `$$DIMENSION`, matches the `3dbenchy.stl` bounding box |
| Layer count | 1600 | **A** | `$$LAYERS/001600` = 48.000/0.030 ✓ |
| Selected layer | z = 51.000 mm (layer 1600) | **A** | Substrate start height 3.000 + part height 48.000 |
| Section size | 6.40 × 6.41 mm, 94 infill + 119 contour | **A** | Measured from `Path.txt` |

Nothing in this block is a compromise — all of it is reproducible from the files.

---

## 2. Path and timing (`CLI_Trans.py`)

| Parameter | Value | Grade | Source / compromise |
|---|---|---|---|
| Infill scan speed | 0.8 m/s | — | **Study variable**, not an unknown |
| Contour speed | 0.6 m/s | **C** | Convention is 0.5–0.8× infill. True value is in the Netfabb build strategy — **not looked up** |
| Jump speed | 5.0 m/s | **C** | Typical LPBF magnitude. True value is in the machine configuration — **not looked up** |
| Fixed jump delay | 2.0e-4 s | **C** | Lumps jump delay + laser on/off delay + mark/polygon delay into one value. Same source |
| Jump time model | `t = distance/speed + delay` | **B** | The old version hard-coded 1e-7 s (instantaneous). 3DThesis Mode 1 is "teleport + dwell"; neither deposits heat, so matching duration is equivalent |

**Quantified cost of the compromise:** varying all three C-grade parameters over
their plausible ranges moves single-layer time between 380 and 453 ms (±9%).
This affects inter-track cooling only — a second-order effect on melt-pool size.

---

## 3. Heat source (`Beam.txt`)

| Parameter | Value | Grade | Source / issue |
|---|---|---|---|
| `Width_X` / `Width_Y` | 6.124e-5 | **B** | = 0.6124 × D4σ(100 µm). Conversion: kernel `exp(-3·dx²/φx)` against Gaussian `exp(-dx²/2σ²)` gives `Width = √6·σ`, consistent with the 3DThesis README |
| `Power` | 200 | **A** | Vendor datasheet |
| `Depth_Z` | 3.5e-5 | **D** | Literature default |
| `Efficiency` | 0.35 | **D** | 3DThesis README's "typical LPBF value" |

### Why these two can only be calibrated

`Depth_Z` is an *effective* depth that lumps together vapour-cavity formation and
melt-pool flow. Real laser absorption depth in metal is tens of nanometres —
energy is deposited at the surface. What carries it downward is powder multiple
reflection, surface depression, and the keyhole. 3DThesis resolves none of them.

The 3DThesis authors say so directly (Coleman, Knapp, Stump, Rolchigo, Kincaid,
Plotkowski, *Additive Manufacturing* 95 (2024) 104531):

> "the dimensions and effective absorption of the volumetric heat source **are
> calibrated to reproduce melt pool dimensions observed in metallographic cross
> sections taken from single-track experiments on bare plate**"

Absorptivity is the same story. Trapp et al. (*Appl. Mater. Today* 9 (2017) 341)
measured it varying **by a factor of two** with power and velocity, driven mainly
by powder scattering and pool morphology. It is not a tabulated material constant.

### Sensitivity (measured)

| `Depth_Z` | Peak T | Width | Depth | Aspect |
|---|---|---|---|---|
| 20 µm | 6218 K | 95.2 µm | 24.8 µm | 0.26 |
| **35 µm (current)** | 4770 K | 86.0 µm | 29.3 µm | 0.34 |
| 60 µm | 3389 K | 71.6 µm | 35.0 µm | 0.49 |

**Melt depth swings 41%; peak temperature halves.** The `(t/D)²` term alone owns
40% of the LOF budget — this is the single largest uncertainty in the chain.

### Order-of-magnitude anchors (not "correct values")

| Source | Width | Depth_Z | Ratio |
|---|---|---|---|
| 3DThesis official examples | 10 µm | 10 µm | 1.00 (pure demo) |
| Repository `316L_Test/316l-1~3` | 75 µm | 45 µm | 0.60 |
| This project | 61.24 µm | 35 µm | 0.57 |

---

## 4. Material properties (`Material.txt`)

Primary source: Pichler, Simonds, Sowards, Pottlacher, *J. Mater. Sci.* 55 (2020)
4081, [doi:10.1007/s10853-019-04261-6](https://doi.org/10.1007/s10853-019-04261-6).
They measured **NIST SRM 1155a** — the same 316L batch.

| Parameter | Value | Grade | Source |
|---|---|---|---|
| `T_L` | 1708 K | **A** | Table 2: T_liquidus = 1708 ± 30 K |
| `c` | 789 J/(kg·K) | **B** | From measured specific enthalpy: 298 K → fully liquid needs 1112 kJ/kg. `c_eff = 1112e3/1409.85 = 789`. **Latent heat (measured 290 kJ/kg) is automatically included** |
| `ρ` | 7904 kg/m³ | **A** | Weighed machined cylinders, 7904 ± 25 |
| `k` | 24 W/(m·K) | **C** | They measured *resistivity*, not conductivity. Two independent estimates: linear Kirchhoff fit gives 23.8; Wiedemann–Franz gives ≈20 (low — ignores phonons). Took 24, ±15% |
| `T_0` | 298.15 K | **C** | Assumes no preheat. True value is the interlayer temperature — **not looked up** |
| CET block | removed | **B** | The original n=3 and a=1.25e6 were not from the same Gäumann fit; `eqFrac` was meaningless |

**Why ρ is the room-temperature value, not a high-temperature mean.** The specific
enthalpy 1112 kJ/kg is per unit mass, and ρ and c appear in the model only as a
product (β = q/(ρc), α = k/(ρc)). Using the reference-state density is what makes
energy balance:

```
7904 × 789 × 1409.85 = 8.79e9 J/m³ = 1112 kJ/kg × 7904 kg/m³   ✓
```

> At one point this was changed to 7600 (attempting a "mean of room temperature
> and near-liquidus"). That was wrong — it breaks the relation above. The
> original 7900 was closer.

**How much T_0 costs** (measured): going from 25 to 200 °C gives +6% width,
+6% depth, but **+15% pool tail length**. Small effect on fusion judgement;
not negligible for G/V and solidification-structure analysis.

**Recorded but not used:** solidus T_S = 1675 ± 15 K; liquid cp = 847 and
end-of-solid cp = 714 J/(kg·K); liquid density at liquidus ≈ 6936 kg/m³.

---

## 5. Domain and numerics

| Item | A (snapshots) | B (solidification) | C (calibration) | Grade |
|---|---|---|---|---|
| Domain type | local window 1.9×1.8×0.08 mm | custom points + regular XY | structured 1.6×0.5×0.12 mm | **B** |
| XY resolution | 10 µm | 25 µm (point spacing) | 5 µm | **B** |
| Z resolution | 5 µm | 15 µm | 2.5 µm | **B** |
| `Timestep` | — | 1.0e-5 s | 2.5e-6 s | **B** |
| `Tracking` | None | **None (required)** | **Surface (required)** | **A** |
| Snapshot count | 20 | — | — | **B** |
| `Buffer` | 0.003 m | 0.0015 m | 0.001 m | **C** |

Resolution criterion: melt width ≈ 86 µm needs ≥ 6 cells; melt depth ≈ 45 µm.
Timestep criterion: the spot advances 8 / 2 µm per step.

### Why B's domain has both `Custom` and a regular XY block

3DThesis stops at `Custom` and computes only the point file. But `lpbf_score`'s
`generate_domain_xy()` needs `X.Min/Max` to build the target-region reference
grid — without it, coverage cannot be computed. Both blocks coexist; verified
that 3DThesis still outputs exactly one row per point.

---

## 6. Coverage points and scoring

| Parameter | Value | Grade | Source |
|---|---|---|---|
| Hatch spacing (fed to scorer) | 100 µm **explicit** | **A** | **Never use inference** — with 119 contour segments, `infer_hatch_spacing()` collapses to 21.33 µm |
| In-plane sampling | 25 µm | **B** | Width 86 µm and hatch 100 µm leave an unmelted ridge only ~12 µm wide; a 50 µm grid takes 2 points per 100 µm period — not enough to resolve it |
| z planes | 0 / −30 µm | **B** | Top surface / interlayer interface |
| Corridor half-width | 50 µm | **B** | = 0.5 × max(hatch, grid) |
| Point count | 80,498 | — | 40,249 XY × 2 planes |
| LOF criterion | `(h/W)² + (t/D)² ≤ 1` | **A** | Tang, Pistorius & Beuth, *Addit. Manuf.* (2017) |
| Keyhole criterion | aspect 0.5 / 0.8 | **A** | Common in the field; Cunningham et al., *Science* (2019) |
| Seven weights | 20/30/10/10/10/12/8 | **C** | **Expert judgement, no data.** Sensitivity study required |
| Pool-consistency scale | 0.35 | **C** | Engineering default |
| Thermal-consistency scale | 0.35 | **C** | Engineering default |

### Algorithmic bugs found and fixed

1. `infer_hatch_spacing()` was polluted by contour segments → added a
   length-weighted dominant-direction filter: 21.33 → 99.63 µm.
2. `match_cli_solid_geometry()` failed to match → added `infill_segments()` which
   groups by speed to isolate infill: 213 → exactly 94 segments. Now matches
   layer 1600 and excludes 2,372 cavity points.
3. `snapshot_geometry()` used the wrong grid (B's 25/15 µm instead of the
   snapshots' 10/5 µm) → now inferred from the data. LOF moved from
   "crosses the boundary, review" to "clearly exceeds" (nominal 0.810 → 1.071).
4. The same wrong-grid bug survived in `minimum_cells_across_p10_pool_dimension`
   because only one of the two call sites was fixed: 3.08 → 9.24 cells. The
   grid derivation is now hoisted so both sites share it.
5. Melt depth was quantised to the z grid (`z.max() - z.min() + dz`) → replaced
   with sub-grid linear interpolation of the liquidus isosurface, per column.
   Verified on real data: all 20 snapshots interpolate 100% of columns.
6. `generate_coverage_points.py` returned **exit code 1 on success** — the
   wrapper was `raise SystemExit(main())` and `main()` returns a statistics dict.
   Any caller treated a successful run as a failure.

27 unit tests cover all of the above.

---

## 7. One-page summary

### Exact (grade A — safe for a methods section)

Power, spot, layer thickness, hatch spacing, layer rotation, part geometry,
selected layer, scan path — all reproducible from vendor data or the files.
`T_L` = 1708 K and `ρ` = 7904 kg/m³ — measured on the same material batch.
LOF and keyhole criteria — primary literature.
Three source-level constraints — located at specific line numbers.

### Derived (grade B — the rule is defensible)

`Width_X` = 0.6124 × D4σ — conversion verified against the kernel.
`c` = 789 J/(kg·K) — from NIST measured enthalpy, latent heat included.
Grid and timestep — explicit resolution criteria (≥ 6 cells across the pool).
Coverage sampling — explicit physical reason (12 µm unmelted ridge).

### Compromise (grade C — findable, not yet looked up)

`T_0` (→ print job settings) · contour speed (→ build strategy) ·
jump speed and delay (→ machine configuration).
`k` = 24 W/(m·K) (two estimates differ by 15%, no same-source measurement).
The seven scoring weights (expert judgement).

**Half a day of desk work moves the first three from C to A.**

### Uncalibrated (grade D — current values not trustworthy)

`Efficiency` = 0.35 · `Depth_Z` = 3.5e-5

**Only these two — but they determine melt depth, and melt depth determines the
entire LOF verdict.**

Quantified: the model puts the critical hatch at 81 µm; *Metals* 11 (2021) 832
measured 99.8% density at hatch **0.12 mm** under otherwise identical parameters
(200 W / 800 mm/s / 0.03 mm layer). **That 1.5× gap is the magnitude of the
uncalibrated error.**

---

## 8. Confidence tiers

**Trustworthy** (insensitive to the D-grade parameters):

- Top-surface coverage ~100%, no track breaks
- Unmelted material concentrated 25–50 µm off the track centreline, deep
  (100% → 99.8% → 75.0% → 79.8% with depth)
- Layer thickness 50 → 30 µm raises interface fusion from 43% to 89%
- Pool consistency CV = 0.033, remelt controlled

**Not trustworthy** (directly dependent on the D-grade parameters):

- The FAIL verdict itself (88.8% vs the 95% threshold)
- "Hatch should be 70 µm"
- Any absolute temperature (peak 6442 K vs boiling point 3090 K)

### External validity check

Published work on the same model class (heat conduction, 316L bare plate, 200+
single tracks — [doi:10.1007/s40192-024-00382-2](https://doi.org/10.1007/s40192-024-00382-2))
establishes a validity indicator: **up to a simulated peak temperature of about
6000 °C, melt-depth error stays below 20%; beyond it the model under-predicts
depth, with the error growing linearly.**

This layer's 20 snapshots peak at 5852–6169 °C (median 5975). **Nine of twenty
exceed the threshold** — the case sits exactly on the validity boundary, and the
error direction is *under-prediction*, meaning the real pool is deeper.

Flipping FAIL to PASS at h = 100 µm requires melt depth to be **27% larger**.
The published in-range error bound is 20% and growing past the threshold.
**The FAIL verdict therefore lies within the model's own documented error band**
and cannot be treated as a physical conclusion.

---

## 9. Order of work

1. **Today:** re-run case B (~19 min) so the new material properties take effect.
2. **Half a day, no cost:** look up `T_0`, contour speed and jump parameters —
   three parameters move from C to A.
3. **Zero-cost external validation:** fix 200 W / 0.8 m/s / 30 µm layer and sweep
   hatch ∈ {70, 80, 100, 120} µm. If the curve passes near 120 µm the model
   agrees with the literature; if it needs 70 µm, that is direct evidence the
   model is conservative.
4. **Then:** 3–5 single-track coupons + cross-section metallography → fit the two
   D-grade parameters with `C_calib`.

Step 3 comes before the experiments because it is cheap and its outcome
determines how the metallography should be interpreted.
