# LPBF Agent

This package screens one LPBF layer from 3DThesis output. One 3DThesis case is
one physical build layer; CSV `z` values are simulated depths below that layer,
not multiple build layers.

The repository-level [README](../../README.md) is the canonical installation and
command reference.

## Structure

```text
LPBF_Agent/
├── run_all.py                  scoring entry point
├── generate_coverage_points.py pre-simulation Custom-point generator
├── lpbf_score/                 production scoring package
├── config/                     thresholds, weights, and output templates
├── docs/                       formulas, evidence, and model boundaries
├── tests/                      automated tests
└── results/                    assessment outputs
```

`lpbf_score/scorer.py` performs the following work:

- reads final Solidification, Snapshot, RDF, and optional temperature-history data;
- reconstructs process parameters from Material, Domain, Beam, Path, and CLI files;
- measures liquidus-isotherm width, length, depth, connectivity, and resolution;
- separates top coverage, interface fusion, and LOF geometry to avoid duplicate penalties;
- reports nominal, half-cell, and conservative LOF values;
- reports data completeness separately from evidence adequacy;
- matches the executed Path to a CLI layer so holes, cavities, and exterior voids
  are excluded while unscanned solid remains in the coverage denominator;
- classifies unmelted points by interior, endpoint, boundary, and centerline distance;
- includes reproducible G/V, G×V, cooling-rate, solidification-event, and remelt
  diagnostics without treating them as release criteria;
- writes human-readable summaries, diagnostics, plots, CSV tables, and JSON;
- supports intentional large offsets with `--stress-test` while preserving hard
  physical criteria and cross-case provenance checks.

## Coverage-point generation

Build `Path.txt` and a regular X/Y/Z Domain first, then run:

```bash
python3 generate_coverage_points.py \
  --case-dir /path/to/current_layer_solidification \
  --domain-file /path/to/regular/Domain.txt \
  --source-cli /path/to/original_slice.cli \
  --layer-thickness-um 30 \
  --hatch-spacing-um 100 \
  --output /path/to/current_layer_solidification/coverage_target_points.txt
```

`--domain-file` must preserve the regular X/Y/Z bounds and resolution. The
physical powder-layer thickness must come from `--source-cli` or
`--layer-thickness-um`; `Z.Res` is only numerical grid spacing. The generated
file is headerless, space-delimited, and expressed in millimetres for a Custom
3DThesis domain:

```text
Custom
{
    File coverage_target_points.txt
}
```

With a source CLI, the complete even-odd solid cross-section is the denominator.
The path corridor is retained as a missing-scan diagnostic, so a solid region
with no path cannot disappear from the assessment.

## Scoring

```bash
python3 run_all.py \
  --solidification-dir /path/to/316L_solidification \
  --snapshots-dir /path/to/316L_snapshots \
  --source-cli /path/to/original_slice.cli \
  --layer-id 316L-layer-001 \
  --output-dir results/316L_layer_001
```

Use explicit values when the source CLI is unavailable:

```bash
python3 run_all.py \
  --solidification-dir /path/to/316L_solidification \
  --snapshots-dir /path/to/316L_snapshots \
  --layer-thickness-um 30 \
  --hatch-spacing-um 100 \
  --stress-test \
  --output-dir results/316L_extreme_assessment
```

Stress-test mode accepts intentional coordinate and parameter offsets as valid
test inputs. It does not disable coverage, LOF, aspect-ratio, or provenance
checks, and it never combines physically inconsistent cases into one score.

## Output

Open `Dashboard.html` first. `00_READ_ME_FIRST.md` is the concise conclusion and
`Assessment_Report.md` is the technical appendix.

```text
results/316L_layer_001/
├── Dashboard.html
├── 00_READ_ME_FIRST.md
├── Layer_Assessment.csv
├── Action_Plan.csv
├── Assessment_Report.md
├── Problem_Diagnosis.csv
├── Coverage_Diagnostics.csv
├── Score_Breakdown.csv
├── assessment.json
├── feature_summary.json
├── snapshot_features.csv
└── coverage_target_points.txt
```

`PASS` means the simulation passed this screening model; it is not CT,
metallography, density, mechanical-property, or fatigue certification.

## Tests

```bash
python3 -m unittest discover -s tests -v
```

The tests cover good geometry, LOF failure, missing evidence, authoritative
coverage failure, resolution sensitivity, CLI hole masking, and extreme Path
line parsing.
