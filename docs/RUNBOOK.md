# Runbook

Every command from compilation to scoring. **Only the three variables in step 0
need editing** — everything after uses them, so no path is ever typed twice.

---

## 0. Set three variables

Re-run these in every new shell (or put them in `~/.bashrc`).

```bash
export PROJ=$HOME/lpbf                           # project root
export THESIS=$PROJ/3DThesis/build/bin/3DThesis  # built in step 2
export NCORES=$(nproc)                           # macOS: sysctl -n hw.physicalcpu
```

Check them:

```bash
echo "$PROJ" && ls "$PROJ/lpbf.py" && echo "cores: $NCORES"
```

If `ls` fails, `PROJ` is wrong. Fix it before continuing.

---

## 1. Transfer (run on your workstation)

```bash
rsync -av \
  --exclude='Data/' --exclude='__pycache__/' --exclude='*.pyc' \
  --exclude='venv/' --exclude='.venv/' --exclude='build/' \
  ./ user@server:~/lpbf/
```

`Data/` holds simulation output that regenerates; `build/` is a host-specific
binary that must be recompiled on the target.

---

## 2. Build 3DThesis

```bash
cd "$PROJ/3DThesis"
rm -rf build && mkdir build && cd build
cmake -D CMAKE_BUILD_TYPE=Release ..
make -j"$NCORES"
```

Verify:

```bash
ls -la "$THESIS" && "$THESIS" 2>&1 | head -5
```

Running it without arguments prints usage or complains about a missing
`ParamInput` — **anything other than "No such file" means the build worked.**

<details>
<summary>cmake cannot find OpenMP</summary>

```bash
sudo apt install libomp-dev      # Debian / Ubuntu
sudo yum install libgomp         # RHEL / CentOS
```

OpenMP is the only hard dependency. **Do not install MPI** — `Run.cpp:690`
disables melt-pool statistics under MPI domain decomposition, and the
calibration case needs them.
</details>

---

## 3. Python environment

```bash
python3 --version                                  # 3.10+ required
pip install numpy pandas pyyaml matplotlib
python3 -c "import numpy, pandas, yaml; print('deps OK')"
```

`matplotlib` is optional — without it plots are skipped and nothing else changes.
For CJK labels: `sudo apt install fonts-noto-cjk`. Without a CJK font the plot
falls back to English labels automatically rather than rendering empty boxes.

---

## 4. Self-test (seconds)

```bash
cd "$PROJ" && python3 lpbf.py test
```

Expect `Ran 27 tests ... OK`.

---

## 5. Validate inputs (seconds — do not skip)

```bash
cd "$PROJ"
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 \
  --dry-run
```

**You must see:**

```
Same layer, same section, same contours across all inputs; only hatch differs
-- controlled-comparison precondition holds.
```

This check has already caught two real export errors: a hatch angle off by one
layer rotation, and a hatch distance entered as 1.0 mm instead of 0.1 mm.
**Never bypass it with `--skip-validation`** — the hours of simulation that
follow would produce data that cannot be compared.

Steps 4 and 5 take under ten seconds combined. Run them **before** submitting a
queued job; they eliminate most wasted queue time.

---

## 6. One point first (~40 min)

3DThesis has not yet consumed the `Mode.txt` and `Domain.txt` this pipeline
generates. Validate that with a single point rather than betting hours on it.

```bash
cd "$PROJ"
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --layer-number -1 \
  --thesis-bin "$THESIS" \
  --max-threads "$NCORES" \
  --resume
```

Then check three things exist:

```bash
ls "$PROJ/source/Test/316L/sweep_hatch/h70/A_snapshots/Data" | wc -l   # 20
ls "$PROJ/source/Test/316L/sweep_hatch/h70/B_solidification/Data"      # Solidification.Final
head -40 "$PROJ/source/LPBF_Agent/results/sweep_hatch/h70/Assessment_Report.md"
```

All three present → the chain works end to end.

---

## 7. Full run (~2 h for three points)

```bash
cd "$PROJ"
nohup python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 \
  --thesis-bin "$THESIS" \
  --max-threads "$NCORES" \
  --resume \
  > sweep.log 2>&1 &

tail -f sweep.log     # Ctrl-C stops watching, not the job
```

`h=70` is already done; `--resume` skips it and only fills in 80 and 120.

<details>
<summary>Slurm</summary>

```bash
cat > sweep.sbatch <<'EOF'
#!/bin/bash
#SBATCH --job-name=lpbf_sweep
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=06:00:00
#SBATCH --mem=32G

cd "$SLURM_SUBMIT_DIR"
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 \
  --thesis-bin "$PWD/3DThesis/build/bin/3DThesis" \
  --max-threads "$SLURM_CPUS_PER_TASK" \
  --resume
EOF

sbatch sweep.sbatch
```

One node is enough — 3DThesis does not use MPI here. If the job is killed on
timeout, resubmit the same file; `--resume` picks up where it stopped.
</details>

**Always pass `--resume`.** Completed paths, coverage points and 3DThesis output
are never recomputed.

**Always pass `--max-threads`.** Without it the case templates keep their
hard-coded thread counts (8 and 13) and the rest of the machine sits idle.

---

## 8. Results

```bash
cd "$PROJ/source/LPBF_Agent/results/sweep_hatch" && ls
cat Hatch_Sweep_Report.md
```

| File | Contents |
|---|---|
| `Hatch_Sweep_Report.md` | **Read this.** §3 gives the critical hatch and its ratio to the literature's 120 µm |
| `hatch_sweep_lof.png` | LOF vs. hatch, with threshold line and literature marker |
| `hatch_sweep.csv` | Raw data |
| `h70/` `h80/` `h120/` | Per-point scoring output |

Re-render the report without re-simulating:

```bash
cd "$PROJ" && python3 lpbf.py sweep --cli 70=... --cli 80=... --cli 120=... --report-only
```

---

## 9. Adding the fourth point

Export the 100 µm slice with hatch distance **0.1 mm** and nothing else changed
(see [`SLICING.md`](SLICING.md)), then add it to the same command:

```bash
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 100=source/Test/CLI/100µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 \
  --thesis-bin "$THESIS" --max-threads "$NCORES" --resume
```

`--resume` skips the three finished points.

---

## 10. Other layers

```bash
cd "$PROJ"

# Pick a layer
python3 lpbf.py layers --cli source/Test/CLI/100µm_3dbenchy.cli --top 10

# Estimate cost -- seconds, do not skip
python3 lpbf.py probe --cli source/Test/CLI/100µm_3dbenchy.cli --layer-z-mm 12.750

# Run. Add --region-mm 6 above ~500k points.
python3 lpbf.py layer \
  --cli source/Test/CLI/100µm_3dbenchy.cli \
  --layer-z-mm 12.750 \
  --region-mm 6 \
  --thesis-bin "$THESIS" --max-threads "$NCORES" --resume
```

The largest sections in this part are near z ≈ 12.6–13.6 mm: 2.04 million
coverage points and roughly 8 hours for the whole layer, versus ~116,000 points
and ~27 minutes with `--region-mm 6` — at the same in-plane resolution.

---

## Quick reference

```bash
export PROJ=$HOME/lpbf
export THESIS=$PROJ/3DThesis/build/bin/3DThesis
export NCORES=$(nproc)

cd "$PROJ/3DThesis" && rm -rf build && mkdir build && cd build
cmake -D CMAKE_BUILD_TYPE=Release .. && make -j"$NCORES"

pip install numpy pandas pyyaml matplotlib

cd "$PROJ"
python3 lpbf.py test
python3 lpbf.py sweep --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
                      --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
                      --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
                      --layer-number -1 --dry-run
python3 lpbf.py sweep --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
                      --layer-number -1 --thesis-bin "$THESIS" \
                      --max-threads "$NCORES" --resume
python3 lpbf.py sweep --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
                      --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
                      --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
                      --layer-number -1 --thesis-bin "$THESIS" \
                      --max-threads "$NCORES" --resume
cat "$PROJ/source/LPBF_Agent/results/sweep_hatch/Hatch_Sweep_Report.md"
```
