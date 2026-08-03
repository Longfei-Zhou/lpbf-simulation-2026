# Cluster deployment

The full-chain driver is `source/LPBF_Agent/sweep_hatch.py`, invoked through
`lpbf.py sweep` or `lpbf.py layer`. Everything after slicing lives in it.
Netfabb slicing is the only step that stays on your workstation.

Step-by-step commands: [`RUNBOOK.md`](RUNBOOK.md). This document covers what is
specific to a shared machine.

---

## 1. 3DThesis must be rebuilt

A binary built on macOS is Mach-O arm64 and will not run on a Linux server. Ship
the source and compile there:

```bash
cd 3DThesis && mkdir -p build && cd build
cmake -D CMAKE_BUILD_TYPE=Release .. && make -j
```

| Dependency | Version | Required |
|---|---|---|
| CMake | 3.9+ | no — a plain makefile also ships |
| **OpenMP** | — | **yes** — the only hard dependency |

g++ on Linux normally bundles OpenMP. If cmake cannot find it:

```bash
sudo apt install libomp-dev     # Debian / Ubuntu
sudo yum install libgomp        # RHEL / CentOS
```

**Do not enable MPI.** `CMakeLists.txt` has `find_package(MPI QUIET)`, so it is
picked up automatically if present — but `Run.cpp:690` disables melt-pool
statistics output under MPI domain decomposition, and the calibration case needs
them. Single-node OpenMP is enough for this workload.

---

## 2. Python

**Python 3.10+** is required (the code uses `X | None`). Three third-party
packages:

```bash
pip install numpy pandas pyyaml
pip install matplotlib          # optional, plots only
```

Everything else is standard library.

---

## 3. Thread count — the easiest thing to get wrong

The case templates hard-code `MaxThreads` to the workstation's core count
(8 for case A, 13 for case B). **Leaving that unchanged wastes every extra core
on the server** — 3DThesis parallelises with OpenMP and `MaxThreads` caps the
thread pool directly.

```bash
python3 lpbf.py sweep ... --max-threads 64
```

This rewrites `Settings.txt` while assembling each case. Check the machine first:

```bash
nproc
lscpu | grep -E "^CPU\(s\)|Core|Socket"
```

Use the **physical** core count rather than the SMT count. For memory-bandwidth
bound work like this, hyperthreading rarely helps.

---

## 4. What to transfer

| Content | Size | Note |
|---|---|---|
| `source/Test/CLI/*.cli` | **200 MB** | Input slices — the bulk of the transfer |
| Case configuration | < 1 MB | Templates. **Do not copy `Data/`** |
| `source/LPBF_Agent/` | a few MB | Scoring code |
| `source/CLI_Trans.py`, `lpbf.py` | small | |
| `3DThesis/` source | a few MB | Recompiled on arrival |

```bash
rsync -av --exclude='Data/' --exclude='__pycache__/' --exclude='*.pyc' \
      --exclude='venv/' --exclude='.venv/' --exclude='build/' \
      ./ user@server:~/lpbf/
```

`Data/` in an existing workspace runs to ~570 MB and regenerates completely.

---

## 5. Disk

Per case, measured on the reference layer:

| Case | Output |
|---|---|
| A — snapshots (20 frames, 587,707 points each) | ~**360 MB** |
| B — solidification (80,498 points, 16 intermediate dumps) | ~**200 MB** |

Three hatch points ≈ **1.7 GB**; four ≈ **2.2 GB**. Add the same again per
extra layer. Reserving 20 GB is comfortable.

To trade history for space, raise `OutputFrequency` in case B's `Mode.txt` (e.g.
`2000` → `100000`) so only the Final file is written — about 1/16 the volume.
The cost is losing the evolution of coverage as the scan proceeds.

---

## 6. Pre-flight, in order

```bash
# 1. Build
cd 3DThesis/build && cmake -D CMAKE_BUILD_TYPE=Release .. && make -j
ls -la bin/3DThesis

# 2. Python
python3 --version                       # 3.10+
python3 -c "import numpy, pandas, yaml; print('ok')"

# 3. Self-test (27 tests, seconds)
cd ~/lpbf && python3 lpbf.py test

# 4. Input validation and case assembly (seconds, no simulation)
python3 lpbf.py sweep \
  --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --cli 80=source/Test/CLI/80µm_3dbenchy.cli \
  --cli 120=source/Test/CLI/120µm_3dbenchy.cli \
  --layer-number -1 --dry-run

# 5. One point first (~40 min) — confirms 3DThesis accepts the generated config
python3 lpbf.py sweep --cli 70=source/Test/CLI/70µm_3dbenchy.cli \
  --layer-number -1 --thesis-bin 3DThesis/build/bin/3DThesis \
  --max-threads <cores> --resume

# 6. Then the rest
```

Steps 1–4 are seconds to minutes. **Run them before submitting a queued job.**

---

## 7. Background execution

```bash
nohup python3 lpbf.py sweep ... --resume > sweep.log 2>&1 &
tail -f sweep.log
```

Slurm:

```bash
#!/bin/bash
#SBATCH --job-name=lpbf_sweep
#SBATCH --nodes=1                 # no MPI — one node
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
```

**Always pass `--resume`.** After a timeout kill, resubmitting the same script
continues from where it stopped — finished paths, coverage points and 3DThesis
output are not recomputed.

---

## 8. Two notes

**Plot labels.** Servers rarely have CJK fonts installed. The plotting code
probes for one and falls back to English labels if none is found, so you get a
readable figure rather than boxes. To keep the localised labels:

```bash
sudo apt install fonts-noto-cjk
```

**Never use `--skip-validation`.** Step 4 is the only thing standing between you
and a set of slices that are not actually comparable. It has already caught two
real export errors. Skipping it means hours of compute producing data you cannot
use.
