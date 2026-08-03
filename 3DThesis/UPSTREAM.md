# Upstream Version

This directory is a distributable copy of
[ORNL-MDF/3DThesis](https://github.com/ORNL-MDF/3DThesis) with its `.git`
directory removed.

| Field | Value |
|---|---|
| Source | <https://github.com/ORNL-MDF/3DThesis.git> |
| Commit | `2de7fc6d8cfa5de78b111df97b1a4d9156a8cf60` |
| Description | `4.0.0-6-g2de7fc6` |
| Date | 2026-06-30 |
| Upstream change | Merge pull request #38 from ORNL-MDF/BeamRadiusFix |

To use a Git submodule instead:

```bash
rm -rf 3DThesis
git submodule add https://github.com/ORNL-MDF/3DThesis.git 3DThesis
cd 3DThesis
git checkout 2de7fc6d8cfa5de78b111df97b1a4d9156a8cf60
```

## Source-level constraints used by this project

These are implementation constraints, not project preferences:

1. The `depth` column requires `Tracking Surface`; it is unavailable for a
   Custom point domain.
2. `MP_Stats` is incompatible with a Custom domain because the implementation
   requires structured-grid `ijk_to_p` indexing.
3. `MP_depth` converts the Z-direction cell count using `xres` rather than
   `zres`. On an anisotropic grid it is wrong by `xres/zres`; use `depth`.
4. Melt-pool statistics are unavailable with MPI domain decomposition, so this
   project builds without MPI and uses single-node OpenMP parallelism.
