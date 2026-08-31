# Peritumoral Distance-Resolved qMRI Analysis

Computes geodesic distance-from-tumor maps from diffusion MRI, then models each
quantitative MRI (qMRI) metric as a smooth function of that distance, per patient, to
locate a data-driven landmark distance at which the metric recovers toward
normal-appearing white matter (NAWM) values. Supports both single-timepoint analysis
and longitudinal comparison across timepoints (e.g. pre- vs. post-treatment), with plotting scripts to
visualize per-patient fits, cohort-level trends, and group comparisons.

## Structure

One shared first step, then two independent analysis paths — use whichever matches
your data (or both):

```
distmap_script.py
        |
        +-- gam_landmark_analysis.R --------- gam_figures.R          (single timepoint)
        +-- gam_landmark_analysis_long.R --- gam_figures_long.R (longitudinal)
```

Each of these three pieces is independent in its own right: regenerate distance maps
with different settings without touching the R scripts, re-run either R analysis
script across metrics/map variants without recomputing distance maps, and re-run a
plotting script without redoing the analysis.

### Distance map

**`distmap_script.py`** — Computes a geodesic distance-from-tumor map per patient
(Hamiltonian Fast Marching). Isotropic/anisotropic, weighted/unweighted variants.
Feeds into both analysis paths below.
- In: seed mask, co-registered (+ FOD peaks and axial-diffusivity map for anisotropic/weighted variants)
- Out: one distance map (NIfTI) per patient per map variant

### Single-timepoint analysis

**`gam_landmark_analysis.R`** — Per patient, fits a smooth curve of metric vs.
distance, locates the landmark peak/trough nearest NAWM, bootstraps its location,
compares against simpler models. Loops over metrics and map variants.
- In: metric map, distance map (from the shared step), BraTS-style tumor segmentation (necrosis / edema / enhancing), white-matter mask — all co-registered
- Out: one CSV per (metric, map variant), one row per patient: landmark distance, bootstrap CI, curve-fit diagnostics, model comparisons

**`gam_figures.R`** — Builds figures for one metric and one map variant at a
time: a per-patient grid (fitted curve, NAWM/T2H reference lines, landmark + CI), a
population overview (all patients' curves plus the cohort median and IQR), and cohort
boxplots of landmark distance and gradient, optionally split by a subgroup column.
Reads the results CSV above, and re-extracts voxels from the same NIfTI files only to
refit each curve for display — no analysis is redone.
- In: results CSV above, plus the same metric map/distance map/tumor segmentation/white-matter mask used to produce it
- Out: PDF figures (patient grid, population overview, boxplots)

### Longitudinal comparison

**`gam_landmark_analysis_long.R`** — For patients with two co-registered scans (a
follow-up visit, before/after any intervention, or any two timepoints), jointly fits a smooth curve
per timepoint, matches the second landmark to the first, and bootstraps both. Use this
instead of the single-timepoint script when comparing the same patient across two
scans. Loops over metrics and map variants.
- In: metric map, distance map, tumor segmentation, white-matter mask — all at BOTH timepoints, co-registered (distance map/segmentation/mask are typically computed once at the first timepoint and reused for the second)
- Out: two CSVs per (metric, map variant): one row per patient per timepoint (curve fit + landmark), and one row per patient (change between timepoints)

**`gam_figures_long.R`** — The two-timepoint counterpart to the single-timepoint
plotting script: a per-patient grid with both timepoints' curves overlaid (landmark +
CI for each, mismatched correspondences flagged), a population overview across
timepoints, a paired shift plot (each patient's landmark at each timepoint as a
connected line), and cohort boxplots of the landmark and gradient shift, optionally
split by a subgroup column. Reads the two results CSVs above, and re-extracts voxels
from the same NIfTI files (at both timepoints) only to refit each curve for display.
- In: per-timepoint and delta results CSVs above, plus the same NIfTI inputs used to produce them, at both timepoints
- Out: PDF figures (patient grid, population overview, paired shift plot, boxplots)

All scripts also need a participant list file. Exact filenames/paths are placeholders
in each script's settings — edit to match your data. See each script's header for full
detail on inputs, settings, and outputs.

## Requirements

* Python: `numpy`, `nibabel`, plus [HamiltonFastMarching](https://github.com/Mirebeau/HamiltonFastMarching) (external toolbox, installed locally).
* R: `install.packages(c("RNifti", "data.table", "mgcv", "furrr", "future"))` for the analysis scripts; add `ggplot2` for the plotting scripts.

## Citation

If you use this pipeline, please cite [paper/preprint] and [HamiltonFastMarching](https://github.com/Mirebeau/HamiltonFastMarching).
