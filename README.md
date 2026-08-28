# Peritumoral Distance-Resolved qMRI Analysis

Computes geodesic distance-from-tumor maps from diffusion MRI, then models each
quantitative MRI (qMRI) metric as a smooth function of that distance, per patient, to
locate a data-driven landmark distance at which the metric recovers toward
normal-appearing white matter (NAWM) values. Supports both single-timepoint analysis
and longitudinal (pre- vs. post-treatment) comparison, with a plotting script to
visualize per-patient fits, cohort-level trends, and group comparisons.

The pipeline steps are independent: regenerate distance maps with different settings
without touching the R scripts, or re-run either R analysis script across metrics/map
variants without recomputing distance maps.

## Pipeline

1. **`distmap_script.py`** — Computes a geodesic distance-from-tumor map per patient
   (Hamiltonian Fast Marching). Isotropic/anisotropic, weighted/unweighted variants.
   - In: seed mask, co-registered (+ FOD peaks and axial-diffusivity map for anisotropic/weighted variants)
   - Out: one distance map (NIfTI) per patient per map variant

2. **`gam_landmark_pipeline.R`** — Per patient, fits a smooth curve of metric vs.
   distance, locates the landmark peak/trough nearest NAWM, bootstraps its location,
   compares against simpler models. Loops over metrics and map variants.
   - In: metric map, distance map (from step 1), BraTS-style tumor segmentation (necrosis / edema / enhancing), white-matter mask — all co-registered
   - Out: one CSV per (metric, map variant), one row per patient: landmark distance, bootstrap CI, curve-fit diagnostics, model comparisons

3. **`gam_landmark_prepost_pipeline.R`** — For patients with both a pre- and a
   post-treatment scan, jointly fits a smooth curve per timepoint, matches the
   post-treatment landmark to the pre-treatment one, and bootstraps both. Use this
   instead of step 2 when comparing the same patient across two timepoints. Loops over
   metrics and map variants.
   - In: metric map, distance map, tumor segmentation, white-matter mask — all at BOTH timepoints, co-registered
   - Out: two CSVs per (metric, map variant): one row per patient per timepoint (curve fit + landmark), and one row per patient (pre-to-post landmark shift)

4. **`gam_landmark_plots.R`** — Builds figures for one metric and one map variant at a
   time: a per-patient grid (fitted curve, NAWM/T2H reference lines, landmark + CI), a
   population overview (all patients' curves plus the cohort median and IQR), and
   cohort boxplots of landmark distance and gradient, optionally split by a subgroup
   column. Reads the results CSV from step 2, and re-extracts voxels from the same
   NIfTI files only to refit each curve for display — no analysis is redone.
   - In: results CSV (from step 2), plus the same metric map/distance map/tumor segmentation/white-matter mask used to produce it
   - Out: PDF figures (patient grid, population overview, boxplots)

All scripts also need a participant list file. Exact filenames/paths are placeholders
in each script's settings — edit to match your data. See each script's header for full
detail on inputs, settings, and outputs.

## Requirements

* Python: `numpy`, `nibabel`, plus [HamiltonFastMarching](https://github.com/Mirebeau/HamiltonFastMarching) (external toolbox, installed locally).
* R: `install.packages(c("RNifti", "data.table", "mgcv", "furrr", "future"))` for the analysis scripts (2, 3); add `ggplot2` for the plotting script (4).

## Citation

If you use this pipeline, please cite [paper/preprint] and [HamiltonFastMarching](https://github.com/Mirebeau/HamiltonFastMarching).
