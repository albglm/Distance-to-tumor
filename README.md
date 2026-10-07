# Distance-to-tumor qMRI profiles in glioblastoma

Code for modeling quantitative MRI (qMRI) metrics as continuous functions of
geodesic distance from the contrast-enhancing tumor (CET), as described in the
accompanying manuscript.

Patient data cannot be shared for privacy reasons; the code is provided so the
analysis can be inspected and applied to other datasets.

## Pipeline

| Step | Script | Output |
|---|---|---|
| 1. Distance maps (diffusion space) | `distance_maps.py` | `derivatives/distancemaps/sub-*/dwi/*_distance.nii.gz` |
| 2. GAM profiles, landmark, local gradient | `gam_distance_profiles.R` | `derivatives/distance_profiles/*.csv` |
| 3. Figures | `plot_profiles.R` | population profiles, participant profiles, example participants |
| 4. Summary tables | `summarize_profiles.R` | fit diagnostics and landmark/gradient summary per metric and map |
| 5. Sensitivity analyses and contralateral control (optional) | `compare_sensitivity.R` | landmark agreement of sensitivity runs with the main analysis; ipsilateral vs contralateral profiles and fit statistics |

Each qMRI map must be in the same space as its distance map and masks; any
registration between diffusion and other image spaces is done by the user
beforehand (see Input).

`distance_maps.py` solves the eikonal equation with the Hamiltonian Fast
Marching library from all CET voxels, within NAWM plus the lesion and
excluding deep gray matter, for the four map variants (isotropic, isotropic
weighted by axial diffusivity, anisotropic from FOD peaks, anisotropic
weighted), the contralateral control (mirrored CET) and the CET-dilation
sensitivity map. Its header lists the equations and the input files.

`gam_distance_profiles.R` fits, for every participant, qMRI metric and distance
map, a generalized additive model (GAM) of metric value against distance, and
extracts the landmark distance and local gradient. It also computes landmark
uncertainty (spatial block bootstrap), fit diagnostics and population curves.

## Requirements

- Python 3 with `numpy`, `nibabel`, `scipy`, and the
  [HamiltonFastMarching](https://github.com/Mirebeau/HamiltonFastMarching)
  library compiled locally (step 1)
- R ≥ 4.4 with `RNifti`, `data.table`, `mgcv`, `future`, `furrr` (step 2) and
  `ggplot2`, `patchwork` (steps 3-5)
- Preprocessed images per participant (see Input)

```
pip install numpy nibabel scipy
Rscript -e 'install.packages(c("RNifti", "data.table", "mgcv", "future", "furrr", "ggplot2", "patchwork"))'
```

## Input

The scripts expect a BIDS-derivatives layout. All images of one participant
must be in the space of the metric map: GRE-derived metrics in GRE space
(`anat`, `space-gre`), ADC and FA in diffusion space (`dwi`, `space-dwi`), so
masks and distance maps are needed in both spaces. `distance_maps.py` writes
the distance maps in diffusion space; bringing them (and the masks) into the
space of the other qMRI maps is left to the user, with any registration tool.
Segmentations and masks must be resampled with nearest-neighbour
interpolation.

```
<bids_root>/
├── participants.tsv                      participant_id, include (1 = analysed)
└── derivatives/
    ├── qmri/sub-P001/
    │   ├── anat/sub-P001_space-gre_R2starmap.nii.gz
    │   │        sub-P001_space-gre_Chimap.nii.gz                 (QSM)
    │   │        sub-P001_space-gre_desc-dia_Chimap.nii.gz        (Xdia)
    │   │        sub-P001_space-gre_desc-para_Chimap.nii.gz       (Xpara)
    │   └── dwi/ sub-P001_space-dwi_model-tensor_param-adc_dwimap.nii.gz
    │            sub-P001_space-dwi_model-tensor_param-fa_dwimap.nii.gz
    │            sub-P001_space-dwi_model-tensor_param-ad_dwimap.nii.gz    axial diffusivity (weighted maps)
    │            sub-P001_space-dwi_model-csd_param-peaks_dwimap.nii.gz    two largest FOD peaks (anisotropic maps)
    ├── masks/sub-P001/{anat,dwi}/
    │        sub-P001_space-<gre|dwi>_desc-tumor_dseg.nii.gz       tumor segmentation: 1 necrosis, 2 T2H, 3 CET
    │        sub-P001_space-<gre|dwi>_label-NAWM_mask.nii.gz       thresholded FAST white matter
    │        sub-P001_space-<gre|dwi>_label-tumorhemi_mask.nii.gz  hemisphere of the tumor
    │        sub-P001_space-dwi_label-deepGM_mask.nii.gz          deep gray matter, excluded (optional)
    │        sub-P001_space-dwi_label-CETmirrored_mask.nii.gz     mirrored CET (contralateral control)
    │        sub-P001_space-dwi_label-contrahemi_mask.nii.gz      contralateral hemisphere
    └── distancemaps/sub-P001/{anat,dwi}/
             sub-P001_space-<gre|dwi>_desc-<map>_distance.nii.gz
```

The tumor segmentation uses the BraTS 2023 labels (1 necrosis, 2 T2H/edema,
3 CET); other label values are set in `TUMOR_LABELS` (`distance_maps.py`) and
`tumor_labels` (`gam_distance_profiles.R`), e.g. `cet = 4` for BraTS 2021. The
seeds of the distance maps are CET + necrosis, and the T2H is the lesion
outside them.

Distance maps (`desc-` label): `iso` (isotropic), `isoweighted` (isotropic
weighted, primary), `aniso` (anisotropic), `anisoweighted` (anisotropic
weighted), `contraisoweighted` (isotropic weighted from the mirrored CET,
contralateral control) and `isoweightedcetdil` (isotropic weighted from the CET
dilated by two voxels). The same labels appear in the output file names.

Other maps can be tried without changing the rest of the pipeline:

- **Other seeds or weighting:** add one line to `MAPS` in `distance_maps.py`,
  e.g. `isoweightedlesion` (seeds = whole lesion, i.e. distance from the
  T2H–NAWM boundary), then use its name in `maps` of
  `gam_distance_profiles.R`. With lesion seeds the T2H voxels have distance 0
  and are left out, so the profile covers the NAWM only and the landmark is
  searched over the whole distance range.
- **Other propagation model:** edit `speed()` (weighting by axial diffusivity)
  or `anisotropy()` (speed along and across the fiber direction from the FOD
  peak ratio) in the block PROPAGATION MODEL of `distance_maps.py`. These
  functions apply to all maps, so add a new entry to `MAPS` (e.g. a copy of
  `aniso` named `anisotest`) and run only that one (`MAPS_TO_RUN`); existing
  maps are kept.

The minimum for one analysis is a qMRI map, the tumor segmentation, a NAWM
mask and (with `restrict_to_tumor_hemisphere <- TRUE`) the tumor-hemisphere
mask, plus the distance map. The plain `iso` map needs no diffusion data; the
weighted maps need the axial diffusivity and the anisotropic maps the FOD peaks.

The qMRI metrics are listed in `metric_table` at the top of
`gam_distance_profiles.R`: one row per metric with its file name, image space,
landmark type (peak or trough) and NAWM tolerance. To analyse another metric,
add a row; `metrics` then selects which rows are run. All file names are
defined in the block **PATHS AND INPUT FILES**, so only this block needs to
change for a different folder structure.

## Example data

`make_example_data.py` writes a small synthetic dataset (three participants,
every input file, and distance maps under all map names) to check the
installation before using real data:

```
python make_example_data.py /path/to/example_bids
```

Then set `bids_root` in `gam_distance_profiles.R` and `results_dir` in
`plot_profiles.R` to that folder. For a quick run (under a minute) also set
`maps <- "isoweighted"` and `n_bootstrap <- 0`:

```
Rscript gam_distance_profiles.R
Rscript plot_profiles.R
```

The landmarks should be at about 5 mm in all three participants. Step 1 needs
the HFM library; set `RECOMPUTE = True` in `distance_maps.py` to test it on
these data.

## Running

```
python distance_maps.py                      # step 1, all included participants
Rscript gam_distance_profiles.R              # step 2, all included participants
Rscript plot_profiles.R                      # step 3
Rscript summarize_profiles.R                 # step 4
Rscript compare_sensitivity.R                # step 5, after the sensitivity / contralateral runs

python distance_maps.py sub-P001             # one or more participants only
Rscript gam_distance_profiles.R sub-P001
```

Without arguments, both scripts process every participant with `include = 1`
in `participants.tsv`; participant IDs given on the command line override the
table. Both scripts keep what is already computed, so an interrupted run
restarts where it stopped, and a missing file or failed fit in one participant
is reported without stopping the run. Set `RECOMPUTE = True`
(`distance_maps.py`) or `recompute <- TRUE` (`gam_distance_profiles.R`) to
redo existing results, e.g. after changing input files or masks.
The summary tables of `gam_distance_profiles.R` (`*_GAM.csv`, `*_curves.csv`,
`*_popcurve.csv`) contain the participants of the last call; running it again
without arguments rebuilds them for the whole cohort from the cache.

The default settings reproduce the **main analysis** (six metrics, four
distance maps, block bootstrap with B = 100).

The other analyses of the paper are obtained by changing the settings below
and running the scripts again:

| Analysis | Settings |
|---|---|
| Main analysis | defaults |
| Contralateral control | step 1 with `MAPS_TO_RUN = ["contraisoweighted"]`; step 2 with `maps <- "contraisoweighted"`, `contralateral <- TRUE` |
| Smoothing sensitivity | `smoothing_penalty <- 1.0`, `n_bootstrap <- 0` |
| Landmark-tolerance sensitivity | `tol_scale <- 0.5` or `2`, `n_bootstrap <- 0` |
| T2H-boundary sensitivity | `t2h_perturb_voxels <- -1L` or `1L`, `n_bootstrap <- 0` |
| CET-boundary sensitivity | step 1 with `MAPS_TO_RUN = ["isoweightedcetdil"]`; step 2 with `maps <- "isoweightedcetdil"`, `n_bootstrap <- 0` |

Each run is saved under its own `run_tag`, so runs do not overwrite each other.
`compare_sensitivity.R` then compares the landmarks of the chosen runs (setting
`compare`) with the main analysis: Spearman ρ with bootstrap CI, median absolute
landmark shift and the share of participants whose landmark changed.
With `run_contralateral <- TRUE` it also compares the main run with the
contralateral run: a figure of the population profiles of both sides and their
paired difference, and a table of deviance explained, AIC gain per 1,000 voxels
and mean white-matter value (paired Wilcoxon test, FDR across metrics).
`ipsilateral_voxels <- "nawm"` refits the ipsilateral profiles on NAWM voxels
only (tissue-matched control) from the voxel cache. Participants listed in
`excluded_participants` (e.g. tumor spread into the contralateral hemisphere)
are left out of the tests and summarized separately.

## Output

For each metric × distance map × run, in `<bids_root>/derivatives/distance_profiles/`:

- `*_GAM.csv` — one row per participant: voxel counts, k, k-index, edf,
  deviance explained, AIC gain, Moran's I of the residuals, landmark
  (`first_peak` from the original fit; `fp_boot_median` and 95% CI from the
  block bootstrap), selection step of the landmark, local gradient
  (`mean_deriv_to_peak`).
- `*_curves.csv` — fitted curve of each participant on a normalized 0–1
  distance grid (raw and mean-centered).
- `*_popcurve.csv` — population curve (equal-weight mean of participant curves)
  with participant-bootstrap 95% CI.

## Method summary

- **Sampling domain:** NAWM plus the whole T2H (taken irrespective of tissue
  class), excluding the CET; metric values above the 99th percentile and
  distances outside the 1st–95th percentile are removed per participant.
- **GAM:** `metric ~ s(distance, bs = "cr", k)`, REML, γ = 1.5;
  k = min(30, max(10, ⌊n/15⌋)), raised to 45 if the k-index is < 0.8 and
  edf/k > 0.9.
- **Landmark:** first peak (trough for ADC) within the T2H distance range
  (all extrema if none lies there) whose fitted value is within a tolerance of
  the NAWM mean; the tolerance is widened in up to five steps, with the
  extremum closest to the NAWM mean as final fallback.
- **Local gradient:** mean derivative from the smallest fitted distance to the
  landmark.
- **Uncertainty:** spatial block bootstrap (5 × 5 × 5 mm cubes, B = 100),
  conditional on the segmentation and distance map.

## Citation

If you use this code, please cite the accompanying article:

> [Authors]. [Title]. [Journal] [Year]. doi:[DOI]

and the Hamiltonian Fast Marching library used for the distance maps:

- J.-M. Mirebeau, J. Portegies. Hamiltonian Fast Marching: a numerical solver
  for anisotropic and non-holonomic eikonal PDEs. *Image Processing On Line*
  2019;9:47–93. doi:[10.5201/ipol.2019.227](https://doi.org/10.5201/ipol.2019.227)
