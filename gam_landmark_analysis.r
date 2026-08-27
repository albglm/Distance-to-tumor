################################################################################
# GAM LANDMARK ANALYSIS OF DISTANCE-FROM-TUMOR PROFILES
################################################################################
# WHAT IT DOES
#   For each requested MRI metric (e.g. R2*, QSM, ...) and each requested
#   distance-map variant, and for every participant:
#     1. Loads the metric map, the geodesic distance-from-tumor map, and a
#        tumor segmentation, and builds a table of the valid tissue voxels:
#        how far each is from the tumor, and its metric value.
#     2. Fits a smooth curve (GAM) of the metric vs. distance, locates the
#        first landmark peak/trough near the normal-appearing white matter (NAWM) reference
#        value, bootstraps a 95% CI for its location, and compares the
#        curve against simpler models.
#     3. Saves one summary row per participant to a CSV, plus prints
#        group-level summary tables.
##
# DEPENDENCIES
#   install.packages(c("RNifti", "data.table", "mgcv", "furrr", "future"))
#   RNifti reads/writes NIfTI (neuroimaging) files; swap in oro.nifti if
#   you prefer it, adjusting the RNifti calls in extract_voxel_data().
#
# REQUIRED INPUT FILES (per subject, NIfTI)
#   - 1 metric map 
#   - 1 geodesic distance map
#   - 1 BraTS-style tumor segmentation (multi-label: necrosis, T2/FLAIR hyperintensity [T2H], contrast enhancing tumor -- see https://github.com/BrainLesion/BraTS)
#   - 1 white matter mask to calculate reference mean values used throughout the landmark selection
#   - 1 participant list file.
#
# OPTIONAL INPUT FILES
#   - 1 additional exclusion mask can be added. Defaults to unset (skipped)
#
# OUTPUT
#   <output_dir>/<measure>_<map_label>_GAM.csv, one row per participant,
#   for every (measure, map_label) combination requested below.
#
# HOW TO RUN
#   1. Edit the USER SETTINGS block below to match your paths/metrics.
#   2. Run: Rscript gam_landmark_pipeline.R  (or Source in RStudio).
################################################################################

library(RNifti)     # reads/writes NIfTI (neuroimaging) files
library(data.table)
library(mgcv)
library(furrr)
library(future)

options(future.globals.maxSize = 16 * 1024^3)

# ============================== USER SETTINGS ================================

## --- parallelization -----------------------------------------------------------
# Number of parallel workers for the GAM-fitting step (the slow part).
#   n_workers <- 1   -- run sequentially, one participant at a time.
#   n_workers <- 4    -- run 4 participants at a time, etc.
n_workers <- min(8, future::availableCores())

## --- study / paths -----------------------------------------------------------
base_dir         <- "/data/project/TP1"  # root folder for this timepoint
output_dir       <- file.path(base_dir, "res")              # where per-participant result CSVs are saved
participant_list <- fread("/data/project/part.csv", sep = "\t")
participant_ids  <- participant_list[["NP"]]

scalar_measures     <- c("R2s", "ChiDia", "ChiPara", "QSM")  # MRI metrics to analyze; one output file per measure
map_labels          <- c("iso")   # distance-map variant tag(s), used to build the distance-map filename below

## --- optional additional mask ---------------------------------------------------
# Defaults to NULL (skipped). To use it, set it to a function that takes a
# subject ID and returns a file path to a binary NIfTI mask.
excluded_region_mask_fn <- NULL
# Example: excluded_region_mask_fn <- function(subject) file.path(base_dir, "derivatives/GRE", subject, "reg/t1_first_seg_all_fast_firstseg_bin.nii.gz")  # TRUE = deep grey matter, removed

## --- voxel extraction settings ------------------------------------------------
scalar_clip_percentiles   <- c(0, 99)   # metric values outside this percentile range (within valid tissue) are dropped
distance_clip_percentiles <- c(1, 95)   # distance values outside this percentile range are dropped
save_qc_masks <- FALSE   # optional: also save each participant's "valid voxel" mask as a NIfTI file, for visual QC

## --- landmark / GAM settings --------------------------------------------------
landmark_type <- "peak"      # "peak" (curve turns from rising to falling) or
                              # "trough" (curve turns from falling to rising) --
                              # which kind of turning point to search for

n_bootstrap         <- 100   # bootstrap resamples per patient (uncertainty estimate)
n_grid_points       <- 200   # points along the distance axis used to evaluate the fitted curve
distance_quantile   <- 1     # if < 1, keep only the nearest quantile of voxel distances per patient
bootstrap_voxel_cap <- NULL  # optional cap on voxels resampled per bootstrap iteration

# --- curve-flexibility (k) settings -----------------------------------------
k_min                       <- 10   # minimum allowed smoothness parameter k
k_max_default               <- 30   # maximum k tried on the first (default) fit
k_max_fallback              <- 45   # larger k tried only if the first fit looks inadequate
voxels_per_k                <- 15   # voxels required to "earn" one extra unit of k
smoothing_penalty           <- 1.5  # gamma; >1 favors smoother/simpler curves
k_index_adequate_threshold  <- 0.8  # k-index below this = fit may be too rigid
edf_k_ratio_flag_threshold  <- 0.9  # edf/k above this = curve using up most of its flexibility

# --- NAWM-proximity tolerance for landmark selection ------------------------
# Strategy: take the FIRST peak (nearest to tumor), whose fitted value is within nawm_tol_start of
# mean_nawm. If none qualifies, widen by nawm_tol_step and retry, up to
# nawm_tol_max_steps times, before falling back to the peak globally
# closest to mean_nawm. Units are the metric's own scalar units, so
# re-tune per metric:
#   metric    tol_start   tol_step
#   adc       0.00002     0.00002
#   fa        0.05        0.01
#   r2s       1.5         0.5
#   QSM       0.002       0.001
#   ChiPara   0.002       0.001
#   ChiDia    0.002       0.001
nawm_tol_start     <- 0.002
nawm_tol_step      <- 0.001
nawm_tol_max_steps <- 5

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (n_workers <= 1) {
  plan(sequential)
} else {
  plan(multisession, workers = n_workers)
}
tryCatch(
  future::nbrOfWorkers(),  # forces workers to actually start now, so failures surface here, not mid-run
  error = function(e) {
    message("multisession workers failed to start (", conditionMessage(e), "); falling back to sequential.")
    plan(sequential)
  }
)

# ============================ INPUT FILE PATHS =================================
# One function per required file, each taking a subject ID and returning
# its file path. Edit these directly to match wherever your files actually
# are.
 
metric_file_fn <- function(subject) file.path(base_dir, "derivatives/GRE", subject, paste0(measure, ".nii"))
  # the quantitative MRI metric map being analyzed (e.g. R2*, QSM...)
 
distance_map_file_fn <- function(subject) file.path(base_dir, "derivatives/DWI", subject,
  sprintf("distancemap_core_%s_dwi.nii.gz", map_label))
  # geodesic distance-from-tumor map
 
tumor_segmentation_file_fn <- function(subject) file.path(base_dir, "derivatives/segm", subject, "segmentation_brats.nii.gz")
  # BraTS-style multi-label tumor segmentation (necrosis / T2H / contrast ehnancing)
 
white_matter_mask_file_fn <- function(subject) file.path(base_dir, "derivatives/DWI", subject, "reg/t1_fast_pve_2_dwi_thr.nii.gz")
  # TRUE = white matter (e.g. thresholded white-matter partial volume
  # estimate ["pve"] from tissue segmentation).

# =============================== VOXEL EXTRACTION ===============================
# Builds a per-voxel table (one row per valid tissue voxel) for one subject,
# one metric, and one distance-map variant. This feeds straight into the GAM
# fitting step further below
extract_voxel_data <- function(subject, measure, map_label) {
  metric_img   <- readNifti(metric_file_fn(subject))   
  metric_map   <- as.array(metric_img)
  distance_map <- as.array(readNifti(distance_map_file_fn(subject)))
  tumor_seg    <- as.array(readNifti(tumor_segmentation_file_fn(subject)))
  white_matter_mask <- as.array(readNifti(white_matter_mask_file_fn(subject))) > 0

  # From the tumor segmentation: exclude necrosis and contrast enhancing tumor from the tissue analysis;
  # T2H is kept, and used below both to label NAWM voxels and to widen the allowed
  # region
  tumor_exclusion_mask <- !(tumor_seg %in% c(1, 3)) # 1 = necrosis; 2 = T2H; 3 = contrast enhancing tumor
  t2h_mask <- tumor_seg == 2
  allowed_region <- white_matter_mask | t2h_mask

  # Voxels kept for analysis: finite, at a positive distance from the
  # tumor, not necrosis/hyperintensity, and white matter or T2H.
  valid <- is.finite(metric_map) & is.finite(distance_map) & (distance_map > 0) &
    tumor_exclusion_mask & allowed_region

  # Optional further restriction, only applied if configured above.
  if (!is.null(excluded_region_mask_fn)) {
    excluded_region_mask <- as.array(readNifti(excluded_region_mask_fn(subject))) > 0
    valid <- valid & !excluded_region_mask
  }

  # Trim extreme outliers (in both the metric and the distance) before use.
  scalar_bounds   <- quantile(metric_map[valid], scalar_clip_percentiles / 100, na.rm = TRUE)
  distance_bounds <- quantile(distance_map[valid], distance_clip_percentiles / 100, na.rm = TRUE)

  valid_mask <- valid &
    metric_map > scalar_bounds[1] & metric_map < scalar_bounds[2] &
    distance_map > distance_bounds[1] & distance_map < distance_bounds[2]

  if (save_qc_masks) {
    qc_file <- file.path(dirname(metric_file_fn(subject)), sprintf("%s_valid_mask.nii.gz", subject))
    valid_mask_image <- asNifti(array(as.integer(valid_mask), dim = dim(valid_mask)), reference = metric_img)
    writeNifti(valid_mask_image, qc_file)
    message("  saved QC mask: ", qc_file)
  }

  # Distance, rescaled to 0-1 across this subject's valid voxels.
  dist_values   <- distance_map[valid_mask]
  d_min <- min(dist_values); d_max <- max(dist_values)
  distance_norm <- pmin(pmax((dist_values - d_min) / (d_max - d_min), 0), 1)

  # Metric value, z-scored across this subject's valid voxels.
  scalar_values <- metric_map[valid_mask]
  scalar_norm   <- (scalar_values - mean(scalar_values)) / sd(scalar_values)

  # Any valid voxel that isn't T2H is NAWM
  region <- ifelse(tumor_seg[valid_mask] == 2, "t2h", "nawm")

  data.table(
    participant   = subject, region = region,
    distance_raw  = dist_values, distance_norm = distance_norm,
    scalar_raw    = scalar_values, scalar_norm = scalar_norm,
    min_distance  = d_min, max_distance = d_max
  )
}

# =============================== GAM HELPERS ====================================

# Drops patients with too little data (<30 voxels); optionally trims to the
# nearest `distance_quantile` fraction of distances.
prep_voxel_data <- function(voxel_data, quantile_cutoff = distance_quantile) {
  if (is.null(voxel_data) || nrow(voxel_data) < 30L) return(NULL)
  if (quantile_cutoff < 1.0) {
    dmax <- quantile(voxel_data$distance_raw, quantile_cutoff, na.rm = TRUE)
    voxel_data <- voxel_data[distance_raw <= dmax]
  }
  if (nrow(voxel_data) < 30L) NULL else voxel_data
}

# Fits scalar_raw ~ smooth(distance_raw) with the given smoothness k.
# Returns NULL (instead of crashing) if the fit fails.
fit_gam_safe <- function(voxel_data, k_value, method = "REML", gamma = smoothing_penalty) {
  tryCatch(
    gam(scalar_raw ~ s(distance_raw, bs = "cr", k = k_value),
        data = voxel_data, method = method, gamma = gamma),
    error = function(e) NULL
  )
}

# Numerically estimates the slope of the fitted curve at each point on
# `distance_grid`, by comparing predictions a tiny step apart.
get_slope <- function(gam_model, distance_grid) {
  tryCatch({
    eps <- diff(range(distance_grid)) / 1e5
    lo  <- predict(gam_model, newdata = data.frame(distance_raw = distance_grid))
    hi  <- predict(gam_model, newdata = data.frame(distance_raw = distance_grid + eps))
    (hi - lo) / eps
  }, error = function(e) NULL)
}

# Finds every location where the slope crosses zero in the requested
# direction, at positive distances only:
#   type = "peak"   -> curve turns from rising to falling (slope + to -)
#   type = "trough" -> curve turns from falling to rising (slope - to +)
# Controlled by the `landmark_type` setting above.
find_all_peaks <- function(slope, distance_grid, type = landmark_type) {
  if (is.null(slope) || length(slope) < 2L) return(numeric(0))
  n <- length(slope)
  crossing_idx <- switch(type,
    peak   = which(slope[-n] > 0 & slope[-1L] <= 0 & distance_grid[-n] > 0),
    trough = which(slope[-n] < 0 & slope[-1L] >= 0 & distance_grid[-n] > 0),
    stop('landmark_type must be "peak" or "trough", got: ', type)
  )
  if (!length(crossing_idx)) return(numeric(0))
  vapply(crossing_idx, function(i)
    distance_grid[i] - slope[i] * (distance_grid[i + 1L] - distance_grid[i]) /
      (slope[i + 1L] - slope[i]), numeric(1))
}

# Average slope of the curve from the tumor up to `peak_loc` -- summarizes
# how steeply the curve rose on its way to that landmark.
mean_slope_to_peak <- function(slope, distance_grid, peak_loc) {
  if (is.na(peak_loc) || is.null(slope)) return(NA_real_)
  idx <- distance_grid <= peak_loc
  if (sum(idx) < 2L) return(NA_real_)
  mean(slope[idx], na.rm = TRUE)
}

# Average measured value in this patient's "nawm" voxels -- the reference value.
get_mean_nawm <- function(voxel_data) {
  if (!"region" %in% names(voxel_data)) return(NA_real_)
  vals <- voxel_data$scalar_raw[voxel_data$region == "nawm"]
  if (!length(vals)) return(NA_real_)
  mean(vals, na.rm = TRUE)
}

# Farthest distance from the tumor covered by this patient's T2H voxels -- used to restrict where peaks are first searched for.
get_t2h_max_dist <- function(voxel_data) {
  if (!"region" %in% names(voxel_data)) return(Inf)
  vals <- voxel_data$distance_raw[voxel_data$region == "t2h"]
  if (!length(vals)) return(Inf)
  max(vals, na.rm = TRUE)
}

# Implements the peak-selection strategy described under "NAWM-proximity
# tolerance" above. Returns the chosen location, the number of candidate
# peaks considered, a text code for which rule matched, and the tolerance
# ultimately used.
locate_landmark_near_ref <- function(gam_model, slope, distance_grid, mean_nawm,
                                      search_max = Inf,
                                      tol_start = nawm_tol_start,
                                      tol_step = nawm_tol_step,
                                      max_steps = nawm_tol_max_steps) {
  peaks <- find_all_peaks(slope, distance_grid)
  n_peaks_total <- length(peaks)
  if (!n_peaks_total) return(list(loc = NA_real_, n_peaks = 0L, selection = NA_character_, tol_used = NA_real_))

  peaks_in_t2h_range <- peaks[peaks <= search_max]
  used_fallback_range  <- length(peaks_in_t2h_range) == 0L
  candidates <- if (used_fallback_range) peaks else peaks_in_t2h_range
  n_peaks <- length(candidates)

  if (is.na(mean_nawm)) {
    sel <- if (used_fallback_range) "first_no_nawm_ref_no_t2h_range" else "first_no_nawm_ref"
    return(list(loc = candidates[1L], n_peaks = n_peaks, selection = sel, tol_used = NA_real_))
  }

  peak_values <- tryCatch(
    predict(gam_model, newdata = data.frame(distance_raw = candidates)),
    error = function(e) rep(NA_real_, n_peaks)
  )
  if (all(is.na(peak_values))) {
    sel <- if (used_fallback_range) "first_predict_failed_no_t2h_range" else "first_predict_failed"
    return(list(loc = candidates[1L], n_peaks = n_peaks, selection = sel, tol_used = NA_real_))
  }

  dist_to_reference <- abs(peak_values - mean_nawm)

  # Progressively widen the tolerance; at each step, take the FIRST
  # (nearest-tumor) candidate within tolerance -- so widening only rescues
  # a near-tumor peak that just missed a tighter window, never jumping
  # straight to a far peak while a tighter match might still exist.
  for (step in 0:max_steps) {
    tolerance <- tol_start + step * tol_step
    within_tolerance <- which(dist_to_reference <= tolerance)
    if (length(within_tolerance)) {
      idx <- within_tolerance[1L]
      sel <- if (idx == 1L) "first_within_nawm_tol" else "later_within_nawm_tol"
      if (step > 0L) sel <- paste0(sel, "_widened", step)
      sel <- if (used_fallback_range) paste0(sel, "_no_t2h_range") else sel
      return(list(loc = candidates[idx], n_peaks = n_peaks, selection = sel, tol_used = tolerance))
    }
  }

  # No match at any tolerance -- fall back to whichever candidate is
  # globally closest to mean_nawm, flagged distinctly.
  closest_idx <- which.min(dist_to_reference)
  sel <- if (closest_idx == 1L) "closest_to_nawm_and_first_no_tol_match"
         else "closest_to_nawm_not_first_no_tol_match"
  sel <- if (used_fallback_range) paste0(sel, "_no_t2h_range") else sel
  list(loc = candidates[closest_idx], n_peaks = n_peaks, selection = sel, tol_used = NA_real_)
}

# Runs mgcv's built-in check of whether k was large enough for the curve
# to be considered adequately flexible.
check_k_adequacy <- function(gam_model) {
  tryCatch({
    kc <- k.check(gam_model)
    list(k_index = kc[1L, "k-index"], k_pval = kc[1L, "p-value"])
  }, error = function(e) list(k_index = NA_real_, k_pval = NA_real_))
}

# Chooses an initial k from sample size, fits the curve, refits with a
# larger k if the adequacy check flags it (see USER SETTINGS above),
# computes the slope, and locates the landmark peak.
select_and_fit <- function(voxel_data) {
  distance_grid <- seq(min(voxel_data$distance_raw), max(voxel_data$distance_raw), length.out = n_grid_points)
  n_voxels <- nrow(voxel_data)

  k_value <- min(k_max_default, max(k_min, floor(n_voxels / voxels_per_k)))

  gam_model <- fit_gam_safe(voxel_data, k_value, gamma = smoothing_penalty)
  if (is.null(gam_model)) return(list(mod = NULL, k = NA_integer_, k_reason = "fit_failed"))

  k_check <- check_k_adequacy(gam_model)
  edf <- sum(gam_model$edf)
  edf_k_ratio <- edf / k_value
  k_reason <- "fixed"

  if (!is.na(k_check$k_index) && k_check$k_index < k_index_adequate_threshold &&
      !is.na(edf_k_ratio) && edf_k_ratio > edf_k_ratio_flag_threshold &&
      k_value < k_max_fallback) {

    k_value_fallback <- k_max_fallback
    gam_model_fallback <- fit_gam_safe(voxel_data, k_value_fallback, gamma = smoothing_penalty)

    if (!is.null(gam_model_fallback)) {
      gam_model   <- gam_model_fallback
      k_value     <- k_value_fallback
      k_check     <- check_k_adequacy(gam_model)
      edf         <- sum(gam_model$edf)
      edf_k_ratio <- edf / k_value
      k_reason    <- "flagged_low_kindex_increased"
    } else {
      k_reason <- "flagged_low_kindex_refit_failed"
    }
  }

  slope <- get_slope(gam_model, distance_grid)
  mean_nawm <- get_mean_nawm(voxel_data)
  t2h_max_dist <- get_t2h_max_dist(voxel_data)
  landmark <- locate_landmark_near_ref(gam_model, slope, distance_grid, mean_nawm, t2h_max_dist)

  list(
    mod = gam_model, k = k_value, k_index = k_check$k_index, k_pval = k_check$k_pval,
    edf = edf, n_peaks = landmark$n_peaks,
    fp = landmark$loc, landmark_selection = landmark$selection, tol_used = landmark$tol_used,
    mean_nawm = mean_nawm, t2h_max_dist = t2h_max_dist,
    dv = slope, grid = distance_grid, gamma = smoothing_penalty, k_reason = k_reason
  )
}

# Resamples this patient's voxels with replacement `n_boot` times, refits
# the curve and re-locates the landmark peak each time, building up a
# distribution of plausible peak locations (used for the CI).
bootstrap_patient <- function(voxel_data, k_value, distance_grid, gamma = smoothing_penalty,
                               mean_nawm = NA_real_, search_max = Inf, n_boot = n_bootstrap) {
  n_voxels <- nrow(voxel_data)
  sample_n <- if (!is.null(bootstrap_voxel_cap) && bootstrap_voxel_cap < n_voxels) bootstrap_voxel_cap else n_voxels
  peak_boot <- rep(NA_real_, n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample(n_voxels, sample_n, replace = (sample_n == n_voxels))
    gam_model <- fit_gam_safe(voxel_data[idx], k_value, gamma = gamma)
    if (is.null(gam_model)) next
    slope <- get_slope(gam_model, distance_grid)
    if (is.null(slope)) next
    peak_boot[b] <- locate_landmark_near_ref(gam_model, slope, distance_grid, mean_nawm, search_max)$loc
  }
  peak_boot
}

# Summarizes the bootstrap distribution of peak locations into a median
# and 95% CI (2.5th / 97.5th percentiles).
summarize_bootstrap <- function(peak_boot) {
  ci <- quantile(peak_boot, c(0.025, 0.975), na.rm = TRUE)
  list(median = round(median(peak_boot, na.rm = TRUE), 4),
       ci_lo = round(ci[1L], 4), ci_hi = round(ci[2L], 4))
}

# Tests whether the fitted smooth curve explains the data meaningfully
# better than just the patient's overall mean value (a "flat line" model).
compare_smooth_to_mean <- function(voxel_data, k_value, distance_grid, gam_model) {
  mean_val <- mean(voxel_data$scalar_raw, na.rm = TRUE)

  out <- list(
    mean_val = round(mean_val, 4),
    aic_null = NA_real_, aic_smooth = NA_real_, delta_aic = NA_real_,
    lrt_p = NA_real_, l2_curve_vs_mean = NA_real_, curve_range = NA_real_,
    smooth_better = NA
  )

  tryCatch({
    null_model <- lm(scalar_raw ~ 1, data = voxel_data)
    out$aic_null   <- round(AIC(null_model), 2)
    out$aic_smooth <- round(AIC(gam_model), 2)
    out$delta_aic  <- round(out$aic_null - out$aic_smooth, 2)
  }, error = function(e) NULL)

  tryCatch({
    null_ml <- gam(scalar_raw ~ 1, data = voxel_data, method = "ML")
    full_ml <- gam(scalar_raw ~ s(distance_raw, bs = "cr", k = k_value),
                    data = voxel_data, method = "ML", gamma = smoothing_penalty)
    lrt_result <- anova(null_ml, full_ml, test = "Chisq")
    p_col <- grep("Pr\\(", names(lrt_result), value = TRUE)
    out$lrt_p <- round(lrt_result[[p_col]][2L], 6)
  }, error = function(e) NULL)

  tryCatch({
    grid_pred <- predict(gam_model, newdata = data.frame(distance_raw = distance_grid))
    out$l2_curve_vs_mean <- round(mean((grid_pred - mean_val)^2, na.rm = TRUE), 6)
    out$curve_range      <- round(diff(range(grid_pred, na.rm = TRUE)), 4)
  }, error = function(e) NULL)

  out$smooth_better <- isTRUE(!is.na(out$lrt_p) && out$lrt_p < 0.05 &&
                                 sum(gam_model$edf) > 1.5)
  out
}

# Tests whether adding the smooth distance curve explains the data
# meaningfully better than discrete tissue-region labels alone (needs >=2
# distinct regions present).
compare_smooth_to_region <- function(voxel_data, k_value) {
  out <- list(
    n_regions = NA_integer_, aic_region = NA_real_, aic_region_smooth = NA_real_,
    delta_aic_region = NA_real_, lrt_p_region = NA_real_, smooth_better_region = NA
  )

  if (!"region" %in% names(voxel_data)) return(out)

  region_factor <- factor(voxel_data$region)
  out$n_regions <- nlevels(region_factor)
  if (out$n_regions < 2L) return(out)

  voxel_data2 <- copy(voxel_data)
  voxel_data2[, region := region_factor]

  tryCatch({
    null_ml <- gam(scalar_raw ~ factor(region), data = voxel_data2, method = "ML")
    full_ml <- gam(scalar_raw ~ factor(region) + s(distance_raw, bs = "cr", k = k_value),
                    data = voxel_data2, method = "ML", gamma = smoothing_penalty)

    out$aic_region        <- round(AIC(null_ml), 2)
    out$aic_region_smooth <- round(AIC(full_ml), 2)
    out$delta_aic_region  <- round(out$aic_region - out$aic_region_smooth, 2)

    lrt_result <- anova(null_ml, full_ml, test = "Chisq")
    p_col <- grep("Pr\\(", names(lrt_result), value = TRUE)
    out$lrt_p_region <- round(lrt_result[[p_col]][2L], 6)

    out$smooth_better_region <- isTRUE(!is.na(out$lrt_p_region) && out$lrt_p_region < 0.05)
  }, error = function(e) NULL)

  out
}

# Full per-patient pipeline: clean data, fit the curve and pick k, locate
# & bootstrap the landmark peak, run the model comparisons, and package
# everything into one output row. Patients with too little data or a
# failed fit get a minimal row flagging why.
compute_patient <- function(patient_id, metric_table, measure, map_label) {

  voxel_data <- prep_voxel_data(metric_table[.(patient_id)])
  if (is.null(voxel_data))
    return(data.table(participant = patient_id, metric = measure, map = map_label,
                      selected_k = NA_integer_, k_reason = "insufficient_data",
                      has_peak = FALSE))

  fit <- select_and_fit(voxel_data)
  if (is.null(fit$mod))
    return(data.table(participant = patient_id, metric = measure, map = map_label,
                      selected_k = NA_integer_, k_reason = fit$k_reason,
                      has_peak = FALSE))

  boot_summary <- summarize_bootstrap(
    bootstrap_patient(voxel_data, fit$k, fit$grid, fit$gamma, fit$mean_nawm, fit$t2h_max_dist))
  mean_slope <- mean_slope_to_peak(fit$dv, fit$grid, fit$fp)
  mean_cmp   <- compare_smooth_to_mean(voxel_data, fit$k, fit$grid, fit$mod)
  region_cmp <- compare_smooth_to_region(voxel_data, fit$k)

  data.table(
    participant        = patient_id, metric = measure, map = map_label,
    selected_k         = fit$k, k_reason = fit$k_reason,
    gamma              = fit$gamma,
    deviance_explained = round(summary(fit$mod)$dev.expl, 4),
    edf                = round(fit$edf, 2),
    edf_k_ratio        = round(fit$edf / fit$k, 3),
    k_index            = round(fit$k_index, 3),
    k_pval             = round(fit$k_pval, 4),
    first_peak         = round(fit$fp, 4),
    landmark_selection = fit$landmark_selection,
    tol_used           = round(fit$tol_used, 4),
    mean_nawm          = round(fit$mean_nawm, 4),
    t2h_max_dist     = round(fit$t2h_max_dist, 4),
    has_peak           = !is.na(fit$fp),
    n_peaks            = fit$n_peaks,
    fp_boot_median     = boot_summary$median,
    fp_boot_ci_lo      = boot_summary$ci_lo,
    fp_boot_ci_hi      = boot_summary$ci_hi,
    mean_deriv_to_peak = round(mean_slope, 6),
    mean_val           = mean_cmp$mean_val,
    aic_null           = mean_cmp$aic_null,
    aic_smooth         = mean_cmp$aic_smooth,
    delta_aic          = mean_cmp$delta_aic,
    lrt_p              = mean_cmp$lrt_p,
    l2_curve_vs_mean   = mean_cmp$l2_curve_vs_mean,
    curve_range        = mean_cmp$curve_range,
    smooth_better      = mean_cmp$smooth_better,
    n_regions            = region_cmp$n_regions,
    aic_region           = region_cmp$aic_region,
    aic_region_smooth    = region_cmp$aic_region_smooth,
    delta_aic_region     = region_cmp$delta_aic_region,
    lrt_p_region         = region_cmp$lrt_p_region,
    smooth_better_region = region_cmp$smooth_better_region)
}

# ================================ MAIN LOOP ==================================
# For each distance-map variant and each metric: extract voxel data for
# every participant directly from the NIfTI files, fit the GAM landmark
# analysis in parallel, print summary tables, and save one CSV of results.

for (map_label in map_labels) {
  for (measure in scalar_measures) {

    message(sprintf("=== metric: %s | distance-map variant: %s ===", measure, map_label))

    message("extracting voxel data from NIfTI files...")
    voxel_tables <- vector("list", length(participant_ids))
    for (i in seq_along(participant_ids)) {
      subject <- participant_ids[i]
      voxel_tables[[i]] <- tryCatch(
        extract_voxel_data(subject, measure, map_label),
        error = function(e) {
          message(sprintf("  %s failed during extraction: %s", subject, e$message))
          NULL
        })
    }
    metric_table <- rbindlist(Filter(Negate(is.null), voxel_tables), fill = TRUE)
    setkey(metric_table, participant)
    message(sprintf("  %d voxels across %d participants", nrow(metric_table), uniqueN(metric_table$participant)))
    print(head(metric_table, 10))
    rm(voxel_tables)

    message("running gam fits across patients...")
    t0 <- proc.time()

    results_all <- future_map(
      participant_ids,
      function(patient_id)
        tryCatch(
          compute_patient(patient_id, metric_table, measure, map_label),
          error = function(e) {
            message(sprintf("  %s failed: %s", patient_id, e$message))
            NULL
          }),
      .options  = furrr_options(seed = TRUE, globals = TRUE),
      .progress = TRUE)
    names(results_all) <- participant_ids

    message(sprintf("done in %.1f min", (proc.time() - t0)["elapsed"] / 60))

    n_failed <- sum(vapply(results_all, is.null, logical(1)))
    if (n_failed) message(sprintf("%d / %d patients failed", n_failed, length(participant_ids)))

    results_per_map <- rbindlist(Filter(Negate(is.null), results_all), fill = TRUE)
    rm(metric_table, results_all); gc()

    message(sprintf("%d rows", nrow(results_per_map)))

    # Curve fit quality across all patients.
    results_per_map[!is.na(selected_k), .(
      median_dev = round(median(deviance_explained, na.rm = TRUE), 3),
      q25_dev    = round(quantile(deviance_explained, 0.25, na.rm = TRUE), 3),
      q75_dev    = round(quantile(deviance_explained, 0.75, na.rm = TRUE), 3),
      median_edf = round(median(edf, na.rm = TRUE), 2),
      median_edf_k_ratio = round(median(edf_k_ratio, na.rm = TRUE), 3),
      n_fixed    = sum(k_reason == "fixed", na.rm = TRUE),
      n_flagged_increased = sum(k_reason == "flagged_low_kindex_increased", na.rm = TRUE),
      n_flagged_failed     = sum(k_reason == "flagged_low_kindex_refit_failed", na.rm = TRUE)
    )] |> print()

    # Peak/landmark detection rate, location, and tolerance-widening frequency.
    results_per_map[, .(
      prop_has_peak   = round(mean(has_peak, na.rm = TRUE), 3),
      median_peak     = round(median(first_peak, na.rm = TRUE), 3),
      median_boot_med = round(median(fp_boot_median, na.rm = TRUE), 3),
      median_ci_width = round(median(fp_boot_ci_hi - fp_boot_ci_lo, na.rm = TRUE), 3),
      median_n_peaks  = round(median(n_peaks, na.rm = TRUE), 1),
      prop_first_at_start_tol = round(mean(landmark_selection == "first_within_nawm_tol", na.rm = TRUE), 3),
      prop_needed_widening    = round(mean(grepl("_widened", landmark_selection), na.rm = TRUE), 3),
      prop_no_tol_match       = round(mean(grepl("no_tol_match", landmark_selection), na.rm = TRUE), 3),
      median_tol_used = round(median(tol_used, na.rm = TRUE), 4)
    )] |> print()

    # How often the smooth curve beats a flat-mean model.
    results_per_map[!is.na(selected_k), .(
      prop_smooth_better = round(mean(smooth_better, na.rm = TRUE), 3),
      median_delta_aic    = round(median(delta_aic, na.rm = TRUE), 2),
      median_lrt_p        = round(median(lrt_p, na.rm = TRUE), 4),
      median_l2           = round(median(l2_curve_vs_mean, na.rm = TRUE), 6),
      median_curve_range  = round(median(curve_range, na.rm = TRUE), 4)
    )] |> print()

    # How often the smooth curve beats a region-only model.
    results_per_map[!is.na(selected_k), .(
      n_with_2plus_regions = sum(n_regions >= 2, na.rm = TRUE),
      prop_smooth_better_region = round(mean(smooth_better_region, na.rm = TRUE), 3),
      median_delta_aic_region   = round(median(delta_aic_region, na.rm = TRUE), 2),
      median_lrt_p_region       = round(median(lrt_p_region, na.rm = TRUE), 4)
    )] |> print()

    # Per-patient landmark locations and their bootstrap CIs.
    results_per_map[!is.na(first_peak),
                   .(participant, first_peak, landmark_selection, n_peaks, fp_boot_median, fp_boot_ci_lo, fp_boot_ci_hi)] |> print()

    out_file <- file.path(output_dir, sprintf("%s_%s_GAM.csv", measure, map_label))
    fwrite(results_per_map, out_file)
    message("saved: ", out_file)

    rm(results_per_map); gc()
  }
}