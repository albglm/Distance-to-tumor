################################################################################
# GAM LANDMARK ANALYSIS: PRE- VS POST-TREATMENT COMPARISON
################################################################################
# WHAT IT DOES
#   For each requested MRI metric and each requested distance-map variant,
#   and for every participant who has BOTH a pre-treatment and a
#   post-treatment scan:
#     1. Loads the metric map, the geodesic distance-from-tumor map, and a
#        tumor segmentation at EACH timepoint, and builds one combined
#        table of valid tissue voxels tagged "Pre" or "Post".
#     2. Fits a single smooth curve per timepoint (metric vs. distance,
#        jointly modeled with a time factor), finds the first landmark
#        peak/trough nearest the tumor in each curve, and matches the Post
#        landmark to the Pre one (flagging cases where they land
#        suspiciously far apart -- likely a mismatched correspondence
#        rather than genuine change).
#     3. Bootstraps both landmark locations together (resampling voxels
#        from both timepoints, refitting, and re-matching each time) for
#        95% CIs, and computes the Pre-to-Post shift.
#     4. Saves one row per participant per timepoint (curve fit &
#        landmark details) and one row per participant summarizing the
#        Pre-to-Post change, plus prints group-level summary tables.
#
#   This script simply takes the first peak/trough nearest the
#   tumor at Pre, then finds whichever Post peak/trough best corresponds
#   to it.
#
# DEPENDENCIES
#   install.packages(c("RNifti", "data.table", "mgcv", "furrr", "future"))
#   RNifti reads/writes NIfTI (neuroimaging) files; swap in oro.nifti if
#   you prefer it, adjusting the RNifti calls in extract_voxel_data().
#
# REQUIRED INPUT FILES (per subject, per timepoint, NIfTI)
#   - 1 metric map
#   - 1 geodesic distance map
#   - 1 BraTS-style tumor segmentation (multi-label: necrosis, T2/FLAIR
#     hyperintensity [T2H], contrast-enhancing tumor)
#   - 1 white matter mask (defines which voxels are analyzed, together
#     with the segmentation's T2H label)
#   - 1 participant list file (participants without both timepoints are
#     skipped automatically)
#
#   Pre and Post must be co-registered (same voxel grid, same space) --
#   this script does voxel-by-voxel comparisons across timepoints with no
#   registration step of its own.
#
# OPTIONAL INPUT FILES
#   - 1 additional exclusion mask can be added. Defaults to unset (skipped)
#
# OUTPUT, per (measure, map_label) combination
#   - <output_dir>/<measure>_<map_label>_GAM_prepost_per_timepoint.csv
#     one row per participant per timepoint (curve fit + landmark details)
#   - <output_dir>/<measure>_<map_label>_GAM_prepost_delta.csv
#     one row per participant (Pre-to-Post landmark shift)
#
# HOW TO RUN
#   1. Edit the USER SETTINGS block below to match your paths/metrics.
#   2. Run: Rscript gam_landmark_prepost_pipeline.R  (or Source in RStudio).
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
base_dir_pre     <- "/data/project/TP1"  # root folder for the pre-treatment timepoint
base_dir_post    <- "/data/project/TP2"  # root folder for the post-treatment timepoint
output_dir       <- file.path(base_dir_post, "res")   # where result CSVs are saved
participant_list <- fread("/data/project/part.csv", sep = "\t")
participant_ids  <- participant_list[["NP"]]

scalar_measures <- c("R2s", "ChiDia", "ChiPara", "QSM")  # MRI metrics to analyze; one pair of output files per measure
map_labels      <- c("iso")   # distance-map variant tag(s), used to build the distance-map filename below

## --- optional additional mask ---------------------------------------------------
# Defaults to NULL (skipped). To use it, set it to a function that takes a
# subject ID and a timepoint ("Pre"/"Post") and returns a file path to a
# binary NIfTI mask.
excluded_region_mask_fn <- NULL
# Example: excluded_region_mask_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/GRE", subject, "reg/t1_first_seg_all_fast_firstseg_bin.nii.gz")  # TRUE = deep grey matter, removed

## --- voxel extraction settings ------------------------------------------------
scalar_clip_percentiles   <- c(0, 99)   # metric values outside this percentile range (within valid tissue) are dropped
distance_clip_percentiles <- c(1, 95)   # distance values outside this percentile range are dropped
save_qc_masks <- FALSE   # optional: also save each participant's "valid voxel" mask as a NIfTI file, for visual QC

## --- landmark / GAM settings --------------------------------------------------
landmark_type <- "peak"      # "peak" (curve turns from rising to falling) or
                              # "trough" (curve turns from falling to rising) --
                              # which kind of turning point to search for

n_bootstrap         <- 100   # bootstrap resamples per patient (uncertainty estimate)
n_grid_points       <- 200   # points along the distance axis used to evaluate the fitted curves
bootstrap_voxel_cap <- NULL  # optional cap on voxels resampled per bootstrap iteration

# --- curve-flexibility (k) settings -----------------------------------------
# Voxel-count heuristic: pick an initial k from sample size,
# refit once with a larger k only if the fit
# looks inadequate (see select_k_and_fit below).
k_min                       <- 10   # minimum allowed smoothness parameter k
k_max_default               <- 30   # maximum k tried on the first (default) fit
k_max_fallback              <- 45   # larger k tried only if the first fit looks inadequate
voxels_per_k                <- 15   # voxels (pooled across both timepoints) required to "earn" one extra unit of k
k_index_adequate_threshold  <- 0.8  # k-index below this = fit may be too rigid, for either curve
edf_k_ratio_flag_threshold  <- 0.9  # edf/k above this = curve using up most of its flexibility

# --- Pre/Post peak-matching settings -----------------------------------------
# How far apart (in distance units) the Post landmark can land from the
# Pre landmark before this is flagged as a likely mismatch rather than
# genuine change, and the search window used when re-locating the Post
# landmark during bootstrap resampling (anchored near the Pre location so
# resampled Post landmarks stay correspondence-consistent).
peak_match_window <- 50

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
# One function per required file, each taking a subject ID and a timepoint
# ("Pre"/"Post") and returning its file path. Edit these directly to match
# wherever your files actually are. If the distance map/segmentation/white matter mask are only
# computed once, at baseline, symlink the expected Post-timepoint path to
# the baseline file.

time_base_dir <- function(time) if (time == "Pre") base_dir_pre else base_dir_post

metric_file_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/GRE", subject, paste0(measure, ".nii"))
  # the quantitative MRI metric map being analyzed (e.g. R2*, QSM...)

distance_map_file_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/DWI", subject,
  sprintf("distancemap_core_%s_dwi.nii.gz", map_label))
  # geodesic distance-from-tumor map

tumor_segmentation_file_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/segm", subject, "segmentation_brats.nii.gz")
  # BraTS-style multi-label tumor segmentation (necrosis / T2H / contrast enhancing)

white_matter_mask_file_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/DWI", subject, "reg/t1_fast_pve_2_dwi_thr.nii.gz")
  # TRUE = white matter (e.g. thresholded white-matter partial volume
  # estimate ["pve"] from tissue segmentation).

# =============================== VOXEL EXTRACTION ===============================
# Builds a per-voxel table (one row per valid tissue voxel) for one subject,
# one timepoint, one metric, and one distance-map variant. This feeds
# straight into the GAM fitting step further below.
extract_voxel_data <- function(subject, time, measure, map_label) {
  metric_img   <- readNifti(metric_file_fn(subject, time))
  metric_map   <- as.array(metric_img)
  distance_map <- as.array(readNifti(distance_map_file_fn(subject, time)))
  tumor_seg    <- as.array(readNifti(tumor_segmentation_file_fn(subject, time)))
  white_matter_mask <- as.array(readNifti(white_matter_mask_file_fn(subject, time))) > 0

  # From the tumor segmentation: exclude necrosis and contrast enhancing
  # tumor from the tissue analysis; T2H is kept, and widens the allowed
  # region.
  tumor_exclusion_mask <- !(tumor_seg %in% c(1, 3))  # 1 = necrosis; 2 = T2H; 3 = contrast enhancing tumor
  t2h_mask <- tumor_seg == 2
  allowed_region <- white_matter_mask | t2h_mask

  # Voxels kept for analysis: finite, at a positive distance from the
  # tumor, not necrosis/enhancing tumor, and white matter or T2H.
  valid <- is.finite(metric_map) & is.finite(distance_map) & (distance_map > 0) &
    tumor_exclusion_mask & allowed_region

  # Optional further restriction, only applied if configured above.
  if (!is.null(excluded_region_mask_fn)) {
    excluded_region_mask <- as.array(readNifti(excluded_region_mask_fn(subject, time))) > 0
    valid <- valid & !excluded_region_mask
  }

  # Trim extreme outliers (in both the metric and the distance) before use.
  scalar_bounds   <- quantile(metric_map[valid], scalar_clip_percentiles / 100, na.rm = TRUE)
  distance_bounds <- quantile(distance_map[valid], distance_clip_percentiles / 100, na.rm = TRUE)

  valid_mask <- valid &
    metric_map > scalar_bounds[1] & metric_map < scalar_bounds[2] &
    distance_map > distance_bounds[1] & distance_map < distance_bounds[2]

  if (save_qc_masks) {
    qc_file <- file.path(dirname(metric_file_fn(subject, time)), sprintf("%s_%s_valid_mask.nii.gz", subject, time))
    valid_mask_image <- asNifti(array(as.integer(valid_mask), dim = dim(valid_mask)), reference = metric_img)
    writeNifti(valid_mask_image, qc_file)
    message("  saved QC mask: ", qc_file)
  }

  dist_values   <- distance_map[valid_mask]
  scalar_values <- metric_map[valid_mask]
  region        <- ifelse(tumor_seg[valid_mask] == 2, "t2h", "nawm")

  data.table(
    participant  = subject, time = time, region = region,
    distance_raw = dist_values, scalar_raw = scalar_values
  )
}

# =============================== GAM HELPERS ====================================

# Drops patients missing either timepoint, or with too little data overall (<30 voxels).
prep_voxel_data <- function(voxel_data) {
  if (is.null(voxel_data) || nrow(voxel_data) < 30L) return(NULL)
  if (uniqueN(voxel_data$time) < 2L) return(NULL)
  voxel_data
}

# Fits one smooth curve per timepoint jointly: scalar_raw ~ time +
# smooth(distance_raw), with a separate smooth for "Pre" and "Post".
fit_gam_safe <- function(voxel_data, k_value) {
  tryCatch(
    gam(scalar_raw ~ time + s(distance_raw, by = time, bs = "cr", k = k_value),
        data = voxel_data, method = "REML"),
    error = function(e) NULL
  )
}

# Numerically estimates the slope of the fitted curve for one timepoint
# ("Pre" or "Post"), at each point on `distance_grid`.
get_slope <- function(gam_model, distance_grid, time_level) {
  tryCatch({
    eps <- diff(range(distance_grid)) / 1e5
    new_lo <- data.frame(distance_raw = distance_grid, time = factor(time_level, levels = c("Pre", "Post")))
    new_hi <- data.frame(distance_raw = distance_grid + eps, time = factor(time_level, levels = c("Pre", "Post")))
    (predict(gam_model, newdata = new_hi) - predict(gam_model, newdata = new_lo)) / eps
  }, error = function(e) NULL)
}

# Finds every location where the slope crosses zero in the requested
# direction, at positive distances only (see `landmark_type` above).
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

# Pairs up the Pre and Post landmarks: takes the first (nearest-tumor) Pre
# landmark as the anchor, then the Post landmark closest to it. Flags the
# pair if they're farther apart than `peak_match_window` -- likely a
# mismatched correspondence (e.g. a different curve feature entirely)
# rather than genuine Pre-to-Post change.
match_landmarks <- function(slope_pre, slope_post, distance_grid, max_shift = peak_match_window) {
  landmarks_pre  <- find_all_peaks(slope_pre,  distance_grid)
  landmarks_post <- find_all_peaks(slope_post, distance_grid)

  if (!length(landmarks_pre))
    return(list(loc_pre = NA_real_, loc_post = NA_real_,
                n_pre = 0L, n_post = length(landmarks_post),
                shift = NA_real_, shift_flag = TRUE))

  loc_pre <- landmarks_pre[1L]

  if (!length(landmarks_post))
    return(list(loc_pre = loc_pre, loc_post = NA_real_,
                n_pre = length(landmarks_pre), n_post = 0L,
                shift = NA_real_, shift_flag = TRUE))

  loc_post <- landmarks_post[which.min(abs(landmarks_post - loc_pre))]
  shift    <- abs(loc_post - loc_pre)

  list(loc_pre = loc_pre, loc_post = loc_post,
       n_pre = length(landmarks_pre), n_post = length(landmarks_post),
       shift = round(shift, 4), shift_flag = shift > max_shift)
}

# During bootstrap, searches for the Post landmark near a given anchor
# (rather than blindly taking the nearest-tumor one) -- keeps resampled
# Post landmarks correspondence-consistent with the Pre anchor.
find_landmark_near <- function(slope, distance_grid, anchor = NA_real_, window = peak_match_window, type = landmark_type) {
  if (is.null(slope) || length(slope) < 2L) return(NA_real_)
  n <- length(slope)
  near_anchor <- if (!is.na(anchor)) distance_grid >= (anchor - window) & distance_grid <= (anchor + window) else rep(TRUE, n)
  if (sum(near_anchor) < 2L) near_anchor <- rep(TRUE, n)

  crossing_idx <- switch(type,
    peak   = which(slope[-n] > 0 & slope[-1L] <= 0 & near_anchor[-n]),
    trough = which(slope[-n] < 0 & slope[-1L] >= 0 & near_anchor[-n])
  )
  crossing_idx <- crossing_idx[distance_grid[crossing_idx] > 0]
  if (!length(crossing_idx)) return(NA_real_)
  i <- crossing_idx[1L]
  distance_grid[i] - slope[i] * (distance_grid[i + 1L] - distance_grid[i]) / (slope[i + 1L] - slope[i])
}

# Average slope of the curve from the tumor up to `landmark_loc`.
mean_slope_to_landmark <- function(slope, distance_grid, landmark_loc) {
  if (is.na(landmark_loc) || is.null(slope)) return(NA_real_)
  idx <- distance_grid <= landmark_loc
  if (sum(idx) < 2L) return(NA_real_)
  mean(slope[idx], na.rm = TRUE)
}

# Runs mgcv's built-in check of whether k was large enough, combined
# conservatively across both the Pre and Post smooth terms (worst k-index,
# largest p-value of the two).
check_k_adequacy <- function(gam_model) {
  tryCatch({
    kc <- k.check(gam_model)
    list(k_index = min(kc[, "k-index"], na.rm = TRUE),
         k_pval  = max(kc[, "p-value"], na.rm = TRUE))
  }, error = function(e) list(k_index = NA_real_, k_pval = NA_real_))
}

# Chooses an initial k from sample size (pooled across both timepoints),
# fits the joint curve, and refits once with a larger k if the fit looks inadequate
select_k_and_fit <- function(voxel_data) {
  n_voxels <- nrow(voxel_data)
  k_value <- min(k_max_default, max(k_min, floor(n_voxels / voxels_per_k)))

  gam_model <- fit_gam_safe(voxel_data, k_value)
  if (is.null(gam_model)) return(list(mod = NULL, k = NA_integer_, k_reason = "fit_failed"))

  k_check <- check_k_adequacy(gam_model)
  edf <- sum(gam_model$edf)
  edf_k_ratio <- edf / k_value
  k_reason <- "fixed"

  if (!is.na(k_check$k_index) && k_check$k_index < k_index_adequate_threshold &&
      !is.na(edf_k_ratio) && edf_k_ratio > edf_k_ratio_flag_threshold &&
      k_value < k_max_fallback) {

    k_value_fallback <- k_max_fallback
    gam_model_fallback <- fit_gam_safe(voxel_data, k_value_fallback)

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

  list(mod = gam_model, k = k_value, k_index = k_check$k_index, k_pval = k_check$k_pval,
       edf = edf, k_reason = k_reason)
}

# Resamples this patient's voxels (pooled across both timepoints, so each
# draw keeps its original Pre/Post label) with replacement `n_boot` times,
# refits the joint curve, and re-locates both landmarks each time --
# Post anchored near the Pre landmark from that same draw -- building up
# paired distributions of plausible Pre and Post landmark locations.
bootstrap_patient <- function(voxel_data, k_value, distance_grid, n_boot = n_bootstrap) {
  n_voxels <- nrow(voxel_data)
  sample_n <- if (!is.null(bootstrap_voxel_cap) && bootstrap_voxel_cap < n_voxels) bootstrap_voxel_cap else n_voxels
  loc_pre_boot  <- rep(NA_real_, n_boot)
  loc_post_boot <- rep(NA_real_, n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample(n_voxels, sample_n, replace = TRUE)
    gam_model <- fit_gam_safe(voxel_data[idx], k_value)
    if (is.null(gam_model)) next
    slope_pre  <- get_slope(gam_model, distance_grid, "Pre")
    slope_post <- get_slope(gam_model, distance_grid, "Post")
    landmarks_pre <- find_all_peaks(slope_pre, distance_grid)
    loc_pre_boot[b]  <- if (length(landmarks_pre)) landmarks_pre[1L] else NA_real_
    loc_post_boot[b] <- find_landmark_near(slope_post, distance_grid, anchor = loc_pre_boot[b])
  }
  list(loc_pre = loc_pre_boot, loc_post = loc_post_boot)
}

# Summarizes a bootstrap distribution into a median and 95% CI.
summarize_bootstrap <- function(boot_values) {
  ci <- quantile(boot_values, c(0.025, 0.975), na.rm = TRUE)
  list(median = round(median(boot_values, na.rm = TRUE), 4),
       ci_lo = round(ci[1L], 4), ci_hi = round(ci[2L], 4))
}

# Full per-patient pipeline: combine both timepoints' voxel data, select k,
# fit the joint curve, match Pre/Post landmarks, bootstrap both, and
# package a two-row (per-timepoint) table plus a one-row (delta) table.
compute_patient <- function(patient_id, metric_table, measure, map_label) {

  voxel_data <- prep_voxel_data(metric_table[.(patient_id)])
  if (is.null(voxel_data))
    return(list(
      per_timepoint = data.table(participant = patient_id, metric = measure, map = map_label,
                                  timepoint = NA_character_, selected_k = NA_integer_,
                                  k_reason = "insufficient_data", has_peak = FALSE),
      delta = data.table()))

  distance_grid <- seq(min(voxel_data$distance_raw), max(voxel_data$distance_raw), length.out = n_grid_points)

  fit <- select_k_and_fit(voxel_data)
  if (is.null(fit$mod))
    return(list(
      per_timepoint = data.table(participant = patient_id, metric = measure, map = map_label,
                                  timepoint = NA_character_, selected_k = NA_integer_,
                                  k_reason = fit$k_reason, has_peak = FALSE),
      delta = data.table()))

  gam_model  <- fit$mod
  dev_expl   <- round(summary(gam_model)$dev.expl, 6)
  slope_pre  <- get_slope(gam_model, distance_grid, "Pre")
  slope_post <- get_slope(gam_model, distance_grid, "Post")

  matched <- match_landmarks(slope_pre, slope_post, distance_grid)
  boot    <- bootstrap_patient(voxel_data, fit$k, distance_grid)
  boot_pre  <- summarize_bootstrap(boot$loc_pre)
  boot_post <- summarize_bootstrap(boot$loc_post)

  per_timepoint <- rbindlist(lapply(c("Pre", "Post"), function(tp) {
    loc   <- if (tp == "Pre") matched$loc_pre else matched$loc_post
    slope <- if (tp == "Pre") slope_pre       else slope_post
    boot_summary <- if (tp == "Pre") boot_pre else boot_post
    data.table(
      participant         = patient_id, metric = measure, map = map_label, timepoint = tp,
      selected_k          = fit$k, k_reason = fit$k_reason,
      deviance_explained  = dev_expl,
      edf                 = round(fit$edf, 2),
      edf_k_ratio         = round(fit$edf / fit$k, 3),
      k_index             = round(fit$k_index, 3),
      k_pval              = round(fit$k_pval, 4),
      landmark            = round(loc, 4),
      has_peak            = !is.na(loc),
      n_landmarks_pre     = matched$n_pre,
      n_landmarks_post    = matched$n_post,
      landmark_shift      = matched$shift,
      landmark_shift_flag = matched$shift_flag,
      boot_median         = boot_summary$median,
      boot_ci_lo          = boot_summary$ci_lo,
      boot_ci_hi          = boot_summary$ci_hi,
      mean_slope_to_landmark = round(mean_slope_to_landmark(slope, distance_grid, loc), 6))
  }))

  delta <- data.table(
    participant           = patient_id, metric = measure, map = map_label,
    landmark_pre           = round(matched$loc_pre,  4),
    landmark_post          = round(matched$loc_post, 4),
    landmark_diff          = round(matched$loc_post - matched$loc_pre, 4),
    landmark_diff_norm     = round((matched$loc_post - matched$loc_pre) / max(voxel_data$distance_raw), 4),
    landmark_shift         = matched$shift,
    landmark_shift_flag    = matched$shift_flag,
    landmark_diff_boot_median = round(boot_post$median - boot_pre$median, 4),
    rate_to_landmark_pre   = round(mean_slope_to_landmark(slope_pre,  distance_grid, matched$loc_pre),  6),
    rate_to_landmark_post  = round(mean_slope_to_landmark(slope_post, distance_grid, matched$loc_post), 6),
    rate_to_landmark_diff  = round(mean_slope_to_landmark(slope_post, distance_grid, matched$loc_post) -
                                      mean_slope_to_landmark(slope_pre, distance_grid, matched$loc_pre), 6))

  list(per_timepoint = per_timepoint, delta = delta)
}

# ================================ MAIN LOOP ==================================
# For each distance-map variant and each metric: extract voxel data (both
# timepoints) for every participant directly from the NIfTI files, fit the
# GAM pre/post comparison in parallel, print summary tables, and save two
# CSVs of results.

for (map_label in map_labels) {
  for (measure in scalar_measures) {

    message(sprintf("=== metric: %s | distance-map variant: %s ===", measure, map_label))

    message("extracting voxel data from NIfTI files (Pre and Post)...")
    voxel_tables <- vector("list", length(participant_ids) * 2L)
    idx <- 1L
    for (i in seq_along(participant_ids)) {
      subject <- participant_ids[i]
      for (time in c("Pre", "Post")) {
        voxel_tables[[idx]] <- tryCatch(
          extract_voxel_data(subject, time, measure, map_label),
          error = function(e) {
            message(sprintf("  %s (%s) failed during extraction: %s", subject, time, e$message))
            NULL
          })
        idx <- idx + 1L
      }
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

    results_ok        <- Filter(Negate(is.null), results_all)
    results_per_tp    <- rbindlist(lapply(results_ok, `[[`, "per_timepoint"), fill = TRUE)
    results_delta     <- rbindlist(lapply(results_ok, `[[`, "delta"), fill = TRUE)
    rm(metric_table, results_all, results_ok); gc()

    message(sprintf("per-timepoint rows: %d | delta rows: %d", nrow(results_per_tp), nrow(results_delta)))

    # Curve fit quality, by timepoint.
    results_per_tp[!is.na(selected_k), .(
      median_dev = round(median(deviance_explained, na.rm = TRUE), 3),
      median_edf = round(median(edf, na.rm = TRUE), 2)
    ), by = timepoint] |> print()

    # Landmark detection rate and location, by timepoint.
    results_per_tp[, .(
      prop_has_peak   = round(mean(has_peak, na.rm = TRUE), 3),
      median_landmark = round(median(landmark, na.rm = TRUE), 3),
      median_boot_med = round(median(boot_median, na.rm = TRUE), 3),
      median_ci_width = round(median(boot_ci_hi - boot_ci_lo, na.rm = TRUE), 3)
    ), by = timepoint] |> print()

    # Pre-to-Post change summary.
    results_delta[, .(
      n                 = .N,
      n_shift_flagged   = sum(landmark_shift_flag, na.rm = TRUE),
      median_shift      = round(median(landmark_shift, na.rm = TRUE), 2),
      median_landmark_diff = round(median(landmark_diff, na.rm = TRUE), 4),
      median_rate_diff  = round(median(rate_to_landmark_diff, na.rm = TRUE), 6)
    )] |> print()

    flagged <- results_delta[landmark_shift_flag == TRUE]
    if (nrow(flagged)) {
      message(sprintf("%d participant(s) flagged for a large/likely-mismatched Pre-Post shift:", nrow(flagged)))
      flagged[, .(participant, landmark_pre, landmark_post, landmark_shift)] |> print()
    } else {
      message("no flagged landmark correspondences")
    }

    out_base <- file.path(output_dir, sprintf("%s_%s_GAM_prepost", measure, map_label))
    fwrite(results_per_tp, paste0(out_base, "_per_timepoint.csv"))
    fwrite(results_delta,  paste0(out_base, "_delta.csv"))
    message("saved: ", out_base, "_per_timepoint.csv, ", out_base, "_delta.csv")

    rm(results_per_tp, results_delta); gc()
  }
}
