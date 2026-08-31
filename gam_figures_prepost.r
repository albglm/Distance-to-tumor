################################################################################
# GAM LANDMARK VISUALIZATION: PRE- VS POST-TREATMENT COMPARISON
################################################################################
# WHAT IT DOES
#   Four figures for ONE metric and ONE distance-map variant at a time (set
#   both at the top of USER SETTINGS -- change and re-run for another):
#     1. PATIENT GRID: one small panel per patient, with their Pre and Post
#        fitted curves overlaid, each with its own landmark point and
#        bootstrap 95% CI. Patients flagged for a likely mismatched
#        Pre/Post correspondence are marked distinctly.
#     2. POPULATION OVERVIEW: every patient's Pre and Post curves faintly in
#        the background, with the cohort median curve and IQR ribbon for
#        each timepoint on top -- shows whether there's a systematic
#        before/after shift at the group level.
#     3. PAIRED SHIFT PLOT: each patient's Pre and Post landmark shown as
#        two connected points -- the individual-level view of who moved,
#        how far, and in which direction (a boxplot alone hides this).
#     4. COHORT BOXPLOTS: the Pre-to-Post landmark shift and gradient shift
#        across the full cohort, split by an optional subgroup label.
#
#   This does NOT re-run the GAM fitting -- it reads the results CSVs produced by
#   gam_landmark_analysis_prepost.R (for the saved k, landmarks, CIs, and
#   shift values) and re-extracts + refits each patient's joint Pre/Post
#   curve from the same NIfTI files purely for display.
#
# DEPENDENCIES
#   install.packages(c("RNifti", "data.table", "mgcv", "ggplot2"))
#
# REQUIRED INPUT
#   - The two results CSVs already produced by gam_landmark_analysis_prepost.R
#     for this (measure, map_label) combination (per-timepoint and delta).
#   - The same NIfTI inputs that pipeline used to produce them, at BOTH
#     timepoints (metric map, distance map, tumor segmentation, white
#     matter mask). Pre and Post must be co-registered.
#
# HOW TO RUN
#   1. Edit USER SETTINGS below (measure, map_label, paths, patient lists).
#   2. Run: Rscript gam_landmark_prepost_plots.R  (or Source in RStudio).
################################################################################

library(RNifti)
library(data.table)
library(mgcv)
library(ggplot2)

# ============================== USER SETTINGS ================================

## --- which metric and distance-map variant to plot ---------------------------
measure   <- "R2s"
map_label <- "iso"

## --- paths ---------
base_dir_pre  <- "/data/project/TP1"
base_dir_post <- "/data/project/TP2"
output_dir    <- file.path(base_dir_post, "res")     # where gam_landmark_analysis_prepost.R saved its results CSVs
figures_dir   <- file.path(base_dir_post, "figures")  # where this script saves its plots
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

participant_list <- fread("/data/project/part.csv", sep = "\t")
participant_ids  <- participant_list[["NP"]]   # full cohort -- used for the population overview and boxplots

# Patients shown in the per-patient grid (Figure 1). Defaults to the full
# cohort above; edit to a subset for a curated main-text figure. Order
# follows the vector given here.
patient_ids_for_grid <- participant_ids
# patient_ids_for_grid <- c("sub-P003", "sub-P012", "sub-P027", "sub-P041")

excluded_participants <- character(0)  # participant IDs to drop from all figures; character(0) = exclude no one
# excluded_participants <- c("sub-P037", "sub-P031", "sub-P015", "sub-P009")

## --- optional subgroup labels, for the boxplots (Figure 4) -------------------
split_by_group <- TRUE     # TRUE: split by `group_column` if that column exists; FALSE: always one box per feature
group_column   <- "group"  # name of the grouping column in participant_list

## --- must match the settings used to PRODUCE the results CSVs ----------------
n_grid_points             <- 200
scalar_clip_percentiles   <- c(0, 99)
distance_clip_percentiles <- c(1, 95)
excluded_region_mask_fn   <- NULL
# excluded_region_mask_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/GRE", subject, "reg/t1_first_seg_all_fast_firstseg_bin.nii.gz")

# ============================ INPUT FILE PATHS =================================
# Identical to gam_landmark_analysis_prepost.R -- edit both if your layout changes.
# If the distance map/segmentation/white matter mask are only
# computed once, at baseline, symlink the expected Post-timepoint path to
# the baseline file.

time_base_dir <- function(time) if (time == "Pre") base_dir_pre else base_dir_post

metric_file_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/GRE", subject, paste0(measure, ".nii"))
distance_map_file_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/DWI", subject,
  sprintf("distancemap_core_%s_dwi.nii.gz", map_label))
tumor_segmentation_file_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/segm", subject, "segmentation_brats.nii.gz")
white_matter_mask_file_fn <- function(subject, time) file.path(time_base_dir(time), "derivatives/DWI", subject, "reg/t1_fast_pve_2_dwi_thr.nii.gz")

# =============================== VOXEL EXTRACTION ===============================
# Identical to gam_landmark_analysis_prepost.R's extract_voxel_data().
extract_voxel_data <- function(subject, time, measure, map_label) {
  metric_map   <- as.array(readNifti(metric_file_fn(subject, time)))
  distance_map <- as.array(readNifti(distance_map_file_fn(subject, time)))
  tumor_seg    <- as.array(readNifti(tumor_segmentation_file_fn(subject, time)))
  white_matter_mask <- as.array(readNifti(white_matter_mask_file_fn(subject, time))) > 0

  tumor_exclusion_mask <- !(tumor_seg %in% c(1, 3))  # 1 = necrosis; 2 = T2H; 3 = contrast enhancing tumor
  t2h_mask <- tumor_seg == 2
  allowed_region <- white_matter_mask | t2h_mask

  valid <- is.finite(metric_map) & is.finite(distance_map) & (distance_map > 0) &
    tumor_exclusion_mask & allowed_region

  if (!is.null(excluded_region_mask_fn)) {
    excluded_region_mask <- as.array(readNifti(excluded_region_mask_fn(subject, time))) > 0
    valid <- valid & !excluded_region_mask
  }

  scalar_bounds   <- quantile(metric_map[valid], scalar_clip_percentiles / 100, na.rm = TRUE)
  distance_bounds <- quantile(distance_map[valid], distance_clip_percentiles / 100, na.rm = TRUE)

  valid_mask <- valid &
    metric_map > scalar_bounds[1] & metric_map < scalar_bounds[2] &
    distance_map > distance_bounds[1] & distance_map < distance_bounds[2]

  data.table(participant = subject, time = time,
             distance_raw = distance_map[valid_mask], scalar_raw = metric_map[valid_mask])
}

# Loads and combines a patient's Pre and Post voxel data, then refits the
# joint curve (using the k already selected by the pipeline) and returns
# predicted values for both timepoints over each timepoint's own 1st-99th
# percentile distance range.
refit_patient_curves <- function(subject, k_value) {
  voxel_data <- rbindlist(lapply(c("Pre", "Post"), function(tp)
    tryCatch(extract_voxel_data(subject, tp, measure, map_label), error = function(e) NULL)))
  if (is.null(voxel_data) || nrow(voxel_data) < 30L || uniqueN(voxel_data$time) < 2L) return(NULL)

  gam_model <- tryCatch(
    gam(scalar_raw ~ time + s(distance_raw, by = time, bs = "cr", k = k_value), data = voxel_data, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(gam_model)) return(NULL)

  rbindlist(lapply(c("Pre", "Post"), function(tp) {
    d_range <- quantile(voxel_data[time == tp]$distance_raw, c(0.01, 0.99), na.rm = TRUE)
    grid <- seq(d_range[1], d_range[2], length.out = n_grid_points)
    data.table(timepoint = tp, distance_raw = grid,
               fitted = predict(gam_model, newdata = data.frame(distance_raw = grid, time = factor(tp, levels = c("Pre", "Post")))))
  }))
}

# Reads the two results CSVs that gam_landmark_analysis_prepost.R produced
# for this (measure, map_label) combination.
load_per_timepoint_results <- function() {
  f <- file.path(output_dir, sprintf("%s_%s_GAM_prepost_per_timepoint.csv", measure, map_label))
  if (!file.exists(f)) stop("Per-timepoint results file not found: ", f, " -- run gam_landmark_analysis_prepost.R first.")
  dt <- fread(f)
  setkey(dt, participant)
  dt
}

load_delta_results <- function() {
  f <- file.path(output_dir, sprintf("%s_%s_GAM_prepost_delta.csv", measure, map_label))
  if (!file.exists(f)) stop("Delta results file not found: ", f, " -- run gam_landmark_analysis_prepost.R first.")
  dt <- fread(f)
  setkey(dt, participant)
  dt
}

# ============================== LABELS & STYLE ================================

theme_pub <- function(base_size = 16, strip_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title   = element_text(size = base_size + 3, face = "bold"),
      axis.title   = element_text(size = base_size),
      axis.text    = element_text(size = base_size - 3),
      strip.text   = element_text(size = strip_size, face = "bold"),
      legend.text  = element_text(size = base_size - 2),
      legend.title = element_text(size = base_size - 1)
    )
}

# Pre/Post colors (Okabe-Ito, colorblind-safe), distinguished by linetype
# too, so they hold up in grayscale.
PRE_COLOR   <- "#0072B2"  # Okabe-Ito blue
POST_COLOR  <- "#D55E00"  # Okabe-Ito vermillion
PRE_LINETYPE  <- "solid"
POST_LINETYPE <- "dashed"

# =============================================================================
# FIGURE 1: per-patient grid -- Pre and Post curves, landmarks, and CIs
# =============================================================================

build_figure_grid <- function() {
  per_tp <- load_per_timepoint_results()
  ids <- setdiff(patient_ids_for_grid, excluded_participants)

  curve_list <- list()
  landmark_list <- list()

  for (pid in ids) {
    rows <- per_tp[.(pid)]
    if (nrow(rows) < 2L || any(is.na(rows$selected_k))) next

    curve <- refit_patient_curves(pid, rows$selected_k[1L])
    if (is.null(curve)) next

    curve_list[[pid]] <- curve[, participant := pid]

    landmark_rows <- rows[, .(participant = pid, timepoint, landmark, boot_ci_lo, boot_ci_hi, landmark_shift_flag)]
    landmark_rows[, fitted_at_landmark := {
      sub <- curve[timepoint == .BY$timepoint]
      if (nrow(sub) == 0L || is.na(landmark)) NA_real_
      else approx(sub$distance_raw, sub$fitted, xout = landmark, rule = 2)$y
    }, by = timepoint]
    landmark_list[[pid]] <- landmark_rows
  }

  if (!length(curve_list)) return(NULL)
  curves_dt    <- rbindlist(curve_list)
  landmarks_dt <- rbindlist(landmark_list)

  panel_order <- ids[ids %in% unique(curves_dt$participant)]
  curves_dt[, facet_label    := factor(participant, levels = panel_order)]
  landmarks_dt[, facet_label := factor(participant, levels = panel_order)]

  n_col <- min(4L, length(panel_order))

  ggplot(curves_dt, aes(x = distance_raw, y = fitted, color = timepoint, linetype = timepoint)) +
    geom_rect(data = landmarks_dt, aes(xmin = boot_ci_lo, xmax = boot_ci_hi, ymin = -Inf, ymax = Inf, fill = timepoint),
              inherit.aes = FALSE, alpha = 0.12) +
    geom_line(linewidth = 0.8) +
    geom_point(data = landmarks_dt, aes(x = landmark, y = fitted_at_landmark, color = timepoint,
                                          shape = landmark_shift_flag),
               inherit.aes = FALSE, size = 2.6) +
    scale_color_manual(values = c(Pre = PRE_COLOR, Post = POST_COLOR)) +
    scale_fill_manual(values  = c(Pre = PRE_COLOR, Post = POST_COLOR)) +
    scale_linetype_manual(values = c(Pre = PRE_LINETYPE, Post = POST_LINETYPE)) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c(`FALSE` = "matched", `TRUE` = "flagged"), name = "correspondence") +
    facet_wrap(~ facet_label, scales = "free", ncol = n_col) +
    labs(title = paste(measure, "-", map_label, "- Pre vs Post"),
         x = "Geodesic distance from CET boundary", y = measure, color = "timepoint", linetype = "timepoint") +
    theme_pub(base_size = 13, strip_size = 11)
}

# =============================================================================
# FIGURE 2: population overview -- Pre vs Post, faint individuals + IQR
# =============================================================================

build_figure_population <- function() {
  per_tp <- load_per_timepoint_results()
  ids <- setdiff(participant_ids, excluded_participants)

  curve_list <- list()

  for (pid in ids) {
    rows <- per_tp[.(pid)]
    if (nrow(rows) < 2L || any(is.na(rows$selected_k))) next
    curve <- refit_patient_curves(pid, rows$selected_k[1L])
    if (is.null(curve)) next
    curve_list[[pid]] <- curve[, participant := pid]
  }

  if (!length(curve_list)) return(NULL)
  all_curves <- rbindlist(curve_list)

  # normalize each patient's distance to 0-1 within its own timepoint range
  all_curves[, rel_distance := (distance_raw - min(distance_raw)) / (max(distance_raw) - min(distance_raw)),
             by = .(participant, timepoint)]

  common_grid <- seq(0, 1, length.out = n_grid_points)
  interp_list <- list()
  for (pid in unique(all_curves$participant)) {
    for (tp in c("Pre", "Post")) {
      sub <- all_curves[participant == pid & timepoint == tp]
      if (nrow(sub) < 2L) next
      interp_list[[paste(pid, tp)]] <- data.table(
        participant = pid, timepoint = tp, rel_distance = common_grid,
        fitted = approx(sub$rel_distance, sub$fitted, xout = common_grid, rule = 2)$y)
    }
  }
  interp_dt <- rbindlist(interp_list)

  summary_dt <- interp_dt[, .(
    median_fitted = median(fitted, na.rm = TRUE),
    q25 = quantile(fitted, 0.25, na.rm = TRUE),
    q75 = quantile(fitted, 0.75, na.rm = TRUE)
  ), by = .(timepoint, rel_distance)]

  ggplot() +
    geom_line(data = interp_dt, aes(x = rel_distance, y = fitted, group = interaction(participant, timepoint), color = timepoint),
              alpha = 0.25, linewidth = 0.35) +
    geom_ribbon(data = summary_dt, aes(x = rel_distance, ymin = q25, ymax = q75, fill = timepoint), alpha = 0.2) +
    geom_line(data = summary_dt, aes(x = rel_distance, y = median_fitted, color = timepoint, linetype = timepoint), linewidth = 1.3) +
    scale_color_manual(values = c(Pre = PRE_COLOR, Post = POST_COLOR)) +
    scale_fill_manual(values  = c(Pre = PRE_COLOR, Post = POST_COLOR)) +
    scale_linetype_manual(values = c(Pre = PRE_LINETYPE, Post = POST_LINETYPE)) +
    labs(title = paste(measure, "-", map_label, "- Pre vs Post"),
         x = "Normalized distance from CET boundary (0-1)", y = measure, color = "timepoint", fill = "timepoint", linetype = "timepoint") +
    theme_pub()
}

# =============================================================================
# FIGURE 3: paired shift plot -- each patient's Pre -> Post landmark
# =============================================================================

build_figure_paired_shift <- function() {
  per_tp <- load_per_timepoint_results()
  dt <- per_tp[!is.na(landmark) & !(participant %in% excluded_participants),
               .(participant, timepoint, landmark, landmark_shift_flag)]
  if (!nrow(dt)) return(NULL)

  dt[, timepoint := factor(timepoint, levels = c("Pre", "Post"))]
  dt[, flag_label := ifelse(landmark_shift_flag, "flagged", "matched")]

  ggplot(dt, aes(x = timepoint, y = landmark, group = participant, color = flag_label)) +
    geom_line(alpha = 0.6, linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c(matched = "grey40", flagged = "#D55E00"), name = "correspondence") +
    labs(title = paste(measure, "-", map_label, "- paired landmark shift"),
         x = NULL, y = paste(measure, "landmark distance")) +
    theme_pub()
}

# =============================================================================
# FIGURE 4: cohort boxplots -- landmark shift & gradient shift, optional subgroup
# =============================================================================

build_figure_boxplots <- function() {
  delta <- load_delta_results()
  dt <- delta[!(participant %in% excluded_participants),
              .(participant, landmark_shift = landmark_diff, gradient_shift = rate_to_landmark_diff)]

  if (split_by_group && group_column %in% names(participant_list)) {
    subgroup_dt <- participant_list[, .(participant = NP, group = get(group_column))]
    dt <- merge(dt, subgroup_dt, by = "participant", all.x = TRUE)
  } else {
    dt[, group := "All patients"]
  }

  long_dt <- melt(dt, id.vars = c("participant", "group"),
                   measure.vars = c("landmark_shift", "gradient_shift"),
                   variable.name = "feature", value.name = "value")
  long_dt[, feature := factor(feature, levels = c("landmark_shift", "gradient_shift"),
                                labels = c("Landmark shift (Post - Pre)", "Gradient shift (Post - Pre)"))]

  base_group_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#999999")
  n_groups <- uniqueN(long_dt$group)

  p <- ggplot(long_dt, aes(x = group, y = value, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.6, size = 1.8) +
    facet_wrap(~ feature, scales = "free_y")

  if (n_groups <= length(base_group_colors)) {
    p <- p + scale_fill_manual(values = base_group_colors[seq_len(n_groups)])
  }

  p +
    labs(title = paste(measure, "-", map_label, "- Pre-to-Post shift"), x = NULL, y = NULL) +
    theme_pub() +
    theme(legend.position = if (n_groups > 1) "right" else "none")
}

# ================================== RUN =======================================

message("metric: ", measure, " | map: ", map_label)

p_grid <- build_figure_grid()
if (!is.null(p_grid)) {
  n_patients <- length(setdiff(patient_ids_for_grid, excluded_participants))
  n_col <- min(4L, n_patients)
  n_rows <- ceiling(n_patients / n_col)
  out_grid <- file.path(figures_dir, sprintf("%s_%s_prepost_patient_grid.pdf", measure, map_label))
  pdf(out_grid, width = 6 * n_col, height = 4 * n_rows)
  print(p_grid)
  dev.off()
  message("saved: ", out_grid)
}

p_population <- build_figure_population()
if (!is.null(p_population)) {
  out_population <- file.path(figures_dir, sprintf("%s_%s_prepost_population_overview.pdf", measure, map_label))
  pdf(out_population, width = 10, height = 7)
  print(p_population)
  dev.off()
  message("saved: ", out_population)
}

p_paired <- build_figure_paired_shift()
if (!is.null(p_paired)) {
  out_paired <- file.path(figures_dir, sprintf("%s_%s_prepost_paired_shift.pdf", measure, map_label))
  pdf(out_paired, width = 6, height = 6)
  print(p_paired)
  dev.off()
  message("saved: ", out_paired)
}

p_boxplots <- build_figure_boxplots()
if (!is.null(p_boxplots)) {
  out_boxplots <- file.path(figures_dir, sprintf("%s_%s_prepost_boxplots.pdf", measure, map_label))
  pdf(out_boxplots, width = 8, height = 5)
  print(p_boxplots)
  dev.off()
  message("saved: ", out_boxplots)
}
