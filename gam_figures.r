################################################################################
# GAM LANDMARK VISUALIZATION
################################################################################
# WHAT IT DOES
#   Three figures for ONE metric and ONE distance-map variant at a time (set
#   both at the top of USER SETTINGS -- change and re-run for another):
#     1. PATIENT GRID: one small panel per patient (from a list you choose),
#        each showing that patient's fitted curve, their T2H/NAWM reference
#        lines, and their landmark with its bootstrap 95% CI.
#     2. POPULATION OVERVIEW: every patient's fitted curve faintly in the
#        background, with the cohort median curve and IQR ribbon on top,
#        plus cohort-average T2H/NAWM reference lines.
#     3. COHORT BOXPLOTS: landmark distance and mean gradient-to-landmark
#        across the full cohort, split by an optional subgroup label (e.g.
#        MGMT status) for comparing two groups.
#
#   This script does NOT re-run the GAM fitting pipeline -- it reads the
#   results CSV produced by gam_landmark_pipeline.R (for the saved k,
#   landmark, and bootstrap CI per patient) and re-extracts + refits each
#   patient's curve from the same NIfTI files purely for plotting, since the
#   pipeline doesn't store the full fitted curve or raw voxels in that CSV.
#
# DEPENDENCIES
#   install.packages(c("RNifti", "data.table", "mgcv", "ggplot2"))
#
# REQUIRED INPUT
#   - The GAM results CSV already produced by gam_landmark_pipeline.R for
#     this (measure, map_label) combination.
#   - The same NIfTI inputs gam_landmark_pipeline.R used to produce it
#     (metric map, distance map, tumor segmentation, white matter mask).
#
# HOW TO RUN
#   1. Edit USER SETTINGS below (measure, map_label, paths, patient lists).
#   2. Run: Rscript gam_landmark_plots.R  (or Source in RStudio).
################################################################################

library(RNifti)
library(data.table)
library(mgcv)
library(ggplot2)

# ============================== USER SETTINGS ================================

## --- which metric and distance-map variant to plot ---------------------------
# One at a time -- change these and re-run for a different metric/map.
measure   <- "R2s"
map_label <- "iso"

## --- paths (copy from your gam_landmark_pipeline.R settings) -----------------
base_dir   <- "/data/project/TP1"
output_dir <- file.path(base_dir, "res")       # where gam_landmark_pipeline.R saved its results CSV
figures_dir <- file.path(base_dir, "figures")  # where this script saves its plots
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

participant_list <- fread("/data/project/part.csv", sep = "\t")
participant_ids  <- participant_list[["NP"]]   # full cohort -- used for the population overview and boxplots

# Patients shown in the per-patient grid (Figure 1). Defaults to the full
# cohort above; edit to a subset or exclude for a curated main-text figure. Order
# follows the vector given here.
patient_ids_for_grid <- participant_ids
# patient_ids_for_grid <- c("sub-P003", "sub-P012", "sub-P027", "sub-P041")   # a hand-picked subset
# patient_ids_for_grid <- participant_list[some_column == TRUE, NP]           # or filtered from a column in participant_list
excluded_participants <- character(0)  # participant IDs to drop from all figures; character(0) = exclude no one
# excluded_participants <- c("sub-P037", "sub-P031", "sub-P015", "sub-P009")

## --- optional subgroup labels, for the boxplots (Figure 3) -------------------
# Add a column to participant_list above (e.g. MGMT status, tumor subtype,
# or any categorical label -- 2 or more levels are both fine) to enable
# splitting the boxplots by it.
split_by_group <- TRUE     # TRUE: split by `group_column` if that column exists; FALSE: always one box per feature
group_column   <- "group"  # name of the grouping column in participant_list

## --- must match the settings used to PRODUCE the results CSV -----------------
# Curves are refit here purely for plotting; only k, the landmark, and its
# CI are stored in the results CSV, not the fitted curve itself.
smoothing_penalty         <- 1.5
n_grid_points             <- 200
scalar_clip_percentiles   <- c(0, 99)
distance_clip_percentiles <- c(1, 95)
excluded_region_mask_fn   <- NULL
# excluded_region_mask_fn <- function(subject) file.path(base_dir, "derivatives/GRE", subject, "reg/t1_first_seg_all_fast_firstseg_bin.nii.gz")

# ============================ INPUT FILE PATHS =================================
# Identical to gam_landmark_pipeline.R -- edit both if your layout changes.

metric_file_fn <- function(subject) file.path(base_dir, "derivatives/GRE", subject, paste0(measure, ".nii"))
distance_map_file_fn <- function(subject) file.path(base_dir, "derivatives/DWI", subject,
  sprintf("distancemap_core_%s_dwi.nii.gz", map_label))
tumor_segmentation_file_fn <- function(subject) file.path(base_dir, "derivatives/segm", subject, "segmentation_brats.nii.gz")
white_matter_mask_file_fn <- function(subject) file.path(base_dir, "derivatives/DWI", subject, "reg/t1_fast_pve_2_dwi_thr.nii.gz")

# =============================== VOXEL EXTRACTION ===============================
# Identical to gam_landmark_pipeline.R's extract_voxel_data().
extract_voxel_data <- function(subject, measure, map_label) {
  metric_map   <- as.array(readNifti(metric_file_fn(subject)))
  distance_map <- as.array(readNifti(distance_map_file_fn(subject)))
  tumor_seg    <- as.array(readNifti(tumor_segmentation_file_fn(subject)))
  white_matter_mask <- as.array(readNifti(white_matter_mask_file_fn(subject))) > 0

  tumor_exclusion_mask <- !(tumor_seg %in% c(1, 3))  # 1 = necrosis; 2 = T2H; 3 = contrast enhancing tumor
  t2h_mask <- tumor_seg == 2
  allowed_region <- white_matter_mask | t2h_mask

  valid <- is.finite(metric_map) & is.finite(distance_map) & (distance_map > 0) &
    tumor_exclusion_mask & allowed_region

  if (!is.null(excluded_region_mask_fn)) {
    excluded_region_mask <- as.array(readNifti(excluded_region_mask_fn(subject))) > 0
    valid <- valid & !excluded_region_mask
  }

  scalar_bounds   <- quantile(metric_map[valid], scalar_clip_percentiles / 100, na.rm = TRUE)
  distance_bounds <- quantile(distance_map[valid], distance_clip_percentiles / 100, na.rm = TRUE)

  valid_mask <- valid &
    metric_map > scalar_bounds[1] & metric_map < scalar_bounds[2] &
    distance_map > distance_bounds[1] & distance_map < distance_bounds[2]

  dist_values   <- distance_map[valid_mask]
  scalar_values <- metric_map[valid_mask]
  region        <- ifelse(tumor_seg[valid_mask] == 2, "t2h", "nawm")

  data.table(participant = subject, region = region, distance_raw = dist_values, scalar_raw = scalar_values)
}

# Refits one patient's curve using their already-selected k (from the
# results CSV), over the 1st-99th percentile range of their distances.
refit_patient_curve <- function(voxel_data, k_value) {
  d_lo <- quantile(voxel_data$distance_raw, 0.01, na.rm = TRUE)
  d_hi <- quantile(voxel_data$distance_raw, 0.99, na.rm = TRUE)

  gam_model <- tryCatch(
    gam(scalar_raw ~ s(distance_raw, bs = "cr", k = k_value), data = voxel_data,
        gamma = smoothing_penalty, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(gam_model)) return(NULL)

  grid <- seq(d_lo, d_hi, length.out = n_grid_points)
  data.table(distance_raw = grid, fitted = predict(gam_model, newdata = data.frame(distance_raw = grid)))
}

# Reads the results CSV that gam_landmark_pipeline.R produced for this
# (measure, map_label) combination.
load_gam_results <- function() {
  f <- file.path(output_dir, sprintf("%s_%s_GAM.csv", measure, map_label))
  if (!file.exists(f)) stop("GAM results file not found: ", f, " -- run gam_landmark_pipeline.R first.")
  dt <- fread(f)
  setkey(dt, participant)
  dt
}

# ============================== LABELS & STYLE ================================

# One place to control font sizing across all figures.
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

# Reference-line style, shared across figures. Colors are from the
# Okabe-Ito palette (verified to stay distinguishable under protanopia/
# deuteranopia/tritanopia); linetypes differ too, so the T2H/NAWM
# distinction survives grayscale printing as well.
T2H_COLOR      <- "#D55E00"  # Okabe-Ito vermillion
NAWM_COLOR     <- "#009E73"  # Okabe-Ito bluish green
T2H_LINETYPE   <- "dashed"
NAWM_LINETYPE  <- "dotted"
REFLINE_WIDTH  <- 0.9

# =============================================================================
# FIGURE 1: per-patient grid -- fitted curve + reference lines + landmark + CI
# =============================================================================

build_figure_grid <- function() {
  gam_results <- load_gam_results()
  ids <- setdiff(patient_ids_for_grid, excluded_participants)

  curve_list <- list()
  landmark_list <- list()

  for (pid in ids) {
    row <- gam_results[.(pid)]
    if (nrow(row) == 0L || is.na(row$selected_k)) next

    voxel_data <- tryCatch(extract_voxel_data(pid, measure, map_label), error = function(e) NULL)
    if (is.null(voxel_data) || nrow(voxel_data) < 30L) next

    curve <- refit_patient_curve(voxel_data, row$selected_k)
    if (is.null(curve)) next

    fitted_at_landmark <- if (!is.na(row$fp_boot_median))
      approx(curve$distance_raw, curve$fitted, xout = row$fp_boot_median, rule = 2)$y else NA_real_

    curve_list[[pid]] <- curve[, participant := pid]
    landmark_list[[pid]] <- data.table(
      participant = pid,
      mean_t2h  = mean(voxel_data$scalar_raw[voxel_data$region == "t2h"],  na.rm = TRUE),
      mean_nawm = mean(voxel_data$scalar_raw[voxel_data$region == "nawm"], na.rm = TRUE),
      fp_boot_median = row$fp_boot_median, fp_boot_ci_lo = row$fp_boot_ci_lo, fp_boot_ci_hi = row$fp_boot_ci_hi,
      fitted_at_landmark = fitted_at_landmark)
  }

  if (!length(curve_list)) return(NULL)
  curves_dt    <- rbindlist(curve_list)
  landmarks_dt <- rbindlist(landmark_list)

  panel_order <- ids[ids %in% unique(curves_dt$participant)]  # preserves the input list's order
  curves_dt[, facet_label    := factor(participant, levels = panel_order)]
  landmarks_dt[, facet_label := factor(participant, levels = panel_order)]

  n_col <- min(4L, length(panel_order))

  ggplot(curves_dt, aes(x = distance_raw, y = fitted)) +
    geom_rect(data = landmarks_dt, aes(xmin = fp_boot_ci_lo, xmax = fp_boot_ci_hi, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "steelblue", alpha = 0.15) +
    geom_line(linewidth = 0.7, color = "black") +
    geom_hline(data = landmarks_dt, aes(yintercept = mean_t2h),  linetype = T2H_LINETYPE,  color = T2H_COLOR,  linewidth = REFLINE_WIDTH) +
    geom_hline(data = landmarks_dt, aes(yintercept = mean_nawm), linetype = NAWM_LINETYPE, color = NAWM_COLOR, linewidth = REFLINE_WIDTH) +
    geom_point(data = landmarks_dt, aes(x = fp_boot_median, y = fitted_at_landmark),
               inherit.aes = FALSE, color = "steelblue4", size = 2.4) +
    facet_wrap(~ facet_label, scales = "free", ncol = n_col) +
    labs(title = paste(measure, "-", map_label),
         x = "Geodesic distance from CET boundary", y = measure) +
    theme_pub(base_size = 13, strip_size = 11)
}

# =============================================================================
# FIGURE 2: population overview -- all patients faint, median + IQR on top
# =============================================================================

build_figure_population <- function() {
  gam_results <- load_gam_results()
  ids <- setdiff(participant_ids, excluded_participants)

  curve_list <- list()
  region_list <- list()

  for (pid in ids) {
    row <- gam_results[.(pid)]
    if (nrow(row) == 0L || is.na(row$selected_k)) next

    voxel_data <- tryCatch(extract_voxel_data(pid, measure, map_label), error = function(e) NULL)
    if (is.null(voxel_data) || nrow(voxel_data) < 30L) next

    curve <- refit_patient_curve(voxel_data, row$selected_k)
    if (is.null(curve)) next

    d_lo <- min(curve$distance_raw); d_hi <- max(curve$distance_raw)
    curve[, rel_distance := (distance_raw - d_lo) / (d_hi - d_lo)]
    curve_list[[pid]] <- curve[, participant := pid]

    region_list[[pid]] <- data.table(
      participant = pid,
      mean_t2h  = mean(voxel_data$scalar_raw[voxel_data$region == "t2h"],  na.rm = TRUE),
      mean_nawm = mean(voxel_data$scalar_raw[voxel_data$region == "nawm"], na.rm = TRUE))
  }

  if (!length(curve_list)) return(NULL)
  all_curves <- rbindlist(curve_list)
  region_dt  <- rbindlist(region_list)

  common_grid <- seq(0, 1, length.out = n_grid_points)
  interp_list <- list()
  for (pid in unique(all_curves$participant)) {
    sub <- all_curves[participant == pid]
    if (nrow(sub) < 2L) next
    interp_list[[pid]] <- data.table(participant = pid, rel_distance = common_grid,
                                      fitted = approx(sub$rel_distance, sub$fitted, xout = common_grid, rule = 2)$y)
  }
  interp_dt <- rbindlist(interp_list)

  summary_dt <- interp_dt[, .(
    median_fitted = median(fitted, na.rm = TRUE),
    q25 = quantile(fitted, 0.25, na.rm = TRUE),
    q75 = quantile(fitted, 0.75, na.rm = TRUE)
  ), by = rel_distance]

  cohort_t2h  <- mean(region_dt$mean_t2h,  na.rm = TRUE)
  cohort_nawm <- mean(region_dt$mean_nawm, na.rm = TRUE)

  ggplot() +
    geom_line(data = interp_dt, aes(x = rel_distance, y = fitted, group = participant),
              color = "grey70", alpha = 0.5, linewidth = 0.4) +
    geom_ribbon(data = summary_dt, aes(x = rel_distance, ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.25) +
    geom_line(data = summary_dt, aes(x = rel_distance, y = median_fitted), color = "steelblue4", linewidth = 1.3) +
    geom_hline(yintercept = cohort_t2h,  linetype = T2H_LINETYPE,  color = T2H_COLOR,  linewidth = REFLINE_WIDTH) +
    geom_hline(yintercept = cohort_nawm, linetype = NAWM_LINETYPE, color = NAWM_COLOR, linewidth = REFLINE_WIDTH) +
    labs(title = paste(measure, "-", map_label),
         x = "Normalized distance from CET boundary (0-1)", y = measure) +
    theme_pub()
}

# =============================================================================
# FIGURE 3: cohort boxplots -- landmark distance & gradient, optional subgroup
# =============================================================================

build_figure_boxplots <- function() {
  gam_results <- load_gam_results()
  dt <- gam_results[!is.na(selected_k) & !(participant %in% excluded_participants),
                     .(participant, landmark_distance = fp_boot_median, gradient = mean_deriv_to_peak)]

  if (split_by_group && group_column %in% names(participant_list)) {
    subgroup_dt <- participant_list[, .(participant = NP, group = get(group_column))]
    dt <- merge(dt, subgroup_dt, by = "participant", all.x = TRUE)
  } else {
    dt[, group := "All patients"]
  }

  long_dt <- melt(dt, id.vars = c("participant", "group"),
                   measure.vars = c("landmark_distance", "gradient"),
                   variable.name = "feature", value.name = "value")
  long_dt[, feature := factor(feature, levels = c("landmark_distance", "gradient"),
                                labels = c("Landmark distance", "Mean gradient to landmark"))]

  # Okabe-Ito-based palette, colorblind-safe up to 8 groups; beyond that,
  # falls back to ggplot's default discrete hue scale rather than
  # recycling colors onto multiple groups.
  base_group_colors <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#999999")
  n_groups <- uniqueN(long_dt$group)

  p <- ggplot(long_dt, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.6, size = 1.8) +
    facet_wrap(~ feature, scales = "free_y")

  if (n_groups <= length(base_group_colors)) {
    p <- p + scale_fill_manual(values = base_group_colors[seq_len(n_groups)])
  }

  p +
    labs(title = paste(measure, "-", map_label), x = NULL, y = NULL) +
    theme_pub() +
    theme(legend.position = if (uniqueN(long_dt$group) > 1) "right" else "none")
}

# ================================== RUN =======================================

message("metric: ", measure, " | map: ", map_label)

p_grid <- build_figure_grid()
if (!is.null(p_grid)) {
  n_patients <- length(setdiff(patient_ids_for_grid, excluded_participants))
  n_col <- min(4L, n_patients)
  n_rows <- ceiling(n_patients / n_col)
  out_grid <- file.path(figures_dir, sprintf("%s_%s_patient_grid.pdf", measure, map_label))
  pdf(out_grid, width = 6 * n_col, height = 4 * n_rows)
  print(p_grid)
  dev.off()
  message("saved: ", out_grid)
}

p_population <- build_figure_population()
if (!is.null(p_population)) {
  out_population <- file.path(figures_dir, sprintf("%s_%s_population_overview.pdf", measure, map_label))
  pdf(out_population, width = 10, height = 7)
  print(p_population)
  dev.off()
  message("saved: ", out_population)
}

p_boxplots <- build_figure_boxplots()
if (!is.null(p_boxplots)) {
  out_boxplots <- file.path(figures_dir, sprintf("%s_%s_boxplots.pdf", measure, map_label))
  pdf(out_boxplots, width = 8, height = 5)
  print(p_boxplots)
  dev.off()
  message("saved: ", out_boxplots)
}