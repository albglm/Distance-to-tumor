################################################################################
# Summary tables of the distance-to-tumor qMRI profiles
################################################################################
#
# Reads the *_GAM.csv files of gam_distance_profiles.R (one run) and writes two
# tables, per metric and distance map (median [IQR] across participants):
#
#   fit_summary.csv       voxels, deviance explained, edf, Moran's I of the
#                         residuals, k-index, landmark CI width, and how the
#                         landmark was selected (initial tolerance, widened,
#                         fallback, full distance range)
#   feature_summary.csv   landmark distance and local gradient
#
# HOW TO RUN
#   1. Set results_dir, run_tag, metrics and maps below.
#   2. Rscript summarize_profiles.R
#
# Required R package: data.table
################################################################################

library(data.table)

# ================================ SETTINGS ====================================

results_dir <- "/path/to/bids/derivatives/distance_profiles"   # output_dir of gam_distance_profiles.R
run_tag     <- "g1.5_tol1_blk5mm_t2h0"
metrics     <- c("R2s", "QSM", "ChiDia", "ChiPara", "adc", "fa")
maps        <- c("iso", "isoweighted", "aniso", "anisoweighted")

# ================================== LOAD ======================================

gam <- rbindlist(lapply(metrics, function(m) rbindlist(lapply(maps, function(mp) {
  f <- file.path(results_dir, sprintf("%s_%s_%s_GAM.csv", m, mp, run_tag))
  if (file.exists(f)) fread(f)[status == "ok"] else { message("missing: ", f); NULL }
}))), fill = TRUE)
if (!nrow(gam)) stop("no results found in ", results_dir)
gam[, `:=`(metric = factor(metric, metrics), map = factor(map, maps))]

med_iqr <- function(x, digits = 2) {
  q <- quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE)
  sprintf("%s [%s, %s]", signif(q[1], digits + 1), signif(q[2], digits + 1), signif(q[3], digits + 1))
}
pct <- function(x) sprintf("%.1f", 100 * mean(x))

# =============================== FIT SUMMARY ==================================

fit <- gam[, .(
  n                       = .N,
  voxels                  = med_iqr(n_voxels, 3),
  deviance_explained_pct  = med_iqr(100 * deviance_explained),
  edf                     = med_iqr(edf),
  morans_i                = med_iqr(moran_i_resid),
  k_index_min             = round(min(k_index, na.rm = TRUE), 2),
  landmark_ci_width       = med_iqr(fp_boot_ci_width_norm),          # fraction of the distance range
  initial_tolerance_pct   = pct(grepl("^initial_tol", landmark_selection)),
  widened_tolerance_pct   = pct(grepl("^widened", landmark_selection)),
  fallback_pct            = pct(grepl("^closest", landmark_selection)),
  full_range_search_pct   = pct(grepl("no_t2h_range", landmark_selection)),
  no_extremum             = sum(landmark_selection == "no_extremum")
), keyby = .(metric, map)]

# ============================= FEATURE SUMMARY ================================

features <- gam[, .(
  n                 = .N,
  landmark_distance = med_iqr(fcoalesce(as.numeric(fp_boot_median), first_peak)),   # bootstrap median, or the fit's landmark without bootstrap; index units (mm for iso)
  local_gradient    = med_iqr(mean_deriv_to_peak)   # metric units per distance unit
), keyby = .(metric, map)]

# ================================== SAVE ======================================

fwrite(fit,      file.path(results_dir, sprintf("fit_summary_%s.csv", run_tag)))
fwrite(features, file.path(results_dir, sprintf("feature_summary_%s.csv", run_tag)))
print(fit); print(features)
message("saved to ", results_dir)
