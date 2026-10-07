################################################################################
# Distance-to-tumor qMRI profiles: GAM fit, landmark and local gradient
################################################################################
#
# For every participant, qMRI metric and distance map, this script
#   1. collects the voxels of the sampling domain (T2H + NAWM) with their
#      geodesic distance from the contrast-enhancing tumor (CET) and metric value,
#   2. fits  metric ~ s(distance)  with a GAM,
#   3. finds the LANDMARK (first peak, or trough for ADC, whose fitted value is
#      close to the NAWM mean) and the LOCAL GRADIENT (mean slope up to it),
#   4. estimates landmark uncertainty with a spatial block bootstrap,
#   5. reports fit diagnostics (deviance explained, edf, k-index, AIC gain,
#      Moran's I of the residuals),
#   6. builds population curves (equal-weight mean of patient curves on a
#      normalized 0-1 distance grid, with a participant-bootstrap 95% CI).
#
# HOW TO RUN
#   1. Edit the two blocks below: "PATHS AND INPUT FILES" and "ANALYSIS".
#   2. Rscript gam_distance_profiles.R                    all included participants
#      Rscript gam_distance_profiles.R sub-P001 sub-P002  only these participants
#   The default settings reproduce the main analysis of the paper. The other
#   analyses (contralateral control, sensitivity analyses) are obtained by
#   changing the settings marked [VARIANT]; see README.md.
#
# OUTPUT (in output_dir), one set of files per metric x map x run:
#   <metric>_<map>_<run_tag>_GAM.csv       one row per participant
#   <metric>_<map>_<run_tag>_curves.csv    participant curves on the 0-1 grid
#   <metric>_<map>_<run_tag>_popcurve.csv  population curve with 95% CI
#
# Required R packages: RNifti, data.table, mgcv, future, furrr
################################################################################

library(RNifti)
library(data.table)
library(mgcv)
library(future)
library(furrr)

# ============================ PATHS AND INPUT FILES ===========================
# Input follows a BIDS-derivatives layout (see README.md):
#   <bids_root>/participants.tsv
#   <bids_root>/derivatives/qmri/sub-<id>/{anat,dwi}/          qMRI maps
#   <bids_root>/derivatives/masks/sub-<id>/{anat,dwi}/         masks
#   <bids_root>/derivatives/distancemaps/sub-<id>/{anat,dwi}/  distance maps

bids_root  <- "/path/to/bids"
output_dir <- file.path(bids_root, "derivatives", "distance_profiles")

# participants: IDs given on the command line, otherwise all rows of
# participants.tsv with include = 1 (columns participant_id, include)
participant_ids <- commandArgs(trailingOnly = TRUE)
if (!length(participant_ids))
  participant_ids <- fread(file.path(bids_root, "participants.tsv"), sep = "\t")[include == 1, participant_id]

# Each metric is sampled in its own space: GRE-derived metrics in GRE space
# (anat), ADC and FA in diffusion space (dwi). All masks and distance maps must
# exist in both spaces.
image_space <- function(metric) if (metric %in% c("adc", "fa")) "dwi" else "gre"

# file name of each qMRI metric (BIDS suffixes; desc-/param- entities)
metric_file <- c(
  R2s     = "space-gre_R2starmap",
  QSM     = "space-gre_Chimap",
  ChiDia  = "space-gre_desc-dia_Chimap",
  ChiPara = "space-gre_desc-para_Chimap",
  adc     = "space-dwi_model-tensor_param-adc_dwimap",
  fa      = "space-dwi_model-tensor_param-fa_dwimap")

input_files <- function(subject, metric, map) {
  space  <- image_space(metric)
  folder <- if (space == "dwi") "dwi" else "anat"
  deriv  <- function(pipeline, name)
    file.path(bids_root, "derivatives", pipeline, subject, folder, sprintf("%s_%s.nii.gz", subject, name))
  list(
    metric   = deriv("qmri",         metric_file[[metric]]),
    distance = deriv("distancemaps", sprintf("space-%s_desc-%s_distance", space, map)),
    tumor    = deriv("masks",        sprintf("space-%s_desc-tumor_dseg", space)),       # tumor segmentation
    wm       = deriv("masks",        sprintf("space-%s_label-NAWM_mask", space)),       # thresholded FAST WM
    hemi     = deriv("masks",        sprintf("space-%s_label-tumorhemi_mask", space))   # tumor hemisphere
  )
}

# labels of the tumor segmentation (BraTS 2023 convention; BraTS 2021 uses cet = 4)
tumor_labels <- c(necrosis = 1, t2h = 2, cet = 3)

# ================================== ANALYSIS ==================================

metrics <- c("R2s", "QSM", "ChiDia", "ChiPara", "adc", "fa")
# distance maps (desc- label of the file); isoweighted is the primary map
maps    <- c("iso", "isoweighted", "aniso", "anisoweighted")

# [VARIANT] contralateral control: maps <- "contraisoweighted" (map from the
# mirrored CET) and contralateral <- TRUE
contralateral <- FALSE

# restrict the sampling domain to the tumor hemisphere (main analysis only)
restrict_to_tumor_hemisphere <- TRUE

# landmark: "trough" for ADC, "peak" for all other metrics
landmark_type <- function(metric) if (metric == "adc") "trough" else "peak"

# initial NAWM tolerance and widening step, in metric units
nawm_tol_start <- c(adc = 2e-5, fa = 0.05, R2s = 1.5, QSM = 0.002, ChiDia = 0.002, ChiPara = 0.002)
nawm_tol_step  <- c(adc = 2e-5, fa = 0.01, R2s = 0.5, QSM = 0.001, ChiDia = 0.001, ChiPara = 0.001)
nawm_tol_max_steps <- 5

# [VARIANT] sensitivity analyses (run with n_bootstrap <- 0)
smoothing_penalty   <- 1.5   # GAM gamma; sensitivity: 1.0
tol_scale           <- 1     # multiplies tolerance start and step; sensitivity: 0.5 and 2
t2h_perturb_voxels  <- 0L    # T2H boundary: -1 = erode, +1 = dilate, 0 = none
# CET boundary sensitivity: use the distance maps recomputed from the dilated CET as `maps`

# voxel trimming (percentiles, computed per participant)
metric_clip_pct   <- c(0, 99)
distance_clip_pct <- c(1, 95)

# GAM basis: k = min(30, max(10, floor(n / 15))); raised to 45 if the k-index
# is < 0.8 and edf/k > 0.9
k_min <- 10; k_max <- 30; k_fallback <- 45; voxels_per_k <- 15
n_grid <- 200                 # grid for the derivative and landmark search

# spatial block bootstrap of the landmark
block_size_mm <- 5
n_bootstrap   <- 100          # 0 = skip

# population curves
pop_grid_n <- 100
n_pop_boot <- 1000

n_workers <- min(8, future::availableCores())   # parallel participants

# Results are cached per participant, so an interrupted run restarts where it
# stopped. Set to TRUE after changing input files, the sampling domain or the
# trimming, otherwise the cached voxels and fits are reused.
recompute <- FALSE

# ================================ SETUP =======================================

run_tag <- sprintf("g%s_tol%s_blk%smm_t2h%d", smoothing_penalty, tol_scale, block_size_mm, t2h_perturb_voxels)
voxel_cache <- file.path(output_dir, "voxel_cache")
fit_cache   <- file.path(output_dir, sprintf("fit_cache_%s_B%d", run_tag, n_bootstrap))
dir.create(voxel_cache, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_cache,   recursive = TRUE, showWarnings = FALSE)
options(future.globals.maxSize = 16 * 1024^3)
if (n_workers > 1) plan(multisession, workers = n_workers) else plan(sequential)

read_mask <- function(f) as.array(readNifti(f)) > 0

# 6-connected dilation / erosion of a binary 3D mask
shift_array <- function(a, d) {
  dims <- dim(a); out <- array(FALSE, dims)
  src <- lapply(1:3, function(ax) if (d[ax] >= 0) seq_len(dims[ax] - d[ax]) else seq.int(1 - d[ax], dims[ax]))
  dst <- lapply(1:3, function(ax) src[[ax]] + d[ax])
  out[dst[[1]], dst[[2]], dst[[3]]] <- a[src[[1]], src[[2]], src[[3]]]
  out
}
morph6 <- function(mask, n_iter, op = c("dilate", "erode")) {
  op <- match.arg(op)
  for (it in seq_len(n_iter)) {
    m <- mask
    for (o in list(c(1,0,0), c(-1,0,0), c(0,1,0), c(0,-1,0), c(0,0,1), c(0,0,-1)))
      m <- if (op == "dilate") m | shift_array(mask, o) else m & shift_array(mask, o)
    mask <- m
  }
  mask
}

# ============================ 1. VOXEL EXTRACTION =============================

extract_voxels <- function(subject, metric, map) {
  f <- input_files(subject, metric, map)
  needed <- c("metric", "distance", "wm", "tumor",
              if (!contralateral && restrict_to_tumor_hemisphere) "hemi")
  missing <- unlist(f[needed])[!file.exists(unlist(f[needed]))]
  if (length(missing)) stop("missing: ", paste(missing, collapse = "; "))

  img      <- readNifti(f$metric)
  value    <- as.array(img)
  distance <- as.array(readNifti(f$distance))
  wm       <- read_mask(f$wm)
  tumor    <- round(as.array(readNifti(f$tumor)))
  cet      <- array(tumor %in% tumor_labels[c("necrosis", "cet")], dim(tumor))   # CET + necrosis
  lesion   <- array(tumor %in% tumor_labels, dim(tumor))                         # CET + necrosis + T2H

  if (contralateral) {
    # distance map from the mirrored CET covers the contralateral hemisphere;
    # the real lesion (dilated by 2 voxels) is excluded; white matter only
    lesion <- morph6(lesion, 2L, "dilate")
    t2h    <- array(FALSE, dim(value))
    keep   <- wm & !lesion
  } else {
    if (t2h_perturb_voxels > 0L) lesion <- morph6(lesion,  t2h_perturb_voxels, "dilate")
    if (t2h_perturb_voxels < 0L) lesion <- morph6(lesion, -t2h_perturb_voxels, "erode")
    t2h  <- lesion & !cet
    keep <- !cet & (wm | t2h)                       # T2H taken whole, NAWM from the WM mask
    if (restrict_to_tumor_hemisphere) keep <- keep & read_mask(f$hemi)
  }
  valid <- keep & is.finite(value) & is.finite(distance) & distance > 0

  # trim metric outliers and extreme distances
  vq <- quantile(value[valid],    metric_clip_pct   / 100)
  dq <- quantile(distance[valid], distance_clip_pct / 100)
  valid <- valid & value > vq[1] & value < vq[2] & distance > dq[1] & distance < dq[2]

  ijk <- which(valid, arr.ind = TRUE)
  vox <- pixdim(img)[1:3]
  data.table(participant = subject,
             region       = ifelse(t2h[valid], "t2h", "nawm"),
             distance_raw = distance[valid],
             scalar_raw   = value[valid],
             vox_i = ijk[, 1], vox_j = ijk[, 2], vox_k = ijk[, 3],
             vox_dx = vox[1], vox_dy = vox[2], vox_dz = vox[3])
}

extract_voxels_cached <- function(subject, metric, map) {
  f <- file.path(voxel_cache, sprintf("%s_%s_t2h%d_%s.rds", metric, map, t2h_perturb_voxels, subject))
  if (!recompute && file.exists(f)) return(readRDS(f))
  dat <- tryCatch(extract_voxels(subject, metric, map), error = function(e) {
    message("  ", subject, ": extraction failed (", e$message, ")"); NULL })
  if (!is.null(dat) && nrow(dat) > 0L) saveRDS(dat, f)
  dat
}

# ============================ 2-3. GAM AND LANDMARK ===========================

fit_gam <- function(dat, k)
  tryCatch(gam(scalar_raw ~ s(distance_raw, bs = "cr", k = k), data = dat,
               method = "REML", gamma = smoothing_penalty), error = function(e) NULL)

# first derivative of the fitted curve (forward difference)
slope_on_grid <- function(mod, grid) {
  eps <- diff(range(grid)) / 1e5
  tryCatch(as.numeric((predict(mod, data.frame(distance_raw = grid + eps)) -
                       predict(mod, data.frame(distance_raw = grid))) / eps),
           error = function(e) NULL)
}

# locations where the slope changes sign (linear interpolation between grid points)
find_extrema <- function(slope, grid, type) {
  if (is.null(slope)) return(numeric(0))
  n <- length(slope)
  i <- if (type == "peak") which(slope[-n] > 0 & slope[-1L] <= 0 & grid[-n] > 0)
       else                which(slope[-n] < 0 & slope[-1L] >= 0 & grid[-n] > 0)
  grid[i] - slope[i] * (grid[i + 1L] - grid[i]) / (slope[i + 1L] - slope[i])
}

# Landmark = first extremum within the T2H distance range (all extrema if none
# lies there) whose fitted value is within the tolerance of the NAWM mean. The
# tolerance is widened stepwise; final fallback = extremum closest to the NAWM mean.
locate_landmark <- function(mod, grid, mean_nawm, t2h_max, tol, type) {
  ext <- find_extrema(slope_on_grid(mod, grid), grid, type)
  if (!length(ext)) return(list(loc = NA_real_, n_ext = 0L, selection = "no_extremum"))

  in_t2h <- ext[ext <= t2h_max]
  suffix <- if (length(in_t2h)) "" else "_no_t2h_range"
  cand   <- if (length(in_t2h)) in_t2h else ext
  dev    <- abs(as.numeric(predict(mod, data.frame(distance_raw = cand))) - mean_nawm)

  for (step in 0:nawm_tol_max_steps) {
    hit <- which(dev <= tol$start + step * tol$step)
    if (length(hit)) {
      label <- if (step == 0L) "initial_tol" else paste0("widened", step)
      return(list(loc = cand[hit[1L]], n_ext = length(cand), selection = paste0(label, suffix)))
    }
  }
  list(loc = cand[which.min(dev)], n_ext = length(cand), selection = paste0("closest_to_nawm", suffix))
}

fit_participant <- function(dat, tol, type) {
  grid <- seq(min(dat$distance_raw), max(dat$distance_raw), length.out = n_grid)
  k    <- min(k_max, max(k_min, floor(nrow(dat) / voxels_per_k)))
  mod  <- fit_gam(dat, k)
  if (is.null(mod)) return(NULL)

  # basis check: refit with a larger basis if k appears too small
  k_index <- tryCatch(k.check(mod)[1L, "k-index"], error = function(e) NA_real_)
  k_increased <- FALSE
  if (!is.na(k_index) && k_index < 0.8 && sum(mod$edf) / k > 0.9) {
    mod2 <- fit_gam(dat, k_fallback)
    if (!is.null(mod2)) {
      mod <- mod2; k <- k_fallback; k_increased <- TRUE
      k_index <- tryCatch(k.check(mod)[1L, "k-index"], error = function(e) NA_real_)
    }
  }

  mean_nawm <- mean(dat$scalar_raw[dat$region == "nawm"])
  t2h_max   <- if (any(dat$region == "t2h")) max(dat$distance_raw[dat$region == "t2h"]) else Inf
  landmark  <- locate_landmark(mod, grid, mean_nawm, t2h_max, tol, type)
  slope     <- slope_on_grid(mod, grid)
  local_gradient <- if (is.na(landmark$loc)) NA_real_ else mean(slope[grid <= landmark$loc])

  list(mod = mod, k = k, k_index = k_index, k_increased = k_increased, grid = grid,
       mean_nawm = mean_nawm, t2h_max = t2h_max, landmark = landmark, local_gradient = local_gradient)
}

# ============================ 4. BLOCK BOOTSTRAP ==============================

# resample whole cubes of block_size_mm, refit with the same k and grid
block_bootstrap <- function(dat, fit, tol, type) {
  bi <- floor((dat$vox_i - 1) * dat$vox_dx[1] / block_size_mm)
  bj <- floor((dat$vox_j - 1) * dat$vox_dy[1] / block_size_mm)
  bk <- floor((dat$vox_k - 1) * dat$vox_dz[1] / block_size_mm)
  blocks <- split(seq_len(nrow(dat)), as.integer(factor(bi + 1e4 * bj + 1e8 * bk)))
  locs <- rep(NA_real_, n_bootstrap)
  for (b in seq_len(n_bootstrap)) {
    idx <- unlist(blocks[sample.int(length(blocks), length(blocks), replace = TRUE)], use.names = FALSE)
    mod <- fit_gam(dat[idx], fit$k)
    if (!is.null(mod)) locs[b] <- locate_landmark(mod, fit$grid, fit$mean_nawm, fit$t2h_max, tol, type)$loc
  }
  list(locs = locs, n_blocks = length(blocks))
}

# ============================ 5. RESIDUAL DIAGNOSTICS =========================

# Moran's I of GAM residuals over face-adjacent (6-connected) voxel pairs
morans_i <- function(dat, resid) {
  ni <- max(dat$vox_i) + 2; nj <- max(dat$vox_j) + 2
  key  <- function(i, j, k) i + ni * (j + nj * k)
  keys <- key(dat$vox_i, dat$vox_j, dat$vox_k)
  z <- resid - mean(resid)
  a <- b <- integer(0)
  for (d in list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))) {
    nb <- match(key(dat$vox_i + d[1], dat$vox_j + d[2], dat$vox_k + d[3]), keys)
    ok <- which(!is.na(nb)); a <- c(a, ok); b <- c(b, nb[ok])
  }
  if (length(a) < 10L) return(NA_real_)
  (length(z) / length(a)) * sum(z[a] * z[b]) / sum(z^2)
}

# ============================ ONE PARTICIPANT =================================

analyse_participant <- function(subject, metric, map) {
  tol  <- list(start = nawm_tol_start[[metric]] * tol_scale, step = nawm_tol_step[[metric]] * tol_scale)
  type <- landmark_type(metric)
  id   <- data.table(participant = subject, metric = metric, map = map, run_tag = run_tag)

  dat <- extract_voxels_cached(subject, metric, map)
  if (is.null(dat) || nrow(dat) < 30L) return(list(row = id[, status := "insufficient_data"]))
  fit <- fit_participant(dat, tol, type)
  if (is.null(fit)) return(list(row = id[, status := "fit_failed"]))

  d_range <- diff(range(dat$distance_raw))
  bs <- if (n_bootstrap > 0L) block_bootstrap(dat, fit, tol, type) else list(locs = NA_real_, n_blocks = NA_integer_)
  ci <- if (all(is.na(bs$locs))) c(NA_real_, NA_real_) else quantile(bs$locs, c(0.025, 0.975), na.rm = TRUE)

  row <- cbind(id, data.table(
    status             = "ok",
    landmark_type      = type,
    n_voxels           = nrow(dat),
    n_voxels_nawm      = sum(dat$region == "nawm"),
    n_voxels_t2h       = sum(dat$region == "t2h"),
    distance_range     = d_range,
    selected_k         = fit$k,
    k_increased        = fit$k_increased,
    k_index            = fit$k_index,
    edf                = sum(fit$mod$edf),
    deviance_explained = summary(fit$mod)$dev.expl,
    delta_aic          = AIC(lm(scalar_raw ~ 1, data = dat)) - AIC(fit$mod),
    moran_i_resid      = morans_i(dat, residuals(fit$mod, type = "response")),
    mean_nawm          = fit$mean_nawm,
    first_peak         = fit$landmark$loc,          # landmark of the non-resampled fit
    has_peak           = !is.na(fit$landmark$loc),
    n_extrema          = fit$landmark$n_ext,
    landmark_selection = fit$landmark$selection,
    mean_deriv_to_peak = fit$local_gradient,        # local gradient
    n_blocks              = bs$n_blocks,
    n_boot_valid          = sum(!is.na(bs$locs)),
    fp_boot_median        = median(bs$locs, na.rm = TRUE),   # landmark estimate
    fp_boot_ci_lo         = ci[[1L]],
    fp_boot_ci_hi         = ci[[2L]],
    fp_boot_ci_width_norm = (ci[[2L]] - ci[[1L]]) / d_range
  ))

  # participant curve on the normalized 0-1 distance grid
  u <- seq(0, 1, length.out = pop_grid_n)
  d <- min(dat$distance_raw) + u * d_range
  fitted <- as.numeric(predict(fit$mod, data.frame(distance_raw = d)))
  curve <- data.table(participant = subject, u_idx = seq_len(pop_grid_n), u = u, distance_raw = d,
                      fitted = fitted, fitted_centered = fitted - mean(dat$scalar_raw))
  list(row = row, curve = curve)
}

analyse_participant_cached <- function(subject, metric, map) {
  f <- file.path(fit_cache, sprintf("%s_%s_%s.rds", metric, map, subject))
  if (!recompute && file.exists(f)) return(readRDS(f))
  res <- tryCatch(analyse_participant(subject, metric, map), error = function(e) {
    message("  ", subject, ": failed (", e$message, ")"); NULL })
  if (!is.null(res) && identical(res$row$status, "ok")) saveRDS(res, f)
  res
}

# ============================ 6. POPULATION CURVE =============================

population_curve <- function(curves, value_col) {
  M <- as.matrix(dcast(curves, participant ~ u_idx, value.var = value_col)[, -1, with = FALSE])
  n <- nrow(M)
  boot <- vapply(seq_len(n_pop_boot),
                 function(b) colMeans(M[sample.int(n, n, replace = TRUE), , drop = FALSE], na.rm = TRUE),
                 numeric(ncol(M)))
  ci <- apply(boot, 1L, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  data.table(scale = value_col, u = seq(0, 1, length.out = ncol(M)), mean = colMeans(M, na.rm = TRUE),
             ci_lo = ci[1L, ], ci_hi = ci[2L, ], n_participants = n)
}

# ================================== MAIN ======================================

message("run_tag: ", run_tag, " | bootstrap: ", n_bootstrap, " | contralateral: ", contralateral)

for (map in maps) {
  for (metric in metrics) {
    message(sprintf("=== %s | %s ===", metric, map))

    res <- future_map(participant_ids, function(p) analyse_participant_cached(p, metric, map),
                      .options = furrr_options(seed = TRUE), .progress = TRUE)
    res <- Filter(Negate(is.null), res)
    results <- rbindlist(lapply(res, `[[`, "row"), fill = TRUE)
    curves  <- rbindlist(lapply(res, `[[`, "curve"))
    if (!nrow(results)) next

    # short summary on screen
    ok <- results[status == "ok"]
    print(ok[, .(n = .N,
                 median_dev_expl  = round(median(deviance_explained), 3),
                 median_moran_i   = round(median(moran_i_resid, na.rm = TRUE), 3),
                 prop_initial_tol = round(mean(grepl("^initial_tol", landmark_selection)), 3),
                 prop_widened     = round(mean(grepl("^widened", landmark_selection)), 3),
                 prop_fallback    = round(mean(grepl("^closest", landmark_selection)), 3),
                 prop_full_range  = round(mean(grepl("no_t2h_range", landmark_selection)), 3))])

    stem <- file.path(output_dir, sprintf("%s_%s_%s", metric, map, run_tag))
    fwrite(results, paste0(stem, "_GAM.csv"))
    if (nrow(curves)) {
      fwrite(curves, paste0(stem, "_curves.csv"))
      set.seed(1L)
      pop <- rbind(population_curve(curves, "fitted"), population_curve(curves, "fitted_centered"))
      fwrite(pop[, `:=`(metric = metric, map = map, run_tag = run_tag)], paste0(stem, "_popcurve.csv"))
    }
  }
}
