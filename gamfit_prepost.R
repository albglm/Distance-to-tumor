library(data.table)
library(mgcv)
library(furrr)

plan(multisession, workers = 8)
options(future.globals.maxSize = 16 * 1024^3)

# --- settings --------------------------------------------------------------

folder_distmap_pre  <- "/data/project/TP1/res/"
folder_distmap_post <- "/data/project/TP2/res/"
df_part             <- fread("/data/project/part.csv", sep = "\t")
participant_ids     <- df_part[NP]

metri <- "R2s"
label <- "eigenv"

B_boot       <- 100
GRID_N       <- 200
K_CANDIDATES <- c(5L, 6L, 7L, 8L, 10L, 12L, 15L)
BOOT_VOXEL_N <- NULL  # set to e.g. 1000 to subsample voxels per bootstrap draw

message("loading data...")
dist_file_pre  <- paste0(folder_distmap_pre,  metri, "_distancemap_core_", label, "_r2s_ero_file.csv")
dist_file_post <- paste0(folder_distmap_post, metri, "_distancemap_core_", label, "_r2s_ero_file.csv")

pre  <- fread(dist_file_pre);  pre[,  time := "Pre"]
post <- fread(dist_file_post); post[, time := "Post"]
metrics_dt <- rbind(pre, post)
metrics_dt[, time := factor(time, levels = c("Pre", "Post"))]
setkey(metrics_dt, participant)
rm(pre, post)

# --- helpers -----------------------------------------------------------

prep_dat <- function(dt) if (is.null(dt) || nrow(dt) < 30L) NULL else dt

fit_gam_safe <- function(dat, k_val) {
  tryCatch(
    gam(scalar_raw ~ time + s(distance_raw, by = time, bs = "cr", k = k_val),
        data = dat, method = "REML"),
    error = function(e) NULL
  )
}

# derivative via finite differences, evaluated at a given timepoint level
get_deriv <- function(mod, grid, t_level) {
  tryCatch({
    eps <- diff(range(grid)) / 1e5
    lo  <- predict(mod, newdata = data.frame(distance_raw = grid,
                                             time = factor(t_level, levels = c("Pre", "Post"))))
    hi  <- predict(mod, newdata = data.frame(distance_raw = grid + eps,
                                             time = factor(t_level, levels = c("Pre", "Post"))))
    (hi - lo) / eps
  }, error = function(e) NULL)
}

# locate all +/- crossings of the derivative past distance 0 -> peaks of the curve
all_peaks <- function(dv, grid) {
  if (is.null(dv) || length(dv) < 2L) return(numeric(0))
  n     <- length(dv)
  cross <- which(dv[-n] > 0 & dv[-1L] <= 0 & grid[-n] > 0)
  if (!length(cross)) return(numeric(0))
  vapply(cross, function(i)
    grid[i] - dv[i] * (grid[i + 1L] - grid[i]) / (dv[i + 1L] - dv[i]), numeric(1))
}

# pair the Post peak closest to the first Pre peak; flag a large jump
match_peaks <- function(dv_pre, dv_post, grid, max_shift = 50) {
  peaks_pre  <- all_peaks(dv_pre,  grid)
  peaks_post <- all_peaks(dv_post, grid)
  
  if (!length(peaks_pre))
    return(list(fp_pre = NA_real_, fp_post = NA_real_,
                n_peaks_pre = 0L, n_peaks_post = length(peaks_post),
                peak_shift = NA_real_, shift_flag = TRUE))
  
  fp_pre <- peaks_pre[1L]
  
  if (!length(peaks_post))
    return(list(fp_pre = fp_pre, fp_post = NA_real_,
                n_peaks_pre = length(peaks_pre), n_peaks_post = 0L,
                peak_shift = NA_real_, shift_flag = TRUE))
  
  fp_post <- peaks_post[which.min(abs(peaks_post - fp_pre))]
  shift   <- abs(fp_post - fp_pre)
  
  list(fp_pre = fp_pre, fp_post = fp_post,
       n_peaks_pre = length(peaks_pre), n_peaks_post = length(peaks_post),
       peak_shift = round(shift, 4), shift_flag = shift > max_shift)
}

# during bootstrap, search for the Post peak near the Pre anchor rather than
# blindly taking the first one — keeps resampled Post peaks correspondence-consistent
first_peak_near <- function(dv, grid, anchor = NA_real_, window = 50) {
  if (is.null(dv) || length(dv) < 2L) return(NA_real_)
  n     <- length(dv)
  valid <- if (!is.na(anchor)) grid >= (anchor - window) & grid <= (anchor + window) else rep(TRUE, n)
  if (sum(valid) < 2L) valid <- rep(TRUE, n)
  
  cross <- which(dv[-n] > 0 & dv[-1L] <= 0 & valid[-n])
  cross <- cross[grid[cross] > 0]
  if (!length(cross)) return(NA_real_)
  i <- cross[1L]
  grid[i] - dv[i] * (grid[i + 1L] - grid[i]) / (dv[i + 1L] - dv[i])
}

mean_deriv_to <- function(dv, grid, loc) {
  if (is.na(loc) || is.null(dv)) return(NA_real_)
  idx <- grid <= loc
  if (sum(idx) < 2L) return(NA_real_)
  mean(dv[idx], na.rm = TRUE)
}

k_index_from_mod <- function(mod) {
  tryCatch({
    kc <- k.check(mod)
    list(k_index = min(kc[, "k-index"], na.rm = TRUE),
         k_pval  = max(kc[, "p-value"], na.rm = TRUE))
  }, error = function(e) list(k_index = NA_real_, k_pval = NA_real_))
}

# fit across candidate k's, using peak counts from both Pre and Post curves
# as a diagnostic; keep adequate-k-index candidates if any, then smallest k
# among those within 2 REML units of the best
select_k <- function(dat, grid) {
  candidates <- list()
  
  for (k_val in K_CANDIDATES) {
    mod <- fit_gam_safe(dat, k_val)
    if (is.null(mod)) next
    ki      <- k_index_from_mod(mod)
    dv_pre  <- get_deriv(mod, grid, "Pre")
    dv_post <- get_deriv(mod, grid, "Post")
    peaks_pre  <- all_peaks(dv_pre,  grid)
    peaks_post <- all_peaks(dv_post, grid)
    candidates[[length(candidates) + 1L]] <- list(
      k = k_val, k_index = ki$k_index, k_pval = ki$k_pval,
      edf = sum(mod$edf), reml = mod$gcv.ubre,
      n_peaks = length(peaks_pre) + length(peaks_post),
      fp_pre  = if (length(peaks_pre))  peaks_pre[1L]  else NA_real_,
      fp_post = if (length(peaks_post)) peaks_post[1L] else NA_real_)
  }
  
  if (!length(candidates)) return(NULL)
  
  adequate <- Filter(function(x) !is.na(x$k_index) && x$k_index >= 0.8, candidates)
  pool     <- if (length(adequate)) adequate else candidates
  reason   <- if (length(adequate)) "adequate" else "best_available_not_adequate"
  
  reml_vals <- vapply(pool, `[[`, numeric(1), "reml")
  close     <- which(reml_vals <= min(reml_vals, na.rm = TRUE) + 2)
  k_close   <- vapply(pool[close], `[[`, integer(1), "k")
  chosen    <- pool[[close[which.min(k_close)]]]
  chosen$k_reason <- reason
  chosen
}

bootstrap_peak <- function(dat, k_val, grid, fp_pre_anchor, B = B_boot) {
  n_vox  <- nrow(dat)
  samp_n <- if (!is.null(BOOT_VOXEL_N) && BOOT_VOXEL_N < n_vox) BOOT_VOXEL_N else n_vox
  out <- list(fp_pre = rep(NA_real_, B), fp_post = rep(NA_real_, B))
  
  for (b in seq_len(B)) {
    idx <- sample(n_vox, samp_n, replace = TRUE)
    mod <- fit_gam_safe(dat[idx, ], k_val)
    if (is.null(mod)) next
    dv_pre  <- get_deriv(mod, grid, "Pre")
    dv_post <- get_deriv(mod, grid, "Post")
    peaks_pre <- all_peaks(dv_pre, grid)
    out$fp_pre[b]  <- if (length(peaks_pre)) peaks_pre[1L] else NA_real_
    out$fp_post[b] <- first_peak_near(dv_post, grid, anchor = fp_pre_anchor, window = 50)
  }
  out
}

boot_stats <- function(bvec) {
  bci <- quantile(bvec, c(0.025, 0.975), na.rm = TRUE)
  list(median = round(median(bvec, na.rm = TRUE), 4),
       ci_lo = round(bci[1L], 4), ci_hi = round(bci[2L], 4))
}

# --- per patient -----------------------------------------------------------

compute_patient <- function(pid, dt) {
  
  dat <- prep_dat(dt[.(pid)])
  if (is.null(dat) || length(unique(dat$time)) < 2L)
    return(list(
      per_map = data.table(participant = pid, metric = metri, map = label,
                           timepoint = NA_character_, selected_k = NA_integer_,
                           k_reason = "insufficient_data", has_peak = FALSE),
      delta = data.table()))
  
  maxdist <- max(dat$distance_raw, na.rm = TRUE)
  grid    <- seq(min(dat$distance_raw), maxdist, length.out = GRID_N)
  
  k_sel <- select_k(dat, grid)
  if (is.null(k_sel))
    return(list(
      per_map = data.table(participant = pid, metric = metri, map = label,
                           timepoint = NA_character_, selected_k = NA_integer_,
                           k_reason = "all_fits_failed", has_peak = FALSE),
      delta = data.table()))
  
  mod <- fit_gam_safe(dat, k_sel$k)
  if (is.null(mod)) return(list(per_map = data.table(), delta = data.table()))
  
  ki       <- k_index_from_mod(mod)
  dev_expl <- round(summary(mod)$dev.expl, 6)
  dv_pre   <- get_deriv(mod, grid, "Pre")
  dv_post  <- get_deriv(mod, grid, "Post")
  
  matched <- match_peaks(dv_pre, dv_post, grid, max_shift = 50)
  
  boot    <- bootstrap_peak(dat, k_sel$k, grid, fp_pre_anchor = matched$fp_pre)
  bs_pre  <- boot_stats(boot$fp_pre)
  bs_post <- boot_stats(boot$fp_post)
  
  per_map <- rbindlist(lapply(c("Pre", "Post"), function(tp) {
    fp <- if (tp == "Pre") matched$fp_pre  else matched$fp_post
    dv <- if (tp == "Pre") dv_pre          else dv_post
    bs <- if (tp == "Pre") bs_pre          else bs_post
    data.table(
      participant        = pid, metric = metri, map = label, timepoint = tp,
      selected_k         = k_sel$k, k_reason = k_sel$k_reason,
      deviance_explained = dev_expl,
      edf                = round(sum(mod$edf), 2),
      edf_k_ratio        = round(sum(mod$edf) / k_sel$k, 3),
      k_index            = round(ki$k_index, 3),
      k_pval             = round(ki$k_pval, 4),
      first_peak         = round(fp, 4),
      has_peak           = !is.na(fp),
      n_peaks_pre        = matched$n_peaks_pre,
      n_peaks_post       = matched$n_peaks_post,
      peak_shift         = matched$peak_shift,
      peak_shift_flag    = matched$shift_flag,
      fp_boot_median     = bs$median,
      fp_boot_ci_lo      = bs$ci_lo,
      fp_boot_ci_hi      = bs$ci_hi,
      mean_deriv_to_peak = round(mean_deriv_to(dv, grid, fp), 6))
  }))
  
  delta <- data.table(
    participant          = pid, metric = metri, map = label,
    first_peak_pre       = round(matched$fp_pre,  4),
    first_peak_post      = round(matched$fp_post, 4),
    first_peak_diff      = round(matched$fp_post - matched$fp_pre, 4),
    first_peak_diff_norm = round((matched$fp_post - matched$fp_pre) / maxdist, 4),
    peak_shift           = matched$peak_shift,
    peak_shift_flag      = matched$shift_flag,
    fp_diff_boot_median  = round(bs_post$median - bs_pre$median, 4),
    rate_to_peak_pre     = round(mean_deriv_to(dv_pre,  grid, matched$fp_pre),  6),
    rate_to_peak_post    = round(mean_deriv_to(dv_post, grid, matched$fp_post), 6),
    rate_to_peak_diff    = round(mean_deriv_to(dv_post, grid, matched$fp_post) -
                                   mean_deriv_to(dv_pre,  grid, matched$fp_pre),  6))
  
  list(per_map = per_map, delta = delta)
}

# --- run ---------------------------------------------------------------

message("running gam fits across patients...")
t0 <- proc.time()

results_all <- future_map(
  participant_ids,
  function(pid)
    tryCatch(
      compute_patient(pid, metrics_dt),
      error = function(e) {
        message(sprintf("  %s failed: %s", pid, e$message))
        NULL
      }),
  .options  = furrr_options(seed = TRUE, globals = TRUE),
  .progress = TRUE)
names(results_all) <- participant_ids

message(sprintf("done in %.1f min", (proc.time() - t0)["elapsed"] / 60))

n_failed <- sum(vapply(results_all, is.null, logical(1)))
if (n_failed) message(sprintf("%d / %d patients failed", n_failed, length(participant_ids)))

results_ok     <- Filter(Negate(is.null), results_all)
results_permap <- rbindlist(lapply(results_ok, `[[`, "per_map"), fill = TRUE)
results_delta  <- rbindlist(lapply(results_ok, `[[`, "delta"),   fill = TRUE)

rm(metrics_dt, results_all, results_ok); gc()

message(sprintf("per-map rows: %d | delta rows: %d", nrow(results_permap), nrow(results_delta)))

# --- quick look ----------------------------------------------------------

results_permap[!is.na(selected_k), .(
  median_dev = round(median(deviance_explained, na.rm = TRUE), 3),
  median_edf = round(median(edf, na.rm = TRUE), 2)
), by = timepoint] |> print()

results_permap[, .(
  prop_has_peak   = round(mean(has_peak, na.rm = TRUE), 3),
  median_peak     = round(median(first_peak, na.rm = TRUE), 3),
  median_boot_med = round(median(fp_boot_median, na.rm = TRUE), 3),
  median_ci_width = round(median(fp_boot_ci_hi - fp_boot_ci_lo, na.rm = TRUE), 3)
), by = timepoint] |> print()

results_delta[, .(
  n                = .N,
  n_shift_flagged  = sum(peak_shift_flag, na.rm = TRUE),
  median_shift     = round(median(peak_shift, na.rm = TRUE), 2),
  median_peak_diff = round(median(first_peak_diff, na.rm = TRUE), 4),
  median_rate_diff = round(median(rate_to_peak_diff, na.rm = TRUE), 6)
)] |> print()

flagged <- results_delta[peak_shift_flag == TRUE]
if (nrow(flagged)) {
  flagged[, .(participant, first_peak_pre, first_peak_post, peak_shift)] |> print()
} else {
  message("no flagged peak correspondences")
}

# --- save ----------------------------------------------------------------

out_base <- paste0(folder_distmap_post, metri, "_", label, "_GAM_prepost")
fwrite(results_permap, paste0(out_base, "_per_map.csv"))
fwrite(results_delta,  paste0(out_base, "_delta.csv"))
message("saved: ", out_base, "_per_map.csv, ", out_base, "_delta.csv")