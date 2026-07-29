library(data.table)
library(mgcv)
library(furrr)

plan(multisession, workers = 8)
options(future.globals.maxSize = 16 * 1024^3)

# --- settings --------------------------------------------------------------

folder_distmap  <- "/data/project/TP1/res/"
df_part         <- fread("/data/project/part.csv", sep = "\t")
participant_ids <- df_part[NP]

metri <- "R2s"
label <- "eigenv"

B_boot       <- 100
GRID_N       <- 200
K_CANDIDATES <- c(5L, 6L, 7L, 8L, 10L, 12L, 15L)
DIST_PCT     <- 1
BOOT_VOXEL_N <- NULL  # set to e.g. 1000 to subsample voxels per bootstrap draw

message("loading data...")
dist_file  <- paste0(folder_distmap, metri, "_distancemap_core_", label, "_r2s_ero_file.csv")
metrics_dt <- fread(dist_file)
setkey(metrics_dt, participant)

# --- helpers -----------------------------------------------------------

# trim to 1st-99th percentile of distance, optionally cut further at DIST_PCT
prep_dat <- function(dt, dist_pct = DIST_PCT) {
  if (is.null(dt) || nrow(dt) < 30L) return(NULL)
  q   <- quantile(dt$distance_raw, c(0.01, 0.99), na.rm = TRUE)
  out <- dt[distance_raw >= q[1L] & distance_raw <= q[2L]]
  if (dist_pct < 1.0) {
    dmax <- quantile(out$distance_raw, dist_pct, na.rm = TRUE)
    out  <- out[distance_raw <= dmax]
  }
  if (nrow(out) < 30L) NULL else out
}

fit_gam_safe <- function(dat, k_val) {
  tryCatch(
    gam(scalar_raw ~ s(distance_raw, bs = "cr", k = k_val), data = dat, method = "REML"),
    error = function(e) NULL
  )
}

# derivative via finite differences on the prediction grid
get_deriv <- function(mod, grid) {
  tryCatch({
    eps   <- diff(range(grid)) / 1e5
    lo    <- predict(mod, newdata = data.frame(distance_raw = grid))
    hi    <- predict(mod, newdata = data.frame(distance_raw = grid + eps))
    (hi - lo) / eps
  }, error = function(e) NULL)
}

# locate first +/- crossing of the derivative past distance 0 -> first peak
all_peaks <- function(dv, grid) {
  if (is.null(dv) || length(dv) < 2L) return(numeric(0))
  n     <- length(dv)
  cross <- which(dv[-n] > 0 & dv[-1L] <= 0 & grid[-n] > 0)
  if (!length(cross)) return(numeric(0))
  vapply(cross, function(i)
    grid[i] - dv[i] * (grid[i + 1L] - grid[i]) / (dv[i + 1L] - dv[i]), numeric(1))
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
    list(k_index = kc[1L, "k-index"], k_pval = kc[1L, "p-value"])
  }, error = function(e) list(k_index = NA_real_, k_pval = NA_real_))
}

# fit across candidate k's, keep the ones with an adequate k-index (>=0.8) if
# any exist, then pick smallest k among those within 2 REML units of the best
select_and_fit <- function(dat) {
  grid <- seq(min(dat$distance_raw), max(dat$distance_raw), length.out = GRID_N)
  candidates <- list()
  
  for (k_val in K_CANDIDATES) {
    mod <- fit_gam_safe(dat, k_val)
    if (is.null(mod)) next
    ki    <- k_index_from_mod(mod)
    dv    <- get_deriv(mod, grid)
    peaks <- all_peaks(dv, grid)
    candidates[[length(candidates) + 1L]] <- list(
      mod = mod, k = k_val, k_index = ki$k_index, k_pval = ki$k_pval,
      edf = sum(mod$edf), n_peaks = length(peaks),
      fp = if (length(peaks)) peaks[1L] else NA_real_,
      dv = dv, grid = grid, reml = mod$gcv.ubre)
  }
  
  if (!length(candidates))
    return(list(mod = NULL, k = NA_integer_, k_reason = "all_fits_failed"))
  
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

bootstrap_patient <- function(dat, k_val, grid, B = B_boot) {
  n_vox  <- nrow(dat)
  samp_n <- if (!is.null(BOOT_VOXEL_N) && BOOT_VOXEL_N < n_vox) BOOT_VOXEL_N else n_vox
  fp_boot <- rep(NA_real_, B)
  
  for (b in seq_len(B)) {
    idx <- sample(n_vox, samp_n, replace = (samp_n == n_vox))
    mod <- fit_gam_safe(dat[idx], k_val)
    if (is.null(mod)) next
    dv <- get_deriv(mod, grid)
    if (is.null(dv)) next
    peaks <- all_peaks(dv, grid)
    fp_boot[b] <- if (length(peaks)) peaks[1L] else NA_real_
  }
  fp_boot
}

boot_stats <- function(bvec) {
  bci <- quantile(bvec, c(0.025, 0.975), na.rm = TRUE)
  list(median = round(median(bvec, na.rm = TRUE), 4),
       ci_lo = round(bci[1L], 4), ci_hi = round(bci[2L], 4))
}

# --- per patient -----------------------------------------------------------

compute_patient <- function(pid, dt) {
  
  dat <- prep_dat(dt[.(pid)])
  if (is.null(dat))
    return(data.table(participant = pid, metric = metri, map = label,
                      selected_k = NA_integer_, k_reason = "insufficient_data",
                      has_peak = FALSE))
  
  fit <- select_and_fit(dat)
  if (is.null(fit$mod))
    return(data.table(participant = pid, metric = metri, map = label,
                      selected_k = NA_integer_, k_reason = fit$k_reason,
                      has_peak = FALSE))
  
  fp_bs <- boot_stats(bootstrap_patient(dat, fit$k, fit$grid))
  md    <- mean_deriv_to(fit$dv, fit$grid, fit$fp)
  
  data.table(
    participant        = pid, metric = metri, map = label,
    selected_k         = fit$k, k_reason = fit$k_reason,
    deviance_explained = round(summary(fit$mod)$dev.expl, 4),
    edf                = round(fit$edf, 2),
    edf_k_ratio        = round(fit$edf / fit$k, 3),
    k_index            = round(fit$k_index, 3),
    k_pval             = round(fit$k_pval, 4),
    first_peak         = round(fit$fp, 4),
    has_peak           = !is.na(fit$fp),
    n_peaks            = fit$n_peaks,
    fp_boot_median     = fp_bs$median,
    fp_boot_ci_lo      = fp_bs$ci_lo,
    fp_boot_ci_hi      = fp_bs$ci_hi,
    mean_deriv_to_peak = round(md, 6))
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

results_permap <- rbindlist(Filter(Negate(is.null), results_all), fill = TRUE)
rm(metrics_dt, results_all); gc()

message(sprintf("%d rows", nrow(results_permap)))

# --- quick look ----------------------------------------------------------

results_permap[!is.na(selected_k), .(
  median_dev = round(median(deviance_explained, na.rm = TRUE), 3),
  q25_dev    = round(quantile(deviance_explained, 0.25, na.rm = TRUE), 3),
  q75_dev    = round(quantile(deviance_explained, 0.75, na.rm = TRUE), 3),
  median_edf = round(median(edf, na.rm = TRUE), 2),
  n_adequate = sum(k_reason == "adequate", na.rm = TRUE),
  n_fallback = sum(k_reason == "best_available_not_adequate", na.rm = TRUE)
)] |> print()

results_permap[, .(
  prop_has_peak   = round(mean(has_peak, na.rm = TRUE), 3),
  median_peak     = round(median(first_peak, na.rm = TRUE), 3),
  median_boot_med = round(median(fp_boot_median, na.rm = TRUE), 3),
  median_ci_width = round(median(fp_boot_ci_hi - fp_boot_ci_lo, na.rm = TRUE), 3)
)] |> print()

results_permap[!is.na(first_peak),
               .(participant, first_peak, fp_boot_median, fp_boot_ci_lo, fp_boot_ci_hi)] |> print()

out_file <- paste0(folder_distmap, metri, "_", label, "_GAM_results.csv")
fwrite(results_permap, out_file)
message("saved: ", out_file)