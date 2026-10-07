################################################################################
# Sensitivity analyses and contralateral control
################################################################################
#
# Compares runs of gam_distance_profiles.R with the main analysis. Nothing in
# the main pipeline is refit.
#
# PART 1  Landmark agreement. Each chosen sensitivity run is compared with the
#         main run, per metric, paired by participant:
#           - Spearman rho with 95% CI (participant bootstrap)
#           - median absolute landmark shift (distance units of the map)
#           - % of participants whose landmark changed
#         The landmarks of the non-resampled fits are compared, so sensitivity
#         runs can be made without bootstrap (n_bootstrap = 0).
#
# PART 2  Contralateral control. The main run (CET seeds, tumor hemisphere) is
#         compared with the contralateral run (mirrored CET seeds, contralateral
#         hemisphere):
#           - figure: A population profiles of both sides, B paired difference
#             (ipsilateral - contralateral); equal-weight mean of participant
#             curves with participant-bootstrap 95% CI
#           - table: deviance explained (%), AIC gain per 1,000 voxels and mean
#             white-matter value; median [IQR] per side, median paired
#             difference with bootstrap 95% CI, two-sided Wilcoxon signed-rank
#             test, FDR (Benjamini-Hochberg) across metrics within each endpoint
#         With ipsilateral_voxels = "nawm" the ipsilateral profiles are refit on
#         NAWM voxels only (tissue-matched control; T2H voxels dropped, distance
#         still from the real CET), from the voxel cache of the main run.
#
# HOW TO RUN
#   1. Make the runs with gam_distance_profiles.R (see README).
#   2. Set results_dir, choose the runs in `compare`, and set PART 2 below.
#   3. Rscript compare_sensitivity.R
#
# OUTPUT (results_dir)
#   sensitivity_comparison.csv
#   contralateral_<all|nawm>_fit_comparison.csv      (+ _excluded.csv)
#   contralateral_<all|nawm>_curves.csv
#   figures/contralateral_<all|nawm>.pdf / .tiff
#
# Required R packages: data.table, ggplot2, patchwork; mgcv for "nawm"
################################################################################

library(data.table)
library(ggplot2)
library(patchwork)

# ================================ SETTINGS ====================================

results_dir <- "/path/to/bids/derivatives/distance_profiles"   # output_dir of gam_distance_profiles.R
metrics     <- c("R2s", "QSM", "ChiDia", "ChiPara", "adc", "fa")

main <- list(map = "isoweighted", run_tag = "g1.5_tol1_blk5mm_t2h0")

# ---- PART 1: sensitivity runs (map and run_tag each run was saved under) -----
sensitivity_runs <- list(
  smoothing_1.0   = list(map = "isoweighted",       run_tag = "g1_tol1_blk5mm_t2h0"),
  tolerance_x0.5  = list(map = "isoweighted",       run_tag = "g1.5_tol0.5_blk5mm_t2h0"),
  tolerance_x2    = list(map = "isoweighted",       run_tag = "g1.5_tol2_blk5mm_t2h0"),
  t2h_eroded      = list(map = "isoweighted",       run_tag = "g1.5_tol1_blk5mm_t2h-1"),
  t2h_dilated     = list(map = "isoweighted",       run_tag = "g1.5_tol1_blk5mm_t2h1"),
  cet_dilated     = list(map = "isoweightedcetdil", run_tag = "g1.5_tol1_blk5mm_t2h0"))

compare <- names(sensitivity_runs)        # or e.g. c("smoothing_1.0", "cet_dilated"); character(0) to skip

# ---- PART 2: contralateral control -------------------------------------------
run_contralateral  <- TRUE
contra             <- list(map = "contraisoweighted", run_tag = main$run_tag)
ipsilateral_voxels <- "all"               # "all" (NAWM + T2H) or "nawm" (tissue-matched)
profile_value      <- "fitted_centered"   # panel A: "fitted_centered" or "fitted" (absolute values)
smoothing_penalty  <- 1.5                 # for the "nawm" refit; as in the main run

# participants left out of the tests and summarized separately, e.g. tumor
# spread into the contralateral hemisphere: c("sub-P003", "sub-P004")
excluded_participants <- character(0)

n_boot     <- 10000   # CIs of rho and of median paired differences
n_pop_boot <- 1000    # CI of the population curves

figures_dir <- file.path(results_dir, "figures")

# ================================ HELPERS =====================================

metric_label <- c(R2s = "bold(R2^'*')~(s^-1)", QSM = "bold(QSM)~(ppm)", ChiDia = "bold(chi[dia])~(ppm)",
                  ChiPara = "bold(chi[para])~(ppm)", adc = "bold(ADC)~(10^-3~mm^2/s)", fa = "bold(FA)")

result_file <- function(metric, run, type)
  file.path(results_dir, sprintf("%s_%s_%s_%s.csv", metric, run$map, run$run_tag, type))

read_result <- function(metric, run, type) {
  f <- result_file(metric, run, type)
  if (!file.exists(f)) { message("missing: ", f); return(NULL) }
  fread(f)
}

boot_ci <- function(x, stat) {
  b <- replicate(n_boot, { i <- sample.int(length(x), replace = TRUE); suppressWarnings(stat(i)) })
  quantile(b, c(0.025, 0.975), na.rm = TRUE)
}

med_iqr <- function(x) {
  q <- sapply(signif(quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE), 3), format, scientific = FALSE)
  sprintf("%s [%s, %s]", q[1], q[2], q[3])
}

set.seed(1)

# ======================== PART 1: LANDMARK AGREEMENT ==========================

landmarks <- function(metric, run) {
  d <- read_result(metric, run, "GAM")
  if (!is.null(d)) d[status == "ok", .(participant, landmark = first_peak)]
}

agreement <- rbindlist(lapply(metrics, function(m) {
  ref <- landmarks(m, main)
  if (is.null(ref)) return(NULL)
  rbindlist(lapply(compare, function(name) {
    alt <- landmarks(m, sensitivity_runs[[name]])
    if (is.null(alt)) return(NULL)
    d <- merge(ref, alt, by = "participant", suffixes = c("_main", "_sens"))[is.finite(landmark_main) & is.finite(landmark_sens)]
    if (nrow(d) < 5L) return(NULL)
    rho   <- cor(d$landmark_main, d$landmark_sens, method = "spearman")
    ci    <- boot_ci(seq_len(nrow(d)), function(i) cor(d$landmark_main[i], d$landmark_sens[i], method = "spearman"))
    shift <- abs(d$landmark_sens - d$landmark_main)
    data.table(metric = m, analysis = name, n = nrow(d),
               rho = round(rho, 2), rho_ci_lo = round(ci[[1]], 2), rho_ci_hi = round(ci[[2]], 2),
               median_abs_shift = signif(median(shift), 3),
               changed_pct = round(100 * mean(shift > 1e-6), 1))
  }))
}))

if (length(compare)) {
  if (!nrow(agreement)) stop("no sensitivity comparisons possible; check results_dir and run tags")
  fwrite(agreement, file.path(results_dir, "sensitivity_comparison.csv"))
  print(dcast(agreement, analysis ~ metric, value.var = "rho"))
  message("saved ", file.path(results_dir, "sensitivity_comparison.csv"))
}

# ======================= PART 2: CONTRALATERAL CONTROL ========================

if (run_contralateral) {

  # ---- ipsilateral fits on NAWM voxels only (ipsilateral_voxels = "nawm") ----
  # same GAM and basis-size rules as gam_distance_profiles.R, no bootstrap
  refit_nawm <- function(metric) {
    files <- list.files(file.path(results_dir, "voxel_cache"),
                        sprintf("^%s_%s_t2h0_.*\\.rds$", metric, main$map), full.names = TRUE)
    if (!length(files)) { message("no cached voxels for ", metric); return(NULL) }
    fits <- lapply(files, function(f) {
      dat <- readRDS(f)[region == "nawm"]
      if (nrow(dat) < 30L) return(NULL)
      fit <- function(k) tryCatch(mgcv::gam(scalar_raw ~ s(distance_raw, bs = "cr", k = k), data = dat,
                                            method = "REML", gamma = smoothing_penalty), error = function(e) NULL)
      k <- min(30, max(10, floor(nrow(dat) / 15))); mod <- fit(k)
      if (is.null(mod)) return(NULL)
      k_index <- tryCatch(mgcv::k.check(mod)[1L, "k-index"], error = function(e) NA_real_)
      if (!is.na(k_index) && k_index < 0.8 && sum(mod$edf) / k > 0.9) { mod2 <- fit(45); if (!is.null(mod2)) mod <- mod2 }
      u <- seq(0, 1, length.out = 100)
      d <- min(dat$distance_raw) + u * diff(range(dat$distance_raw))
      fitted <- as.numeric(predict(mod, data.frame(distance_raw = d)))
      list(gam = data.table(participant = dat$participant[1], n_voxels = nrow(dat),
                            deviance_explained = summary(mod)$dev.expl,
                            delta_aic = AIC(lm(scalar_raw ~ 1, data = dat)) - AIC(mod),
                            mean_nawm = mean(dat$scalar_raw)),
           curves = data.table(participant = dat$participant[1], u_idx = seq_along(u), u = u,
                               fitted = fitted, fitted_centered = fitted - mean(dat$scalar_raw)))
    })
    fits <- Filter(Negate(is.null), fits)
    list(gam = rbindlist(lapply(fits, `[[`, "gam")), curves = rbindlist(lapply(fits, `[[`, "curves")))
  }

  # ---- per-participant fit statistics and curves of both sides ---------------
  ok   <- function(d) if (!is.null(d)) d[status == "ok"]
  keep <- function(d) if (!is.null(d)) d[!participant %in% excluded_participants]
  ipsi <- lapply(setNames(metrics, metrics), function(m)
    if (ipsilateral_voxels == "nawm") refit_nawm(m)
    else list(gam = ok(read_result(m, main, "GAM")), curves = read_result(m, main, "curves")))
  ctrl <- lapply(setNames(metrics, metrics), function(m)
    list(gam = ok(read_result(m, contra, "GAM")), curves = read_result(m, contra, "curves")))

  # ---- population curves: A profiles of both sides, B paired difference ------
  population_curve <- function(cv, value_col) {
    M <- as.matrix(dcast(cv, participant ~ u_idx, value.var = value_col)[, -1, with = FALSE])
    n <- nrow(M)
    set.seed(1L)
    boot <- vapply(seq_len(n_pop_boot),
                   function(b) colMeans(M[sample.int(n, n, replace = TRUE), , drop = FALSE], na.rm = TRUE),
                   numeric(ncol(M)))
    ci <- apply(boot, 1L, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    data.table(u = seq(0, 1, length.out = ncol(M)), mean = colMeans(M, na.rm = TRUE),
               ci_lo = ci[1L, ], ci_hi = ci[2L, ], n = n)
  }

  curves <- rbindlist(lapply(metrics, function(m) {
    a <- keep(ipsi[[m]]$curves); b <- keep(ctrl[[m]]$curves)
    if (is.null(a) || is.null(b) || !nrow(a) || !nrow(b)) return(NULL)
    paired <- merge(a[, .(participant, u_idx, ipsi = fitted)], b[, .(participant, u_idx, contra = fitted)],
                    by = c("participant", "u_idx"))[, difference := ipsi - contra]
    rbind(population_curve(a, profile_value)[, curve := "Ipsilateral"],
          population_curve(b, profile_value)[, curve := "Contralateral"],
          population_curve(paired, "difference")[, curve := "Difference"])[, metric := m]
  }))
  if (!nrow(curves)) stop("no curves found for the contralateral control; check contra$map")

  # ---- paired fit statistics -------------------------------------------------
  endpoints <- c(dev = "Deviance explained, %", aic_per_kvox = "AIC gain per 1,000 voxels",
                 mean_nawm = "Mean white-matter value", n_voxels = "Voxels")
  stats_of <- function(d, side) if (!is.null(d) && nrow(d))
    d[, .(participant, side, dev = 100 * deviance_explained, aic_per_kvox = 1000 * delta_aic / n_voxels,
          mean_nawm, n_voxels = as.numeric(n_voxels))]
  fits <- rbindlist(lapply(metrics, function(m)
    rbind(stats_of(ipsi[[m]]$gam, "ipsi"), stats_of(ctrl[[m]]$gam, "contra"))[, metric := m]))
  paired <- dcast(melt(fits, id.vars = c("participant", "metric", "side"), variable.name = "endpoint"),
                  participant + metric + endpoint ~ side)[is.finite(ipsi) & is.finite(contra)]

  summarize_pairs <- function(d, test) d[, {
    diff <- ipsi - contra
    do_test <- test && endpoint[1] != "n_voxels"
    ci <- if (do_test) boot_ci(diff, function(i) median(diff[i])) else c(NA, NA)
    .(n = .N, ipsilateral = med_iqr(ipsi), contralateral = med_iqr(contra),
      diff_median = signif(median(diff), 3), diff_ci_lo = signif(ci[[1]], 3), diff_ci_hi = signif(ci[[2]], 3),
      p = if (do_test) suppressWarnings(wilcox.test(ipsi, contra, paired = TRUE, exact = FALSE)$p.value) else NA_real_)
  }, keyby = .(metric = factor(metric, metrics), endpoint = factor(endpoint, names(endpoints)))]

  tests <- summarize_pairs(paired[!participant %in% excluded_participants], test = TRUE)
  tests[, p_fdr := p.adjust(p, "BH"), by = endpoint]
  tests[, endpoint := endpoints[as.character(endpoint)]]

  # ---- save ------------------------------------------------------------------
  out <- file.path(results_dir, sprintf("contralateral_%s", ipsilateral_voxels))
  fwrite(tests,  paste0(out, "_fit_comparison.csv"))
  fwrite(curves, paste0(out, "_curves.csv"))
  print(tests, digits = 3)
  if (length(excluded_participants)) {
    excl <- summarize_pairs(paired[participant %in% excluded_participants], test = FALSE)
    fwrite(excl[, endpoint := endpoints[as.character(endpoint)]], paste0(out, "_fit_comparison_excluded.csv"))
  }

  # ---- figure ----------------------------------------------------------------
  pl <- copy(curves)
  pl[metric == "adc", c("mean", "ci_lo", "ci_hi") := .(mean * 1e3, ci_lo * 1e3, ci_hi * 1e3)]   # 10^-3 mm²/s
  pl[, metric := factor(metric_label[metric], levels = metric_label[metrics])]
  side_lab <- c(Ipsilateral = sprintf("Ipsilateral%s (tumor seeds)", if (ipsilateral_voxels == "nawm") " NAWM" else ""),
                Contralateral = "Contralateral (mirrored seeds)")
  pl[curve != "Difference", side := factor(side_lab[curve], levels = side_lab)]
  cols <- setNames(c("#D55E00", "#0072B2"), side_lab)

  zero_line   <- geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3)
  panel_style <- list(
    facet_wrap(~ metric, ncol = 1, scales = "free_y", labeller = label_parsed),
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(mult = 0.01)),
    theme_classic(base_size = 8),
    theme(plot.title = element_text(face = "bold"), strip.background = element_blank(),
          strip.text = element_text(face = "bold", hjust = 0), panel.spacing = unit(0.8, "lines")))

  pA <- ggplot(pl[curve != "Difference"], aes(u)) + panel_style +
    (if (profile_value == "fitted_centered") zero_line) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi, fill = side), alpha = 0.2) +
    geom_line(aes(y = mean, colour = side, linetype = side), linewidth = 0.6) +
    scale_colour_manual(values = cols, name = NULL) + scale_fill_manual(values = cols, name = NULL) +
    scale_linetype_manual(values = setNames(c("solid", "22"), side_lab), name = NULL) +
    labs(title = "A  Profiles", x = "Normalized distance from seed boundary",
         y = if (profile_value == "fitted") "Fitted value" else "Mean-centered value")

  pB <- ggplot(pl[curve == "Difference"], aes(u)) + panel_style + zero_line +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), fill = "grey50", alpha = 0.25) +
    geom_line(aes(y = mean), linewidth = 0.6) +
    labs(title = "B  Ipsilateral - contralateral", x = "Normalized distance from seed boundary", y = "Difference")

  p <- (pA | pB) + plot_layout(guides = "collect") & theme(legend.position = "top")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(figures_dir, sprintf("contralateral_%s", ipsilateral_voxels))
  ggsave(paste0(f, ".pdf"),  p, width = 180, height = 170, units = "mm", device = cairo_pdf)
  ggsave(paste0(f, ".tiff"), p, width = 180, height = 170, units = "mm", dpi = 600, compression = "lzw")
  message("saved ", out, "_*.csv and ", f, ".pdf/.tiff")
}
