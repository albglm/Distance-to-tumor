################################################################################
# Plots of the distance-to-tumor qMRI profiles
################################################################################
#
# Reads the output of gam_distance_profiles.R and draws, for one distance map:
#   1. population profiles: all metrics in one grid; mean-centered participant
#      curves (grey), equal-weight population mean with 95% CI (blue), mean
#      T2H and NAWM values (orange dashed, green dotted)
#   2. participant profiles: one panel per participant for one metric, with the
#      landmark (blue point) and the participant's T2H and NAWM means
#   3. example participants: fitted profile, landmark and T2H/NAWM means for a
#      few participants, with the voxel distributions of metric value (left)
#      and distance (bottom)
#
# HOW TO RUN
#   1. Set results_dir and run_tag below (as in gam_distance_profiles.R).
#   2. Rscript plot_profiles.R
#
# Required R packages: data.table, ggplot2, patchwork
################################################################################

library(data.table)
library(ggplot2)
library(patchwork)

# ================================ SETTINGS ====================================

results_dir <- "/path/to/bids/derivatives/distance_profiles"   # output_dir of gam_distance_profiles.R
run_tag     <- "g1.5_tol1_blk5mm_t2h0"                          # main analysis
map         <- "isoweighted"
metrics     <- c("adc", "fa", "R2s", "QSM", "ChiDia", "ChiPara")

participant_metric <- "ChiDia"   # metric for plots 2 and 3
n_examples         <- 4L         # participants in plot 3, spread over the T2H-NAWM difference

figures_dir <- file.path(results_dir, "figures")

# ============================ LABELS AND STYLE ================================

metric_name <- list(R2s = "R2*", QSM = "QSM", ChiDia = quote(chi[dia]), ChiPara = quote(chi[para]),
                    adc = "ADC", fa = "FA")
metric_unit <- list(R2s = expression(s^-1), QSM = "ppm", ChiDia = "ppm", ChiPara = "ppm",
                    adc = expression(mm^2/s), fa = "unitless")
# metrics without an entry above are labeled with their code
name_of <- function(m) if (is.null(metric_name[[m]])) m else metric_name[[m]]
distance_unit <- if (map == "iso") "mm" else "index units"

T2H_COLOR <- "#D55E00"; NAWM_COLOR <- "#009E73"; MEAN_COLOR <- "#0072B2"   # Okabe-Ito
region_colors <- c(T2H = T2H_COLOR, NAWM = NAWM_COLOR)

ref_lines <- function(t2h, nawm) list(
  geom_hline(yintercept = t2h,  colour = T2H_COLOR,  linetype = "dashed", linewidth = 0.9),
  geom_hline(yintercept = nawm, colour = NAWM_COLOR, linetype = "dotted", linewidth = 0.9))

# ================================== DATA ======================================

result_file <- function(metric, type)
  file.path(results_dir, sprintf("%s_%s_%s_%s.csv", metric, map, run_tag, type))

# voxels from the cache of gam_distance_profiles.R, and T2H / NAWM means per participant
load_voxels <- function(metric) {
  files <- list.files(file.path(results_dir, "voxel_cache"),
                      sprintf("^%s_%s_t2h0_.*\\.rds$", metric, map), full.names = TRUE)
  if (!length(files)) stop("no cached voxels for ", metric, " / ", map)
  v <- rbindlist(lapply(files, function(f) readRDS(f)[, .(participant, distance_raw, scalar_raw, region)]))
  v[, region := fifelse(region == "nawm", "NAWM", "T2H")]
}

region_means <- function(v) {
  r <- dcast(v[, .(m = mean(scalar_raw)), by = .(participant, region)], participant ~ region, value.var = "m")
  r <- merge(r, v[, .(all_mean = mean(scalar_raw)), by = participant], by = "participant")
  r[, t2h_minus_nawm := T2H - NAWM][order(t2h_minus_nawm)]
}

landmarks <- function(metric, curves) {
  d <- fread(result_file(metric, "GAM"))[status == "ok"]
  d <- d[, .(participant, landmark = fcoalesce(as.numeric(fp_boot_median), first_peak))][is.finite(landmark)]
  d[, fitted := { s <- curves[participant == .BY$participant]
                  approx(s$distance_raw, s$fitted, xout = landmark, rule = 2)$y }, by = participant]
}

# ======================= 1. POPULATION PROFILES ===============================

population_panel <- function(metric) {
  curves <- fread(result_file(metric, "curves"))
  pop    <- fread(result_file(metric, "popcurve"))[scale == "fitted_centered"]
  r      <- region_means(load_voxels(metric))
  ggplot() +
    geom_line(data = curves, aes(u, fitted_centered, group = participant), colour = "grey70", linewidth = 0.3) +
    geom_ribbon(data = pop, aes(u, ymin = ci_lo, ymax = ci_hi), fill = MEAN_COLOR, alpha = 0.25) +
    geom_line(data = pop, aes(u, mean), colour = MEAN_COLOR, linewidth = 1.1) +
    ref_lines(mean(r$T2H - r$all_mean), mean(r$NAWM - r$all_mean)) +
    labs(title = name_of(metric), subtitle = metric_unit[[metric]], x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 11, colour = "grey30"))
}

plot_population <- function() {
  panels <- wrap_plots(lapply(metrics, population_panel), ncol = 3)
  y_title <- wrap_elements(full = grid::textGrob("Metric value (mean-centered)", rot = 90, gp = grid::gpar(fontsize = 14)))
  x_title <- wrap_elements(full = grid::textGrob("Normalized distance from seed boundary", gp = grid::gpar(fontsize = 14)))
  ((y_title | panels) + plot_layout(widths = c(0.03, 1))) / x_title + plot_layout(heights = c(1, 0.04))
}

# ======================= 2. PARTICIPANT PROFILES ==============================

plot_participants <- function(metric) {
  r      <- region_means(load_voxels(metric))
  curves <- merge(fread(result_file(metric, "curves")), r, by = "participant")
  lm     <- landmarks(metric, curves)
  order_ <- function(d) d[, participant := factor(participant, levels = r$participant)]
  ggplot(order_(curves), aes(distance_raw, fitted)) +
    geom_line(linewidth = 0.6) +
    geom_hline(aes(yintercept = T2H),  colour = T2H_COLOR,  linetype = "dashed", linewidth = 0.9) +
    geom_hline(aes(yintercept = NAWM), colour = NAWM_COLOR, linetype = "dotted", linewidth = 0.9) +
    geom_point(data = order_(lm), aes(landmark, fitted), colour = MEAN_COLOR, size = 2) +
    facet_wrap(~ participant, scales = "free", ncol = 5) +
    labs(title = name_of(metric), x = sprintf("Distance from seed boundary (%s)", distance_unit),
         y = metric_unit[[metric]]) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
}

# ======================= 3. EXAMPLE PARTICIPANTS ==============================

plot_examples <- function(metric, n) {
  v      <- load_voxels(metric)
  r      <- region_means(v)
  r      <- r[unique(round(seq(1, nrow(r), length.out = min(n, nrow(r)))))]
  curves <- fread(result_file(metric, "curves"))
  lm     <- landmarks(metric, curves)

  panel <- function(pid) {
    vx <- v[participant == pid]; cv <- curves[participant == pid]; rr <- r[participant == pid]
    xr <- range(cv$distance_raw); yr <- range(vx$scalar_raw)
    main <- ggplot(cv, aes(distance_raw, fitted)) + geom_line(linewidth = 0.8) +
      ref_lines(rr$T2H, rr$NAWM) +
      geom_point(data = lm[participant == pid], aes(landmark, fitted), colour = MEAN_COLOR, size = 2.6) +
      coord_cartesian(xlim = xr, ylim = yr) +
      labs(title = pid, x = sprintf("Distance from seed boundary (%s)", distance_unit), y = metric_unit[[metric]]) +
      theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold"))
    density_plot <- function(var) ggplot(vx, aes(.data[[var]], colour = region, fill = region)) +
      geom_density(alpha = 0.25, linewidth = 0.5) + scale_y_reverse() +
      scale_colour_manual(values = region_colors) + scale_fill_manual(values = region_colors) +
      theme_void() + theme(legend.position = "none")
    left   <- density_plot("scalar_raw") + coord_flip(xlim = yr)
    bottom <- density_plot("distance_raw") + coord_cartesian(xlim = xr)
    (left + main + plot_layout(widths = c(1, 3))) /
      (plot_spacer() + bottom + plot_layout(widths = c(1, 3))) + plot_layout(heights = c(3, 1))
  }
  wrap_plots(lapply(r$participant, panel), ncol = 2)
}

# ================================== RUN =======================================

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
save_figure <- function(p, name, width, height) {
  f <- file.path(figures_dir, sprintf("%s_%s", name, map))
  ggsave(paste0(f, ".pdf"),  p, width = width, height = height, device = cairo_pdf)
  ggsave(paste0(f, ".tiff"), p, width = width, height = height, dpi = 600, compression = "lzw")
  message("saved ", f, ".pdf/.tiff")
}

save_figure(plot_population(), "population_profiles", 13, 8)
save_figure(plot_participants(participant_metric), sprintf("participant_profiles_%s", participant_metric), 15, 14)
save_figure(plot_examples(participant_metric, n_examples), sprintf("example_profiles_%s", participant_metric), 12, 8)
