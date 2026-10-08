################################################################################
# Exploratory survival analysis of the distance-profile features
################################################################################
#
# Relates the landmark distance and local gradient of each qMRI metric (output
# of gam_distance_profiles.R, one distance map) to overall survival with Cox
# proportional-hazards models. All analyses are exploratory; AIC and
# concordance are in-sample.
#
#   1. Cohort: Kaplan-Meier median overall survival, reverse Kaplan-Meier median
#      follow-up, availability of the clinical variables.
#   2. Primary analysis: each z-scored feature added to the null model
#      (likelihood-ratio test, FDR across metric x feature), hazard ratio, AIC
#      gain, concordance, Schoenfeld test, influence diagnostics.
#   3. Sensitivity analyses: the same models with one clinical covariate at a
#      time, and stratified by MGMT status; events per variable reported.
#   4. Added information beyond regional summaries (T2H mean, T2H - NAWM
#      difference of the same metric, from the voxels entering the GAM fits).
#   5. MGMT status: log-rank test, proportional-hazards diagnostics, imaging
#      variables by MGMT group (Wilcoxon rank-sum, Hodges-Lehmann shift) and a
#      Kaplan-Meier figure with numbers at risk.
#
# INPUT
#   <bids_root>/participants.tsv   participant_id, include, and the clinical
#                                  columns named in SETTINGS (missing: n/a)
#   <results_dir>                  *_GAM.csv and voxel_cache/ of
#                                  gam_distance_profiles.R
#
# HOW TO RUN
#   1. Set the paths and the clinical column names below.
#   2. Rscript survival_analysis.R
#
# OUTPUT (<results_dir>/survival/)
#   clinical_availability.csv, clinical_covariates_vs_null.csv
#   survival_primary_unadjusted.csv, survival_sensitivity_adjusted.csv,
#   survival_regional_incremental.csv, survival_influence.csv
#   mgmt_ph_diagnostics.csv, mgmt_group_comparison.csv,
#   mgmt_kaplan_meier.pdf / .tiff / _risk_table.csv
#
# Required R packages: data.table, survival, ggplot2, patchwork, scales
################################################################################

library(data.table)
library(survival)
library(ggplot2)
library(patchwork)
library(scales)

# ================================ SETTINGS ====================================

bids_root   <- "/path/to/bids"
results_dir <- file.path(bids_root, "derivatives", "distance_profiles")   # output_dir of gam_distance_profiles.R
run_tag     <- "g1.5_tol1_blk5mm_t2h0"
map         <- "isoweighted"                                             # primary distance map
metrics     <- c("R2s", "QSM", "ChiDia", "ChiPara", "adc", "fa")

# columns of participants.tsv
time_col  <- "os_days"     # overall survival time (days from surgery)
event_col <- "os_event"    # 1 = death, 0 = censored
mgmt_col  <- "mgmt"        # 1 / methylated, 0 / unmethylated; n/a = unknown

# clinical covariates of the sensitivity analysis, one per model:
# label = column in participants.tsv. Numeric columns are z-scored, others are
# used as factors; columns not in the file are skipped and reported.
covariates <- c(age                   = "age",
                tumor_volume_fraction = "tumor_volume_fraction",
                kps                   = "kps",
                extent_of_resection   = "residual_ce")      # 1 = residual contrast-enhancing tumor

min_events_per_variable <- 10
relevant_delta_aic      <- 2
relevant_delta_c        <- 0.03

out_dir <- file.path(results_dir, "survival")

# ================================== DATA ======================================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
metric_label <- c(R2s = "R2*", QSM = "QSM", ChiDia = "Xdia", ChiPara = "Xpara", adc = "ADC", fa = "FA")
label_of <- function(m) ifelse(m %in% names(metric_label), metric_label[m], m)
features <- c(landmark = "fp_boot_median", gradient = "mean_deriv_to_peak")

pt <- fread(file.path(bids_root, "participants.tsv"), sep = "\t", na.strings = c("n/a", "NA", ""))
if ("include" %in% names(pt)) pt <- pt[include == 1]

clinical <- data.table(participant = pt$participant_id,
                       os_time = as.numeric(pt[[time_col]]), os_event = as.numeric(pt[[event_col]]))
clinical <- clinical[!is.na(os_time) & !is.na(os_event)]
pt <- pt[match(clinical$participant, pt$participant_id)]

mgmt_value <- tolower(trimws(as.character(if (mgmt_col %in% names(pt)) pt[[mgmt_col]] else NA)))
clinical[, mgmt := fifelse(mgmt_value %in% c("1", "methylated"), 1L,
                   fifelse(mgmt_value %in% c("0", "unmethylated"), 0L, NA_integer_))]

# covariate terms: z-scored if numeric, factor otherwise
present <- covariates[covariates %in% names(pt)]
terms <- character(0)
for (a in names(present)) {
  x <- pt[[present[[a]]]]
  if (is.numeric(x)) { clinical[, (a) := as.numeric(scale(x))]; terms[a] <- a }
  else               { clinical[, (a) := factor(x)];            terms[a] <- a }
}

message(sprintf("clinical data: %d patients, %d events", nrow(clinical), sum(clinical$os_event)))
availability <- data.table(variable = c(names(covariates), "mgmt"),
                           column = c(unname(covariates), mgmt_col),
                           in_file = c(names(covariates) %in% names(present), mgmt_col %in% names(pt)))
availability[, n_available := sapply(variable, function(v) if (v %in% names(clinical)) sum(!is.na(clinical[[v]])) else 0L)]
print(availability)
fwrite(availability, file.path(out_dir, "clinical_availability.csv"))

# imaging features (GAM output) and regional summaries (voxel cache)
load_features <- function(m) {
  f <- file.path(results_dir, sprintf("%s_%s_%s_GAM.csv", m, map, run_tag))
  if (!file.exists(f)) { message("missing: ", f); return(NULL) }
  d <- fread(f)[status == "ok"]
  # bootstrap-median landmark; the landmark of the main fit if run without bootstrap
  d[, fp_boot_median := fcoalesce(as.numeric(fp_boot_median), as.numeric(first_peak))]
  d[, c("participant", unname(features)), with = FALSE]
}
load_regional <- function(m) {
  files <- file.path(results_dir, "voxel_cache", sprintf("%s_%s_t2h0_%s.rds", m, map, clinical$participant))
  files <- files[file.exists(files)]
  if (!length(files)) { message("no voxel cache for ", m); return(NULL) }
  rbindlist(lapply(files, function(f) {
    v <- readRDS(f)
    if (!all(c("t2h", "nawm") %in% v$region)) return(NULL)
    data.table(participant = v$participant[1],
               t2h_mean    = mean(v$scalar_raw[v$region == "t2h"]),
               region_diff = mean(v$scalar_raw[v$region == "t2h"]) - mean(v$scalar_raw[v$region == "nawm"]))
  }))
}

# ================================= HELPERS ====================================

cox <- function(d, covars, strata = NULL) {
  rhs <- c(if (length(covars)) covars else "1", if (!is.null(strata)) sprintf("strata(%s)", strata))
  f <- as.formula(paste("Surv(os_time, os_event) ~", paste(rhs, collapse = " + ")))
  tryCatch(coxph(f, data = d, model = TRUE, x = TRUE), error = function(e) NULL)
}
lrt_p <- function(m0, m1) {
  a <- if (is.null(m0) || is.null(m1)) NULL else tryCatch(anova(m0, m1, test = "Chisq"), error = function(e) NULL)
  if (is.null(a)) NA_real_ else a[[grep("^P", names(a), value = TRUE)[1]]][2L]
}
delta_aic   <- function(m0, m1) if (is.null(m0) || is.null(m1)) NA_real_ else AIC(m0) - AIC(m1)
concordance <- function(m) if (is.null(m) || !length(coef(m))) NA_real_ else unname(m$concordance["concordance"])
schoenfeld  <- function(m, term) {
  z <- if (is.null(m) || !length(coef(m))) NULL else tryCatch(cox.zph(m)$table, error = function(e) NULL)
  if (!is.null(z) && term %in% rownames(z)) z[term, "p"] else NA_real_
}
hazard_ratio <- function(m, term) {
  if (is.null(m) || !(term %in% names(coef(m)))) return(list(hr = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_))
  s <- summary(m)
  list(hr = s$conf.int[term, "exp(coef)"], ci_lo = s$conf.int[term, "lower .95"],
       ci_hi = s$conf.int[term, "upper .95"], p = s$coefficients[term, "Pr(>|z|)"])
}
skewness   <- function(x) { x <- x[is.finite(x)]; mean((x - mean(x))^3) / sd(x)^3 }
n_outliers <- function(x) { q <- quantile(x, c(.25, .75), na.rm = TRUE); sum(x < q[1] - 3 * diff(q) | x > q[2] + 3 * diff(q), na.rm = TRUE) }
km_table   <- function(sf) { t <- summary(sf)$table; cols <- c("records", "events", "median", "0.95LCL", "0.95UCL")
                             if (is.matrix(t)) t[, cols, drop = FALSE] else t[cols] }

# =========================== 1. COHORT SURVIVAL ===============================

cat("\n== Overall survival (Kaplan-Meier) ==\n");          print(km_table(survfit(Surv(os_time, os_event) ~ 1, data = clinical)))
cat("\n== Median follow-up (reverse Kaplan-Meier) ==\n");  print(km_table(survfit(Surv(os_time, 1 - os_event) ~ 1, data = clinical)))

# each clinical covariate alone vs the null model (descriptive)
clin_desc <- rbindlist(lapply(names(terms), function(a) {
  d <- clinical[!is.na(get(a))]; m0 <- cox(d, character(0)); m1 <- cox(d, a)
  data.table(covariate = a, n = nrow(d), events = sum(d$os_event), df = length(coef(m1)),
             delta_aic_vs_null = delta_aic(m0, m1), lrt_p = lrt_p(m0, m1), concordance = concordance(m1),
             schoenfeld_p = tryCatch(cox.zph(m1)$table["GLOBAL", "p"], error = function(e) NA_real_))
}))
cat("\n== Clinical covariates vs null (descriptive) ==\n"); print(clin_desc)
fwrite(clin_desc, file.path(out_dir, "clinical_covariates_vs_null.csv"))

# ================ 2-4. PRIMARY, SENSITIVITY AND REGIONAL MODELS ===============

primary <- list(); adjusted <- list(); regional <- list(); influence <- list()

for (m in metrics) {
  fd <- load_features(m); if (is.null(fd)) next
  rd <- load_regional(m)
  base <- merge(clinical, fd, by = "participant")
  if (!is.null(rd)) base <- merge(base, rd, by = "participant", all.x = TRUE)

  for (fl in names(features)) {
    d <- base[!is.na(get(features[[fl]]))]
    d[, feat_z := as.numeric(scale(get(features[[fl]])))]

    # ---- primary: unadjusted ----------------------------------------------
    m0 <- cox(d, character(0)); m1 <- cox(d, "feat_z")
    primary[[length(primary) + 1L]] <- data.table(
      metric = m, feature = fl, n = nrow(d), events = sum(d$os_event), as.data.table(hazard_ratio(m1, "feat_z")),
      lrt_p = lrt_p(m0, m1), delta_aic = delta_aic(m0, m1), c_apparent = concordance(m1),
      schoenfeld_p = schoenfeld(m1, "feat_z"),
      skewness = skewness(d[[features[[fl]]]]), n_outliers_3iqr = n_outliers(d[[features[[fl]]]]))
    db <- as.numeric(residuals(m1, type = "dfbeta"))
    influence[[length(influence) + 1L]] <- data.table(
      metric = m, feature = fl, participant = d$participant, dfbeta = db,
      dfbeta_flag = abs(db) > 2 / sqrt(nrow(d)), deviance_resid = residuals(m1, type = "deviance"))

    # ---- sensitivity: one clinical covariate at a time ---------------------
    for (a in names(terms)) {
      da <- d[!is.na(get(a))]; if (nrow(da) < 10L) next
      ma0 <- cox(da, a); ma1 <- cox(da, c(a, "feat_z"))
      epv <- sum(da$os_event) / if (is.null(ma1)) NA_integer_ else length(coef(ma1))
      adjusted[[length(adjusted) + 1L]] <- data.table(
        metric = m, feature = fl, adjustment = a, n = nrow(da), events = sum(da$os_event),
        epv = epv, epv_adequate = epv >= min_events_per_variable, as.data.table(hazard_ratio(ma1, "feat_z")),
        lrt_p = lrt_p(ma0, ma1), delta_aic = delta_aic(ma0, ma1),
        cor_feat_cov = if (is.factor(da[[a]])) NA_real_ else cor(da$feat_z, da[[a]]),
        schoenfeld_p = schoenfeld(ma1, "feat_z"))
    }

    # ---- sensitivity: stratified by MGMT -----------------------------------
    ds <- d[!is.na(mgmt)]
    if (nrow(ds) >= 10L) {
      ms0 <- cox(ds, character(0), "mgmt"); ms1 <- cox(ds, "feat_z", "mgmt")
      adjusted[[length(adjusted) + 1L]] <- data.table(
        metric = m, feature = fl, adjustment = "mgmt_stratified", n = nrow(ds), events = sum(ds$os_event),
        epv = sum(ds$os_event), epv_adequate = TRUE, as.data.table(hazard_ratio(ms1, "feat_z")),
        lrt_p = lrt_p(ms0, ms1), delta_aic = delta_aic(ms0, ms1), cor_feat_cov = NA_real_,
        schoenfeld_p = schoenfeld(ms1, "feat_z"))
    }

    # ---- added information beyond regional summaries -----------------------
    if (!is.null(rd)) {
      dr <- d[complete.cases(d[, .(t2h_mean, region_diff)])]
      dr[, `:=`(t2h_mean_z = as.numeric(scale(t2h_mean)), t2h_minus_nawm_z = as.numeric(scale(region_diff)))]
      for (bl in c("t2h_mean", "t2h_minus_nawm")) {
        bt <- paste0(bl, "_z"); mr0 <- cox(dr, bt); mr1 <- cox(dr, c(bt, "feat_z")); hb <- hazard_ratio(mr0, bt)
        regional[[length(regional) + 1L]] <- data.table(
          metric = m, feature = fl, baseline = bl, n = nrow(dr), events = sum(dr$os_event),
          baseline_hr = hb$hr, baseline_ci_lo = hb$ci_lo, baseline_ci_hi = hb$ci_hi, baseline_p = hb$p,
          as.data.table(hazard_ratio(mr1, "feat_z")), lrt_p = lrt_p(mr0, mr1), delta_aic = delta_aic(mr0, mr1),
          delta_c_apparent = concordance(mr1) - concordance(mr0), cor_feat_baseline = cor(dr$feat_z, dr[[bt]]),
          schoenfeld_p = schoenfeld(mr1, "feat_z"))
      }
    }
  }
}

# ---- FDR, labels, output ----------------------------------------------------
primary <- rbindlist(primary)
if (!nrow(primary)) stop("no imaging features found in ", results_dir)
primary[, `:=`(lrt_p_fdr = p.adjust(lrt_p, "BH"), aic_favors = delta_aic >= relevant_delta_aic,
               c_favors = c_apparent - 0.5 >= relevant_delta_c)]
adjusted <- rbindlist(adjusted, fill = TRUE)
if (nrow(adjusted)) adjusted[, lrt_p_fdr := p.adjust(lrt_p, "BH"), by = adjustment]
regional <- rbindlist(regional, fill = TRUE)
if (nrow(regional)) regional[, lrt_p_fdr := p.adjust(lrt_p, "BH"), by = baseline]
influence <- rbindlist(influence)

order_rows <- function(x) { if (!nrow(x)) return(x)
  x[, metric := factor(label_of(metric), levels = label_of(metrics))]; setorder(x, metric, feature)[] }
primary <- order_rows(primary); adjusted <- order_rows(adjusted); regional <- order_rows(regional)

cat("\n== Primary analysis (unadjusted, exploratory) ==\n")
print(primary[, .(metric, feature, n, events, hr, ci_lo, ci_hi, p, lrt_p, lrt_p_fdr, delta_aic, c_apparent, schoenfeld_p)], digits = 3)
cat(sprintf("\n%d / %d features with pFDR < 0.05; %d with dAIC >= %g; %d with dC >= %g\n",
            sum(primary$lrt_p_fdr < 0.05, na.rm = TRUE), nrow(primary), sum(primary$aic_favors, na.rm = TRUE),
            relevant_delta_aic, sum(primary$c_favors, na.rm = TRUE), relevant_delta_c))

if (nrow(adjusted)) {
  cat("\n== Sensitivity analyses (one clinical covariate, or MGMT strata) ==\n")
  print(adjusted[, .(min_pfdr = min(lrt_p_fdr, na.rm = TRUE), n_sig = sum(lrt_p_fdr < 0.05, na.rm = TRUE),
                     n_models = .N, n = sprintf("%d-%d", min(n), max(n)), epv_min = round(min(epv), 1)), by = adjustment])
  hr_change <- merge(primary[, .(metric, feature, hr_unadj = hr)], adjusted[, .(metric, feature, adjustment, hr)],
                     by = c("metric", "feature"))
  print(hr_change[, .(max_abs_log_hr_change = round(max(abs(log(hr) - log(hr_unadj)), na.rm = TRUE), 3),
                      same_direction = sum(sign(log(hr)) == sign(log(hr_unadj)), na.rm = TRUE), n = .N), by = adjustment])
}
if (nrow(regional)) {
  cat("\n== Added information beyond regional summaries ==\n")
  print(regional[, .(min_pfdr = min(lrt_p_fdr, na.rm = TRUE), n_sig = sum(lrt_p_fdr < 0.05, na.rm = TRUE)), by = baseline])
}

fwrite(primary,   file.path(out_dir, "survival_primary_unadjusted.csv"))
fwrite(adjusted,  file.path(out_dir, "survival_sensitivity_adjusted.csv"))
fwrite(regional,  file.path(out_dir, "survival_regional_incremental.csv"))
fwrite(influence, file.path(out_dir, "survival_influence.csv"))

# ================================ 5. MGMT =====================================

dm <- clinical[!is.na(mgmt)]
if (nrow(dm) >= 10L && length(unique(dm$mgmt)) == 2L) {

  # ---- survival by MGMT status and proportional-hazards diagnostics --------
  cat("\n== MGMT: Kaplan-Meier medians ==\n"); print(km_table(survfit(Surv(os_time, os_event) ~ mgmt, data = dm)))
  lr   <- survdiff(Surv(os_time, os_event) ~ mgmt, data = dm)
  m_mg <- coxph(Surv(os_time, os_event) ~ mgmt, data = dm)
  zt   <- cox.zph(m_mg, transform = "km")
  m_tt <- tryCatch(coxph(Surv(os_time, os_event) ~ mgmt + tt(mgmt), data = dm, tt = function(x, t, ...) x * log(t)),
                   error = function(e) NULL)
  tt_s <- if (is.null(m_tt)) NULL else summary(m_tt)$coefficients
  mgmt_ph <- data.table(
    n = nrow(dm), events = sum(dm$os_event), n_methylated = sum(dm$mgmt == 1),
    logrank_chisq = lr$chisq, logrank_p = pchisq(lr$chisq, 1, lower.tail = FALSE), cox_hr = exp(coef(m_mg)),
    schoenfeld_chisq = zt$table["mgmt", "chisq"], schoenfeld_p = zt$table["mgmt", "p"],
    tt_coef_logt = if (is.null(tt_s)) NA_real_ else tt_s["tt(mgmt)", "coef"],
    tt_p         = if (is.null(tt_s)) NA_real_ else tt_s["tt(mgmt)", "Pr(>|z|)"])
  cat("\n== MGMT: log-rank and proportional-hazards diagnostics ==\n"); print(mgmt_ph)
  fwrite(mgmt_ph, file.path(out_dir, "mgmt_ph_diagnostics.csv"))

  # ---- imaging variables by MGMT group ---------------------------------------
  # Wilcoxon rank-sum, Hodges-Lehmann shift (methylated - unmethylated);
  # FDR within landmark/gradient and within regional summaries
  rows <- list()
  for (m in metrics) {
    fd <- load_features(m); if (is.null(fd)) next
    rd <- load_regional(m)
    d  <- merge(dm[, .(participant, mgmt)], fd, by = "participant")
    if (!is.null(rd)) d <- merge(d, rd, by = "participant", all.x = TRUE)
    vars <- c(features, t2h_mean = "t2h_mean", t2h_minus_nawm = "region_diff")
    for (vl in names(vars)) {
      if (!vars[[vl]] %in% names(d)) next
      x1 <- d[mgmt == 1, get(vars[[vl]])]; x0 <- d[mgmt == 0, get(vars[[vl]])]
      x1 <- x1[is.finite(x1)]; x0 <- x0[is.finite(x0)]
      if (length(x1) < 3L || length(x0) < 3L) next
      w <- suppressWarnings(wilcox.test(x1, x0, conf.int = TRUE, exact = FALSE))
      rows[[length(rows) + 1L]] <- data.table(
        metric = m, variable = vl, family = if (vl %in% names(features)) "landmark_gradient" else "regional",
        n_methylated = length(x1), n_unmethylated = length(x0),
        median_methylated = median(x1), median_unmethylated = median(x0),
        hl_diff = unname(w$estimate), hl_ci_lo = w$conf.int[1], hl_ci_hi = w$conf.int[2], p = w$p.value)
    }
  }
  mgmt_cmp <- rbindlist(rows)
  if (nrow(mgmt_cmp)) {
    mgmt_cmp[, p_fdr := p.adjust(p, "BH"), by = family]
    mgmt_cmp[, metric := factor(label_of(metric), levels = label_of(metrics))]
    setorder(mgmt_cmp, family, metric, variable)
    cat("\n== Imaging variables by MGMT status ==\n")
    print(mgmt_cmp[, .(min_p = min(p), min_pfdr = min(p_fdr), n_sig = sum(p_fdr < 0.05)), by = family])
    fwrite(mgmt_cmp, file.path(out_dir, "mgmt_group_comparison.csv"))
  }

  # ---- Kaplan-Meier figure with numbers at risk -------------------------------
  group_label <- c(`1` = "Methylated", `0` = "Unmethylated")
  group_color <- c(Methylated = "#E69F00", Unmethylated = "#009E73")
  risk_step   <- 250                                          # days between risk-table columns

  dm[, group := factor(group_label[as.character(mgmt)], levels = group_label)]
  sf <- survfit(Surv(os_time, os_event) ~ group, data = dm)
  strata <- rep(sub("^group=", "", names(sf$strata)), sf$strata)
  curve <- rbind(data.table(group = levels(dm$group), time = 0, surv = 1, n_censor = 0),
                 data.table(group = strata, time = sf$time, surv = sf$surv, n_censor = sf$n.censor))
  curve[, group := factor(group, levels = group_label)]
  breaks <- seq(0, ceiling(max(dm$os_time) / risk_step) * risk_step, by = risk_step)
  rk <- summary(sf, times = breaks, extend = TRUE)
  risk <- data.table(group = factor(sub("^group=", "", as.character(rk$strata)), levels = group_label),
                     time = rk$time, n_risk = rk$n.risk)
  tab <- summary(sf)$table
  legend_label <- setNames(sprintf("%s (n = %d; median %s d)", group_label, as.integer(tab[paste0("group=", group_label), "records"]),
                                   format(tab[paste0("group=", group_label), "median"], trim = TRUE, drop0trailing = TRUE)), group_label)
  x_scale <- scale_x_continuous(breaks = breaks, limits = range(breaks), expand = expansion(mult = c(0.01, 0.02)))

  p_km <- ggplot(curve, aes(time, surv, colour = group)) +
    geom_step(linewidth = 0.5) +
    geom_point(data = curve[n_censor > 0], shape = 3, size = 1.4, stroke = 0.4, show.legend = FALSE) +
    annotate("text", x = 0, y = 0.04, hjust = 0, size = 8 / .pt, label = sprintf("Log-rank p = %.3f", mgmt_ph$logrank_p)) +
    scale_colour_manual(values = group_color, labels = legend_label) + x_scale +
    scale_y_continuous(labels = label_percent(), limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Overall survival", colour = NULL) +
    theme_classic(base_size = 8) +
    theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1), legend.background = element_blank())
  p_risk <- ggplot(risk, aes(time, group, label = n_risk, colour = group)) +
    geom_text(size = 7.5 / .pt, show.legend = FALSE) +
    scale_colour_manual(values = group_color) + x_scale + scale_y_discrete(limits = rev(group_label)) +
    labs(x = "Time since surgery (days)", y = NULL, title = "Number at risk") +
    theme_classic(base_size = 8) +
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y = element_text(colour = "black", face = "bold"),
          plot.title = element_text(size = 8, face = "bold"))

  f <- file.path(out_dir, "mgmt_kaplan_meier")
  fig <- p_km / p_risk + plot_layout(heights = c(4, 1))
  ggsave(paste0(f, ".pdf"),  fig, width = 120, height = 95, units = "mm", device = cairo_pdf)
  ggsave(paste0(f, ".tiff"), fig, width = 120, height = 95, units = "mm", dpi = 600, compression = "lzw")
  fwrite(dcast(risk, group ~ time, value.var = "n_risk"), paste0(f, "_risk_table.csv"))
}

message("saved to ", out_dir)
