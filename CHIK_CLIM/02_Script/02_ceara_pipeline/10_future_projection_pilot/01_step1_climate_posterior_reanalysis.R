# STEP 1 -- Climate-transmission inference audit.
#
# Re-analyses the ALREADY-FITTED, accepted climate-forced v4.9 canary
# posterior (03_Output/.../ce_climate_forced_canary_q0.05.rds) and the
# accepted non-climate v4.9 baseline posterior (q0.05). NO Stan refit is
# performed anywhere in this script -- both inputs already exist on disk.
#
# Extends (does not duplicate) 03_Output/02_ceara_pipeline/reports/
# climate_forced_v4_9/CLIMATE_V49_CANARY_ASSESSMENT.md: that report already
# covers HMC, checkpoint S_prop stability, confounding correlations, PPC and
# the beta_T/beta_P summary at 5 checkpoints. This script (a) independently
# re-derives the same core numbers directly from the posterior draws as a
# cross-check, (b) extends the correlation table to year-specific effects
# and seasonal-amplitude parameters, (c) adds an explicit numerical
# reconstruction check of the log_R0(t) equation against the paper's stated
# form, (d) adds a full 573-week (not just 5-checkpoint) S_prop comparison,
# and (e) adds a per-draw variance decomposition of log_R0(t) into
# seasonal / year / climate components.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot/00_config_future_pilot.R"))

message("[Step1] Loading accepted climate-forced and baseline fits (read-only; no refit)...")
climate  <- readRDS(ACCEPTED_CLIMATE_FIT)
baseline <- readRDS(ACCEPTED_BASELINE_FIT)
stopifnot(identical(as.Date(climate$weekly_data$week_start), as.Date(baseline$weekly_data$week_start)))
dates <- as.Date(climate$weekly_data$week_start)
N <- length(dates)
sd_clim <- climate$stan_data
year_id <- sd_clim$year_id
Y <- sd_clim$Y

# ================================================================
# 1.1 Reconstruct and VERIFY the exact implemented log_R0(t) expression
# ================================================================
message("[Step1 1.1] Extracting transformed-parameter series for equation verification...")
draws <- rstan::extract(climate$fit, pars = c(
  "alpha_R", "year_effect", "A_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
  "sigma_season_year", "beta_T", "beta_P", "log_R0", "R0_t", "climate_multiplier",
  "S_prop", "immune_prop", "phi_obs"
), permuted = TRUE)
n_draws <- length(draws$alpha_R)

seasonal_sin  <- sd_clim$seasonal_sin;  seasonal_cos  <- sd_clim$seasonal_cos
seasonal_sin2 <- sd_clim$seasonal_sin2; seasonal_cos2 <- sd_clim$seasonal_cos2
z_T_anom <- sd_clim$z_T_anom; z_P_anom <- sd_clim$z_P_anom

check_idx <- sample.int(n_draws, min(50L, n_draws))
max_abs_err <- 0
for (d in check_idx) {
  seasonal_t <- draws$A_year[d, year_id] * (draws$beta_sin1[d] * seasonal_sin + draws$beta_cos1[d] * seasonal_cos) +
    draws$beta_sin2[d] * seasonal_sin2 + draws$beta_cos2[d] * seasonal_cos2
  year_t <- draws$year_effect[d, year_id]
  climate_t <- draws$beta_T[d] * z_T_anom + draws$beta_P[d] * z_P_anom
  reconstructed <- draws$alpha_R[d] + year_t + seasonal_t + climate_t
  max_abs_err <- max(max_abs_err, max(abs(reconstructed - draws$log_R0[d, ])))
}
equation_check <- tibble(
  check = "log_R0(t) = alpha_R + year_effect[year_id(t)] + A_year[year_id(t)]*(beta_sin1*sin1+beta_cos1*cos1) + beta_sin2*sin2+beta_cos2*cos2 + beta_T*z_T_anom(t) + beta_P*z_P_anom(t)",
  n_draws_checked = length(check_idx), max_abs_error = max_abs_err,
  pass = max_abs_err < 1e-8
)
message(sprintf("[Step1 1.1] Equation reconstruction check: max abs error = %.3e across %d draws -- %s",
                 max_abs_err, length(check_idx), if (equation_check$pass) "PASS" else "FAIL"))
write_csv(equation_check, file.path(DIR_TABLES, "STEP1_equation_reconstruction_check.csv"))
if (!equation_check$pass) stop("Equation reconstruction check FAILED -- do not proceed with downstream Step 1 interpretation until resolved.")

# Document exact lag specification / climatology / standardisation actually used by the fit.
scaling_meta <- read_csv(ACCEPTED_SCALING_CSV, show_col_types = FALSE) |> dplyr::filter(lag_spec == "primary")
lag_spec_doc <- tibble(
  item = c("temperature lag window", "precipitation lag window", "temperature transform", "precipitation transform",
           "climatology reference", "standardisation", "boundary-week NA handling", "fixed q", "generation interval",
           "seasonal period (weeks)"),
  value = c("mean of weeks (t-2, t-1, t): primary lag spec, NOT re-derived here", "sum of weeks (t-6...t-2): primary lag spec, NOT re-derived here",
            "raw mean temperature (deg C)", "log1p(summed precipitation, mm)",
            "mean lagged value by (state, ISO week-of-year) across all available years, state=CE",
            sprintf("pooled mean/SD across CE/BA/RJ/MT (T_anom_mean=%.4f, T_anom_sd=%.4f, P_anom_mean=%.4f, P_anom_sd=%.4f)",
                    scaling_meta$T_anom_mean, scaling_meta$T_anom_sd, scaling_meta$P_anom_mean, scaling_meta$P_anom_sd),
            "set to 0 (assume average climate) for the earliest weeks in 2015-01/02 lacking a full lag window -- disclosed, not imputed",
            "0.05 (fixed scenario input, not estimated in this or any other current-pipeline fit)",
            "G=8 weeks, Gamma(shape=4, rate=2), discretised and renormalised to sum to 1",
            "52.1775 (365.25/7)")
)
write_csv(lag_spec_doc, file.path(DIR_TABLES, "STEP1_lag_and_standardisation_specification.csv"))

# ================================================================
# 1.2.A-C Posterior distributions of beta_T, beta_P: posterior vs prior N(0,0.15)
# ================================================================
message("[Step1 1.2] Posterior summaries for beta_T, beta_P...")
PRIOR_SD_BETA <- 0.15
beta_summary <- tibble(
  parameter = c("beta_T", "beta_P"),
  posterior_median = c(median(draws$beta_T), median(draws$beta_P)),
  posterior_mean = c(mean(draws$beta_T), mean(draws$beta_P)),
  posterior_sd = c(sd(draws$beta_T), sd(draws$beta_P)),
  lo95 = c(quantile(draws$beta_T, .025), quantile(draws$beta_P, .025)),
  hi95 = c(quantile(draws$beta_T, .975), quantile(draws$beta_P, .975)),
  prob_positive = c(mean(draws$beta_T > 0), mean(draws$beta_P > 0)),
  prob_negative = c(mean(draws$beta_T < 0), mean(draws$beta_P < 0)),
  prior_mean = 0, prior_sd = PRIOR_SD_BETA,
  posterior_sd_over_prior_sd = c(sd(draws$beta_T) / PRIOR_SD_BETA, sd(draws$beta_P) / PRIOR_SD_BETA)
)
message("Posterior/prior SD ratio < ~0.9 would indicate the data meaningfully update a weakly-informative prior; ~1.0 means the posterior is prior-dominated.")
print(as.data.frame(beta_summary), digits = 4)
write_csv(beta_summary, file.path(DIR_TABLES, "STEP1_beta_posterior_vs_prior_summary.csv"))

# Save raw draws (compact) for FIG01's prior/posterior overlay + FIG02's time series.
saveRDS(list(beta_T = draws$beta_T, beta_P = draws$beta_P, climate_multiplier = draws$climate_multiplier,
             dates = dates, prior_sd = PRIOR_SD_BETA),
        file.path(DIR_RESULTS, "step1_beta_and_multiplier_draws.rds"))

# ================================================================
# 1.2.D Posterior correlations -- EXTENDED to year-specific effects and
# seasonal-amplitude parameters (not just their scalar hyperparameters)
# ================================================================
message("[Step1 1.2.D] Posterior correlation table (extended)...")
year_effect_cors_T <- sapply(seq_len(Y), function(y) cor(draws$beta_T, draws$year_effect[, y]))
year_effect_cors_P <- sapply(seq_len(Y), function(y) cor(draws$beta_P, draws$year_effect[, y]))
A_year_mean_draw <- rowMeans(draws$A_year)
final_S_prop <- draws$S_prop[, N]; final_immune_prop <- draws$immune_prop[, N]

cor_tbl <- tibble(
  with_quantity = c("beta_P (cross-term)", "alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                     "sigma_season_year (seasonal-amplitude SD)", "mean(A_year) (seasonal-amplitude realisation)",
                     "max|cor| with any single year_effect[y]", "mean|cor| across all year_effect[y]",
                     "phi_obs (observation dispersion)", "S_prop at final week (terminal susceptibility)",
                     "immune_prop at final week (terminal immunity)"),
  cor_beta_T = c(cor(draws$beta_T, draws$beta_P), cor(draws$beta_T, draws$alpha_R), cor(draws$beta_T, draws$beta_sin1),
                 cor(draws$beta_T, draws$beta_cos1), cor(draws$beta_T, draws$beta_sin2), cor(draws$beta_T, draws$beta_cos2),
                 cor(draws$beta_T, draws$sigma_season_year), cor(draws$beta_T, A_year_mean_draw),
                 max(abs(year_effect_cors_T)), mean(abs(year_effect_cors_T)),
                 cor(draws$beta_T, draws$phi_obs), cor(draws$beta_T, final_S_prop), cor(draws$beta_T, final_immune_prop)),
  cor_beta_P = c(NA_real_, cor(draws$beta_P, draws$alpha_R), cor(draws$beta_P, draws$beta_sin1),
                 cor(draws$beta_P, draws$beta_cos1), cor(draws$beta_P, draws$beta_sin2), cor(draws$beta_P, draws$beta_cos2),
                 cor(draws$beta_P, draws$sigma_season_year), cor(draws$beta_P, A_year_mean_draw),
                 max(abs(year_effect_cors_P)), mean(abs(year_effect_cors_P)),
                 cor(draws$beta_P, draws$phi_obs), cor(draws$beta_P, final_S_prop), cor(draws$beta_P, final_immune_prop))
)
print(as.data.frame(cor_tbl), digits = 3)
write_csv(cor_tbl, file.path(DIR_TABLES, "STEP1_posterior_correlations_extended.csv"))

# Per-year breakdown (the aggregated max/mean above hides WHICH year drives
# the largest cor(beta, year_effect[y]) -- report it explicitly since it
# is materially larger than every other correlation in the table above).
calendar_years <- sort(unique(as.integer(format(dates, "%Y"))))
year_cor_tbl <- tibble(year = calendar_years, cor_beta_T_year_effect = year_effect_cors_T, cor_beta_P_year_effect = year_effect_cors_P)
print(as.data.frame(year_cor_tbl), digits = 3)
write_csv(year_cor_tbl, file.path(DIR_TABLES, "STEP1_beta_vs_year_effect_per_year_correlations.csv"))
message(sprintf("[Step1 1.2.D] Largest single-year correlation: beta_T vs year_effect[%d] = %.3f; beta_P vs year_effect[%d] = %.3f",
                 calendar_years[which.max(abs(year_effect_cors_T))], year_effect_cors_T[which.max(abs(year_effect_cors_T))],
                 calendar_years[which.max(abs(year_effect_cors_P))], year_effect_cors_P[which.max(abs(year_effect_cors_P))]))

# ================================================================
# 1.2.E-F Historical climate multiplier and its posterior uncertainty
# ================================================================
message("[Step1 1.2.E-F] Climate multiplier exp(beta_T*zT + beta_P*zP) over 2015-2025...")
mult_summary <- tibble(
  week_start = dates,
  median = apply(draws$climate_multiplier, 2, median),
  lo95 = apply(draws$climate_multiplier, 2, quantile, .025),
  hi95 = apply(draws$climate_multiplier, 2, quantile, .975)
)
write_csv(mult_summary, file.path(DIR_TABLES, "STEP1_climate_multiplier_weekly.csv"))
message(sprintf("Climate multiplier posterior median range: [%.3f, %.3f]; widest 95%% CI half-width: %.3f",
                 min(mult_summary$median), max(mult_summary$median), max((mult_summary$hi95 - mult_summary$lo95) / 2)))

# ================================================================
# 1.2.G Sensitivity of inferred S(t) to including climate -- FULL series
# ================================================================
message("[Step1 1.2.G] Full 573-week S_prop(t) comparison, baseline vs climate-forced...")
draws_base <- rstan::extract(baseline$fit, pars = c("S_prop", "immune_prop"), permuted = TRUE)
S_compare <- tibble(
  week_start = dates,
  S_prop_baseline_median = apply(draws_base$S_prop, 2, median),
  S_prop_climate_median = apply(draws$S_prop, 2, median)
) |> mutate(abs_diff = abs(S_prop_climate_median - S_prop_baseline_median),
            rel_diff_pct = 100 * abs_diff / S_prop_baseline_median)
write_csv(S_compare, file.path(DIR_TABLES, "STEP1_S_prop_full_series_comparison.csv"))
message(sprintf("Full-series S_prop relative difference: max=%.3f%%, mean=%.4f%% (5-checkpoint figure in the existing canary assessment reported max 0.27%%)",
                 max(S_compare$rel_diff_pct), mean(S_compare$rel_diff_pct)))

# ================================================================
# 1.3 Climate vs seasonality: per-draw variance decomposition of log_R0(t)
# (descriptive only -- NOT a causal attribution; components are only
# exactly additive in variance if mutually uncorrelated, which the
# Phase-1 orthogonality check found to be ~true for climate vs.
# harmonics/year but is verified again here directly)
# ================================================================
message("[Step1 1.3] Per-draw variance decomposition of log_R0(t) - alpha_R into seasonal / year / climate components...")
decomp_idx <- sample.int(n_draws, min(400L, n_draws))
decomp_rows <- vector("list", length(decomp_idx))
for (i in seq_along(decomp_idx)) {
  d <- decomp_idx[i]
  seasonal_t <- draws$A_year[d, year_id] * (draws$beta_sin1[d] * seasonal_sin + draws$beta_cos1[d] * seasonal_cos) +
    draws$beta_sin2[d] * seasonal_sin2 + draws$beta_cos2[d] * seasonal_cos2
  year_t <- draws$year_effect[d, year_id]
  climate_t <- draws$beta_T[d] * z_T_anom + draws$beta_P[d] * z_P_anom
  total_t <- seasonal_t + year_t + climate_t # = log_R0(t) - alpha_R
  v_total <- var(total_t); v_seasonal <- var(seasonal_t); v_year <- var(year_t); v_climate <- var(climate_t)
  decomp_rows[[i]] <- tibble(draw = d, var_total = v_total, var_seasonal = v_seasonal, var_year = v_year, var_climate = v_climate,
                             share_seasonal = v_seasonal / v_total, share_year = v_year / v_total, share_climate = v_climate / v_total,
                             covariance_remainder_share = 1 - (v_seasonal + v_year + v_climate) / v_total)
}
decomp_draws <- bind_rows(decomp_rows)
decomp_summary <- decomp_draws |>
  summarise(across(c(share_seasonal, share_year, share_climate, covariance_remainder_share),
                    list(median = median, lo95 = ~quantile(.x, .025), hi95 = ~quantile(.x, .975))))
message("Variance share summary (median [95% CI] across draws, descriptive decomposition, NOT causal attribution):")
print(as.data.frame(decomp_summary), digits = 3)
write_csv(decomp_draws, file.path(DIR_RESULTS, "STEP1_variance_decomposition_per_draw.csv"))
write_csv(decomp_summary, file.path(DIR_TABLES, "STEP1_variance_decomposition_summary.csv"))

# Posterior-median weekly component series for FIG04.
seasonal_med <- draws$A_year[which.min(abs(draws$alpha_R - median(draws$alpha_R))), year_id] # representative draw near posterior median alpha_R
med_draw <- which.min(abs(draws$alpha_R - median(draws$alpha_R)))
weekly_components <- tibble(
  week_start = dates, year_id = year_id,
  seasonal = draws$A_year[med_draw, year_id] * (draws$beta_sin1[med_draw] * seasonal_sin + draws$beta_cos1[med_draw] * seasonal_cos) +
    draws$beta_sin2[med_draw] * seasonal_sin2 + draws$beta_cos2[med_draw] * seasonal_cos2,
  year_effect_component = draws$year_effect[med_draw, year_id],
  climate_component = draws$beta_T[med_draw] * z_T_anom + draws$beta_P[med_draw] * z_P_anom,
  log_R0_minus_alpha = draws$log_R0[med_draw, ] - draws$alpha_R[med_draw]
)
write_csv(weekly_components, file.path(DIR_TABLES, "STEP1_log_R0_weekly_component_decomposition.csv"))

# ================================================================
# 1.4 Existing comparator / HMC diagnostics / q and lag-spec disclosures
# ================================================================
message("[Step1 1.4] Citing existing Stage-03 HMC/PPC comparator tables (not recomputed -- already validated)...")
existing_hmc_path <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/CE_canary_HMC_comparison.csv")
if (file.exists(existing_hmc_path)) {
  file.copy(existing_hmc_path, file.path(DIR_TABLES, "STEP1_HMC_comparison_from_stage03.csv"), overwrite = TRUE)
} else {
  message("  (existing HMC comparison table not found on disk -- recomputing minimal version from saved hmc gate objects)")
  write_csv(tibble(model = c("baseline", "climate"), divergences = c(baseline$hmc$divergences, climate$hmc$divergences),
                    max_rhat = c(baseline$hmc$maximum_rhat, climate$hmc$maximum_rhat)),
            file.path(DIR_TABLES, "STEP1_HMC_comparison_from_stage03.csv"))
}

disclosures <- tibble(
  item = c("q (ascertainment fraction)", "lag specification", "multi-state pooling", "q-sensitivity sweep availability"),
  status = c(
    "FIXED at 0.05 in every current-pipeline fit (this one included); not estimated, not uniquely identified by this or any other fit",
    "Only the 'primary' lag spec (T: mean t-2..t; P: sum t-6..t-2) was ever fit to the climate-forced Stan model. 'sensitivity_A' (T: t-1..t; P: t-4..t-1) and 'sensitivity_B' (T: t-3..t-1; P: t-8..t-3) anomaly series are already computed (climate_anomaly_covariates.csv) but have NEVER been fit -- this is a required, not-yet-performed sensitivity check, not an expensive refit done here per instruction",
    "No multi-state/pooled climate coefficient model exists yet (CLIMATE_V49_CANARY_ASSESSMENT.md explicitly stops before this 'Phase 4' step); beta_T/beta_P are CE-specific, single-state estimates only",
    "The existing q in {0.05,0.10,0.15,0.20,0.25,0.30} sweep (v4.7/v4.8/v4.9 baseline) contains NO climate forcing at all -- it is a q-identifiability check for the NON-climate model and must not be read as a sensitivity analysis of beta_T/beta_P"
  )
)
write_csv(disclosures, file.path(DIR_TABLES, "STEP1_disclosures_q_and_lag_spec.csv"))

message("\n[Step1] Posterior reanalysis complete. Tables written to: ", DIR_TABLES)
