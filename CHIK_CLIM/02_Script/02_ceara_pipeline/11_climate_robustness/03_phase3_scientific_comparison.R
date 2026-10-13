# PHASE 3 -- Scientific comparison across M0-M4.
#
# Loads the accepted M0 (read-only) and the newly fit M1-M4 bundles
# (produced by 02_phase2_fit_M1_M4.R -- run that first). No Stan fitting
# here. Reports effect direction/magnitude/identifiability separately;
# does not interpret posterior stability as evidence of a causal effect.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/11_climate_robustness/00_config_climate_robustness.R"))

message("[Phase3] Loading M0 (accepted, read-only) + M1-M4 (newly fit) bundles...")
bundles <- list(M0 = readRDS(ACCEPTED_CLIMATE_FIT_M0))
for (m in c("M1", "M2", "M3", "M4")) {
  p <- file.path(DIR_FITS_CR, sprintf("ce_climate_robustness_%s.rds", m))
  if (!file.exists(p)) stop("Missing fit for ", m, ": ", p, " -- run 02_phase2_fit_M1_M4.R first.")
  bundles[[m]] <- readRDS(p)
}
dates_hist <- as.Date(bundles$M0$weekly_data$week_start)
N_hist <- length(dates_hist)
obs_cases <- bundles$M0$weekly_data$cases
years_hist <- as.integer(format(dates_hist, "%Y"))
hist_year_labels <- sort(unique(years_hist))

PRIOR_SD_BETA <- 0.15 # unchanged across all 5 models, per instruction

scaling_meta <- read_csv(ACCEPTED_SCALING_CSV, show_col_types = FALSE)
lag_spec_of <- c(M0 = "primary", M1 = "sensitivity_A", M2 = "sensitivity_B", M3 = "primary", M4 = "primary")
year_sd_of  <- c(M0 = 0.40, M1 = 0.40, M2 = 0.40, M3 = 0.20, M4 = 0.60)

draws_all <- lapply(names(bundles), function(m) {
  d <- rstan::extract(bundles[[m]]$fit, pars = c(
    "alpha_R", "year_effect", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year",
    "beta_T", "beta_P", "phi_obs", "S_prop", "immune_prop", "R0_t", "climate_multiplier", "C_pred"
  ), permuted = TRUE)
  d
})
names(draws_all) <- names(bundles)

# ================================================================
# 1-3. beta_T/beta_P posteriors, P(>0), prior-to-posterior updating
# ================================================================
message("[Phase3] 1-3: beta_T/beta_P posterior summaries across M0-M4...")
beta_summary <- bind_rows(lapply(names(bundles), function(m) {
  d <- draws_all[[m]]
  tibble(
    model = m, label = MODEL_SPECS[[m]]$label, lag_spec = lag_spec_of[[m]], year_effect_prior_sd = year_sd_of[[m]],
    parameter = c("beta_T", "beta_P"),
    posterior_median = c(median(d$beta_T), median(d$beta_P)), posterior_mean = c(mean(d$beta_T), mean(d$beta_P)),
    posterior_sd = c(sd(d$beta_T), sd(d$beta_P)), lo95 = c(quantile(d$beta_T, .025), quantile(d$beta_P, .025)),
    hi95 = c(quantile(d$beta_T, .975), quantile(d$beta_P, .975)),
    prob_positive = c(mean(d$beta_T > 0), mean(d$beta_P > 0)),
    posterior_sd_over_prior_sd = c(sd(d$beta_T) / PRIOR_SD_BETA, sd(d$beta_P) / PRIOR_SD_BETA)
  )
}))
print(as.data.frame(beta_summary), digits = 4)
write_csv(beta_summary, file.path(DIR_TABLES_CR, "PHASE3_beta_posterior_summary_M0_M4.csv"))

# ---- Physical-scale comparison (per instruction: do not compare raw
# z-scale coefficients across different standardisation denominators as
# though they share identical physical units). beta_per_physical_unit =
# beta_on_z_scale / covariate_sd_for_that_lag_spec -- converts to "effect
# per 1 degC of raw lagged-mean-temperature anomaly" (T) and "per 1 unit
# of log1p(lagged summed precipitation mm) anomaly" (P; precipitation's
# OWN transform is already log-scale, so a fully linear physical unit is
# not available for P -- this is the closest comparable quantity and is
# flagged as such, not presented as mm of rainfall). ----
physical_scale <- beta_summary |> left_join(scaling_meta |> select(lag_spec, T_anom_sd, P_anom_sd), by = "lag_spec") |>
  mutate(covariate_sd = ifelse(parameter == "beta_T", T_anom_sd, P_anom_sd),
         physical_unit = ifelse(parameter == "beta_T", "per 1 degC of lagged mean-temperature anomaly", "per 1 unit of log1p(lagged summed precipitation mm) anomaly (NOT raw mm -- precipitation's own transform is already log-scale)"),
         beta_per_physical_unit_median = posterior_median / covariate_sd,
         beta_per_physical_unit_lo95 = lo95 / covariate_sd, beta_per_physical_unit_hi95 = hi95 / covariate_sd) |>
  select(model, parameter, lag_spec, covariate_sd, physical_unit, posterior_median, beta_per_physical_unit_median, beta_per_physical_unit_lo95, beta_per_physical_unit_hi95)
print(as.data.frame(physical_scale), digits = 4)
write_csv(physical_scale, file.path(DIR_TABLES_CR, "PHASE3_beta_physical_scale_comparison.csv"))

# ================================================================
# 4. Climate multiplier exp(beta_T*z_T + beta_P*z_P) per model
# ================================================================
message("[Phase3] 4: Weekly climate multiplier per model...")
multiplier_weekly <- bind_rows(lapply(names(bundles), function(m) {
  mult <- draws_all[[m]]$climate_multiplier
  tibble(model = m, week_start = dates_hist, median = apply(mult, 2, median), lo95 = apply(mult, 2, quantile, .025), hi95 = apply(mult, 2, quantile, .975))
}))
write_csv(multiplier_weekly, file.path(DIR_TABLES_CR, "PHASE3_climate_multiplier_weekly_M0_M4.csv"))
multiplier_range <- multiplier_weekly |> group_by(model) |> summarise(min_median = min(median), max_median = max(median), max_width_95 = max(hi95 - lo95), .groups = "drop")
print(as.data.frame(multiplier_range), digits = 4)
write_csv(multiplier_range, file.path(DIR_TABLES_CR, "PHASE3_climate_multiplier_range_summary.csv"))

# ================================================================
# 5. Posterior correlations with year effects, especially 2022
# ================================================================
message("[Phase3] 5: Posterior correlations of beta_T/beta_P with individual year effects...")
year_cor <- bind_rows(lapply(names(bundles), function(m) {
  d <- draws_all[[m]]
  Y_m <- ncol(d$year_effect)
  yrs <- hist_year_labels # same 11 calendar years for all 5 models (same 2015-2025 window)
  stopifnot(length(yrs) == Y_m)
  tibble(model = m, year = yrs,
         cor_beta_T = sapply(seq_len(Y_m), function(y) cor(d$beta_T, d$year_effect[, y])),
         cor_beta_P = sapply(seq_len(Y_m), function(y) cor(d$beta_P, d$year_effect[, y])))
}))
write_csv(year_cor, file.path(DIR_TABLES_CR, "PHASE3_beta_vs_year_effect_correlations_M0_M4.csv"))
cor_2022 <- year_cor |> dplyr::filter(year == 2022)
message("cor(beta_T, year_effect[2022]) and cor(beta_P, year_effect[2022]) by model:")
print(as.data.frame(cor_2022), digits = 3)
write_csv(cor_2022, file.path(DIR_TABLES_CR, "PHASE3_beta_vs_2022_year_effect_correlation.csv"))

max_abs_cor <- year_cor |> group_by(model) |> summarise(max_abs_cor_T = max(abs(cor_beta_T)), year_of_max_T = year[which.max(abs(cor_beta_T))],
                                                          max_abs_cor_P = max(abs(cor_beta_P)), year_of_max_P = year[which.max(abs(cor_beta_P))], .groups = "drop")
print(as.data.frame(max_abs_cor), digits = 3)
write_csv(max_abs_cor, file.path(DIR_TABLES_CR, "PHASE3_max_abs_year_correlation_by_model.csv"))

# ================================================================
# 6. Historical R0(t) and S(t)/N(t) per model
# ================================================================
message("[Phase3] 6: Historical R0(t) and S_prop(t) medians per model...")
recon_weekly <- bind_rows(lapply(names(bundles), function(m) {
  d <- draws_all[[m]]
  tibble(model = m, week_start = dates_hist, R0_median = apply(d$R0_t, 2, median), S_prop_median = apply(d$S_prop, 2, median))
}))
write_csv(recon_weekly, file.path(DIR_TABLES_CR, "PHASE3_R0_Sprop_weekly_M0_M4.csv"))

recon_checkpoints <- recon_weekly |> dplyr::filter(week_start %in% as.Date(c("2016-12-31", "2018-12-31", "2020-12-31", "2022-12-31", "2025-12-21")))
# nearest-date match since exact dates may not land on a fitted week boundary
checkpoint_dates <- as.Date(c("2016-12-31", "2018-12-31", "2020-12-31", "2022-12-31", "2025-12-21"))
checkpoint_idx <- sapply(checkpoint_dates, function(cd) which.min(abs(dates_hist - cd)))
recon_checkpoints <- bind_rows(lapply(names(bundles), function(m) {
  d <- draws_all[[m]]
  tibble(model = m, checkpoint = dates_hist[checkpoint_idx], S_prop_median = apply(d$S_prop, 2, median)[checkpoint_idx], R0_median = apply(d$R0_t, 2, median)[checkpoint_idx])
}))
print(as.data.frame(recon_checkpoints), digits = 4)
write_csv(recon_checkpoints, file.path(DIR_TABLES_CR, "PHASE3_reconstruction_checkpoints_M0_M4.csv"))

S_vs_M0 <- recon_weekly |> select(model, week_start, S_prop_median) |> pivot_wider(names_from = model, values_from = S_prop_median) |>
  pivot_longer(c(M1, M2, M3, M4), names_to = "model", values_to = "S_prop") |> mutate(abs_diff_vs_M0 = abs(S_prop - M0), rel_diff_pct = 100 * abs_diff_vs_M0 / M0)
S_vs_M0_summary <- S_vs_M0 |> group_by(model) |> summarise(max_rel_diff_pct = max(rel_diff_pct), mean_rel_diff_pct = mean(rel_diff_pct), .groups = "drop")
print(as.data.frame(S_vs_M0_summary), digits = 4)
write_csv(S_vs_M0_summary, file.path(DIR_TABLES_CR, "PHASE3_S_prop_relative_difference_vs_M0.csv"))

# ================================================================
# 7. PPC and reconstruction error per model
# ================================================================
message("[Phase3] 7: Posterior predictive checks per model...")
ppc_summary <- bind_rows(lapply(names(bundles), function(m) {
  pred <- draws_all[[m]]$C_pred
  coverage <- mean(obs_cases >= apply(pred, 2, quantile, .025) & obs_cases <= apply(pred, 2, quantile, .975))
  rmse <- sqrt(mean((obs_cases - apply(pred, 2, median))^2))
  tibble(model = m, weekly_coverage_95_pct = 100 * coverage, weekly_RMSE = rmse)
}))
print(as.data.frame(ppc_summary), digits = 4)
write_csv(ppc_summary, file.path(DIR_TABLES_CR, "PHASE3_PPC_summary_M0_M4.csv"))

annual_ppc <- bind_rows(lapply(names(bundles), function(m) {
  pred <- draws_all[[m]]$C_pred
  annual_obs <- tapply(obs_cases, years_hist, sum)
  annual_pred <- apply(t(apply(pred, 1, function(r) tapply(r, years_hist, sum))), 2, median)
  tibble(model = m, year = as.integer(names(annual_obs)), observed = as.numeric(annual_obs), predicted_median = annual_pred, ratio = annual_pred / as.numeric(annual_obs))
}))
write_csv(annual_ppc, file.path(DIR_TABLES_CR, "PHASE3_annual_PPC_M0_M4.csv"))

# ================================================================
# HMC diagnostics for ALL 5 (M0 recomputed here for a side-by-side table; M1-M4 already saved by Phase2)
# ================================================================
hmc_all <- bind_rows(lapply(names(bundles), function(m) {
  h <- bundles[[m]]$hmc
  tibble(model = m, divergences = h$divergences, max_treedepth_hits = h$max_treedepth_hits, maximum_rhat = h$maximum_rhat,
         minimum_bulk_ess = h$minimum_bulk_ess, hmc_pass = h$hmc_pass)
}))
print(as.data.frame(hmc_all), digits = 4)
write_csv(hmc_all, file.path(DIR_TABLES_CR, "PHASE3_HMC_diagnostics_all_M0_M4.csv"))

message("\n[Phase3] Scientific comparison complete. Tables -> ", DIR_TABLES_CR)
