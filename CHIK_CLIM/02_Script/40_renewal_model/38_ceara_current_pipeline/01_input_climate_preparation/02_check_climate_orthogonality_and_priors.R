# Climate-forced v4.9 canary -- pre-fit checks requested before Phase 2:
#   1. Prior-predictive climate_multiplier using CE's ACTUAL observed
#      (z_T_anom, z_P_anom) weekly pairs (not marginal/hypothetical quantiles).
#   2. Seasonal mean of z_T_anom / z_P_anom by week_of_year (should be ~0).
#   3. cor(z_T_anom, z_P_anom).
#   4. Orthogonality of the anomalies to the EXISTING harmonic regressors
#      (seasonal_sin/cos/sin2/cos2, exactly as used in the Stan model) --
#      the critical identifiability check: if anomalies correlate strongly
#      with the harmonics already in log_R0[t], beta_T/beta_P would compete
#      with beta_sin1/beta_cos1/beta_sin2/beta_cos2 rather than adding new
#      information.
#
# No Stan fit is run here -- this is entirely R-side, using the ALREADY
# BUILT anomaly table (Phase 1) and the ALREADY FITTED CE q=0.05 bundle's
# stan_data (to get the EXACT harmonic values/week alignment used in
# fitting, not a re-derivation).

required_packages <- c("here", "dplyr", "readr", "tibble", "ggplot2", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(ggplot2); library(tidyr) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

# ---- Load CE q=0.05's exact fitting-window harmonic values + dates ----
ce_bundle <- readRDS(file.path(root, "02_Script/40_renewal_model/28_v4_9_hierarchical_seasonality/outputs/q0.05/renewal_ceara_v4_9_fit_q0.05.rds"))
ce_dates <- as.Date(ce_bundle$weekly_data$week_start)
harmonics <- tibble(week_start = ce_dates,
                     seasonal_sin = ce_bundle$stan_data$seasonal_sin, seasonal_cos = ce_bundle$stan_data$seasonal_cos,
                     seasonal_sin2 = ce_bundle$stan_data$seasonal_sin2, seasonal_cos2 = ce_bundle$stan_data$seasonal_cos2)
message(sprintf("CE fitting window: %s to %s (N=%d)", min(ce_dates), max(ce_dates), length(ce_dates)))

# ---- Load Phase-1 anomalies, primary spec, CE only, aligned to the fitting window ----
anom <- read_csv(file.path(table_dir, "climate_anomaly_covariates.csv"), show_col_types = FALSE) |>
  filter(state == "CE", lag_spec == "primary") |> mutate(week_start = as.Date(week_start))

ce <- harmonics |> inner_join(anom, by = "week_start")
message(sprintf("Matched %d / %d fitting weeks to Phase-1 anomaly rows (%d dropped, expected only for incomplete-lag-window edge weeks)",
                 nrow(ce), length(ce_dates), length(ce_dates) - nrow(ce)))
if (anyNA(ce$z_T_anom) || anyNA(ce$z_P_anom)) {
  message("NOTE: ", sum(is.na(ce$z_T_anom) | is.na(ce$z_P_anom)), " rows have NA z-anomaly (incomplete lag window at series start) -- excluded from checks below.")
  ce <- ce |> filter(!is.na(z_T_anom), !is.na(z_P_anom))
}

# ================================================================
# Check 1: seasonal mean of anomalies by week_of_year (should be ~0)
# ================================================================
seasonal_mean_check <- ce |> group_by(week_of_year) |>
  summarise(mean_T_anom = mean(T_anom), mean_P_anom = mean(P_anom), n = n(), .groups = "drop")
message("\n=== Check 1: mean anomaly by week-of-year (CE, primary spec) ===")
message(sprintf("Range of per-week-of-year mean T_anom: [%.4f, %.4f] deg C (should hug 0)", min(seasonal_mean_check$mean_T_anom), max(seasonal_mean_check$mean_T_anom)))
message(sprintf("Range of per-week-of-year mean P_anom (log1p scale): [%.4f, %.4f] (should hug 0)", min(seasonal_mean_check$mean_P_anom), max(seasonal_mean_check$mean_P_anom)))
message(sprintf("Overall mean T_anom = %.2e | overall mean P_anom = %.2e (exactly 0 up to floating point, by construction)", mean(ce$T_anom), mean(ce$P_anom)))
write_csv(seasonal_mean_check, file.path(table_dir, "CE_check1_seasonal_mean_by_week.csv"))

# ================================================================
# Check 2: cor(z_T_anom, z_P_anom)
# ================================================================
cor_TP <- cor(ce$z_T_anom, ce$z_P_anom)
message(sprintf("\n=== Check 2: cor(z_T_anom, z_P_anom) = %.4f ===", cor_TP))

# ================================================================
# Check 3 (CRITICAL): orthogonality of anomalies to the harmonic regressors
# already in log_R0[t]
# ================================================================
orthogonality <- tibble(
  regressor = c("seasonal_sin", "seasonal_cos", "seasonal_sin2", "seasonal_cos2"),
  cor_with_z_T_anom = c(cor(ce$z_T_anom, ce$seasonal_sin), cor(ce$z_T_anom, ce$seasonal_cos),
                         cor(ce$z_T_anom, ce$seasonal_sin2), cor(ce$z_T_anom, ce$seasonal_cos2)),
  cor_with_z_P_anom = c(cor(ce$z_P_anom, ce$seasonal_sin), cor(ce$z_P_anom, ce$seasonal_cos),
                         cor(ce$z_P_anom, ce$seasonal_sin2), cor(ce$z_P_anom, ce$seasonal_cos2))
)
message("\n=== Check 3 (CRITICAL): cor(anomaly, existing harmonic regressor) ===")
print(as.data.frame(orthogonality), digits = 3)
max_abs_cor <- max(abs(c(orthogonality$cor_with_z_T_anom, orthogonality$cor_with_z_P_anom)))
message(sprintf("\nMax |correlation| with any harmonic regressor = %.4f -- %s",
                 max_abs_cor, if (max_abs_cor < 0.15) "small, anomalies are close to orthogonal to the existing harmonics (good)" else "NON-TRIVIAL -- flag possible confounding with beta_sin/beta_cos before fitting"))
write_csv(orthogonality, file.path(table_dir, "CE_check3_anomaly_harmonic_orthogonality.csv"))

checks_summary <- tibble(check = c("cor(z_T_anom, z_P_anom)", "max |cor(anomaly, harmonic)|", "mean T_anom (overall)", "mean P_anom (overall)"),
                          value = c(cor_TP, max_abs_cor, mean(ce$T_anom), mean(ce$P_anom)))
write_csv(checks_summary, file.path(table_dir, "CE_identifiability_checks_summary.csv"))

# ================================================================
# Prior-predictive climate_multiplier using ACTUAL observed (z_T, z_P) pairs
# ================================================================
set.seed(20260916L)
N_DRAWS <- 4000L
beta_T_draws <- rnorm(N_DRAWS, 0, 0.15)
beta_P_draws <- rnorm(N_DRAWS, 0, 0.15)

pp <- bind_rows(lapply(seq_len(nrow(ce)), function(i) {
  mult <- exp(beta_T_draws * ce$z_T_anom[i] + beta_P_draws * ce$z_P_anom[i])
  tibble(week_start = ce$week_start[i], z_T_anom = ce$z_T_anom[i], z_P_anom = ce$z_P_anom[i],
         mult_median = median(mult), mult_lo95 = quantile(mult, .025), mult_hi95 = quantile(mult, .975))
}))
message(sprintf("\n=== Prior-predictive climate_multiplier over CE's %d actual observed weeks ===", nrow(pp)))
message(sprintf("Median multiplier range across all real weeks: [%.3f, %.3f]", min(pp$mult_median), max(pp$mult_median)))
message(sprintf("Widest 95%% prior interval observed at any real week: [%.3f, %.3f]", min(pp$mult_lo95), max(pp$mult_hi95)))
write_csv(pp, file.path(table_dir, "CE_prior_predictive_climate_multiplier.csv"))

p_pp_time <- ggplot(pp, aes(week_start)) +
  geom_ribbon(aes(ymin = mult_lo95, ymax = mult_hi95), fill = "#E69F00", alpha = 0.25) +
  geom_line(aes(y = mult_median), colour = "#E69F00", linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  labs(title = "CE prior-predictive climate_multiplier at ACTUAL observed weekly (z_T_anom, z_P_anom)",
       subtitle = sprintf("beta_T, beta_P ~ Normal(0, 0.15); %d prior draws per week; no fitting performed", N_DRAWS),
       x = NULL, y = "climate_multiplier = exp(beta_T*zT + beta_P*zP)") + theme_v4
ggsave(file.path(figure_dir, "CE_prior_predictive_climate_multiplier_timeseries.png"), p_pp_time, width = 200, height = 110, units = "mm", dpi = 300, bg = "white")

p_pp_scatter <- ggplot(pp, aes(z_T_anom, z_P_anom, colour = mult_median)) +
  geom_point(size = 1.3) +
  scale_colour_viridis_c(name = "Prior median\nclimate_multiplier", option = "plasma") +
  labs(title = "CE observed (z_T_anom, z_P_anom) pairs, coloured by prior-predictive median multiplier",
       subtitle = "Every point is a REAL observed week (primary lag spec) -- this is the actual joint distribution the prior will be evaluated over",
       x = "z_T_anom", y = "z_P_anom") + theme_v4
ggsave(file.path(figure_dir, "CE_prior_predictive_observed_zT_zP_scatter.png"), p_pp_scatter, width = 170, height = 140, units = "mm", dpi = 300, bg = "white")

message("\n[saved] CE_prior_predictive_climate_multiplier_timeseries.png, CE_prior_predictive_observed_zT_zP_scatter.png")
message("[saved] CE_check1_seasonal_mean_by_week.csv, CE_check3_anomaly_harmonic_orthogonality.csv, CE_identifiability_checks_summary.csv, CE_prior_predictive_climate_multiplier.csv")
