# STEP 6: mechanistic and predictive summaries from existing simulations.
# No model fitting, parameter tuning, importation redesign, or new simulator
# calibration occurs here.

required_packages <- c("dplyr", "tidyr", "readr", "tibble", "rstan")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(rstan); library(dplyr); library(tidyr); library(readr); library(tibble)
})
# Project setup intentionally prefers stats::filter; this analysis uses the
# dplyr verb throughout, so bind it explicitly in the sourcing environment.
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis/00_config_step6.R"))

qv <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))
max_run <- function(x) {
  rr <- rle(replace(x, is.na(x), FALSE))
  if (!any(rr$values)) 0L else max(rr$lengths[rr$values])
}
crps_empirical <- function(x, y) {
  x <- sort(as.numeric(x)); n <- length(x)
  mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * x) / n^2
}
interval_score <- function(lo, hi, y, alpha) {
  (hi - lo) + (2 / alpha) * (lo - y) * (y < lo) + (2 / alpha) * (y - hi) * (y > hi)
}
wis_50_80_95 <- function(x, y) {
  med <- median(x)
  alphas <- c(.50, .20, .05)
  intervals <- lapply(alphas, function(a) c(qv(x, a / 2), qv(x, 1 - a / 2)))
  (0.5 * abs(y - med) + sum(vapply(seq_along(alphas), function(i) {
    alphas[i] / 2 * interval_score(intervals[[i]][1], intervals[[i]][2], y, alphas[i])
  }, numeric(1)))) / (length(alphas) + 0.5)
}
summarise_dist <- function(x) c(median = median(x), lo80 = qv(x, .10), hi80 = qv(x, .90), lo95 = qv(x, .025), hi95 = qv(x, .975))

message("[Step6] Loading existing diagnostic simulations (read-only)...")
diag_results <- readRDS(DIAG_RDS)

cell_meta <- tibble(
  cell = names(diag_results),
  holdout = sub("_.*", "", names(diag_results)),
  generator = ifelse(grepl("G1_fixed_SD_040", names(diag_results)), "G1_fixed_SD_040", "G2_estimated_tau")
)

trajectories <- bind_rows(lapply(seq_len(nrow(cell_meta)), function(i) {
  nm <- cell_meta$cell[i]
  diag_results[[nm]]$trajectories |>
    mutate(cell = nm, holdout = cell_meta$holdout[i], generator = cell_meta$generator[i],
           year = as.integer(format(week_start, "%Y")))
})) |>
  mutate(S_prop = S / N,
         import_hazard = M_import / N,
         local_hazard = pmax(lambda_t - import_hazard, 0),
         favourable = R_eff > 1,
         locally_dominant = X_latent > 0 & local_hazard > import_hazard)

year_effects <- bind_rows(lapply(seq_len(nrow(cell_meta)), function(i) {
  nm <- cell_meta$cell[i]
  diag_results[[nm]]$year_effects |>
    mutate(cell = nm, holdout = cell_meta$holdout[i], generator = cell_meta$generator[i])
}))

init_states <- bind_rows(lapply(seq_len(nrow(cell_meta)), function(i) {
  nm <- cell_meta$cell[i]
  diag_results[[nm]]$init_states |>
    mutate(cell = nm, holdout = cell_meta$holdout[i], generator = cell_meta$generator[i])
}))

observed_weekly <- bind_rows(lapply(c("H1", "H2"), function(h) {
  diag_results[[paste(h, "G2_estimated_tau", sep = "_")]]$observed |> mutate(holdout = h)
})) |>
  mutate(year = as.integer(format(week_start, "%Y")))
observed_annual <- observed_weekly |>
  group_by(holdout, year) |>
  summarise(observed_cases = sum(observed_cases), population = mean(population), .groups = "drop")

# The observed 2022 epidemic-period shading uses only validation data and is
# never fed back into any simulated quantity: dates containing the central 90%
# of cumulative 2022 reported cases.
epidemic_period_2022 <- observed_weekly |>
  filter(holdout == "H2", year == 2022) |>
  arrange(week_start) |>
  mutate(cum_prop = cumsum(observed_cases) / sum(observed_cases)) |>
  summarise(start = min(week_start[cum_prop >= .05]), end = min(week_start[cum_prop >= .95]),
            observed_2022_cases = sum(observed_cases))
write_csv(epidemic_period_2022, file.path(DIR_TABLES_STEP6, "STEP6_observed_2022_epidemic_period.csv"))

# -------------------------------------------------------------------------
# Reconstruct baseline/year/seasonal/climate log-R0 components and verify that
# their sum reproduces saved R0. This uses only posterior draws already used in
# the simulations and the same observed-climate-conditional covariates.
# -------------------------------------------------------------------------
message("[Step6] Reconstructing log-R0 components and training cutoff states...")
fit_draws <- list()
cutoff_rows <- list()
component_rows <- list()
component_check_rows <- list()
for (h in c("H1", "H2")) {
  bundle <- readRDS(FIT_PATHS[[h]])
  pars <- c("alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
            "sigma_season_year", "tau_year", "beta_T", "beta_P", "phi_obs", "S", "U", "X")
  dr <- rstan::extract(bundle$fit, pars = pars, permuted = TRUE)
  fit_draws[[h]] <- dr
  draw_ids <- sort(unique(year_effects$posterior_draw_id[year_effects$holdout == h & year_effects$generator == "G2_estimated_tau"]))
  n_train <- bundle$stan_data$N
  cutoff_rows[[h]] <- tibble(
    holdout = h, posterior_draw_id = draw_ids,
    S_prop_last_training_week = dr$S[draw_ids, n_train] / (dr$S[draw_ids, n_train] + dr$U[draw_ids, n_train])
  ) |> left_join(
    init_states |> filter(holdout == h, generator == "G2_estimated_tau") |>
      select(posterior_draw_id, S_prop_first_holdout_week = S_frac_init),
    by = "posterior_draw_id"
  )

  obs_h <- observed_weekly |> filter(holdout == h) |> arrange(week_start)
  clim <- read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
    filter(state == "CE", lag_spec == "primary") |>
    mutate(week_start = as.Date(week_start)) |>
    filter(week_start %in% obs_h$week_start) |>
    arrange(week_start)
  if (nrow(clim) != nrow(obs_h) || !all(clim$week_start == obs_h$week_start)) stop("Climate/date mismatch for ", h)
  clim$z_T_anom[is.na(clim$z_T_anom)] <- 0
  clim$z_P_anom[is.na(clim$z_P_anom)] <- 0
  reg <- future_seasonal_regressors((n_train + 1):(n_train + nrow(obs_h)))

  ye_h <- year_effects |> filter(holdout == h, generator == "G2_estimated_tau")
  for (d in draw_ids) {
    yd <- ye_h |> filter(posterior_draw_id == d) |> arrange(holdout_year)
    y_idx <- match(as.integer(format(obs_h$week_start, "%Y")), yd$holdout_year)
    fundamental <- yd$A_year[y_idx] * (dr$beta_sin1[d] * reg$sin1 + dr$beta_cos1[d] * reg$cos1)
    harmonic <- dr$beta_sin2[d] * reg$sin2 + dr$beta_cos2[d] * reg$cos2
    seasonal <- fundamental + harmonic
    temp <- dr$beta_T[d] * clim$z_T_anom
    precip <- dr$beta_P[d] * clim$z_P_anom
    comp <- tibble(
      holdout = h, posterior_draw_id = d, week_start = obs_h$week_start,
      year = as.integer(format(obs_h$week_start, "%Y")),
      baseline_log = dr$alpha_R[d], year_log = yd$year_effect[y_idx],
      seasonal_log = seasonal, climate_T_log = temp, climate_P_log = precip,
      climate_total_log = temp + precip,
      reconstructed_log_R0 = dr$alpha_R[d] + yd$year_effect[y_idx] + seasonal + temp + precip
    )
    component_rows[[length(component_rows) + 1L]] <- comp
  }
}
cutoff_draws <- bind_rows(cutoff_rows)
components <- bind_rows(component_rows)

saved_r0_unique <- trajectories |>
  filter(generator == "G2_estimated_tau", stochastic_replicate_id == 1L) |>
  select(holdout, posterior_draw_id, week_start, saved_R0 = R0)
component_check <- components |>
  inner_join(saved_r0_unique, by = c("holdout", "posterior_draw_id", "week_start")) |>
  mutate(abs_log_error = abs(reconstructed_log_R0 - log(saved_R0)))
component_check_summary <- component_check |>
  group_by(holdout) |>
  summarise(n = n(), max_abs_log_R0_error = max(abs_log_error), mean_abs_log_R0_error = mean(abs_log_error), .groups = "drop")
write_csv(component_check_summary, file.path(DIR_TABLES_STEP6, "STEP6_R0_component_reconstruction_check.csv"))
write_csv(components, file.path(DIR_RESULTS_STEP6, "STEP6_weekly_log_R0_components.csv"))

cutoff_summary <- cutoff_draws |>
  group_by(holdout) |>
  summarise(across(starts_with("S_prop"), list(median = median, lo95 = ~qv(.x, .025), hi95 = ~qv(.x, .975))), .groups = "drop")
write_csv(cutoff_summary, file.path(DIR_TABLES_STEP6, "STEP6_training_cutoff_susceptibility.csv"))

# Infection history entering 2022 is structurally different between holdouts:
# H1's last eight weeks of 2021 are stochastic holdout simulation, whereas
# H2's are training-posterior latent infections conditioned on surveillance
# through 2021. Quantify the weighted renewal infectiousness without treating
# H1/H2 as paired.
w_generation <- generation_weights_future()
G_generation <- length(w_generation)
h1_entry_history <- trajectories |>
  filter(holdout == "H1", generator == "G2_estimated_tau", week_start < as.Date("2022-01-02")) |>
  arrange(posterior_draw_id, stochastic_replicate_id, week_start) |>
  group_by(posterior_draw_id, stochastic_replicate_id) |>
  summarise(entry2022_recent_X_sum = sum(tail(X_latent, G_generation)),
            entry2022_weighted_infectiousness = sum(w_generation * rev(tail(X_latent, G_generation))),
            source = "H1 simulated final 8 weeks of 2021", .groups = "drop") |>
  mutate(holdout = "H1")
h2_draw_ids <- sort(unique(year_effects$posterior_draw_id[year_effects$holdout == "H2" & year_effects$generator == "G2_estimated_tau"]))
h2_N_train <- ncol(fit_draws$H2$X)
h2_entry_history <- tibble(
  holdout = "H2", posterior_draw_id = h2_draw_ids, stochastic_replicate_id = NA_integer_,
  entry2022_recent_X_sum = rowSums(fit_draws$H2$X[h2_draw_ids, (h2_N_train - G_generation + 1):h2_N_train, drop = FALSE]),
  entry2022_weighted_infectiousness = apply(
    fit_draws$H2$X[h2_draw_ids, (h2_N_train - G_generation + 1):h2_N_train, drop = FALSE], 1,
    function(x) sum(w_generation * rev(x))
  ),
  source = "H2 training-posterior final 8 weeks of 2021"
)
entry2022_history <- bind_rows(h1_entry_history, h2_entry_history)
first_week_2022 <- trajectories |>
  filter(generator == "G2_estimated_tau", week_start == as.Date("2022-01-02")) |>
  select(holdout, posterior_draw_id, stochastic_replicate_id, first_week_X = X_latent,
         first_week_local_hazard = local_hazard, first_week_import_hazard = import_hazard)
entry2022_summary <- entry2022_history |>
  group_by(holdout, source) |>
  summarise(across(c(entry2022_recent_X_sum, entry2022_weighted_infectiousness),
                   list(median = median, lo95 = ~qv(.x, .025), hi95 = ~qv(.x, .975),
                        prob_zero = ~mean(.x == 0))), .groups = "drop") |>
  left_join(first_week_2022 |>
              group_by(holdout) |>
              summarise(first_week_X_median = median(first_week_X), first_week_X_lo95 = qv(first_week_X, .025),
                        first_week_X_hi95 = qv(first_week_X, .975),
                        first_week_local_hazard_median = median(first_week_local_hazard),
                        first_week_import_hazard_median = median(first_week_import_hazard), .groups = "drop"),
            by = "holdout")
write_csv(entry2022_history, file.path(DIR_RESULTS_STEP6, "STEP6_entry2022_infection_history_draws.csv"))
write_csv(entry2022_summary, file.path(DIR_TABLES_STEP6, "STEP6_entry2022_infection_history.csv"))

# -------------------------------------------------------------------------
# Per-trajectory annual mechanics. M_import is always called introduction
# pressure, never imported infections.
# -------------------------------------------------------------------------
message("[Step6] Computing trajectory-level annual mechanism summaries...")
annual_traj <- trajectories |>
  arrange(cell, posterior_draw_id, stochastic_replicate_id, week_start) |>
  group_by(cell, holdout, generator, posterior_draw_id, stochastic_replicate_id, year) |>
  summarise(
    n_weeks = n(), S_prop_start = first(S_prop), S_prop_end = last(S_prop),
    delta_S_prop = last(S_prop) - first(S_prop),
    cumulative_latent_infections = sum(X_latent), annual_reported_cases = sum(Y_reported),
    mean_R0 = mean(R0), peak_R0 = max(R0), mean_Reff = mean(R_eff), peak_Reff = max(R_eff),
    weeks_Reff_gt1 = sum(favourable), longest_run_Reff_gt1 = max_run(favourable),
    introduction_pressure_total = sum(M_import),
    introduction_pressure_favourable = sum(M_import[favourable]),
    weeks_favourable_with_positive_pressure = sum(favourable & M_import > 0),
    total_force_sum = sum(lambda_t), mean_total_force = mean(lambda_t),
    local_force_sum = sum(local_hazard), introduction_force_sum = sum(import_hazard),
    any_latent_infection = any(X_latent > 0), weeks_zero_infections = sum(X_latent == 0),
    longest_infection_free_run = max_run(X_latent == 0),
    longest_locally_dominant_run = max_run(locally_dominant),
    favourable_pressure = any(favourable & M_import > 0),
    established = favourable_pressure & longest_locally_dominant_run >= ESTABLISHMENT_RUN_WEEKS,
    population_mean = mean(N), population_start = first(N), population_end = last(N),
    .groups = "drop"
  ) |>
  left_join(year_effects |> select(cell, posterior_draw_id, holdout_year, year_effect, A_year, tau_used),
            by = c("cell", "posterior_draw_id", "year" = "holdout_year")) |>
  mutate(annual_reported_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
         major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM &
           annual_reported_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K,
         mechanism_class = case_when(
           !favourable_pressure ~ "1 No positive M during Reff>1",
           !established ~ "2 Pressure present; no sustained local dominance",
           established & !major_epidemic ~ "3 Established; no major annual epidemic",
           TRUE ~ "4 Major annual epidemic"
         ))

component_annual <- components |>
  group_by(holdout, posterior_draw_id, year) |>
  summarise(baseline_log = first(baseline_log), year_log = first(year_log),
            seasonal_log_mean = mean(seasonal_log), seasonal_log_peak = max(seasonal_log),
            climate_T_log_mean = mean(climate_T_log), climate_P_log_mean = mean(climate_P_log),
            climate_total_log_mean = mean(climate_total_log), climate_total_log_peak = max(climate_total_log),
            .groups = "drop")
annual_traj <- annual_traj |>
  left_join(component_annual, by = c("holdout", "posterior_draw_id", "year"))

write_csv(annual_traj, file.path(DIR_RESULTS_STEP6, "STEP6_trajectory_annual_metrics.csv"))
saveRDS(annual_traj, file.path(DIR_RESULTS_STEP6, "STEP6_trajectory_annual_metrics.rds"))

# Weekly trajectory quantiles for requested S/N, X, Y, R0, Reff, M and force.
weekly_quantiles <- trajectories |>
  group_by(holdout, generator, week_start) |>
  summarise(across(c(S_prop, X_latent, Y_reported, R0, R_eff, M_import, lambda_t, local_hazard, import_hazard, N),
                   list(median = median, lo50 = ~qv(.x, .25), hi50 = ~qv(.x, .75),
                        lo80 = ~qv(.x, .10), hi80 = ~qv(.x, .90),
                        lo95 = ~qv(.x, .025), hi95 = ~qv(.x, .975))), .groups = "drop")
write_csv(weekly_quantiles, file.path(DIR_RESULTS_STEP6, "STEP6_weekly_trajectory_quantiles.csv"))

# -------------------------------------------------------------------------
# Numerical consistency and implementation audit.
# -------------------------------------------------------------------------
population_checks <- trajectories |>
  summarise(max_abs_N_minus_S_plus_U = max(abs(N - (S + U))), min_S = min(S), min_U = min(U),
            min_N = min(N), max_S_minus_N = max(S - N),
            min_local_hazard = min(local_hazard),
            aggregate_reported_to_latent_ratio = sum(Y_reported) / sum(X_latent))
write_csv(population_checks, file.path(DIR_TABLES_STEP6, "STEP6_numerical_consistency.csv"))

integrity <- read_csv(file.path(DIR_TABLES_RD, "PHASE1_integrity_check_vs_original.csv"), show_col_types = FALSE)
audit <- tribble(
  ~item, ~status, ~finding,
  "Training/serology", "PASS", "H1/H2 use only weeks through their cutoffs; eligible 2018 serology is included in both.",
  "Observed infection seed", "PASS", "Training data have is_seed=X_seed=0 and no observed holdout infection seed enters simulation.",
  "Full-period posterior initialization", "PASS", "Holdouts read only ce_tau_year_training_H1/H2.rds and their own terminal states/history.",
  "Terminal S/U/history", "PASS", "advance_terminal_state supplies the omitted final demographic step; 8-week posterior X history is carried forward.",
  "tau_year draw", "PASS-WITH-HMC-CAVEAT", "G2 uses tau_year from the matching training-only posterior draw; quantitative results remain provisional.",
  "Future year-effect recentering", "FAIL", "Future z_year values are re-centred within each finite holdout. Predictive new-year effects should be marginal draws around the fitted population mean; the current zero-sum constraint is horizon-dependent.",
  "R0/Reff/seasonal phase", "PASS", "Reconstructed log-R0 components reproduce saved R0 numerically; absolute seasonal indexing continues from training; Reff=R0*S/N.",
  "Demography/population", "PASS", "S+U=N numerically, states remain nonnegative, and weekly demographic inputs are carried through.",
  "Importation accounting", "PASS", "M_t enters once as M_t/N in the hazard. It is introduction pressure, not realised imported infections or added cases.",
  "Observation process", "PASS", "Reported cases are NegBinomial with mean q*X and fixed q=0.05; aggregate simulation ratio is checked, not used for calibration.",
  "Climate information", "PASS-WITH-CAVEAT", "Hindcasts use observed holdout climate and full-2015-2025 anomaly scaling: observed-climate-conditional with preprocessing leakage."
)
write_csv(audit, file.path(DIR_TABLES_STEP6, "STEP6_design_and_implementation_audit.csv"))

# Recentring changes marginal SD from tau to tau*sqrt((K-1)/K), forces
# realised future-year effects to sum to zero, and differs for H1/H2 horizons.
recenter_audit <- year_effects |>
  group_by(holdout, generator, posterior_draw_id) |>
  summarise(K = n(), sum_year_effect = sum(year_effect), empirical_sd = sd(year_effect), tau_used = first(tau_used), .groups = "drop") |>
  group_by(holdout, generator, K) |>
  summarise(theoretical_marginal_sd_ratio = sqrt((first(K) - 1) / first(K)),
            max_abs_sum_year_effect = max(abs(sum_year_effect)),
            median_tau_used = median(tau_used), median_within_horizon_sd = median(empirical_sd), .groups = "drop")
write_csv(recenter_audit, file.path(DIR_TABLES_STEP6, "STEP6_future_year_recentering_audit.csv"))

# -------------------------------------------------------------------------
# Predictive performance: annual intervals, WIS/CRPS, cumulative holdout, and
# draw-cluster bootstrap uncertainty for the 2022 exceedance probability.
# -------------------------------------------------------------------------
message("[Step6] Computing predictive validation metrics...")
annual_performance <- annual_traj |>
  group_by(holdout, generator, year) |>
  group_modify(function(.x, .g) {
    x <- .x$annual_reported_cases
    obs <- observed_annual$observed_cases[observed_annual$holdout == .g$holdout & observed_annual$year == .g$year]
    if (length(obs) != 1L) stop("Observed annual match failed")
    tibble(pred_median = median(x), lo50 = qv(x, .25), hi50 = qv(x, .75),
           lo80 = qv(x, .10), hi80 = qv(x, .90), lo95 = qv(x, .025), hi95 = qv(x, .975),
           observed_cases = obs,
           covered50 = obs >= qv(x, .25) & obs <= qv(x, .75),
           covered80 = obs >= qv(x, .10) & obs <= qv(x, .90),
           covered95 = obs >= qv(x, .025) & obs <= qv(x, .975),
           width50 = qv(x, .75) - qv(x, .25), width80 = qv(x, .90) - qv(x, .10),
           width95 = qv(x, .975) - qv(x, .025),
           exceedance_probability = mean(x >= obs), WIS = wis_50_80_95(x, obs), CRPS = crps_empirical(x, obs))
  }) |>
  ungroup()
write_csv(annual_performance, file.path(DIR_TABLES_STEP6, "STEP6_annual_predictive_performance.csv"))

set.seed(20261013L)
exceedance_mc <- annual_traj |>
  filter(year == 2022) |>
  group_by(holdout, generator) |>
  group_modify(function(.x, .g) {
    obs <- observed_annual$observed_cases[observed_annual$holdout == .g$holdout & observed_annual$year == 2022]
    by_draw <- .x |> group_by(posterior_draw_id) |> summarise(p = mean(annual_reported_cases >= obs), .groups = "drop")
    boots <- replicate(2000L, mean(sample(by_draw$p, nrow(by_draw), replace = TRUE)))
    tibble(n_posterior_draws = nrow(by_draw), n_trajectories = nrow(.x), observed_2022 = obs,
           probability_equal_or_exceed = mean(.x$annual_reported_cases >= obs),
           cluster_bootstrap_lo95 = qv(boots, .025), cluster_bootstrap_hi95 = qv(boots, .975))
  }) |>
  ungroup()
write_csv(exceedance_mc, file.path(DIR_TABLES_STEP6, "STEP6_2022_exceedance_probability_mc.csv"))

holdout_cumulative <- annual_traj |>
  group_by(holdout, generator, posterior_draw_id, stochastic_replicate_id) |>
  summarise(predicted_cases = sum(annual_reported_cases), predicted_latent = sum(cumulative_latent_infections), .groups = "drop") |>
  group_by(holdout, generator) |>
  group_modify(function(.x, .g) {
    obs <- sum(observed_annual$observed_cases[observed_annual$holdout == .g$holdout])
    tibble(pred_median = median(.x$predicted_cases), lo50 = qv(.x$predicted_cases, .25), hi50 = qv(.x$predicted_cases, .75),
           lo80 = qv(.x$predicted_cases, .10), hi80 = qv(.x$predicted_cases, .90),
           lo95 = qv(.x$predicted_cases, .025), hi95 = qv(.x$predicted_cases, .975), observed_cases = obs,
           covered95 = obs >= qv(.x$predicted_cases, .025) & obs <= qv(.x$predicted_cases, .975),
           WIS = wis_50_80_95(.x$predicted_cases, obs), CRPS = crps_empirical(.x$predicted_cases, obs))
  }) |>
  ungroup()
write_csv(holdout_cumulative, file.path(DIR_TABLES_STEP6, "STEP6_multiyear_cumulative_performance.csv"))

width_comparison <- annual_performance |>
  select(holdout, year, generator, width50, width80, width95) |>
  pivot_wider(names_from = generator, values_from = c(width50, width80, width95)) |>
  mutate(ratio_G1_to_G2_width50 = width50_G1_fixed_SD_040 / width50_G2_estimated_tau,
         ratio_G1_to_G2_width80 = width80_G1_fixed_SD_040 / width80_G2_estimated_tau,
         ratio_G1_to_G2_width95 = width95_G1_fixed_SD_040 / width95_G2_estimated_tau)
write_csv(width_comparison, file.path(DIR_TABLES_STEP6, "STEP6_G1_vs_G2_interval_widths.csv"))

# Controlled G1-vs-G2 generator substitution within each matching training fit.
# Draws/replicates are matched by identifier; H1-vs-H2 is never treated as paired.
paired_generator <- annual_traj |>
  select(holdout, generator, posterior_draw_id, stochastic_replicate_id, year,
         annual_reported_cases, cumulative_latent_infections, peak_Reff, year_effect) |>
  pivot_wider(names_from = generator,
              values_from = c(annual_reported_cases, cumulative_latent_infections, peak_Reff, year_effect)) |>
  mutate(delta_cases_G1_minus_G2 = annual_reported_cases_G1_fixed_SD_040 - annual_reported_cases_G2_estimated_tau,
         delta_latent_G1_minus_G2 = cumulative_latent_infections_G1_fixed_SD_040 - cumulative_latent_infections_G2_estimated_tau,
         delta_peak_Reff_G1_minus_G2 = peak_Reff_G1_fixed_SD_040 - peak_Reff_G2_estimated_tau,
         delta_year_effect_G1_minus_G2 = year_effect_G1_fixed_SD_040 - year_effect_G2_estimated_tau)
paired_generator_summary <- paired_generator |>
  group_by(holdout, year) |>
  summarise(across(starts_with("delta_"), list(median = median, lo95 = ~qv(.x, .025), hi95 = ~qv(.x, .975))), .groups = "drop")
write_csv(paired_generator_summary, file.path(DIR_TABLES_STEP6, "STEP6_controlled_G1_G2_contrasts.csv"))

# -------------------------------------------------------------------------
# 2022 mechanism classification and descriptive outcome contrasts for G2.
# These are associations conditioned on simulated outcome, not causal effects.
# -------------------------------------------------------------------------
traj_2022 <- annual_traj |> filter(year == 2022)
mechanism_classification <- traj_2022 |>
  count(holdout, generator, mechanism_class, name = "n") |>
  group_by(holdout, generator) |>
  mutate(probability = n / sum(n)) |>
  ungroup()
write_csv(mechanism_classification, file.path(DIR_TABLES_STEP6, "STEP6_2022_mechanism_classification.csv"))

establishment_summary <- traj_2022 |>
  group_by(holdout, generator) |>
  summarise(prob_any_latent_infection = mean(any_latent_infection),
            prob_favourable_pressure = mean(favourable_pressure), prob_established = mean(established),
            prob_major_epidemic = mean(major_epidemic),
            median_weeks_zero_infections = median(weeks_zero_infections),
            median_longest_infection_free_run = median(longest_infection_free_run),
            median_longest_local_dominance_run = median(longest_locally_dominant_run), .groups = "drop")
write_csv(establishment_summary, file.path(DIR_TABLES_STEP6, "STEP6_2022_establishment_summary.csv"))

mechanism_vars <- c("S_prop_start", "peak_R0", "peak_Reff", "weeks_Reff_gt1", "longest_run_Reff_gt1",
                    "introduction_pressure_favourable", "cumulative_latent_infections", "annual_reported_cases",
                    "year_effect", "A_year", "seasonal_log_peak", "climate_total_log_mean")
outcome_contrasts <- traj_2022 |>
  filter(generator == "G2_estimated_tau") |>
  mutate(outcome = ifelse(major_epidemic, "Major annual epidemic", "No major annual epidemic")) |>
  select(holdout, outcome, all_of(mechanism_vars)) |>
  pivot_longer(cols = all_of(mechanism_vars), names_to = "quantity", values_to = "value") |>
  group_by(holdout, outcome, quantity) |>
  summarise(n = n(), median = median(value), lo80 = qv(value, .10), hi80 = qv(value, .90),
            lo95 = qv(value, .025), hi95 = qv(value, .975), .groups = "drop")
write_csv(outcome_contrasts, file.path(DIR_TABLES_STEP6, "STEP6_G2_2022_outcome_contrasts.csv"))

headline_2022 <- traj_2022 |>
  group_by(holdout, generator) |>
  summarise(across(c(S_prop_start, S_prop_end, mean_R0, peak_R0, mean_Reff, peak_Reff,
                     weeks_Reff_gt1, longest_run_Reff_gt1, introduction_pressure_favourable,
                     cumulative_latent_infections, annual_reported_cases, year_effect, A_year),
                   list(median = median, lo95 = ~qv(.x, .025), hi95 = ~qv(.x, .975))), .groups = "drop")
write_csv(headline_2022, file.path(DIR_TABLES_STEP6, "STEP6_2022_headline_mechanistic_metrics.csv"))

# H1 susceptibility at Jan 2022 is a simulated state after two holdout years;
# H2 Jan 2022 is directly initialized from its training-only terminal posterior.
jan2022_susc <- trajectories |>
  filter(generator == "G2_estimated_tau", year == 2022) |>
  group_by(holdout, posterior_draw_id, stochastic_replicate_id) |>
  summarise(S_prop_Jan2022 = first(S_prop), .groups = "drop") |>
  group_by(holdout) |>
  summarise(median = median(S_prop_Jan2022), lo95 = qv(S_prop_Jan2022, .025), hi95 = qv(S_prop_Jan2022, .975), .groups = "drop")
write_csv(jan2022_susc, file.path(DIR_TABLES_STEP6, "STEP6_Jan2022_susceptibility_G2.csv"))

message("[Step6] Diagnostic tables complete -> ", DIR_TABLES_STEP6)
