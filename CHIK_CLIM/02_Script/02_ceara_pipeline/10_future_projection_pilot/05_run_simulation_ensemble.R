# STEP 2 -- Stochastic future epidemic generator: simulation ensemble runner.
#
# Produces, for each requested scenario (year-effect arm x climate arm x
# importation level x demography arm), an ensemble of posterior-draw x
# stochastic-replicate trajectories over 2026-2035. Run 04_run_unit_tests.R
# FIRST -- this script assumes the engine already passed validation.
#
# Illustrative scenario simulation; no SSP forcing or vaccination applied.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot/00_config_future_pilot.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/10_future_projection_pilot/03_future_renewal_functions.R"))

message("[Step2 ensemble] Loading accepted climate-forced fit (read-only)...")
climate <- readRDS(ACCEPTED_CLIMATE_FIT)
sd_clim <- climate$stan_data
dates_hist <- as.Date(climate$weekly_data$week_start)
N_hist <- length(dates_hist)
G <- sd_clim$G
w <- generation_weights_future()

draws <- rstan::extract(climate$fit, pars = c(
  "alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year",
  "beta_T", "beta_P", "phi_obs", "S", "U", "X"
), permuted = TRUE)
n_draws_available <- length(draws$alpha_R)

future_dates <- build_future_dates()
n_future_weeks <- length(future_dates)
t_index_future <- (N_hist + 1):(N_hist + n_future_weeks)
years_future <- as.integer(format(future_dates, "%Y"))
future_year_labels <- sort(unique(years_future))
n_future_years <- length(future_year_labels)
year_id_future <- match(years_future, future_year_labels)
weeks_per_future_year <- as.integer(table(factor(years_future, levels = future_year_labels)))

# ---- Historical climate-by-year lookup, for C1/C2 whole-year block resampling ----
hist_years <- as.integer(format(dates_hist, "%Y"))
hist_climate_by_year <- split(tibble(z_T_anom = sd_clim$z_T_anom, z_P_anom = sd_clim$z_P_anom), hist_years)

# ---- Future demography arms (built ONCE -- not posterior-draw-dependent) ----
reconciliation_last_hist <- sd_clim$N_end[N_hist] - sd_clim$N_start[N_hist] - sd_clim$births[N_hist] + sd_clim$deaths[N_hist]
N_end_anchor <- sd_clim$N_end[N_hist]
ibge_uf <- readRDS(ACCEPTED_IBGE_UF_PROJ) |> dplyr::filter(uf_code == UF_CODE_CEARA, year >= min(future_year_labels) - 1, year <= max(future_year_labels))
demography_D1 <- build_future_demography("D1_ibge_projection", future_dates, N_end_anchor, ibge_uf_annual = ibge_uf)
demography_D0 <- build_future_demography("D0_stationary_2025_levels", future_dates, N_end_anchor, hist_weekly = climate$weekly_data |> tail(52))
demography_tables <- list(D1_ibge_projection = demography_D1, D0_stationary_2025_levels = demography_D0)
write_csv(demography_D1, file.path(DIR_TABLES, "STEP2_future_demography_D1_ibge_projection.csv"))
write_csv(demography_D0, file.path(DIR_TABLES, "STEP2_future_demography_D0_stationary.csv"))
message(sprintf("[Step2 ensemble] Future population: D1(IBGE-rate) 2035 end = %.0f (vs. 2026 start %.0f); D0(stationary) 2035 end = %.0f",
                 tail(demography_D1$N_end, 1), demography_D1$N_start[1], tail(demography_D0$N_end, 1)))

# ================================================================
# Deterministic sub-seed derivation: stable, reproducible, and distinct
# per (scenario, draw-slot, climate-sequence, year-effect, replicate).
# ================================================================
sub_seed <- function(...) {
  key <- paste(..., sep = "|")
  (BASE_SEED + sum(utf8ToInt(key)) * 2654435761) %% .Machine$integer.max
}

#' Run one fully-specified scenario: n_draws posterior draws x n_replicates
#' stochastic replicates each. Returns a list(trajectories, summary).
run_scenario <- function(scenario_id, year_arm, climate_arm, lambda_import, demography_arm, n_draws, n_replicates) {
  started <- Sys.time()
  demog <- demography_tables[[demography_arm]]
  set.seed(sub_seed("drawselect", scenario_id))
  draw_ids <- sample.int(n_draws_available, n_draws)
  climate_build_arm <- if (climate_arm == "C2_resampled_weather_zero_beta") "C1_resampled_weather_variability" else climate_arm
  beta_forced_zero <- climate_arm == "C2_resampled_weather_zero_beta"

  traj_list <- vector("list", n_draws * n_replicates)
  wave_list <- vector("list", n_draws * n_replicates)
  slot <- 0L
  for (di in seq_along(draw_ids)) {
    d <- draw_ids[di]
    adv <- advance_terminal_state(draws$S[d, N_hist], draws$U[d, N_hist], draws$X[d, N_hist], sd_clim$births[N_hist], sd_clim$deaths[N_hist], reconciliation_last_hist)
    X_recent <- draws$X[d, (N_hist - G + 1):N_hist]

    set.seed(sub_seed("yeareffect", scenario_id, d))
    year_struct <- sample_future_year_structure(year_arm, n_future_years, year_effect_prior_sd = 0.40, sigma_season_year_draw = draws$sigma_season_year[d], rng = function(n) rnorm(n))

    set.seed(sub_seed("climateblock", scenario_id, d))
    climate_block <- sample_future_climate_block(climate_build_arm, future_year_labels, weeks_per_future_year, hist_climate_by_year,
                                                  rng_discrete = function(n, choices) sample(choices, n, replace = TRUE))

    beta_T_eff <- if (beta_forced_zero) 0 else draws$beta_T[d]
    beta_P_eff <- if (beta_forced_zero) 0 else draws$beta_P[d]
    r0_future <- compute_R0_future(t_index_future, year_id_future, draws$alpha_R[d], year_struct$year_effect, year_struct$A_year,
                                    draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d],
                                    beta_T_eff, beta_P_eff, climate_block$z_T, climate_block$z_P)

    for (r in seq_len(n_replicates)) {
      slot <- slot + 1L
      set.seed(sub_seed("replicate", scenario_id, d, r))
      sim <- simulate_renewal_core(
        R0_vec = r0_future$R0, w = w, S0 = adv$S, U0 = adv$U,
        N_start_vec = demog$N_start, N_end_vec = demog$N_end, births_vec = demog$births, deaths_vec = demog$deaths,
        X_recent_history = X_recent, lambda_import = lambda_import, stochastic = TRUE,
        rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p)
      )
      Y_rep <- simulate_reported_cases(sim$X, q = sd_clim$q, phi_obs = draws$phi_obs[d], stochastic = TRUE,
                                        rng_nb = function(n, mu, phi) rnbinom(n, mu = mu, size = phi))
      traj_list[[slot]] <- tibble(
        scenario_id = scenario_id, posterior_draw_id = d, climate_sequence_id = paste(climate_block$source_years, collapse = ";"),
        year_effect_id = sub_seed("yeareffect", scenario_id, d), importation_scenario_id = lambda_import, stochastic_replicate_id = r,
        week_start = future_dates, S = sim$S, U = sim$U, N = sim$N, X_latent = sim$X, Y_reported = Y_rep,
        R0 = sim$R0, R_eff = sim$R_eff, lambda = sim$lambda, M_import = sim$M
      )
      waves <- detect_simulated_waves(future_dates, Y_rep, sim$N)
      if (nrow(waves)) wave_list[[slot]] <- waves |> mutate(scenario_id = scenario_id, posterior_draw_id = d, stochastic_replicate_id = r, .before = 1)
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  message(sprintf("[Step2 ensemble] %-45s %3d draws x %2d reps = %5d trajectories | %.1fs", scenario_id, n_draws, n_replicates, n_draws * n_replicates, elapsed))
  list(trajectories = bind_rows(traj_list), waves = bind_rows(wave_list))
}

#' Weekly -> annual / cumulative / S-proportion summary for one scenario's trajectories.
summarise_scenario <- function(traj) {
  traj <- traj |> mutate(year = as.integer(format(week_start, "%Y")))
  annual <- traj |> group_by(scenario_id, posterior_draw_id, stochastic_replicate_id, year) |>
    summarise(annual_latent_infections = sum(X_latent), annual_reported_cases = sum(Y_reported), mean_S_over_N = mean(S / N), .groups = "drop")
  cumulative <- traj |> group_by(scenario_id, posterior_draw_id, stochastic_replicate_id) |>
    summarise(
      cum_5y = sum(X_latent[year <= FUTURE_YEAR_START + 4]), cum_10y = sum(X_latent),
      incidence_per_100k_5y = 1e5 * sum(X_latent[year <= FUTURE_YEAR_START + 4]) / mean(N[year <= FUTURE_YEAR_START + 4]),
      incidence_per_100k_10y = 1e5 * sum(X_latent) / mean(N), final_S_over_N = S[which.max(week_start)] / N[which.max(week_start)],
      .groups = "drop"
    )
  list(annual = annual, cumulative = cumulative)
}

# ================================================================
# Scenario definitions
# ================================================================
ENSEMBLE_MODE <- Sys.getenv("FUTURE_PILOT_ENSEMBLE_MODE", "smoke") # "smoke" or "pilot"
settings_smoke <- ENSEMBLE_SETTINGS$smoke
settings_pilot <- ENSEMBLE_SETTINGS$pilot

reference_scenario <- list(id = "REF_Y1_C1_D1_import1", year_arm = "Y1_resampled_distributional",
                            climate_arm = "C1_resampled_weather_variability", lambda_import = IMPORTATION_REFERENCE, demography_arm = "D1_ibge_projection")

sensitivity_grid <- tidyr::expand_grid(
  year_arm = YEAR_EFFECT_ARMS, climate_arm = CLIMATE_ARMS, lambda_import = IMPORTATION_GRID
) |> mutate(demography_arm = "D1_ibge_projection",
            id = sprintf("SENS_%s_%s_import%s", sub("_.*", "", year_arm), sub("_.*", "", climate_arm), lambda_import))
demography_sensitivity <- tibble(id = c("DEMOG_D1", "DEMOG_D0"), year_arm = "Y1_resampled_distributional",
                                  climate_arm = "C1_resampled_weather_variability", lambda_import = IMPORTATION_REFERENCE,
                                  demography_arm = c("D1_ibge_projection", "D0_stationary_2025_levels"))

message(sprintf("[Step2 ensemble] MODE=%s | reference scenario at %s scale; %d sensitivity scenarios + 2 demography-arm scenarios at smoke scale",
                 ENSEMBLE_MODE, ENSEMBLE_MODE, nrow(sensitivity_grid)))

# ---- Run the reference scenario at the requested ensemble scale ----
ref_settings <- if (ENSEMBLE_MODE == "pilot") settings_pilot else settings_smoke
ref_result <- run_scenario(reference_scenario$id, reference_scenario$year_arm, reference_scenario$climate_arm,
                            reference_scenario$lambda_import, reference_scenario$demography_arm, ref_settings$n_draws, ref_settings$n_replicates)
saveRDS(ref_result$trajectories, file.path(DIR_RESULTS, "STEP2_reference_scenario_trajectories.rds"))
ref_summary <- summarise_scenario(ref_result$trajectories)
write_csv(ref_summary$annual, file.path(DIR_RESULTS, "STEP2_reference_scenario_annual_summary.csv"))
write_csv(ref_summary$cumulative, file.path(DIR_RESULTS, "STEP2_reference_scenario_cumulative_summary.csv"))
write_csv(ref_result$waves, file.path(DIR_RESULTS, "STEP2_reference_scenario_waves.csv"))

# ---- Run the structural sensitivity grid + demography comparison, always at smoke scale ----
all_sens <- bind_rows(sensitivity_grid, demography_sensitivity)
sens_cumulative_list <- vector("list", nrow(all_sens))
for (i in seq_len(nrow(all_sens))) {
  s <- all_sens[i, ]
  res <- run_scenario(s$id, s$year_arm, s$climate_arm, s$lambda_import, s$demography_arm, settings_smoke$n_draws, settings_smoke$n_replicates)
  summ <- summarise_scenario(res$trajectories)
  sens_cumulative_list[[i]] <- summ$cumulative |> mutate(year_arm = s$year_arm, climate_arm = s$climate_arm, lambda_import = s$lambda_import, demography_arm = s$demography_arm)
}
sensitivity_cumulative <- bind_rows(sens_cumulative_list)
write_csv(sensitivity_cumulative, file.path(DIR_RESULTS, "STEP2_sensitivity_grid_cumulative_summary.csv"))

message("\n[Step2 ensemble] Complete. Reference scenario + ", nrow(all_sens), " sensitivity scenarios written to ", DIR_RESULTS)
