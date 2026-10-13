# PHASE 2 -- Posterior-predictive future year effects: G0 vs. G1 vs. G2.
#
# G0: Neutral (year_effect=0, A_year=1) -- diagnostic floor only.
# G1: Existing fixed-SD=0.40 year_effect, WITH the resampled seasonal-
#     amplitude hierarchy retained (A_year varies) -- this is the
#     original Y1 generator's year-effect mechanism exactly.
# G2: Posterior-predictive year_effect using the SAME draw's ESTIMATED
#     tau_year (from the Phase 1 fit using prior P2 by default), with
#     the SAME seasonal-amplitude mechanism retained.
#
# PAIRED design: G0/G1/G2 use the SAME posterior draw d (alpha_R, beta_sin1,
# ..., sigma_season_year, beta_T, beta_P, terminal state ALL identical
# across the three arms for a given d) and the SAME fresh future z_year/
# z_season innovations -- only the year-effect SCALE differs (0 / fixed
# 0.40 / that draw's own tau_year[d]). This isolates the scale-choice
# question from ordinary posterior and Monte Carlo noise, exactly as
# Phase 2's instructions require ("common stochastic inputs where a valid
# paired comparison is possible").
#
# Historical sum-to-zero vs. future generation: the historical
# `year_effect = tau_year*(z_year-mean(z_year))` identification uses
# EXACTLY the Y=11 FITTED years' z_year to compute `mean(z_year)` --
# that mean is a property of the fitted past, not re-estimated here. For
# FUTURE years there is no such fitted set: each future year's deviation
# is a FRESH draw from the same N(0,1) data-generating assumption the
# model places on z_year, demeaned only across THAT future block (10
# future years), analogous to -- but not numerically tied to -- the
# historical demeaning. This mirrors exactly how 03_future_renewal_functions.R's
# `sample_future_year_structure()` already treats Y1 (reused unmodified
# below); it is NOT automatically re-centred against a combined
# historical+future set, and the future demeaning window is the
# arbitrary 10-year projection horizon chosen for this exercise, not a
# principled infinite-population mean -- disclosed explicitly here and
# in the final report.
#
# alpha_R consistency: alpha_R is the SAME fitted intercept in every arm
# (drawn from the SAME posterior draw), so it continues to represent
# "typical" log R0 net of the (now-zero-mean-within-window, by
# construction) year effect -- its interpretation does not change
# between arms because only the year_effect TERM varies, never alpha_R itself.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/12_future_year_effect_validation/00_config_year_effect_validation.R"))

TAU_PRIOR_FOR_G2 <- "P2" # default; see PHASE1_tau_year_posterior_summary.csv -- revisit if Phase 1 strongly favours another prior
tau_fit_path <- file.path(DIR_FITS_YV, sprintf("ce_tau_year_%s.rds", TAU_PRIOR_FOR_G2))
if (!file.exists(tau_fit_path)) stop("Missing Phase 1 tau_year fit: ", tau_fit_path, " -- run 02_phase1_fit_tau_year.R first.")

message("[Phase2] Loading the tau_year (", TAU_PRIOR_FOR_G2, ") full-period fit (read-only)...")
tau_bundle <- readRDS(tau_fit_path)
sd_tau <- tau_bundle$stan_data
dates_hist <- as.Date(tau_bundle$weekly_data$week_start)
N_hist <- length(dates_hist)
G <- sd_tau$G
w <- generation_weights_future()

draws <- rstan::extract(tau_bundle$fit, pars = c(
  "alpha_R", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2", "sigma_season_year", "tau_year",
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
z_zero <- rep(0, n_future_weeks) # C0 climate-neutral throughout -- isolates year-effect mechanisms, consistent with the earlier year_effect_diagnostic design

reconciliation_last_hist <- sd_tau$N_end[N_hist] - sd_tau$N_start[N_hist] - sd_tau$births[N_hist] + sd_tau$deaths[N_hist]
N_end_anchor <- sd_tau$N_end[N_hist]
ibge_uf <- readRDS(ACCEPTED_IBGE_UF_PROJ) |> dplyr::filter(uf_code == UF_CODE_CEARA, year >= min(future_year_labels) - 1, year <= max(future_year_labels))
demog <- build_future_demography("D1_ibge_projection", future_dates, N_end_anchor, ibge_uf_annual = ibge_uf)

GENERATORS <- c("G0_neutral", "G1_fixed_SD_040", "G2_estimated_tau")
GENERATOR_LABELS <- c(G0_neutral = "G0: neutral (diagnostic)", G1_fixed_SD_040 = "G1: fixed SD=0.40 (existing Y1)", G2_estimated_tau = "G2: estimated tau_year")

build_three_generators <- function(z_year_f, z_season_f, sigma_season_year_draw, tau_year_draw) {
  A_year_varying <- exp(sigma_season_year_draw * z_season_f) # SAME amplitude draw reused across G0/G1/G2 below (G0 overrides to 1)
  list(
    G0_neutral        = list(year_effect = rep(0, n_future_years), A_year = rep(1, n_future_years)),
    G1_fixed_SD_040   = list(year_effect = 0.40 * (z_year_f - mean(z_year_f)), A_year = A_year_varying),
    G2_estimated_tau  = list(year_effect = tau_year_draw * (z_year_f - mean(z_year_f)), A_year = A_year_varying)
  )
}

run_generator_comparison <- function(n_draws, n_replicates, label) {
  started <- Sys.time()
  set.seed(BASE_SEED_YV)
  draw_ids <- sample.int(n_draws_available, n_draws)

  traj_list <- vector("list", n_draws * n_replicates * length(GENERATORS)); slot <- 0L
  for (di in seq_along(draw_ids)) {
    d <- draw_ids[di]
    adv <- advance_terminal_state(draws$S[d, N_hist], draws$U[d, N_hist], draws$X[d, N_hist], sd_tau$births[N_hist], sd_tau$deaths[N_hist], reconciliation_last_hist)
    X_recent <- draws$X[d, (N_hist - G + 1):N_hist]
    set.seed(BASE_SEED_YV + d)
    z_year_f <- rnorm(n_future_years); z_season_f <- rnorm(n_future_years)
    gens <- build_three_generators(z_year_f, z_season_f, draws$sigma_season_year[d], draws$tau_year[d])
    r0_by_gen <- lapply(GENERATORS, function(g_name) {
      a <- gens[[g_name]]
      compute_R0_future(t_index_future, year_id_future, draws$alpha_R[d], a$year_effect, a$A_year,
                         draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d],
                         draws$beta_T[d], draws$beta_P[d], z_zero, z_zero)
    }); names(r0_by_gen) <- GENERATORS
    for (r in seq_len(n_replicates)) {
      common_seed <- BASE_SEED_YV + d * 1000L + r
      for (g_name in GENERATORS) {
        slot <- slot + 1L
        set.seed(common_seed)
        sim <- simulate_renewal_core(R0_vec = r0_by_gen[[g_name]]$R0, w = w, S0 = adv$S, U0 = adv$U,
                                      N_start_vec = demog$N_start, N_end_vec = demog$N_end, births_vec = demog$births, deaths_vec = demog$deaths,
                                      X_recent_history = X_recent, lambda_import = IMPORTATION_REFERENCE, stochastic = TRUE,
                                      rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p))
        Y_rep <- simulate_reported_cases(sim$X, q = QVAL, phi_obs = draws$phi_obs[d], stochastic = TRUE, rng_nb = function(n, mu, phi) rnbinom(n, mu = mu, size = phi))
        traj_list[[slot]] <- tibble(generator = g_name, posterior_draw_id = d, stochastic_replicate_id = r, week_start = future_dates,
                                     S = sim$S, N = sim$N, X_latent = sim$X, Y_reported = Y_rep, R0 = sim$R0, R_eff = sim$R_eff)
      }
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  message(sprintf("[Phase2] %-6s %3d draws x %2d reps x %d generators = %6d trajectories | %.1fs", label, n_draws, n_replicates, length(GENERATORS), n_draws * n_replicates * length(GENERATORS), elapsed))
  bind_rows(traj_list)
}

summarise_generator_trajectories <- function(traj) {
  traj <- traj |> mutate(year = as.integer(format(week_start, "%Y")))
  annual_R0 <- traj |> group_by(generator, posterior_draw_id, year) |> summarise(mean_R0 = mean(R0), peak_R0 = max(R0), .groups = "drop")
  per_traj <- traj |> group_by(generator, posterior_draw_id, stochastic_replicate_id) |>
    summarise(cum_10y = sum(X_latent), incidence_per_100k_10y = 1e5 * sum(X_latent) / mean(N),
              frac_weeks_Reff_gt1 = mean(R_eff > 1), n_Reff_episodes = { rr <- rle(R_eff > 1); sum(rr$values, na.rm = TRUE) },
              min_S_over_N = min(S / N), final_S_over_N = S[which.max(week_start)] / N[which.max(week_start)], .groups = "drop")
  major <- traj |> group_by(generator, posterior_draw_id, stochastic_replicate_id) |>
    group_modify(~ { waves <- detect_simulated_waves(.x$week_start, .x$Y_reported, .x$N); tibble(has_major_epidemic = if (nrow(waves)) any(waves$major_epidemic) else FALSE) }) |> ungroup()
  list(annual_R0 = annual_R0, per_trajectory = per_traj |> left_join(major, by = c("generator", "posterior_draw_id", "stochastic_replicate_id")))
}

message("[Phase2] Smoke test...")
smoke_traj <- run_generator_comparison(ENSEMBLE_RETRO$smoke$n_draws, ENSEMBLE_RETRO$smoke$n_replicates, "SMOKE")
smoke_summary <- summarise_generator_trajectories(smoke_traj)
print(as.data.frame(smoke_summary$per_trajectory |> group_by(generator) |> summarise(median_cum10y = median(cum_10y), .groups = "drop")))

message("\n[Phase2] Full-scale paired run...")
full_traj <- run_generator_comparison(ENSEMBLE_RETRO$full$n_draws, ENSEMBLE_RETRO$full$n_replicates, "FULL")
saveRDS(full_traj, file.path(DIR_RESULTS_YV, "PHASE2_generator_trajectories.rds"))

# Individual future year_effect realisations for G1/G2 (for FIG18 -- NOT
# derivable from the trajectory summaries above, which only keep R0(t)/X(t)).
message("[Phase2] Saving individual future year_effect realisations (G1/G2) for FIG18...")
set.seed(BASE_SEED_YV)
draw_ids_ye <- sample.int(n_draws_available, ENSEMBLE_RETRO$full$n_draws)
future_ye <- bind_rows(lapply(draw_ids_ye, function(d) {
  set.seed(BASE_SEED_YV + d)
  z_year_f <- rnorm(n_future_years); invisible(rnorm(n_future_years)) # consume z_season_f identically to build_three_generators()'s own stream, for consistency (not otherwise used here)
  bind_rows(
    tibble(generator = "G1_fixed_SD_040", draw = d, future_year_index = seq_len(n_future_years), year_effect = 0.40 * (z_year_f - mean(z_year_f))),
    tibble(generator = "G2_estimated_tau", draw = d, future_year_index = seq_len(n_future_years), year_effect = draws$tau_year[d] * (z_year_f - mean(z_year_f)))
  )
}))
write_csv(future_ye, file.path(DIR_RESULTS_YV, "PHASE2_future_year_effect_samples.csv"))
full_summary <- summarise_generator_trajectories(full_traj)
write_csv(full_summary$per_trajectory, file.path(DIR_RESULTS_YV, "PHASE2_per_trajectory_summary.csv"))
write_csv(full_summary$annual_R0, file.path(DIR_RESULTS_YV, "PHASE2_annual_R0_by_generator.csv"))

gen_summary <- full_summary$per_trajectory |> group_by(generator) |>
  summarise(n = n(), median_cum10y = median(cum_10y), mean_cum10y = mean(cum_10y),
            median_min_S_over_N = median(min_S_over_N), proportion_major_epidemic = mean(has_major_epidemic),
            median_frac_weeks_Reff_gt1 = median(frac_weeks_Reff_gt1), .groups = "drop")
message("\n=== Phase2 generator summary (full scale) ===")
print(as.data.frame(gen_summary), digits = 4)
write_csv(gen_summary, file.path(DIR_TABLES_YV, "PHASE2_generator_summary.csv"))

# Extreme annual peak R0 frequency, matched to the climate_robustness/year_effect_diagnostic convention
peak_extreme <- full_summary$annual_R0 |> group_by(generator) |>
  summarise(median_peak_R0 = median(peak_R0), p90_peak_R0 = quantile(peak_R0, .90), p99_peak_R0 = quantile(peak_R0, .99), max_peak_R0 = max(peak_R0), .groups = "drop")
print(as.data.frame(peak_extreme), digits = 4)
write_csv(peak_extreme, file.path(DIR_TABLES_YV, "PHASE2_peak_R0_extremes_by_generator.csv"))

# Paired contrasts G2-G1 and G1-G0 on cum_10y (same draw x replicate)
wide <- full_summary$per_trajectory |> select(posterior_draw_id, stochastic_replicate_id, generator, cum_10y) |>
  pivot_wider(names_from = generator, values_from = cum_10y)
paired <- wide |> mutate(diff_G1_minus_G0 = G1_fixed_SD_040 - G0_neutral, diff_G2_minus_G0 = G2_estimated_tau - G0_neutral, diff_G2_minus_G1 = G2_estimated_tau - G1_fixed_SD_040)
paired_summary <- tibble(
  contrast = c("G1 - G0 (fixed-SD year effect vs. neutral)", "G2 - G0 (estimated-tau year effect vs. neutral)", "G2 - G1 (estimated-tau vs. fixed-SD)"),
  median_diff = c(median(paired$diff_G1_minus_G0), median(paired$diff_G2_minus_G0), median(paired$diff_G2_minus_G1)),
  lo95 = c(quantile(paired$diff_G1_minus_G0, .025), quantile(paired$diff_G2_minus_G0, .025), quantile(paired$diff_G2_minus_G1, .025)),
  hi95 = c(quantile(paired$diff_G1_minus_G0, .975), quantile(paired$diff_G2_minus_G0, .975), quantile(paired$diff_G2_minus_G1, .975))
)
print(as.data.frame(paired_summary), digits = 4)
write_csv(paired_summary, file.path(DIR_TABLES_YV, "PHASE2_paired_contrasts.csv"))
write_csv(paired, file.path(DIR_RESULTS_YV, "PHASE2_paired_differences_cum10y.csv"))

message("\n[Phase2] Generator comparison complete. Tables -> ", DIR_TABLES_YV)
