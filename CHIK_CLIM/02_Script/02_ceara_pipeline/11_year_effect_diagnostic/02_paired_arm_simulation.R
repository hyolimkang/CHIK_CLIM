# Year-effect diagnostic -- Part 3/4: paired 4-arm decomposition.
#
# FIXES the confound in 05_run_simulation_ensemble.R's run_scenario():
# there, draw_ids and the future year-effect innovations are both seeded
# by `sub_seed(..., scenario_id)`, so two different scenarios (e.g. Y0 vs
# Y1) use DIFFERENT posterior draws AND different stochastic innovations
# -- any outcome contrast is confounded with pure Monte Carlo sampling
# noise, not just the structural arm difference.
#
# This script instead uses a PAIRED design: for each (posterior draw,
# stochastic replicate), the SAME draw, the SAME future z_year/z_season
# innovations, the SAME climate input (C0, zero, per instruction), the
# SAME demography, and the SAME starting RNG seed are used across all 4
# arms -- only the year_effect/A_year CONSTRUCTION differs between arms.
# This isolates the structural contrast from sampling noise.
#
# Does NOT modify 03_future_renewal_functions.R -- reuses its functions
# unchanged (simulate_renewal_core, compute_R0_future, advance_terminal_state,
# build_future_demography, simulate_reported_cases, detect_simulated_waves).

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/11_year_effect_diagnostic/00_config_year_effect_diagnostic.R"))

message("[Diag 2] Loading accepted climate-forced fit (read-only; no refit)...")
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
z_zero <- rep(0, n_future_weeks) # C0 climate-neutral, shared across all arms

reconciliation_last_hist <- sd_clim$N_end[N_hist] - sd_clim$N_start[N_hist] - sd_clim$births[N_hist] + sd_clim$deaths[N_hist]
N_end_anchor <- sd_clim$N_end[N_hist]
ibge_uf <- readRDS(ACCEPTED_IBGE_UF_PROJ) |> dplyr::filter(uf_code == UF_CODE_CEARA, year >= min(future_year_labels) - 1, year <= max(future_year_labels))
demog <- build_future_demography("D1_ibge_projection", future_dates, N_end_anchor, ibge_uf_annual = ibge_uf) # single, shared demography arm (instruction: hold common)

YEAR_EFFECT_PRIOR_SD <- 0.40

#' Construct the 4 arms' (year_effect, A_year) from ONE shared pair of
#' future innovations -- this is what makes the design paired.
build_four_arms <- function(z_year_f, z_season_f, sigma_season_year_draw) {
  year_effect_y1 <- YEAR_EFFECT_PRIOR_SD * (z_year_f - mean(z_year_f))
  A_year_y1 <- exp(sigma_season_year_draw * z_season_f)
  list(
    A_neutral        = list(year_effect = rep(0, n_future_years), A_year = rep(1, n_future_years)),
    B_year_only      = list(year_effect = year_effect_y1,          A_year = rep(1, n_future_years)),
    C_amplitude_only = list(year_effect = rep(0, n_future_years),  A_year = A_year_y1),
    D_both           = list(year_effect = year_effect_y1,          A_year = A_year_y1)
  )
}

#' Run the paired design for n_draws x n_replicates. Returns one long
#' tibble: (draw, replicate, arm, weekly trajectory columns).
run_paired_diagnostic <- function(n_draws, n_replicates, label) {
  started <- Sys.time()
  set.seed(DIAG_BASE_SEED) # ONE common draw-selection seed -- shared by construction (not per-arm)
  draw_ids <- sample.int(n_draws_available, n_draws)

  traj_list <- vector("list", n_draws * n_replicates * length(DECOMP_ARMS))
  slot <- 0L
  for (di in seq_along(draw_ids)) {
    d <- draw_ids[di]
    adv <- advance_terminal_state(draws$S[d, N_hist], draws$U[d, N_hist], draws$X[d, N_hist], sd_clim$births[N_hist], sd_clim$deaths[N_hist], reconciliation_last_hist)
    X_recent <- draws$X[d, (N_hist - G + 1):N_hist]

    set.seed(DIAG_BASE_SEED + d) # ONE common future-innovation seed per draw -- shared across all 4 arms
    z_year_f <- rnorm(n_future_years); z_season_f <- rnorm(n_future_years)
    arms <- build_four_arms(z_year_f, z_season_f, draws$sigma_season_year[d])

    r0_by_arm <- lapply(DECOMP_ARMS, function(arm_name) {
      a <- arms[[arm_name]]
      compute_R0_future(t_index_future, year_id_future, draws$alpha_R[d], a$year_effect, a$A_year,
                         draws$beta_sin1[d], draws$beta_cos1[d], draws$beta_sin2[d], draws$beta_cos2[d],
                         draws$beta_T[d], draws$beta_P[d], z_zero, z_zero) # beta coefficients irrelevant since z=0 (C0); kept for consistency with compute_R0_future's signature
    })
    names(r0_by_arm) <- DECOMP_ARMS

    for (r in seq_len(n_replicates)) {
      common_seed <- DIAG_BASE_SEED + d * 1000L + r # ONE common starting seed per (draw, replicate) -- reset before EACH arm below (common random numbers, "where feasible": true for the first random draws consumed each week; decoheres after arms' trajectories diverge enough to change rbinom's internal algorithm path -- disclosed in the report)
      for (arm_name in DECOMP_ARMS) {
        slot <- slot + 1L
        set.seed(common_seed) # reset to the SAME seed for every arm
        sim <- simulate_renewal_core(
          R0_vec = r0_by_arm[[arm_name]]$R0, w = w, S0 = adv$S, U0 = adv$U,
          N_start_vec = demog$N_start, N_end_vec = demog$N_end, births_vec = demog$births, deaths_vec = demog$deaths,
          X_recent_history = X_recent, lambda_import = DIAG_LAMBDA_IMPORT, stochastic = TRUE,
          rng_pois = function(n, lam) rpois(n, lam), rng_binom = function(n, size, p) rbinom(n, size, p)
        )
        Y_rep <- simulate_reported_cases(sim$X, q = sd_clim$q, phi_obs = draws$phi_obs[d], stochastic = TRUE,
                                          rng_nb = function(n, mu, phi) rnbinom(n, mu = mu, size = phi))
        traj_list[[slot]] <- tibble(
          posterior_draw_id = d, stochastic_replicate_id = r, arm = arm_name,
          week_start = future_dates, S = sim$S, N = sim$N, X_latent = sim$X, Y_reported = Y_rep, R0 = sim$R0, R_eff = sim$R_eff
        )
      }
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  message(sprintf("[Diag 2] %-8s %3d draws x %2d reps x %d arms = %6d trajectories | %.1fs",
                   label, n_draws, n_replicates, length(DECOMP_ARMS), n_draws * n_replicates * length(DECOMP_ARMS), elapsed))
  bind_rows(traj_list)
}

# ================================================================
# Per-trajectory summary metrics (Section 4)
# ================================================================
summarise_trajectories <- function(traj) {
  traj <- traj |> mutate(year = as.integer(format(week_start, "%Y")))
  per_traj <- traj |> group_by(posterior_draw_id, stochastic_replicate_id, arm) |>
    summarise(
      cum_10y = sum(X_latent), incidence_per_100k_10y = 1e5 * sum(X_latent) / mean(N),
      weeks_Reff_gt1 = sum(R_eff > 1), frac_weeks_Reff_gt1 = mean(R_eff > 1),
      min_S_over_N = min(S / N), final_S_over_N = S[which.max(week_start)] / N[which.max(week_start)],
      .groups = "drop"
    )
  # number of distinct R_eff>1 "episodes" (runs), per trajectory
  episodes <- traj |> group_by(posterior_draw_id, stochastic_replicate_id, arm) |>
    summarise(n_Reff_episodes = { r <- rle(R_eff > 1); sum(r$values, na.rm = TRUE) }, .groups = "drop")
  annual_R0 <- traj |> group_by(posterior_draw_id, arm, year) |>
    summarise(mean_R0 = mean(R0), peak_R0 = max(R0), .groups = "drop")
  major <- traj |> group_by(posterior_draw_id, stochastic_replicate_id, arm) |>
    group_modify(~ {
      waves <- detect_simulated_waves(.x$week_start, .x$Y_reported, .x$N)
      tibble(has_major_epidemic = if (nrow(waves)) any(waves$major_epidemic) else FALSE, n_major_epidemics = if (nrow(waves)) sum(waves$major_epidemic) else 0L)
    }) |> ungroup()
  list(per_trajectory = per_traj |> left_join(episodes, by = c("posterior_draw_id", "stochastic_replicate_id", "arm")) |> left_join(major, by = c("posterior_draw_id", "stochastic_replicate_id", "arm")),
       annual_R0 = annual_R0)
}

# ================================================================
# Run: smoke test first, then full scale if feasible
# ================================================================
message("[Diag 2] Smoke test...")
smoke_traj <- run_paired_diagnostic(DIAG_ENSEMBLE$smoke$n_draws, DIAG_ENSEMBLE$smoke$n_replicates, "SMOKE")
smoke_summary <- summarise_trajectories(smoke_traj)
write_csv(smoke_summary$per_trajectory, file.path(DIR_RESULTS_DIAG, "DIAG_smoke_per_trajectory_summary.csv"))
message("[Diag 2] Smoke test per-arm median cum_10y:")
print(as.data.frame(smoke_summary$per_trajectory |> group_by(arm) |> summarise(median_cum10y = median(cum_10y), .groups = "drop")))

message("\n[Diag 2] Full-scale paired run...")
full_traj <- run_paired_diagnostic(DIAG_ENSEMBLE$full$n_draws, DIAG_ENSEMBLE$full$n_replicates, "FULL")
saveRDS(full_traj, file.path(DIR_RESULTS_DIAG, "DIAG_full_trajectories.rds")) # gitignored (03_Output/**/*.rds), regenerable by re-running this script
full_summary <- summarise_trajectories(full_traj)
write_csv(full_summary$per_trajectory, file.path(DIR_RESULTS_DIAG, "DIAG_full_per_trajectory_summary.csv"))
write_csv(full_summary$annual_R0, file.path(DIR_RESULTS_DIAG, "DIAG_full_annual_R0_by_arm.csv"))

arm_summary <- full_summary$per_trajectory |> group_by(arm) |>
  summarise(n = n(), median_cum10y = median(cum_10y), mean_cum10y = mean(cum_10y),
            median_min_S_over_N = median(min_S_over_N), median_final_S_over_N = median(final_S_over_N),
            median_frac_weeks_Reff_gt1 = median(frac_weeks_Reff_gt1), median_n_Reff_episodes = median(n_Reff_episodes),
            proportion_major_epidemic = mean(has_major_epidemic), .groups = "drop")
message("\n=== Full-scale per-arm summary ===")
print(as.data.frame(arm_summary), digits = 4)
write_csv(arm_summary, file.path(DIR_TABLES_DIAG, "DIAG_arm_summary.csv"))

# ================================================================
# Paired contrasts (at the SAME draw x replicate across arms) and the
# 2x2 factorial decomposition: D-A = (B-A) + (C-A) + interaction
# ================================================================
wide <- full_summary$per_trajectory |> select(posterior_draw_id, stochastic_replicate_id, arm, cum_10y) |>
  pivot_wider(names_from = arm, values_from = cum_10y)
paired_diff <- wide |> mutate(
  diff_B_minus_A = B_year_only - A_neutral, diff_C_minus_A = C_amplitude_only - A_neutral,
  diff_D_minus_A = D_both - A_neutral, interaction = D_both - B_year_only - C_amplitude_only + A_neutral
)
write_csv(paired_diff, file.path(DIR_RESULTS_DIAG, "DIAG_paired_differences_cum10y.csv"))

paired_summary <- tibble(
  contrast = c("B - A (year-level effect)", "C - A (amplitude effect)", "D - A (both, total)", "Interaction (D-B-C+A)"),
  median_diff = c(median(paired_diff$diff_B_minus_A), median(paired_diff$diff_C_minus_A), median(paired_diff$diff_D_minus_A), median(paired_diff$interaction)),
  lo95 = c(quantile(paired_diff$diff_B_minus_A, .025), quantile(paired_diff$diff_C_minus_A, .025), quantile(paired_diff$diff_D_minus_A, .025), quantile(paired_diff$interaction, .025)),
  hi95 = c(quantile(paired_diff$diff_B_minus_A, .975), quantile(paired_diff$diff_C_minus_A, .975), quantile(paired_diff$diff_D_minus_A, .975), quantile(paired_diff$interaction, .975)),
  fraction_positive = c(mean(paired_diff$diff_B_minus_A > 0), mean(paired_diff$diff_C_minus_A > 0), mean(paired_diff$diff_D_minus_A > 0), mean(paired_diff$interaction > 0))
)
message("\n=== Paired contrasts on cum_10y (same draw x replicate across arms) ===")
print(as.data.frame(paired_summary), digits = 4)
write_csv(paired_summary, file.path(DIR_TABLES_DIAG, "DIAG_paired_contrast_summary.csv"))

# Also on log1p scale (more appropriate given extreme right-skew/population bound)
wide_log <- wide |> mutate(across(c(A_neutral, B_year_only, C_amplitude_only, D_both), ~ log1p(.x)))
paired_diff_log <- wide_log |> mutate(
  diff_B_minus_A = B_year_only - A_neutral, diff_C_minus_A = C_amplitude_only - A_neutral,
  diff_D_minus_A = D_both - A_neutral, interaction = D_both - B_year_only - C_amplitude_only + A_neutral
)
paired_summary_log <- tibble(
  contrast = c("B - A (year-level effect)", "C - A (amplitude effect)", "D - A (both, total)", "Interaction (D-B-C+A)"),
  median_diff_log1p = c(median(paired_diff_log$diff_B_minus_A), median(paired_diff_log$diff_C_minus_A), median(paired_diff_log$diff_D_minus_A), median(paired_diff_log$interaction))
)
message("\n=== Same contrasts on log1p(cum_10y) scale (variance-decomposition-appropriate for skewed data) ===")
print(as.data.frame(paired_summary_log), digits = 4)
write_csv(paired_summary_log, file.path(DIR_TABLES_DIAG, "DIAG_paired_contrast_summary_log1p.csv"))

message("\n[Diag 2] Paired 4-arm simulation complete. Tables -> ", DIR_TABLES_DIAG, " | Results -> ", DIR_RESULTS_DIAG)
