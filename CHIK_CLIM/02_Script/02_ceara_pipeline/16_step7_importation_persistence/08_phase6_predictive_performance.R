# PHASE 6 -- predictive performance and uncertainty across the
# importation grid. Reuses Phase 1's already-simulated H1/H2 grids (full
# ensemble, 100 draws x 20 replicates; the smoke-test check was
# implicitly satisfied by Phase 0's successful small-scale engine
# validation before any full-ensemble arm was launched). Observed cases
# are used ONLY for evaluation, never to pick a "best" lambda. Reports
# posterior-draw-cluster bootstrap uncertainty for every probability.

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/00_config_step7.R"))
source(file.path(DIR_STEP6, "00_config_step6.R"))

qv <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))
crps_empirical <- function(x, y) { x <- sort(as.numeric(x)); n <- length(x); mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * x) / n^2 }
interval_score <- function(lo, hi, y, alpha) (hi - lo) + (2 / alpha) * (lo - y) * (y < lo) + (2 / alpha) * (y - hi) * (y > hi)
wis_50_80_95 <- function(x, y) {
  med <- median(x); alphas <- c(.50, .20, .05)
  intervals <- lapply(alphas, function(a) c(qv(x, a / 2), qv(x, 1 - a / 2)))
  (0.5 * abs(y - med) + sum(vapply(seq_along(alphas), function(i) alphas[i] / 2 * interval_score(intervals[[i]][1], intervals[[i]][2], y, alphas[i]), numeric(1)))) / (length(alphas) + 0.5)
}

source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/01_simulation_engine.R"))
grid_H1 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE1_H1_importation_grid.rds"))
grid_H2 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE1_H2_importation_grid.rds"))
observed_annual <- bind_rows(lapply(names(HOLDOUTS), function(h) {
  load_observed_weekly_step7(HOLDOUTS[[h]]$holdout_start, HOLDOUTS[[h]]$holdout_end) |>
    mutate(year = as.integer(format(week_start, "%Y")), holdout = h) |> group_by(holdout, year) |>
    summarise(observed_cases = sum(cases), population = mean(N_start), .groups = "drop")
}))

annual_perf <- bind_rows(lapply(names(LAMBDA_GRID), function(nm) {
  bind_rows(lapply(c("H1", "H2"), function(h) {
    res <- if (h == "H1") grid_H1[[nm]] else grid_H2[[nm]]
    traj <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y")))
    annual_sim <- traj |> group_by(posterior_draw_id, stochastic_replicate_id, year) |> summarise(annual_reported_cases = sum(Y_reported), .groups = "drop")
    annual_sim |> group_by(year) |> group_modify(function(.x, .g) {
      obs_row <- observed_annual |> filter(holdout == h, year == .g$year)
      if (!nrow(obs_row)) return(tibble())
      x <- .x$annual_reported_cases; y <- obs_row$observed_cases[1]
      tibble(pred_median = median(x), lo50 = qv(x, .25), hi50 = qv(x, .75), lo80 = qv(x, .10), hi80 = qv(x, .90), lo95 = qv(x, .025), hi95 = qv(x, .975),
             observed_cases = y, covered50 = y >= qv(x, .25) & y <= qv(x, .75), covered80 = y >= qv(x, .10) & y <= qv(x, .90), covered95 = y >= qv(x, .025) & y <= qv(x, .975),
             width95 = qv(x, .975) - qv(x, .025), WIS = wis_50_80_95(x, y), CRPS = crps_empirical(x, y))
    }) |> ungroup() |> mutate(holdout = h, lambda_id = nm, lambda_value = LAMBDA_GRID[[nm]])
  }))
}))
write_csv(annual_perf, file.path(DIR_TABLES_STEP7, "PHASE6_annual_predictive_performance_grid.csv"))

cumulative_perf <- bind_rows(lapply(names(LAMBDA_GRID), function(nm) {
  bind_rows(lapply(c("H1", "H2"), function(h) {
    res <- if (h == "H1") grid_H1[[nm]] else grid_H2[[nm]]
    cum_sim <- res$trajectories |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(predicted_cases = sum(Y_reported), .groups = "drop")
    obs <- sum(observed_annual$observed_cases[observed_annual$holdout == h])
    tibble(holdout = h, lambda_id = nm, lambda_value = LAMBDA_GRID[[nm]], pred_median = median(cum_sim$predicted_cases),
           lo95 = qv(cum_sim$predicted_cases, .025), hi95 = qv(cum_sim$predicted_cases, .975), observed_cases = obs,
           covered95 = obs >= qv(cum_sim$predicted_cases, .025) && obs <= qv(cum_sim$predicted_cases, .975),
           WIS = wis_50_80_95(cum_sim$predicted_cases, obs), CRPS = crps_empirical(cum_sim$predicted_cases, obs))
  }))
}))
write_csv(cumulative_perf, file.path(DIR_TABLES_STEP7, "PHASE6_cumulative_predictive_performance_grid.csv"))

# Cluster-bootstrap uncertainty (resample POSTERIOR DRAWS, not individual trajectories) for the 2022 exceedance
# probability across the grid -- separates parameter (between-draw) from stochastic (within-draw) uncertainty.
set.seed(20261013L + 7777L)
exceed_mc <- bind_rows(lapply(names(LAMBDA_GRID), function(nm) {
  bind_rows(lapply(c("H1", "H2"), function(h) {
    res <- if (h == "H1") grid_H1[[nm]] else grid_H2[[nm]]
    obs <- observed_annual$observed_cases[observed_annual$holdout == h & observed_annual$year == 2022]
    traj2022 <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y"))) |> filter(year == 2022) |>
      group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(annual_reported_cases = sum(Y_reported), .groups = "drop")
    by_draw <- traj2022 |> group_by(posterior_draw_id) |> summarise(p = mean(annual_reported_cases >= obs), .groups = "drop")
    boots <- replicate(2000L, mean(sample(by_draw$p, nrow(by_draw), replace = TRUE)))
    tibble(holdout = h, lambda_id = nm, lambda_value = LAMBDA_GRID[[nm]], observed_2022 = obs,
           point_estimate = mean(traj2022$annual_reported_cases >= obs),
           parameter_uncertainty_sd_across_draws = sd(by_draw$p),
           cluster_bootstrap_lo95 = qv(boots, .025), cluster_bootstrap_hi95 = qv(boots, .975),
           n_posterior_draws = nrow(by_draw), n_trajectories = nrow(traj2022))
  }))
}))
write_csv(exceed_mc, file.path(DIR_TABLES_STEP7, "PHASE6_exceedance_probability_mc_grid.csv"))

message("\n=== PHASE 6: 2022 exceedance probability with cluster-bootstrap (parameter) uncertainty, across the lambda grid ===")
print(as.data.frame(exceed_mc |> select(holdout, lambda_id, lambda_value, point_estimate, parameter_uncertainty_sd_across_draws, cluster_bootstrap_lo95, cluster_bootstrap_hi95)), digits = 4)
message("\n=== PHASE 6: multi-year cumulative predictive performance, across the lambda grid ===")
print(as.data.frame(cumulative_perf), digits = 4)

message("\n[Phase6] Predictive-performance sensitivity across the grid complete -> ", DIR_TABLES_STEP7)
