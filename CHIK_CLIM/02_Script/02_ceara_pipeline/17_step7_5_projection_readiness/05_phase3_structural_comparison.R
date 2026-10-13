# PHASE 3 -- does episodic introduction (same mean pressure) produce
# more realistic persistence/amplification, or merely occasional
# artificial outbreaks? Pure analysis of Phase 2's already-simulated
# Model A/B trajectories -- no new simulation. Uses a REVISED
# establishment proxy that does NOT require a same-year import (Phase 0
# item 5): `chain_active` = sustained local transmission regardless of
# whether a fresh import coincided with that specific year.

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/00_config_step7_5.R"))

max_run <- function(x) { rr <- rle(replace(x, is.na(x), FALSE)); if (!any(rr$values)) 0L else max(rr$lengths[rr$values]) }
qv <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))
crps_empirical <- function(x, y) { x <- sort(as.numeric(x)); n <- length(x); mean(abs(x - y)) - sum((2 * seq_len(n) - n - 1) * x) / n^2 }
interval_score <- function(lo, hi, y, alpha) (hi - lo) + (2 / alpha) * (lo - y) * (y < lo) + (2 / alpha) * (y - hi) * (y > hi)
wis_50_80_95 <- function(x, y) {
  med <- median(x); alphas <- c(.50, .20, .05)
  intervals <- lapply(alphas, function(a) c(qv(x, a / 2), qv(x, 1 - a / 2)))
  (0.5 * abs(y - med) + sum(vapply(seq_along(alphas), function(i) alphas[i] / 2 * interval_score(intervals[[i]][1], intervals[[i]][2], y, alphas[i]), numeric(1)))) / (length(alphas) + 0.5)
}

arms <- readRDS(file.path(DIR_RESULTS_STEP75, "PHASE2_structural_arms.rds"))
observed_weekly <- load_observed_weekly_step7(HOLDOUTS$H1$holdout_start, HOLDOUTS$H1$holdout_end) |> mutate(year = as.integer(format(week_start, "%Y")))
observed_annual <- observed_weekly |> group_by(year) |> summarise(observed_cases = sum(cases), .groups = "drop")
observed_2022 <- observed_annual$observed_cases[observed_annual$year == 2022]
observed_full_period <- sum(observed_annual$observed_cases)

ARM_LABS <- c(ModelA_constant = "Model A: constant Poisson(lambda=1)", ModelB_p10 = "Model B: episodic p=0.10 (mean=1)", ModelB_p25 = "Model B: episodic p=0.25 (mean=1)")

decompose <- function(res, arm_id) {
  traj <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y")),
                                      import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                                      favourable = R_eff > 1,
                                      # introduction EVENT: a week with M_import>0 (regardless of Reff)
                                      introduction_event = M_import > 0,
                                      # successful LOCAL CHAIN (revised, year-agnostic establishment proxy): local
                                      # transmission dominates the hazard for several consecutive weeks, with NO
                                      # requirement that a fresh import coincide with THIS specific calendar year.
                                      locally_dominant = X_latent > 0 & local_hazard > import_hazard)

  entry <- traj |> filter(week_start < as.Date("2022-01-02")) |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |>
    group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(entry_weighted_infectiousness = sum(W_GEN * rev(tail(X_latent, G_GEN))), .groups = "drop")

  # "Chain active" across the WHOLE pre-2022 period (2020-2021), year-agnostic: a run of >=4 consecutive
  # locally-dominant weeks ANYWHERE in 2020-2021 -- does not require the triggering import to fall in 2022.
  chain_status <- traj |> filter(year %in% 2020:2021) |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |>
    group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(any_introduction_event = any(introduction_event), longest_chain_run = max_run(locally_dominant),
              chain_active = max_run(locally_dominant) >= ESTABLISHMENT_RUN_WEEKS, .groups = "drop")

  annual_2022 <- traj |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(weeks_Reff_gt1 = sum(favourable), longest_run_Reff_gt1 = max_run(favourable),
              annual_reported_cases = sum(Y_reported), annual_latent = sum(X_latent),
              first_sustained_dominance_week = { idx <- which(locally_dominant); if (!length(idx)) NA_integer_ else idx[1] },
              population_mean = mean(N), .groups = "drop") |>
    mutate(annual_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
           major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM & annual_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K,
           reaches_observed_2022 = annual_reported_cases >= observed_2022)

  full_period <- traj |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(cumulative_reported = sum(Y_reported), .groups = "drop")

  annual_all <- traj |> group_by(posterior_draw_id, stochastic_replicate_id, year) |> summarise(annual_reported_cases = sum(Y_reported), .groups = "drop")
  wis_by_year <- annual_all |> group_by(year) |> group_modify(function(.x, .g) {
    obs <- observed_annual$observed_cases[observed_annual$year == .g$year]; if (!length(obs)) return(tibble())
    tibble(WIS = wis_50_80_95(.x$annual_reported_cases, obs), CRPS = crps_empirical(.x$annual_reported_cases, obs), observed = obs, pred_median = median(.x$annual_reported_cases))
  }) |> ungroup() |> mutate(arm = arm_id)

  tibble(
    arm = arm_id, arm_label = ARM_LABS[[arm_id]],
    median_entry_weighted_infectiousness = median(entry$entry_weighted_infectiousness),
    lo95_entry = qv(entry$entry_weighted_infectiousness, .025), hi95_entry = qv(entry$entry_weighted_infectiousness, .975),
    median_local_to_import_hazard_ratio_favourable = median((traj |> filter(favourable, year == 2022))$local_hazard / pmax((traj |> filter(favourable, year == 2022))$import_hazard, 1e-12)),
    prob_any_introduction_event_2020_2021 = mean(chain_status$any_introduction_event),
    prob_chain_active_2020_2021 = mean(chain_status$chain_active),
    median_longest_chain_run_2020_2021 = median(chain_status$longest_chain_run),
    median_weeks_Reff_gt1_2022 = median(annual_2022$weeks_Reff_gt1),
    median_longest_run_Reff_gt1_2022 = median(annual_2022$longest_run_Reff_gt1),
    prob_major_epidemic_2022 = mean(annual_2022$major_epidemic),
    prob_reaches_observed_2022 = mean(annual_2022$reaches_observed_2022),
    median_annual_reported_2022 = median(annual_2022$annual_reported_cases),
    median_cumulative_reported_full_period = median(full_period$cumulative_reported),
    lo95_cumulative = qv(full_period$cumulative_reported, .025), hi95_cumulative = qv(full_period$cumulative_reported, .975),
    width95_cumulative = qv(full_period$cumulative_reported, .975) - qv(full_period$cumulative_reported, .025),
    WIS_cumulative_full_period = wis_50_80_95(full_period$cumulative_reported, observed_full_period),
    CRPS_cumulative_full_period = crps_empirical(full_period$cumulative_reported, observed_full_period),
    n_trajectories = nrow(annual_2022)
  ) -> arm_summary
  list(summary = arm_summary, wis_by_year = wis_by_year)
}

message("[Phase3] Decomposing Model A vs Model B (p=0.10, p=0.25) outcomes...")
results <- lapply(names(arms), function(nm) decompose(arms[[nm]], nm))
names(results) <- names(arms)
structural_summary <- bind_rows(lapply(results, `[[`, "summary"))
wis_by_year_all <- bind_rows(lapply(results, `[[`, "wis_by_year"))

write_csv(structural_summary, file.path(DIR_TABLES_STEP75, "PHASE3_structural_comparison_summary.csv"))
write_csv(wis_by_year_all, file.path(DIR_TABLES_STEP75, "PHASE3_WIS_by_year_structural_arms.csv"))

message("\n=== PHASE 3: Model A vs Model B structural comparison (H1, full 2020-2025 holdout) ===")
print(as.data.frame(structural_summary |> select(arm, median_entry_weighted_infectiousness, prob_chain_active_2020_2021,
                                                    prob_major_epidemic_2022, prob_reaches_observed_2022, median_annual_reported_2022,
                                                    median_cumulative_reported_full_period, width95_cumulative, WIS_cumulative_full_period)), digits = 4)

# "Occasional artificial outbreak" check: distribution of PEAK annual reported cases
# in years OTHER than 2022, across arms -- are episodic bursts creating unrealistic spikes elsewhere?
other_year_peaks <- bind_rows(lapply(names(arms), function(nm) {
  arms[[nm]]$trajectories |> mutate(year = as.integer(format(week_start, "%Y"))) |> filter(year != 2022) |>
    group_by(posterior_draw_id, stochastic_replicate_id, year) |> summarise(annual_reported_cases = sum(Y_reported), .groups = "drop") |>
    group_by(year) |> summarise(median_annual = median(annual_reported_cases), p95_annual = qv(annual_reported_cases, .95), p99_annual = qv(annual_reported_cases, .99), .groups = "drop") |>
    mutate(arm = nm)
}))
write_csv(other_year_peaks, file.path(DIR_TABLES_STEP75, "PHASE3_other_year_peaks_artificial_outbreak_check.csv"))
message("\n=== Phase3: non-2022 year peaks, by arm (checking for artificial outbreaks outside the target year) ===")
print(as.data.frame(other_year_peaks |> pivot_wider(names_from = arm, values_from = c(median_annual, p95_annual, p99_annual))), digits = 4)

message("\n[Phase3] Structural comparison complete -> ", DIR_TABLES_STEP75)
