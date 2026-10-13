# PHASE 1 -- why does H1's free-running generator produce so little
# infection pressure entering 2022? Full 2020-2021 analysis, reusing
# Step 7's already-simulated, already-validated lambda=1 baseline
# (NO new simulation). Observed surveillance used ONLY for evaluation.

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/00_config_step7_5.R"))

max_run <- function(x) { rr <- rle(replace(x, is.na(x), FALSE)); if (!any(rr$values)) 0L else max(rr$lengths[rr$values]) }
qv <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))

modelA <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))$H1
traj <- modelA$trajectories |> mutate(year = as.integer(format(week_start, "%Y")),
                                       import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                                       favourable = R_eff > 1, locally_dominant = X_latent > 0 & local_hazard > import_hazard)

observed_weekly <- load_observed_weekly_step7(HOLDOUTS$H1$holdout_start, HOLDOUTS$H1$holdout_end) |> mutate(year = as.integer(format(week_start, "%Y")))

# ================================================================
# Weekly predictive distribution (2020-2021 only) vs. observed.
# ================================================================
weekly_q <- traj |> filter(year %in% 2020:2021) |> group_by(week_start) |>
  summarise(across(c(X_latent, Y_reported, R0, R_eff, M_import, local_hazard, import_hazard),
                   list(median = median, lo50 = ~qv(.x, .25), hi50 = ~qv(.x, .75), lo95 = ~qv(.x, .025), hi95 = ~qv(.x, .975))), .groups = "drop")
weekly_q <- weekly_q |> left_join(observed_weekly |> select(week_start, observed_cases = cases, N_start), by = "week_start")
write_csv(weekly_q, file.path(DIR_RESULTS_STEP75, "PHASE1_weekly_2020_2021_predictive_vs_observed.csv"))

# ================================================================
# Annual predictive distribution + coverage, 2020-2021.
# ================================================================
annual_sim <- traj |> filter(year %in% 2020:2021) |> group_by(posterior_draw_id, stochastic_replicate_id, year) |> summarise(annual_reported_cases = sum(Y_reported), .groups = "drop")
observed_annual <- observed_weekly |> filter(year %in% 2020:2021) |> group_by(year) |> summarise(observed_cases = sum(cases), .groups = "drop")
annual_perf <- annual_sim |> group_by(year) |> group_modify(function(.x, .g) {
  obs <- observed_annual$observed_cases[observed_annual$year == .g$year]; x <- .x$annual_reported_cases
  tibble(pred_median = median(x), lo50 = qv(x, .25), hi50 = qv(x, .75), lo95 = qv(x, .025), hi95 = qv(x, .975),
         observed_cases = obs, covered50 = obs >= qv(x, .25) & obs <= qv(x, .75), covered95 = obs >= qv(x, .025) & obs <= qv(x, .975))
}) |> ungroup()
write_csv(annual_perf, file.path(DIR_TABLES_STEP75, "PHASE1_annual_2020_2021_coverage.csv"))
message("=== Phase1: annual 2020-2021 predictive coverage (H1/G2, lambda=1) ==="); print(as.data.frame(annual_perf), digits = 4)

# ================================================================
# Sustained local transmission, infection-free periods, and the A/B/C/D
# decomposition for WHEN the shortfall develops.
# ================================================================
per_traj_2020 <- traj |> filter(year == 2020) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
  summarise(weeks_zero_X = sum(X_latent == 0), longest_zero_run = max_run(X_latent == 0), weeks_Reff_gt1 = sum(favourable),
            longest_locally_dominant_run = max_run(locally_dominant), cumulative_X = sum(X_latent), .groups = "drop") |> mutate(year = 2020L)
per_traj_2021 <- traj |> filter(year == 2021) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
  summarise(weeks_zero_X = sum(X_latent == 0), longest_zero_run = max_run(X_latent == 0), weeks_Reff_gt1 = sum(favourable),
            longest_locally_dominant_run = max_run(locally_dominant), cumulative_X = sum(X_latent), .groups = "drop") |> mutate(year = 2021L)
year_summary <- bind_rows(per_traj_2020, per_traj_2021) |> group_by(year) |>
  summarise(median_weeks_zero_X = median(weeks_zero_X), median_longest_zero_run = median(longest_zero_run),
            median_weeks_Reff_gt1 = median(weeks_Reff_gt1), median_longest_dominant_run = median(longest_locally_dominant_run),
            median_cumulative_X = median(cumulative_X), prob_near_total_fadeout = mean(weeks_zero_X >= 45), .groups = "drop")
write_csv(year_summary, file.path(DIR_TABLES_STEP75, "PHASE1_annual_persistence_summary.csv"))
message("\n=== Phase1: year-by-year persistence summary (where does the shortfall develop?) ==="); print(as.data.frame(year_summary), digits = 4)

# Entry to January 2022 (reused from Step7's own definition for consistency)
entry <- traj |> filter(week_start < as.Date("2022-01-02")) |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |>
  group_by(posterior_draw_id, stochastic_replicate_id) |>
  summarise(entry_weighted_infectiousness = sum(W_GEN * rev(tail(X_latent, G_GEN))), .groups = "drop")
H2_ref <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))$H2
# H2's entry reference, computed exactly as in Step7 Phase1/Phase0 (training-posterior X, labelled diagnostic reference ONLY)
training_H2 <- readRDS(FIT_PATHS[["H2"]])
dr_H2 <- rstan::extract(training_H2$fit, pars = "X", permuted = TRUE); n_train_H2 <- training_H2$stan_data$N
draw_ids_H2 <- sort(unique(H2_ref$year_effects$posterior_draw_id)); X_tail_H2 <- dr_H2$X[draw_ids_H2, (n_train_H2 - G_GEN + 1):n_train_H2, drop = FALSE]
H2_entry_median <- median(apply(X_tail_H2, 1, function(x) sum(W_GEN * rev(x))))
entry_summary <- tibble(source = c("H1 (simulated, this module)", "H2 (training-posterior, DIAGNOSTIC REFERENCE ONLY -- never used to initialise H1)"),
                         median_entry_weighted_infectiousness = c(median(entry$entry_weighted_infectiousness), H2_entry_median),
                         lo95 = c(qv(entry$entry_weighted_infectiousness, .025), NA), hi95 = c(qv(entry$entry_weighted_infectiousness, .975), NA))
write_csv(entry_summary, file.path(DIR_TABLES_STEP75, "PHASE1_entry2022_H1_vs_H2reference.csv"))
message("\n=== Phase1: entry-2022 infectiousness, H1 (simulated) vs H2 (diagnostic reference) ==="); print(as.data.frame(entry_summary), digits = 4)

# A/B/C/D classification evidence (descriptive synthesis of the above, numeric support only)
abcd <- tibble(
  explanation = c("A. Insufficient local transmission growth", "B. Repeated stochastic fadeout", "C. Insufficient introduction pressure", "D. State-propagation/observation-model errors"),
  supported = c(
    sprintf("NOT primarily -- Reff>1 occurs in most trajectories most years (median %.0f weeks/yr in 2020, %.0f in 2021); local/import hazard ratio is large once underway (Step7 Phase4).", year_summary$median_weeks_Reff_gt1[1], year_summary$median_weeks_Reff_gt1[2]),
    sprintf("NOT PRIMARILY, but not irrelevant to the tail -- at the MEDIAN, zero weeks with zero latent infections occur (%.0f of ~52) in either 2020 or 2021, and the median longest locally-dominant run spans nearly the full year (%.0f weeks) -- the typical trajectory sustains persistently LOW, not fully extinct, transmission; P(near-total fadeout, >=45 zero weeks)=%.1f%% (2020), %.1f%% (2021). Fadeout risk is a real but secondary/tail phenomenon, not the median story.", year_summary$median_weeks_zero_X[1], year_summary$median_longest_dominant_run[1], 100*year_summary$prob_near_total_fadeout[1], 100*year_summary$prob_near_total_fadeout[2]),
    sprintf("PRIMARY -- lambda=1 is a small constant rate relative to the gap (median entry infectiousness %.2f vs. H2's reference %.1f, a ~%.0fx gap); see Step7 Phase1-2 for the controlled lambda/pre-post contrasts quantifying this directly.", median(entry$entry_weighted_infectiousness), H2_entry_median, H2_entry_median/median(entry$entry_weighted_infectiousness)),
    "NOT SUPPORTED -- R0 reconstruction matches to 4.4e-16, S+U=N exactly, q reproduced exactly (Step6/Step7 Phase0 numerical-consistency checks); no evidence of a propagation or observation-model error."
  ))
write_csv(abcd, file.path(DIR_TABLES_STEP75, "PHASE1_ABCD_classification.csv"))
message("\n=== Phase1: A/B/C/D decomposition ===")
for (i in seq_len(nrow(abcd))) message(sprintf("%s\n  -> %s\n", abcd$explanation[i], abcd$supported[i]))

message("\n[Phase1] Pre-resurgence persistence analysis complete -> ", DIR_TABLES_STEP75)
