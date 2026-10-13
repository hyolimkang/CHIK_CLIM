# PHASE 1 -- importation intensity sensitivity: lambda_import in
# {0, 0.25, 1, 4}, constant throughout the holdout. PRIMARY: H1/G2
# (2020-2025). SECONDARY: H2/G2 (2022-2025). All four arms per holdout
# share the SAME seed_offset (hence the SAME sampled posterior draws
# and the SAME z_year/z_season/A_year/climate realisations per draw --
# only the lambda-dependent, shared-uniform-coupled M_t differs), per
# the task's common-random-number requirement. lambda=1 reuses Phase 0's
# already-run baseline (same seed) rather than resimulating it.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/00_config_step7.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/01_simulation_engine.R"))
source(file.path(DIR_STEP6, "00_config_step6.R"))

max_run <- function(x) { rr <- rle(replace(x, is.na(x), FALSE)); if (!any(rr$values)) 0L else max(rr$lengths[rr$values]) }
qv <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))

training_bundles <- list(H1 = readRDS(FIT_PATHS[["H1"]]), H2 = readRDS(FIT_PATHS[["H2"]]))
phase0 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))

SEED_H1 <- BASE_SEED_STEP7 + 1L  # SAME as Phase 0's H1 baseline -- common random numbers across the lambda grid
SEED_H2 <- BASE_SEED_STEP7 + 2L  # SAME as Phase 0's H2 baseline

run_grid <- function(holdout_id, seed_offset, baseline_res) {
  res_list <- list()
  for (nm in names(LAMBDA_GRID)) {
    lam <- LAMBDA_GRID[[nm]]
    if (lam == LAMBDA_BASELINE) { res_list[[nm]] <- baseline_res; message(sprintf("  [%s / %s / lambda=%.2f] reusing Phase0 baseline (same seed)", holdout_id, nm, lam)); next }
    message(sprintf("  [%s / %s / lambda=%.2f] simulating (full ensemble, shared-uniform-coupled to the baseline)...", holdout_id, nm, lam))
    res_list[[nm]] <- simulate_holdout_step7(training_bundles[[holdout_id]], HOLDOUTS[[holdout_id]], lambda_mode = "constant", lambda_const = lam,
                                              year_effect_mode = "centered", n_draws = ENSEMBLE_STEP7$full$n_draws, n_replicates = ENSEMBLE_STEP7$full$n_replicates,
                                              seed_offset = seed_offset)
  }
  res_list
}

message("[Phase1] H1/G2 full importation grid (PRIMARY)...")
grid_H1 <- run_grid("H1", SEED_H1, phase0$H1)
message("[Phase1] H2/G2 full importation grid (SECONDARY)...")
grid_H2 <- run_grid("H2", SEED_H2, phase0$H2)
saveRDS(grid_H1, file.path(DIR_RESULTS_STEP7, "PHASE1_H1_importation_grid.rds"))
saveRDS(grid_H2, file.path(DIR_RESULTS_STEP7, "PHASE1_H2_importation_grid.rds"))

# ================================================================
# Outcome extraction, common to both holdouts.
# ================================================================
summarise_arm <- function(res, holdout_id, lambda_id, lambda_val, training_bundle) {
  traj <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y")),
                                      S_prop = S / N, import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                                      favourable = R_eff > 1, locally_dominant = X_latent > 0 & local_hazard > import_hazard)

  # A. Pre-2022 persistence (H1 only: weeks 2020-2021 exist in the holdout trajectory)
  pre2022 <- traj |> filter(week_start < as.Date("2022-01-02"))
  persistence <- if (nrow(pre2022)) pre2022 |> group_by(posterior_draw_id, stochastic_replicate_id) |>
      summarise(cumulative_latent_2020_2021 = sum(X_latent), weeks_zero_infections = sum(X_latent == 0),
                longest_infection_free_run = max_run(X_latent == 0), frac_weeks_zero = mean(X_latent == 0), .groups = "drop") else NULL

  # B. Entry-2022 infectiousness
  entry <- if (holdout_id == "H1") {
    pre2022 |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
      summarise(entry_weighted_infectiousness = sum(W_GEN * rev(tail(X_latent, G_GEN))), .groups = "drop")
  } else {
    dr <- rstan::extract(training_bundle$fit, pars = "X", permuted = TRUE); n_train <- training_bundle$stan_data$N
    draw_ids <- sort(unique(res$year_effects$posterior_draw_id)); X_tail <- dr$X[draw_ids, (n_train - G_GEN + 1):n_train, drop = FALSE]
    tibble(posterior_draw_id = draw_ids, stochastic_replicate_id = NA_integer_, entry_weighted_infectiousness = apply(X_tail, 1, function(x) sum(W_GEN * rev(x))))
  }

  # C. 2022 resurgence
  annual_2022 <- traj |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(S_prop_start = first(S_prop), peak_R0 = max(R0), peak_Reff = max(R_eff), weeks_Reff_gt1 = sum(favourable),
              longest_run_Reff_gt1 = max_run(favourable), cumulative_latent_infections = sum(X_latent), annual_reported_cases = sum(Y_reported),
              favourable_pressure = any(favourable & M_import > 0), longest_locally_dominant_run = max_run(locally_dominant),
              population_mean = mean(N), .groups = "drop") |>
    mutate(established = favourable_pressure & longest_locally_dominant_run >= ESTABLISHMENT_RUN_WEEKS,
           annual_reported_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
           major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM & annual_reported_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K)

  observed_2022 <- 50986
  cell_summary <- tibble(
    holdout = holdout_id, lambda_id = lambda_id, lambda_value = lambda_val,
    median_entry_weighted_infectiousness = median(entry$entry_weighted_infectiousness),
    lo95_entry = qv(entry$entry_weighted_infectiousness, .025), hi95_entry = qv(entry$entry_weighted_infectiousness, .975),
    prob_entry_near_zero = mean(entry$entry_weighted_infectiousness < 1),
    median_cumulative_latent_2020_2021 = if (!is.null(persistence)) median(persistence$cumulative_latent_2020_2021) else NA_real_,
    median_frac_weeks_zero_2020_2021 = if (!is.null(persistence)) median(persistence$frac_weeks_zero) else NA_real_,
    median_peak_R0 = median(annual_2022$peak_R0), median_peak_Reff = median(annual_2022$peak_Reff),
    median_weeks_Reff_gt1 = median(annual_2022$weeks_Reff_gt1),
    prob_favourable_pressure = mean(annual_2022$favourable_pressure), prob_established = mean(annual_2022$established),
    prob_major_epidemic = mean(annual_2022$major_epidemic),
    prob_exceed_observed_2022 = mean(annual_2022$annual_reported_cases >= observed_2022),
    median_annual_reported_cases_2022 = median(annual_2022$annual_reported_cases),
    n_trajectories = nrow(annual_2022)
  )
  list(cell_summary = cell_summary, entry = entry |> mutate(holdout = holdout_id, lambda_id = lambda_id),
       annual_2022 = annual_2022 |> mutate(holdout = holdout_id, lambda_id = lambda_id),
       persistence = if (!is.null(persistence)) persistence |> mutate(holdout = holdout_id, lambda_id = lambda_id) else NULL)
}

message("[Phase1] Summarising outcomes across the grid...")
all_summaries <- list()
all_entry <- list(); all_annual <- list(); all_persistence <- list()
for (nm in names(LAMBDA_GRID)) {
  s1 <- summarise_arm(grid_H1[[nm]], "H1", nm, LAMBDA_GRID[[nm]], training_bundles$H1)
  s2 <- summarise_arm(grid_H2[[nm]], "H2", nm, LAMBDA_GRID[[nm]], training_bundles$H2)
  all_summaries[[paste0("H1_", nm)]] <- s1$cell_summary; all_summaries[[paste0("H2_", nm)]] <- s2$cell_summary
  all_entry[[paste0("H1_", nm)]] <- s1$entry; all_entry[[paste0("H2_", nm)]] <- s2$entry
  all_annual[[paste0("H1_", nm)]] <- s1$annual_2022; all_annual[[paste0("H2_", nm)]] <- s2$annual_2022
  all_persistence[[paste0("H1_", nm)]] <- s1$persistence
}
cell_summary_all <- bind_rows(all_summaries)
write_csv(cell_summary_all, file.path(DIR_TABLES_STEP7, "PHASE1_importation_grid_summary.csv"))
write_csv(bind_rows(all_entry), file.path(DIR_RESULTS_STEP7, "PHASE1_entry_infectiousness_draws.csv"))
write_csv(bind_rows(all_annual), file.path(DIR_RESULTS_STEP7, "PHASE1_annual_2022_outcomes.csv"))
write_csv(bind_rows(all_persistence), file.path(DIR_RESULTS_STEP7, "PHASE1_pre2022_persistence_H1.csv"))

message("\n=== PHASE 1: importation grid summary (H1 primary, H2 secondary) ===")
print(as.data.frame(cell_summary_all |> select(holdout, lambda_id, lambda_value, median_entry_weighted_infectiousness,
                                                 median_cumulative_latent_2020_2021, prob_established, prob_major_epidemic,
                                                 prob_exceed_observed_2022, median_annual_reported_cases_2022)), digits = 4)

message("\n[Phase1] Importation grid complete -> ", DIR_TABLES_STEP7)
