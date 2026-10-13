# PHASE 3 -- introduction TIMING at fixed total pressure. H1/G2 only,
# lambda_import=1 baseline. T0 = the ALREADY-REALISED M_t sequence from
# Phase 0/1's lambda=1 baseline (captured, not redrawn). T1 = the SAME
# per-(draw,replicate) M_t VALUES, randomly permuted WITHIN each
# calendar-year block (so each year's TOTAL introduction pressure is
# unchanged, only which week gets which value changes), replayed
# through the identical renewal mechanics with the SAME rng_binom
# stream position as T0 (the fixed-M rng_pois closure still consumes
# one runif() per week, so the infection-generation stream stays
# aligned -- see 01_simulation_engine.R's documented coupling). This is
# an intervention-style diagnostic: it does NOT claim real importations
# are uniformly distributed or that permuted timing is epidemiologically
# correct -- it isolates whether WHEN (not how much) introduction
# pressure arrives matters, holding the total fixed.

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

training_H1 <- readRDS(FIT_PATHS[["H1"]])
phase0 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE0_baseline_lambda1_simulations.rds"))
SEED_H1 <- BASE_SEED_STEP7 + 1L

message("[Phase3] Building T0 (original realised timing) lookup from the Phase0 lambda=1 baseline...")
T0_traj <- phase0$H1$trajectories |> mutate(year = as.integer(format(week_start, "%Y"))) |> arrange(posterior_draw_id, stochastic_replicate_id, week_start)
week_order <- T0_traj |> distinct(week_start, year) |> arrange(week_start)

# Build, per (draw,replicate), the PERMUTED M_t vector: shuffle within each calendar year.
set.seed(SEED_H1 + 99999L) # ONE fixed permutation per (draw,replicate), deterministic and documented, not re-drawn per analysis
permute_within_year <- function(m_vec, years) {
  out <- m_vec
  for (y in unique(years)) { idx <- which(years == y); out[idx] <- sample(m_vec[idx]) }
  out
}
fixed_M_table <- T0_traj |> group_by(posterior_draw_id, stochastic_replicate_id) |>
  summarise(M_original = list(M_import), .groups = "drop") |>
  rowwise() |> mutate(M_permuted = list(permute_within_year(M_original, week_order$year))) |> ungroup()
fixed_M_lookup_fn <- function(d, r) {
  row <- fixed_M_table[fixed_M_table$posterior_draw_id == d & fixed_M_table$stochastic_replicate_id == r, ]
  if (!nrow(row)) stop("No T0 M_t found for draw ", d, " replicate ", r)
  row$M_permuted[[1]]
}

message("[Phase3] Simulating T1 (permuted timing, same total/year, same posterior draws/replicates)...")
T1 <- simulate_holdout_step7(training_H1, HOLDOUTS$H1, lambda_mode = "fixed_M", fixed_M_lookup = fixed_M_lookup_fn,
                              year_effect_mode = "centered", n_draws = ENSEMBLE_STEP7$full$n_draws, n_replicates = ENSEMBLE_STEP7$full$n_replicates,
                              seed_offset = SEED_H1)
saveRDS(T1, file.path(DIR_RESULTS_STEP7, "PHASE3_T1_permuted_timing_trajectories.rds"))

# Sanity check: per-(draw,replicate) yearly M_t totals must be IDENTICAL between T0 and T1 (permutation preserves the total).
yearly_totals_T0 <- T0_traj |> group_by(posterior_draw_id, stochastic_replicate_id, year) |> summarise(total_M = sum(M_import), .groups = "drop")
yearly_totals_T1 <- T1$trajectories |> mutate(year = as.integer(format(week_start, "%Y"))) |> group_by(posterior_draw_id, stochastic_replicate_id, year) |>
  summarise(total_M = sum(M_import), .groups = "drop")
total_check <- inner_join(yearly_totals_T0, yearly_totals_T1, by = c("posterior_draw_id", "stochastic_replicate_id", "year"), suffix = c("_T0", "_T1")) |>
  mutate(diff = total_M_T0 - total_M_T1)
write_csv(total_check |> summarise(max_abs_diff = max(abs(diff)), n = n()), file.path(DIR_TABLES_STEP7, "PHASE3_timing_permutation_total_preservation_check.csv"))
message(sprintf("[Phase3] Yearly M_t total preservation check: max abs diff across %d (draw,replicate,year) cells = %g (should be 0)", nrow(total_check), max(abs(total_check$diff))))

# ================================================================
# Outcome comparison, T0 vs T1
# ================================================================
summarise_timing_arm <- function(traj, label) {
  traj <- traj |> mutate(year = as.integer(format(week_start, "%Y")), import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                          favourable = R_eff > 1, locally_dominant = X_latent > 0 & local_hazard > import_hazard)
  pre2022 <- traj |> filter(week_start < as.Date("2022-01-02"))
  entry <- pre2022 |> arrange(posterior_draw_id, stochastic_replicate_id, week_start) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(entry_weighted_infectiousness = sum(W_GEN * rev(tail(X_latent, G_GEN))), .groups = "drop")
  pressure_timing <- traj |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(frac_import_during_favourable = sum(M_import[favourable]) / pmax(sum(M_import), 1e-9), .groups = "drop")
  annual_2022 <- traj |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(annual_reported_cases = sum(Y_reported), cumulative_latent_infections = sum(X_latent),
              favourable_pressure = any(favourable & M_import > 0), longest_locally_dominant_run = max_run(locally_dominant),
              population_mean = mean(N), .groups = "drop") |>
    mutate(established = favourable_pressure & longest_locally_dominant_run >= ESTABLISHMENT_RUN_WEEKS,
           annual_reported_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
           major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM & annual_reported_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K)
  tibble(arm = label, median_entry_weighted_infectiousness = median(entry$entry_weighted_infectiousness),
         median_frac_import_during_favourable = median(pressure_timing$frac_import_during_favourable),
         median_cumulative_latent_2022 = median(annual_2022$cumulative_latent_infections),
         median_annual_reported_cases_2022 = median(annual_2022$annual_reported_cases),
         prob_established = mean(annual_2022$established), prob_major_epidemic = mean(annual_2022$major_epidemic),
         prob_exceed_observed_2022 = mean(annual_2022$annual_reported_cases >= 50986), n = nrow(annual_2022))
}
timing_summary <- bind_rows(summarise_timing_arm(phase0$H1$trajectories, "T0_original_timing"), summarise_timing_arm(T1$trajectories, "T1_permuted_timing"))
write_csv(timing_summary, file.path(DIR_TABLES_STEP7, "PHASE3_timing_permutation_summary.csv"))
message("\n=== PHASE 3: introduction timing at fixed total pressure (H1/G2) ==="); print(as.data.frame(timing_summary), digits = 4)

# Paired per-(draw,replicate) deltas (same draws/replicates in both arms -- a valid paired contrast)
paired <- inner_join(
  phase0$H1$trajectories |> mutate(year = as.integer(format(week_start, "%Y"))) |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cases_T0 = sum(Y_reported), .groups = "drop"),
  T1$trajectories |> mutate(year = as.integer(format(week_start, "%Y"))) |> filter(year == 2022) |> group_by(posterior_draw_id, stochastic_replicate_id) |> summarise(cases_T1 = sum(Y_reported), .groups = "drop"),
  by = c("posterior_draw_id", "stochastic_replicate_id")
) |> mutate(delta_T1_minus_T0 = cases_T1 - cases_T0)
write_csv(tibble(median_delta = median(paired$delta_T1_minus_T0), lo95_delta = qv(paired$delta_T1_minus_T0, .025), hi95_delta = qv(paired$delta_T1_minus_T0, .975), n = nrow(paired)),
          file.path(DIR_TABLES_STEP7, "PHASE3_paired_timing_delta_2022_cases.csv"))
message(sprintf("Paired T1-T0 delta in 2022 reported cases: median=%.1f [%.1f, %.1f] (n=%d paired draws x replicates)",
                 median(paired$delta_T1_minus_T0), qv(paired$delta_T1_minus_T0, .025), qv(paired$delta_T1_minus_T0, .975), nrow(paired)))

message("\n[Phase3] Timing-permutation experiment complete -> ", DIR_TABLES_STEP7)
