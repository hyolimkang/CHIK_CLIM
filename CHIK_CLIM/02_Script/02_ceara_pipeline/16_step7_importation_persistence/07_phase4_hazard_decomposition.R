# PHASE 4 -- epidemic amplification: decompose the weekly force of
# infection into local_hazard = R0_t*weighted_infectiousness_t/N_t and
# introduction_hazard = M_t/N_t (an explicit competing-hazard split of
# the already-saved lambda_t = local_hazard + import_hazard; NO
# attribution of individual infection COUNTS to one hazard or the other
# is attempted -- X(t) is generated from the COMBINED hazard, and any
# count-level split would require an arbitrary assumption not supported
# by the renewal model's own structure). Pure analysis of Phase 1's
# already-simulated H1/G2 grid (I0/I2/I3); no new simulation. Outcome-
# conditioned comparisons (large vs. no epidemic) are DESCRIPTIVE, not causal.

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/16_step7_importation_persistence/00_config_step7.R"))
source(file.path(DIR_STEP6, "00_config_step6.R"))

max_run <- function(x) { rr <- rle(replace(x, is.na(x), FALSE)); if (!any(rr$values)) 0L else max(rr$lengths[rr$values]) }
qv <- function(x, p) as.numeric(stats::quantile(x, p, na.rm = TRUE, names = FALSE, type = 8))

grid_H1 <- readRDS(file.path(DIR_RESULTS_STEP7, "PHASE1_H1_importation_grid.rds"))

decompose_arm <- function(res, lambda_id) {
  traj <- res$trajectories |> mutate(year = as.integer(format(week_start, "%Y")),
                                      import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0),
                                      ratio_local_to_import = local_hazard / pmax(import_hazard, 1e-12),
                                      favourable = R_eff > 1)
  traj_2022 <- traj |> filter(year == 2022) |> arrange(posterior_draw_id, stochastic_replicate_id, week_start)

  # Growth dynamics per trajectory: weeks from the FIRST favourable+pressure week to the ANNUAL peak latent
  # incidence; longest run of week-on-week INCREASING X_latent (a simple amplification-duration proxy);
  # time (in weeks from Jan-1) to cross latent-infection thresholds of 100 / 1000; extinction before the
  # seasonal favourable window closes (last favourable week with zero subsequent infections that year).
  growth <- traj_2022 |> group_by(posterior_draw_id, stochastic_replicate_id) |>
    summarise(
      first_favourable_week = { w <- which(favourable & M_import > 0); if (!length(w)) NA_integer_ else w[1] },
      peak_week = which.max(X_latent), peak_X = max(X_latent),
      weeks_favourable_to_peak = { w <- which(favourable & M_import > 0); if (!length(w) || peak_week < w[1]) NA_integer_ else peak_week - w[1] },
      longest_increasing_run = max_run(c(FALSE, diff(X_latent) > 0)),
      week_to_100 = { cs <- cumsum(X_latent); w <- which(cs >= 100); if (!length(w)) NA_integer_ else w[1] },
      week_to_1000 = { cs <- cumsum(X_latent); w <- which(cs >= 1000); if (!length(w)) NA_integer_ else w[1] },
      extinct_after_favourable = { w <- which(favourable); if (!length(w)) NA else all(X_latent[(max(w) + 1):n()] == 0) },
      cumulative_latent = sum(X_latent), annual_reported_cases = sum(Y_reported),
      mean_ratio_local_to_import_favourable = mean(ratio_local_to_import[favourable], na.rm = TRUE),
      population_mean = mean(N), .groups = "drop"
    ) |>
    mutate(annual_reported_incidence_per_100k = annual_reported_cases / population_mean * 1e5,
           major_epidemic = annual_reported_cases >= MAJOR_ABS_CASE_MINIMUM & annual_reported_incidence_per_100k >= MAJOR_INCIDENCE_MINIMUM_PER_100K,
           outcome = ifelse(major_epidemic, "Major annual epidemic", "No major annual epidemic"), lambda_id = lambda_id)
  growth
}

message("[Phase4] Decomposing hazards and growth dynamics across the H1/G2 lambda grid (I0, I2, I3)...")
growth_all <- bind_rows(lapply(c("I0", "I2", "I3"), function(nm) decompose_arm(grid_H1[[nm]], nm)))
write_csv(growth_all, file.path(DIR_RESULTS_STEP7, "PHASE4_growth_dynamics_per_trajectory.csv"))

growth_summary <- growth_all |> group_by(lambda_id, outcome) |>
  summarise(n = n(), median_weeks_favourable_to_peak = median(weeks_favourable_to_peak, na.rm = TRUE),
            median_longest_increasing_run = median(longest_increasing_run),
            median_week_to_100 = median(week_to_100, na.rm = TRUE), prob_never_reach_100 = mean(is.na(week_to_100)),
            median_week_to_1000 = median(week_to_1000, na.rm = TRUE), prob_never_reach_1000 = mean(is.na(week_to_1000)),
            prob_extinct_after_favourable = mean(extinct_after_favourable, na.rm = TRUE),
            median_ratio_local_to_import_favourable = median(mean_ratio_local_to_import_favourable, na.rm = TRUE), .groups = "drop")
write_csv(growth_summary, file.path(DIR_TABLES_STEP7, "PHASE4_growth_dynamics_summary.csv"))
message("\n=== PHASE 4: growth dynamics, by lambda arm and 2022 outcome (H1/G2; DESCRIPTIVE) ==="); print(as.data.frame(growth_summary), digits = 4)

# Weekly hazard-ratio trajectory (median local/import hazard ratio among FAVOURABLE weeks, by lambda arm)
hazard_weekly <- bind_rows(lapply(c("I0", "I2", "I3"), function(nm) {
  grid_H1[[nm]]$trajectories |> mutate(import_hazard = M_import / N, local_hazard = pmax(lambda_t - M_import / N, 0), favourable = R_eff > 1, lambda_id = nm) |>
    filter(favourable) |> group_by(lambda_id, week_start) |>
    summarise(median_local_hazard = median(local_hazard), median_import_hazard = median(import_hazard), n = n(), .groups = "drop")
}))
write_csv(hazard_weekly, file.path(DIR_RESULTS_STEP7, "PHASE4_weekly_hazard_components_favourable_weeks.csv"))

message("\n[Phase4] Hazard decomposition and amplification analysis complete -> ", DIR_TABLES_STEP7)
