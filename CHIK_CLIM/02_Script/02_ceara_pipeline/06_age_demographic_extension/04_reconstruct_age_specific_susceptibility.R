# Age-structured vaccine extension -- dedicated diagnostic: age-specific
# susceptible fraction / infection incidence under a COMMON force of
# infection (no age-specific transmission, no vaccination). Runs the
# already-validated (Phase 4 identity-pass) age-cohort simulator for many
# posterior draws and characterises how S_prop_age evolves by age and
# time purely through births/ageing/differential cumulative exposure.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/model_fits/ceara/climate_forced_v4_9")
table_dir <- file.path(root, "03_Output/tables/climate_forced_v4_9")
figure_dir <- file.path(root, "03_Output/figures/climate_forced_v4_9")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
source(file.path(root, "02_Script/02_ceara_pipeline/06_age_demographic_extension/02_age_cohort_simulator.R")) # moved from 36/.../07_age_cohort_simulator.R

MAX_AGE <- 100L; A <- MAX_AGE + 1L
ages <- 0:MAX_AGE

# ================================================================
# Setup (identical inputs to Phase 4)
# ================================================================
climate <- readRDS(file.path(base_dir, "outputs/ce_canary/ce_climate_forced_canary_q0.05.rds"))
fit <- climate$fit
weekly <- climate$weekly_data
dates <- as.Date(weekly$week_start)
years <- as.integer(format(dates, "%Y"))
N <- length(dates)
age_shares <- read_csv(file.path(table_dir, "CE_age_population_shares.csv"), show_col_types = FALSE)

sd_data <- climate$stan_data
w <- sd_data$w; imports_per_week <- sd_data$imports_per_week
is_seed <- sd_data$is_seed; X_seed <- sd_data$X_seed
births <- weekly$births; deaths <- weekly$all_cause_deaths
reconciliation <- weekly$net_population_reconciliation; N_start <- weekly$N_start

# Actual CE wave dates (national wave census) -- NOT guessed.
CE_2017_ONSET <- as.Date("2017-02-12"); CE_2017_END <- as.Date("2018-12-16") # CE_wave_03, 106,695 cases
CE_2022_SEED_START <- min(climate$seed_dates); CE_2022_WAVE_END <- as.Date("2023-12-17") # CE_wave_06

AGE_GROUP_BREAKS <- c(0, 5, 10, 20, 40, 60, MAX_AGE + 1)
AGE_GROUP_LABELS <- c("0-4", "5-9", "10-19", "20-39", "40-59", "60+")
age_group <- cut(ages, breaks = AGE_GROUP_BREAKS, labels = AGE_GROUP_LABELS, right = FALSE)

draws_v49 <- rstan::extract(fit, pars = c("R0_t", "S", "U", "X"), permuted = TRUE)
n_draws_total <- nrow(draws_v49$R0_t)
N_DRAWS <- 80L
set.seed(20260918L)
test_draw_idx <- sample.int(n_draws_total, N_DRAWS)
message(sprintf("Running age-specific diagnostic for %d / %d posterior draws...", N_DRAWS, n_draws_total))

# ================================================================
# Run simulator for all draws, accumulate single-year S_prop/immune_prop
# arrays [draw, t, age] for the heatmap, and age-GROUP absolute S/Uinf/N/Inew
# arrays [draw, t, group] for trajectories/incidence.
# ================================================================
S_prop_arr <- array(NA_real_, dim = c(N_DRAWS, N, A))
immune_prop_arr <- array(NA_real_, dim = c(N_DRAWS, N, A))
S_group_arr <- array(NA_real_, dim = c(N_DRAWS, N, length(AGE_GROUP_LABELS)))
Uinf_group_arr <- array(NA_real_, dim = c(N_DRAWS, N, length(AGE_GROUP_LABELS)))
N_group_arr <- array(NA_real_, dim = c(N_DRAWS, N, length(AGE_GROUP_LABELS)))
Inew_group_arr <- array(NA_real_, dim = c(N_DRAWS, N, length(AGE_GROUP_LABELS)))
hazard_mat <- matrix(NA_real_, N_DRAWS, N)
validation_rows <- vector("list", N_DRAWS)

# Cross-section checkpoints (Section 6) -- actual model dates
checkpoint_dates <- as.Date(c("2015-01-04", "2017-02-05", "2018-12-16", "2021-12-26", "2023-12-17", "2025-12-21"))
checkpoint_labels <- c("Baseline (start)", "Before 2017 epidemic", "After 2017 epidemic",
                        "Before 2022 seed/epidemic", "After 2022 epidemic", "End of simulation")
checkpoint_idx <- sapply(checkpoint_dates, function(cd) which.min(abs(dates - cd)))
S_prop_checkpoint_arr <- array(NA_real_, dim = c(N_DRAWS, length(checkpoint_dates), A))
immune_prop_checkpoint_arr <- array(NA_real_, dim = c(N_DRAWS, length(checkpoint_dates), A))

for (i in seq_len(N_DRAWS)) {
  d <- test_draw_idx[i]
  sim <- simulate_age_cohort(
    R0_t_draw = draws_v49$R0_t[d, ], births = births, deaths = deaths, reconciliation = reconciliation,
    N_start = N_start, w = w, imports_per_week = imports_per_week, is_seed = is_seed, X_seed = X_seed,
    age_shares = age_shares, years = years, max_age = MAX_AGE, vaccination = NULL
  )
  N_age <- sim$S_age + sim$Uinf_age + sim$Uvac_age
  S_prop_arr[i, , ] <- sim$S_age / N_age
  immune_prop_arr[i, , ] <- sim$Uinf_age / N_age
  hazard_mat[i, ] <- sim$hazard

  for (g in seq_along(AGE_GROUP_LABELS)) {
    cols <- which(age_group == AGE_GROUP_LABELS[g])
    S_group_arr[i, , g] <- rowSums(sim$S_age[, cols, drop = FALSE])
    Uinf_group_arr[i, , g] <- rowSums(sim$Uinf_age[, cols, drop = FALSE])
    N_group_arr[i, , g] <- rowSums(N_age[, cols, drop = FALSE])
    Inew_group_arr[i, , g] <- rowSums(sim$Inew_age[, cols, drop = FALSE])
  }
  S_prop_checkpoint_arr[i, , ] <- sim$S_age[checkpoint_idx, ] / N_age[checkpoint_idx, ]
  immune_prop_checkpoint_arr[i, , ] <- sim$Uinf_age[checkpoint_idx, ] / N_age[checkpoint_idx, ]

  # Section 8 critical validation (re-confirmed here, not assumed)
  validation_rows[[i]] <- tibble(
    draw = d,
    max_abs_diff_S = max(abs(rowSums(sim$S_age) - draws_v49$S[d, ])),
    max_abs_diff_U = max(abs(rowSums(sim$Uinf_age) - draws_v49$U[d, ])),
    max_abs_diff_X = max(abs(sim$X_total - draws_v49$X[d, ])),
    max_abs_diff_N = max(abs(rowSums(N_age) - N_start))
  )
  if (i %% 20 == 0) message(sprintf("  ...%d / %d draws done", i, N_DRAWS))
}
validation_tbl <- bind_rows(validation_rows)
message("\n=== Section 8: critical validation (re-confirmed at single-year age resolution) ===")
print(as.data.frame(summarise(validation_tbl, across(starts_with("max_abs_diff"), max))), digits = 4)
write_csv(validation_tbl, file.path(table_dir, "AGE_DIAGNOSTIC_validation.csv"))
identity_still_holds <- all(validation_tbl$max_abs_diff_S < 1e-3, validation_tbl$max_abs_diff_U < 1e-3,
                             validation_tbl$max_abs_diff_X < 1e-3, validation_tbl$max_abs_diff_N < 1e-3)
message(sprintf("Aggregate identity still holds at single-year resolution: %s", identity_still_holds))
if (!identity_still_holds) stop("STOP: aggregate identity broke during age-specific diagnostic computation -- do not proceed.")

# ---- Explicit numerical check: p_infection_given_susceptible identical across ages ----
message("\n=== Section 1: verifying p_infection_given_susceptible is identical across ages ===")
spotcheck_t <- c(50, 120, 165, 400, 550) # avoid seed weeks for this specific check
spotcheck_rows <- lapply(spotcheck_t, function(t_idx) {
  d <- test_draw_idx[1]
  sim <- simulate_age_cohort(draws_v49$R0_t[d, ], births, deaths, reconciliation, N_start, w, imports_per_week,
                              is_seed, X_seed, age_shares, years, MAX_AGE, vaccination = NULL)
  ratio_by_age <- sim$Inew_age[t_idx, ] / pmax(sim$S_age[t_idx, ], 1e-9)
  tibble(t = t_idx, week_start = dates[t_idx], hazard = sim$hazard[t_idx],
         max_abs_dev_across_ages = max(abs(ratio_by_age - sim$hazard[t_idx])))
})
spotcheck_tbl <- bind_rows(spotcheck_rows)
print(as.data.frame(spotcheck_tbl), digits = 10)
write_csv(spotcheck_tbl, file.path(table_dir, "AGE_DIAGNOSTIC_hazard_uniformity_check.csv"))
message(sprintf("Max deviation of age-specific infection ratio from the common hazard, across %d spot-check weeks: %.2e (should be ~0)",
                 length(spotcheck_t), max(spotcheck_tbl$max_abs_dev_across_ages)))

message("\n[saved] AGE_DIAGNOSTIC_validation.csv, AGE_DIAGNOSTIC_hazard_uniformity_check.csv")
saveRDS(list(S_prop_arr = S_prop_arr, immune_prop_arr = immune_prop_arr,
             S_group_arr = S_group_arr, Uinf_group_arr = Uinf_group_arr, N_group_arr = N_group_arr, Inew_group_arr = Inew_group_arr,
             S_prop_checkpoint_arr = S_prop_checkpoint_arr, immune_prop_checkpoint_arr = immune_prop_checkpoint_arr,
             checkpoint_dates = checkpoint_dates, checkpoint_labels = checkpoint_labels, checkpoint_idx = checkpoint_idx,
             dates = dates, ages = ages, age_group = age_group, age_group_labels = AGE_GROUP_LABELS, test_draw_idx = test_draw_idx),
        file.path(base_dir, "outputs/ce_canary/age_diagnostic_arrays.rds"))
message("[saved] age_diagnostic_arrays.rds (intermediate, for plotting script)")
