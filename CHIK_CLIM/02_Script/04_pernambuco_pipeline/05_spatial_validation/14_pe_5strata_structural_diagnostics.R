# Pernambuco Spatial Model 1 -- Sections 15-16: biological/structural
# diagnostics per stratum, and a municipality-level consistency check
# against the earlier crude case-implied attack-fraction diagnostic.

required_packages <- c("rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
pe_root <- file.path(root, "03_Output/model_fits/pernambuco/v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_model")
turnover_table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")

b <- readRDS(file.path(pe_root, "outputs/spatial5/pe_5strata_pilot.rds"))
fit <- b$fit
stratum_levels <- b$stratum_levels
R <- b$stan_data$R
draws <- rstan::extract(fit, permuted = TRUE)

weeks <- readRDS(file.path(turnover_table_dir, "pe_municipality_week_panel.rds")) |>
  pull(week_start) |> unique() |> sort()
waves <- read_csv(file.path(turnover_table_dir, "pe_wave_definitions.csv"), show_col_types = FALSE) |> arrange(start_week)
wave_contrib <- read_csv(file.path(table_dir, "PE_5strata_wave_contributions.csv"), show_col_types = FALSE)

message("=== Section 15: per-stratum biological/structural diagnostics ===")
struct_tbl <- bind_rows(lapply(seq_len(R), function(r) {
  S_prop_r <- draws$S_prop[, , r]
  R0_r <- draws$R0_t[, , r]
  U_prop_r <- draws$immune_prop[, , r]
  tibble(
    stratum = stratum_levels[r],
    min_S_prop_median = median(apply(S_prop_r, 1, min)),
    min_S_prop_lo95 = quantile(apply(S_prop_r, 1, min), .025),
    max_R0_median = median(apply(R0_r, 1, max)),
    max_R0_hi95 = quantile(apply(R0_r, 1, max), .975),
    cumulative_infected_2025_median = median(U_prop_r[, ncol(U_prop_r)]),
    cumulative_infected_2025_hi95 = quantile(U_prop_r[, ncol(U_prop_r)], .975)
  )
}))
print(as.data.frame(struct_tbl), digits = 3)

# Flag: implausible near-total infection (>0.95) or extreme R0 (>10) needed
# merely to reproduce later recurrence.
struct_tbl <- struct_tbl |> mutate(
  flag_near_total_infection = cumulative_infected_2025_hi95 > 0.95,
  flag_extreme_R0 = max_R0_hi95 > 10
)
message("\n=== Flags ===")
print(as.data.frame(struct_tbl |> select(stratum, flag_near_total_infection, flag_extreme_R0)))

# ---- Observed vs predicted wave totals and peak timing, by stratum --------
C_pred <- draws$C_pred # [iter, N, R]
wave_ppc <- bind_rows(lapply(seq_len(nrow(waves)), function(i) {
  w <- waves[i, ]
  idx <- which(weeks >= w$start_week & weeks <= w$end_week)
  bind_rows(lapply(seq_len(R), function(r) {
    obs <- b$stan_data$C[idx, r]
    pred_draws <- rowSums(C_pred[, idx, r, drop = FALSE])
    obs_peak_week <- weeks[idx][which.max(obs)]
    pred_peak_idx <- idx[apply(C_pred[, idx, r, drop = FALSE], 1, which.max)]
    tibble(wave_id = w$wave_id, stratum = stratum_levels[r],
           observed_total = sum(obs), predicted_median = median(pred_draws),
           predicted_lo95 = quantile(pred_draws, .025), predicted_hi95 = quantile(pred_draws, .975),
           observed_peak_week = obs_peak_week,
           predicted_peak_week_median = weeks[round(median(pred_peak_idx))])
  }))
}))
message("\n=== Observed vs predicted wave totals + peak timing, by stratum (first 15 rows) ===")
print(as.data.frame(head(wave_ppc, 15)), digits = 3)
write_csv(wave_ppc, file.path(table_dir, "PE_5strata_wave_ppc_by_stratum.csv"))
message("[saved] ", file.path(table_dir, "PE_5strata_wave_ppc_by_stratum.csv"))

write_csv(struct_tbl, file.path(table_dir, "PE_5strata_structural_diagnostics.csv"))
message("\n[saved] ", file.path(table_dir, "PE_5strata_structural_diagnostics.csv"))

# ---- Section 16: municipality consistency check ----------------------------
# Compare posterior regional cumulative-infection burden with the earlier
# crude municipality-level case-implied attack fraction (spatial-turnover
# diagnostic script 27), aggregated to the SAME 5 strata using the SAME
# stratum lookup. NOT a re-fit -- q_municipality is NOT estimated here.
implied <- read_csv(file.path(turnover_table_dir, "pe_municipality_implied_attack_by_q.csv"), show_col_types = FALSE,
                     col_types = cols(muni6 = col_character()))
stratum_lookup <- read_csv(file.path(table_dir, "PE_spatial_stratum_lookup.csv"), show_col_types = FALSE,
                            col_types = cols(muni6 = col_character()))

crude_by_stratum <- implied |>
  dplyr::filter(checkpoint == as.Date("2025-12-21"), q_label == "global-q median") |>
  left_join(stratum_lookup |> select(muni6, final_model_stratum), by = "muni6") |>
  dplyr::filter(!is.na(final_model_stratum)) |>
  group_by(final_model_stratum) |>
  summarise(crude_pop_weighted_A = weighted.mean(A_i, w = population_at_checkpoint, na.rm = TRUE),
            crude_median_A = median(A_i), crude_max_A = max(A_i), .groups = "drop") |>
  rename(stratum = final_model_stratum)

consistency_tbl <- struct_tbl |>
  select(stratum, model_cumulative_infected_2025_median = cumulative_infected_2025_median,
         model_cumulative_infected_2025_hi95 = cumulative_infected_2025_hi95) |>
  left_join(crude_by_stratum, by = "stratum") |>
  mutate(ratio_model_over_crude = model_cumulative_infected_2025_median / crude_pop_weighted_A)

message("\n=== Section 16: model regional burden vs municipality-level crude implied attack (q=homogeneous global-q median, for reference) ===")
print(as.data.frame(consistency_tbl), digits = 3)
write_csv(consistency_tbl, file.path(table_dir, "PE_5strata_municipality_consistency_check.csv"))
message("[saved] ", file.path(table_dir, "PE_5strata_municipality_consistency_check.csv"))
message("\nNote: crude comparison uses the HOMOGENEOUS model's q (0.0098) since the crude municipality diagnostic was built before this spatial model existed. No municipality-specific q is fit here (Section 16 prohibits this).")
