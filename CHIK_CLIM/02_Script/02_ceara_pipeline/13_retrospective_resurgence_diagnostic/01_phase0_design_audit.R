# PHASE 0 -- Verify the existing H1/H2 retrospective design for leakage
# or implementation errors BEFORE any new mechanism diagnosis. Read-only
# code audit (10 checklist items) + one quantitative check: whether the
# climate anomaly standardisation (climatology mean/SD) computed over
# the FULL 2015-2025 period leaks holdout-period information into the
# "observed climate" retrospective inputs used by simulate_holdout().

required_packages <- c("here", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/13_retrospective_resurgence_diagnostic/00_config_retro_diagnostic.R"))

# ================================================================
# Checklist 1-10: code-audit findings (verified by reading the actual
# source referenced below; PASS/PASS-WITH-CAVEAT/FAIL).
# ================================================================
audit <- tribble(
  ~item, ~check, ~status, ~evidence,
  1, "Training-only parameter estimation",
    "PASS", "04_phase3_fit_holdout_training.R fits ce_tau_year_training_{H1,H2}.rds using ONLY build_training_only_tau_data(train_start,train_end) -- the ACCEPTED_CLIMATE_FIT_M0 / P1-P3 full-period fits are never read by the training step.",
  2, "No holdout cases used during fitting",
    "PASS", "build_training_only_tau_data() (03_build_training_only_data.R) filters week_start<=train_end BEFORE calling make_v4_3_data(); fit_index=seq_len(n) covers only training weeks; N_fit=n (training length).",
  3, "No 2022 observed-data-conditioned seed in holdout",
    "PASS", "build_training_only_tau_data() sets is_seed=integer(n) (all 0) and X_seed=numeric(n) (all 0) by construction -- never calls the accepted seed-requiring make_v4_8/v4_9 builder. simulate_holdout() (12_future_year_effect_validation/06_phase3_4_retrospective_validation.R L121) passes NO is_seed/X_seed into simulate_renewal_core() for the holdout period either.",
  4, "No full-period posterior draws used in holdout",
    "PASS", "simulate_holdout() L64-65 extracts draws exclusively from training_bundle$fit (ce_tau_year_training_H1/H2.rds); ACCEPTED_CLIMATE_FIT_M0 and ce_tau_year_{P1,P2,P3}.rds are not referenced anywhere in 06_phase3_4_retrospective_validation.R.",
  5, "2018 serology correctly included",
    "PASS", "build_training_only_tau_data() L29-32 computes find_sero_week_index() and explicitly stop()s if the serology week postdates train_end; both H1 (cutoff 2019-12-29) and H2 (cutoff 2021-12-26) pass (serology midpoint 2018-09-16 precedes both).",
  6, "Climate climatology/scaling uses only admissible historical information",
    "PASS-WITH-CAVEAT", "01_build_climate_anomaly_covariates.R computes climatology (T_clim/logP_clim by week_of_year) and the standardisation scale (T_anom_mean/sd, P_anom_mean/sd) by pooling ACROSS THE FULL 2015-2025 period with NO date filter -- i.e. the holdout years' own weather contributes to the long-run mean/SD used to standardise ALL years including the training period. Quantified below: this is a genuine but small normalisation-level leakage, strictly distinct from (and much smaller than) using holdout CASES. The retrospective arm is already explicitly labelled 'observed-climate-conditional hindcasting, not prospective forecasting' for the separate, larger reason that it uses the holdout period's ACTUAL recorded weather (unavoidable for any retrospective exercise) -- this item concerns only the standardisation CONSTANTS, not the weather values themselves.",
  7, "Correct transfer of S/U at training cutoff",
    "PASS", "advance_terminal_state() (03_future_renewal_functions.R L115-120) completes the one demographic update the Stan model's `if (t<N)` guard skips at the final fitted week; identical function used for both the original future-pilot horizon and simulate_holdout(); unit-tested in 10_future_projection_pilot/04_run_unit_tests.R (test 2.1, not rerun here).",
  8, "Correct continuation of seasonal phase",
    "PASS", "build_training_only_tau_data() builds seasonal_sin/cos/sin2/cos2 via seq_len(n) starting at t=1 on the SAME train_start=2015-01-04 as the full/historical builder; simulate_holdout() L70 continues t_index_holdout=(N_train+1):(N_train+n_holdout) -- phase is continuous across the training/holdout boundary by construction (both series share the same absolute week-1 anchor).",
  9, "No accidental recentering/rescaling of future year effects",
    "PASS", "simulate_holdout() L106-108 draws z_year_h over n_holdout_years ONLY (not training years) and recentres year_effect_h = tau*(z_year_h-mean(z_year_h)) within that holdout-years set alone -- structurally identical to how the historical fit recentres year_effect within its own fitted years; no cross-contamination of the two year sets.",
  10, "Correct stochastic importation and observation-model implementation",
    "PASS", "simulate_renewal_core() (03_future_renewal_functions.R L211-213) adds M_t~Poisson(lambda_import) as a hazard addition (lambda_t = R0*infectiousness/N + M_t/N), never double-counted as a case; simulate_reported_cases() applies Y~NegBinomial(mu=q*X,phi_obs) with q=0.05 fixed, phi_obs drawn from the SAME training-only posterior -- consistent with the accepted model's observation layer."
)
write_csv(audit, file.path(DIR_TABLES_RD, "PHASE0_design_audit.csv"))
message("\n=== PHASE 0 design audit (10-point checklist) ===")
print(as.data.frame(audit[, c("item","check","status")]), row.names = FALSE)

# ================================================================
# Quantify item 6: full-period vs training-only climatology/scaling,
# and how much it would shift the actual z_T_anom/z_P_anom values used
# as "observed climate" inputs during each holdout.
# ================================================================
PRIMARY_STATES <- c("CE", "BA", "RJ", "MT") # pooled standardisation set, byte-identical to 01_build_climate_anomaly_covariates.R
raw_path <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_lagged_raw.csv")
scale_path_full <- file.path(ROOT, "03_Output/02_ceara_pipeline/tables/climate_forced_v4_9/climate_anomaly_scaling_metadata.csv")
lagged <- read_csv(raw_path, show_col_types = FALSE) |> dplyr::filter(lag_spec == "primary")
scaling_full <- read_csv(scale_path_full, show_col_types = FALSE) |> dplyr::filter(lag_spec == "primary")

recompute_scaling <- function(cutoff) {
  lag_cut <- lagged |> dplyr::filter(week_start <= cutoff)
  clim_cut <- lag_cut |> dplyr::filter(!is.na(lagged_temp), !is.na(log_precip)) |>
    group_by(state, week_of_year) |> summarise(T_clim = mean(lagged_temp), logP_clim = mean(log_precip), .groups = "drop")
  anom_cut <- lag_cut |> left_join(clim_cut, by = c("state", "week_of_year")) |>
    mutate(T_anom = lagged_temp - T_clim, P_anom = log_precip - logP_clim)
  anom_cut |> dplyr::filter(state %in% PRIMARY_STATES, !is.na(T_anom), !is.na(P_anom)) |>
    summarise(T_anom_mean = mean(T_anom), T_anom_sd = sd(T_anom), P_anom_mean = mean(P_anom), P_anom_sd = sd(P_anom), n_obs_pooled = n())
}

scaling_H1 <- recompute_scaling(HOLDOUTS$H1$train_end) |> mutate(spec = "H1_training_only_le_2019_12_29")
scaling_H2 <- recompute_scaling(HOLDOUTS$H2$train_end) |> mutate(spec = "H2_training_only_le_2021_12_26")
scaling_cmp <- bind_rows(
  scaling_full |> select(T_anom_mean, T_anom_sd, P_anom_mean, P_anom_sd, n_obs_pooled = n_obs_pooled) |> mutate(spec = "FULL_2015_2025_as_used_in_simulate_holdout"),
  scaling_H1, scaling_H2
) |> relocate(spec)
write_csv(scaling_cmp, file.path(DIR_TABLES_RD, "PHASE0_climatology_leakage_check.csv"))

pct_shift <- scaling_cmp |> dplyr::filter(spec != "FULL_2015_2025_as_used_in_simulate_holdout") |>
  mutate(pct_diff_T_sd = 100 * (T_anom_sd - scaling_full$T_anom_sd) / scaling_full$T_anom_sd,
         pct_diff_P_sd = 100 * (P_anom_sd - scaling_full$P_anom_sd) / scaling_full$P_anom_sd,
         pct_diff_T_mean_in_T_sd_units = (T_anom_mean - scaling_full$T_anom_mean) / scaling_full$T_anom_sd,
         pct_diff_P_mean_in_P_sd_units = (P_anom_mean - scaling_full$P_anom_mean) / scaling_full$P_anom_sd)
write_csv(pct_shift, file.path(DIR_TABLES_RD, "PHASE0_climatology_leakage_magnitude.csv"))

message("\n=== PHASE 0 item 6 quantification: full-period vs training-only climate standardisation ===")
print(as.data.frame(scaling_cmp), digits = 4)
message("\nPercent shift in SD, and mean shift expressed in full-period SD units (i.e. how many z-score units a given week's anomaly would move if re-standardised training-only):")
print(as.data.frame(pct_shift[, c("spec","pct_diff_T_sd","pct_diff_P_sd","pct_diff_T_mean_in_T_sd_units","pct_diff_P_mean_in_P_sd_units")]), digits = 3)

gate_status <- tibble(
  check = "PHASE0_gate",
  all_items_pass_or_caveat = all(audit$status %in% c("PASS", "PASS-WITH-CAVEAT")),
  serious_leakage_requiring_fix_before_proceeding = FALSE,
  note = "Item 6 (climate standardisation) carries a disclosed, quantified normalisation-level leakage, LARGER than initially assumed: the full-period (2015-2025) pooled SD used to build z_T_anom/z_P_anom is 13.1%/9.5% larger than a true training-only SD for H1 and 7.6%/5.2% larger for H2 (PHASE0_climatology_leakage_magnitude.csv) -- i.e. the full period's own more extreme years (including the 2022 resurgence period itself) inflate the long-run SD used to standardise EVERY week, training and holdout alike. Direction matters: a LARGER divisor SD means every week's |z_T_anom|/|z_P_anom| -- training and holdout -- is systematically SMALLER in magnitude than it would be under a true leak-free training-only standardisation. Since beta_T/beta_P were themselves estimated (in the training-only Stan fit) against these same deflated z-scores, the fitted coefficients and the holdout's climate inputs are at least INTERNALLY CONSISTENT (same standardisation applied on both sides of the training/holdout boundary) -- the leak does not let the model see holdout-magnitude information it could not otherwise infer. If anything, a smaller-than-true |z| would UNDERSTATE rather than inflate climate-driven transmission risk, so this leakage cannot explain H1's under-prediction of 2022 and works, if at all, in the opposite direction. NOT fixed here: correcting it would require rebuilding climate_anomaly_covariates.csv with a training-only standardisation and refitting the accepted M0 climate model and all three training-only tau_year fits -- out of scope under this task's 'NO NEW STAN FITTING' / 'do not modify accepted models' constraints. Disclosed, not corrected, and flagged as a candidate root-cause check for future work (Phase 4/Unresolved issues)."
)
write_csv(gate_status, file.path(DIR_TABLES_RD, "PHASE0_gate_status.csv"))
message("\n[Phase0] Design audit complete -> ", DIR_TABLES_RD)
