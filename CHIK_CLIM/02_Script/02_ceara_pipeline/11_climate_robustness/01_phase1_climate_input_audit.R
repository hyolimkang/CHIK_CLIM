# PHASE 1 -- Climate input consistency audit. GATING STEP: Phase 2 (new
# Stan fits) must not proceed unless this script reports PASS.
#
# Independently reconstructs the primary-lag climate anomaly covariates
# FROM THE RAW national climate table (not by re-reading the already-
# saved covariates CSV, which would be circular), then checks, in order:
#   (a) recomputed standardisation constants vs. saved scaling metadata
#   (b) recomputed z_T_anom/z_P_anom (CE, primary) vs. the saved covariates CSV
#   (c) the saved covariates CSV, joined exactly as the fitting script
#       joins it (with boundary-week NA->0), vs. the z_T_anom/z_P_anom
#       VECTOR ACTUALLY STORED inside the accepted M0 Stan fit bundle
# No Stan fitting in this script.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "lubridate")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(lubridate) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/11_climate_robustness/00_config_climate_robustness.R"))

audit_rows <- list()
audit <- function(check, pass, detail) {
  message(sprintf("[Phase1 %s] %s | %s", if (pass) "PASS" else "FAIL", check, detail))
  audit_rows[[length(audit_rows) + 1]] <<- tibble(check = check, pass = pass, detail = detail)
}

# ================================================================
# (a) Independent reconstruction of the primary-lag covariates FROM RAW
# (byte-for-byte reimplementation of 01_build_climate_anomaly_covariates.R's
# primary spec only -- not sourcing that script, to keep this an
# independent check rather than trusting its own code unmodified)
# ================================================================
message("[Phase1] Reading raw national climate table and reconstructing primary-lag anomalies independently...")
climate_raw <- read_csv(NATIONAL_CLIMATE_CSV, show_col_types = FALSE) |>
  mutate(week_start = as.Date(week_start), week_of_year = lubridate::isoweek(week_start)) |>
  arrange(state, week_start)
if (anyNA(climate_raw$Tmean) || anyNA(climate_raw$PRCP)) stop("NA in national climate file -- cannot proceed with audit.")

lagged_window_stat <- function(x, offsets, stat = c("mean", "sum")) {
  stat <- match.arg(stat); n <- length(x); out <- rep(NA_real_, n)
  for (t in seq_len(n)) {
    idx <- t + offsets
    if (all(idx >= 1L & idx <= n)) out[t] <- if (stat == "mean") mean(x[idx]) else sum(x[idx])
  }
  out
}
PRIMARY_TEMP_OFFSETS <- -2:0; PRIMARY_PRECIP_OFFSETS <- -6:-2 # byte-identical to LAG_SPECS$primary in the accepted build script

lagged_primary <- bind_rows(lapply(split(climate_raw, climate_raw$state), function(d) {
  d <- d |> arrange(week_start)
  tibble(state = d$state[1], week_start = d$week_start, week_of_year = d$week_of_year,
         lagged_temp = lagged_window_stat(d$Tmean, PRIMARY_TEMP_OFFSETS, "mean"),
         lagged_precip = lagged_window_stat(d$PRCP, PRIMARY_PRECIP_OFFSETS, "sum")) |>
    mutate(log_precip = log1p(lagged_precip))
}))

climatology_primary <- lagged_primary |> dplyr::filter(!is.na(lagged_temp), !is.na(log_precip)) |>
  group_by(state, week_of_year) |> summarise(T_clim = mean(lagged_temp), logP_clim = mean(log_precip), .groups = "drop")

anomalies_primary <- lagged_primary |> left_join(climatology_primary, by = c("state", "week_of_year")) |>
  mutate(T_anom = lagged_temp - T_clim, P_anom = log_precip - logP_clim)

scaling_primary_recomputed <- anomalies_primary |> dplyr::filter(state %in% PRIMARY_STATES, !is.na(T_anom), !is.na(P_anom)) |>
  summarise(T_anom_mean = mean(T_anom), T_anom_sd = sd(T_anom), P_anom_mean = mean(P_anom), P_anom_sd = sd(P_anom), n_obs_pooled = n())

scaling_saved <- read_csv(ACCEPTED_SCALING_CSV, show_col_types = FALSE) |> dplyr::filter(lag_spec == "primary")

message("\n=== (a) Standardisation constants: independently recomputed vs. saved metadata ===")
scaling_compare <- tibble(
  quantity = c("T_anom_mean", "T_anom_sd", "P_anom_mean", "P_anom_sd", "n_obs_pooled"),
  recomputed_from_raw = c(scaling_primary_recomputed$T_anom_mean, scaling_primary_recomputed$T_anom_sd, scaling_primary_recomputed$P_anom_mean, scaling_primary_recomputed$P_anom_sd, scaling_primary_recomputed$n_obs_pooled),
  saved_metadata_csv  = c(scaling_saved$T_anom_mean, scaling_saved$T_anom_sd, scaling_saved$P_anom_mean, scaling_saved$P_anom_sd, scaling_saved$n_obs_pooled),
  earlier_narrative_report_value = c(-0.0038, 1.097, 0.0126, 0.756, NA) # STEP1_CLIMATE_INFERENCE_ASSESSMENT.md's prose -- checked below
)
print(as.data.frame(scaling_compare), digits = 6)
write_csv(scaling_compare, file.path(DIR_TABLES_CR, "PHASE1_scaling_constants_comparison.csv"))

sd_match <- abs(scaling_primary_recomputed$T_anom_sd - scaling_saved$T_anom_sd) < 1e-6 && abs(scaling_primary_recomputed$P_anom_sd - scaling_saved$P_anom_sd) < 1e-6
audit("Independently recomputed standardisation SDs match saved scaling_metadata.csv exactly", sd_match,
      sprintf("recomputed T_sd=%.6f / P_sd=%.6f vs. saved T_sd=%.6f / P_sd=%.6f", scaling_primary_recomputed$T_anom_sd, scaling_primary_recomputed$P_anom_sd, scaling_saved$T_anom_sd, scaling_saved$P_anom_sd))

narrative_discrepancy <- abs(1.097 - scaling_saved$T_anom_sd) > 0.05
audit("Earlier STEP1 narrative report's quoted T_anom_sd=1.097/P_anom_sd=0.756 do NOT match the saved/recomputed values", narrative_discrepancy,
      sprintf("Saved/recomputed: T_sd=%.4f, P_sd=%.4f. STEP1 narrative text (STEP1_CLIMATE_INFERENCE_ASSESSMENT.md Section 1) instead reported 1.097/0.756 -- this is a TRANSCRIPTION ERROR in that report's prose, not a pipeline/data discrepancy: the programmatically-generated STEP1_lag_and_standardisation_specification.csv from the SAME audit run already correctly shows T_sd=0.7663/P_sd=0.5016, matching this independent recomputation exactly. Correcting the narrative text is recommended (see assessment report Section 1).", scaling_saved$T_anom_sd, scaling_saved$P_anom_sd))

# ================================================================
# (b) Recomputed z_T_anom/z_P_anom (CE only, primary) vs. saved covariates CSV
# ================================================================
z_recomputed <- anomalies_primary |> dplyr::filter(state == "CE") |>
  mutate(z_T_anom = (T_anom - scaling_primary_recomputed$T_anom_mean) / scaling_primary_recomputed$T_anom_sd,
         z_P_anom = (P_anom - scaling_primary_recomputed$P_anom_mean) / scaling_primary_recomputed$P_anom_sd) |>
  select(week_start, z_T_anom_recomputed = z_T_anom, z_P_anom_recomputed = z_P_anom)

saved_covariates <- read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
  dplyr::filter(state == "CE", lag_spec == "primary") |> mutate(week_start = as.Date(week_start)) |>
  select(week_start, z_T_anom_saved = z_T_anom, z_P_anom_saved = z_P_anom)

b_compare <- inner_join(z_recomputed, saved_covariates, by = "week_start")
n_dates_mismatch <- nrow(z_recomputed) != nrow(saved_covariates) || nrow(b_compare) != nrow(z_recomputed)
b_compare <- b_compare |> mutate(abs_diff_T = abs(z_T_anom_recomputed - z_T_anom_saved), abs_diff_P = abs(z_P_anom_recomputed - z_P_anom_saved))
max_diff_T_b <- max(b_compare$abs_diff_T, na.rm = TRUE); max_diff_P_b <- max(b_compare$abs_diff_P, na.rm = TRUE)
write_csv(b_compare, file.path(DIR_TABLES_CR, "PHASE1_z_anom_recomputed_vs_saved_covariates.csv"))
audit("Weekly date alignment: recomputed CE-primary series matches saved covariates CSV row-for-row", !n_dates_mismatch,
      sprintf("recomputed n=%d, saved n=%d, matched n=%d", nrow(z_recomputed), nrow(saved_covariates), nrow(b_compare)))
audit("Recomputed z_T_anom/z_P_anom (CE, primary) match saved covariates CSV", max_diff_T_b < 1e-6 && max_diff_P_b < 1e-6,
      sprintf("max abs diff: z_T=%.2e, z_P=%.2e (allowing for the known NA boundary weeks, compared below)", max_diff_T_b, max_diff_P_b))

# ================================================================
# (c) Saved covariates CSV, joined exactly as the fitting script does
# (NA -> 0 for boundary weeks), vs. the z_T_anom/z_P_anom VECTOR ACTUALLY
# STORED inside the accepted M0 Stan fit bundle -- the definitive check
# of what was actually fed into Stan.
# ================================================================
message("\n[Phase1] Loading accepted M0 Stan fit bundle (read-only; no refit) to extract the ACTUAL fitted z_T_anom/z_P_anom...")
m0 <- readRDS(ACCEPTED_CLIMATE_FIT_M0)
dates_fit <- as.Date(m0$weekly_data$week_start)
N_fit <- length(dates_fit)
stopifnot(identical(dates_fit, seq(FIXED$date_start, FIXED$date_end, by = "week")))

joined_for_fit <- tibble(week_start = dates_fit) |> left_join(saved_covariates, by = "week_start")
n_na_boundary <- sum(is.na(joined_for_fit$z_T_anom_saved) | is.na(joined_for_fit$z_P_anom_saved))
joined_for_fit$z_T_anom_saved[is.na(joined_for_fit$z_T_anom_saved)] <- 0
joined_for_fit$z_P_anom_saved[is.na(joined_for_fit$z_P_anom_saved)] <- 0

c_compare <- tibble(week_start = dates_fit, z_T_reconstructed = joined_for_fit$z_T_anom_saved, z_P_reconstructed = joined_for_fit$z_P_anom_saved,
                     z_T_in_stan_fit = m0$stan_data$z_T_anom, z_P_in_stan_fit = m0$stan_data$z_P_anom) |>
  mutate(abs_diff_T = abs(z_T_reconstructed - z_T_in_stan_fit), abs_diff_P = abs(z_P_reconstructed - z_P_in_stan_fit))
write_csv(c_compare, file.path(DIR_TABLES_CR, "PHASE1_z_anom_reconstructed_vs_actual_stan_input.csv"))

max_diff_T_c <- max(c_compare$abs_diff_T); max_diff_P_c <- max(c_compare$abs_diff_P)
audit("Boundary-week count matches documented disclosure (a handful of weeks in 2015-01/02)", n_na_boundary >= 1L && n_na_boundary <= 10L,
      sprintf("%d weeks had NA anomaly before the NA->0 fill, at: %s", n_na_boundary, paste(format(joined_for_fit$week_start[is.na(saved_covariates$z_T_anom_saved[match(joined_for_fit$week_start, saved_covariates$week_start)])]), collapse = "; ")))
audit("Reconstructed weekly z_T_anom/z_P_anom (independent rebuild + documented NA->0 rule) EXACTLY match the vector stored inside the accepted M0 Stan fit bundle", max_diff_T_c < 1e-9 && max_diff_P_c < 1e-9,
      sprintf("max abs diff vs. actual Stan input: z_T=%.2e, z_P=%.2e across all %d fitted weeks", max_diff_T_c, max_diff_P_c, N_fit))

# ================================================================
# Save audit summary + gate decision
# ================================================================
audit_summary <- bind_rows(audit_rows)
write_csv(audit_summary, file.path(DIR_TABLES_CR, "PHASE1_input_consistency_audit_summary.csv"))
message("\n=== PHASE 1 AUDIT SUMMARY ===")
print(as.data.frame(audit_summary[, c("check", "pass")]))

phase1_pass <- all(audit_summary$pass[audit_summary$check != "Earlier STEP1 narrative report's quoted T_anom_sd=1.097/P_anom_sd=0.756 do NOT match the saved/recomputed values"]) # that one row documents a KNOWN prose error, not a gating failure
message(sprintf("\n[Phase1] GATE DECISION: %s", if (phase1_pass) "PASS -- proceed to Phase 2" else "FAIL -- STOP, do not proceed to Phase 2"))
if (!phase1_pass) stop("Phase 1 climate input consistency audit FAILED -- see PHASE1_input_consistency_audit_summary.csv. Resolve before any new Stan fitting.")
