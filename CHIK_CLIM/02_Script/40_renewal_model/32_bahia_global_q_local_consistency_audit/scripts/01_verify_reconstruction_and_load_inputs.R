# Bahia global-q local-consistency audit -- Section 1 & 2.
#
# Post-hoc plausibility/consistency audit only. Does NOT modify Stan, does
# NOT refit q, does NOT estimate municipality-specific latent infections.
# Loads and cross-checks all validated inputs before any arithmetic.

required_packages <- c("rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/bahia_global_q_local_consistency_audit")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# ---- A/B: municipality-week panel (reuse the already-validated one from ----
#           the spatial-turnover diagnostic -- same source, same window)
turnover_table_dir <- file.path(root, "03_Output/tables/bahia_spatial_turnover_diagnostic")
bahia_panel <- readRDS(file.path(turnover_table_dir, "bahia_municipality_week_panel.rds"))

# ---- D: the current full Bahia global-q + geo-adjusted-serology fit --------
fit_bundle <- readRDS(file.path(root, "02_Script/40_renewal_model/30_bahia_v4_9_replication",
                                 "outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds"))
fit <- fit_bundle$fit
weekly_data <- fit_bundle$weekly_data

# ---- Section 2: verify municipality -> state reconstruction ---------------
state_from_muni <- bahia_panel |>
  group_by(week_start) |>
  summarise(cases_muni_agg = sum(cases_confirmed, na.rm = TRUE), .groups = "drop")

state_from_stan <- tibble(week_start = as.Date(weekly_data$week_start), cases_stan = weekly_data$cases)

reconciled <- inner_join(state_from_stan, state_from_muni, by = "week_start")

n_weeks <- nrow(reconciled)
total_stan <- sum(reconciled$cases_stan)
total_muni <- sum(reconciled$cases_muni_agg)
abs_discrepancy <- abs(total_stan - total_muni)
rel_discrepancy <- abs_discrepancy / total_stan
weekly_diff <- abs(reconciled$cases_stan - reconciled$cases_muni_agg)
max_weekly_discrepancy <- max(weekly_diff)

recon_report <- tibble(
  n_weeks = n_weeks, total_state_cases_stan_input = total_stan,
  total_cases_municipality_aggregation = total_muni,
  absolute_discrepancy = abs_discrepancy, relative_discrepancy = rel_discrepancy,
  max_weekly_discrepancy = max_weekly_discrepancy,
  n_weeks_with_any_discrepancy = sum(weekly_diff > 0)
)

message("=== Section 2: municipality -> state reconstruction check ===")
print(as.data.frame(recon_report))
write_csv(recon_report, file.path(table_dir, "bahia_reconstruction_check.csv"))

if (rel_discrepancy > 0.001) {
  stop("STOP: municipality aggregation does not reproduce the Stan input state series closely enough (relative discrepancy > 0.1%). Investigate before proceeding.")
}
message("\nCONFIRMED: municipality aggregation reproduces the exact Stan-input Bahia state case series (relative discrepancy = ",
        signif(rel_discrepancy, 4), ").")

# ---- D continued: extract required posterior quantities -------------------
message("\n=== Posterior quantities available in the fit ===")
q_draws <- as.vector(rstan::extract(fit, "q")$q)
X_draws <- rstan::extract(fit, "X")$X       # [draws x N]
S_draws <- rstan::extract(fit, "S")$S       # [draws x N]
U_draws <- rstan::extract(fit, "U")$U       # [draws x N]
message(sprintf("q draws: %d | X/S/U dims: %d x %d weeks", length(q_draws), nrow(X_draws), ncol(X_draws)))

# ---- E: the six Bahia serosurveys (already validated window audit) --------
source(file.path(root, "02_Script/40_renewal_model/30_bahia_v4_9_replication/scripts/06_build_bahia_serology_windows.R"))
sero_windows <- build_bahia_sero_stan_fields(as.Date(weekly_data$week_start))
message("\n=== Six Bahia serosurveys (existing, validated) ===")
print(as.data.frame(sero_windows$audit |> select(sero_id, location, n_positive, n_tested, observed_prevalence, start_idx, n_weeks)))
message("\nNote: U19-U22 are four distinct COMMUNITIES within Salvador municipality -- their municipality-level reported-case CONTEXT (Salvador) is identical or nearly identical by construction. Municipality-level incidence cannot resolve within-Salvador community heterogeneity (Section 1 caveat).")

# ---- Population source check -----------------------------------------------
message("\n=== Population source ===")
message("Municipality population: chik_dlnm_panel_muni_week_2015_2025.rds 'population' column (weekly, per-municipality; same source as the renewal pipeline's muni-week panel).")
message("State population (for the crude arithmetic checks): weekly_data$N_start / N_end from the frozen model's IBGE population-projection-revision-2024 series (same demographic vintage used throughout this project).")

saveRDS(list(bahia_panel = bahia_panel, fit_bundle = fit_bundle, sero_windows = sero_windows,
             q_draws = q_draws, X_draws = X_draws, S_draws = S_draws, U_draws = U_draws),
        file.path(table_dir, "audit_inputs_bundle.rds"))
message("\n[saved] ", file.path(table_dir, "audit_inputs_bundle.rds"), " (shared inputs for subsequent audit scripts)")
