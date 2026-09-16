# Pernambuco Spatial Model 1 -- Sections 4-11: build the 5-strata Stan data.
#
# Constructs per-stratum weekly case series, demographic bookkeeping
# (N_start/N_end/births/deaths, closed by a residual reconciliation term
# exactly as the state-level pipeline already does), population-scaled
# importation, and Recife-only U14 serology linkage. NO Stan fitting here.

required_packages <- c("dplyr", "readr", "tidyr", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(stringr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/40_renewal_model/17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R"))
source(file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication/scripts/02_build_pe_serology_window.R"))

spatial_table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_model")
turnover_table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")

STRATUM_LEVELS <- c("1_Recife", "2_Metropolitana_remainder", "3_Agreste", "4_Sertao", "5_Vale_Sao_Francisco_Araripe")
RECIFE_INDEX <- 1L

make_year_contrast_basis <- function(Y) {
  qr_input <- cbind(rep(1, Y), diag(Y)[, seq_len(Y - 1L), drop = FALSE])
  basis <- qr.Q(qr(qr_input))[, 2:Y, drop = FALSE]
  if (max(abs(colSums(basis))) > 1e-10 || max(abs(crossprod(basis) - diag(Y - 1L))) > 1e-10) {
    stop("Invalid sum-to-zero contrast basis.")
  }
  basis
}

build_pe_5strata_stan_data <- function() {
  stratum_lookup <- read_csv(file.path(spatial_table_dir, "PE_spatial_stratum_lookup.csv"), show_col_types = FALSE,
                              col_types = cols(muni6 = col_character()))
  pe_panel <- readRDS(file.path(turnover_table_dir, "pe_municipality_week_panel.rds")) |> mutate(muni6 = as.character(muni6))
  weekly_state <- readRDS(file.path(root, "03_Output/tables/pernambuco_v4_9_replication/pernambuco_weekly_input.rds"))

  panel_stratum <- pe_panel |> left_join(stratum_lookup |> select(muni6, final_model_stratum), by = "muni6")
  if (any(is.na(panel_stratum$final_model_stratum))) stop("STOP: unassigned municipality in the case panel.")

  weeks <- sort(unique(panel_stratum$week_start))
  N <- length(weeks)
  R <- length(STRATUM_LEVELS)
  if (!identical(weeks, sort(unique(weekly_state$week_start)))) stop("STOP: week index mismatch between muni panel and state series.")

  # ---- Cases: sum by stratum x week ------------------------------------
  case_wide <- panel_stratum |>
    group_by(final_model_stratum, week_start) |>
    summarise(cases = sum(cases_confirmed, na.rm = TRUE), .groups = "drop") |>
    mutate(final_model_stratum = factor(final_model_stratum, levels = STRATUM_LEVELS)) |>
    arrange(week_start, final_model_stratum) |>
    pivot_wider(names_from = final_model_stratum, values_from = cases)
  if (!identical(case_wide$week_start, weeks)) stop("STOP: case_wide week order does not match the panel week index.")
  C_mat <- as.matrix(case_wide[, STRATUM_LEVELS])
  storage.mode(C_mat) <- "integer"
  if (any(is.na(C_mat))) stop("STOP: NA in the 5-strata case matrix -- a stratum x week combination is missing.")
  case_diff <- as.numeric(rowSums(C_mat)) - as.numeric(weekly_state$cases)
  if (max(abs(case_diff)) > 0) {
    message("Max abs diff: ", max(abs(case_diff)), " at weeks: ", paste(weeks[abs(case_diff) > 0], collapse = ", "))
    stop("STOP: 5-strata case aggregation does not exactly reproduce the state series.")
  }
  message("CONFIRMED: 5-strata case matrix sums exactly to the state case series (", N, " weeks, ", R, " strata).")

  # ---- Population by stratum x week (N_start), rescaled to close exactly ----
  # to the already-validated state-level N_start series (the muni-panel
  # population and the state demography pipeline are independently built;
  # rescale stratum shares -- NOT the state total -- so accounting closes
  # exactly, preserving each stratum's relative population share).
  pop_wide <- panel_stratum |>
    group_by(final_model_stratum, week_start) |>
    summarise(population = sum(population, na.rm = TRUE), .groups = "drop") |>
    mutate(final_model_stratum = factor(final_model_stratum, levels = STRATUM_LEVELS)) |>
    arrange(week_start, final_model_stratum) |>
    pivot_wider(names_from = final_model_stratum, values_from = population)
  pop_mat_raw <- as.matrix(pop_wide[, STRATUM_LEVELS])
  pop_state_from_munis <- rowSums(pop_mat_raw)
  scale_factor <- weekly_state$N_start / pop_state_from_munis
  message(sprintf("Muni-panel-summed population vs state N_start: max relative diff = %.6f%% (rescaled to close exactly)",
                   100 * max(abs(scale_factor - 1))))
  N_start_mat <- pop_mat_raw * scale_factor # each row now sums exactly to weekly_state$N_start

  # N_end[t,r] = N_start[t+1,r] for t<N (same convention as the state
  # pipeline, confirmed N_end[t]==N_start[t+1] there); final week scaled by
  # the state's own last-week growth ratio (boundary approximation, does not
  # propagate since t=N has no t=N+1 recursion step).
  N_end_mat <- rbind(N_start_mat[-1, , drop = FALSE],
                      N_start_mat[N, , drop = FALSE] * (weekly_state$N_end[N] / weekly_state$N_start[N]))
  stopifnot(max(abs(rowSums(N_end_mat) - weekly_state$N_end)) < 1e-3)

  # ---- Births / deaths: allocate state totals by population share; ----
  #      reconciliation absorbs the stratum-specific residual (identical in
  #      spirit to the state pipeline's own net_population_reconciliation
  #      catch-all for net migration/data revision).
  pop_share <- N_start_mat / rowSums(N_start_mat)
  births_mat <- pop_share * weekly_state$births
  deaths_mat <- pop_share * weekly_state$all_cause_deaths

  # ---- Importation: population-scaled allocation of the statewide rate ----
  imports_r_mat <- pop_share * 1 # imports_per_week = 1 statewide, unchanged from the homogeneous model
  stopifnot(max(abs(rowSums(imports_r_mat) - 1)) < 1e-9)

  # ---- Year / seasonal terms (shared across strata, time-indexed only) ----
  dates <- weeks
  years <- sort(unique(as.integer(format(dates, "%Y"))))
  Y <- length(years)
  year_id <- match(as.integer(format(dates, "%Y")), years)
  seasonal_sin <- sin(2 * pi * seq_len(N) / 52.1775)
  seasonal_cos <- cos(2 * pi * seq_len(N) / 52.1775)
  seasonal_sin2 <- sin(4 * pi * seq_len(N) / 52.1775)
  seasonal_cos2 <- cos(4 * pi * seq_len(N) / 52.1775)

  # ---- Recife-only U14 serology (Section 9: direct link, no eta_geo) ----
  sero_window <- build_pe_sero_window_index(dates, pe_primary_serosurvey$survey_start[1], pe_primary_serosurvey$survey_end[1])

  stan_data <- list(
    N = N, R = R, G = 8L,
    C = C_mat, w = generation_weights(),
    N_start = N_start_mat, N_end = N_end_mat, births = births_mat, deaths = deaths_mat,
    imports_r = imports_r_mat,
    Y = Y, year_id = year_id,
    seasonal_sin = seasonal_sin, seasonal_cos = seasonal_cos,
    seasonal_sin2 = seasonal_sin2, seasonal_cos2 = seasonal_cos2,
    year_effect_prior_sd = 0.40,
    year_contrast_basis = make_year_contrast_basis(Y),
    region_contrast_basis = make_year_contrast_basis(R),
    is_seed = matrix(0L, N, R), X_seed = matrix(0, N, R),
    RECIFE_INDEX = RECIFE_INDEX,
    J_sero = 1L,
    sero_n_positive = as.array(pe_primary_serosurvey$n_positive[1]),
    sero_n_tested = as.array(pe_primary_serosurvey$n_tested[1]),
    sero_window_start_idx = as.array(sero_window$start_idx),
    sero_window_n_weeks = as.array(sero_window$n_weeks),
    kappa_sero = 50,
    logit_q_prior_mean = qlogis(0.10), logit_q_prior_sd = 1.0
  )

  list(stan_data = stan_data, weeks = weeks, stratum_levels = STRATUM_LEVELS, recife_index = RECIFE_INDEX,
       sero_window = sero_window, scale_factor_summary = summary(as.vector(scale_factor)))
}

if (sys.nframe() == 0L) {
  built <- build_pe_5strata_stan_data()
  message("\n=== Stan data dimensions ===")
  message("N=", built$stan_data$N, " R=", built$stan_data$R, " Y=", built$stan_data$Y, " J_sero=", built$stan_data$J_sero)
  message("Strata (in matrix-column order): ", paste(built$stratum_levels, collapse = ", "))
  message("Recife serology window: start_idx=", built$sero_window$start_idx, " n_weeks=", built$sero_window$n_weeks)
  saveRDS(built, file.path(spatial_table_dir, "PE_5strata_stan_data.rds"))
  message("[saved] ", file.path(spatial_table_dir, "PE_5strata_stan_data.rds"))
}
