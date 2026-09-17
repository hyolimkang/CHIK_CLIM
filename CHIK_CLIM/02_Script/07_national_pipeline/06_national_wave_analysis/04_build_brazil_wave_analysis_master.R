# Integrate frozen V2 wave definitions with nationwide early-Re and
# retrospective susceptibility outputs.  No climate or downstream association
# model is fitted here.

required_master_packages <- c("here")
missing_master_packages <- required_master_packages[
  !vapply(required_master_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_master_packages)) stop("Missing package(s): ", paste(missing_master_packages, collapse = ", "))

master_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(master_root(), "02_Script", "07_national_pipeline", "06_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

read_required_csv <- function(path, instruction) {
  if (!file.exists(path)) stop(instruction)
  read_csv(path, show_col_types = FALSE)
}

run_national_wave_master <- function() {
  paths <- ensure_national_output_dirs()
  census <- read_national_wave_census() |> dplyr::filter(major_epidemic_primary)
  early_re <- read_required_csv(
    file.path(paths$table, "brazil_chik_major_wave_early_re.csv"),
    "Run 01_fit_brazil_major_wave_early_re.R first."
  )
  onset_s <- read_required_csv(
    file.path(paths$table, "brazil_chik_wave_onset_susceptibility.csv"),
    "Run 03_reconstruct_brazil_weekly_susceptibility.R first."
  )
  master <- census |>
    select(state, wave_id, wave_order, is_recurrence, onset_week, peak_week, end_week,
           total_cases, cumulative_incidence_per_100k) |>
    left_join(
      early_re |>
        select(state, wave_id, weeks_since_previous_wave, previous_wave_total_cases,
               preceding_trough_ratio, re6_prefit_eligible, re6_fit_status, re6_fit_pass, re6_analysis_usable,
               Re_early_4_median, Re_early_4_q025, Re_early_4_q975,
               Re_early_6_median, Re_early_6_q025, Re_early_6_q975,
               Re_early_8_median, Re_early_8_q025, Re_early_8_q975),
      by = c("state", "wave_id")
    ) |>
    left_join(onset_s, by = c("state", "wave_id")) |>
    mutate(calendar_year = year(as.Date(onset_week))) |>
    arrange(state, onset_week)
  write_csv(master, file.path(paths$table, "brazil_chik_wave_analysis_master.csv"))

  annual_hmc <- read_required_csv(file.path(paths$table, "brazil_chik_annual_shape_hmc_gate.csv"),
                                  "Annual SHAPE HMC output is missing.")
  q_implied <- read_required_csv(file.path(paths$table, "brazil_chik_annual_shape_q_implied.csv"),
                                 "Annual SHAPE q_implied output is missing.")
  annual_s <- read_required_csv(file.path(paths$table, "brazil_chik_annual_susceptibility_summary.csv"),
                                "Annual susceptibility summary is missing.")
  boundary <- read_required_csv(file.path(paths$table, "brazil_chik_weekly_susceptibility_year_boundary_check.csv"),
                                "Weekly boundary check is missing.")
  diagnostics <- bind_rows(
    tibble(component = "early_Re", metric = c("total_major_waves", "Re6_prefit_eligible", "Re6_fits_attempted", "Re6_HMC_pass", "Re6_HMC_fail", "primary_usable_Re6"),
           value = c(nrow(master), sum(master$re6_prefit_eligible, na.rm = TRUE),
                     sum(master$re6_fit_status == "fitted", na.rm = TRUE), sum(master$re6_fit_pass, na.rm = TRUE),
                     sum(master$re6_prefit_eligible & !master$re6_fit_pass, na.rm = TRUE), sum(master$re6_analysis_usable, na.rm = TRUE))),
    tibble(component = "susceptibility", metric = c("states_fitted", "states_passing_HMC", "states_failing_HMC", "annual_lambda_min_q025", "annual_lambda_max_q975", "median_end_2025_S_prop", "maximum_weekly_annual_S_discrepancy_prop"),
           value = c(sum(annual_hmc$fit_status == "fitted"), sum(annual_hmc$hmc_pass, na.rm = TRUE), sum(!annual_hmc$hmc_pass, na.rm = TRUE),
                     min(annual_s$lambda_q025, na.rm = TRUE), max(annual_s$lambda_q975, na.rm = TRUE),
                     median(as.numeric(sub("%$", "", dplyr::filter(annual_s, year == 2025)$S_end_prop_median)) / 100, na.rm = TRUE), # S_end_prop_median is stored as a "84.6%"-style string; converted to a 0-1 proportion here only so this diagnostic median is numeric -- no change to the underlying value
                     max(boundary$max_abs_discrepancy_prop, na.rm = TRUE))),
    tibble(component = "integration", metric = c("major_waves_with_S_onset", "major_waves_with_usable_Re6", "major_waves_with_both", "recurrent_waves_with_both"),
           value = c(sum(!is.na(master$S_onset_median)), sum(master$re6_analysis_usable, na.rm = TRUE),
                     sum(!is.na(master$S_onset_median) & master$re6_analysis_usable, na.rm = TRUE),
                     sum(master$is_recurrence & !is.na(master$S_onset_median) & master$re6_analysis_usable, na.rm = TRUE)))
  )
  write_csv(diagnostics, file.path(paths$table, "brazil_chik_wave_analysis_diagnostics.csv"))
  write_csv(q_implied, file.path(paths$table, "brazil_chik_state_q_implied_diagnostic.csv"))
  invisible(list(master = master, diagnostics = diagnostics))
}

if (sys.nframe() == 0L) run_national_wave_master()
