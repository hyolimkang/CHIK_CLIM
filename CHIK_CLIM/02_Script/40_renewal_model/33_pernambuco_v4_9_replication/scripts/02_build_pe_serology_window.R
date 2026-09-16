# Pernambuco external replication -- Section 4/5 primary serology (U14).
#
# U14 / source P10: Recife, household residents, stratified random
# household survey by socioeconomic level, IgM or IgG, Aug 2018-Feb 2019,
# n_tested=2070, n_positive=770 (raw prevalence 0.372, reported 95% CI
# 0.340-0.404). Recife (municipality) -- NOT Pernambuco-state-representative
# -- hence the SAME geographic-offset transport model already validated in
# Bahia is reused unchanged (no new sigma_geo, no tightened prior).
#
# U07 (Fulni-o) and U08 (Truka) are explicitly EXCLUDED from this primary
# model (Section 6) -- restricted, non-state-representative Indigenous
# cohorts, reserved for the later Model C local-context sensitivity only.

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(readr) })

pe_primary_serosurvey <- tibble::tribble(
  ~sero_id, ~source_id, ~location, ~survey_start, ~survey_end, ~endpoint, ~n_tested, ~n_positive,
  "U14", "P10", "Recife", "2018-08-01", "2019-02-28", "IgM or IgG", 2070L, 770L
) |>
  mutate(survey_start = as.Date(survey_start), survey_end = as.Date(survey_end),
         observed_prevalence = n_positive / n_tested)

build_pe_sero_window_index <- function(dates, survey_start, survey_end) {
  idx <- which(dates >= survey_start & dates <= survey_end)
  if (length(idx) == 0L) stop("No model weeks found in survey window ", survey_start, " to ", survey_end)
  if (!all(diff(idx) == 1L)) stop("Survey window weeks are not contiguous")
  list(start_idx = min(idx), n_weeks = length(idx), week_starts = dates[idx])
}

build_pe_sero_stan_fields <- function(dates, surveys = pe_primary_serosurvey) {
  windows <- lapply(seq_len(nrow(surveys)), function(j) {
    build_pe_sero_window_index(dates, surveys$survey_start[j], surveys$survey_end[j])
  })
  audit <- bind_rows(lapply(seq_along(windows), function(j) {
    tibble(sero_id = surveys$sero_id[j], location = surveys$location[j],
           survey_start = surveys$survey_start[j], survey_end = surveys$survey_end[j],
           n_tested = surveys$n_tested[j], n_positive = surveys$n_positive[j],
           observed_prevalence = surveys$observed_prevalence[j],
           start_idx = windows[[j]]$start_idx, n_weeks = windows[[j]]$n_weeks,
           first_week = min(windows[[j]]$week_starts), last_week = max(windows[[j]]$week_starts))
  }))
  # as.array() is required here (not just for J_sero>1 like Bahia): with
  # J_sero=1, a plain length-1 R vector is silently passed to rstan as a
  # scalar, and Stan's array[J_sero] declarations then fail with
  # "mismatch in number dimensions declared and found ... dims found=()".
  list(
    stan_fields = list(
      J_sero = nrow(surveys),
      sero_n_positive = as.array(surveys$n_positive),
      sero_n_tested = as.array(surveys$n_tested),
      sero_window_start_idx = as.array(vapply(windows, function(w) w$start_idx, integer(1))),
      sero_window_n_weeks = as.array(vapply(windows, function(w) w$n_weeks, integer(1)))
    ),
    audit = audit
  )
}

if (sys.nframe() == 0L) {
  root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
  weekly <- readRDS(file.path(root, "03_Output/tables/pernambuco_v4_9_replication/pernambuco_weekly_input.rds"))
  dates <- as.Date(weekly$week_start)
  built <- build_pe_sero_stan_fields(dates)
  print(as.data.frame(built$audit))
  out_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication")
  write_csv(built$audit, file.path(out_dir, "pernambuco_serology_window_audit.csv"))
  message("[saved] ", file.path(out_dir, "pernambuco_serology_window_audit.csv"))
}
