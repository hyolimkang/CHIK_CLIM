# Bahia PRIMARY serology surveys (6 observations) for the geographically-
# adjusted integrated-serology extension of the frozen v4.9 model. Source:
# Brazil_CHIKV_serology_inventory_2026-09-10.xlsx (aggregate counts only,
# no individual-level data). Excludes U23 (pooled duplicate of U19-U22),
# U01/U02 (hotspots), U09/U24 (restricted cohorts), U28 (inadequately
# retrieved), U36 (uncertain sampling frame/window), U13 (no numerator --
# reserved for a later sensitivity analysis using its reported prevalence
# + 95% CI only, NOT implemented here).
#
# Builds and AUDITS the exact model weeks each survey's collection window
# covers (uniform weighting over model weeks whose week_start falls within
# [survey_start, survey_end]) -- done here in R, not in Stan. Confirms each
# survey's covered weeks form one contiguous run (true by construction: a
# closed date interval over a strictly ordered weekly series), so Stan can
# use a simple (start_idx, n_weeks) pair and average over that contiguous
# range, exactly equivalent to a uniform weight vector summing to 1.

suppressPackageStartupMessages({ library(dplyr); library(tibble); library(readr) })

bahia_primary_serosurveys <- tibble::tribble(
  ~sero_id, ~location, ~survey_start, ~survey_end, ~endpoint, ~n_tested, ~n_positive,
  "U03", "Chapada rural district, Riachao do Jacuipe", "2016-04-01", "2016-04-30", "IgM or IgG", 120L, 24L,
  "U04", "Pau da Lima, Salvador", "2016-11-01", "2017-02-28", "IgG", 1772L, 209L,
  "U19", "Alto do Cabrito, Salvador", "2018-03-01", "2018-10-31", "IgG", 376L, 23L,
  "U20", "Marechal Rondon, Salvador", "2018-03-01", "2018-10-31", "IgG", 337L, 16L,
  "U21", "Nova Constituinte, Salvador", "2018-03-01", "2018-10-31", "IgG", 305L, 69L,
  "U22", "Rio Sena, Salvador", "2018-03-01", "2018-10-31", "IgG", 298L, 13L
) |>
  mutate(survey_start = as.Date(survey_start), survey_end = as.Date(survey_end),
         observed_prevalence = n_positive / n_tested)

build_sero_window_index <- function(dates, survey_start, survey_end) {
  idx <- which(dates >= survey_start & dates <= survey_end)
  if (length(idx) == 0L) stop("No model weeks found in survey window ", survey_start, " to ", survey_end)
  if (!all(diff(idx) == 1L)) stop("Survey window weeks are not contiguous -- unexpected for a closed date interval on a sorted weekly series")
  list(start_idx = min(idx), n_weeks = length(idx), week_starts = dates[idx])
}

build_bahia_sero_stan_fields <- function(dates, surveys = bahia_primary_serosurveys) {
  windows <- lapply(seq_len(nrow(surveys)), function(j) {
    build_sero_window_index(dates, surveys$survey_start[j], surveys$survey_end[j])
  })
  audit <- bind_rows(lapply(seq_along(windows), function(j) {
    tibble(sero_id = surveys$sero_id[j], location = surveys$location[j],
           survey_start = surveys$survey_start[j], survey_end = surveys$survey_end[j],
           n_tested = surveys$n_tested[j], n_positive = surveys$n_positive[j],
           observed_prevalence = surveys$observed_prevalence[j],
           start_idx = windows[[j]]$start_idx, n_weeks = windows[[j]]$n_weeks,
           first_week = min(windows[[j]]$week_starts), last_week = max(windows[[j]]$week_starts))
  }))
  list(
    stan_fields = list(
      J_sero = nrow(surveys),
      sero_n_positive = surveys$n_positive,
      sero_n_tested = surveys$n_tested,
      sero_window_start_idx = vapply(windows, function(w) w$start_idx, integer(1)),
      sero_window_n_weeks = vapply(windows, function(w) w$n_weeks, integer(1))
    ),
    audit = audit
  )
}

if (sys.nframe() == 0L) {
  root <- ROOT # from 00_project_setup.R (sourced by .Rprofile) -- no scientific change
  b <- readRDS(file.path(root, "03_Output/03_bahia_pipeline/model_fits/v4_9_replication/outputs/q0.05/renewal_bahia_v4_9_fit_q0.05.rds"))
  dates <- as.Date(b$weekly_data$week_start)
  built <- build_bahia_sero_stan_fields(dates)
  print(as.data.frame(built$audit))
  out_dir <- file.path(root, "03_Output/03_bahia_pipeline/tables/renewal_bahia_v4_9_replication")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(built$audit, file.path(out_dir, "bahia_serology_window_audit.csv"))
  message("[saved] ", file.path(out_dir, "bahia_serology_window_audit.csv"))
}
