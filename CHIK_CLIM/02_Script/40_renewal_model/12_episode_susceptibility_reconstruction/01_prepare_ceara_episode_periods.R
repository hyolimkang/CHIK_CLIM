# Prepare the Ceará period-level input for episode-sequential susceptibility.
#
# The surveillance series is deliberately aligned to the existing epidemic-wave
# audit: all notified SINAN chikungunya records, assigned to onset week when
# available and notification week otherwise. Confirmed-only counts are not used
# here because they would no longer match the audit-defined episode boundaries.

required_packages <- c("here", "dplyr", "tibble", "readr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
})

# This repository has an outer Git root and an inner R-project root.  `here()`
# can resolve to the outer root when this script is executed from a shell, so
# resolve the R-project root explicitly and use it for every data/output path.
project_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM R-project root.")
}

project_path <- function(...) file.path(project_root(), ...)

pilot_settings <- list(
  state_code = "23",
  state_abbreviation = "CE",
  analysis_start = as.Date("2015-01-04"),
  # The audited episode table currently ends in 2024; retain its coverage.
  analysis_end = as.Date("2024-12-29")
)

read_ceara_notified_weekly <- function(settings) {
  case_path <- project_path("01_Data", "chik_brazil_muni_week_2015_2025.csv")
  weekly <- read.csv(case_path, stringsAsFactors = FALSE) |>
    mutate(week_start = as.Date(week_start)) |>
    filter(substr(muni6, 1, 2) == settings$state_code) |>
    group_by(week_start) |>
    summarise(cases = sum(cases_notified), .groups = "drop") |>
    filter(week_start >= settings$analysis_start, week_start <= settings$analysis_end) |>
    arrange(week_start)

  expected_dates <- seq(settings$analysis_start, settings$analysis_end, by = "week")
  if (!identical(weekly$week_start, expected_dates)) {
    stop("Notified case series does not provide complete weekly coverage for the pilot period.")
  }
  list(weekly = weekly, case_path = case_path)
}

read_ceara_audit_episodes <- function(settings) {
  audit_path <- project_path("03_Output", "tables", "chik_state_epidemic_wave_audit.csv")
  episodes <- read.csv(audit_path, stringsAsFactors = FALSE) |>
    filter(uf == settings$state_abbreviation) |>
    transmute(
      state = "Ceara",
      episode_id = sprintf("CE_wave_%02d", wave_id),
      episode_start = as.Date(start_week),
      episode_end = as.Date(end_week),
      reported_cases_episode_audit = as.integer(total_cases),
      duration_weeks_audit = as.integer(duration_weeks)
    ) |>
    arrange(episode_start)

  if (nrow(episodes) == 0L) stop("No Ceará episodes found in the audit table.")
  if (any(episodes$episode_start > episodes$episode_end)) stop("Episode end precedes episode start.")
  if (any(episodes$episode_start[-1] <= episodes$episode_end[-nrow(episodes)])) {
    stop("Audit episodes overlap; period timeline cannot be partitioned.")
  }
  list(episodes = episodes, audit_path = audit_path)
}

build_sequential_periods <- function(episodes, settings) {
  periods <- list()
  cursor <- settings$analysis_start
  period_number <- 1L

  add_period <- function(type, start, end, episode_id = NA_character_) {
    if (start > end) return(invisible(NULL))
    periods[[length(periods) + 1L]] <<- tibble(
      state = "Ceara",
      period_id = sprintf("CE_period_%02d", period_number),
      period_type = type,
      episode_id = episode_id,
      period_start = start,
      period_end = end
    )
    period_number <<- period_number + 1L
  }

  for (episode_number in seq_len(nrow(episodes))) {
    episode <- episodes[episode_number, ]
    add_period("inter_episode", cursor, episode$episode_start - 7)
    add_period("epidemic", episode$episode_start, episode$episode_end, episode$episode_id)
    cursor <- episode$episode_end + 7
  }
  add_period("inter_episode", cursor, settings$analysis_end)

  bind_rows(periods)
}

aggregate_period_inputs <- function(periods, weekly_cases, settings) {
  demography_path <- project_path("01_Data", "ceara_weekly_demography_2015_2025.rds")
  demography <- readRDS(demography_path) |>
    mutate(week_start = as.Date(week_start)) |>
    filter(week_start >= settings$analysis_start, week_start <= settings$analysis_end)

  period_rows <- lapply(seq_len(nrow(periods)), function(period_number) {
    period <- periods[period_number, ]
    case_slice <- weekly_cases |>
      filter(week_start >= period$period_start, week_start <= period$period_end)
    demography_slice <- demography |>
      filter(week_start >= period$period_start, week_start <= period$period_end)

    expected_dates <- seq(period$period_start, period$period_end, by = "week")
    if (nrow(case_slice) != length(expected_dates) || nrow(demography_slice) != length(expected_dates)) {
      stop("A period does not have complete cases and demography coverage: ", period$period_id)
    }

    accounting_error <-
      demography_slice$N_end[nrow(demography_slice)] - demography_slice$N_start[1] -
      sum(demography_slice$births) + sum(demography_slice$all_cause_deaths) -
      sum(demography_slice$net_population_reconciliation)
    if (abs(accounting_error) > 1e-5) {
      stop("Period demographic accounting failed: ", period$period_id)
    }

    midpoint <- period$period_start + floor((period$period_end - period$period_start) / 2)
    bind_cols(
      period,
      tibble(
        n_weeks = length(expected_dates),
        reported_cases_period = as.integer(sum(case_slice$cases)),
        N_start = unname(demography_slice$N_start[1]),
        N_end = unname(demography_slice$N_end[nrow(demography_slice)]),
        births = sum(demography_slice$births),
        deaths = sum(demography_slice$all_cause_deaths),
        reconciliation = sum(demography_slice$net_population_reconciliation),
        year_centered = as.numeric(midpoint - as.Date("2020-07-01")) / 365.25,
        demographic_accounting_error = accounting_error
      )
    )
  })

  list(
    periods = bind_rows(period_rows),
    demography_path = demography_path
  )
}

make_serology_metadata <- function() {
  # Both surveys are municipality-scale; they are not state-representative and
  # are intentionally excluded from the primary Ceará state likelihood.
  tibble(
    study_id = c("juazeiro_2018", "quixada_2018_2019"),
    state = "Ceara",
    municipality = c("Juazeiro do Norte", "Quixada"),
    survey_start = as.Date(c("2018-06-03", "2018-06-03")),
    survey_end = as.Date(c("2018-12-30", "2019-12-29")),
    n_tested = c(404L, 409L),
    n_positive = c(103L, 289L),
    geographic_scale = "municipality",
    source = "Survey specifications supplied for renewal v2.1; publication provenance not stored in this repository.",
    primary_state_likelihood = FALSE
  )
}

prepare_ceara_episode_periods <- function() {
  case_input <- read_ceara_notified_weekly(pilot_settings)
  episode_input <- read_ceara_audit_episodes(pilot_settings)
  periods <- build_sequential_periods(episode_input$episodes, pilot_settings)
  prepared <- aggregate_period_inputs(periods, case_input$weekly, pilot_settings)
  serology_metadata <- make_serology_metadata()

  # unlist() strips the Date class. Reconstruct it explicitly before comparing
  # with the expected weekly sequence.
  all_period_weeks <- as.Date(unlist(Map(
    seq,
    prepared$periods$period_start,
    prepared$periods$period_end,
    MoreArgs = list(by = "week")
  )), origin = "1970-01-01")
  expected_weeks <- seq(pilot_settings$analysis_start, pilot_settings$analysis_end, by = "week")
  if (!identical(all_period_weeks, expected_weeks)) {
    stop("Periods do not form one complete, non-overlapping weekly timeline.")
  }
  if (sum(prepared$periods$reported_cases_period) != sum(case_input$weekly$cases)) {
    stop("Period reported cases do not sum to the original state total.")
  }

  list(
    periods = prepared$periods,
    episodes = episode_input$episodes,
    weekly_cases = case_input$weekly,
    serology_metadata = serology_metadata,
    config = list(
      state = "Ceara",
      case_definition = "Notified SINAN chikungunya records; residence municipality preferred; symptom-onset week preferred and notification week used when onset date is missing.",
      audit_source = normalizePath(episode_input$audit_path, winslash = "/", mustWork = TRUE),
      case_source = normalizePath(case_input$case_path, winslash = "/", mustWork = TRUE),
      demography_source = normalizePath(prepared$demography_path, winslash = "/", mustWork = TRUE),
      analysis_start = pilot_settings$analysis_start,
      analysis_end = pilot_settings$analysis_end
    )
  )
}


if (sys.nframe() == 0L) {
  output_directory <- project_path(
    "03_Output", "tables", "episode_susceptibility", "ceara_pilot_v1"
  )
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  prepared <- prepare_ceara_episode_periods()
  saveRDS(
    prepared,
    project_path("02_Script", "stan", "ceara_episode_susceptibility_periods.rds")
  )
  write.csv(prepared$periods, file.path(output_directory, "ceara_period_definition.csv"), row.names = FALSE)
  write.csv(prepared$episodes, file.path(output_directory, "ceara_episode_definition.csv"), row.names = FALSE)
  write.csv(prepared$serology_metadata, file.path(output_directory, "serology_metadata.csv"), row.names = FALSE)
}
