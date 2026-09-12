# Ceara major-outbreak early-phase renewal model (fixed six-week window).
#
# Purpose
# -------
# This analysis deliberately separates outbreak definition from transmission
# estimation. Candidate onsets come from the existing state epidemic-wave audit;
# only candidate episodes with >= 4,000 reported cases are retained as major
# outbreaks. Each retained outbreak has one constant Re estimated from exactly
# the first six weeks after its audit-defined onset.
#
# This is a reported-incidence renewal model. Re is interpretable as an
# effective reproduction number when reporting probability and reporting delay
# are approximately stable inside a six-week onset window.

required_packages <- c("here", "rstan", "dplyr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
  library(tibble)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

settings <- list(
  state = "CE",
  major_total_cases_min = 4000L,
  fit_window_weeks = 6L,
  generation_interval_weeks = 8L,
  generation_interval_weights = {
    weights <- diff(pgamma(0:8, shape = 4, rate = 2))
    weights / sum(weights)
  },
  control = list(adapt_delta = 0.90, max_treedepth = 11)
)

read_positive_integer_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)

  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < 1L) stop(name, " must be a positive integer.")
  value
}

read_major_episode_definition <- function(settings) {
  audit_path <- here::here("03_Output", "tables", "chik_state_epidemic_wave_audit.csv")
  audit <- read.csv(audit_path, stringsAsFactors = FALSE)

  candidates <- audit |>
    filter(uf == settings$state) |>
    transmute(
      episode_id = sprintf("CE_major_wave_%02d", wave_id),
      wave_id = as.integer(wave_id),
      onset = as.Date(start_week),
      peak_week = as.Date(peak_week),
      peak_cases = as.integer(peak_cases),
      total_episode_cases = as.integer(total_cases),
      duration_weeks = as.integer(duration_weeks),
      passes_major_threshold = total_episode_cases >= settings$major_total_cases_min
    ) |>
    arrange(onset)

  major_episodes <- candidates |>
    filter(passes_major_threshold)
  if (nrow(major_episodes) == 0L) {
    stop("No candidate episodes satisfy the major-outbreak threshold.")
  }

  list(candidates = candidates, major_episodes = major_episodes, audit_path = audit_path)
}

read_weekly_cases <- function() {
  weekly_path <- here::here("01_Data", "ce_weekly_2014_2025.rds")
  weekly <- readRDS(weekly_path) |>
    transmute(week_start = as.Date(week_start), cases = as.integer(cases)) |>
    arrange(week_start)

  expected_dates <- seq(min(weekly$week_start), max(weekly$week_start), by = "week")
  if (!identical(weekly$week_start, expected_dates) || anyNA(weekly$cases) || any(weekly$cases < 0L)) {
    stop("Weekly case input must be a complete non-negative integer time series.")
  }

  list(weekly = weekly, weekly_path = weekly_path)
}

renewal_lambda <- function(cases, target_index, weights) {
  history_index <- target_index - seq_along(weights)
  if (any(history_index < 1L)) stop("An onset has insufficient renewal history.")
  sum(weights * cases[history_index])
}

build_major_episode_data <- function(weekly, episodes, settings) {
  W <- settings$fit_window_weeks
  G <- settings$generation_interval_weeks

  bind_rows(lapply(seq_len(nrow(episodes)), function(episode_number) {
    onset_index <- match(episodes$onset[episode_number], weekly$week_start)
    target_index <- onset_index + seq_len(W) - 1L
    first_eight_index <- onset_index + seq_len(8L) - 1L

    if (is.na(onset_index) || onset_index <= G || max(first_eight_index) > nrow(weekly)) {
      stop("Episode ", episodes$episode_id[episode_number], " is outside the available incidence history.")
    }

    lambda <- vapply(
      target_index,
      renewal_lambda,
      numeric(1),
      cases = weekly$cases,
      weights = settings$generation_interval_weights
    )
    if (any(lambda <= 0)) {
      stop("Episode ", episodes$episode_id[episode_number], " has zero renewal infectiousness.")
    }

    tibble(
      episode_id = episodes$episode_id[episode_number],
      wave_id = episodes$wave_id[episode_number],
      onset = episodes$onset[episode_number],
      peak_week = episodes$peak_week[episode_number],
      peak_cases = episodes$peak_cases[episode_number],
      total_episode_cases = episodes$total_episode_cases[episode_number],
      duration_weeks = episodes$duration_weeks[episode_number],
      cases_first_4wk = sum(weekly$cases[onset_index + 0:3]),
      cases_first_8wk = sum(weekly$cases[first_eight_index]),
      week_in_episode = seq_len(W),
      week_start = weekly$week_start[target_index],
      I_obs = weekly$cases[target_index],
      Lambda = lambda
    )
  }))
}

make_stan_data <- function(episode_data) {
  episode_ids <- unique(episode_data$episode_id)
  list(
    K = length(episode_ids),
    W = settings$fit_window_weeks,
    I_obs = matrix(episode_data$I_obs, nrow = length(episode_ids), byrow = TRUE),
    Lambda = matrix(episode_data$Lambda, nrow = length(episode_ids), byrow = TRUE)
  )
}

run_major_episode_fit <- function() {
  iterations <- read_positive_integer_env("RENEWAL_MAJOR_EPISODE_ITER", 2000L)
  warmup <- read_positive_integer_env("RENEWAL_MAJOR_EPISODE_WARMUP", 1000L)
  chains <- read_positive_integer_env("RENEWAL_MAJOR_EPISODE_CHAINS", 4L)
  if (warmup >= iterations) stop("Warmup must be smaller than iterations.")

  episode_definition <- read_major_episode_definition(settings)
  incidence_input <- read_weekly_cases()
  episode_data <- build_major_episode_data(
    incidence_input$weekly, episode_definition$major_episodes, settings
  )
  stan_data <- make_stan_data(episode_data)

  stan_source <- here::here("02_Script", "stan", "renewal_ceara_episode_re.stan")
  output_path <- here::here("02_Script", "stan", "renewal_ceara_major_episode_re6_fit.rds")
  table_directory <- here::here(
    "03_Output", "tables", "renewal_major_episode_re", "ceara_major_outbreaks_6wk"
  )
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

  message(sprintf(
    "[fit] %d major episodes | onset threshold from audit | total cases >= %s | fixed %d-week window",
    stan_data$K,
    format(settings$major_total_cases_min, big.mark = ","),
    settings$fit_window_weeks
  ))
  model <- rstan::stan_model(stan_source)
  fit <- rstan::sampling(
    model,
    data = stan_data,
    chains = chains,
    iter = iterations,
    warmup = warmup,
    seed = 24052024L,
    control = settings$control,
    refresh = max(1L, floor(iterations / 10L))
  )

  write.csv(
    episode_definition$candidates,
    file.path(table_directory, "candidate_episode_audit.csv"),
    row.names = FALSE
  )
  write.csv(
    episode_data,
    file.path(table_directory, "major_episode_six_week_windows.csv"),
    row.names = FALSE
  )
  saveRDS(
    list(
      fit = fit,
      stan_data = stan_data,
      episodes = episode_definition$major_episodes,
      episode_data = episode_data,
      weekly_cases = incidence_input$weekly,
      config = list(
        model_version = "major_episode_re6",
        formulation = "I_obs[t] ~ NB2(Re_episode * Lambda[t], phi)",
        onset_rule = "Existing audit: 3-week smoothed incidence >= 1 per 100,000/week; bridge dips <= 3 weeks; candidate duration >= 4 weeks and total cases >= 50.",
        major_rule = "Candidate episode total reported cases >= 4,000.",
        major_total_cases_min = settings$major_total_cases_min,
        fit_window_weeks = settings$fit_window_weeks,
        audit_source = normalizePath(episode_definition$audit_path, winslash = "/", mustWork = TRUE),
        weekly_source = normalizePath(incidence_input$weekly_path, winslash = "/", mustWork = TRUE),
        generation_interval_weights = settings$generation_interval_weights,
        iterations = iterations,
        warmup = warmup,
        chains = chains,
        stan_source = normalizePath(stan_source, winslash = "/", mustWork = TRUE)
      )
    ),
    output_path
  )
  message("[save] ", output_path)
}


if (sys.nframe() == 0L) {
  run_major_episode_fit()
}
