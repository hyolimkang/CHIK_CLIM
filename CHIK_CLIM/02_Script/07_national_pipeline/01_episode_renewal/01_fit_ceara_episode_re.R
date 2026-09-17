# Ceara outbreak-level early-phase renewal model.
#
# Scientific target
# -----------------
# Estimate one effective reproduction number (Re) for each pre-defined Ceara
# outbreak, using only its first 4, 6, or 8 observed weeks. This intentionally
# does not fit a ten-year latent R0(t), susceptible trajectory, reporting
# trajectory, or demographic process.
#
# Incidence is reported SINAN case incidence. Thus Re has the usual renewal
# interpretation only if reporting probability and reporting delay are roughly
# stable within an individual short fitting window.
#
# Inputs
# ------
# - 01_Data/ce_weekly_2014_2025.rds
# - 03_Output/07_national_pipeline/tables/chik_state_epidemic_wave_audit.csv
#
# Outputs
# -------
# - 03_Output/07_national_pipeline/model_fits/episode_renewal/renewal_ceara_episode_re_fit.rds
# - 03_Output/07_national_pipeline/tables/renewal_episode_re/ceara_early_phase_re/episode_windows.csv
#
# Environment overrides: RENEWAL_EPISODE_ITER, RENEWAL_EPISODE_WARMUP,
# RENEWAL_EPISODE_CHAINS. Defaults are 2,000 / 1,000 / 4.

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


# -----------------------------------------------------------------------------
# Fixed analytical choices
# -----------------------------------------------------------------------------

episode_settings <- list(
  state = "CE",
  generation_interval_weeks = 8L,
  # Discretised Gamma(shape = 4, rate = 2), identical support to v2.x models.
  generation_interval_weights = {
    weights <- diff(pgamma(0:8, shape = 4, rate = 2))
    weights / sum(weights)
  },
  sensitivity_windows = c(4L, 6L, 8L),
  primary_window = 8L,
  control = list(adapt_delta = 0.90, max_treedepth = 11)
)

read_positive_integer_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)

  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < 1L) stop(name, " must be a positive integer.")
  value
}


# -----------------------------------------------------------------------------
# Episode definition and renewal input construction
# -----------------------------------------------------------------------------

read_ceara_episode_starts <- function(settings) {
  audit_path <- here::here("03_Output", "07_national_pipeline", "tables", "chik_state_epidemic_wave_audit.csv")
  audit <- read.csv(audit_path, stringsAsFactors = FALSE)

  episodes <- audit |>
    dplyr::filter(uf == settings$state) |>
    transmute(
      episode_id = sprintf("CE_wave_%02d", wave_id),
      wave_id = as.integer(wave_id),
      start_week = as.Date(start_week),
      audit_end_week = as.Date(end_week),
      audit_peak_week = as.Date(peak_week),
      audit_total_cases = total_cases,
      audit_peak_cases = peak_cases
    ) |>
    arrange(start_week)

  if (nrow(episodes) == 0L) {
    stop("No Ceara outbreak episodes were found in the wave audit.")
  }

  list(episodes = episodes, audit_path = audit_path)
}

read_ceara_weekly_cases <- function() {
  weekly_path <- here::here("01_Data", "ce_weekly_2014_2025.rds")
  weekly <- readRDS(weekly_path) |>
    transmute(week_start = as.Date(week_start), cases = as.integer(cases)) |>
    arrange(week_start)

  if (anyNA(weekly$cases) || any(weekly$cases < 0L)) {
    stop("Weekly reported cases must be non-missing non-negative integers.")
  }
  if (!identical(weekly$week_start, seq(min(weekly$week_start), max(weekly$week_start), by = "week"))) {
    stop("Ceara weekly incidence is incomplete or not in chronological order.")
  }

  list(weekly = weekly, weekly_path = weekly_path)
}

calculate_renewal_lambda <- function(incidence, target_index, weights) {
  generation_lags <- seq_along(weights)
  history_index <- target_index - generation_lags
  if (any(history_index < 1L)) {
    stop("An episode has insufficient pre-onset incidence history.")
  }
  sum(weights * incidence[history_index])
}

build_episode_data <- function(weekly, episodes, window_weeks, settings) {
  G <- settings$generation_interval_weeks
  rows <- lapply(seq_len(nrow(episodes)), function(episode_number) {
    start_index <- match(episodes$start_week[episode_number], weekly$week_start)
    target_index <- start_index + seq_len(window_weeks) - 1L

    if (is.na(start_index) || max(target_index) > nrow(weekly)) {
      stop("Episode ", episodes$episode_id[episode_number], " is outside the weekly incidence panel.")
    }
    if (start_index <= G) {
      stop("Episode ", episodes$episode_id[episode_number], " lacks ", G, " weeks of renewal history.")
    }

    lambda <- vapply(
      target_index,
      calculate_renewal_lambda,
      numeric(1),
      incidence = weekly$cases,
      weights = settings$generation_interval_weights
    )
    if (any(lambda <= 0)) {
      stop(
        "Episode ", episodes$episode_id[episode_number],
        " has zero renewal infectiousness in its target window; this no-importation model is undefined."
      )
    }

    tibble(
      episode_id = episodes$episode_id[episode_number],
      wave_id = episodes$wave_id[episode_number],
      window_weeks = window_weeks,
      week_in_episode = seq_len(window_weeks),
      week_start = weekly$week_start[target_index],
      I_obs = weekly$cases[target_index],
      Lambda = lambda,
      audit_peak_week = episodes$audit_peak_week[episode_number],
      audit_total_cases = episodes$audit_total_cases[episode_number]
    )
  })

  bind_rows(rows)
}

make_stan_data <- function(episode_data) {
  episode_ids <- unique(episode_data$episode_id)
  window_weeks <- unique(episode_data$window_weeks)
  stopifnot(length(window_weeks) == 1L)

  list(
    K = length(episode_ids),
    W = window_weeks,
    I_obs = matrix(
      episode_data$I_obs,
      nrow = length(episode_ids),
      byrow = TRUE
    ),
    Lambda = matrix(
      episode_data$Lambda,
      nrow = length(episode_ids),
      byrow = TRUE
    )
  )
}


# -----------------------------------------------------------------------------
# Fitting workflow
# -----------------------------------------------------------------------------

fit_one_window <- function(model, episode_data, iterations, warmup, chains, settings) {
  stan_data <- make_stan_data(episode_data)
  fit <- rstan::sampling(
    object = model,
    data = stan_data,
    chains = chains,
    iter = iterations,
    warmup = warmup,
    seed = 24052018L + stan_data$W,
    control = settings$control,
    refresh = max(1L, floor(iterations / 10L))
  )

  list(fit = fit, stan_data = stan_data, episode_data = episode_data)
}

run_episode_renewal_fit <- function() {
  iterations <- read_positive_integer_env("RENEWAL_EPISODE_ITER", 2000L)
  warmup <- read_positive_integer_env("RENEWAL_EPISODE_WARMUP", 1000L)
  chains <- read_positive_integer_env("RENEWAL_EPISODE_CHAINS", 4L)
  if (warmup >= iterations) stop("Warmup must be smaller than iterations.")

  episode_input <- read_ceara_episode_starts(episode_settings)
  incidence_input <- read_ceara_weekly_cases()
  stan_source <- here::here("02_Script", "00_shared", "stan", "current", "national", "renewal_ceara_episode_re.stan")
  output_path <- here::here("03_Output", "07_national_pipeline", "model_fits", "episode_renewal", "renewal_ceara_episode_re_fit.rds")
  table_directory <- here::here(
    "03_Output", "07_national_pipeline", "tables", "renewal_episode_re", "ceara_early_phase_re"
  )
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

  model <- rstan::stan_model(stan_source)
  fits <- list()
  episode_windows <- list()

  for (window_weeks in episode_settings$sensitivity_windows) {
    episode_data <- build_episode_data(
      weekly = incidence_input$weekly,
      episodes = episode_input$episodes,
      window_weeks = window_weeks,
      settings = episode_settings
    )
    message("[fit] first ", window_weeks, " weeks; ", nrow(episode_input$episodes), " episodes")
    fits[[as.character(window_weeks)]] <- fit_one_window(
      model, episode_data, iterations, warmup, chains, episode_settings
    )
    episode_windows[[as.character(window_weeks)]] <- episode_data
  }

  write.csv(
    bind_rows(episode_windows),
    file.path(table_directory, "episode_windows.csv"),
    row.names = FALSE
  )

  saveRDS(
    list(
      fits = fits,
      episodes = episode_input$episodes,
      weekly_cases = incidence_input$weekly,
      config = list(
        model_version = "episode_renewal_re",
        formulation = "I_obs[t] ~ NB2(Re_episode * Lambda[t], phi)",
        reported_incidence_assumption = "Reporting probability and delay are approximately stable within each short window.",
        audit_source = normalizePath(episode_input$audit_path, winslash = "/", mustWork = TRUE),
        weekly_source = normalizePath(incidence_input$weekly_path, winslash = "/", mustWork = TRUE),
        generation_interval_weights = episode_settings$generation_interval_weights,
        sensitivity_windows = episode_settings$sensitivity_windows,
        primary_window = episode_settings$primary_window,
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
  run_episode_renewal_fit()
}
