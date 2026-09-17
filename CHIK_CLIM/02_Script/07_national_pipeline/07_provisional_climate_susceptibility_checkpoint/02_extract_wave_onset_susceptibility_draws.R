# Raw per-draw susceptible proportion at each major epidemic wave's onset
# week, for the provisional Re x climate x susceptibility checkpoint.
#
# This is a NEW, additive script -- it does not modify
# 03_reconstruct_brazil_weekly_susceptibility.R. It reuses that script's
# exact weekly-placement recursion, but (a) saves the RAW per-draw onset
# value instead of collapsing immediately to quantiles, and (b) supports an
# optional FOI_SCALE_FACTOR to rescale the long-term FOI assumption as a
# documented sensitivity scenario (Section 14 of the task spec) -- purely as
# post-processing of the already-fitted annual SHAPE posterior (lambda), no
# Stan refit. Because attack_prob = 1-exp(-lambda) is nonlinear and the
# annual S/U recursion compounds year to year, scaling requires re-running
# the FULL annual recursion (mirroring chik_dynamic_annual_foi_shape_v2.stan's
# transformed-parameters block exactly, in R) before re-running the weekly
# placement -- not just rescaling the already-extracted X/S_end in place.

required_extract_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "lubridate")
missing_extract_packages <- required_extract_packages[
  !vapply(required_extract_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_extract_packages)) stop("Missing package(s): ", paste(missing_extract_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble); library(lubridate) })

extract_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(extract_root(), "02_Script", "07_national_pipeline", "06_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

extract_settings <- list(
  years = 2015:2025,
  foi_scale_factor = as.numeric(Sys.getenv("FOI_SCALE_FACTOR", "1.0")),
  scale_tag = Sys.getenv("FOI_SCALE_TAG", "1_0")
)

# Mirrors chik_dynamic_annual_foi_shape_v2.stan's transformed-parameters
# recursion exactly (S_start[1]=N_start[1], U_start[1]=0; composition-based
# S/U split of deaths+reconciliation), but vectorised across posterior draws
# in plain R, with an explicit lambda scale factor. At foi_scale=1.0 this
# must reproduce the stored draws$X / draws$S_end within floating-point
# tolerance -- that reproduction IS the correctness check for this function.
recompute_annual_scaled <- function(lambda, annual, foi_scale) {
  n_draws <- nrow(lambda); Y <- ncol(lambda)
  S_start <- matrix(NA_real_, n_draws, Y); U_start <- matrix(NA_real_, n_draws, Y)
  S_end <- matrix(NA_real_, n_draws, Y); U_end <- matrix(NA_real_, n_draws, Y)
  X <- matrix(NA_real_, n_draws, Y)
  S_start[, 1] <- annual$N_start[1]; U_start[, 1] <- 0
  for (y in seq_len(Y)) {
    attack_prob <- -expm1(-(lambda[, y] * foi_scale))
    X[, y] <- S_start[, y] * attack_prob
    s_after <- S_start[, y] - X[, y]
    u_after <- U_start[, y] + X[, y]
    frac_s <- s_after / annual$N_start[y]
    S_end[, y] <- s_after + annual$births[y] - annual$deaths[y] * frac_s + annual$reconciliation[y] * frac_s
    U_end[, y] <- u_after - annual$deaths[y] * (1 - frac_s) + annual$reconciliation[y] * (1 - frac_s)
    if (y < Y) { S_start[, y + 1] <- S_end[, y]; U_start[, y + 1] <- U_end[, y] }
  }
  list(X = X, S_start = S_start, S_end = S_end, U_start = U_start, U_end = U_end)
}

# Same weekly-placement day-loop as 03_reconstruct_brazil_weekly_susceptibility.R,
# but taking (X, S_start, U_start) as arguments rather than extracting them
# directly from the stanfit, so it works identically for scaled or unscaled
# input.
place_weekly_and_capture_onset <- function(state, annual, X, S_start, U_start, state_week, settings, major_onsets) {
  weekly <- state_week |> dplyr::filter(state == !!state) |> arrange(week_start) |>
    mutate(year = year(week_start), smoothed_cases = pmax(centered_mean_3(reported_cases), 0))
  onset_rows <- list()

  for (y_index in seq_along(settings$years)) {
    annual_year <- settings$years[y_index]
    year_weeks <- weekly |> dplyr::filter(year == annual_year)
    if (!nrow(year_weeks)) next
    timing_unidentified <- sum(year_weeks$reported_cases) == 0
    weights <- if (timing_unidentified) rep(NA_real_, nrow(year_weeks)) else year_weeks$smoothed_cases / sum(year_weeks$smoothed_cases)

    start_date <- as.Date(sprintf("%d-01-01", annual_year)); end_date <- as.Date(sprintf("%d-12-31", annual_year))
    days <- seq(start_date, end_date, by = "day"); n_days <- length(days)
    birth_day <- annual$births[y_index] / n_days
    death_day <- annual$deaths[y_index] / n_days
    reconciliation_day <- annual$reconciliation[y_index] / n_days
    S <- S_start[, y_index]; U <- U_start[, y_index]
    n_draws <- length(S)
    start_by_week <- matrix(NA_real_, nrow = n_draws, ncol = nrow(year_weeks))
    infection_week <- if (timing_unidentified) matrix(NA_real_, nrow = n_draws, ncol = nrow(year_weeks)) else outer(X[, y_index], weights)
    date_index <- match(year_weeks$week_start, days)

    for (d in seq_along(days)) {
      week_hit <- which(date_index == d)
      if (length(week_hit)) {
        start_by_week[, week_hit] <- S
        if (!timing_unidentified) S <- S - rowSums(infection_week[, week_hit, drop = FALSE])
        U <- (S + U) - S
      }
      population_before_flow <- S + U
      susceptible_fraction <- S / population_before_flow
      S <- S + birth_day - death_day * susceptible_fraction + reconciliation_day * susceptible_fraction
      U <- U - death_day * (1 - susceptible_fraction) + reconciliation_day * (1 - susceptible_fraction)
    }
    if (timing_unidentified) next

    this_year_onsets <- major_onsets |> dplyr::filter(state == !!state, year(onset_week) == annual_year)
    if (!nrow(this_year_onsets)) next
    population_week_start <- vapply(seq_len(nrow(year_weeks)), function(j) {
      interpolate_uf_population(read_national_demography() |> dplyr::filter(state == !!state), year_weeks$week_start[j])
    }, numeric(1))
    for (i in seq_len(nrow(this_year_onsets))) {
      idx <- match(this_year_onsets$onset_week[i], year_weeks$week_start)
      if (is.na(idx)) next
      s_prop <- start_by_week[, idx] / population_week_start[idx]
      onset_rows[[length(onset_rows) + 1L]] <- tibble(
        state = state, wave_id = this_year_onsets$wave_id[i], draw = seq_len(n_draws), S_onset = s_prop
      )
    }
  }
  bind_rows(onset_rows)
}

run_extract_wave_onset_susceptibility_draws <- function() {
  paths <- ensure_national_output_dirs()
  fit_path <- file.path(paths$fit, "brazil_chik_annual_shape_fits.rds")
  if (!file.exists(fit_path)) stop("Run 02_fit_brazil_annual_shape_susceptibility.R first.")
  fit_bundle <- readRDS(fit_path)
  state_week <- read_national_state_week()
  census <- read_national_wave_census()
  major_onsets <- census |> dplyr::filter(major_epidemic_primary) |> select(state, wave_id, onset_week)

  states <- sort(names(fit_bundle$fits))
  message("[wave-onset-S-draws] foi_scale_factor = ", extract_settings$foi_scale_factor, ", ", length(states), " states")

  validation_rows <- list()
  all_onsets <- vector("list", length(states)); names(all_onsets) <- states
  for (state in states) {
    entry <- fit_bundle$fits[[state]]
    if (!isTRUE(entry$hmc$hmc_pass)) { message("  [skip] ", state, " (annual SHAPE HMC failed)"); next }
    message("  ", state)
    annual <- entry$input$annual
    draws <- rstan::extract(entry$fit, pars = c("lambda", "X", "S_end"), permuted = TRUE)
    recomputed <- recompute_annual_scaled(draws$lambda, annual, extract_settings$foi_scale_factor)

    if (extract_settings$foi_scale_factor == 1.0) {
      max_x_diff <- max(abs(recomputed$X - draws$X))
      max_send_diff <- max(abs(recomputed$S_end - draws$S_end))
      validation_rows[[state]] <- tibble(state = state, max_abs_X_diff = max_x_diff, max_abs_S_end_diff = max_send_diff)
      if (max_x_diff > 1e-6 * max(abs(draws$X)) || max_send_diff > 1e-6 * max(abs(draws$S_end))) {
        stop("Re-derived annual recursion does not reproduce the stored Stan draws for ", state,
             " -- max_abs_X_diff=", max_x_diff, ", max_abs_S_end_diff=", max_send_diff)
      }
    }

    all_onsets[[state]] <- place_weekly_and_capture_onset(
      state, annual, recomputed$X, recomputed$S_start, recomputed$U_start, state_week, extract_settings, major_onsets
    )
  }

  onset_draws <- bind_rows(all_onsets)
  if (!nrow(onset_draws)) stop("No wave-onset susceptibility draws were produced.")

  if (length(validation_rows)) {
    validation <- bind_rows(validation_rows)
    write_csv(validation, file.path(paths$table, "brazil_chik_wave_onset_susceptibility_draws_validation.csv"))
    message("[wave-onset-S-draws] scale=1.0 validation: max|X diff|=", max(validation$max_abs_X_diff),
            ", max|S_end diff|=", max(validation$max_abs_S_end_diff), " (should be ~0)")
  }

  rds_path <- file.path(paths$table, sprintf("brazil_chik_wave_onset_susceptibility_draws_scale%s.rds", extract_settings$scale_tag))
  saveRDS(onset_draws, rds_path, compress = "xz")

  summary_table <- onset_draws |> group_by(state, wave_id) |>
    summarise(
      S_onset_q025 = quantile(S_onset, .025), S_onset_q25 = quantile(S_onset, .25),
      S_onset_median = median(S_onset), S_onset_q75 = quantile(S_onset, .75), S_onset_q975 = quantile(S_onset, .975),
      n_draws = n(), .groups = "drop"
    )
  csv_path <- file.path(paths$table, sprintf("brazil_chik_wave_onset_susceptibility_summary_scale%s.csv", extract_settings$scale_tag))
  write_csv(summary_table, csv_path)

  message("[wave-onset-S-draws] saved: ", rds_path)
  message("[wave-onset-S-draws] saved: ", csv_path)
  message("[wave-onset-S-draws] ", n_distinct(onset_draws$wave_id), " waves x ", n_distinct(onset_draws$draw), " draws")
  invisible(list(draws = onset_draws, summary = summary_table))
}

if (sys.nframe() == 0L) run_extract_wave_onset_susceptibility_draws()
