# Retrospective weekly placement of annual SHAPE-model infections.
#
# Annual posterior X is fixed before this step.  Three-week-smoothed confirmed
# cases determine only within-year timing weights.  No weekly Stan model and no
# re-estimation of infection magnitude are used here.

required_weekly_s_packages <- c("here", "rstan")
missing_weekly_s_packages <- required_weekly_s_packages[
  !vapply(required_weekly_s_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_weekly_s_packages)) stop("Missing package(s): ", paste(missing_weekly_s_packages, collapse = ", "))

weekly_s_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(weekly_s_root(), "02_Script", "40_renewal_model", "15_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

weekly_s_settings <- list(years = 2015:2025, discrepancy_tolerance_prop = 1e-4)

matrix_quantiles <- function(x, probabilities = c(.025, .25, .5, .75, .975)) {
  result <- t(apply(x, 2, quantile, probs = probabilities, na.rm = TRUE))
  colnames(result) <- c("q025", "q25", "median", "q75", "q975")
  as_tibble(result)
}

reconstruct_one_state <- function(state, fit_bundle, state_week, demography, settings) {
  stored <- fit_bundle$fits[[state]]
  if (is.null(stored) || !isTRUE(stored$hmc$hmc_pass)) {
    return(list(status = "annual_shape_hmc_failed"))
  }
  annual <- stored$input$annual
  draws <- rstan::extract(stored$fit, pars = c("X", "S_start", "S_end", "U_start"), permuted = TRUE)
  n_draws <- nrow(draws$X)
  weekly <- state_week |> filter(state == !!state) |> arrange(week_start) |>
    mutate(year = year(week_start), smoothed_cases = pmax(centered_mean_3(reported_cases), 0))
  weekly_rows <- list(); onset_stock <- list(); discrepancy <- list(); timing_audit <- list()

  for (y_index in seq_along(settings$years)) {
    annual_year <- settings$years[y_index]
    year_weeks <- weekly |> filter(year == annual_year)
    if (!nrow(year_weeks)) stop("No weekly surveillance rows for ", state, " in ", annual_year)
    timing_unidentified <- sum(year_weeks$reported_cases) == 0
    weights <- if (timing_unidentified) rep(NA_real_, nrow(year_weeks)) else year_weeks$smoothed_cases / sum(year_weeks$smoothed_cases)
    if (!timing_unidentified && abs(sum(weights) - 1) > 1e-12) stop("Timing weights do not sum to one.")

    start_date <- as.Date(sprintf("%d-01-01", annual_year))
    end_date <- as.Date(sprintf("%d-12-31", annual_year))
    days <- seq(start_date, end_date, by = "day")
    n_days <- length(days)
    birth_day <- annual$births[y_index] / n_days
    death_day <- annual$deaths[y_index] / n_days
    reconciliation_day <- annual$reconciliation[y_index] / n_days
    S <- draws$S_start[, y_index]
    U <- draws$U_start[, y_index]
    start_by_week <- matrix(NA_real_, nrow = n_draws, ncol = nrow(year_weeks))
    infection_week <- if (timing_unidentified) matrix(NA_real_, nrow = n_draws, ncol = nrow(year_weeks)) else outer(draws$X[, y_index], weights)
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
    if (any(S < -1e-6 | U < -1e-6)) stop("Negative reconstructed stock in ", state, " ", annual_year)
    if (!timing_unidentified && max(abs(rowSums(infection_week) - draws$X[, y_index])) > 1e-6) {
      stop("Annual-to-weekly infection allocation failed in ", state, " ", annual_year)
    }
    end_discrepancy <- if (timing_unidentified) rep(NA_real_, n_draws) else S - draws$S_end[, y_index]
    discrepancy[[as.character(annual_year)]] <- tibble(
      state = state, year = annual_year,
      S_end_weekly_minus_annual_q025 = quantile(end_discrepancy, .025, na.rm = TRUE),
      S_end_weekly_minus_annual_median = median(end_discrepancy, na.rm = TRUE),
      S_end_weekly_minus_annual_q975 = quantile(end_discrepancy, .975, na.rm = TRUE),
      max_abs_discrepancy_prop = if (timing_unidentified) NA_real_ else max(abs(end_discrepancy) / annual$N_end[y_index]),
      weekly_timing_unidentified = timing_unidentified
    )
    timing_audit[[as.character(annual_year)]] <- tibble(
      state = state, year = annual_year, weekly_timing_unidentified = timing_unidentified,
      raw_cases_year = sum(year_weeks$reported_cases), smoothed_weight_sum = sum(weights, na.rm = TRUE),
      annual_infection_allocation_max_abs_error = if (timing_unidentified) NA_real_ else max(abs(rowSums(infection_week) - draws$X[, y_index]))
    )
    if (!timing_unidentified) {
      # Denominator is the contemporaneous interpolated population stock,
      # matching the annual input's 1-July population interpolation.
      population_week_start <- sapply(seq_len(nrow(year_weeks)), function(j) {
        # S+U is reconstructed implicitly; annual population interpolation is
        # used only for presenting a proportion at each weekly start.
        interpolate_uf_population(demography |> filter(state == !!state), year_weeks$week_start[j])
      })
      stock_summary <- matrix_quantiles(start_by_week / rep(population_week_start, each = n_draws))
      weekly_rows[[as.character(annual_year)]] <- bind_cols(
        year_weeks |> select(state, week_start, year, reported_cases, smoothed_cases),
        tibble(weekly_timing_unidentified = FALSE), stock_summary
      )
      onset_stock[[as.character(annual_year)]] <- list(weeks = year_weeks$week_start, stock = start_by_week,
                                                        population = population_week_start)
    } else {
      weekly_rows[[as.character(annual_year)]] <- year_weeks |>
        transmute(state, week_start, year, reported_cases, smoothed_cases,
                  weekly_timing_unidentified = TRUE, q025 = NA_real_, q25 = NA_real_, median = NA_real_, q75 = NA_real_, q975 = NA_real_)
    }
  }
  list(status = "ok", weekly = bind_rows(weekly_rows), onset_stock = onset_stock,
       discrepancy = bind_rows(discrepancy), timing_audit = bind_rows(timing_audit))
}

extract_wave_onset_s <- function(census, reconstructions) {
  major <- census |> filter(major_epidemic_primary)
  bind_rows(lapply(seq_len(nrow(major)), function(i) {
    wave <- major[i, ]
    reconstructed <- reconstructions[[wave$state]]
    if (is.null(reconstructed) || !identical(reconstructed$status, "ok")) {
      return(tibble(state = wave$state, wave_id = wave$wave_id, susceptibility_timing_status = "annual_shape_hmc_failed",
                    S_onset_q025 = NA_real_, S_onset_q25 = NA_real_, S_onset_median = NA_real_, S_onset_q75 = NA_real_, S_onset_q975 = NA_real_))
    }
    wave_year <- as.character(year(wave$onset_week))
    values <- reconstructed$onset_stock[[wave_year]]
    index <- match(wave$onset_week, values$weeks)
    if (is.na(index)) {
      return(tibble(state = wave$state, wave_id = wave$wave_id, susceptibility_timing_status = "weekly_timing_unidentified_or_date_missing",
                    S_onset_q025 = NA_real_, S_onset_q25 = NA_real_, S_onset_median = NA_real_, S_onset_q75 = NA_real_, S_onset_q975 = NA_real_))
    }
    s_prop <- values$stock[, index] / values$population[index]
    bind_cols(tibble(state = wave$state, wave_id = wave$wave_id, susceptibility_timing_status = "retrospective_smoothed"),
              summarise_draws(s_prop, "S_onset"),
              tibble(S_onset_q25 = quantile(s_prop, .25), S_onset_q75 = quantile(s_prop, .75))) |>
      select(state, wave_id, susceptibility_timing_status, S_onset_q025, S_onset_q25, S_onset_median, S_onset_q75, S_onset_q975)
  }))
}

run_national_weekly_susceptibility <- function() {
  paths <- ensure_national_output_dirs()
  fit_path <- file.path(paths$fit, "brazil_chik_annual_shape_fits.rds")
  if (!file.exists(fit_path)) stop("Run 02_fit_brazil_annual_shape_susceptibility.R first.")
  fit_bundle <- readRDS(fit_path)
  state_week <- read_national_state_week()
  demography <- read_national_demography()
  census <- read_national_wave_census()
  reconstructions <- lapply(sort(names(fit_bundle$fits)), function(state) {
    message("[weekly-S] ", state)
    reconstruct_one_state(state, fit_bundle, state_week, demography, weekly_s_settings)
  })
  names(reconstructions) <- sort(names(fit_bundle$fits))
  weekly_summary <- bind_rows(lapply(reconstructions, `[[`, "weekly"))
  timing_audit <- bind_rows(lapply(reconstructions, `[[`, "timing_audit"))
  discrepancy <- bind_rows(lapply(reconstructions, `[[`, "discrepancy"))
  onset <- extract_wave_onset_s(census, reconstructions)
  write_csv(weekly_summary, file.path(paths$table, "brazil_chik_weekly_susceptibility_summary.csv"))
  write_csv(timing_audit, file.path(paths$table, "brazil_chik_weekly_susceptibility_timing_audit.csv"))
  write_csv(discrepancy, file.path(paths$table, "brazil_chik_weekly_susceptibility_year_boundary_check.csv"))
  write_csv(onset, file.path(paths$table, "brazil_chik_wave_onset_susceptibility.csv"))
  write_csv(tibble(
    status = "retrospective_smoothed_susceptibility",
    timing_curve = "three-week centred mean of confirmed cases; no lag shift",
    surveillance_date_definition = "event_date: symptom onset when available, notification date otherwise",
    annual_anchor = "State-specific posterior annual infections from validated SHAPE model",
    initial_condition = "S_start[2015] = N_start[2015]; U_start[2015] = 0",
    demographic_flow = "IBGE annual births, deaths, and exact reconciliation; constant daily-rate allocation within calendar year; proportional S/U allocation after infections",
    boundary_rule = "Each calendar year is re-initialised from its annual posterior S_start; no discrepancy is carried forward."
  ), file.path(paths$table, "brazil_chik_weekly_susceptibility_metadata.csv"))
  invisible(list(onset = onset, discrepancy = discrepancy))
}

if (sys.nframe() == 0L) run_national_weekly_susceptibility()
