# Load the accepted climate-forced v4.9 fit (READ-ONLY) and truncate every
# input to the 2015-2021 primary analysis window. No refitting.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble) })

load_accepted_inputs <- function() {
  climate <- readRDS(ACCEPTED_CLIMATE_FIT)
  weekly_full <- climate$weekly_data
  dates_full <- as.Date(weekly_full$week_start)
  keep <- dates_full >= WINDOW_START & dates_full <= WINDOW_END
  if (!any(dates_full == WINDOW_END)) stop("WINDOW_END is not an exact week boundary in the fitted series -- check WINDOW_END.")

  sd_full <- climate$stan_data
  N_keep <- sum(keep)
  message(sprintf("[load] Truncating %d fitted weeks -> %d weeks (%s to %s)", length(dates_full), N_keep, WINDOW_START, WINDOW_END))

  if (any(sd_full$is_seed[keep] == 1)) stop("STOP: the 2022 seed window overlaps the requested 2015-2021 window -- check WINDOW_END.")

  age_shares <- read_csv(ACCEPTED_AGE_SHARES, show_col_types = FALSE)
  years_keep <- as.integer(format(dates_full[keep], "%Y"))
  if (!all(unique(years_keep) %in% age_shares$year)) stop("STOP: age-share coverage does not span the requested window.")

  draws <- rstan::extract(climate$fit, pars = "R0_t", permuted = TRUE)$R0_t
  n_draws_total <- nrow(draws)

  list(
    dates = dates_full[keep], years = years_keep, N = N_keep,
    R0_t_all_draws = draws[, keep, drop = FALSE], n_draws_total = n_draws_total,
    w = sd_full$w, imports_per_week = sd_full$imports_per_week,
    is_seed = integer(N_keep), X_seed = numeric(N_keep), # all-zero: no seed window within 2015-2021
    births = weekly_full$births[keep], deaths = weekly_full$all_cause_deaths[keep],
    reconciliation = weekly_full$net_population_reconciliation[keep], N_start = weekly_full$N_start[keep],
    age_shares = age_shares |> dplyr::filter(year %in% years_keep),
    source_fit_path = ACCEPTED_CLIMATE_FIT, source_age_shares_path = ACCEPTED_AGE_SHARES
  )
}

# NOTE: this script only DEFINES load_accepted_inputs(); it does not call
# it or source 01_config.R itself. Driver scripts source 01_config.R then
# this file, then call load_accepted_inputs().
