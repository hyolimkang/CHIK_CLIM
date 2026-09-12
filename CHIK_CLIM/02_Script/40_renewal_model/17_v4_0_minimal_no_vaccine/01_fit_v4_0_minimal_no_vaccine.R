# Ceará v4.0 all-age, no-vaccination renewal benchmark (2015-2025).
#
# The model deliberately estimates only a transmission intercept, 11 centred
# annual deviations, two seasonal terms, and NB observation dispersion. q and
# recurrent importation are fixed scenario inputs, not fitted parameters.

required_packages <- c("here", "rstan", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
})
rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

v4_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

# When launched by Rscript from the outer Git repository, `here::here()`
# resolves there rather than at the inner R project. Resolve the project first.
setwd(v4_root())
here::i_am("CHIK_CLIM.Rproj")
source(file.path(
  v4_root(), "02_Script", "40_renewal_model", "00_shared",
  "01_build_ceara_state_weekly.R"
))

env_integer <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < 1L) stop(name, " must be a positive integer.")
  value
}

env_number <- function(name, default, lower = -Inf, upper = Inf) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  value <- suppressWarnings(as.numeric(value))
  if (!is.finite(value) || value < lower || value > upper) {
    stop(name, " must lie in [", lower, ", ", upper, "].")
  }
  value
}

generation_weights <- function() {
  weights <- diff(pgamma(0:8, shape = 4, rate = 2))
  weights / sum(weights)
}

v4_paths <- function(root = v4_root()) {
  run_tag <- Sys.getenv("RENEWAL_V4_RUN_TAG", unset = "base")
  if (!grepl("^[A-Za-z0-9_-]+$", run_tag)) stop("RENEWAL_V4_RUN_TAG may use only letters, numbers, _ and -.")
  suffix <- if (identical(run_tag, "base")) "" else paste0("_", run_tag)
  list(
    root = root,
    stan = file.path(root, "02_Script", "stan", "renewal_ceara_v4_0_minimal_no_vaccine.stan"),
    fit = file.path(root, "02_Script", "stan", paste0("renewal_ceara_v4_0_minimal_no_vaccine_fit", suffix, ".rds")),
    tables = file.path(root, "03_Output", "tables", paste0("renewal_v4_0_minimal_no_vaccine", suffix)),
    figures = file.path(root, "03_Output", "figures", paste0("renewal_v4_0_minimal_no_vaccine", suffix)),
    run_tag = run_tag
  )
}

make_v4_data <- function(weekly, q_fixed, imports_per_week, year_effect_prior_sd) {
  dates <- as.Date(weekly$week_start)
  years <- sort(unique(as.integer(format(dates, "%Y"))))
  expected_dates <- seq(min(dates), max(dates), by = "week")
  if (!identical(dates, expected_dates)) stop("Weekly data must be complete and consecutive.")
  accounting_error <- max(abs(
    weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation
  ))
  if (accounting_error > 1e-6) stop("Demographic accounting failed before Stan fitting.")
  list(
    stan_data = list(
      N = nrow(weekly), G = 8L, C = as.integer(weekly$cases), w = generation_weights(),
      N_start = as.numeric(weekly$N_start), N_end = as.numeric(weekly$N_end),
      births = as.numeric(weekly$births), deaths = as.numeric(weekly$all_cause_deaths),
      Y = length(years), year_id = match(as.integer(format(dates, "%Y")), years),
      seasonal_sin = sin(2 * pi * seq_len(nrow(weekly)) / 52.1775),
      seasonal_cos = cos(2 * pi * seq_len(nrow(weekly)) / 52.1775),
      q_fixed = q_fixed, imports_per_week = imports_per_week,
      year_effect_prior_sd = year_effect_prior_sd
    ),
    years = years,
    demographic_accounting_error = accounting_error
  )
}

compute_hmc_gate <- function(fit, pars, max_treedepth) {
  summary <- rstan::summary(fit, pars = pars)$summary
  bfmi <- rstan::get_bfmi(fit)
  result <- list(
    divergences = rstan::get_num_divergent(fit),
    max_treedepth_hits = rstan::get_num_max_treedepth(fit),
    maximum_rhat = max(summary[, "Rhat"], na.rm = TRUE),
    minimum_bulk_ess = min(summary[, "n_eff"], na.rm = TRUE),
    bfmi = bfmi,
    max_treedepth = max_treedepth
  )
  result$hmc_pass <- result$divergences == 0L && result$max_treedepth_hits == 0L &&
    result$maximum_rhat <= 1.01 && result$minimum_bulk_ess >= 100 && all(result$bfmi >= .3)
  result
}

# RStan's default unconstrained draws in [-2, 2] can imply implausibly large
# initial R0 values for this deterministic epidemic recurrence.  These values
# are only sampler initials (not priors): all chains start near the declared
# prior centre with small, deterministic between-chain offsets.
v4_initial_values <- function(stan_data) {
  function(chain_id = 1L) {
    offsets <- c(-0.06, -0.02, 0.02, 0.06)
    offset <- offsets[(as.integer(chain_id) - 1L) %% length(offsets) + 1L]
    list(
      alpha_R = log(1.2) + offset,
      z_year = rep(offset / stan_data$year_effect_prior_sd, stan_data$Y),
      beta_sin = offset / 2,
      beta_cos = -offset / 2,
      phi_obs = 20
    )
  }
}

run_v4_0_minimal_no_vaccine <- function() {
  root <- v4_root()
  paths <- v4_paths(root)
  dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)

  settings <- list(
    date_start = as.Date("2015-01-04"),
    date_end = as.Date("2025-12-21"),
    q_fixed = env_number("RENEWAL_V4_Q_FIXED", .10, lower = 1e-6, upper = 1 - 1e-6),
    imports_per_week = env_number("RENEWAL_V4_IMPORTS_PER_WEEK", 1, lower = 0),
    year_effect_prior_sd = env_number("RENEWAL_V4_YEAR_EFFECT_PRIOR_SD", .40, lower = 1e-6),
    iter = env_integer("RENEWAL_V4_ITER", 2000L),
    warmup = env_integer("RENEWAL_V4_WARMUP", 1000L),
    chains = env_integer("RENEWAL_V4_CHAINS", 4L),
    adapt_delta = .95,
    max_treedepth = env_integer("RENEWAL_V4_MAX_TREEDEPTH", 12L)
  )
  if (settings$warmup >= settings$iter) stop("Warmup must be smaller than iterations.")

  panel <- readRDS(file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |>
    filter(week_start >= settings$date_start, week_start <= settings$date_end)
  demography_path <- file.path(root, "01_Data", "ceara_weekly_demography_2015_2025.rds")
  weekly_demography <- readRDS(demography_path) |>
    filter(week_start >= settings$date_start, week_start <= settings$date_end)
  weekly <- weekly_cases |>
    select(-population) |>
    inner_join(weekly_demography, by = "week_start") |>
    arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) {
    stop("Case and demographic panels do not have the same full weekly sequence.")
  }

  prepared <- make_v4_data(
    weekly, settings$q_fixed, settings$imports_per_week, settings$year_effect_prior_sd
  )
  message(sprintf(
    "[v4.0] %s to %s | %d weeks | %s reported cases | q fixed=%.3f | imports/week=%.2f",
    settings$date_start, settings$date_end, prepared$stan_data$N,
    format(sum(prepared$stan_data$C), big.mark = ","), settings$q_fixed, settings$imports_per_week
  ))
  started <- Sys.time()
  fit <- rstan::stan(
    file = paths$stan, data = prepared$stan_data,
    chains = settings$chains, iter = settings$iter, warmup = settings$warmup,
    seed = 20260911L,
    init = v4_initial_values(prepared$stan_data),
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth),
    refresh = max(1L, floor(settings$iter / 10L))
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs"), settings$max_treedepth)
  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly,
    config = c(settings, list(
      model_version = "v4_0_minimal_no_vaccine", run_tag = paths$run_tag, years = prepared$years,
      demographic_accounting_error = prepared$demographic_accounting_error,
      demography_source = normalizePath(demography_path, winslash = "/", mustWork = TRUE),
      elapsed_seconds = elapsed_seconds, stan_source = normalizePath(paths$stan, winslash = "/", mustWork = TRUE),
      q_role = "fixed observation-scale sensitivity input; not estimated and not used to calibrate susceptibility",
      imports_role = "fixed weekly external infection pressure; not estimated",
      serology_role = "external posterior predictive consistency check only"
    )),
    hmc = hmc
  )
  saveRDS(bundle, paths$fit)
  message("[v4.0] HMC gate: ", if (hmc$hmc_pass) "PASS" else "FAIL", "; saved: ", paths$fit)
  invisible(bundle)
}

if (sys.nframe() == 0L) run_v4_0_minimal_no_vaccine()
