# Ceará v4.1 all-age, no-vaccination renewal benchmark (2015-2025).
#
# Same scientific backbone as the HMC-passing v4.0 (max_treedepth=14), with
# exactly two substantive changes:
#
# A. Annual transmission-deviation term: v4.0 sampled Y independent z_year
#    and recentred them (z_year - mean(z_year)), a sum-to-zero effect that
#    still retains one likelihood-redundant direction. v4.1 samples only the
#    Y-1 identified degrees of freedom via an explicit orthonormal sum-to-
#    zero basis (a normalised Helmert contrast matrix) and learns a single
#    annual-heterogeneity scale sigma_year ~ half-normal(0, 0.35).
#
# B. Reporting fraction: v4.0 fixed q_fixed = 0.10 with zero uncertainty.
#    v4.1 estimates ONE constant q with an informative prior approximating
#    the induced distribution of Beta(30,28) [p_symptomatic] *
#    Beta(20,60) [p_detect_given_symptomatic] -- Beta(16.2, 108.7), confirmed
#    adequate by Monte Carlo + KS check in 00_check_q_prior_approximation.R.
#    p_symptomatic and p_detect_given_symptomatic are NOT estimated
#    separately. The prior's (a, b) are exposed as env vars so the
#    prespecified q-prior sensitivity runs (narrower/wider, same centre) can
#    reuse this exact script without touching the primary specification.
#
# imports_per_week remains a fixed scenario input, not a fitted parameter.

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

v4_1_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

# When launched by Rscript from the outer Git repository, `here::here()`
# resolves there rather than at the inner R project. Resolve the project first.
setwd(v4_1_root())
here::i_am("CHIK_CLIM.Rproj")
source(file.path(
  v4_1_root(), "02_Script", "00_shared", "functions", "01_build_ceara_state_weekly.R"
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

# Y x (Y-1) orthonormal sum-to-zero contrast basis: normalised Helmert
# contrasts. Columns are orthogonal to each other (t(B) %*% B = I) and to the
# all-ones vector (colSums(B) = 0), verified numerically in
# 04_check_orthonormal_year_basis.R.
build_orthonormal_year_basis <- function(Y) {
  H <- stats::contr.helmert(Y)
  sweep(H, 2, sqrt(colSums(H^2)), "/")
}

v4_1_paths <- function(root = v4_1_root()) {
  run_tag <- Sys.getenv("RENEWAL_V4_1_RUN_TAG", unset = "base")
  if (!grepl("^[A-Za-z0-9_-]+$", run_tag)) stop("RENEWAL_V4_1_RUN_TAG may use only letters, numbers, _ and -.")
  suffix <- if (identical(run_tag, "base")) "" else paste0("_", run_tag)
  list(
    root = root,
    stan = file.path(root, "02_Script", "stan", "renewal_ceara_v4_1_minimal_no_vaccine.stan"),
    fit = file.path(root, "02_Script", "stan", paste0("renewal_ceara_v4_1_minimal_no_vaccine_fit", suffix, ".rds")),
    tables = file.path(root, "03_Output", "tables", paste0("renewal_v4_1_minimal_no_vaccine", suffix)),
    figures = file.path(root, "03_Output", "figures", paste0("renewal_v4_1_minimal_no_vaccine", suffix)),
    progress = file.path(root, "02_Script", "stan", "sampling_progress", paste0("v4_1", suffix)),
    run_tag = run_tag
  )
}

make_v4_1_data <- function(weekly, q_prior_a, q_prior_b, imports_per_week) {
  dates <- as.Date(weekly$week_start)
  years <- sort(unique(as.integer(format(dates, "%Y"))))
  expected_dates <- seq(min(dates), max(dates), by = "week")
  if (!identical(dates, expected_dates)) stop("Weekly data must be complete and consecutive.")
  accounting_error <- max(abs(
    weekly$N_end - weekly$N_start - weekly$births + weekly$all_cause_deaths - weekly$net_population_reconciliation
  ))
  if (accounting_error > 1e-6) stop("Demographic accounting failed before Stan fitting.")
  B_year <- build_orthonormal_year_basis(length(years))
  list(
    stan_data = list(
      N = nrow(weekly), G = 8L, C = as.integer(weekly$cases), w = generation_weights(),
      N_start = as.numeric(weekly$N_start), N_end = as.numeric(weekly$N_end),
      births = as.numeric(weekly$births), deaths = as.numeric(weekly$all_cause_deaths),
      Y = length(years), year_id = match(as.integer(format(dates, "%Y")), years),
      seasonal_sin = sin(2 * pi * seq_len(nrow(weekly)) / 52.1775),
      seasonal_cos = cos(2 * pi * seq_len(nrow(weekly)) / 52.1775),
      imports_per_week = imports_per_week,
      B_year = B_year, q_prior_a = q_prior_a, q_prior_b = q_prior_b
    ),
    years = years,
    B_year = B_year,
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

# All chains start near the declared prior centre with small, deterministic
# between-chain offsets -- these are sampler initials, not priors.
v4_1_initial_values <- function(stan_data) {
  q_prior_mean <- stan_data$q_prior_a / (stan_data$q_prior_a + stan_data$q_prior_b)
  function(chain_id = 1L) {
    offsets <- c(-0.06, -0.02, 0.02, 0.06)
    offset <- offsets[(as.integer(chain_id) - 1L) %% length(offsets) + 1L]
    list(
      alpha_R = log(1.2) + offset,
      z_year_free = rep(offset, stan_data$Y - 1L),
      sigma_year = 0.3,
      beta_sin = offset / 2,
      beta_cos = -offset / 2,
      phi_obs = 20,
      q = min(max(q_prior_mean + offset / 4, 1e-3), 1 - 1e-3)
    )
  }
}

run_v4_1_minimal_no_vaccine <- function() {
  root <- v4_1_root()
  paths <- v4_1_paths(root)
  dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(paths$progress)) unlink(paths$progress, recursive = TRUE)
  dir.create(paths$progress, recursive = TRUE, showWarnings = FALSE)

  settings <- list(
    date_start = as.Date("2015-01-04"),
    date_end = as.Date("2025-12-21"),
    # Primary prior: Beta(16.2, 108.7), confirmed by Monte Carlo (n=2e6) + KS
    # check (KS=0.0096) to adequately approximate the induced distribution of
    # Beta(30,28) * Beta(20,60) -- see 00_check_q_prior_approximation.R.
    # Sensitivity run_tags (q_prior_narrow / q_prior_wide) override these via
    # env vars while keeping approximately the same prior mean (~0.13).
    q_prior_a = env_number("RENEWAL_V4_1_Q_PRIOR_A", 16.2, lower = 1e-6),
    q_prior_b = env_number("RENEWAL_V4_1_Q_PRIOR_B", 108.7, lower = 1e-6),
    imports_per_week = env_number("RENEWAL_V4_1_IMPORTS_PER_WEEK", 1, lower = 0),
    iter = env_integer("RENEWAL_V4_1_ITER", 2000L),
    warmup = env_integer("RENEWAL_V4_1_WARMUP", 1000L),
    chains = env_integer("RENEWAL_V4_1_CHAINS", 4L),
    adapt_delta = .95,
    max_treedepth = env_integer("RENEWAL_V4_1_MAX_TREEDEPTH", 14L)
  )
  if (settings$warmup >= settings$iter) stop("Warmup must be smaller than iterations.")
  if ((settings$iter - settings$warmup) < 1000L) stop("Post-warmup draws per chain must be >= 1000.")

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

  prepared <- make_v4_1_data(weekly, settings$q_prior_a, settings$q_prior_b, settings$imports_per_week)
  q_prior_mean <- settings$q_prior_a / (settings$q_prior_a + settings$q_prior_b)
  message(sprintf(
    "[v4.1] run_tag=%s | %s to %s | %d weeks | %s reported cases | q ~ Beta(%.3g, %.3g) [mean=%.3f] | imports/week=%.2f | Y=%d",
    paths$run_tag, settings$date_start, settings$date_end, prepared$stan_data$N,
    format(sum(prepared$stan_data$C), big.mark = ","), settings$q_prior_a, settings$q_prior_b, q_prior_mean,
    settings$imports_per_week, prepared$stan_data$Y
  ))
  # A single base sample_file path lets rstan auto-suffix "_<chain_id>.csv"
  # per chain (rstan's stan() rejects a pre-built vector of per-chain paths
  # here -- passing one triggers "the condition has length > 1" inside each
  # worker, confirmed by smoke test). These per-chain CSVs let
  # 07_check_sampling_progress.R report live iteration counts and treedepth/
  # divergence by counting/parsing written rows -- rstan's Windows PSOCK
  # chain parallelisation does not stream refresh() console messages back to
  # the parent process, so this file-based row count is the only real-time
  # progress signal available while chains are running.
  sample_file <- file.path(paths$progress, "chain")
  message("[v4.1] per-chain progress CSVs: ", paths$progress, "/chain_<id>.csv")
  started <- Sys.time()
  fit <- rstan::stan(
    file = paths$stan, data = prepared$stan_data,
    chains = settings$chains, iter = settings$iter, warmup = settings$warmup,
    seed = 20260911L,
    init = v4_1_initial_values(prepared$stan_data),
    control = list(adapt_delta = settings$adapt_delta, max_treedepth = settings$max_treedepth),
    refresh = max(1L, floor(settings$iter / 10L)),
    sample_file = sample_file
  )
  elapsed_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  hmc <- compute_hmc_gate(fit, c("alpha_R", "z_year_free", "sigma_year", "beta_sin", "beta_cos", "phi_obs", "q"), settings$max_treedepth)
  bundle <- list(
    fit = fit, stan_data = prepared$stan_data, weekly_data = weekly, B_year = prepared$B_year,
    config = c(settings, list(
      model_version = "v4_1_minimal_no_vaccine", run_tag = paths$run_tag, years = prepared$years,
      demographic_accounting_error = prepared$demographic_accounting_error,
      demography_source = normalizePath(demography_path, winslash = "/", mustWork = TRUE),
      elapsed_seconds = elapsed_seconds, stan_source = normalizePath(paths$stan, winslash = "/", mustWork = TRUE),
      q_role = "estimated constant reporting fraction; informative prior from external epidemiological information, not a weak/uniform prior; not weekly, not yearly, not time-varying",
      imports_role = "fixed weekly external infection pressure; not estimated",
      serology_role = "external posterior predictive consistency check only",
      annual_effect_parameterisation = "explicit (Y-1)-dimensional orthonormal sum-to-zero basis (normalised Helmert contrasts); sigma_year ~ half-normal(0, 0.35) learned, not fixed"
    )),
    hmc = hmc
  )
  saveRDS(bundle, paths$fit)
  message("[v4.1] HMC gate: ", if (hmc$hmc_pass) "PASS" else "FAIL", "; saved: ", paths$fit)
  invisible(bundle)
}

if (sys.nframe() == 0L) run_v4_1_minimal_no_vaccine()
