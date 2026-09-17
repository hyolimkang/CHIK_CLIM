# Fit the period-level, episode-sequential susceptibility pilot for Ceara.
#
# The primary output is FILTERED susceptibility at each episode onset. For an
# episode, the filtered fit ends in the period immediately before its onset;
# it cannot use the target episode or any later observations. A single fit over
# every period is retained separately as a SMOOTHED retrospective sensitivity.

required_packages <- c("here", "rstan", "dplyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

# See the preparation script for why the inner R-project root is resolved
# explicitly rather than assuming that here::here() already points to it.
project_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM R-project root.")
}

project_path <- function(...) file.path(project_root(), ...)

read_positive_integer_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < 1L) stop(name, " must be a positive integer.")
  value
}

default_prior_settings <- list(
  # The pilot uses one state, so state offsets are deliberately absent. The
  # intercept/trend are the requested low-dimensional detection process.
  q_alpha_prior_mean = qlogis(0.02),
  q_alpha_prior_sd = 1.0,
  q_beta_prior_mean = 0.0,
  q_beta_prior_sd = 0.15,
  # Medians: epidemic attack probability approx. 14%; inter-episode approx. 1%.
  log_lambda_epi_prior_mean = log(0.15),
  log_lambda_epi_prior_sd = 1.0,
  log_lambda_inter_prior_mean = log(0.01),
  log_lambda_inter_prior_sd = 0.75,
  log_phi_prior_mean = log(20),
  log_phi_prior_sd = 1.0
)

make_stan_data <- function(periods, priors = default_prior_settings) {
  if (nrow(periods) < 1L) stop("At least one period is required.")
  if (any(periods$N_start <= 0 | periods$N_end <= 0)) stop("Population must be positive.")
  stan_vector <- function(x) array(as.numeric(x), dim = length(x))

  c(
    list(
      P = nrow(periods),
      # Preserve a one-dimensional array for the first filtered cut fit
      # (P = 1). RStan otherwise serializes length-one integer vectors as
      # scalars, which does not match Stan's array[P] declarations.
      C = array(as.integer(periods$reported_cases_period), dim = nrow(periods)),
      is_epidemic = array(as.integer(periods$period_type == "epidemic"), dim = nrow(periods)),
      year_centered = stan_vector(periods$year_centered),
      N_start = stan_vector(periods$N_start),
      N_end = stan_vector(periods$N_end),
      births = stan_vector(periods$births),
      deaths = stan_vector(periods$deaths),
      reconciliation = stan_vector(periods$reconciliation)
    ),
    priors
  )
}

fit_one_reconstruction <- function(model, periods, iterations, warmup, chains, seed, label) {
  message(sprintf("[fit] %s: %d sequential periods, %d chains, %d iterations", label, nrow(periods), chains, iterations))
  started <- Sys.time()
  fit <- rstan::sampling(
    object = model,
    data = make_stan_data(periods),
    chains = chains,
    iter = iterations,
    warmup = warmup,
    seed = seed,
    cores = if (nrow(periods) == 1L) 1L else chains,
    control = list(adapt_delta = 0.99, max_treedepth = 12),
    refresh = max(1L, floor(iterations / 10L))
  )
  list(
    fit = fit,
    periods = periods,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    label = label
  )
}

run_ceara_episode_susceptibility_fit <- function() {
  # Keep all sampling controls configurable for a short computational check.
  iterations <- read_positive_integer_env("EPISODE_SUSCEPTIBILITY_ITER", 2000L)
  warmup <- read_positive_integer_env("EPISODE_SUSCEPTIBILITY_WARMUP", 1000L)
  chains <- read_positive_integer_env("EPISODE_SUSCEPTIBILITY_CHAINS", 4L)
  if (warmup >= iterations) stop("Warmup must be smaller than total iterations.")

  source(project_path(
    "02_Script", "07_national_pipeline", "03_episode_susceptibility_reconstruction",
    "01_prepare_ceara_episode_periods.R"
  ))
  prepared <- prepare_ceara_episode_periods()
  periods <- prepared$periods

  stan_source <- project_path("02_Script", "stan", "episode_sequential_susceptibility.stan")
  output_path <- project_path("02_Script", "stan", "ceara_episode_susceptibility_pilot_v1_fit.rds")
  table_directory <- project_path("03_Output", "tables", "episode_susceptibility", "ceara_pilot_v1")
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)

  model <- rstan::stan_model(stan_source)
  smoothed <- fit_one_reconstruction(
    model, periods, iterations, warmup, chains, seed = 12012026L, label = "smoothed full timeline"
  )

  # The first period before every audit-defined epidemic contains all and only
  # observations before the episode onset. Its S_end is the target S_at_onset.
  episode_period_index <- which(periods$period_type == "epidemic")
  filtered <- lapply(seq_along(episode_period_index), function(index) {
    target_period <- episode_period_index[index]
    cutoff_period <- target_period - 1L
    target_episode <- periods[target_period, ]
    result <- fit_one_reconstruction(
      model = model,
      periods = periods[seq_len(cutoff_period), , drop = FALSE],
      iterations = iterations,
      warmup = warmup,
      chains = chains,
      seed = 12012026L + index,
      label = paste0("filtered through ", periods$period_id[cutoff_period], " for ", target_episode$episode_id)
    )
    result$target_episode <- target_episode
    result$cutoff_period_id <- periods$period_id[cutoff_period]
    result
  })
  names(filtered) <- periods$episode_id[episode_period_index]

  bundle <- list(
    model_version = "episode_sequential_susceptibility_pilot_v1",
    formulation = paste(
      "Period-level lambda; attack = 1-exp(-lambda); C ~ NB2(q * X, phi).",
      "No weekly transmission, R0, Reff, introductions, or climate covariates."
    ),
    prepared = prepared,
    priors = default_prior_settings,
    sampling = list(iterations = iterations, warmup = warmup, chains = chains),
    smoothed = smoothed,
    filtered = filtered,
    stan_source = normalizePath(stan_source, winslash = "/", mustWork = TRUE)
  )
  saveRDS(bundle, output_path)
  write.csv(periods, file.path(table_directory, "ceara_period_definition.csv"), row.names = FALSE)
  write.csv(prepared$episodes, file.path(table_directory, "ceara_episode_definition.csv"), row.names = FALSE)
  write.csv(prepared$serology_metadata, file.path(table_directory, "serology_metadata.csv"), row.names = FALSE)
  message("[save] ", output_path)
  invisible(bundle)
}


if (sys.nframe() == 0L) {
  run_ceara_episode_susceptibility_fit()
}
