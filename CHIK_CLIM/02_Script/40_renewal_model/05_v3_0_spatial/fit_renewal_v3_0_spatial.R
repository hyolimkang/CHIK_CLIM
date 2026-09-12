# =============================================================================
# fit_renewal_v3_0_spatial.R
#
# Compile and, only when explicitly requested, run HMC for the v3.0
# municipality-level renewal model. The default action is data validation
# plus Stan compilation; it never launches sampling unless one of the two
# selection switches below is set.
#
# Two mutually exclusive selection modes:
#   RENEWAL_V3_0_RUN_PILOT=true   ad-hoc debug subset (required serology/PPC
#                                  municipalities + top-5 + bottom-5 by cases)
#   RENEWAL_V3_0_MIN_CASES=250    every municipality with 2015-2019 total
#                                  cases above this threshold, plus the
#                                  required serology/PPC municipalities
#                                  regardless of their own case count. Case
#                                  counts are extremely concentrated in
#                                  Ceara (56/184 municipalities already
#                                  cover 94% of all reported cases), so this
#                                  is the intended route to a fit that
#                                  shares transmission/reporting structure
#                                  but estimates a real trajectory only
#                                  where the data can support one.
#
# Environment switches:
#   RENEWAL_V3_0_PILOT_ITER=100       total iterations (minimum 20)
#   RENEWAL_V3_0_PILOT_WARMUP=50      warmup iterations
#   RENEWAL_V3_0_PILOT_CHAINS=1       number of chains
#   RENEWAL_V3_0_ADAPT_DELTA=0.90     NUTS target acceptance
#   RENEWAL_V3_0_MAX_TREEDEPTH=12     NUTS max treedepth
# =============================================================================

required_packages <- c("here", "rstan")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))

source(here::here(
  "02_Script/40_renewal_model/05_v3_0_spatial/prepare_renewal_v3_0_spatial_data.R"
))

as_flag <- function(value) {
  tolower(value) %in% c("1", "true", "yes")
}

v3_0_parameter_count <- function(stan_data) {
  # Scalars: sigma_R, sigma_muni_R, p_symp, rho_sym_global, phi_obs.
  5L + (stan_data$N - stan_data$seed_weeks) + stan_data$M +
    stan_data$M * stan_data$seed_weeks
}

pilot_diagnostics <- function(fit) {
  list(
    divergences = rstan::get_num_divergent(fit),
    max_treedepth_hits = rstan::get_num_max_treedepth(fit),
    min_bfmi = min(rstan::get_bfmi(fit))
  )
}

write_v3_0_implementation_report <- function(prepared, compiled, pilot = NULL, sampling_chains = 1L) {
  report_path <- here::here(
    "03_Output/tables/renewal_v3_0/renewal_v3_0_spatial_implementation_report.md"
  )
  dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

  d <- prepared$diagnostics
  s <- prepared$serology
  p <- prepared$stan_data
  selection_label <- if (isTRUE(prepared$config$debug_subset)) {
    "ad-hoc debug subset"
  } else if (!is.null(prepared$config$min_total_cases)) {
    sprintf("case-count threshold subset (>%g total cases)", prepared$config$min_total_cases)
  } else {
    "all Ceara municipalities"
  }
  pilot_lines <- if (is.null(pilot)) {
    c(
      "- Pilot HMC: not run. The default workflow deliberately stops after ",
      "  deterministic validation and compilation."
    )
  } else {
    c(
      sprintf(
        "- HMC run on the %s (%d municipalities, %d chain(s)): divergences %d; maximum-treedepth hits %d; minimum E-BFMI %.3f.",
        selection_label, p$M, sampling_chains, pilot$divergences, pilot$max_treedepth_hits, pilot$min_bfmi
      ),
      if (sampling_chains == 1L) {
        "- One short chain is a computational smoke test and has no convergence claim (Rhat is undefined with a single chain)."
      } else {
        "- Multi-chain run; consult the diagnostics table for per-parameter Rhat/ESS before treating this as converged."
      }
    )
  }

  lines <- c(
    "# v3.0 spatial renewal implementation report",
    "",
    "## Files created",
    "",
    "- 02_Script/stan/renewal_ceara_v3_0_spatial.stan",
    "- 02_Script/40_renewal_model/05_v3_0_spatial/prepare_renewal_v3_0_spatial_data.R",
    "- 02_Script/40_renewal_model/05_v3_0_spatial/fit_renewal_v3_0_spatial.R",
    "- 02_Script/40_renewal_model/05_v3_0_spatial/diagnose_renewal_v3_0_spatial.R",
    "",
    "## Differences from v2.2",
    "",
    "- Replaces the single state S/U/X trajectory with one trajectory per municipality.",
    "- Applies the v2.2 demographic recursion independently to every municipality.",
    "- Replaces state-to-site serology offsets with direct municipality-window immune proportions.",
    "- Uses one shared weekly log-R0 random walk plus non-centred municipality R0 offsets.",
    "- Uses one common reporting level (p_symp * rho_sym_global); v2.2's Cear<U+00E1> reporting offset and trend are omitted.",
    "",
    "## Data and deterministic checks",
    "",
    sprintf(
      "- Fitting period: %s to %s; %d complete weeks.",
      prepared$config$date_start, prepared$config$date_end, d$weeks
    ),
    sprintf(
      "- Cear<U+00E1> panel: %d municipality-week rows, %d municipalities. Stan data: %d municipalities.",
      d$panel_rows, d$municipalities_all_ceara, d$municipalities_in_stan_data
    ),
    sprintf(
      "- Maximum demographic identity error |S + U - N_start| in deterministic recursion: %.3g.",
      d$maximum_accounting_error
    ),
    sprintf(
      "- Maximum absolute difference between summed municipality population and v2.2 state stock: %.3f%%.",
      100 * d$maximum_population_relative_gap
    ),
    sprintf("- Births: %s.", d$births_source),
    sprintf("- Deaths: %s.", d$deaths_source),
    "",
    "## Serology mapping",
    "",
    sprintf(
      "- %s (IBGE municipality code %s): %d/%d, %s to %s.",
      s$site[1], s$muni6[1], s$positive[1], s$n[1], s$window_start[1], s$window_end[1]
    ),
    sprintf(
      "- %s (IBGE municipality code %s): %d/%d, %s to %s.",
      s$site[2], s$muni6[2], s$positive[2], s$n[2], s$window_start[2], s$window_end[2]
    ),
    "- Fortaleza is retained only for an external posterior predictive comparison and is not in the likelihood.",
    "",
    "## Model size and compilation",
    "",
    sprintf(
      "- Estimated parameter dimension: %d (%d municipality seed hazards, %d shared-week R0 states, %d municipality R0 standard-normal effects, and 5 scalar parameters).",
      v3_0_parameter_count(p), p$M * p$seed_weeks,
      p$N - p$seed_weeks, p$M
    ),
    sprintf("- Stan compilation: %s.", if (compiled) "successful" else "not completed"),
    "",
    "## Pilot HMC",
    "",
    pilot_lines,
    "",
    "## Remaining warnings",
    "",
    "- Municipal all-cause deaths are a transparent population-share allocation of the documented Cear<U+00E1> annual total, not municipality-specific mortality observations.",
    "- Municipality annual population stocks are linearly interpolated between 1-January annual anchors; their summed stock is checked against the v2.2 state demographic series.",
    "- Reporting, initial infections, and susceptibility remain potentially confounded. Compilation or a short pilot is not scientific validation.",
    "- A future all-Cear<U+00E1> 2015-2019 production fit must be assessed with multi-chain convergence, posterior predictive, serology, and sensitivity diagnostics before interpretation."
  )
  writeLines(lines, report_path, useBytes = TRUE)
  report_path
}

run_v3_0_spatial <- function() {
  run_pilot <- as_flag(Sys.getenv("RENEWAL_V3_0_RUN_PILOT", "false"))
  min_cases_env <- Sys.getenv("RENEWAL_V3_0_MIN_CASES", "")
  min_total_cases <- if (nzchar(min_cases_env)) as.numeric(min_cases_env) else NULL
  if (run_pilot && !is.null(min_total_cases)) {
    stop(
      "RENEWAL_V3_0_RUN_PILOT (ad-hoc debug subset) and RENEWAL_V3_0_MIN_CASES ",
      "(case-count threshold subset) are mutually exclusive"
    )
  }
  run_sampling <- run_pilot || !is.null(min_total_cases)

  pilot_iter <- as.integer(Sys.getenv("RENEWAL_V3_0_PILOT_ITER", "100"))
  pilot_warmup <- as.integer(Sys.getenv("RENEWAL_V3_0_PILOT_WARMUP", "50"))
  pilot_chains <- as.integer(Sys.getenv("RENEWAL_V3_0_PILOT_CHAINS", "1"))
  adapt_delta <- as.numeric(Sys.getenv("RENEWAL_V3_0_ADAPT_DELTA", "0.90"))
  max_treedepth <- as.integer(Sys.getenv("RENEWAL_V3_0_MAX_TREEDEPTH", "12"))
  if (is.na(pilot_iter) || pilot_iter < 20L || is.na(pilot_warmup) ||
    pilot_warmup < 1L || pilot_warmup >= pilot_iter) {
    stop("Iteration settings require iter >= 20 and 1 <= warmup < iter")
  }
  if (is.na(pilot_chains) || pilot_chains < 1L) {
    stop("RENEWAL_V3_0_PILOT_CHAINS must be a positive integer")
  }

  # The full all-Ceara data are prepared and checked by default, but no
  # sampling is started unless a selection mode above requests it.
  prepared <- prepare_renewal_v3_0_spatial_data(
    debug_subset = run_pilot, min_total_cases = min_total_cases
  )
  stan_file <- here::here("02_Script/stan/renewal_ceara_v3_0_spatial.stan")
  model <- rstan::stan_model(file = stan_file)
  message("[compile] successful: ", stan_file)

  pilot <- NULL
  if (run_sampling) {
    data <- prepared$stan_data
    message(sprintf(
      "[sampling] M=%d municipalities, %d chain(s), iter=%d, warmup=%d, adapt_delta=%.2f, max_treedepth=%d",
      data$M, pilot_chains, pilot_iter, pilot_warmup, adapt_delta, max_treedepth
    ))
    init_fn <- function() {
      list(
        sigma_R = 0.05,
        mu_R = rep(log(1.1), data$N - data$seed_weeks),
        sigma_muni_R = 0.05,
        z_muni_R = rep(0, data$M),
        log_seed_hazard = log(10) - log(data$N_start[, seq_len(data$seed_weeks), drop = FALSE]),
        p_symp = 0.52,
        rho_sym_global = 0.25,
        phi_obs = 20
      )
    }
    fit <- rstan::sampling(
      object = model,
      data = data,
      iter = pilot_iter,
      warmup = pilot_warmup,
      chains = pilot_chains,
      seed = 30052018,
      init = init_fn,
      control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
      refresh = max(1L, floor(pilot_iter / 5L))
    )
    pilot <- pilot_diagnostics(fit)
    pilot_suffix <- if (run_pilot) {
      "_debug_pilot"
    } else {
      sprintf("_min%d_cases", min_total_cases)
    }
    pilot_path <- here::here(
      "02_Script/stan", paste0("renewal_ceara_v3_0_spatial_fit", pilot_suffix, ".rds")
    )
    saveRDS(list(fit = fit, prepared = prepared, pilot_diagnostics = pilot), pilot_path)
    message("[pilot save] ", pilot_path)
  }

  report_path <- write_v3_0_implementation_report(
    prepared,
    compiled = TRUE, pilot = pilot, sampling_chains = pilot_chains
  )
  message("[report] ", report_path)
  invisible(list(model = model, prepared = prepared, pilot = pilot, report_path = report_path))
}

if (sys.nframe() == 0) {
  run_v3_0_spatial()
}
