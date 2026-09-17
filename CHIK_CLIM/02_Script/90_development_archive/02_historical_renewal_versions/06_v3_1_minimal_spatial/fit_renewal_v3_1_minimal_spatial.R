# =============================================================================
# fit_renewal_v3_1_minimal_spatial.R
#
# Default action: build and validate the six-unit 2015--2019 data, then
# compile Stan. A short four-chain pilot is opt-in and is the only sampling
# mode implemented here; no production HMC is launched automatically.
#
# Environment switches:
#   RENEWAL_V3_1_RUN_PILOT=true
#   RENEWAL_V3_1_PILOT_ITER=100
#   RENEWAL_V3_1_PILOT_WARMUP=50
#   RENEWAL_V3_1_HIGH_BURDEN_MIN=250
#   RENEWAL_V3_1_MODERATE_BURDEN_MIN=25
# =============================================================================

required_packages <- c("here", "rstan", "posterior")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(posterior)
})

rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))

source(here::here(
  "02_Script/90_development_archive/02_historical_renewal_versions/06_v3_1_minimal_spatial/prepare_renewal_v3_1_minimal_spatial_data.R"
))

as_flag <- function(value) tolower(value) %in% c("1", "true", "yes")

v3_1_parameter_count <- function(data) {
  # sigma_R, sigma_unit, p_symp, rho_sym, phi_obs plus vector/matrix terms.
  5L + (data$N - data$seed_weeks) + data$K + data$K * data$seed_weeks
}

v3_1_hmc_diagnostics <- function(fit) {
  draws_array <- posterior::as_draws_array(
    rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE)
  )
  convergence <- posterior::summarise_draws(
    draws_array,
    posterior::rhat,
    posterior::ess_bulk,
    posterior::ess_tail
  )
  convergence <- as.data.frame(convergence)
  names(convergence) <- sub("^posterior::", "", names(convergence))
  if (".variable" %in% names(convergence)) {
    names(convergence)[names(convergence) == ".variable"] <- "variable"
  }
  finite_rhat <- convergence$rhat[is.finite(convergence$rhat)]
  worst <- convergence[order(convergence$rhat,
    decreasing = TRUE,
    na.last = NA
  ), , drop = FALSE]
  list(
    divergences = rstan::get_num_divergent(fit),
    max_treedepth_hits = rstan::get_num_max_treedepth(fit),
    bfmi = rstan::get_bfmi(fit),
    max_rhat = if (length(finite_rhat)) max(finite_rhat) else NA_real_,
    min_bulk_ess = min(convergence$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(convergence$ess_tail, na.rm = TRUE),
    convergence = convergence,
    worst_mixing = utils::head(worst, 10)
  )
}

write_v3_1_implementation_report <- function(prepared, compiled, pilot = NULL,
                                             elapsed_seconds = NULL) {
  report_path <- here::here(
    "03_Output/tables/renewal_v3_1/renewal_v3_1_minimal_spatial_implementation_report.md"
  )
  dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

  d <- prepared$diagnostics
  x <- prepared$stan_data
  pilot_lines <- if (is.null(pilot)) {
    c("- Four-chain pilot: not run; default execution stops after compilation.")
  } else {
    c(
      sprintf("- Four-chain pilot elapsed time: %.1f seconds.", elapsed_seconds),
      sprintf(
        "- Divergences: %d; treedepth hits: %d; maximum Rhat: %.3f.",
        pilot$divergences, pilot$max_treedepth_hits, pilot$max_rhat
      ),
      sprintf(
        "- Minimum bulk ESS: %.1f; minimum tail ESS: %.1f; minimum BFMI: %.3f.",
        pilot$min_bulk_ess, pilot$min_tail_ess, min(pilot$bfmi)
      ),
      "- This short pilot is a computational diagnostic, not a converged production analysis."
    )
  }

  lines <- c(
    "# v3.1 minimal spatial renewal implementation report",
    "",
    "## Model structure",
    "",
    "- Six units: Fortaleza, Juazeiro do Norte, Quixada, and fixed high/moderate/low burden pooled strata.",
    sprintf(
      "- Fixed thresholds: high >= %g total cases; moderate %g to < %g; low < %g.",
      prepared$config$high_burden_min, prepared$config$moderate_burden_min,
      prepared$config$high_burden_min, prepared$config$moderate_burden_min
    ),
    "- The named municipalities are excluded from pooled strata before threshold assignment.",
    "- Unit-specific S/U/X states use the v2.2 demographic recursion; statewide susceptibility is population-weighted and derived.",
    "- Serology likelihood includes only Juazeiro and Quixada. Fortaleza is external only.",
    "- Reporting is common across units: p_symp * rho_sym.",
    "",
    "## Checks",
    "",
    sprintf(
      "- %d municipalities collapsed to %d units across %d weeks.",
      d$municipalities_all_ceara, d$units, d$weeks
    ),
    sprintf(
      "- Maximum deterministic S + U accounting error: %.3g.",
      d$maximum_accounting_error
    ),
    sprintf(
      "- Maximum municipality-sum versus v2.2 state population gap: %.3f%%.",
      100 * d$maximum_population_relative_gap
    ),
    sprintf("- Unit assignment: %s.", prepared$assignment_path),
    sprintf("- Stan parameter dimension: %d.", v3_1_parameter_count(x)),
    sprintf("- Stan compilation: %s.", if (compiled) "successful" else "not completed"),
    "",
    "## Pilot diagnostics",
    "",
    pilot_lines,
    "",
    "## Remaining limitations",
    "",
    "- Municipality deaths remain an explicit population-share allocation of the Cear<U+00E1> annual total before aggregation, not observed municipality mortality.",
    "- No serology anchor is created for Fortaleza or pooled units; their uncertainty should remain data-driven.",
    "- Reporting, seeding, and susceptibility can remain confounded. The intended future test is robustness of the Cear<U+00E1> trajectory to refined pooling, not a preferred point estimate."
  )
  writeLines(lines, report_path, useBytes = TRUE)
  report_path
}

run_v3_1_minimal_spatial <- function() {
  run_pilot <- as_flag(Sys.getenv("RENEWAL_V3_1_RUN_PILOT", "false"))
  iter <- as.integer(Sys.getenv("RENEWAL_V3_1_PILOT_ITER", "100"))
  warmup <- as.integer(Sys.getenv("RENEWAL_V3_1_PILOT_WARMUP", "50"))
  high <- as.numeric(Sys.getenv("RENEWAL_V3_1_HIGH_BURDEN_MIN", "250"))
  moderate <- as.numeric(Sys.getenv("RENEWAL_V3_1_MODERATE_BURDEN_MIN", "25"))
  if (is.na(iter) || is.na(warmup) || iter < 20L || warmup < 1L || warmup >= iter) {
    stop("Pilot settings require iter >= 20 and 1 <= warmup < iter")
  }

  prepared <- prepare_renewal_v3_1_minimal_spatial_data(
    high_burden_min = high,
    moderate_burden_min = moderate
  )
  stan_file <- here::here("02_Script/stan/renewal_ceara_v3_1_minimal_spatial.stan")
  model <- rstan::stan_model(file = stan_file)
  message("[compile] successful: ", stan_file)

  pilot <- NULL
  elapsed_seconds <- NULL
  if (run_pilot) {
    data <- prepared$stan_data
    init_fn <- function() {
      list(
        sigma_R = 0.05,
        mu_R = stats::rnorm(data$N - data$seed_weeks, log(1.1), 0.02),
        sigma_unit = 0.05,
        z_unit = stats::rnorm(data$K, 0, 0.05),
        log_seed_hazard = log(10) -
          log(data$N_start[, seq_len(data$seed_weeks), drop = FALSE]) +
          matrix(stats::rnorm(data$K * data$seed_weeks, 0, 0.02),
            nrow = data$K
          ),
        p_symp = 0.52,
        rho_sym = 0.25,
        phi_obs = 20
      )
    }
    message(sprintf(
      "[pilot] 4 chains; iter=%d; warmup=%d; K=%d; parameter dimension=%d",
      iter, warmup, data$K, v3_1_parameter_count(data)
    ))
    timing <- system.time({
      fit <- rstan::sampling(
        object = model,
        data = data,
        iter = iter,
        warmup = warmup,
        chains = 4,
        seed = 31052018,
        init = init_fn,
        control = list(adapt_delta = 0.90, max_treedepth = 12),
        refresh = max(1L, floor(iter / 5L))
      )
    })
    elapsed_seconds <- unname(timing["elapsed"])
    pilot <- v3_1_hmc_diagnostics(fit)
    out_path <- here::here(
      "02_Script/stan/renewal_ceara_v3_1_minimal_spatial_pilot.rds"
    )
    saveRDS(
      list(
        fit = fit, prepared = prepared, pilot_diagnostics = pilot,
        elapsed_seconds = elapsed_seconds
      ),
      out_path
    )
    message("[pilot save] ", out_path)
  }

  report_path <- write_v3_1_implementation_report(
    prepared,
    compiled = TRUE, pilot = pilot, elapsed_seconds = elapsed_seconds
  )
  message("[report] ", report_path)
  invisible(list(
    model = model, prepared = prepared, pilot = pilot,
    report_path = report_path
  ))
}

if (sys.nframe() == 0) {
  run_v3_1_minimal_spatial()
}
