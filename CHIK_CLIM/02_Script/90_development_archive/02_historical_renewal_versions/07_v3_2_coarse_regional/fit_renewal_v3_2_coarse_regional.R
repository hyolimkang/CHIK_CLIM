# Compile v3.2 and optionally run only the prescribed diagnostic pilot.
required_packages <- c("here", "rstan", "posterior")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(posterior)
})
rstan_options(auto_write = TRUE)
options(mc.cores = min(4, parallel::detectCores()))
source(here::here("02_Script/90_development_archive/02_historical_renewal_versions/07_v3_2_coarse_regional/prepare_renewal_v3_2_coarse_regional_data.R"))
as_flag <- function(x) tolower(x) %in% c("1", "true", "yes")
param_count_v3_2 <- function(d) 2L + (d$N - 1L) + 1L + d$K + 2L + d$K + 2L
hmc_v3_2 <- function(fit) {
  a <- posterior::as_draws_array(rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE))
  s <- as.data.frame(posterior::summarise_draws(a, posterior::rhat, posterior::ess_bulk, posterior::ess_tail))
  names(s) <- sub("^posterior::", "", names(s))
  if (".variable" %in% names(s)) names(s)[names(s) == ".variable"] <- "variable"
  fr <- s$rhat[is.finite(s$rhat)]
  list(
    divergences = rstan::get_num_divergent(fit), treedepth_hits = rstan::get_num_max_treedepth(fit),
    bfmi = rstan::get_bfmi(fit), max_rhat = if (length(fr)) max(fr) else NA_real_,
    rhat_gt_1_01 = sum(s$rhat > 1.01, na.rm = TRUE),
    min_bulk_ess = min(s$ess_bulk, na.rm = TRUE), min_tail_ess = min(s$ess_tail, na.rm = TRUE),
    convergence = s
  )
}
write_v3_2_report <- function(x, compiled, hmc = NULL, elapsed = NULL) {
  path <- normalizePath(here::here("..", "docs", "v3_2_coarse_regional_implementation.md"), winslash = "/", mustWork = FALSE)
  d <- x$diagnostics
  q <- x$q_prior
  h <- if (is.null(hmc)) {
    c("- Pilot not run.")
  } else {
    c(
      sprintf("- Elapsed: %.1f seconds; divergences: %d; treedepth hits: %d.", elapsed, hmc$divergences, hmc$treedepth_hits),
      sprintf(
        "- max Rhat: %.3f; Rhat > 1.01: %d; min bulk/tail ESS: %.1f / %.1f; BFMI by chain: %s.",
        hmc$max_rhat, hmc$rhat_gt_1_01, hmc$min_bulk_ess, hmc$min_tail_ess, paste(sprintf("%.3f", hmc$bfmi), collapse = ", ")
      ),
      "- Pilot status is assessed by the diagnostic script; no scientific inference follows from this pilot."
    )
  }
  lines <- c(
    "# v3.2 coarse-regional renewal implementation", "", "## Definition", "",
    "- Eight contiguous official-mapping units: Fortaleza, its remainder, Norte/Sobral, Cariri remainder, Juazeiro do Norte, Sert<U+00E3>o Central remainder, Quixad<U+00E1>, and Litoral Leste/Jaguaribe.",
    "- Municipality membership and original official macroregion are saved in the assignment CSV. Juazeiro, Quixad<U+00E1>, and Fortaleza are removed from parent regions before aggregation.",
    sprintf("- Assignment CSV: %s.", x$assignment_path), "", "## Differences", "",
    "- Replaces v3.1 burden strata with official contiguous health-macroregion units.",
    "- Uses direct non-centred shared log-R0 random walk, not free weekly hazards.",
    "- Uses one hierarchical seed intensity per region and a fixed equal eight-week seed profile.",
    "- Replaces p_symp times reporting with one q parameter.", "", "## q prior", "",
    sprintf("- %d simulated products of Beta(30,28) and Beta(20,60), seed %d: mean %.6f, SD %.6f.", q$n_draws, q$seed, q$mean, q$sd),
    sprintf("- Moment-matched q prior: Beta(%.6f, %.6f).", q$a, q$b), "", "## Checks", "",
    sprintf("- %d municipalities assigned exactly once to %d regions; %d weeks.", d$municipalities, d$units, d$weeks),
    sprintf("- Maximum deterministic S+U accounting error %.3g; population-gap tolerance result %.3f%%.", d$accounting_error, 100 * d$population_gap),
    "- Weekly regional cases exactly reproduce the Cear<U+00E1> total; births/deaths are non-negative; survey windows map directly to Juazeiro and Quixad<U+00E1>.",
    sprintf("- Parameter dimension: %d. Stan compilation: %s.", param_count_v3_2(x$stan_data), if (compiled) "successful" else "not completed"),
    "", "## Pilot diagnostics", "", h, "", "## Warnings", "",
    "- Municipality deaths retain the explicit v3.0 population-share allocation of Cear<U+00E1> annual deaths before aggregation.",
    "- Regional R0 is a nuisance device for estimating statewide susceptibility, not the primary target.",
    "- Good case fit alone is not scientific validation."
  )
  writeLines(lines, path)
  path
}
run_v3_2 <- function() {
  run_pilot <- as_flag(Sys.getenv("RENEWAL_V3_2_RUN_PILOT", "false"))
  iter <- as.integer(Sys.getenv("RENEWAL_V3_2_ITER", "500"))
  warmup <- as.integer(Sys.getenv("RENEWAL_V3_2_WARMUP", "250"))
  if (is.na(iter) || is.na(warmup) || warmup < 1 || warmup >= iter) stop("Require 1 <= warmup < iter")
  x <- prepare_renewal_v3_2_coarse_regional_data()
  model <- rstan::stan_model(here::here("02_Script/stan/renewal_ceara_v3_2_coarse_regional.stan"))
  h <- NULL
  elapsed <- NULL
  if (run_pilot) {
    d <- x$stan_data
    init <- function() {
      list(
        alpha_R = log(1.1), sigma_R = .05, z_R = rep(0, d$N - 1),
        sigma_region = .05, z_region = rep(0, d$K), mu_seed = log(10), sigma_seed = .1, z_seed = rep(0, d$K),
        q = d$q_prior_a / (d$q_prior_a + d$q_prior_b), phi_obs = 20
      )
    }
    timing <- system.time(fit <- rstan::sampling(model,
      data = d, chains = 4, iter = iter, warmup = warmup,
      seed = 32052018, init = init, control = list(adapt_delta = .90, max_treedepth = 12), refresh = max(1L, floor(iter / 10))
    ))
    elapsed <- unname(timing["elapsed"])
    h <- hmc_v3_2(fit)
    saveRDS(
      list(fit = fit, prepared = x, pilot_diagnostics = h, elapsed_seconds = elapsed),
      here::here("02_Script/stan/renewal_ceara_v3_2_coarse_regional_pilot.rds")
    )
  }
  report <- write_v3_2_report(x, TRUE, h, elapsed)
  message("[report] ", report)
  invisible(list(prepared = x, diagnostics = h, report = report))
}
if (sys.nframe() == 0) run_v3_2()
