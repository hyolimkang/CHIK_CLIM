# Prior-predictive calibration of kappa_sero=50 (task Section 8), BEFORE
# any fitting. kappa_sero is a FIXED, pragmatic observation-dispersion
# parameter -- not a biological quantity, not estimated. This reports what
# range of OBSERVED Juazeiro do Norte seroprevalences (n=404) the
# beta-binomial(404, p*kappa, (1-p)*kappa) implies for several assumed true
# Ceará STATE prevalences p, so the choice can be judged for plausibility
# before use -- not tuned after seeing how the renewal-model posterior
# behaves.

required_packages <- c("here", "tibble", "readr", "dplyr", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(tibble); library(readr); library(dplyr); library(scales) })

kappa_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

# Beta-binomial(n, alpha, beta) via a Beta-mixed-Binomial Monte Carlo
# (exact and simple; avoids a VGAM/extraDistr dependency).
rbetabinom <- function(n_draws, size, alpha, beta) {
  p <- rbeta(n_draws, alpha, beta)
  rbinom(n_draws, size, p)
}

run_kappa_calibration <- function() {
  root <- kappa_root()
  out_dir <- file.path(root, "02_Script", "40_renewal_model", "21_v4_2_q_juazeiro_serology")
  sero_n <- 404L
  kappa_sero <- 50
  true_state_prevalences <- c(.10, .20, .25, .30, .40)

  set.seed(20260911L)
  n_mc <- 2e5
  results <- bind_rows(lapply(true_state_prevalences, function(p) {
    alpha <- p * kappa_sero; beta <- (1 - p) * kappa_sero
    draws <- rbetabinom(n_mc, sero_n, alpha, beta)
    tibble(
      assumed_true_state_prevalence = p, kappa_sero = kappa_sero, alpha_sero = alpha, beta_sero = beta,
      implied_observed_x_median = median(draws), implied_observed_x_q025 = quantile(draws, .025, names = FALSE),
      implied_observed_x_q975 = quantile(draws, .975, names = FALSE),
      implied_observed_prevalence_q025 = quantile(draws, .025, names = FALSE) / sero_n,
      implied_observed_prevalence_median = median(draws) / sero_n,
      implied_observed_prevalence_q975 = quantile(draws, .975, names = FALSE) / sero_n
    )
  }))
  write_csv(results, file.path(out_dir, "kappa_sero_calibration.csv"))

  message("[kappa calibration] kappa_sero=", kappa_sero, ", sero_n=", sero_n, " (Juazeiro do Norte, 103/404 observed):")
  print(as.data.frame(results))
  message("\n[kappa calibration] Interpretation: at kappa_sero=50, if the TRUE Ceará-state prevalence were exactly ",
          "the OBSERVED Juazeiro value's own rate (25.5%), the beta-binomial's 95% predictive range for a repeat ",
          "404-person Juazeiro-like sample is [", round(results$implied_observed_x_q025[results$assumed_true_state_prevalence == .25]), ", ",
          round(results$implied_observed_x_q975[results$assumed_true_state_prevalence == .25]), "] out of 404 -- i.e. observed ",
          "prevalence could plausibly range from ", scales::percent(results$implied_observed_prevalence_q025[results$assumed_true_state_prevalence == .25], accuracy = .1),
          " to ", scales::percent(results$implied_observed_prevalence_q975[results$assumed_true_state_prevalence == .25], accuracy = .1),
          " even if the true state prevalence were exactly 25%. This is the deliberate representativeness-weakening ",
          "the beta-binomial (vs a naive binomial) is meant to encode: a single municipality's survey should NOT pin ",
          "state prevalence down to within a percentage point or two, as a plain binomial(404, p) would.")
  invisible(results)
}

if (sys.nframe() == 0L) run_kappa_calibration()
