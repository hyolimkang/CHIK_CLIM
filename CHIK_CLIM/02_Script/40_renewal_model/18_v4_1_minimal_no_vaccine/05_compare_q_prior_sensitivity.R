# Prespecified q-prior sensitivity comparison (run only after the primary
# v4.1 fit passes the HMC gate, per the task spec).
#
# Two additional priors are fit with the SAME v4.1 model/data/init strategy,
# preserving the same prior mean (~0.129) as the primary Beta(16.2, 108.7):
#   - narrow: Beta(32.4, 217.4)  (primary (a,b) x2 -- half the prior variance)
#   - wide:   Beta(8.1, 54.35)   (primary (a,b) x0.5 -- double the prior variance)
# These are NOT chosen based on which gives the most attractive S trajectory;
# they are a fixed multiplicative rescaling of the primary prior's (a,b)
# decided before any of the three fits were compared.
#
# Run (after 01_fit has been run once with defaults for the primary tag):
#   RENEWAL_V4_1_RUN_TAG=q_prior_narrow RENEWAL_V4_1_Q_PRIOR_A=32.4 RENEWAL_V4_1_Q_PRIOR_B=217.4 Rscript 01_fit_v4_1_minimal_no_vaccine.R
#   RENEWAL_V4_1_RUN_TAG=q_prior_wide   RENEWAL_V4_1_Q_PRIOR_A=8.1  RENEWAL_V4_1_Q_PRIOR_B=54.35  Rscript 01_fit_v4_1_minimal_no_vaccine.R
# then run this script with default run_tag (reads all three fits directly).

required_packages <- c("here", "rstan", "dplyr", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr) })

sensitivity_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

extract_scenario_summary <- function(root, run_tag, label) {
  suffix <- if (identical(run_tag, "base")) "" else paste0("_", run_tag)
  fit_path <- file.path(root, "02_Script", "stan", paste0("renewal_ceara_v4_1_minimal_no_vaccine_fit", suffix, ".rds"))
  if (!file.exists(fit_path)) {
    message("[q-prior sensitivity] missing fit for scenario '", label, "': ", fit_path)
    return(NULL)
  }
  bundle <- readRDS(fit_path)
  fit <- bundle$fit
  years <- as.integer(format(as.Date(bundle$weekly_data$week_start), "%Y"))
  q_draws <- rstan::extract(fit, pars = "q")$q
  alpha_R_draws <- rstan::extract(fit, pars = "alpha_R")$alpha_R
  year_effect_draws <- rstan::extract(fit, pars = "year_effect")$year_effect
  X_draws <- rstan::extract(fit, pars = "X")$X
  S_draws <- rstan::extract(fit, pars = "S")$S
  sorted_years <- sort(unique(years))
  r0_year <- function(y) {
    i <- match(y, sorted_years)
    exp(alpha_R_draws + year_effect_draws[, i])
  }
  tibble(
    scenario = label, q_prior_a = bundle$stan_data$q_prior_a, q_prior_b = bundle$stan_data$q_prior_b,
    q_prior_mean = bundle$stan_data$q_prior_a / (bundle$stan_data$q_prior_a + bundle$stan_data$q_prior_b),
    hmc_pass = bundle$hmc$hmc_pass,
    q_posterior_median = median(q_draws), q_posterior_q025 = quantile(q_draws, .025), q_posterior_q975 = quantile(q_draws, .975),
    S_2025_prop_median = median(S_draws[, which(years == 2025L)[1]] / bundle$stan_data$N_start[which(years == 2025L)[1]]),
    S_2025_prop_q025 = quantile(S_draws[, which(years == 2025L)[1]] / bundle$stan_data$N_start[which(years == 2025L)[1]], .025),
    S_2025_prop_q975 = quantile(S_draws[, which(years == 2025L)[1]] / bundle$stan_data$N_start[which(years == 2025L)[1]], .975),
    cumulative_infections_median = median(rowSums(X_draws)), cumulative_infections_q025 = quantile(rowSums(X_draws), .025), cumulative_infections_q975 = quantile(rowSums(X_draws), .975),
    R0_2017_median = median(r0_year(2017L)), R0_2017_q025 = quantile(r0_year(2017L), .025), R0_2017_q975 = quantile(r0_year(2017L), .975),
    R0_2022_median = median(r0_year(2022L)), R0_2022_q025 = quantile(r0_year(2022L), .025), R0_2022_q975 = quantile(r0_year(2022L), .975)
  )
}

run_compare_q_prior_sensitivity <- function() {
  root <- sensitivity_root()
  scenarios <- list(c("base", "primary_Beta_16.2_108.7"), c("q_prior_narrow", "narrow_Beta_32.4_217.4"), c("q_prior_wide", "wide_Beta_8.1_54.35"))
  rows <- lapply(scenarios, function(s) extract_scenario_summary(root, s[1], s[2]))
  comparison <- bind_rows(rows)
  if (nrow(comparison) < 2) {
    message("[q-prior sensitivity] fewer than 2 scenarios available -- nothing to compare yet.")
    return(invisible(comparison))
  }

  # Qualitative-stability check: does S_2025 / cumulative infections / R0_2017
  # / R0_2022 stay within the SAME order of magnitude and same qualitative
  # regime (e.g. R0 > 1 or < 1) across all available scenarios?
  stability <- tibble(
    quantity = c("S_2025_prop", "cumulative_infections", "R0_2017", "R0_2022"),
    range_ratio = c(
      max(comparison$S_2025_prop_median) / min(comparison$S_2025_prop_median),
      max(comparison$cumulative_infections_median) / min(comparison$cumulative_infections_median),
      max(comparison$R0_2017_median) / min(comparison$R0_2017_median),
      max(comparison$R0_2022_median) / min(comparison$R0_2022_median)
    )
  ) |> mutate(qualitatively_stable = range_ratio < 1.25) # <25% swing in point estimate across priors

  paths_tables <- file.path(root, "03_Output", "tables", "renewal_v4_1_minimal_no_vaccine")
  dir.create(paths_tables, recursive = TRUE, showWarnings = FALSE)
  write_csv(comparison, file.path(paths_tables, "q_prior_sensitivity_comparison.csv"))
  write_csv(stability, file.path(paths_tables, "q_prior_sensitivity_stability_verdict.csv"))
  message("[q-prior sensitivity] stability verdict:")
  for (i in seq_len(nrow(stability))) message(sprintf("  %s: range_ratio(max/min median)=%.3f -> %s", stability$quantity[i], stability$range_ratio[i], if (stability$qualitatively_stable[i]) "stable" else "NOT stable"))
  invisible(list(comparison = comparison, stability = stability))
}

if (sys.nframe() == 0L) run_compare_q_prior_sensitivity()
