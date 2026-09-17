# Applies the prespecified v4.1 decision rule by reading the outputs already
# written by scripts 01-05. Does not refit or recompute anything itself.
#
# v4.1 is a STABILIZED HISTORICAL NO-VACCINATION BACKBONE only if:
#   1. HMC passes (primary run);
#   2. simulation recovery is acceptable;
#   3. q is not pathologically confounded with transmission parameters;
#   4. major epidemic periods are reasonably reproduced;
#   5. q-prior sensitivity does not make S essentially arbitrary.
# If HMC passes but q/S remain almost entirely prior-determined, the verdict
# is reported as "computationally stable but absolute susceptibility weakly
# identified" rather than a pass/fail label.

required_packages <- c("here", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr) })

decision_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

read_if_exists <- function(path) if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL

run_v4_1_decision_rule <- function() {
  root <- decision_root()
  table_dir <- file.path(root, "03_Output", "tables", "renewal_v4_1_minimal_no_vaccine")

  hmc <- read_if_exists(file.path(table_dir, "hmc_diagnostics.csv"))
  correlations <- read_if_exists(file.path(table_dir, "q_posterior_correlations.csv"))
  recovery_summary <- read_if_exists(file.path(table_dir, "simulation_recovery_summary.csv"))
  mismatch <- read_if_exists(file.path(table_dir, "annual_mismatch_comparison_v4_0_vs_v4_1.csv"))
  stability <- read_if_exists(file.path(table_dir, "q_prior_sensitivity_stability_verdict.csv"))
  q_summary <- read_if_exists(file.path(table_dir, "q_and_sigma_year_prior_vs_posterior.csv"))

  missing <- c(
    if (is.null(hmc)) "hmc_diagnostics.csv (run 02_diagnose)",
    if (is.null(correlations)) "q_posterior_correlations.csv (run 02_diagnose)",
    if (is.null(recovery_summary)) "simulation_recovery_summary.csv (run 03_simulation_recovery)",
    if (is.null(stability)) "q_prior_sensitivity_stability_verdict.csv (run 05_compare_q_prior_sensitivity)"
  )
  if (length(missing)) {
    message("[v4.1 decision] cannot yet render a final verdict -- missing: ", paste(missing, collapse = "; "))
    return(invisible(NULL))
  }

  criterion_1_hmc <- hmc$value[hmc$metric == "divergences"] == 0 &&
    hmc$value[hmc$metric == "max_treedepth_hits"] == 0 &&
    hmc$value[hmc$metric == "maximum_Rhat"] <= 1.01 &&
    hmc$value[hmc$metric == "minimum_bulk_ESS"] >= 100

  criterion_2_recovery <- recovery_summary$value[recovery_summary$metric == "q_and_S_95pct_coverage"] >= 0.8 &&
    recovery_summary$value[recovery_summary$metric == "HMC_gate_pass"] == 1

  criterion_3_confounding <- !any(correlations$pathological_confounding_flag)

  criterion_4_major_epidemics <- if (!is.null(mismatch)) {
    rows <- filter(mismatch, year %in% c(2017L, 2022L))
    all(abs(rows$v4_1_relative_residual) < 0.5)
  } else NA

  criterion_5_sprior_stability <- all(stability$qualitatively_stable)

  q_prior_data_conflict <- if (!is.null(q_summary)) {
    row <- filter(q_summary, parameter == "q")
    abs(row$posterior_minus_prior_mean_in_prior_sd) < 0.5 # posterior barely moved from prior
  } else NA

  all_pass <- isTRUE(criterion_1_hmc) && isTRUE(criterion_2_recovery) && isTRUE(criterion_3_confounding) &&
    isTRUE(criterion_4_major_epidemics) && isTRUE(criterion_5_sprior_stability)

  weakly_identified <- isTRUE(criterion_1_hmc) && isTRUE(q_prior_data_conflict) && !isTRUE(criterion_5_sprior_stability)

  verdict <- if (all_pass) {
    "STABILIZED HISTORICAL NO-VACCINATION BACKBONE: all five criteria satisfied."
  } else if (weakly_identified) {
    "COMPUTATIONALLY STABLE BUT ABSOLUTE SUSCEPTIBILITY WEAKLY IDENTIFIED: HMC passes and q/S are essentially prior-determined -- do not present S(t) as data-driven without this caveat."
  } else {
    "NOT YET A STABILIZED BACKBONE: at least one decision-rule criterion failed -- see criteria table for which."
  }

  criteria_table <- tibble::tibble(
    criterion = c("1_HMC_gate_passes", "2_simulation_recovery_acceptable", "3_q_not_pathologically_confounded",
                  "4_major_epidemics_reasonably_reproduced", "5_q_prior_sensitivity_not_arbitrary"),
    result = c(criterion_1_hmc, criterion_2_recovery, criterion_3_confounding, criterion_4_major_epidemics, criterion_5_sprior_stability)
  )
  write_csv(criteria_table, file.path(table_dir, "v4_1_decision_rule_criteria.csv"))
  writeLines(c(
    "# v4.1 decision-rule verdict", "",
    paste0("**", verdict, "**"), "",
    "## Criteria", "",
    paste0("- HMC gate passes: ", criterion_1_hmc),
    paste0("- Simulation recovery acceptable (q & S coverage >= 80%, HMC pass): ", criterion_2_recovery),
    paste0("- q not pathologically confounded with transmission parameters (|r| <= 0.9): ", criterion_3_confounding),
    paste0("- Major epidemic periods (2017, 2022) reasonably reproduced (|relative residual| < 50%): ", criterion_4_major_epidemics),
    paste0("- q-prior sensitivity does not make S essentially arbitrary (<25% swing across priors): ", criterion_5_sprior_stability),
    paste0("- q posterior essentially unmoved from prior (< 0.5 prior SD shift): ", q_prior_data_conflict)
  ), file.path(table_dir, "v4_1_decision_rule_verdict.md"))
  message("[v4.1 decision] ", verdict)
  invisible(list(verdict = verdict, criteria = criteria_table))
}

if (sys.nframe() == 0L) run_v4_1_decision_rule()
