# ARM C -- routine age-12 vaccination applied, but transmission feedback
# DISABLED (feedback_on = FALSE): lambda(t) is fixed to Arm A's no-vaccine
# trajectory. Vaccinated age-12 individuals are still protected (direct
# effect), but their averted infections do NOT reduce future
# infectiousness/FOI, so non-target age groups should see ~zero benefit.
# NOT a physical policy counterfactual -- a diagnostic decomposition tool
# only (Section 5 of the task specification).
run_arm_C <- function(inputs, R0_draw, coverage, ve_infection, lambda_external) {
  coverage_by_year <- setNames(rep(coverage, length(unique(inputs$years))), sort(unique(inputs$years)))
  vaccination <- list(target_age = TARGET_AGE, coverage_by_year = coverage_by_year, ve_infection = ve_infection)
  simulate_age_vaccine(
    R0_t_draw = R0_draw, births = inputs$births, deaths = inputs$deaths, reconciliation = inputs$reconciliation,
    N_start = inputs$N_start, w = inputs$w, imports_per_week = inputs$imports_per_week,
    is_seed = inputs$is_seed, X_seed = inputs$X_seed, age_shares = inputs$age_shares, years = inputs$years,
    max_age = MAX_AGE, vaccination = vaccination, feedback_on = FALSE, lambda_external = lambda_external
  )
}
