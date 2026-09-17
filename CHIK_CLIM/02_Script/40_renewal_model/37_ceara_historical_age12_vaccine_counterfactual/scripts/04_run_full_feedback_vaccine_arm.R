# ARM B -- routine age-12 vaccination, FULL transmission feedback
# (feedback_on = TRUE): the vaccine arm generates its OWN infectiousness/
# FOI recursively from its OWN (vaccination-reduced) X_total history.
# This is what produces the indirect (herd) effect.
run_arm_B <- function(inputs, R0_draw, coverage, ve_infection) {
  coverage_by_year <- setNames(rep(coverage, length(unique(inputs$years))), sort(unique(inputs$years)))
  vaccination <- list(target_age = TARGET_AGE, coverage_by_year = coverage_by_year, ve_infection = ve_infection)
  simulate_age_vaccine(
    R0_t_draw = R0_draw, births = inputs$births, deaths = inputs$deaths, reconciliation = inputs$reconciliation,
    N_start = inputs$N_start, w = inputs$w, imports_per_week = inputs$imports_per_week,
    is_seed = inputs$is_seed, X_seed = inputs$X_seed, age_shares = inputs$age_shares, years = inputs$years,
    max_age = MAX_AGE, vaccination = vaccination, feedback_on = TRUE, lambda_external = NULL
  )
}
