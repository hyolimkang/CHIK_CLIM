# ARM A -- no vaccine, closed-loop historical replay (feedback_on = TRUE).
run_arm_A <- function(inputs, R0_draw) {
  simulate_age_vaccine(
    R0_t_draw = R0_draw, births = inputs$births, deaths = inputs$deaths, reconciliation = inputs$reconciliation,
    N_start = inputs$N_start, w = inputs$w, imports_per_week = inputs$imports_per_week,
    is_seed = inputs$is_seed, X_seed = inputs$X_seed, age_shares = inputs$age_shares, years = inputs$years,
    max_age = MAX_AGE, vaccination = NULL, feedback_on = TRUE, lambda_external = NULL
  )
}
