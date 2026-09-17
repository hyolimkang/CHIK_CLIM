# Execute Arms A/B/C for N_DRAWS paired posterior draws, aggregate to
# reporting age strata immediately (to keep memory bounded), and save
# compact per-draw results to results/. Does not keep full single-year
# age-resolved matrices in memory across draws.

required_packages <- c("here", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble) })

execute_three_arms_all_draws <- function(inputs) {
  ages <- 0:MAX_AGE
  age_group <- cut(ages, breaks = REPORTING_AGE_BREAKS, labels = REPORTING_AGE_LABELS, right = FALSE)
  n_groups <- length(REPORTING_AGE_LABELS)

  set.seed(RANDOM_SEED)
  n_avail <- inputs$n_draws_total
  n_use <- min(N_DRAWS, n_avail)
  draw_idx <- sample.int(n_avail, n_use)
  message(sprintf("[execute] Using %d / %d available posterior draws (seed=%d)", n_use, n_avail, RANDOM_SEED))

  N <- inputs$N
  stresstest <- SCENARIOS$stresstest

  agg_S <- array(NA_real_, dim = c(n_use, N, n_groups, 3), dimnames = list(NULL, NULL, REPORTING_AGE_LABELS, c("A", "B", "C")))
  agg_N <- array(NA_real_, dim = c(n_use, N, n_groups, 3), dimnames = list(NULL, NULL, REPORTING_AGE_LABELS, c("A", "B", "C")))
  agg_Inew <- array(NA_real_, dim = c(n_use, N, n_groups, 3), dimnames = list(NULL, NULL, REPORTING_AGE_LABELS, c("A", "B", "C")))

  totals <- vector("list", n_use)
  vax_bookkeeping <- vector("list", n_use)
  qa_rows <- vector("list", n_use)

  for (i in seq_len(n_use)) {
    d <- draw_idx[i]
    R0_draw <- inputs$R0_t_all_draws[d, ]

    sim_A <- run_arm_A(inputs, R0_draw)
    sim_B <- run_arm_B(inputs, R0_draw, coverage = stresstest$coverage, ve_infection = stresstest$ve_infection)
    sim_C <- run_arm_C(inputs, R0_draw, coverage = stresstest$coverage, ve_infection = stresstest$ve_infection, lambda_external = sim_A$lambda)

    for (g in seq_len(n_groups)) {
      cols <- which(age_group == REPORTING_AGE_LABELS[g])
      agg_S[i, , g, "A"] <- rowSums(sim_A$S_age[, cols, drop = FALSE]); agg_N[i, , g, "A"] <- rowSums(sim_A$S_age[, cols, drop = FALSE] + sim_A$Uinf_age[, cols, drop = FALSE] + sim_A$Uvac_age[, cols, drop = FALSE]); agg_Inew[i, , g, "A"] <- rowSums(sim_A$Inew_age[, cols, drop = FALSE])
      agg_S[i, , g, "B"] <- rowSums(sim_B$S_age[, cols, drop = FALSE]); agg_N[i, , g, "B"] <- rowSums(sim_B$S_age[, cols, drop = FALSE] + sim_B$Uinf_age[, cols, drop = FALSE] + sim_B$Uvac_age[, cols, drop = FALSE]); agg_Inew[i, , g, "B"] <- rowSums(sim_B$Inew_age[, cols, drop = FALSE])
      agg_S[i, , g, "C"] <- rowSums(sim_C$S_age[, cols, drop = FALSE]); agg_N[i, , g, "C"] <- rowSums(sim_C$S_age[, cols, drop = FALSE] + sim_C$Uinf_age[, cols, drop = FALSE] + sim_C$Uvac_age[, cols, drop = FALSE]); agg_Inew[i, , g, "C"] <- rowSums(sim_C$Inew_age[, cols, drop = FALSE])
    }

    totals[[i]] <- tibble(
      draw = d, t = seq_len(N), week_start = inputs$dates,
      S_A = sim_A$S_total, U_A = sim_A$Uinf_total, X_A = sim_A$X_total, lambda_A = sim_A$lambda, Reff_A = sim_A$R_eff, infectiousness_A = sim_A$infectiousness, N_A = sim_A$N_total,
      S_B = sim_B$S_total, U_B = sim_B$Uinf_total, Uvac_B = sim_B$Uvac_total, X_B = sim_B$X_total, lambda_B = sim_B$lambda, Reff_B = sim_B$R_eff, infectiousness_B = sim_B$infectiousness, N_B = sim_B$N_total,
      S_C = sim_C$S_total, U_C = sim_C$Uinf_total, Uvac_C = sim_C$Uvac_total, X_C = sim_C$X_total, lambda_C = sim_C$lambda, Reff_C = sim_C$R_eff, N_C = sim_C$N_total,
      R0 = R0_draw
    )

    vax_bookkeeping[[i]] <- tibble(
      draw = d, t = seq_len(N), week_start = inputs$dates,
      entrants_total_B = sim_B$entrants_to_target_total, entrants_susceptible_B = sim_B$entrants_to_target_susceptible,
      doses_B = sim_B$doses_administered, effective_B = sim_B$effective_protected,
      doses_C = sim_C$doses_administered, effective_C = sim_C$effective_protected
    )

    qa_rows[[i]] <- tibble(
      draw = d,
      pop_conservation_max_err_A = max(abs(sim_A$N_total - inputs$N_start)),
      pop_conservation_max_err_B = max(abs(sim_B$N_total - inputs$N_start)),
      pop_conservation_max_err_C = max(abs(sim_C$N_total - inputs$N_start)),
      R0_identical_A_B = max(abs(R0_draw - R0_draw)), # trivially 0 -- same input vector; recorded for completeness
      max_diff_lambda_A_C_before_vax = NA_real_ # filled in validation script using per-cohort timing
    )

    if (i %% 25 == 0) message(sprintf("  ...%d / %d draws done", i, n_use))
  }

  saveRDS(list(agg_S = agg_S, agg_N = agg_N, agg_Inew = agg_Inew, age_group_labels = REPORTING_AGE_LABELS,
               totals = bind_rows(totals), vax_bookkeeping = bind_rows(vax_bookkeeping), qa = bind_rows(qa_rows),
               draw_idx = draw_idx, dates = inputs$dates, config = list(n_use = n_use, seed = RANDOM_SEED, window_start = WINDOW_START, window_end = WINDOW_END)),
          file.path(DIR_RESULTS, "three_arm_results.rds"))
  message("[saved] ", file.path(DIR_RESULTS, "three_arm_results.rds"))
  invisible(NULL)
}
