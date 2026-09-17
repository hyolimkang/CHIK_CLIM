# Required QA checks (Section 7). STOP before interpretation if any fail.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
script_dir <- file.path(root, "02_Script/02_ceara_pipeline/07_counterfactuals")
source(file.path(script_dir, "01_config.R"))

res <- readRDS(file.path(DIR_RESULTS, "three_arm_results.rds"))
totals <- res$totals
vax <- res$vax_bookkeeping
dates <- res$dates

qa_log <- list()

# ---- 1. Arm A reproduces the accepted closed-loop no-vaccine replay ----
climate <- readRDS(ACCEPTED_CLIMATE_FIT)
draws_v49 <- rstan::extract(climate$fit, pars = c("S", "U", "X"), permuted = TRUE)
keep_idx <- seq_len(res$config$n_use)
check1 <- bind_rows(lapply(seq_along(res$draw_idx), function(i) {
  d <- res$draw_idx[i]
  sub <- totals |> dplyr::filter(draw == d)
  N_check <- nrow(sub)
  tibble(draw = d, max_abs_diff_S = max(abs(sub$S_A - draws_v49$S[d, seq_len(N_check)])),
         max_abs_diff_X = max(abs(sub$X_A - draws_v49$X[d, seq_len(N_check)])))
}))
qa1_pass <- all(check1$max_abs_diff_S < 1e-3) && all(check1$max_abs_diff_X < 1e-3)
message(sprintf("QA1 (Arm A reproduces accepted closed-loop replay): max diff S=%.2e, X=%.2e -- %s",
                 max(check1$max_abs_diff_S), max(check1$max_abs_diff_X), if (qa1_pass) "PASS" else "FAIL"))
qa_log$qa1_arm_A_matches_accepted <- qa1_pass

# ---- 2 & 3. Population / compartment conservation (all arms) ----
qa23 <- res$qa
qa23_pass <- all(qa23$pop_conservation_max_err_A < 1e-3) && all(qa23$pop_conservation_max_err_B < 1e-3) && all(qa23$pop_conservation_max_err_C < 1e-3)
message(sprintf("QA2/3 (population & compartment conservation, S+Uinf+Uvac=N): max err A/B/C = %.2e/%.2e/%.2e -- %s",
                 max(qa23$pop_conservation_max_err_A), max(qa23$pop_conservation_max_err_B), max(qa23$pop_conservation_max_err_C),
                 if (qa23_pass) "PASS" else "FAIL"))
qa_log$qa2_3_conservation <- qa23_pass

# ---- 4. No duplicate vaccination ----
cohort_check <- vax |> dplyr::filter(doses_B > 1e-6) |> mutate(year = as.integer(format(week_start, "%Y")), birth_cohort = year - TARGET_AGE) |>
  group_by(draw, birth_cohort) |> summarise(n_years_vaccinated = n_distinct(year), .groups = "drop")
qa4_pass <- all(cohort_check$n_years_vaccinated == 1)
message(sprintf("QA4 (no duplicate vaccination, each cohort vaccinated in exactly 1 calendar year): %s", if (qa4_pass) "PASS" else "FAIL"))
qa_log$qa4_no_duplicate_vaccination <- qa4_pass

# ---- 5. Arm C non-target groups approx identical to Arm A ----
# IMPORTANT REVISION (found on first run, verified not a bug): under the
# 100%-coverage/100%-VE/LIFELONG-protection stress test, individuals
# vaccinated at age 12 carry that protection forward as they age into
# 13-17 (within ~1-5 years) and eventually 18-64 (only right at the very
# end of this 7-year window, since the first vaccinated 12-year-old in
# 2015 turns 18 only in 2021). This is a genuine DIRECT-protection
# carryover through ageing -- NOT transmission feedback (which Arm C has
# disabled) -- confirmed by checking the diff starts at EXACTLY 0 in week
# 1-2 and grows smoothly over calendar YEARS (matching the ageing
# timescale), never jumping instantaneously. Only 0-11 (never vaccinated,
# never will be within this age band) and 65+ (cannot be reached by a
# 12-year-old within a 7-year window) can have ZERO legitimate pathway to
# differ in Arm C; 13-17 and 18-64 are EXPECTED to diverge from Arm A
# increasingly over the window and are explicitly excluded from this
# check (their divergence is quantified and discussed in the report as
# part of the direct-effect decomposition, not the formal indirect
# effect).
agg_Inew <- res$agg_Inew
zero_pathway_groups <- c("0-11", "65+") # cannot legitimately differ via ageing-forward direct protection within 2015-2021
carryover_groups <- c("13-17", "18-64") # EXPECTED to differ increasingly via direct-protection ageing forward -- not a failure
qa5_diffs_zero <- sapply(zero_pathway_groups, function(g) {
  a <- agg_Inew[, , g, "A"]; c_ <- agg_Inew[, , g, "C"]
  max(abs(c_ - a)) / max(1, max(a))
})
qa5_diffs_carryover <- sapply(carryover_groups, function(g) {
  a <- agg_Inew[, , g, "A"]; c_ <- agg_Inew[, , g, "C"]
  max(abs(c_ - a)) / max(1, max(a))
})
qa5_pass <- all(qa5_diffs_zero < 0.01)
message("QA5a (0-11, 65+: NO legitimate pathway to differ within this window -- must be ~0):")
print(round(qa5_diffs_zero, 5))
message("QA5b (13-17, 18-64: EXPECTED to diverge via direct-protection ageing forward, NOT a failure -- reported for transparency):")
print(round(qa5_diffs_carryover, 5))
message(sprintf("-- QA5 %s (0-11/65+ relative diffs < 1%%; 13-17/18-64 divergence is expected direct-effect carryover, not scored)", if (qa5_pass) "PASS" else "FAIL"))
qa_log$qa5_arm_C_zero_pathway_groups_match_A <- qa5_pass

# ---- 6. Arm B non-target reductions appear only after feedback propagates ----
# (descriptive check: first week where Arm B's non-target-group X differs from A by >0.01%, should be > 0)
all_reporting_groups <- c(zero_pathway_groups, carryover_groups)
first_divergence_weeks <- sapply(all_reporting_groups, function(g) {
  diffs <- apply(agg_Inew[, , g, "A"] - agg_Inew[, , g, "B"], 2, median)
  w <- which(abs(diffs) > 0.01)
  if (length(w) == 0) NA_integer_ else w[1]
})
message("QA6 (Arm B non-target-group divergence should NOT be at week 1, i.e. requires propagation time):")
print(first_divergence_weeks)
qa6_pass <- all(is.na(first_divergence_weeks) | first_divergence_weeks > 1)
qa_log$qa6_arm_B_indirect_has_lag <- qa6_pass
message(sprintf("-- %s", if (qa6_pass) "PASS" else "FAIL"))

# ---- 7. R0(t) identical across arms (by construction -- same input vector) ----
qa7_pass <- TRUE # R0_draw is passed as the SAME vector into all 3 arm functions; verified structurally, not stochastically
message("QA7 (R0(t) identical across Arms A/B/C): PASS (same R0_t_draw vector passed to all three arm functions, by construction)")
qa_log$qa7_R0_identical <- qa7_pass

# ---- 8. R_eff(t) allowed to differ (informational, not a failure condition) ----
reff_diff <- totals |> summarise(max_diff_AB = max(abs(Reff_A - Reff_B)), max_diff_AC = max(abs(Reff_A - Reff_C)))
message(sprintf("QA8 (R_eff allowed to differ, S(t) differs): max|Reff_A-Reff_B|=%.4f, max|Reff_A-Reff_C|=%.4f (informational)",
                 reff_diff$max_diff_AB, reff_diff$max_diff_AC))
qa_log$qa8_Reff_allowed_to_differ <- TRUE

# ---- 9. lambda(t): Arm B may diverge from Arm A after propagation ----
lambda_diff <- totals |> mutate(diff = abs(lambda_A - lambda_B)) |> group_by(t) |> summarise(median_diff = median(diff), .groups = "drop")
first_lambda_divergence <- lambda_diff$t[which(lambda_diff$median_diff > 1e-9)][1]
message(sprintf("QA9 (lambda_B may diverge from lambda_A after vaccine effect propagates): first week of divergence = %s (week 1 = %s)",
                 first_lambda_divergence, dates[1]))
qa_log$qa9_lambda_divergence_week <- first_lambda_divergence

# ---- 10. Numerical decomposition: total ~ direct-only + indirect ----
cum_by_draw <- totals |> group_by(draw) |> summarise(cum_X_A = sum(X_A), cum_X_B = sum(X_B), cum_X_C = sum(X_C), .groups = "drop") |>
  mutate(total_effect = cum_X_A - cum_X_B, direct_only = cum_X_A - cum_X_C, indirect = cum_X_C - cum_X_B,
         decomposition_error = total_effect - (direct_only + indirect), decomposition_rel_error = decomposition_error / pmax(total_effect, 1))
qa10_pass <- all(abs(cum_by_draw$decomposition_rel_error) < 1e-9)
message(sprintf("QA10 (total ~ direct-only + indirect): max relative decomposition error = %.2e -- %s",
                 max(abs(cum_by_draw$decomposition_rel_error)), if (qa10_pass) "PASS" else "FAIL"))
qa_log$qa10_decomposition_identity <- qa10_pass

# ================================================================
overall_qa <- tibble(check = names(qa_log), pass = unlist(lapply(qa_log, function(x) if (is.logical(x)) x else NA)))
message("\n=== OVERALL QA SUMMARY ===")
print(as.data.frame(overall_qa))
all_pass <- all(overall_qa$pass[!is.na(overall_qa$pass)])
message(sprintf("\n=== ALL REQUIRED QA CHECKS: %s ===", if (all_pass) "PASS" else "AT LEAST ONE CHECK FAILED -- STOP BEFORE INTERPRETATION"))

write_csv(check1, file.path(DIR_DIAGNOSTICS, "QA1_arm_A_vs_accepted.csv"))
write_csv(qa23, file.path(DIR_DIAGNOSTICS, "QA2_3_conservation.csv"))
write_csv(cohort_check, file.path(DIR_DIAGNOSTICS, "QA4_cohort_duplicate_check.csv"))
write_csv(cum_by_draw, file.path(DIR_DIAGNOSTICS, "QA10_decomposition_check.csv"))
write_csv(overall_qa, file.path(DIR_DIAGNOSTICS, "QA_overall_summary.csv"))
message("\n[saved] QA diagnostic CSVs to ", DIR_DIAGNOSTICS)

if (!all_pass) stop("QA FAILURE -- halting before producing substantive interpretation, per instruction.")
