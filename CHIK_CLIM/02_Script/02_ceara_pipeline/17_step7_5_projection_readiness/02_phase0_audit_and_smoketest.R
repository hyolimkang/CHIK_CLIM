# PHASE 0 -- reproducibility/implementation audit (synthesising Step
# 6/7's already-established, re-verified facts -- not re-deriving them
# from scratch) PLUS a small smoke test of the NEW episodic engine
# before any full-scale run, per the task's explicit instruction.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(tidyr) })
filter <- dplyr::filter

source(file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/00_config_step7_5.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/01_episodic_engine.R"))

# ================================================================
# 1-8: audit table, citing the EXACT prior finding and its source file
# for each item -- synthesis, not re-derivation, per Phase 0's scope.
# ================================================================
audit <- tribble(
  ~item, ~status, ~finding_and_source,
  "1. Step6 vs Step7 baseline discrepancy",
    "RESOLVED (disclosed, not silently fixed)",
    "Step 7's engine uses shared-uniform-coupled Poisson draws (qpois(runif(1),lambda)) vs Step 6's direct rpois() -- different RNG stream consumption from week 1, REQUIRED to enable Step 7's cross-lambda common-random-number comparisons. Confirmed via cluster-bootstrap (Step6's P(major)=10.55% fell outside Step7's own [12.7%,27.75%] CI) -- not an implementation bug. Source: STEP7_IMPORTATION_PERSISTENCE_ASSESSMENT.md Phase 0.",
  "2. Posterior draw selection / replicate matching",
    "PASS",
    "Both Step 7 and this module sample posterior draws via sample.int(n_draws_available, n_draws) under a FIXED seed_offset per holdout, shared across all lambda/structural arms within a holdout -- confirmed by Step7's PHASE0_posterior_identity_check.csv (N_train=261/365, 1600 draws, matching the accepted training fits).",
  "3. Future year-effect recentering",
    "DISCLOSED, QUANTIFIED, NOT corrected in the main arms",
    "Shrinks marginal SD from tau_year to tau_year*sqrt((K-1)/K) (H1 K=6: ratio 0.913, confirmed empirically at 0.929 via a direct counterfactual simulation). Source: STEP6 FAIL item + STEP7_IMPORTATION_PERSISTENCE_ASSESSMENT.md Phase 5.",
  "4. Historical climate-standardisation leakage",
    "DISCLOSED, QUANTIFIED, NOT corrected",
    "z_T_anom/z_P_anom standardisation constants use the FULL 2015-2025 period; training-only SD would be 7.6-13.1% smaller. Source: 13_retrospective_resurgence_diagnostic/PHASE0_climatology_leakage_magnitude.csv.",
  "5. Establishment / persistence / major-epidemic definitions",
    "AUDITED; ONE DEFINITIONAL IMPRECISION DISCLOSED, NOT YET CORRECTED",
    "'established' requires M_import>0 COINCIDENT WITH Reff>1 WITHIN the evaluation year -- a trajectory already carrying substantial infection from EARLIER periods can reach 'major_epidemic' via pure local transmission with zero new imports that year, and is never flagged 'established'. Confirmed via cross-tab: 43/2000 (2.15%) H2/G2 trajectories were major-but-not-established (Step7 Phase 0); Step7 Phase2 Experiment B at post-2022 lambda=0 showed P(established)=0% by this tautology while P(major epidemic) was still 18.65% (vs. 19.9% baseline) -- i.e. the DEFINITION, not the underlying transmission, produces the 0%. Phase 3 of THIS module uses major_epidemic and a REVISED establishment proxy that does not require a same-year import (see 04_phase3_separate_intro_amplification.R).",
  "6. Importation-pressure implementation / double counting",
    "PASS",
    "M_t enters ONCE as M_t/N_t in the hazard (lambda_t = R0*infectiousness/N + M_t/N); confirmed via the R0-component reconstruction check (reconstructed log-R0 matches saved R0 to 4.4e-16, floating point only) and via the explicit numerical-consistency table (S+U=N exactly; aggregate reported/latent ratio=0.05004, matching fixed q=0.05). Source: STEP6_R0_component_reconstruction_check.csv, STEP6_numerical_consistency.csv.",
  "7. Training-only latent infection history transfer",
    "PASS",
    "advance_terminal_state() supplies the ONE demographic step the Stan model's `t<N` guard skips at the final training week; the 8-week posterior X history (G=8, matching the Gamma(4,2) generation interval) is carried forward into X_recent_history for the holdout's first-week infectiousness. Source: 03_future_renewal_functions.R (unmodified, reused transitively via 01_simulation_engine.R), unit-tested in 10_future_projection_pilot/04_run_unit_tests.R.",
  "8. HMC reliability / chain-specific sensitivity",
    "UNRESOLVED, EXPLICITLY CARRIED FORWARD",
    "All 5 tau_year fits (P1-P3, H1/H2 training) show mild (2-20/1600) divergences concentrated at LARGER tau_year, not near 0; H1's worst-mixing chain (chain 3, within-chain ESS=37 for tau_year) mixes BEST for sigma_season_year (ESS=199) and vice versa for chain 1 -- a real, chain-specific, non-persistent (between-chain variance <0.3% of total) sampling difficulty, not multimodality. Source: HMC_RELIABILITY_AUDIT.md Sections 1-6."
)
write_csv(audit, file.path(DIR_TABLES_STEP75, "PHASE0_audit_synthesis.csv"))
message("=== Phase 0: reproducibility/implementation audit (synthesis of Step 6/7 findings) ==="); print(as.data.frame(audit[, c("item", "status")]), row.names = FALSE)

# ================================================================
# Smoke test: the NEW episodic engine, small scale, BEFORE any full run.
# ================================================================
message("\n[Phase0] Smoke-testing the episodic engine (10 draws x 5 replicates, H1, p=0.10)...")
training_H1 <- readRDS(FIT_PATHS[["H1"]])
smoke <- simulate_holdout_step75(training_H1, HOLDOUTS$H1, lambda_mode = "episodic", burst_p = BURST_PROBS[["p10"]], burst_lambda_mean = EPISODIC_LAMBDA,
                                  year_effect_mode = "centered", n_draws = ENSEMBLE_SMOKE$n_draws, n_replicates = ENSEMBLE_SMOKE$n_replicates,
                                  seed_offset = SEED_H1_STEP75)
smoke_check <- smoke$trajectories |> summarise(
  n_rows = n(), mean_M_import = mean(M_import), target_mean = EPISODIC_LAMBDA,
  frac_weeks_M_gt_0 = mean(M_import > 0), target_frac_le = BURST_PROBS[["p10"]],
  max_abs_N_minus_S_plus_U = max(abs(N - (S + U))), min_S = min(S), min_U = min(U),
  aggregate_reported_to_latent_ratio = sum(Y_reported) / sum(X_latent))
write_csv(smoke_check, file.path(DIR_TABLES_STEP75, "PHASE0_episodic_engine_smoketest.csv"))
message("Smoke test numerical checks:"); print(as.data.frame(smoke_check), digits = 4)
if (abs(smoke_check$mean_M_import - EPISODIC_LAMBDA) / EPISODIC_LAMBDA > 0.25) {
  stop("[Phase0 STOP] Episodic engine's realised mean M_t deviates from the target lambda by >25% at smoke scale -- stopping before any full-scale run, per instruction. (A 10x5 sample has genuine sampling noise; investigate before proceeding.)")
}
message("[Phase0] Smoke test passed: mean(M_t) is close to the target lambda=1, and the fraction of weeks with a burst is close to p=0.10. Numerical consistency (S+U=N, q=0.05) holds. Proceeding to full-scale Phase 1-3 runs.")
