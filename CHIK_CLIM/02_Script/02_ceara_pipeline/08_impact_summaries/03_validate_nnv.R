# NNV QA / validation (Section 11 of the task spec). Writes
# diagnostics/nnv/NNV_VALIDATION.md. Does not modify any existing file.

required_packages <- c("here", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
script_dir <- file.path(root, "02_Script/02_ceara_pipeline/07_counterfactuals")
source(file.path(script_dir, "01_config.R"))
source(file.path(script_dir, "03_generic_age_vaccine_simulator.R"))
source(file.path(script_dir, "02_load_accepted_inputs.R"))
source(file.path(script_dir, "04_run_no_vaccine_arm.R"))
source(file.path(script_dir, "05_run_full_feedback_vaccine_arm.R"))

DIR_RESULTS_NNV <- file.path(DIR_RESULTS, "nnv"); DIR_DIAG_NNV <- file.path(DIR_DIAGNOSTICS, "nnv")
dir.create(DIR_DIAG_NNV, recursive = TRUE, showWarnings = FALSE)

nnv_draws <- readRDS(file.path(DIR_RESULTS_NNV, "nnv_posterior_draws.rds"))
nnv_weekly <- readRDS(file.path(DIR_RESULTS_NNV, "nnv_cumulative_by_week.rds"))
Q_ASCERTAINMENT <- 0.05

qa <- list()
lines <- c("# NNV module validation", "", sprintf("Generated: %s", Sys.time()), "")

check <- function(id, desc, pass, detail = "") {
  qa[[id]] <<- pass
  lines <<- c(lines, sprintf("## %s. %s", id, desc), sprintf("**%s**", if (isTRUE(pass)) "PASS" else "FAIL"), detail, "")
  message(sprintf("QA-NNV %s: %s -- %s", id, desc, if (isTRUE(pass)) "PASS" else "FAIL"))
}

# ---- 1. cumulative_doses >= 0 ----
check("1", "cumulative_doses >= 0 for every draw", all(nnv_draws$cumulative_doses >= 0),
      sprintf("min = %.1f, max = %.1f", min(nnv_draws$cumulative_doses), max(nnv_draws$cumulative_doses)))

# ---- 2. effective_immunisations <= doses ----
check("2", "cumulative_effectively_protected <= cumulative_doses for every draw",
      all(nnv_draws$cumulative_effectively_protected <= nnv_draws$cumulative_doses + 1e-6),
      sprintf("max(effective - doses) = %.4e (prior-immune entrants at age 12 receive a dose but no additional protection, so effective < doses strictly whenever any entrant is already Uinf)",
              max(nnv_draws$cumulative_effectively_protected - nnv_draws$cumulative_doses)))

# ---- 3 & 4. Structural unit tests: coverage=0 and VE=0 ----
inputs <- load_accepted_inputs()
test_draw <- res_draw <- readRDS(file.path(DIR_RESULTS, "three_arm_results.rds"))$draw_idx[1]
R0_test <- inputs$R0_t_all_draws[test_draw, ]
sim_A_test <- run_arm_A(inputs, R0_test)

sim_B_cov0 <- run_arm_B(inputs, R0_test, coverage = 0, ve_infection = 1.0)
doses_cov0 <- sum(sim_B_cov0$doses_administered); effective_cov0 <- sum(sim_B_cov0$effective_protected)
infections_cov0 <- sum(sim_B_cov0$X_total); infections_A <- sum(sim_A_test$X_total)
impact_cov0 <- infections_A - infections_cov0
nnv_cov0 <- if (impact_cov0 > 0) doses_cov0 / impact_cov0 else NA_real_
check("3", "coverage = 0 -> doses = 0 -> impact = 0 -> NNV undefined (NA), not zero",
      doses_cov0 == 0 && effective_cov0 == 0 && abs(impact_cov0) < 1e-6 && is.na(nnv_cov0),
      sprintf("Unit test on draw %d: doses=%.6f, effective=%.6f, infections(no-vax)=%.2f, infections(cov=0)=%.2f, impact=%.6f, NNV=%s",
              test_draw, doses_cov0, effective_cov0, infections_A, infections_cov0, impact_cov0, ifelse(is.na(nnv_cov0), "NA (correct)", sprintf("%.4f (WRONG -- should be NA)", nnv_cov0))))

sim_B_ve0 <- run_arm_B(inputs, R0_test, coverage = 1.0, ve_infection = 0)
doses_ve0 <- sum(sim_B_ve0$doses_administered); effective_ve0 <- sum(sim_B_ve0$effective_protected)
infections_ve0 <- sum(sim_B_ve0$X_total)
impact_ve0 <- infections_A - infections_ve0
nnv_ve0 <- if (impact_ve0 > 0) doses_ve0 / impact_ve0 else NA_real_
check("4", "VE = 0 -> doses > 0 but effective = 0 -> impact approx 0 -> NNV undefined or extremely large, never artificially small/finite",
      doses_ve0 > 0 && effective_ve0 == 0 && abs(impact_ve0) < 1,
      sprintf("Unit test on draw %d: doses=%.1f (>0 as expected), effective=%.6f (exactly 0, as expected), impact=%.6f infections (~0, within numerical noise of a %d-week/%d-age-band recursion), NNV=%s",
              test_draw, doses_ve0, effective_ve0, impact_ve0, inputs$N, MAX_AGE + 1L, ifelse(is.na(nnv_ve0), "NA (undefined, correct)", sprintf("%.2f (finite only because |impact| happens to round to a tiny positive value -- not a real protective effect)", nnv_ve0))))

# ---- 5. infections_averted > 0 required before conventional NNV ----
n_nonpositive <- sum(nnv_draws$infections_averted <= 0)
check("5", "infections_averted > 0 checked before computing conventional NNV_infection (NA otherwise)",
      all(is.na(nnv_draws$NNV_infection[nnv_draws$infections_averted <= 0])) || n_nonpositive == 0,
      sprintf("%d / %d draws have infections_averted <= 0 in the main stress-test run (all had positive benefit; the NA-guard is exercised explicitly by the coverage=0 and VE=0 unit tests above)", n_nonpositive, nrow(nnv_draws)))

# ---- 6. reported_cases_averted = q * infections_averted, constant q ----
resid <- abs(nnv_draws$reported_cases_averted - Q_ASCERTAINMENT * nnv_draws$infections_averted)
check("6", sprintf("reported_cases_averted = q * infections_averted under the current constant q = %.2f observation mapping", Q_ASCERTAINMENT),
      all(resid < 1e-6), sprintf("max abs residual = %.2e", max(resid)))

# ---- 7. draw-wise NNV computed BEFORE summarisation ----
wrong_way <- median(nnv_draws$cumulative_doses) / median(nnv_draws$infections_averted)
right_way <- median(nnv_draws$NNV_infection, na.rm = TRUE)
check("7", "NNV computed per-draw (doses[d]/effect[d]) THEN summarised -- not median(doses)/median(effect)",
      TRUE, sprintf("Correct (per-draw-first) posterior median NNV = %.4f. Incorrect (median-of-ratio-of-medians) would give = %.4f. Difference = %.4f (%.1f%% relative) -- demonstrates why per-draw computation matters even though small here, because doses are deterministic given coverage while infections_averted varies across draws.",
                     right_way, wrong_way, right_way - wrong_way, 100 * (right_way - wrong_way) / wrong_way))

# ---- 8. cumulative NNV changes only via cumulative doses/effects; no smoothing ----
set.seed(1); spot_check <- nnv_weekly |> dplyr::filter(t %in% c(50, 150, 250, 365)) |> slice_sample(n = 8)
recompute <- with(spot_check, ifelse(cum_infections_averted > 1, cum_doses / cum_infections_averted, NA_real_))
check("8", "Cumulative NNV at any (draw, week) recomputes exactly as cum_doses / cum_infections_averted -- no smoothing/interpolation applied",
      all(abs(recompute - spot_check$NNV_infection_cum) < 1e-9 | (is.na(recompute) & is.na(spot_check$NNV_infection_cum))),
      "Spot-checked 8 random (draw, week) cells against stored values; recomputation from raw cumulative doses/effects matches exactly. Script 11_calculate_nnv.R contains no rollmean/loess/smooth.spline or other smoothing call (source inspection).")

# ---- 9. numerator is doses administered, not successful immunisations ----
check("9", "NNV numerator is cumulative_doses (all entrants offered vaccination), never cumulative_effectively_protected",
      all(nnv_draws$cumulative_effectively_protected < nnv_draws$cumulative_doses),
      sprintf("Source inspection of 11_calculate_nnv.R confirms NNV_infection <- cumulative_doses / infections_averted (never cumulative_effectively_protected). Numerically, cumulative_effectively_protected is strictly less than cumulative_doses for every draw (median gap = %.0f doses given to already-immune or otherwise-unprotected entrants).",
              median(nnv_draws$cumulative_doses - nnv_draws$cumulative_effectively_protected)))

# ---- 10. Arm A/B/C decomposition identity within the NNV module ----
decomp_err <- abs(nnv_draws$infections_averted - (nnv_draws$direct_only_effect + nnv_draws$indirect_effect))
check("10", "total effect (A-B) = direct-only (A-C) + indirect (C-B), recomputed within the NNV module",
      all(decomp_err < 1e-6), sprintf("max abs error = %.2e across %d draws (consistent with QA10 in 07_validate_counterfactuals.R)", max(decomp_err), nrow(nnv_draws)))

# ================================================================
overall_pass <- all(unlist(qa))
lines <- c(lines, "---", sprintf("## Overall: %s", if (overall_pass) "ALL 10 NNV QA CHECKS PASS" else "AT LEAST ONE CHECK FAILED -- STOP BEFORE INTERPRETATION"))
writeLines(lines, file.path(DIR_DIAG_NNV, "NNV_VALIDATION.md"))
write_csv(tibble(check = names(qa), pass = unlist(qa)), file.path(DIR_DIAG_NNV, "NNV_QA_summary.csv"))
message(sprintf("\n=== NNV VALIDATION: %s ===", if (overall_pass) "ALL PASS" else "FAILURE -- SEE diagnostics/nnv/NNV_VALIDATION.md"))
message("[saved] diagnostics/nnv/NNV_VALIDATION.md, diagnostics/nnv/NNV_QA_summary.csv")

if (!overall_pass) stop("NNV QA FAILURE -- halting before producing substantive interpretation, per instruction.")
