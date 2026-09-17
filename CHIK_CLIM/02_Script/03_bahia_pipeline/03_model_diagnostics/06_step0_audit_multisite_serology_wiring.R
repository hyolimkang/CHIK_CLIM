# Step 0 audit (per explicit instruction) for the Bahia multisite-serology
# q=0.05 pilot, BEFORE any global-q model is built: confirm the serology
# likelihood is actually wired into `target`, report full data-side audit
# (J_sero, IDs, n_positive/n_tested, window indices, sero_geographic_sd,
# kappa_sero), report posterior summaries of log_lik_serology / eta_geo /
# p_state_window / p_site_window, and report the FULL, itemised reason the
# q=0.05 pilot was labelled HMC FAIL/marginal (max Rhat, divergences,
# max-treedepth hits, min bulk ESS, min tail ESS, BFMI by chain).

required_packages <- c("rstan", "dplyr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/03_bahia_pipeline/model_fits/v4_9_replication")

audit_one_q <- function(qval) {
  tag <- paste0("q", sprintf("%.2f", qval))
  b <- readRDS(file.path(base_dir, "outputs_multisite_serology", tag, paste0("renewal_bahia_multisite_fit_", tag, ".rds")))
  fit <- b$fit
  sd <- b$stan_data
  audit <- b$sero_audit

  message("\n================ Step 0 audit: ", tag, " ================")

  message("\n--- Data-side wiring (as passed to Stan) ---")
  wiring <- tibble(
    sero_id = audit$sero_id, location = audit$location,
    n_positive = sd$sero_n_positive, n_tested = sd$sero_n_tested,
    window_start_idx = sd$sero_window_start_idx, window_n_weeks = sd$sero_window_n_weeks
  )
  print(as.data.frame(wiring))
  message("J_sero = ", sd$J_sero, " | sero_geographic_sd (fixed) = ", sd$sero_geographic_sd,
          " | kappa_sero (fixed) = ", sd$kappa_sero)

  message("\n--- Confirm serology sampling statement contributes to target ---")
  model_code <- rstan::get_stancode(fit)
  model_block <- sub(".*model \\{", "", model_code)
  model_block <- sub("\\}\\s*generated quantities.*", "", model_block)
  has_sero_sampling <- grepl("sero_n_positive\\[j\\] ~ beta_binomial", model_block)
  message("Serology sampling statement found inside model block: ", has_sero_sampling)
  has_log_lik_gq <- "log_lik_serology" %in% names(rstan::extract(fit, pars = "log_lik_serology"))
  message("log_lik_serology[] present as generated quantity: ", has_log_lik_gq)

  message("\n--- Posterior: log-likelihood contribution of serology ---")
  log_lik_serology <- rstan::extract(fit, "log_lik_serology")$log_lik_serology
  sum_ll <- rowSums(log_lik_serology)
  message(sprintf("sum(log_lik_serology): median=%.2f [%.2f, %.2f]",
                   median(sum_ll), quantile(sum_ll, .025), quantile(sum_ll, .975)))
  ll_table <- tibble(
    sero_id = audit$sero_id,
    log_lik_median = apply(log_lik_serology, 2, median),
    log_lik_lo = apply(log_lik_serology, 2, quantile, .025),
    log_lik_hi = apply(log_lik_serology, 2, quantile, .975)
  )
  print(as.data.frame(ll_table))

  message("\n--- Posterior: eta_geo / p_state_window / p_site_window ---")
  eta_geo <- rstan::extract(fit, "eta_geo")$eta_geo
  p_state_window <- rstan::extract(fit, "p_state_window")$p_state_window
  p_site_window <- rstan::extract(fit, "p_site_window")$p_site_window
  summary_table <- tibble(
    sero_id = audit$sero_id,
    observed_prevalence = audit$observed_prevalence,
    eta_geo_median = apply(eta_geo, 2, median), eta_geo_lo = apply(eta_geo, 2, quantile, .025), eta_geo_hi = apply(eta_geo, 2, quantile, .975),
    p_state_window_median = apply(p_state_window, 2, median),
    p_site_window_median = apply(p_site_window, 2, median)
  )
  print(as.data.frame(summary_table))

  message("\n--- Exact reason for HMC FAIL/marginal label ---")
  h <- b$hmc
  message(sprintf("divergences=%d | max_treedepth_hits=%d | max_rhat=%.4f | min_bulk_ess=%.1f | hmc_pass=%s",
                   h$divergences, h$max_treedepth_hits, h$maximum_rhat, h$minimum_bulk_ess, h$hmc_pass))
  message("BFMI by chain: ", paste(sprintf("chain%d=%.3f", seq_along(h$bfmi), h$bfmi), collapse = " | "))

  message("\n--- Min tail ESS (not part of the original gate function; computed here for completeness) ---")
  pars_gate <- c("alpha_R", "z_year", "beta_sin1", "beta_cos1", "beta_sin2", "beta_cos2",
                  "sigma_season_year", "z_season_year", "phi_obs", "z_geo")
  arr <- as.array(fit, pars = pars_gate)
  tail_ess <- apply(arr, 3, function(x) posterior::ess_tail(x))
  bulk_ess_check <- apply(arr, 3, function(x) posterior::ess_bulk(x))
  message(sprintf("min tail ESS = %.1f (parameter: %s)", min(tail_ess), names(which.min(tail_ess))))
  message(sprintf("min bulk ESS (posterior pkg, cross-check vs rstan summary) = %.1f (parameter: %s)", min(bulk_ess_check), names(which.min(bulk_ess_check))))

  which_rhat_max <- {
    smry <- rstan::summary(fit, pars = pars_gate)$summary
    rownames(smry)[which.max(smry[, "Rhat"])]
  }
  message("Parameter with maximum Rhat: ", which_rhat_max)

  invisible(list(wiring = wiring, ll_table = ll_table, summary_table = summary_table, hmc = h,
                 min_tail_ess = min(tail_ess), max_rhat_param = which_rhat_max))
}

if (!requireNamespace("posterior", quietly = TRUE)) {
  message("[note] 'posterior' package not installed -- installing not attempted automatically; tail-ESS section will fail.")
}

results <- lapply(c(0.05, 0.10), audit_one_q)
