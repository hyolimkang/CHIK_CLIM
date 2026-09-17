# Mode audit: posterior decomposition (Section 2), bridge-sampling mode
# mass comparison (Section 3), epidemiological self-consistency of the
# high-q mode (Section 4), serology contrast (Section 5), and short-vs-full
# period comparison (Section 6). Does NOT modify the model, prior,
# kappa_sero, or sampler. Uses the two CLEAN mode-targeted fits:
# modeB_lowclean (3 chains, all confirmed q~0.05) and modeB_high (4 chains,
# all confirmed q~0.83).

required_packages <- c("here", "rstan", "bridgesampling", "dplyr", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(bridgesampling); library(dplyr); library(tibble); library(readr) })

audit_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

load_mode_fit <- function(root, tag) {
  path <- file.path(root, "03_Output/02_ceara_pipeline/model_fits/baseline/v4_3_short_2015_2019_q_calibration/outputs", tag, paste0("renewal_ceara_v4_3_fit_", tag, ".rds"))
  readRDS(path)
}

# Section 2: decompose lp__ into case log-lik, serology log-lik, eta_q prior,
# and all OTHER prior contributions -- computed post-hoc in R from the
# already-extracted parameter draws (no Stan model change).
decompose_posterior <- function(bundle, label) {
  fit <- bundle$fit
  alpha_R <- rstan::extract(fit, pars = "alpha_R")$alpha_R
  z_year <- rstan::extract(fit, pars = "z_year")$z_year
  beta_sin <- rstan::extract(fit, pars = "beta_sin")$beta_sin
  beta_cos <- rstan::extract(fit, pars = "beta_cos")$beta_cos
  phi_obs <- rstan::extract(fit, pars = "phi_obs")$phi_obs
  q <- rstan::extract(fit, pars = "q")$q
  log_lik_cases <- rstan::extract(fit, pars = "log_lik_cases")$log_lik_cases
  log_lik_serology <- rstan::extract(fit, pars = "log_lik_serology")$log_lik_serology
  log_prior_q <- rstan::extract(fit, pars = "log_prior_q")$log_prior_q
  lp <- rstan::extract(fit, pars = "lp__")$lp__

  log_prior_alpha_R <- dnorm(alpha_R, log(1.2), 0.5, log = TRUE)
  log_prior_z_year <- rowSums(dnorm(z_year, 0, 1, log = TRUE))
  log_prior_beta_sin <- dnorm(beta_sin, 0, 0.25, log = TRUE)
  log_prior_beta_cos <- dnorm(beta_cos, 0, 0.25, log = TRUE)
  log_prior_phi_obs <- dgamma(phi_obs, 2, 0.1, log = TRUE)
  case_ll_sum <- rowSums(log_lik_cases)
  other_priors_sum <- log_prior_alpha_R + log_prior_z_year + log_prior_beta_sin + log_prior_beta_cos + log_prior_phi_obs

  tibble(
    mode = label, draw = seq_along(as.numeric(lp)),
    case_loglik = as.numeric(case_ll_sum), serology_loglik = as.numeric(log_lik_serology),
    q_prior_logdensity = as.numeric(log_prior_q), other_priors_logdensity = as.numeric(other_priors_sum),
    lp__ = as.numeric(lp), q = as.numeric(q), alpha_R = as.numeric(alpha_R)
  )
}

run_mode_audit_analysis <- function() {
  root <- audit_root()
  base_dir <- file.path(root, "03_Output", "02_ceara_pipeline", "model_fits", "baseline", "v4_3_short_2015_2019_q_calibration")

  low_bundle <- load_mode_fit(root, "modeB_lowclean")
  high_bundle <- load_mode_fit(root, "modeB_high")

  # --- Section 2: posterior decomposition ------------------------------------
  decomp <- bind_rows(decompose_posterior(low_bundle, "low_q"), decompose_posterior(high_bundle, "high_q"))
  write_csv(decomp, file.path(base_dir, "mode_audit_posterior_decomposition.csv"))
  decomp_summary <- decomp |> group_by(mode) |> summarise(
    case_loglik_median = median(case_loglik), serology_loglik_median = median(serology_loglik),
    q_prior_logdensity_median = median(q_prior_logdensity), other_priors_median = median(other_priors_logdensity),
    lp___median = median(lp__), .groups = "drop"
  )
  write_csv(decomp_summary, file.path(base_dir, "mode_audit_decomposition_summary.csv"))
  message("[mode audit] posterior decomposition (medians):")
  print(as.data.frame(decomp_summary))
  diffs <- decomp_summary |> summarise(
    d_case_loglik = case_loglik_median[mode == "low_q"] - case_loglik_median[mode == "high_q"],
    d_serology_loglik = serology_loglik_median[mode == "low_q"] - serology_loglik_median[mode == "high_q"],
    d_q_prior = q_prior_logdensity_median[mode == "low_q"] - q_prior_logdensity_median[mode == "high_q"],
    d_lp = lp___median[mode == "low_q"] - lp___median[mode == "high_q"]
  )
  message("[mode audit] low_q MINUS high_q (positive = low_q favoured):")
  print(as.data.frame(diffs))

  # --- Section 3: bridge sampling (rigorous, since the package is available) --
  # A stanfit's compiled-model pointer does not survive saveRDS/readRDS across
  # R sessions (rstan's well-known "model object is not created or not valid"
  # log_prob error) -- bridge_sampler() needs a live log_prob, so reattach a
  # freshly compiled stanmodel (identical source, not a model change) before
  # calling it.
  message("[mode audit] recompiling stanmodel and reattaching a live log_prob instance for bridge sampling...")
  fresh_model <- rstan::stan_model(file.path(root, "02_Script", "00_shared", "stan", "current", "ceara", "renewal_ceara_v4_3_weak_q.stan"))
  reattach_live_log_prob <- function(bundle) {
    throwaway <- rstan::sampling(fresh_model, data = bundle$stan_data, chains = 1, iter = 2, warmup = 1, refresh = 0)
    bundle$fit@stanmodel <- fresh_model
    bundle$fit@.MISC <- throwaway@.MISC
    bundle$fit
  }
  low_bundle$fit <- reattach_live_log_prob(low_bundle)
  high_bundle$fit <- reattach_live_log_prob(high_bundle)
  message("[mode audit] running bridge sampling for each mode (this evaluates the log_prob many times, may take a moment)...")
  bridge_low <- bridgesampling::bridge_sampler(low_bundle$fit, silent = TRUE)
  bridge_high <- bridgesampling::bridge_sampler(high_bundle$fit, silent = TRUE)
  bf_low_over_high <- bridgesampling::bf(bridge_low, bridge_high)
  mode_mass <- tibble(
    mode = c("low_q", "high_q"), log_marginal_likelihood = c(bridge_low$logml, bridge_high$logml)
  ) |> mutate(relative_posterior_mass = exp(log_marginal_likelihood - max(log_marginal_likelihood)))
  mode_mass$relative_posterior_mass <- mode_mass$relative_posterior_mass / sum(mode_mass$relative_posterior_mass)
  write_csv(mode_mass, file.path(base_dir, "mode_audit_bridge_sampling.csv"))
  message("[mode audit] bridge sampling log marginal likelihoods:")
  print(as.data.frame(mode_mass))
  message(sprintf("[mode audit] Bayes factor (low_q vs high_q) = %.4g", as.numeric(bf_low_over_high$bf)))
  message("[mode audit] NOTE: this approximates the RELATIVE mass of the two identified local modes assuming no other modes exist; it is not a full-model marginal likelihood.")

  # --- Section 4: epidemiological self-consistency of each mode --------------
  eco_summary <- function(bundle, label) {
    fit <- bundle$fit
    years <- as.integer(format(as.Date(bundle$weekly_data$week_start), "%Y"))
    R0_t <- rstan::extract(fit, pars = "R0_t")$R0_t
    R_eff_t <- rstan::extract(fit, pars = "R_eff_t")$R_eff_t
    X <- rstan::extract(fit, pars = "X")$X
    R0_t_median <- apply(R0_t, 2, median)
    R_eff_t_median <- apply(R_eff_t, 2, median)
    growth_idx <- which(years %in% c(2016L, 2017L))
    tibble(
      mode = label,
      max_R0_t = max(R0_t_median), frac_weeks_R0_gt1 = mean(R0_t_median > 1),
      R0_growth_median = median(R0_t_median[growth_idx]), Reff_growth_median = median(R_eff_t_median[growth_idx]),
      frac_weeks_Reff_gt1 = mean(R_eff_t_median > 1)
    )
  }
  eco <- bind_rows(eco_summary(low_bundle, "low_q"), eco_summary(high_bundle, "high_q"))
  write_csv(eco, file.path(base_dir, "mode_audit_epidemiological_consistency.csv"))
  message("[mode audit] epidemiological self-consistency:")
  print(as.data.frame(eco))

  # Annual case PPC per mode (Section 4 last bullet).
  annual_ppc <- function(bundle, label) {
    fit <- bundle$fit
    years <- as.integer(format(as.Date(bundle$weekly_data$week_start), "%Y"))
    expected <- rstan::extract(fit, pars = "expected_reported_cases")$expected_reported_cases
    bind_rows(lapply(sort(unique(years)), function(yr) {
      idx <- which(years == yr)
      exp_annual <- rowSums(expected[, idx, drop = FALSE])
      tibble(mode = label, year = yr, observed = sum(bundle$weekly_data$cases[idx]),
             expected_median = median(exp_annual), expected_q025 = quantile(exp_annual, .025), expected_q975 = quantile(exp_annual, .975))
    }))
  }
  ppc <- bind_rows(annual_ppc(low_bundle, "low_q"), annual_ppc(high_bundle, "high_q")) |>
    mutate(ratio_expected_over_observed = expected_median / observed)
  write_csv(ppc, file.path(base_dir, "mode_audit_annual_case_ppc.csv"))
  message("[mode audit] annual case PPC by mode:")
  print(as.data.frame(ppc))

  # --- Section 5: serology contrast -------------------------------------------
  sero_contrast <- function(bundle, label) {
    fit <- bundle$fit
    p_sero <- rstan::extract(fit, pars = "p_state_sero_at_anchor")$p_state_sero_at_anchor
    sero_pred <- rstan::extract(fit, pars = "sero_pred")$sero_pred
    tibble(mode = label, immune_2018_median = median(p_sero), immune_2018_q025 = quantile(p_sero, .025), immune_2018_q975 = quantile(p_sero, .975),
           sero_pred_median = median(sero_pred), sero_pred_q025 = quantile(sero_pred, .025), sero_pred_q975 = quantile(sero_pred, .975),
           observed_juazeiro = 103 / 404)
  }
  sero <- bind_rows(sero_contrast(low_bundle, "low_q"), sero_contrast(high_bundle, "high_q"))
  write_csv(sero, file.path(base_dir, "mode_audit_serology_contrast.csv"))
  message("[mode audit] serology contrast:")
  print(as.data.frame(sero))

  # --- Section 6: short low-q mode vs full-period v4.2 ------------------------
  full_v4_2_path <- file.path(root, "02_Script/90_development_archive/02_historical_renewal_versions/21_v4_2_q_juazeiro_serology/outputs/stageB4/renewal_ceara_v4_2_fit_stageB4.rds")
  if (file.exists(full_v4_2_path)) {
    full_bundle <- readRDS(full_v4_2_path)
    q_full <- rstan::extract(full_bundle$fit, pars = "q")$q
    p_sero_full <- rstan::extract(full_bundle$fit, pars = "p_state_sero_at_anchor")$p_state_sero_at_anchor
    q_low <- rstan::extract(low_bundle$fit, pars = "q")$q
    p_sero_low <- rstan::extract(low_bundle$fit, pars = "p_state_sero_at_anchor")$p_state_sero_at_anchor
    comparison <- tibble(
      period = c("2015-2019 (low-q mode)", "2015-2025 (v4.2 Stage B4)"),
      q_median = c(median(q_low), median(q_full)), q_q025 = c(quantile(q_low, .025), quantile(q_full, .025)), q_q975 = c(quantile(q_low, .975), quantile(q_full, .975)),
      immune_2018_median = c(median(p_sero_low), median(p_sero_full))
    )
    write_csv(comparison, file.path(base_dir, "mode_audit_short_vs_full_period.csv"))
    message("[mode audit] short (low-q mode) vs full-period v4.2:")
    print(as.data.frame(comparison))
  } else {
    message("[mode audit] full-period v4.2 Stage B4 fit not found -- skipping Section 6 comparison.")
  }

  invisible(list(decomp = decomp, mode_mass = mode_mass, eco = eco, ppc = ppc, sero = sero))
}

if (sys.nframe() == 0L) run_mode_audit_analysis()
