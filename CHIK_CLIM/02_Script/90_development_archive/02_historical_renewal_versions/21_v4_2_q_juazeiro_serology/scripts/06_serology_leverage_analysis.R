# Serology leverage/conflict analysis: assembles the comparison table and
# leverage curve across the fixed kappa_sero ladder (50/100/200/500) plus
# the binomial limiting stress test (Section 2-5 of the task spec).
# Does NOT select a kappa or propose a production model -- reports only.

required_packages <- c("here", "rstan", "dplyr", "tibble", "readr", "ggplot2", "patchwork", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2); library(patchwork); library(tidyr) })

lev_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

summarise_one <- function(root, tag, kappa_numeric) {
  base <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "21_v4_2_q_juazeiro_serology")
  if (identical(tag, "kappa50")) {
    # kappa=50 IS the original Stage B4 fit -- reused directly, not refit (Section 2).
    out_dir <- file.path(base, "outputs", "stageB4")
    fit_path <- file.path(out_dir, "renewal_ceara_v4_2_fit_stageB4.rds")
  } else {
    out_dir <- file.path(base, "outputs", paste0("stageB4_", tag))
    fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_2_fit_", tag, ".rds"))
  }
  bundle <- readRDS(fit_path)
  fit <- bundle$fit
  years <- as.integer(format(as.Date(bundle$weekly_data$week_start), "%Y"))

  core_pars <- c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs", "q")
  raw <- rstan::extract(fit, pars = core_pars, permuted = FALSE, inc_warmup = FALSE)
  draws_arr <- posterior::as_draws_array(raw)
  ps <- as.data.frame(posterior::summarise_draws(draws_arr, posterior::rhat, posterior::ess_bulk, posterior::ess_tail))
  names(ps) <- sub("^posterior::", "", names(ps))
  bfmi <- rstan::get_bfmi(fit)
  divergences <- rstan::get_num_divergent(fit)
  td_hits <- rstan::get_num_max_treedepth(fit)
  hmc_pass <- divergences == 0 && td_hits == 0 && max(ps$rhat, na.rm = TRUE) <= 1.01 && min(ps$ess_bulk, na.rm = TRUE) >= 100 && all(bfmi >= .3)

  q_draws <- rstan::extract(fit, pars = "q")$q
  alpha_R_draws <- rstan::extract(fit, pars = "alpha_R")$alpha_R
  p_sero_draws <- rstan::extract(fit, pars = "p_state_sero_at_anchor")$p_state_sero_at_anchor
  S_prop_draws <- rstan::extract(fit, pars = "S_prop")$S_prop
  sero_pred_draws <- rstan::extract(fit, pars = "sero_pred")$sero_pred
  log_lik_cases_draws <- rstan::extract(fit, pars = "log_lik_cases")$log_lik_cases
  log_lik_serology_draws <- rstan::extract(fit, pars = "log_lik_serology")$log_lik_serology
  expected_cases_draws <- rstan::extract(fit, pars = "expected_reported_cases")$expected_reported_cases

  idx_2021 <- which(years == 2021L)[1]
  idx_2022_start <- which(years == 2022L)[1]
  idx_2025_end <- tail(which(years == 2025L), 1)
  annual_expected <- function(year) {
    idx <- which(years == year)
    rowSums(expected_cases_draws[, idx, drop = FALSE])
  }
  observed_annual <- function(year) sum(bundle$weekly_data$cases[years == year])

  tibble(
    tag = tag, kappa_sero = kappa_numeric, hmc_pass = hmc_pass,
    divergences = divergences, treedepth_hits = td_hits,
    max_rhat = max(ps$rhat, na.rm = TRUE), min_bulk_ess = min(ps$ess_bulk, na.rm = TRUE), min_bfmi = min(bfmi),
    q_median = median(q_draws), q_q025 = quantile(q_draws, .025), q_q975 = quantile(q_draws, .975),
    alpha_R_median = median(alpha_R_draws), alpha_R_q025 = quantile(alpha_R_draws, .025), alpha_R_q975 = quantile(alpha_R_draws, .975),
    exp_alpha_R_median = exp(median(alpha_R_draws)),
    immune_prop_2018_median = median(p_sero_draws), immune_prop_2018_q025 = quantile(p_sero_draws, .025), immune_prop_2018_q975 = quantile(p_sero_draws, .975),
    S_over_N_pre2022_median = median(S_prop_draws[, idx_2022_start - 1]),
    S_over_N_end2025_median = median(S_prop_draws[, idx_2025_end]),
    sero_pred_median = median(sero_pred_draws), sero_pred_q025 = quantile(sero_pred_draws, .025), sero_pred_q975 = quantile(sero_pred_draws, .975),
    observed_2017 = observed_annual(2017L), expected_2017_median = median(annual_expected(2017L)),
    expected_2017_q025 = quantile(annual_expected(2017L), .025), expected_2017_q975 = quantile(annual_expected(2017L), .975),
    observed_2022 = observed_annual(2022L), expected_2022_median = median(annual_expected(2022L)),
    expected_2022_q025 = quantile(annual_expected(2022L), .025), expected_2022_q975 = quantile(annual_expected(2022L), .975),
    case_loglik_sum_median = median(rowSums(log_lik_cases_draws)),
    serology_loglik_median = median(log_lik_serology_draws)
  )
}

run_serology_leverage_analysis <- function() {
  root <- lev_root()
  base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "21_v4_2_q_juazeiro_serology")
  ladder <- list(c("kappa50", 50), c("kappa100", 100), c("kappa200", 200), c("kappa500", 500), c("binomial", Inf))
  comparison <- bind_rows(lapply(ladder, function(x) summarise_one(root, x[[1]], as.numeric(x[[2]]))))
  write_csv(comparison, file.path(base_dir, "serology_leverage_comparison.csv"))
  message("[leverage] comparison table:")
  print(as.data.frame(comparison[, c("tag", "hmc_pass", "divergences", "treedepth_hits", "max_rhat", "q_median", "q_q025", "q_q975", "exp_alpha_R_median", "immune_prop_2018_median")]))

  q_prior_mean <- 16.2 / (16.2 + 108.7)
  juazeiro_observed <- 103 / 404
  plot_data <- comparison |> mutate(log_kappa = ifelse(is.infinite(kappa_sero), log10(1e5), log10(kappa_sero)))

  theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))
  p1 <- ggplot(plot_data, aes(log_kappa, q_median)) +
    geom_ribbon(aes(ymin = q_q025, ymax = q_q975), alpha = .15, fill = "#0072B2") +
    geom_line(colour = "#0072B2") + geom_point(colour = "#0072B2", size = 2) +
    geom_hline(yintercept = q_prior_mean, linetype = 2, colour = "#D55E00") +
    annotate("text", x = min(plot_data$log_kappa), y = q_prior_mean, label = "q prior mean", vjust = -.5, hjust = 0, size = 2.5, colour = "#D55E00") +
    labs(title = "A. q posterior vs log10(kappa_sero)", x = "log10(kappa_sero)  [rightmost = binomial limit]", y = "q") + theme_v4
  p2 <- ggplot(plot_data, aes(log_kappa, immune_prop_2018_median)) +
    geom_ribbon(aes(ymin = immune_prop_2018_q025, ymax = immune_prop_2018_q975), alpha = .15, fill = "#009E73") +
    geom_line(colour = "#009E73") + geom_point(colour = "#009E73", size = 2) +
    geom_hline(yintercept = juazeiro_observed, linetype = 2, colour = "#D55E00") +
    annotate("text", x = min(plot_data$log_kappa), y = juazeiro_observed, label = "Juazeiro observed", vjust = 1.5, hjust = 0, size = 2.5, colour = "#D55E00") +
    labs(title = "B. Model-implied 2018 Ceará immune proportion vs log10(kappa_sero)", x = "log10(kappa_sero)", y = "immune proportion") + theme_v4
  p3 <- ggplot(plot_data, aes(log_kappa, S_over_N_end2025_median)) +
    geom_line(colour = "#CC79A7") + geom_point(colour = "#CC79A7", size = 2) +
    labs(title = "C. S/N at end-2025 vs log10(kappa_sero)", x = "log10(kappa_sero)", y = "S/N") + theme_v4
  p4 <- ggplot(plot_data, aes(log_kappa, exp_alpha_R_median)) +
    geom_line(colour = "#E69F00") + geom_point(colour = "#E69F00", size = 2) +
    labs(title = "D. exp(alpha_R) [baseline R0] vs log10(kappa_sero)", x = "log10(kappa_sero)", y = "exp(alpha_R)") + theme_v4
  figure <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "v4.2 serology leverage/conflict analysis")
  figure_dir <- file.path(root, "03_Output", "figures", "renewal_v4_1_q_only") # reuse existing figures root pattern
  figure_dir <- file.path(root, "03_Output", "figures", "renewal_v4_2_q_juazeiro_serology")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figure_dir, "serology_leverage_curve.png"), figure, width = 220, height = 160, units = "mm", dpi = 300)
  message("[leverage] figure saved: ", file.path(figure_dir, "serology_leverage_curve.png"))
  invisible(comparison)
}

if (sys.nframe() == 0L) run_serology_leverage_analysis()
