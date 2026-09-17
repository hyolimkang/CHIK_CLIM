# v4.4 Section 6-8: HMC diagnostics by chain, block-level randomized
# quantile residual ACF (Section 7), and the direct weekly-vs-4-week
# information-content comparison table (Section 8).

required_packages <- c("here", "rstan", "dplyr", "tibble", "readr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr); library(ggplot2); library(patchwork) })

v4_4_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

randomized_residuals_nb <- function(observed, mu_draws, phi_draws, seed = 20260912L) {
  set.seed(seed)
  n_draws <- length(phi_draws); n_units <- length(observed)
  resid <- matrix(NA_real_, n_draws, n_units)
  for (d in seq_len(n_draws)) {
    mu <- mu_draws[d, ]; phi <- phi_draws[d]
    lower <- ifelse(observed == 0, 0, pnbinom(observed - 1, size = phi, mu = mu))
    upper <- pnbinom(observed, size = phi, mu = mu)
    u <- pmin(pmax(runif(n_units, lower, upper), 1e-9), 1 - 1e-9)
    resid[d, ] <- qnorm(u)
  }
  colMeans(resid)
}

run_analysis <- function(fit_id) {
  root <- v4_4_root()
  base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "23_v4_4_four_week_observation")
  tag <- paste0("fit", fit_id)
  bundle <- readRDS(file.path(base_dir, "outputs", tag, paste0("renewal_ceara_v4_4_fit_", tag, ".rds")))
  fit <- bundle$fit

  # --- HMC + per-chain q/alpha_R (Section 6) ---------------------------------
  core_pars <- c("alpha_R", "z_year", "beta_sin", "beta_cos", "phi_obs_block", "eta_q")
  raw <- rstan::extract(fit, pars = core_pars, permuted = FALSE, inc_warmup = FALSE)
  draws_arr <- posterior::as_draws_array(raw)
  ps <- as.data.frame(posterior::summarise_draws(draws_arr, posterior::rhat, posterior::ess_bulk))
  names(ps) <- sub("^posterior::", "", names(ps))
  bfmi <- rstan::get_bfmi(fit)
  divergences <- rstan::get_num_divergent(fit); td_hits <- rstan::get_num_max_treedepth(fit)
  max_rhat <- max(ps$rhat, na.rm = TRUE); min_ess <- min(ps$ess_bulk, na.rm = TRUE)
  hmc_pass <- divergences == 0 && td_hits == 0 && max_rhat <= 1.01 && min_ess >= 100 && all(bfmi >= .3)
  message(sprintf("[%s] HMC GATE: %s | div=%d td_hits=%d max_rhat=%.4f min_ess=%.1f", tag, if (hmc_pass) "PASS" else "FAIL", divergences, td_hits, max_rhat, min_ess))

  q_by_chain <- rstan::extract(fit, pars = "q", permuted = FALSE, inc_warmup = FALSE)
  alpha_R_by_chain <- rstan::extract(fit, pars = "alpha_R", permuted = FALSE, inc_warmup = FALSE)
  chain_summary <- bind_rows(lapply(1:4, function(c) tibble(
    tag = tag, chain = c, q_mean = mean(q_by_chain[, c, 1]), q_sd = sd(q_by_chain[, c, 1]),
    alpha_R_mean = mean(alpha_R_by_chain[, c, 1]), exp_alpha_R_mean = exp(mean(alpha_R_by_chain[, c, 1]))
  )))
  write_csv(chain_summary, file.path(base_dir, paste0("v4_4_", tag, "_chain_summary.csv")))
  message("[", tag, "] per-chain q / alpha_R (chains 1-2 started low-q, 3-4 started high-q):")
  print(as.data.frame(chain_summary))

  # --- Block residual ACF (Section 7) ----------------------------------------
  mu_block <- rstan::extract(fit, pars = "expected_block")$expected_block
  phi_block <- rstan::extract(fit, pars = "phi_obs_block")$phi_obs_block
  resid <- randomized_residuals_nb(bundle$blocks$observed_block, mu_block, phi_block)
  acf_obj <- acf(resid, lag.max = 6, plot = FALSE)
  acf_table <- tibble(tag = tag, lag = as.integer(acf_obj$lag), acf = as.numeric(acf_obj$acf)) |> filter(lag >= 1)
  write_csv(acf_table, file.path(base_dir, paste0("v4_4_", tag, "_block_residual_acf.csv")))
  message("[", tag, "] block residual ACF (lags 1-6):"); print(as.data.frame(acf_table))

  list(tag = tag, hmc_pass = hmc_pass, max_rhat = max_rhat, min_ess = min_ess, chain_summary = chain_summary,
       acf_table = acf_table, resid = resid, q_all = as.numeric(rstan::extract(fit, pars = "q")$q))
}

run_v4_4_analysis <- function() {
  root <- v4_4_root()
  base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "23_v4_4_four_week_observation")
  resA4 <- run_analysis("A4")
  resB4 <- run_analysis("B4")

  # --- Comparison table (Section 8) -------------------------------------------
  weekly_acf <- read_csv(file.path(root, "02_Script/00_shared/legacy_model_functions/22_v4_3_short_2015_2019_q_calibration/weekly_residual_acf.csv"), show_col_types = FALSE)
  weekly_lag1_low <- weekly_acf$acf[weekly_acf$mode == "low_q" & weekly_acf$lag == 1]
  weekly_lag1_high <- weekly_acf$acf[weekly_acf$mode == "high_q" & weekly_acf$lag == 1]

  comparison <- tibble(
    metric = c("n_observations", "lag1_residual_acf", "hmc_gate", "max_rhat", "min_bulk_ess"),
    weekly_model_low_q = c(313, weekly_lag1_low, NA, NA, NA),
    weekly_model_high_q = c(313, weekly_lag1_high, NA, NA, NA),
    fourweek_model_B4 = c(resB4$chain_summary$q_mean[1] |> (\(x) 78)(), mean(resB4$acf_table$acf[resB4$acf_table$lag == 1]),
                            resB4$hmc_pass, resB4$max_rhat, resB4$min_ess)
  )
  write_csv(comparison, file.path(base_dir, "four_week_sensitivity_comparison.csv"))
  message("[v4.4] weekly vs 4-week comparison:"); print(as.data.frame(comparison))

  # --- Figure: block residual ACF + q posterior by fit ------------------------
  theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))
  acf_combined <- bind_rows(resA4$acf_table, resB4$acf_table)
  p1 <- ggplot(acf_combined, aes(lag, acf, fill = tag)) + geom_col(position = "dodge") +
    geom_hline(yintercept = c(-1, 1) * 1.96 / sqrt(78), linetype = 2, colour = "grey40") +
    labs(title = "A. 4-week block residual ACF", x = "lag (blocks)", y = "ACF") + theme_v4
  q_df <- bind_rows(tibble(fit = "A4", q = resA4$q_all), tibble(fit = "B4", q = resB4$q_all))
  p2 <- ggplot(q_df, aes(q, fill = fit)) + geom_histogram(position = "identity", alpha = .5, bins = 60) +
    labs(title = "B. q posterior under 4-week likelihood", x = "q", y = "count") + theme_v4
  figure <- p1 / p2
  figure_dir <- file.path(root, "03_Output", "figures", "renewal_v4_4_four_week_observation")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figure_dir, "block_residual_acf_and_q_posterior.png"), figure, width = 180, height = 180, units = "mm", dpi = 300)
  message("[v4.4] figure saved: ", file.path(figure_dir, "block_residual_acf_and_q_posterior.png"))

  invisible(list(resA4 = resA4, resB4 = resB4, comparison = comparison))
}

if (sys.nframe() == 0L) run_v4_4_analysis()
