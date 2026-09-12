# =============================================================================
# Diagnose the AR(1)-free normalized annual-FOI revision.
#
# This script does not refit anything.  It summarizes the six saved posterior
# fits, applies the prespecified computational gate, and writes figures/tables
# into the v2-specific output directories.
# =============================================================================

required_packages <- c("here", "rstan", "posterior", "dplyr", "tidyr",
                       "readr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(rstan)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

project_root_v2 <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not locate the inner CHIK_CLIM R project.")
}

qfun <- function(x, probability) {
  unname(stats::quantile(x, probability, na.rm = TRUE))
}

summarise_draw <- function(x) {
  c(q025 = qfun(x, .025), median = qfun(x, .5), q975 = qfun(x, .975))
}

hmc_diagnostics_v2 <- function(fit, elapsed_seconds, max_treedepth = 12L) {
  draws_array <- posterior::as_draws_array(rstan::extract(fit, permuted = FALSE))
  draw_summary <- posterior::summarise_draws(
    draws_array, posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  ) |>
    rename(
      rhat = `posterior::rhat`,
      ess_bulk = `posterior::ess_bulk`,
      ess_tail = `posterior::ess_tail`
    ) |>
    filter(!grepl("^(C_pred|lp__)", variable))

  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(sampler, function(x) sum(x[, "divergent__"]), numeric(1)))
  treedepth_hits <- sum(vapply(
    sampler, function(x) sum(x[, "treedepth__"] >= max_treedepth), numeric(1)
  ))
  bfmi <- vapply(
    sampler,
    function(x) mean(diff(x[, "energy__"])^2) / stats::var(x[, "energy__"]),
    numeric(1)
  )
  headline <- tibble(
    hmc_gate = ifelse(
      divergences == 0 && treedepth_hits == 0 &&
        max(draw_summary$rhat, na.rm = TRUE) <= 1.01 &&
        min(draw_summary$ess_bulk, na.rm = TRUE) >= 400 &&
        min(draw_summary$ess_tail, na.rm = TRUE) >= 400 && all(bfmi >= .3),
      "PASS", "FAIL"
    ),
    divergences = divergences,
    max_treedepth_hits = treedepth_hits,
    max_rhat = max(draw_summary$rhat, na.rm = TRUE),
    n_rhat_gt_1_01 = sum(draw_summary$rhat > 1.01, na.rm = TRUE),
    min_bulk_ess = min(draw_summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(draw_summary$ess_tail, na.rm = TRUE),
    min_bfmi = min(bfmi),
    elapsed_seconds = elapsed_seconds
  )
  list(
    headline = headline,
    parameter_summary = draw_summary,
    bfmi = tibble(chain = seq_along(bfmi), bfmi = bfmi)
  )
}

annual_summary_v2 <- function(bundle, fit_id) {
  dr <- rstan::extract(bundle$fit, permuted = TRUE)
  annual <- bundle$prepared$annual
  years <- annual$year
  n_year <- length(years)
  expected_cases <- if (!is.null(dr$expected_cases)) {
    dr$expected_cases
  } else {
    # The conditional-on-total model has no absolute detection parameter.
    # Its fitted expected allocation is total observed cases times case_shape.
    dr$case_shape * sum(annual$observed_cases)
  }

  tibble(
    fit_id = fit_id,
    year = years,
    observed_cases = annual$observed_cases,
    lambda_q025 = apply(dr$lambda, 2, qfun, probability = .025),
    lambda_median = apply(dr$lambda, 2, qfun, probability = .5),
    lambda_q975 = apply(dr$lambda, 2, qfun, probability = .975),
    weight_q025 = apply(dr$weight, 2, qfun, probability = .025),
    weight_median = apply(dr$weight, 2, qfun, probability = .5),
    weight_q975 = apply(dr$weight, 2, qfun, probability = .975),
    attack_prob_q025 = apply(dr$attack_prob, 2, qfun, probability = .025),
    attack_prob_median = apply(dr$attack_prob, 2, qfun, probability = .5),
    attack_prob_q975 = apply(dr$attack_prob, 2, qfun, probability = .975),
    latent_infections_q025 = apply(dr$X, 2, qfun, probability = .025),
    latent_infections_median = apply(dr$X, 2, qfun, probability = .5),
    latent_infections_q975 = apply(dr$X, 2, qfun, probability = .975),
    S_start_prop_q025 = apply(sweep(dr$S_start, 2, annual$N_start, "/"), 2, qfun, probability = .025),
    S_start_prop_median = apply(sweep(dr$S_start, 2, annual$N_start, "/"), 2, qfun, probability = .5),
    S_start_prop_q975 = apply(sweep(dr$S_start, 2, annual$N_start, "/"), 2, qfun, probability = .975),
    S_end_prop_q025 = apply(dr$S_end_prop, 2, qfun, probability = .025),
    S_end_prop_median = apply(dr$S_end_prop, 2, qfun, probability = .5),
    S_end_prop_q975 = apply(dr$S_end_prop, 2, qfun, probability = .975),
    expected_cases_median = apply(expected_cases, 2, qfun, probability = .5),
    predicted_cases_q025 = apply(dr$C_pred, 2, qfun, probability = .025),
    predicted_cases_median = apply(dr$C_pred, 2, qfun, probability = .5),
    predicted_cases_q975 = apply(dr$C_pred, 2, qfun, probability = .975)
  )
}

trace_data_v2 <- function(fit, parameters) {
  sampled <- rstan::extract(fit, pars = parameters, permuted = FALSE,
                            inc_warmup = FALSE)
  bind_rows(lapply(seq_along(parameters), function(i) {
    tidyr::expand_grid(
      iteration = seq_len(dim(sampled)[1]), chain = seq_len(dim(sampled)[2])
    ) |>
      mutate(parameter = parameters[i], value = as.vector(sampled[, , i]))
  }))
}

diagnose_annual_foi_shape_v2 <- function(primary_fit = "M1v2_SHAPE_s0_75") {
  root <- project_root_v2()
  table_dir <- file.path(root, "03_Output", "tables", "annual_foi_shape_v2")
  figure_dir <- file.path(root, "03_Output", "figures", "annual_foi_shape_v2")
  stan_dir <- file.path(root, "02_Script", "stan")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  paths <- list.files(stan_dir, pattern = "^annual_foi_shape_v2_.*\\.rds$",
                      full.names = TRUE)
  if (!length(paths)) stop("No v2 posterior fits were found in ", stan_dir)
  fit_ids <- sub("^annual_foi_shape_v2_(.*)\\.rds$", "\\1", basename(paths))
  bundles <- stats::setNames(lapply(paths, readRDS), fit_ids)

  hmc_all <- list()
  annual_all <- list()
  scalar_all <- list()
  for (fit_id in names(bundles)) {
    bundle <- bundles[[fit_id]]
    hmc <- hmc_diagnostics_v2(bundle$fit, bundle$elapsed_seconds)
    hmc_all[[fit_id]] <- mutate(hmc$headline, fit_id = fit_id, .before = 1)
    write_csv(hmc$bfmi, file.path(table_dir, paste0(fit_id, "_bfmi.csv")))
    write_csv(arrange(hmc$parameter_summary, desc(rhat)) |> slice_head(n = 30),
              file.path(table_dir, paste0(fit_id, "_worst_rhat.csv")))
    write_csv(arrange(hmc$parameter_summary, ess_bulk) |> slice_head(n = 30),
              file.path(table_dir, paste0(fit_id, "_worst_bulk_ess.csv")))
    annual_all[[fit_id]] <- annual_summary_v2(bundle, fit_id)
    dr <- rstan::extract(bundle$fit, permuted = TRUE)
    scalar_all[[fit_id]] <- tibble(
      fit_id = fit_id,
      observation_model = bundle$config$model,
      lambda_bar_q025 = qfun(dr$lambda_bar, .025),
      lambda_bar_median = qfun(dr$lambda_bar, .5),
      lambda_bar_q975 = qfun(dr$lambda_bar, .975),
      sigma_year_q025 = qfun(dr$sigma_year, .025),
      sigma_year_median = qfun(dr$sigma_year, .5),
      sigma_year_q975 = qfun(dr$sigma_year, .975),
      cumulative_infections_median = qfun(dr$cumulative_infections, .5),
      S_end_2025_median = qfun(dr$S_end_prop[, ncol(dr$S_end_prop)], .5),
      q_implied_median = qfun(dr$q_implied, .5),
      q_implied_gt_one_probability = mean(dr$q_implied_gt_one),
      q_median = if (!is.null(dr$q)) qfun(dr$q, .5) else NA_real_,
      phi_median = if (!is.null(dr$phi)) qfun(dr$phi, .5) else NA_real_,
      kappa_median = if (!is.null(dr$kappa)) qfun(dr$kappa, .5) else NA_real_
    )
  }
  hmc_table <- bind_rows(hmc_all)
  annual_table <- bind_rows(annual_all)
  scalar_table <- bind_rows(scalar_all) |> left_join(hmc_table, by = "fit_id")
  write_csv(hmc_table, file.path(table_dir, "annual_foi_shape_v2_hmc_gate.csv"))
  write_csv(scalar_table, file.path(table_dir, "annual_foi_shape_v2_model_comparison.csv"))
  write_csv(annual_table, file.path(table_dir, "annual_foi_shape_v2_annual_posterior_summary.csv"))

  if (!primary_fit %in% names(bundles)) stop("Primary fit not found: ", primary_fit)
  primary_bundle <- bundles[[primary_fit]]
  primary <- filter(annual_table, fit_id == primary_fit)
  dr <- rstan::extract(primary_bundle$fit, permuted = TRUE)
  anchor <- primary_bundle$prepared$foi_anchor$median

  p_cases <- ggplot(primary, aes(year, observed_cases)) +
    geom_col(fill = "#0072B2") +
    geom_ribbon(aes(ymin = predicted_cases_q025, ymax = predicted_cases_q975),
                fill = "grey60", alpha = .45) +
    geom_line(aes(y = predicted_cases_median), colour = "#D55E00", linewidth = .8) +
    labs(title = "A. Annual cases: observed and posterior predictive", y = "Reported cases", x = NULL) +
    theme_classic()
  p_lambda <- ggplot(primary, aes(year, lambda_median)) +
    geom_ribbon(aes(ymin = lambda_q025, ymax = lambda_q975), fill = "#56B4E9", alpha = .35) +
    geom_line(linewidth = .8) + geom_hline(yintercept = anchor, linetype = 2, colour = "#D55E00") +
    scale_y_continuous(trans = "log10") +
    labs(title = "B. Annual FOI", y = "FOI (log scale)", x = NULL) + theme_classic()
  p_weight <- ggplot(primary, aes(year, weight_median)) +
    geom_ribbon(aes(ymin = weight_q025, ymax = weight_q975), fill = "#CC79A7", alpha = .35) +
    geom_line(linewidth = .8) + geom_hline(yintercept = 1, linetype = 2) +
    labs(title = "C. Mean-preserving annual FOI multiplier", y = "Multiplier", x = NULL) + theme_classic()
  p_attack <- ggplot(primary, aes(year, attack_prob_median)) +
    geom_ribbon(aes(ymin = attack_prob_q025, ymax = attack_prob_q975), fill = "#E69F00", alpha = .35) +
    geom_line(linewidth = .8) + labs(title = "D. Annual attack probability", y = "Probability", x = NULL) + theme_classic()
  p_s <- ggplot(primary, aes(year, S_end_prop_median)) +
    geom_ribbon(aes(ymin = S_end_prop_q025, ymax = S_end_prop_q975), fill = "#009E73", alpha = .35) +
    geom_line(linewidth = .8) + coord_cartesian(ylim = c(0, 1)) +
    labs(title = "E. Susceptible proportion at year-end", y = "S/N", x = NULL) + theme_classic()
  p_inf <- ggplot(primary, aes(year, latent_infections_median)) +
    geom_ribbon(aes(ymin = latent_infections_q025, ymax = latent_infections_q975), fill = "#999999", alpha = .4) +
    geom_line(linewidth = .8) + labs(title = "F. Posterior latent infections", y = "Infections", x = NULL) + theme_classic()
  q_data <- tibble(q_implied = dr$q_implied)
  p_q <- ggplot(q_data, aes(q_implied)) + geom_histogram(bins = 35, fill = "#56B4E9") +
    geom_vline(xintercept = 1, colour = "#D55E00", linetype = 2) +
    labs(title = "G. Implied infection-to-case fraction", x = "q implied", y = "Posterior draws") + theme_classic()
  prior_post <- bind_rows(
    tibble(source = "Prior", value = exp(rnorm(4000, primary_bundle$prepared$foi_anchor$log_mean,
                                               primary_bundle$prepared$foi_anchor$log_sd))),
    tibble(source = "Posterior", value = dr$lambda_bar)
  )
  p_anchor <- ggplot(prior_post, aes(value, fill = source)) + geom_density(alpha = .35) +
    scale_x_continuous(trans = "log10") +
    labs(title = "H. Long-term FOI anchor: prior vs posterior", x = "lambda_bar (log scale)", y = "Density") +
    theme_classic() + theme(legend.position = "bottom")
  main_plot <- (p_cases | p_lambda) / (p_weight | p_attack) / (p_s | p_inf) / (p_q | p_anchor)
  ggsave(file.path(figure_dir, "annual_foi_shape_v2_primary_posterior_diagnostics.pdf"),
         main_plot, width = 250, height = 310, units = "mm", device = cairo_pdf)

  sensitivity <- scalar_table |>
    select(fit_id, observation_model, lambda_bar_median, sigma_year_median,
           cumulative_infections_median, S_end_2025_median, q_implied_median,
           hmc_gate) |>
    pivot_longer(c(lambda_bar_median, sigma_year_median, cumulative_infections_median,
                   S_end_2025_median, q_implied_median), names_to = "parameter", values_to = "value")
  p_sens <- ggplot(sensitivity, aes(reorder(fit_id, value), value, colour = hmc_gate)) +
    geom_point(size = 2.5) + coord_flip() + facet_wrap(~parameter, scales = "free_x", ncol = 1) +
    labs(title = "Posterior sensitivity across v2 specifications", x = NULL, y = NULL) + theme_classic()
  ggsave(file.path(figure_dir, "annual_foi_shape_v2_sensitivity.pdf"), p_sens,
         width = 230, height = 240, units = "mm", device = cairo_pdf)

  trace_parameters <- c("lambda_bar", "sigma_year", "kappa", "z_year[1]", "z_year[3]", "z_year[8]")
  trace_plot <- ggplot(trace_data_v2(primary_bundle$fit, trace_parameters),
                       aes(iteration, value, colour = factor(chain))) +
    geom_line(linewidth = .18) + facet_wrap(~parameter, scales = "free_y", ncol = 3) +
    labs(title = "Primary v2 SHAPE model trace plots", colour = "Chain") +
    theme_classic(base_size = 8)
  ggsave(file.path(figure_dir, "annual_foi_shape_v2_primary_traceplots.pdf"), trace_plot,
         width = 250, height = 160, units = "mm", device = cairo_pdf)

  invisible(list(hmc = hmc_table, scalars = scalar_table, annual = annual_table))
}

if (sys.nframe() == 0L) diagnose_annual_foi_shape_v2()
