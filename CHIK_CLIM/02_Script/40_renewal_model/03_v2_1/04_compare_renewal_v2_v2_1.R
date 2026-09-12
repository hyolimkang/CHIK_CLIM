# =============================================================================
# 04_compare_renewal_v2_v2_1.R
#
# Read-only posterior diagnostics comparing saved Ceara renewal-model v2 and
# v2.1 fit bundles. No Stan model is compiled and no sampling is performed.
#
# Environment overrides:
#   RENEWAL_V2_COMPARE_FIT, RENEWAL_V2_1_COMPARE_FIT,
#   RENEWAL_V2_COMPARE_TAG
# =============================================================================

required_packages <- c("here", "rstan", "ggplot2", "dplyr", "tidyr", "tibble", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
})

COLORS <- c("v2" = "#0072B2", "v2.1" = "#D55E00")

fit_path <- function(env_name, default_filename) {
  override <- Sys.getenv(env_name, unset = "")
  if (nzchar(override)) {
    return(override)
  }
  here::here("02_Script/stan", default_filename)
}

read_fit_bundle <- function(path, label) {
  if (!file.exists(path)) stop(label, " fit bundle not found: ", path)
  bundle <- readRDS(path)
  if (!is.list(bundle) || !inherits(bundle$fit, "stanfit") ||
    is.null(bundle$stan_data) || is.null(bundle$weekly_data)) {
    stop(label, " file is not a renewal fit bundle with fit, stan_data, and weekly_data")
  }
  bundle
}

as_draw_matrix <- function(x) {
  # rstan::extract(fit, pars = ...) returns even a scalar parameter as a
  # 1-D array (dim = n_draws, not NULL) rather than a plain vector, so the
  # dimensionality check must look at length(dim(x)), not just is.null().
  if (is.null(dim(x)) || length(dim(x)) == 1L) {
    return(matrix(as.vector(x), ncol = 1L))
  }
  x
}

summarise_draw_matrix <- function(draw_matrix, dates, version, quantity) {
  draw_matrix <- as_draw_matrix(draw_matrix)
  if (ncol(draw_matrix) != length(dates)) {
    stop("Draw matrix and weekly dates have incompatible dimensions for ", quantity)
  }
  quantiles <- t(vapply(seq_len(ncol(draw_matrix)), function(t) {
    stats::quantile(draw_matrix[, t], probs = c(0.025, 0.5, 0.975))
  }, numeric(3)))
  tibble::tibble(
    version = version,
    quantity = quantity,
    week_start = dates,
    q025 = quantiles[, 1],
    median = quantiles[, 2],
    q975 = quantiles[, 3]
  )
}

extract_v2 <- function(bundle) {
  pars <- c(
    "X", "S_prop", "immune_prop", "rho_sym_t", "rho_total_t",
    "R0_t", "R_eff_t", "C_pred", "expected_reported_cases",
    "sero_positive_pred", "sero_geographic_offset", "rho_sym_ceara_mid",
    "log_infection_scale", "sigma_R", "reporting_trend", "phi_obs"
  )
  draws <- rstan::extract(bundle$fit, pars = pars)
  n_draws <- nrow(as_draw_matrix(draws$C_pred))
  N_pop <- bundle$stan_data$N_pop
  if (length(N_pop) != 1L || !is.finite(N_pop) || N_pop <= 0) {
    stop("v2 bundle has an invalid scalar N_pop")
  }
  list(
    label = "v2",
    dates = bundle$weekly_data$week_start,
    observed_cases = bundle$weekly_data$cases,
    seed_weeks = bundle$stan_data$seed_weeks,
    susceptible_prop = draws$S_prop + draws$X / N_pop,
    cumulative_prop = draws$immune_prop,
    reporting = draws$rho_sym_t,
    overall_detection = draws$rho_total_t,
    R0 = draws$R0_t,
    Reff = draws$R_eff_t,
    C_pred = draws$C_pred,
    expected_cases = draws$expected_reported_cases,
    p_symp = rep(bundle$config$p_symp, n_draws),
    p_symp_fixed = TRUE,
    sero_pred_counts = as_draw_matrix(draws$sero_positive_pred),
    sero_sites = "Juazeiro do Norte",
    sero_positive = bundle$stan_data$sero_positive,
    sero_n = bundle$stan_data$sero_tested,
    site_offsets = as_draw_matrix(draws$sero_geographic_offset),
    sigma_geo = rep(bundle$config$sero_geographic_sd, n_draws),
    sigma_geo_fixed = TRUE,
    scale_draws = list(
      log_infection_scale = draws$log_infection_scale,
      sigma_R = draws$sigma_R,
      rho_sym_ceara_mid = draws$rho_sym_ceara_mid,
      reporting_trend = draws$reporting_trend,
      phi_obs = draws$phi_obs,
      sero_geographic_offset = draws$sero_geographic_offset
    ),
    fit = bundle$fit
  )
}

extract_v2_1 <- function(bundle) {
  pars <- c(
    "S", "immune_prop", "rho_sym_t", "overall_detection", "R0_t", "R_eff_t",
    "C_pred", "expected_reported_cases", "p_symp", "sero_positive_pred",
    "sero_geographic_offset", "sigma_geo", "rho_sym_ceara_mid",
    "log_infection_scale", "sigma_R", "reporting_trend", "phi_obs"
  )
  draws <- rstan::extract(bundle$fit, pars = pars)
  N_pop_t <- bundle$stan_data$N_pop_t
  if (length(N_pop_t) != ncol(draws$S) || any(!is.finite(N_pop_t)) || any(N_pop_t <= 0)) {
    stop("v2.1 bundle has an invalid weekly N_pop_t")
  }
  list(
    label = "v2.1",
    dates = bundle$weekly_data$week_start,
    observed_cases = bundle$weekly_data$cases,
    seed_weeks = bundle$stan_data$seed_weeks,
    susceptible_prop = sweep(draws$S, 2, N_pop_t, "/"),
    cumulative_prop = draws$immune_prop,
    reporting = draws$rho_sym_t,
    overall_detection = draws$overall_detection,
    R0 = draws$R0_t,
    Reff = draws$R_eff_t,
    C_pred = draws$C_pred,
    expected_cases = draws$expected_reported_cases,
    p_symp = draws$p_symp,
    p_symp_fixed = FALSE,
    sero_pred_counts = as_draw_matrix(draws$sero_positive_pred),
    sero_sites = bundle$config$serology_sites,
    sero_positive = bundle$config$serology_positive,
    sero_n = bundle$config$serology_n,
    site_offsets = as_draw_matrix(draws$sero_geographic_offset),
    sigma_geo = draws$sigma_geo,
    sigma_geo_fixed = FALSE,
    scale_draws = list(
      log_infection_scale = draws$log_infection_scale,
      sigma_R = draws$sigma_R,
      p_symp = draws$p_symp,
      rho_sym_ceara_mid = draws$rho_sym_ceara_mid,
      reporting_trend = draws$reporting_trend,
      phi_obs = draws$phi_obs,
      sigma_geo = draws$sigma_geo
    ),
    fit = bundle$fit
  )
}

make_serology_draws <- function(model) {
  if (ncol(model$sero_pred_counts) != length(model$sero_sites) ||
    length(model$sero_sites) != length(model$sero_n)) {
    stop(model$label, " serology draws and metadata have incompatible dimensions")
  }
  dplyr::bind_rows(lapply(seq_along(model$sero_sites), function(j) {
    tibble::tibble(
      version = model$label,
      site = model$sero_sites[j],
      predicted_prevalence = model$sero_pred_counts[, j] / model$sero_n[j]
    )
  }))
}

make_serology_observed <- function(model) {
  tibble::tibble(
    version = model$label,
    site = model$sero_sites,
    observed_prevalence = model$sero_positive / model$sero_n,
    lower = vapply(seq_along(model$sero_sites), function(j) {
      stats::binom.test(model$sero_positive[j], model$sero_n[j])$conf.int[1]
    }, numeric(1)),
    upper = vapply(seq_along(model$sero_sites), function(j) {
      stats::binom.test(model$sero_positive[j], model$sero_n[j])$conf.int[2]
    }, numeric(1))
  )
}

make_site_offset_draws <- function(model) {
  dplyr::bind_rows(lapply(seq_along(model$sero_sites), function(j) {
    tibble::tibble(
      version = model$label,
      site = model$sero_sites[j],
      log_odds_offset = model$site_offsets[, j]
    )
  }))
}

parameter_diagnostics <- function(fit, version) {
  summary_matrix <- summary(fit)$summary
  tibble::tibble(
    version = version,
    parameter = rownames(summary_matrix),
    ess = summary_matrix[, "n_eff"],
    rhat = summary_matrix[, "Rhat"]
  )
}

residual_diagnostics <- function(model) {
  phi_obs_mean <- mean(model$scale_draws$phi_obs)
  expected_mean <- colMeans(as_draw_matrix(model$expected_cases))
  keep <- seq_along(model$dates) > model$seed_weeks
  tibble::tibble(
    version = model$label,
    week_start = model$dates[keep],
    observed = model$observed_cases[keep],
    expected = expected_mean[keep]
  ) |>
    dplyr::mutate(
      pearson_resid = (observed - expected) / sqrt(expected + expected^2 / phi_obs_mean)
    )
}

coverage_summary <- function(model) {
  draw_matrix <- as_draw_matrix(model$C_pred)
  quantiles <- t(vapply(seq_len(ncol(draw_matrix)), function(t) {
    stats::quantile(draw_matrix[, t], probs = c(0.025, 0.25, 0.75, 0.975))
  }, numeric(4)))
  keep <- seq_along(model$dates) > model$seed_weeks
  tibble::tibble(
    version = model$label,
    coverage_50 = mean(
      model$observed_cases[keep] >= quantiles[keep, 2] & model$observed_cases[keep] <= quantiles[keep, 3]
    ),
    coverage_95 = mean(
      model$observed_cases[keep] >= quantiles[keep, 1] & model$observed_cases[keep] <= quantiles[keep, 4]
    )
  )
}

sampler_diagnostics <- function(fit, version) {
  parameter_summary <- summary(fit)$summary
  bfmi <- rstan::get_bfmi(fit)
  tibble::tibble(
    version = version,
    max_rhat = max(parameter_summary[, "Rhat"], na.rm = TRUE),
    n_rhat_gt_1.01 = sum(parameter_summary[, "Rhat"] > 1.01, na.rm = TRUE),
    min_ess = min(parameter_summary[, "n_eff"], na.rm = TRUE),
    median_ess = stats::median(parameter_summary[, "n_eff"], na.rm = TRUE),
    divergences = rstan::get_num_divergent(fit),
    max_treedepth_hits = rstan::get_num_max_treedepth(fit),
    min_bfmi = min(bfmi),
    median_bfmi = stats::median(bfmi)
  )
}

correlation_long <- function(scale_draws, version) {
  draw_matrix <- do.call(cbind, scale_draws)
  correlation <- stats::cor(draw_matrix)
  as.data.frame(as.table(correlation), stringsAsFactors = FALSE) |>
    dplyr::rename(parameter_x = Var1, parameter_y = Var2, correlation = Freq) |>
    dplyr::mutate(version = version)
}

plot_trajectory <- function(data, quantity, title, y_label, limits = NULL) {
  plot_data <- dplyr::filter(data, .data$quantity == quantity)
  ggplot(plot_data, aes(week_start, median, colour = version, fill = version)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.13, colour = NA) +
    geom_line(linewidth = 0.6) +
    scale_colour_manual(values = COLORS) +
    scale_fill_manual(values = COLORS) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = limits) +
    labs(title = title, x = NULL, y = y_label, colour = "Model", fill = "Model") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")
}

if (sys.nframe() == 0) {
  v2_path <- fit_path("RENEWAL_V2_COMPARE_FIT", "renewal_ceara_v2_fit.rds")
  v2_1_path <- fit_path("RENEWAL_V2_1_COMPARE_FIT", "renewal_ceara_v2_1_fit.rds")
  v2 <- extract_v2(read_fit_bundle(v2_path, "v2"))
  v2_1 <- extract_v2_1(read_fit_bundle(v2_1_path, "v2.1"))
  models <- list(v2, v2_1)

  trajectories <- dplyr::bind_rows(lapply(models, function(model) {
    dplyr::bind_rows(
      summarise_draw_matrix(model$C_pred, model$dates, model$label, "posterior_predictive_cases"),
      summarise_draw_matrix(model$susceptible_prop, model$dates, model$label, "susceptible_proportion"),
      summarise_draw_matrix(model$cumulative_prop, model$dates, model$label, "cumulative_proportion"),
      summarise_draw_matrix(model$reporting, model$dates, model$label, "symptomatic_reporting"),
      summarise_draw_matrix(model$overall_detection, model$dates, model$label, "overall_detection"),
      summarise_draw_matrix(model$R0, model$dates, model$label, "R0"),
      summarise_draw_matrix(model$Reff, model$dates, model$label, "Reff")
    )
  }))
  observed_cases <- dplyr::bind_rows(lapply(models, function(model) {
    tibble::tibble(version = model$label, week_start = model$dates, observed = model$observed_cases)
  }))

  p_cases <- ggplot(
    dplyr::filter(trajectories, quantity == "posterior_predictive_cases"),
    aes(week_start, median)
  ) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#56B4E9", alpha = 0.25) +
    geom_line(colour = "#0072B2", linewidth = 0.6) +
    geom_point(data = observed_cases, aes(week_start, observed), inherit.aes = FALSE, size = 0.45) +
    facet_wrap(~version, ncol = 1, scales = "free_x") +
    scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
    labs(title = "Observed versus posterior-predicted weekly cases", x = NULL, y = "Cases") +
    theme_classic(base_size = 10)

  p_states <- plot_trajectory(
    trajectories, "susceptible_proportion", "Susceptible proportion", "Proportion", c(0, 1)
  ) + plot_trajectory(
    trajectories, "cumulative_proportion", "Cumulative proportion infected", "Proportion", c(0, 1)
  ) + patchwork::plot_layout(guides = "collect")

  p_reporting <- plot_trajectory(
    trajectories, "symptomatic_reporting", "Symptomatic-case reporting probability", "Probability", c(0, 1)
  ) + plot_trajectory(
    trajectories, "overall_detection", "Overall infection-to-report detection", "Probability", c(0, 1)
  )
  p_symp_draws <- dplyr::bind_rows(lapply(models, function(model) {
    tibble::tibble(version = model$label, p_symp = model$p_symp, fixed = model$p_symp_fixed)
  }))
  p_p_symp <- ggplot(dplyr::filter(p_symp_draws, !fixed), aes(p_symp, colour = version, fill = version)) +
    geom_density(alpha = 0.18, linewidth = 0.6) +
    geom_vline(
      data = dplyr::filter(p_symp_draws, fixed) |> dplyr::distinct(version, p_symp),
      aes(xintercept = p_symp, colour = version), linetype = "dashed", linewidth = 0.7
    ) +
    scale_colour_manual(values = COLORS) +
    scale_fill_manual(values = COLORS) +
    labs(title = "Symptomatic probability", subtitle = "v2 is fixed; v2.1 is estimated", x = expression(p[symp]), y = "Posterior density") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")
  p_case_ascertainment <- (p_reporting | p_p_symp) + patchwork::plot_layout(guides = "collect")

  reproduction <- dplyr::bind_rows(
    dplyr::filter(trajectories, quantity == "R0") |> dplyr::mutate(number = "R0"),
    dplyr::filter(trajectories, quantity == "Reff") |> dplyr::mutate(number = "Reff")
  ) |>
    dplyr::group_by(version, number) |>
    dplyr::mutate(t = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::filter(t > ifelse(version == "v2", v2$seed_weeks, v2_1$seed_weeks))
  p_reproduction <- ggplot(reproduction, aes(week_start, median, colour = version, fill = version)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.13, colour = NA) +
    geom_line(linewidth = 0.6) +
    facet_wrap(~number, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = COLORS) +
    scale_fill_manual(values = COLORS) +
    labs(title = "Basic and effective reproduction numbers", x = NULL, y = "Reproduction number", colour = "Model", fill = "Model") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")

  serology_draws <- dplyr::bind_rows(lapply(models, make_serology_draws))
  serology_observed <- dplyr::bind_rows(lapply(models, make_serology_observed))
  p_serology <- ggplot(serology_draws, aes(predicted_prevalence, colour = version, fill = version)) +
    geom_density(alpha = 0.17, linewidth = 0.55) +
    geom_vline(
      data = serology_observed,
      aes(xintercept = observed_prevalence), colour = "black", linetype = "dashed"
    ) +
    facet_grid(site ~ version, scales = "free_y") +
    scale_colour_manual(values = COLORS) +
    scale_fill_manual(values = COLORS) +
    scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
    labs(title = "Posterior-predicted site seroprevalence", subtitle = "Dashed line: observed proportion", x = "Seroprevalence", y = "Posterior density") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")

  site_offsets <- dplyr::bind_rows(lapply(models, make_site_offset_draws))
  sigma_draws <- dplyr::bind_rows(lapply(models, function(model) {
    tibble::tibble(version = model$label, sigma_geo = model$sigma_geo, fixed = model$sigma_geo_fixed)
  }))
  p_offsets <- ggplot(site_offsets, aes(log_odds_offset, colour = version, fill = version)) +
    geom_density(alpha = 0.17, linewidth = 0.55) +
    facet_wrap(~site, scales = "free_y") +
    scale_colour_manual(values = COLORS) +
    scale_fill_manual(values = COLORS) +
    labs(title = "Site-level geographic log-odds offsets", x = "Site random effect", y = "Posterior density") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")
  p_sigma_geo <- ggplot(dplyr::filter(sigma_draws, !fixed), aes(sigma_geo, colour = version, fill = version)) +
    geom_density(alpha = 0.17, linewidth = 0.55) +
    geom_vline(
      data = dplyr::filter(sigma_draws, fixed) |> dplyr::distinct(version, sigma_geo),
      aes(xintercept = sigma_geo, colour = version), linetype = "dashed", linewidth = 0.7
    ) +
    scale_colour_manual(values = COLORS) +
    scale_fill_manual(values = COLORS) +
    labs(title = expression(sigma[geo]), subtitle = "v2 is fixed; v2.1 is estimated", x = expression(sigma[geo]), y = "Posterior density") +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")
  p_serology_geo <- (p_serology / (p_offsets | p_sigma_geo)) + patchwork::plot_layout(heights = c(1, 1))

  parameter_diag <- dplyr::bind_rows(lapply(models, function(model) {
    parameter_diagnostics(model$fit, model$label)
  }))
  sampler_diag <- dplyr::bind_rows(lapply(models, function(model) {
    sampler_diagnostics(model$fit, model$label)
  }))

  # ---- Posterior predictive checks: residuals and interval coverage -------
  residuals_data <- dplyr::bind_rows(lapply(models, residual_diagnostics))
  coverage_data <- dplyr::bind_rows(lapply(models, coverage_summary))
  message("[PPC coverage] (nominal 50%/95%, post-seed weeks only)")
  print(coverage_data)

  p_residuals_series <- ggplot(residuals_data, aes(week_start, pearson_resid, colour = version)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_point(size = 0.5, alpha = 0.6) +
    facet_wrap(~version, ncol = 1) +
    scale_colour_manual(values = COLORS) +
    labs(
      title = "Posterior predictive check: Pearson residuals", subtitle = "Weeks after model seeding only",
      x = NULL, y = "Pearson residual"
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
  p_residuals_hist <- ggplot(residuals_data, aes(pearson_resid, fill = version)) +
    geom_histogram(bins = 30, alpha = 0.75, position = "identity") +
    geom_vline(xintercept = 0, linetype = "22", colour = "grey40") +
    facet_wrap(~version, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = COLORS) +
    labs(title = "Residual distribution", x = "Pearson residual", y = "Weeks") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
  p_ppc_residuals <- p_residuals_series | p_residuals_hist

  # ---- v2.1-only: p_symp identifiability over time -------------------------
  # p_symp is fixed in v2, so this trade-off against the time-varying
  # reporting quantities only exists to check in v2.1.
  p_symp_draws <- v2_1$scale_draws$p_symp
  reporting_draws <- as_draw_matrix(v2_1$reporting)
  detection_draws <- as_draw_matrix(v2_1$overall_detection)
  p_symp_time_correlation <- tibble::tibble(
    week_start = v2_1$dates,
    p_symp_vs_reporting = vapply(seq_len(ncol(reporting_draws)), function(t) {
      stats::cor(p_symp_draws, reporting_draws[, t])
    }, numeric(1)),
    p_symp_vs_detection = vapply(seq_len(ncol(detection_draws)), function(t) {
      stats::cor(p_symp_draws, detection_draws[, t])
    }, numeric(1))
  ) |>
    tidyr::pivot_longer(-week_start, names_to = "pair", values_to = "correlation")
  p_p_symp_identifiability <- ggplot(p_symp_time_correlation, aes(week_start, correlation, colour = pair)) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_line(linewidth = 0.6) +
    scale_colour_brewer(
      palette = "Dark2",
      labels = c(
        p_symp_vs_reporting = "vs symptomatic-case reporting",
        p_symp_vs_detection = "vs overall detection"
      )
    ) +
    labs(
      title = "v2.1: p_symp identifiability over time",
      subtitle = "Posterior correlation between p_symp and time-varying reporting quantities",
      x = NULL, y = "Correlation", colour = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")
  diagnostic_plot_data <- sampler_diag |>
    dplyr::select(version, max_rhat, min_ess, divergences, max_treedepth_hits, min_bfmi) |>
    tidyr::pivot_longer(-version, names_to = "metric", values_to = "value")
  p_sampler_diag <- ggplot(diagnostic_plot_data, aes(version, value, fill = version)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = format(round(value, 3), trim = TRUE)), vjust = -0.25, size = 2.8) +
    facet_wrap(~metric, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = COLORS) +
    labs(title = "Sampler diagnostics", subtitle = "Full per-parameter Rhat and ESS are saved as a table", x = NULL, y = NULL) +
    theme_classic(base_size = 10)

  correlation_data <- dplyr::bind_rows(lapply(models, function(model) {
    correlation_long(model$scale_draws, model$label)
  }))
  p_correlations <- ggplot(correlation_data, aes(parameter_x, parameter_y, fill = correlation)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.5) +
    facet_wrap(~version, scales = "free") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1)) +
    labs(title = "Posterior correlation of main scale parameters", x = NULL, y = NULL, fill = "Correlation") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

  figure_dir <- here::here("03_Output/figures/renewal_v2_1/compare_v2_vs_v2_1")
  table_dir <- here::here("03_Output/tables/renewal_v2_1/compare_v2_vs_v2_1")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  tag <- Sys.getenv("RENEWAL_V2_COMPARE_TAG", unset = "")
  suffix <- if (nzchar(tag)) paste0("_", tag) else ""
  output_base <- paste0("renewal_v2_vs_v2_1", suffix)

  save_plot <- function(name, plot, width, height) {
    path <- file.path(figure_dir, paste0(output_base, "_", name, ".png"))
    ggplot2::ggsave(path, plot, width = width, height = height, units = "in", dpi = 220, bg = "white")
    message("[save] ", path)
  }
  save_plot("cases", p_cases, 9, 7)
  save_plot("states", p_states, 11, 4.5)
  save_plot("ascertainment", p_case_ascertainment, 12, 4.5)
  save_plot("reproduction", p_reproduction, 10, 7)
  save_plot("serology_geography", p_serology_geo, 11, 10)
  save_plot("sampler_diagnostics", p_sampler_diag, 8, 6)
  save_plot("posterior_correlations", p_correlations, 11, 6)
  save_plot("ppc_residuals", p_ppc_residuals, 10, 6)
  save_plot("v2_1_p_symp_identifiability", p_p_symp_identifiability, 9, 4)

  utils::write.csv(trajectories, file.path(table_dir, paste0(output_base, "_trajectory_summary.csv")), row.names = FALSE)
  utils::write.csv(serology_observed, file.path(table_dir, paste0(output_base, "_serology_observed.csv")), row.names = FALSE)
  utils::write.csv(sampler_diag, file.path(table_dir, paste0(output_base, "_sampler_diagnostics.csv")), row.names = FALSE)
  utils::write.csv(parameter_diag, file.path(table_dir, paste0(output_base, "_parameter_diagnostics.csv")), row.names = FALSE)
  utils::write.csv(correlation_data, file.path(table_dir, paste0(output_base, "_posterior_correlations.csv")), row.names = FALSE)
  utils::write.csv(residuals_data, file.path(table_dir, paste0(output_base, "_ppc_residuals.csv")), row.names = FALSE)
  utils::write.csv(coverage_data, file.path(table_dir, paste0(output_base, "_ppc_coverage.csv")), row.names = FALSE)
  message("[save] diagnostic tables in ", table_dir)
}
