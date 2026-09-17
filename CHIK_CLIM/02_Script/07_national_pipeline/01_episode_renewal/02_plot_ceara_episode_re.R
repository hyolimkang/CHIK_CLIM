# Diagnostics and figures for the Ceara episode-level early-phase renewal fit.
#
# This script is read-only: it never refits the model. It summarises each
# 4-, 6-, and 8-week sensitivity analysis, evaluates standard HMC diagnostics,
# and writes separate overview, early-phase posterior-predictive, and HMC PDFs.

required_packages <- c("here", "rstan", "posterior", "dplyr", "tibble", "tidyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(posterior)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

colour_re <- "#315A7D"
colour_cases <- "#D55E00"
colour_prediction <- "#56B4E9"

publication_theme <- function() {
  theme_classic(base_size = 8.5) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}

summarise_hmc <- function(fit, window_weeks) {
  parameters <- c("log_Re", "phi")
  draws <- posterior::as_draws_array(rstan::extract(
    fit, pars = parameters, permuted = FALSE, inc_warmup = FALSE
  ))
  parameter_summary <- posterior::summarise_draws(
    draws, posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  ) |>
    as.data.frame()
  names(parameter_summary) <- sub("^posterior::", "", names(parameter_summary))

  bfmi <- rstan::get_bfmi(fit)
  tibble(
    window_weeks = window_weeks,
    divergences = rstan::get_num_divergent(fit),
    max_treedepth_hits = rstan::get_num_max_treedepth(fit),
    maximum_rhat = max(parameter_summary$rhat, na.rm = TRUE),
    n_rhat_gt_1_01 = sum(parameter_summary$rhat > 1.01, na.rm = TRUE),
    minimum_bulk_ess = min(parameter_summary$ess_bulk, na.rm = TRUE),
    minimum_tail_ess = min(parameter_summary$ess_tail, na.rm = TRUE),
    minimum_bfmi = min(bfmi),
    gate_pass =
      divergences == 0 &&
      max_treedepth_hits == 0 &&
      maximum_rhat <= 1.01 &&
      minimum_bulk_ess >= 100 &&
      minimum_tail_ess >= 100 &&
      minimum_bfmi >= 0.3
  )
}

summarise_re <- function(fit_object, episodes, window_weeks) {
  re_draws <- rstan::extract(fit_object$fit, pars = "Re")$Re
  tibble(
    episode_id = episodes$episode_id,
    wave_id = episodes$wave_id,
    start_week = episodes$start_week,
    window_weeks = window_weeks,
    Re_q025 = apply(re_draws, 2, quantile, probs = 0.025),
    Re_median = apply(re_draws, 2, median),
    Re_q975 = apply(re_draws, 2, quantile, probs = 0.975),
    phi_q025 = quantile(rstan::extract(fit_object$fit, pars = "phi")$phi, 0.025),
    phi_median = median(rstan::extract(fit_object$fit, pars = "phi")$phi),
    phi_q975 = quantile(rstan::extract(fit_object$fit, pars = "phi")$phi, 0.975)
  )
}

summarise_ppc <- function(fit_object) {
  ppc_draws <- rstan::extract(fit_object$fit, pars = "I_rep")$I_rep
  episode_data <- fit_object$episode_data
  n_episode <- fit_object$stan_data$K
  n_week <- fit_object$stan_data$W

  bind_rows(lapply(seq_len(n_episode), function(episode_index) {
    tibble(
      episode_id = unique(episode_data$episode_id)[episode_index],
      week_in_episode = seq_len(n_week),
      I_rep_q025 = apply(ppc_draws[, episode_index, , drop = FALSE], 3, quantile, probs = 0.025),
      I_rep_median = apply(ppc_draws[, episode_index, , drop = FALSE], 3, median),
      I_rep_q975 = apply(ppc_draws[, episode_index, , drop = FALSE], 3, quantile, probs = 0.975)
    )
  })) |>
    left_join(
      episode_data |>
        select(episode_id, week_in_episode, week_start, I_obs, Lambda),
      by = c("episode_id", "week_in_episode")
    )
}

make_hmc_figure <- function(hmc_table) {
  hmc_long <- hmc_table |>
    select(
      window_weeks, divergences, max_treedepth_hits, maximum_rhat,
      minimum_bulk_ess, minimum_tail_ess, minimum_bfmi
    ) |>
    tidyr::pivot_longer(-window_weeks, names_to = "metric", values_to = "value")

  ggplot(hmc_long, aes(factor(window_weeks), value)) +
    geom_point(size = 2, colour = colour_re) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    labs(
      title = "Episode renewal model: HMC diagnostics",
      subtitle = "All diagnostic gates must pass before interpreting Re estimates.",
      x = "Early-phase fitting window (weeks)", y = NULL
    ) +
    publication_theme()
}

run_episode_renewal_diagnostics <- function() {
  fit_path <- Sys.getenv(
    "RENEWAL_EPISODE_FIT",
    here::here("02_Script", "stan", "renewal_ceara_episode_re_fit.rds")
  )
  if (!file.exists(fit_path)) stop("Fit bundle not found: ", fit_path)

  table_directory <- here::here(
    "03_Output", "tables", "renewal_episode_re", "ceara_early_phase_re"
  )
  figure_directory <- here::here(
    "03_Output", "figures", "renewal_episode_re", "ceara_early_phase_re"
  )
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

  bundle <- readRDS(fit_path)
  window_names <- names(bundle$fits)
  window_weeks <- as.integer(window_names)

  hmc_table <- bind_rows(lapply(window_names, function(window_name) {
    summarise_hmc(bundle$fits[[window_name]]$fit, as.integer(window_name))
  })) |>
    arrange(window_weeks)
  re_summary <- bind_rows(lapply(window_names, function(window_name) {
    summarise_re(bundle$fits[[window_name]], bundle$episodes, as.integer(window_name))
  })) |>
    arrange(start_week, window_weeks)

  primary_window <- bundle$config$primary_window
  primary_fit <- bundle$fits[[as.character(primary_window)]]
  ppc_summary <- summarise_ppc(primary_fit)

  write.csv(hmc_table, file.path(table_directory, "hmc_diagnostics.csv"), row.names = FALSE)
  write.csv(re_summary, file.path(table_directory, "episode_re_summary.csv"), row.names = FALSE)
  write.csv(ppc_summary, file.path(table_directory, "primary_window_ppc.csv"), row.names = FALSE)

  all_gates_pass <- all(hmc_table$gate_pass)
  interpretation_label <- if (all_gates_pass) {
    "All HMC gates pass"
  } else {
    "UNCONVERGED - do not interpret posterior Re estimates"
  }

  primary_windows <- primary_fit$episode_data |>
    group_by(episode_id) |>
    summarise(
      start_week = min(week_start),
      end_week = max(week_start) + 7,
      .groups = "drop"
    )
  state_series_plot <- ggplot(bundle$weekly_cases, aes(week_start, cases)) +
    geom_rect(
      data = primary_windows,
      aes(xmin = start_week, xmax = end_week, ymin = 0, ymax = Inf),
      inherit.aes = FALSE, fill = "#E69F00", alpha = 0.18
    ) +
    geom_line(colour = colour_cases, linewidth = 0.35) +
    # Pseudo-log retains zero weeks while showing intuitive case-count labels.
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      labels = scales::label_number(big.mark = ",")
    ) +
    labs(
      title = "Ceara reported chikungunya incidence",
      subtitle = "Gold bands: eight-week early-phase renewal windows; y-axis uses log(1 + cases).",
      x = NULL, y = "Reported weekly cases"
    ) +
    publication_theme()

  re_plot <- ggplot(re_summary, aes(Re_median, reorder(episode_id, start_week), colour = factor(window_weeks))) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = Re_q025, xmax = Re_q975), height = 0.15) +
    geom_point(size = 1.8) +
    scale_colour_manual(values = c("4" = "#999999", "6" = "#56B4E9", "8" = colour_re)) +
    labs(
      title = "Episode-level effective reproduction numbers",
      subtitle = interpretation_label,
      x = expression(R[e]), y = "Outbreak episode", colour = "Window"
    ) +
    publication_theme()

  overview <- state_series_plot / re_plot +
    plot_annotation(title = "Ceara early-phase renewal analysis")
  ggsave(
    file.path(figure_directory, "ceara_episode_re_overview.pdf"),
    overview, width = 183, height = 180, units = "mm", device = cairo_pdf
  )

  ppc_plot <- ggplot(ppc_summary, aes(week_in_episode, I_rep_median)) +
    geom_ribbon(aes(ymin = I_rep_q025, ymax = I_rep_q975), fill = colour_prediction, alpha = 0.25) +
    geom_line(colour = colour_prediction, linewidth = 0.55) +
    geom_point(aes(y = I_obs), colour = colour_cases, size = 1.1) +
    facet_wrap(~ episode_id, scales = "free_y", ncol = 2) +
    scale_x_continuous(breaks = seq_len(primary_window)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(4),
      labels = scales::label_number(big.mark = ",")
    ) +
    labs(
      title = "Eight-week early-phase posterior predictive checks",
      subtitle = interpretation_label,
      x = "Week since audit-defined outbreak onset",
      y = "Reported cases"
    ) +
    publication_theme()
  ggsave(
    file.path(figure_directory, "ceara_episode_re_primary_window_ppc.pdf"),
    ppc_plot, width = 183, height = 180, units = "mm", device = cairo_pdf
  )

  ggsave(
    file.path(figure_directory, "ceara_episode_re_hmc_diagnostics.pdf"),
    make_hmc_figure(hmc_table), width = 183, height = 145, units = "mm", device = cairo_pdf
  )

  message("[gate] ", interpretation_label)
}


if (sys.nframe() == 0L) {
  run_episode_renewal_diagnostics()
}
