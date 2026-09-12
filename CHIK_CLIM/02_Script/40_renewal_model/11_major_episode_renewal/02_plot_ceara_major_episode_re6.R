# Read-only diagnostics and figures for the fixed six-week major-outbreak fit.

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

colour_cases <- "#D55E00"
colour_prediction <- "#56B4E9"
colour_re <- "#315A7D"

publication_theme <- function() {
  theme_classic(base_size = 8.5) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      strip.background = element_rect(fill = "grey95", colour = NA)
    )
}

summarise_hmc <- function(fit) {
  sampler_draws <- posterior::as_draws_array(rstan::extract(
    fit, pars = c("log_Re", "phi"), permuted = FALSE, inc_warmup = FALSE
  ))
  parameter_summary <- posterior::summarise_draws(
    sampler_draws, posterior::rhat, posterior::ess_bulk, posterior::ess_tail
  ) |>
    as.data.frame()
  names(parameter_summary) <- sub("^posterior::", "", names(parameter_summary))

  bfmi <- rstan::get_bfmi(fit)
  tibble(
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

make_major_episode_table <- function(bundle) {
  re_draws <- rstan::extract(bundle$fit, pars = "Re")$Re
  episode_level_data <- bundle$episode_data |>
    group_by(episode_id, wave_id, onset, peak_week, peak_cases, total_episode_cases, duration_weeks) |>
    summarise(
      cases_first_4wk = first(cases_first_4wk),
      cases_first_8wk = first(cases_first_8wk),
      .groups = "drop"
    ) |>
    arrange(onset)

  episode_level_data |>
    mutate(
      Re6_q025 = apply(re_draws, 2, quantile, probs = 0.025),
      Re6_median = apply(re_draws, 2, median),
      Re6_q975 = apply(re_draws, 2, quantile, probs = 0.975)
    )
}

make_ppc_table <- function(bundle) {
  ppc_draws <- rstan::extract(bundle$fit, pars = "I_rep")$I_rep
  n_episode <- bundle$stan_data$K
  n_week <- bundle$stan_data$W
  episode_ids <- unique(bundle$episode_data$episode_id)

  bind_rows(lapply(seq_len(n_episode), function(episode_index) {
    tibble(
      episode_id = episode_ids[episode_index],
      week_in_episode = seq_len(n_week),
      I_rep_q025 = apply(ppc_draws[, episode_index, , drop = FALSE], 3, quantile, probs = 0.025),
      I_rep_median = apply(ppc_draws[, episode_index, , drop = FALSE], 3, median),
      I_rep_q975 = apply(ppc_draws[, episode_index, , drop = FALSE], 3, quantile, probs = 0.975)
    )
  })) |>
    left_join(
      bundle$episode_data |>
        select(episode_id, week_in_episode, week_start, I_obs, Lambda),
      by = c("episode_id", "week_in_episode")
    )
}

run_major_episode_diagnostics <- function() {
  fit_path <- Sys.getenv(
    "RENEWAL_MAJOR_EPISODE_FIT",
    here::here("02_Script", "stan", "renewal_ceara_major_episode_re6_fit.rds")
  )
  if (!file.exists(fit_path)) stop("Fit bundle not found: ", fit_path)

  table_directory <- here::here(
    "03_Output", "tables", "renewal_major_episode_re", "ceara_major_outbreaks_6wk"
  )
  figure_directory <- here::here(
    "03_Output", "figures", "renewal_major_episode_re", "ceara_major_outbreaks_6wk"
  )
  dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

  bundle <- readRDS(fit_path)
  hmc <- summarise_hmc(bundle$fit)
  episode_table <- make_major_episode_table(bundle)
  ppc_table <- make_ppc_table(bundle)

  write.csv(hmc, file.path(table_directory, "hmc_diagnostics.csv"), row.names = FALSE)
  write.csv(episode_table, file.path(table_directory, "major_episode_re6_summary.csv"), row.names = FALSE)
  write.csv(ppc_table, file.path(table_directory, "major_episode_re6_ppc.csv"), row.names = FALSE)

  label <- if (hmc$gate_pass) {
    "HMC diagnostic gate passed"
  } else {
    "UNCONVERGED - do not interpret posterior Re estimates"
  }

  window_rectangles <- bundle$episode_data |>
    group_by(episode_id) |>
    summarise(
      onset = min(week_start),
      end_week = max(week_start) + 7,
      .groups = "drop"
    )
  state_plot <- ggplot(bundle$weekly_cases, aes(week_start, cases)) +
    geom_rect(
      data = window_rectangles,
      aes(xmin = onset, xmax = end_week, ymin = 0, ymax = Inf),
      inherit.aes = FALSE, fill = "#E69F00", alpha = 0.18
    ) +
    geom_line(colour = colour_cases, linewidth = 0.35) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      breaks = c(0, 1, 10, 100, 1000, 10000),
      labels = scales::label_number(big.mark = ",")
    ) +
    labs(
      title = "Ceara reported chikungunya incidence",
      subtitle = "Gold bands: fixed six-week windows for major outbreaks (candidate total cases >= 4,000).",
      x = NULL, y = "Reported weekly cases"
    ) +
    publication_theme()

  re_plot <- ggplot(episode_table, aes(Re6_median, reorder(episode_id, onset))) +
    geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
    geom_errorbarh(aes(xmin = Re6_q025, xmax = Re6_q975), height = 0.15, colour = colour_re) +
    geom_point(colour = colour_re, size = 2) +
    labs(
      title = "Major-outbreak early-phase effective reproduction numbers",
      subtitle = label,
      x = expression(R[e]~"(first six weeks)"), y = "Major outbreak"
    ) +
    publication_theme()

  overview <- state_plot / re_plot +
    plot_annotation(title = "Ceara major-outbreak renewal analysis")
  ggsave(
    file.path(figure_directory, "ceara_major_episode_re6_overview.pdf"),
    overview, width = 183, height = 180, units = "mm", device = cairo_pdf
  )

  ppc_plot <- ggplot(ppc_table, aes(week_in_episode, I_rep_median)) +
    geom_ribbon(aes(ymin = I_rep_q025, ymax = I_rep_q975), fill = colour_prediction, alpha = 0.25) +
    geom_line(colour = colour_prediction, linewidth = 0.55) +
    geom_point(aes(y = I_obs), colour = colour_cases, size = 1.2) +
    facet_wrap(~ episode_id, scales = "free_y", ncol = 2) +
    scale_x_continuous(breaks = seq_len(bundle$config$fit_window_weeks)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(4),
      labels = scales::label_number(big.mark = ",")
    ) +
    labs(
      title = "Major-outbreak six-week posterior predictive checks",
      subtitle = label,
      x = "Week since audit-defined onset", y = "Reported cases"
    ) +
    publication_theme()
  ggsave(
    file.path(figure_directory, "ceara_major_episode_re6_ppc.pdf"),
    ppc_plot, width = 183, height = 145, units = "mm", device = cairo_pdf
  )

  message("[gate] ", label)
}


if (sys.nframe() == 0L) {
  run_major_episode_diagnostics()
}
