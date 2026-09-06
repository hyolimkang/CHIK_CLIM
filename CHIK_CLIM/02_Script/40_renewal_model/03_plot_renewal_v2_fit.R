# =============================================================================
# 03_plot_renewal_v2_fit.R
#
# Publication-quality trajectory figure for the Ceara renewal v2 model.
# Produces a six-panel figure and a tidy posterior trajectory summary.
# =============================================================================

required_packages <- c("here", "rstan", "ggplot2", "dplyr", "patchwork", "scales", "tibble")
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
  library(patchwork)
})

COL_INK <- "#202124"
COL_MUTED <- "#6B7280"
COL_GRID <- "#E5E7EB"
COL_NAVY <- "#315A7D"
COL_BLUE <- "#56B4E9"
COL_VERMILLION <- "#D55E00"
COL_ORANGE <- "#E69F00"
COL_GREEN <- "#009E73"
COL_PURPLE <- "#76558F"

theme_publication <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_family = "Arial", base_size = base_size) +
    ggplot2::theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, color = COL_MUTED, margin = margin(b = 7)),
      axis.title = element_text(size = base_size, color = COL_INK),
      axis.text = element_text(size = base_size - 0.5, color = COL_MUTED),
      axis.line = element_line(color = COL_INK, linewidth = 0.35),
      axis.ticks = element_line(color = COL_INK, linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.width = grid::unit(13, "pt"),
      legend.key.height = grid::unit(7, "pt"),
      plot.margin = margin(7, 8, 7, 7)
    )
}

summarise_trajectory <- function(draw_matrix, prefix) {
  stopifnot(length(dim(draw_matrix)) == 2L)
  tibble::tibble(
    variable = prefix,
    t = seq_len(ncol(draw_matrix)),
    q025 = apply(draw_matrix, 2, stats::quantile, probs = 0.025),
    q25 = apply(draw_matrix, 2, stats::quantile, probs = 0.25),
    median = apply(draw_matrix, 2, stats::median),
    q75 = apply(draw_matrix, 2, stats::quantile, probs = 0.75),
    q975 = apply(draw_matrix, 2, stats::quantile, probs = 0.975)
  )
}

add_interval_layers <- function(plot, data, colour, fill = colour) {
  plot +
    geom_ribbon(
      data = data, aes(ymin = q025, ymax = q975),
      fill = fill, alpha = 0.13, colour = NA
    ) +
    geom_ribbon(
      data = data, aes(ymin = q25, ymax = q75),
      fill = fill, alpha = 0.24, colour = NA
    ) +
    geom_line(data = data, aes(y = median), colour = colour, linewidth = 0.65)
}

if (sys.nframe() == 0) {
  FIT_TAG <- Sys.getenv("RENEWAL_V2_TAG", unset = "")
  suffix <- if (nzchar(FIT_TAG)) paste0("_", FIT_TAG) else ""
  fit_path <- here::here(
    "02_Script/stan", paste0("renewal_ceara_v2_fit", suffix, ".rds")
  )
  if (!file.exists(fit_path)) stop("Fit bundle not found: ", fit_path)

  bundle <- readRDS(fit_path)
  fit <- bundle$fit
  weekly <- bundle$weekly_data
  dates <- weekly$week_start

  pars <- c(
    "X", "S_prop", "immune_prop", "R0_t", "R_eff_t",
    "rho_sym_t", "rho_total_t", "expected_reported_cases", "C_pred",
    "site_attack_at_serosurvey"
  )
  draws <- rstan::extract(fit, pars = pars, permuted = TRUE)

  summaries <- dplyr::bind_rows(
    summarise_trajectory(draws$X, "latent_infections"),
    summarise_trajectory(draws$S_prop, "susceptible_prop"),
    summarise_trajectory(draws$immune_prop, "immune_prop"),
    summarise_trajectory(draws$R0_t, "R0_t"),
    summarise_trajectory(draws$R_eff_t, "R_eff_t"),
    summarise_trajectory(draws$rho_sym_t, "rho_sym_t"),
    summarise_trajectory(draws$rho_total_t, "rho_total_t"),
    summarise_trajectory(draws$expected_reported_cases, "expected_reported_cases"),
    summarise_trajectory(draws$C_pred, "posterior_predictive_cases")
  ) |>
    dplyr::mutate(week_start = dates[t])

  get_summary <- function(name) dplyr::filter(summaries, variable == name)
  date_breaks <- seq(as.Date("2015-01-01"), as.Date("2019-01-01"), by = "1 year")
  x_scale <- scale_x_date(
    breaks = date_breaks,
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.015))
  )

  expected <- get_summary("expected_reported_cases")
  predictive <- get_summary("posterior_predictive_cases")
  cases_df <- tibble::tibble(week_start = dates, observed = weekly$cases)

  p_cases <- ggplot() +
    geom_ribbon(
      data = predictive,
      aes(week_start, ymin = q025, ymax = q975),
      fill = COL_BLUE, alpha = 0.15
    ) +
    geom_ribbon(
      data = expected,
      aes(week_start, ymin = q25, ymax = q75),
      fill = COL_NAVY, alpha = 0.24
    ) +
    geom_line(
      data = expected, aes(week_start, median, colour = "Model expectation"),
      linewidth = 0.65
    ) +
    geom_point(
      data = cases_df, aes(week_start, observed, colour = "Observed cases"),
      size = 0.55, alpha = 0.72
    ) +
    scale_colour_manual(values = c(
      "Observed cases" = COL_INK,
      "Model expectation" = COL_NAVY
    )) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    x_scale +
    labs(
      title = "Reported cases and posterior prediction",
      subtitle = "Shading: 95% predictive interval; dark band: 50% CrI",
      x = NULL, y = "Weekly reported cases"
    ) +
    theme_publication()

  latent <- get_summary("latent_infections")
  p_infections <- add_interval_layers(
    ggplot(latent, aes(week_start)), latent, COL_VERMILLION
  ) +
    x_scale +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    labs(
      title = "Latent infection incidence",
      subtitle = "Median, 50% and 95% credible intervals",
      x = NULL, y = "True infections per week"
    ) +
    theme_publication()

  immunity <- dplyr::bind_rows(
    get_summary("susceptible_prop") |> mutate(state = "Susceptible"),
    get_summary("immune_prop") |> mutate(state = "Infection-derived immune")
  )
  p_immunity <- ggplot(immunity, aes(week_start, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.7) +
    scale_colour_manual(values = c(
      "Susceptible" = COL_GREEN,
      "Infection-derived immune" = COL_PURPLE
    )) +
    scale_fill_manual(values = c(
      "Susceptible" = COL_GREEN,
      "Infection-derived immune" = COL_PURPLE
    )) +
    scale_y_continuous(
      limits = c(0, 1), labels = scales::label_percent(accuracy = 1),
      breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))
    ) +
    x_scale +
    labs(
      title = "Population susceptibility and immunity",
      subtitle = "Closed population; lifelong immunity",
      x = NULL, y = "Population proportion"
    ) +
    theme_publication()

  reproduction <- dplyr::bind_rows(
    get_summary("R0_t") |>
      filter(t > bundle$stan_data$seed_weeks) |>
      mutate(number = "R0(t)"),
    get_summary("R_eff_t") |>
      filter(t > bundle$stan_data$seed_weeks) |>
      mutate(number = "Reff(t)")
  )
  p_reproduction <- ggplot(reproduction, aes(week_start, colour = number, fill = number)) +
    geom_hline(yintercept = 1, colour = COL_MUTED, linewidth = 0.35, linetype = "22") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.09, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.65) +
    scale_colour_manual(
      values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION),
      labels = c("R0(t)" = expression(R[0](t)), "Reff(t)" = expression(R[eff](t)))
    ) +
    scale_fill_manual(
      values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION),
      labels = c("R0(t)" = expression(R[0](t)), "Reff(t)" = expression(R[eff](t)))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    x_scale +
    labs(
      title = "Time-varying reproduction numbers",
      subtitle = expression(R[eff](t) == R[0](t) %*% S(t-1)/N),
      x = NULL, y = "Reproduction number"
    ) +
    theme_publication()

  reporting <- dplyr::bind_rows(
    get_summary("rho_sym_t") |> mutate(probability = "Symptomatic-case reporting"),
    get_summary("rho_total_t") |> mutate(probability = "Infection-to-report")
  )
  p_reporting <- ggplot(reporting, aes(week_start, colour = probability, fill = probability)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.7) +
    scale_colour_manual(values = c(
      "Symptomatic-case reporting" = COL_PURPLE,
      "Infection-to-report" = COL_ORANGE
    )) +
    scale_fill_manual(values = c(
      "Symptomatic-case reporting" = COL_PURPLE,
      "Infection-to-report" = COL_ORANGE
    )) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(
      title = "Case ascertainment",
      subtitle = sprintf("External symptomatic proportion fixed at %.2f", bundle$config$p_symp),
      x = NULL, y = "Probability"
    ) +
    theme_publication()

  attack <- get_summary("immune_prop")
  sero_date <- as.Date(bundle$config$sero_date)
  sero_observed <- 103 / 404
  sero_interval <- stats::binom.test(103, 404)$conf.int
  site_draws <- draws$site_attack_at_serosurvey
  site_summary <- tibble::tibble(
    week_start = sero_date,
    median = median(site_draws),
    q025 = quantile(site_draws, 0.025),
    q975 = quantile(site_draws, 0.975)
  )
  sero_observed_df <- tibble::tibble(
    week_start = sero_date,
    estimate = sero_observed,
    lower = sero_interval[1],
    upper = sero_interval[2]
  )

  p_attack <- ggplot(attack, aes(week_start)) +
    annotate(
      "rect", xmin = as.Date("2018-06-01"), xmax = as.Date("2018-12-31"),
      ymin = -Inf, ymax = Inf, fill = COL_ORANGE, alpha = 0.055
    ) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_GREEN, alpha = 0.14) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL_GREEN, alpha = 0.24) +
    geom_line(aes(y = median, colour = "Ceara attack rate"), linewidth = 0.7) +
    geom_errorbar(
      data = sero_observed_df,
      aes(week_start, ymin = lower, ymax = upper, colour = "Observed Juazeiro seroprevalence"),
      width = 18, linewidth = 0.55
    ) +
    geom_point(
      data = sero_observed_df,
      aes(week_start, estimate, colour = "Observed Juazeiro seroprevalence"),
      shape = 18, size = 2.3
    ) +
    geom_errorbar(
      data = site_summary,
      aes(week_start, ymin = q025, ymax = q975, colour = "Model-implied Juazeiro attack rate"),
      width = 18, linewidth = 0.55, position = position_nudge(x = 25)
    ) +
    geom_point(
      data = site_summary,
      aes(week_start, median, colour = "Model-implied Juazeiro attack rate"),
      shape = 16, size = 1.8, position = position_nudge(x = 25)
    ) +
    scale_colour_manual(values = c(
      "Ceara attack rate" = COL_GREEN,
      "Observed Juazeiro seroprevalence" = COL_ORANGE,
      "Model-implied Juazeiro attack rate" = COL_PURPLE
    ), labels = c(
      "Ceara attack rate" = "Ceara",
      "Observed Juazeiro seroprevalence" = "Juazeiro observed",
      "Model-implied Juazeiro attack rate" = "Juazeiro model"
    )) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(
      title = "Cumulative infection and serological anchor",
      subtitle = "Yellow field: June-December 2018 survey",
      x = NULL, y = "Cumulative proportion infected"
    ) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_publication(base_size = 8) +
    theme(legend.text = element_text(size = 7.2))

  figure <- (p_cases | p_infections) /
    (p_immunity | p_reproduction) /
    (p_reporting | p_attack) +
    patchwork::plot_annotation(
      title = "Chikungunya transmission and population susceptibility in Cear\u00e1, Brazil",
      subtitle = "Bayesian renewal model, weekly observations from 2015 to 2018",
      tag_levels = "A",
      theme = theme(
        text = element_text(family = "Arial", colour = COL_INK),
        plot.title = element_text(size = 12.5, face = "bold", margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9.5, colour = COL_MUTED, margin = margin(b = 8)),
        plot.tag = element_text(size = 11, face = "bold")
      )
    )

  figure_dir <- here::here("03_Output/figures")
  table_dir <- here::here("03_Output/tables")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  output_base <- paste0("renewal_v2_trajectories_ceara", suffix)
  pdf_path <- file.path(figure_dir, paste0(output_base, ".pdf"))
  tiff_path <- file.path(figure_dir, paste0(output_base, ".tiff"))
  png_path <- file.path(figure_dir, paste0(output_base, ".png"))
  table_path <- file.path(table_dir, paste0("renewal_v2_trajectory_summary", suffix, ".csv"))

  ggsave(pdf_path, figure, width = 183, height = 225, units = "mm", device = cairo_pdf)
  ggsave(tiff_path, figure, width = 183, height = 225, units = "mm", dpi = 600, compression = "lzw")
  ggsave(png_path, figure, width = 183, height = 225, units = "mm", dpi = 300, bg = "white")
  utils::write.csv(summaries, table_path, row.names = FALSE)

  message("[save] ", pdf_path)
  message("[save] ", tiff_path)
  message("[save] ", png_path)
  message("[save] ", table_path)
}
