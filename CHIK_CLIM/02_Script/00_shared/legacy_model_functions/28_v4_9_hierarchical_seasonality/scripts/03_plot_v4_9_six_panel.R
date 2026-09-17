# Six-panel trajectory figure for a single v4.9 fixed-q fit, same format
# as 27_v4_8.../03_plot_v4_8_six_panel.R (2022 seed window shaded grey in
# panels A/B, not scored as predictive). Only addition relative to v4.8's
# script: panel D's subtitle notes the hierarchical seasonal-amplitude
# structure so the changed R0(t) mechanism is visible in the figure itself.
#
# Usage: QVAL=0.10 Rscript 03_plot_v4_9_six_panel.R

required_packages <- c("here", "rstan", "ggplot2", "dplyr", "patchwork", "scales", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(ggplot2); library(dplyr); library(patchwork) })

COL_INK <- "#202124"; COL_MUTED <- "#6B7280"; COL_GRID <- "#E5E7EB"
COL_NAVY <- "#315A7D"; COL_BLUE <- "#56B4E9"; COL_VERMILLION <- "#D55E00"
COL_ORANGE <- "#E69F00"; COL_GREEN <- "#009E73"; COL_PURPLE <- "#76558F"; COL_SEED <- "grey40"

theme_publication <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, color = COL_MUTED, margin = margin(b = 7)),
      axis.title = element_text(size = base_size, color = COL_INK),
      axis.text = element_text(size = base_size - 0.5, color = COL_MUTED),
      axis.line = element_line(color = COL_INK, linewidth = 0.35),
      axis.ticks = element_line(color = COL_INK, linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "top", legend.justification = "left", legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.width = grid::unit(13, "pt"), legend.key.height = grid::unit(7, "pt"),
      plot.margin = margin(7, 8, 7, 7)
    )
}

summarise_trajectory <- function(draw_matrix, prefix) {
  stopifnot(length(dim(draw_matrix)) == 2L)
  tibble::tibble(
    variable = prefix, t = seq_len(ncol(draw_matrix)),
    q025 = apply(draw_matrix, 2, stats::quantile, probs = 0.025),
    q25 = apply(draw_matrix, 2, stats::quantile, probs = 0.25),
    median = apply(draw_matrix, 2, stats::median),
    q75 = apply(draw_matrix, 2, stats::quantile, probs = 0.75),
    q975 = apply(draw_matrix, 2, stats::quantile, probs = 0.975)
  )
}

add_interval_layers <- function(plot, data, colour, fill = colour) {
  plot +
    geom_ribbon(data = data, aes(ymin = q025, ymax = q975), fill = fill, alpha = 0.13, colour = NA) +
    geom_ribbon(data = data, aes(ymin = q25, ymax = q75), fill = fill, alpha = 0.24, colour = NA) +
    geom_line(data = data, aes(y = median), colour = colour, linewidth = 0.65)
}

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "03_Output/02_ceara_pipeline/model_fits/baseline/v4_9_hierarchical_seasonality")

plot_v4_9_six_panel <- function(qval) {
  tag <- paste0("q", sprintf("%.2f", qval))
  bundle <- readRDS(file.path(base_dir, "outputs", tag, paste0("renewal_ceara_v4_9_fit_", tag, ".rds")))
  fit <- bundle$fit
  weekly <- bundle$weekly_data
  dates <- as.Date(weekly$week_start)
  seed_xmin <- min(bundle$seed_dates); seed_xmax <- max(bundle$seed_dates) + 7

  pars <- c("X", "S_prop", "immune_prop", "R0_t", "R_eff_t", "expected_reported_cases", "C_pred")
  draws <- rstan::extract(fit, pars = pars, permuted = TRUE)

  summaries <- dplyr::bind_rows(
    summarise_trajectory(draws$X, "latent_infections"),
    summarise_trajectory(draws$S_prop, "susceptible_prop"),
    summarise_trajectory(draws$immune_prop, "immune_prop"),
    summarise_trajectory(draws$R0_t, "R0_t"),
    summarise_trajectory(draws$R_eff_t, "R_eff_t"),
    summarise_trajectory(draws$expected_reported_cases, "expected_reported_cases"),
    summarise_trajectory(draws$C_pred, "posterior_predictive_cases")
  ) |> dplyr::mutate(week_start = dates[t])

  get_summary <- function(name) dplyr::filter(summaries, variable == name)
  date_breaks <- seq(as.Date("2015-01-01"), as.Date("2026-01-01"), by = "1 year")
  x_scale <- scale_x_date(breaks = date_breaks, date_labels = "%Y", expand = expansion(mult = c(0.01, 0.015)))
  seed_rect <- function(p) p + annotate("rect", xmin = seed_xmin, xmax = seed_xmax, ymin = -Inf, ymax = Inf, fill = COL_SEED, alpha = 0.18)

  expected <- get_summary("expected_reported_cases")
  predictive <- get_summary("posterior_predictive_cases")
  cases_df <- tibble::tibble(week_start = dates, observed = weekly$cases)

  p_cases <- seed_rect(ggplot()) +
    geom_ribbon(data = predictive, aes(week_start, ymin = q025, ymax = q975), fill = COL_BLUE, alpha = 0.15) +
    geom_ribbon(data = expected, aes(week_start, ymin = q25, ymax = q75), fill = COL_NAVY, alpha = 0.24) +
    geom_line(data = expected, aes(week_start, median, colour = "Model expectation"), linewidth = 0.65) +
    geom_point(data = cases_df, aes(week_start, observed, colour = "Observed cases"), size = 0.5, alpha = 0.65) +
    scale_colour_manual(values = c("Observed cases" = COL_INK, "Model expectation" = COL_NAVY)) +
    scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    x_scale +
    labs(title = "Reported cases and posterior prediction",
         subtitle = "Grey band = 2022 seed window (conditioned, excluded from likelihood -- not predictive)",
         x = NULL, y = "Weekly reported cases") +
    theme_publication()

  latent <- get_summary("latent_infections")
  p_infections <- seed_rect(add_interval_layers(ggplot(latent, aes(week_start)), latent, COL_VERMILLION)) +
    x_scale + scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
    labs(title = "Latent infection incidence", subtitle = "Grey band = 2022 seed window (X = C/q, conditioned)",
         x = NULL, y = "True infections per week") +
    theme_publication()

  immunity <- dplyr::bind_rows(
    get_summary("susceptible_prop") |> mutate(state = "Susceptible"),
    get_summary("immune_prop") |> mutate(state = "Infection-derived immune")
  )
  p_immunity <- ggplot(immunity, aes(week_start, colour = state, fill = state)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.7) +
    scale_colour_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune" = COL_PURPLE)) +
    scale_fill_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune" = COL_PURPLE)) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))) +
    x_scale +
    labs(title = "Population susceptibility and immunity", subtitle = "Closed population; LIFELONG immunity; continuous through the 2022 seed",
         x = NULL, y = "Population proportion") +
    theme_publication()

  reproduction <- dplyr::bind_rows(
    get_summary("R0_t") |> mutate(number = "R0(t)"),
    get_summary("R_eff_t") |> mutate(number = "Reff(t)")
  )
  p_reproduction <- ggplot(reproduction, aes(week_start, colour = number, fill = number)) +
    geom_hline(yintercept = 1, colour = COL_MUTED, linewidth = 0.35, linetype = "22") +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.09, colour = NA) +
    geom_line(aes(y = median), linewidth = 0.65) +
    scale_colour_manual(values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION),
                         labels = c("R0(t)" = expression(R[0](t)), "Reff(t)" = expression(R[eff](t)))) +
    scale_fill_manual(values = c("R0(t)" = COL_NAVY, "Reff(t)" = COL_VERMILLION),
                       labels = c("R0(t)" = expression(R[0](t)), "Reff(t)" = expression(R[eff](t)))) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
    x_scale +
    labs(title = "Time-varying reproduction numbers", subtitle = "Hierarchical seasonal amplitude: A_year modulates the first harmonic per year",
         x = NULL, y = "Reproduction number") +
    theme_publication()

  ascertainment_df <- tibble::tibble(week_start = dates, q = qval)
  p_ascertainment <- ggplot(ascertainment_df, aes(week_start, q)) +
    geom_line(colour = COL_ORANGE, linewidth = 0.9) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, max(0.35, qval * 1.3)), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(title = "Case ascertainment", subtitle = sprintf("q FIXED at %.2f (not estimated -- v4.9 hierarchical-seasonality sweep)", qval),
         x = NULL, y = "Ascertainment fraction (q)") +
    theme_publication()

  attack <- get_summary("immune_prop")
  sero_date <- as.Date(bundle$sero_week$matched_week_start)
  sero_observed <- 103 / 404
  sero_interval <- stats::binom.test(103, 404)$conf.int
  sero_observed_df <- tibble::tibble(week_start = sero_date, estimate = sero_observed, lower = sero_interval[1], upper = sero_interval[2])
  sero_idx <- bundle$sero_week$index
  model_at_sero <- tibble::tibble(
    week_start = sero_date, median = median(draws$immune_prop[, sero_idx]),
    q025 = quantile(draws$immune_prop[, sero_idx], 0.025), q975 = quantile(draws$immune_prop[, sero_idx], 0.975)
  )

  p_attack <- ggplot(attack, aes(week_start)) +
    annotate("rect", xmin = as.Date("2018-06-01"), xmax = as.Date("2018-12-31"), ymin = -Inf, ymax = Inf, fill = COL_ORANGE, alpha = 0.055) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_GREEN, alpha = 0.14) +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL_GREEN, alpha = 0.24) +
    geom_line(aes(y = median, colour = "Ceara immune fraction"), linewidth = 0.7) +
    geom_errorbar(data = sero_observed_df, aes(week_start, ymin = lower, ymax = upper, colour = "Observed Juazeiro seroprevalence"), width = 30, linewidth = 0.55) +
    geom_point(data = sero_observed_df, aes(week_start, estimate, colour = "Observed Juazeiro seroprevalence"), shape = 18, size = 2.3) +
    geom_errorbar(data = model_at_sero, aes(week_start, ymin = q025, ymax = q975, colour = "Model Ceara immune fraction (at survey)"), width = 30, linewidth = 0.55, position = position_nudge(x = 40)) +
    geom_point(data = model_at_sero, aes(week_start, median, colour = "Model Ceara immune fraction (at survey)"), shape = 16, size = 1.8, position = position_nudge(x = 40)) +
    scale_colour_manual(values = c(
      "Ceara immune fraction" = COL_GREEN,
      "Observed Juazeiro seroprevalence" = COL_ORANGE,
      "Model Ceara immune fraction (at survey)" = COL_PURPLE
    )) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
    x_scale +
    labs(title = "Cumulative infection and serological anchor", subtitle = "Yellow field: June-December 2018 survey",
         x = NULL, y = "Cumulative proportion infected") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_publication(base_size = 8) + theme(legend.text = element_text(size = 7.2))

  figure <- (p_cases | p_infections) / (p_immunity | p_reproduction) / (p_ascertainment | p_attack) +
    patchwork::plot_annotation(
      title = "Chikungunya transmission and population susceptibility in Cear\u00e1, Brazil (2015-2025)",
      subtitle = sprintf("v4.9 hierarchical-seasonality sweep, q = %.2f (%s) | 2022 seed conditioned, year-varying seasonal amplitude",
                          qval, if (bundle$hmc$hmc_pass) "HMC PASS" else "HMC FAIL/marginal"),
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 12.5, face = "bold", margin = margin(b = 3)),
        plot.subtitle = element_text(size = 9.5, colour = COL_MUTED, margin = margin(b = 8)),
        plot.tag = element_text(size = 11, face = "bold")
      )
    )

  figure_dir <- file.path(root, "03_Output/02_ceara_pipeline/figures/transmission/renewal_v4_9_hierarchical_seasonality")
  table_dir <- file.path(root, "03_Output/02_ceara_pipeline/tables/renewal_v4_9_hierarchical_seasonality")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  output_base <- paste0("renewal_v4_9_trajectories_ceara_", tag)
  png_path <- file.path(figure_dir, paste0(output_base, ".png"))
  pdf_path <- file.path(figure_dir, paste0(output_base, ".pdf"))
  table_path <- file.path(table_dir, paste0("renewal_v4_9_trajectory_summary_", tag, ".csv"))

  ggsave(png_path, figure, width = 183, height = 225, units = "mm", dpi = 300, bg = "white")
  ggsave(pdf_path, figure, width = 183, height = 225, units = "mm", device = cairo_pdf)
  utils::write.csv(summaries, table_path, row.names = FALSE)
  message("[", tag, "] saved: ", png_path)
  invisible(figure)
}

if (sys.nframe() == 0L) {
  qval <- as.numeric(Sys.getenv("QVAL", "0.10"))
  plot_v4_9_six_panel(qval)
}
