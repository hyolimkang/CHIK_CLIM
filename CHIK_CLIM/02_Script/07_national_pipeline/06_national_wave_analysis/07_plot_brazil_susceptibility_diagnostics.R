# State-by-state retrospective susceptibility diagnostics.
# The weekly trajectory is explicitly retrospective/smoothed: annual posterior
# infection burden is allocated within each year using the 3-week case curve.

required_s_plot_packages <- c("here", "readr", "dplyr", "ggplot2", "patchwork", "scales")
missing_s_plot_packages <- required_s_plot_packages[
  !vapply(required_s_plot_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_s_plot_packages)) stop("Missing package(s): ", paste(missing_s_plot_packages, collapse = ", "))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(patchwork); library(scales) })

sus_plot_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}

run_susceptibility_diagnostics <- function() {
  root <- sus_plot_root()
  table_dir <- file.path(root, "03_Output", "tables", "national_wave_analysis")
  figure_dir <- file.path(root, "03_Output", "figures", "national_wave_analysis", "susceptibility")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  weekly <- read_csv(file.path(table_dir, "brazil_chik_weekly_susceptibility_summary.csv"), show_col_types = FALSE) |>
    mutate(week_start = as.Date(week_start))
  annual <- read_csv(file.path(table_dir, "brazil_chik_annual_susceptibility_summary.csv"), show_col_types = FALSE)
  q <- read_csv(file.path(table_dir, "brazil_chik_annual_shape_q_implied.csv"), show_col_types = FALSE)
  hmc <- read_csv(file.path(table_dir, "brazil_chik_annual_shape_hmc_gate.csv"), show_col_types = FALSE)
  boundary <- read_csv(file.path(table_dir, "brazil_chik_weekly_susceptibility_year_boundary_check.csv"), show_col_types = FALSE)
  states <- sort(unique(weekly$state))
  output_path <- file.path(figure_dir, "brazil_chik_weekly_susceptibility_diagnostics.pdf")
  grDevices::pdf(output_path, width = 10.5, height = 7.5, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (state_value in states) {
    w <- dplyr::filter(weekly, state == state_value)
    a <- dplyr::filter(annual, state == state_value)
    p_cases <- ggplot(w, aes(week_start, reported_cases)) +
      geom_line(colour = "grey25", linewidth = .28) +
      scale_y_continuous(trans = pseudo_log_trans(10), labels = label_number()) +
      labs(title = paste0(state_value, ": confirmed weekly cases"), x = NULL, y = "Cases") + theme_classic(base_size = 9)
    p_s <- ggplot(w, aes(week_start, median)) +
      geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#56B4E9", alpha = .25) +
      geom_line(colour = "#0072B2", linewidth = .45) +
      scale_y_continuous(limits = c(0, 1), labels = label_percent()) +
      labs(title = "Retrospective weekly susceptible fraction", subtitle = "Annual posterior infection burden, placed by 3-week-smoothed weekly case timing", x = NULL, y = "S / N") + theme_classic(base_size = 9)
    p_lambda <- ggplot(a, aes(year, lambda_median)) +
      geom_ribbon(aes(ymin = lambda_q025, ymax = lambda_q975), fill = "#E69F00", alpha = .25) +
      geom_line(colour = "#D55E00", linewidth = .5) + geom_point(colour = "#D55E00", size = 1.2) +
      scale_y_continuous(trans = "log10") + scale_x_continuous(breaks = seq(2015, 2025, 2)) +
      labs(title = "Annual FOI posterior", x = NULL, y = expression(lambda)) + theme_classic(base_size = 9)
    p_attack <- ggplot(a, aes(year, infections_median)) +
      geom_ribbon(aes(ymin = infections_q025, ymax = infections_q975), fill = "#009E73", alpha = .25) +
      geom_line(colour = "#009E73", linewidth = .5) + geom_point(colour = "#009E73", size = 1.2) +
      scale_y_continuous(trans = pseudo_log_trans(10), labels = label_number()) + scale_x_continuous(breaks = seq(2015, 2025, 2)) +
      labs(title = "Annual latent infections posterior", x = NULL, y = "Infections") + theme_classic(base_size = 9)
    print((p_cases | p_s) / (p_lambda | p_attack) + plot_annotation(title = paste0("UF ", state_value)))
  }
  p_q <- ggplot(q, aes(reorder(state, q_implied_median), q_implied_median)) +
    geom_linerange(aes(ymin = q_implied_q025, ymax = q_implied_q975), colour = "#0072B2") + geom_point(colour = "#0072B2") +
    coord_flip() + geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    labs(title = "Implied overall confirmed-case detection", subtitle = "Diagnostic only; not fed back into the SHAPE likelihood", x = NULL, y = expression(q[implied])) + theme_classic(base_size = 9)
  p_hmc <- ggplot(hmc, aes(min_bulk_ess, min_tail_ess, colour = hmc_pass)) +
    geom_vline(xintercept = 400, linetype = 2) + geom_hline(yintercept = 400, linetype = 2) + geom_point(size = 2) +
    geom_text(aes(label = state), nudge_y = 100, size = 2.5, show.legend = FALSE) +
    scale_x_log10(labels = label_number()) + scale_y_log10(labels = label_number()) +
    scale_colour_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#D55E00"), name = "Annual HMC pass") +
    labs(title = "Annual SHAPE HMC diagnostics", x = "Minimum bulk ESS", y = "Minimum tail ESS") + theme_classic(base_size = 9)
  p_boundary <- ggplot(boundary, aes(year, max_abs_discrepancy_prop, colour = state)) +
    geom_line(show.legend = FALSE) + scale_y_continuous(trans = "log10") +
    labs(title = "Annual-versus-weekly S end-of-year discrepancy", x = "Year", y = "Maximum absolute discrepancy / annual N") + theme_classic(base_size = 9)
  print((p_q | p_hmc) / p_boundary + plot_annotation(title = "Nationwide retrospective susceptibility diagnostics"))
  invisible(output_path)
}

if (sys.nframe() == 0L) run_susceptibility_diagnostics()
