# Interim nationwide early-Re diagnostics.
# This script is intentionally independent of the annual/weekly S branch so
# that the completed renewal results can be inspected while that batch runs.

required_re_plot_packages <- c("here", "readr", "dplyr", "ggplot2", "patchwork", "scales")
missing_re_plot_packages <- required_re_plot_packages[
  !vapply(required_re_plot_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_re_plot_packages)) stop("Missing package(s): ", paste(missing_re_plot_packages, collapse = ", "))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2); library(patchwork); library(scales) })

plot_re_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}

run_early_re_diagnostics <- function() {
  root <- plot_re_root()
  table_path <- file.path(root, "03_Output", "tables", "national_wave_analysis", "brazil_chik_major_wave_early_re.csv")
  figure_dir <- file.path(root, "03_Output", "figures", "national_wave_analysis", "early_re")
  if (!file.exists(table_path)) stop("Early-Re output is missing; run 01_fit_brazil_major_wave_early_re.R first.")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  waves <- read_csv(table_path, show_col_types = FALSE) |>
    mutate(onset_week = as.Date(onset_week, format = "%Y-%m-%d"),
           label = paste(state, format(onset_week, "%Y-%m"), sep = " ")) |>
    dplyr::filter(re6_fit_status == "fitted") |>
    arrange(Re_early_6_median)
  if (!nrow(waves)) stop("No fitted Re6 waves are available for plotting.")
  waves <- waves |> mutate(plot_id = row_number())

  p_re <- ggplot(waves, aes(plot_id, Re_early_6_median)) +
    geom_linerange(aes(ymin = Re_early_6_q025, ymax = Re_early_6_q975), colour = "#56B4E9", linewidth = .35) +
    geom_point(aes(colour = is_recurrence), size = 1.3) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    scale_colour_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"), labels = c("First major wave", "Recurrent"), name = NULL) +
    labs(title = "Nationwide six-week early Re", subtitle = "Raw weekly confirmed cases; 95% posterior intervals; ordered by median", x = "Major wave", y = expression(R[e])) +
    theme_classic(base_size = 9) + theme(legend.position = "bottom")

  p_qc <- ggplot(waves, aes(re6_bulk_ess, re6_tail_ess, colour = is_recurrence)) +
    geom_vline(xintercept = 400, linetype = 2, colour = "grey50") +
    geom_hline(yintercept = 400, linetype = 2, colour = "grey50") +
    geom_point(size = 1.4, alpha = .8) +
    scale_colour_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"), guide = "none") +
    scale_x_log10(labels = label_number()) + scale_y_log10(labels = label_number()) +
    labs(title = "Wave-specific MCMC ESS", subtitle = "Dashed lines: pre-specified 400 threshold", x = "Bulk ESS", y = "Tail ESS") +
    theme_classic(base_size = 9)

  p_ppc <- ggplot(waves, aes(re6_ppc_all_weeks_95_coverage)) +
    geom_histogram(binwidth = .1, boundary = 0, fill = "#009E73", colour = "white") +
    scale_x_continuous(limits = c(0, 1), labels = label_percent()) +
    labs(title = "Posterior-predictive interval coverage", subtitle = "Fraction of the six observed weeks within its 95% posterior-predictive interval", x = "Coverage", y = "Number of waves") +
    theme_classic(base_size = 9)

  p_rhat <- ggplot(waves, aes(re6_rhat)) +
    geom_histogram(binwidth = .001, fill = "#CC79A7", colour = "white") +
    geom_vline(xintercept = 1.01, linetype = 2, colour = "grey30") +
    labs(title = "Wave-specific Rhat", subtitle = "Dashed line: pre-specified 1.01 threshold", x = "Rhat", y = "Number of waves") +
    theme_classic(base_size = 9)

  plot <- (p_re / (p_qc | p_ppc) / p_rhat) + plot_annotation(
    caption = "Shared-phi K-episode renewal model; primary estimates are used only when the per-wave HMC gate passes."
  )
  output_path <- file.path(figure_dir, "brazil_chik_major_wave_early_re_diagnostics.pdf")
  ggsave(output_path, plot, width = 240, height = 250, units = "mm", device = cairo_pdf)
  message("[save] ", output_path)
  invisible(output_path)
}

if (sys.nframe() == 0L) run_early_re_diagnostics()
