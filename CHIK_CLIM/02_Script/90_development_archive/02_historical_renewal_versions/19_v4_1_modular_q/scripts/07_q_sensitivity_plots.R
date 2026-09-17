# Checks whether five q nodes are sufficient (task Section 13): plots
# S_2025, cumulative infections, R0_2017, and R0_2022 against q. If these
# relationships look smooth and approximately monotonic, five nodes are
# enough for the primary modular approximation -- do NOT automatically add
# more. Only flag for a possible later 7-node sensitivity run if strong
# nonlinear curvature is visible.

required_packages <- c("here", "dplyr", "readr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(patchwork) })

sens_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_q_sensitivity_plots <- function() {
  root <- sens_root()
  base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "19_v4_1_modular_q")
  scalars_path <- file.path(base_dir, "combined_modular_ensemble_scalars.csv")
  if (!file.exists(scalars_path)) stop("Run 05_combine_modular_posterior.R first: ", scalars_path, " not found.")
  scalars <- read_csv(scalars_path, show_col_types = FALSE)

  node_medians <- scalars |> group_by(node_id, q_value) |> summarise(
    S_2025 = median(S_2025), cumulative_infections = median(cumulative_infections),
    R0_2017 = median(R0_2017), R0_2022 = median(R0_2022), .groups = "drop"
  ) |> arrange(q_value)

  theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))
  make_panel <- function(y, ylab, title) {
    ggplot(node_medians, aes(q_value, .data[[y]])) +
      geom_point(size = 2, colour = "#0072B2") + geom_line(colour = "#0072B2", linewidth = .5) +
      labs(title = title, x = "q (node value)", y = ylab) + theme_v4
  }
  p1 <- make_panel("S_2025", "S_2025 (proportion)", "A. S_2025 vs q")
  p2 <- make_panel("cumulative_infections", "Cumulative infections 2015-2025", "B. Cumulative infections vs q")
  p3 <- make_panel("R0_2017", "R0_2017", "C. R0_2017 vs q")
  p4 <- make_panel("R0_2022", "R0_2022", "D. R0_2022 vs q")
  figure <- (p1 | p2) / (p3 | p4) + plot_annotation(title = "v4.1 modular-q sensitivity: key outputs vs q node value")

  dir.create(base_dir, showWarnings = FALSE)
  ggsave(file.path(base_dir, "q_sensitivity_plots.png"), figure, width = 200, height = 180, units = "mm", dpi = 300)

  # Simple smoothness/monotonicity check: sign changes in successive
  # differences (a monotonic relationship has zero sign changes).
  check_monotonic <- function(y) {
    d <- diff(node_medians[[y]])
    sign_changes <- sum(diff(sign(d)) != 0)
    list(sign_changes = sign_changes, monotonic = sign_changes == 0)
  }
  checks <- lapply(c("S_2025", "cumulative_infections", "R0_2017", "R0_2022"), check_monotonic)
  names(checks) <- c("S_2025", "cumulative_infections", "R0_2017", "R0_2022")
  for (nm in names(checks)) {
    message(sprintf("[q-sensitivity] %s: %s (%d sign change(s) across the 5 nodes)", nm, if (checks[[nm]]$monotonic) "monotonic" else "NOT strictly monotonic", checks[[nm]]$sign_changes))
  }
  all_monotonic <- all(vapply(checks, function(x) x$monotonic, logical(1)))
  if (all_monotonic) {
    message("[q-sensitivity] All four outputs are smooth/monotonic across the 5 nodes -- five nodes are sufficient for the primary modular approximation. Do not add more nodes.")
  } else {
    message("[q-sensitivity] At least one output shows non-monotonic behaviour across the 5 nodes -- consider a later 7-node sensitivity analysis (do not do this automatically as part of this run).")
  }
  invisible(list(node_medians = node_medians, checks = checks))
}

if (sys.nframe() == 0L) run_q_sensitivity_plots()
