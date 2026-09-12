# Six-panel diagnostic figure for one completed q-node fit, matching the
# exact panel layout/style used for v4.0 (17_v4_0_minimal_no_vaccine/
# 02_diagnose_v4_0_minimal_no_vaccine.R): A weekly cases + posterior
# prediction, B latent infections, C susceptible/immune stocks, D R0/Reff,
# E annual reported-case validation, F external serology consistency. The
# subtitle reports q_external and the HMC gate result for this node.
#
# Usage: Q_NODE_ID=q50 Rscript 09_plot_q_node_figures.R

required_packages <- c("here", "rstan", "dplyr", "tibble", "ggplot2", "patchwork", "scales", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(ggplot2); library(patchwork); library(scales); library(readr) })

plot_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

summarise_matrix <- function(x, dates, variable) {
  tibble(
    week_start = dates, variable = variable,
    q025 = apply(x, 2, quantile, probs = .025), q25 = apply(x, 2, quantile, probs = .25),
    median = apply(x, 2, median), q75 = apply(x, 2, quantile, probs = .75), q975 = apply(x, 2, quantile, probs = .975)
  )
}

add_band <- function(plot, summary, colour) {
  plot +
    geom_ribbon(data = summary, aes(ymin = q025, ymax = q975), fill = colour, alpha = .13) +
    geom_ribbon(data = summary, aes(ymin = q25, ymax = q75), fill = colour, alpha = .24) +
    geom_line(data = summary, aes(y = median), colour = colour, linewidth = .65)
}

run_plot_q_node_figures <- function(node_id = Sys.getenv("Q_NODE_ID", "q50")) {
  root <- plot_root()
  base_dir <- file.path(root, "02_Script", "40_renewal_model", "19_v4_1_modular_q")
  node_out_dir <- file.path(base_dir, "outputs", node_id)
  fit_path <- file.path(node_out_dir, paste0("renewal_ceara_v4_1_modular_q_fit_", node_id, ".rds"))
  if (!file.exists(fit_path)) stop("No fit found for node '", node_id, "': ", fit_path)
  bundle <- readRDS(fit_path)
  fit <- bundle$fit
  dates <- as.Date(bundle$weekly_data$week_start)
  years <- as.integer(format(dates, "%Y"))

  hmc_path <- file.path(node_out_dir, paste0("hmc_diagnostics_", node_id, ".csv"))
  hmc_pass <- if (file.exists(hmc_path)) {
    h <- read_csv(hmc_path, show_col_types = FALSE)
    h$value[h$metric == "divergences"] == 0 && h$value[h$metric == "max_treedepth_hits"] == 0 &&
      h$value[h$metric == "maximum_Rhat"] <= 1.01 && h$value[h$metric == "minimum_bulk_ESS"] >= 100
  } else isTRUE(bundle$hmc$hmc_pass)

  figure_dir <- file.path(root, "03_Output", "figures", "renewal_v4_1_modular_q", node_id)
  table_dir <- file.path(root, "03_Output", "tables", "renewal_v4_1_modular_q", node_id)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  draws <- rstan::extract(fit, pars = c("C_pred", "X", "S_prop", "immune_prop", "R0_t", "R_eff_t", "expected_reported_cases"))
  summary <- bind_rows(
    summarise_matrix(draws$C_pred, dates, "posterior_predictive_cases"), summarise_matrix(draws$X, dates, "latent_infections"),
    summarise_matrix(draws$S_prop, dates, "susceptible_proportion"), summarise_matrix(draws$immune_prop, dates, "infection_derived_immune_proportion"),
    summarise_matrix(draws$R0_t, dates, "R0"), summarise_matrix(draws$R_eff_t, dates, "Reff"), summarise_matrix(draws$expected_reported_cases, dates, "expected_reported_cases")
  )
  write_csv(summary, file.path(table_dir, paste0("trajectory_summary_", node_id, ".csv")))
  get_summary <- function(variable_name) dplyr::filter(summary, .data$variable == variable_name)
  case_summary <- get_summary("posterior_predictive_cases")

  annual <- bind_rows(lapply(sort(unique(years)), function(year) {
    index <- which(years == year)
    infections <- rowSums(draws$X[, index, drop = FALSE]); expected_cases <- rowSums(draws$expected_reported_cases[, index, drop = FALSE])
    tibble(
      year = year, observed_cases = sum(bundle$weekly_data$cases[index]),
      expected_cases_q025 = quantile(expected_cases, .025), expected_cases_median = median(expected_cases), expected_cases_q975 = quantile(expected_cases, .975),
      infections_median = median(infections)
    )
  }))
  write_csv(annual, file.path(table_dir, paste0("annual_attack_decomposition_", node_id, ".csv")))

  serology_windows <- tibble(
    site = c("Juazeiro do Norte", "Quixada"), observed_positive = c(103L, 289L), observed_n = c(404L, 409L),
    start = as.Date(c("2018-06-03", "2018-06-03")), end = as.Date(c("2018-12-30", "2019-12-29"))
  ) |> rowwise() |> mutate(
    state_immune_q025 = quantile(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE]), .025),
    state_immune_median = median(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE])),
    state_immune_q975 = quantile(rowMeans(draws$immune_prop[, dates >= start & dates <= end, drop = FALSE]), .975),
    observed_proportion = observed_positive / observed_n
  ) |> ungroup()

  theme_v4 <- theme_classic(base_size = 8.5) + theme(legend.position = "top", legend.title = element_blank(), panel.grid.major.y = element_line(colour = "grey90"))
  xscale <- scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2026-01-01"), by = "2 years"), date_labels = "%Y")
  p1 <- add_band(ggplot(case_summary, aes(week_start)), case_summary, "#56B4E9") +
    geom_point(data = tibble(week_start = dates, observed_cases = bundle$weekly_data$cases), aes(week_start, observed_cases), size = .3) +
    xscale + labs(title = "A. Reported cases and posterior prediction", x = NULL, y = "Weekly cases") + theme_v4
  p2_data <- get_summary("latent_infections")
  p2 <- add_band(ggplot(p2_data, aes(week_start)), p2_data, "#D55E00") + xscale + labs(title = "B. Latent infection incidence", x = NULL, y = "Infections") + theme_v4
  p3_data <- bind_rows(mutate(get_summary("susceptible_proportion"), stock = "Susceptible"), mutate(get_summary("infection_derived_immune_proportion"), stock = "Infection-derived immune"))
  p3 <- ggplot(p3_data, aes(week_start, median, colour = stock, fill = stock)) +
    geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .12, colour = NA) + geom_line() +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) + xscale + labs(title = "C. Susceptible and immune stocks", x = NULL, y = "Population proportion") + theme_v4
  p4_data <- bind_rows(mutate(get_summary("R0"), reproduction_number = "R0"), mutate(get_summary("Reff"), reproduction_number = "Reff"))
  p4 <- ggplot(p4_data, aes(week_start, median, colour = reproduction_number, fill = reproduction_number)) +
    geom_hline(yintercept = 1, linetype = 2) + geom_ribbon(aes(ymin = q025, ymax = q975), alpha = .12, colour = NA) + geom_line() +
    xscale + labs(title = "D. Reproduction numbers", x = NULL, y = "R") + theme_v4
  annual_long <- bind_rows(
    transmute(annual, year, type = "Observed reported cases", estimate = observed_cases, lower = observed_cases, upper = observed_cases),
    transmute(annual, year, type = "Posterior expected reported cases", estimate = expected_cases_median, lower = expected_cases_q025, upper = expected_cases_q975)
  )
  p5 <- ggplot(annual_long, aes(year, estimate, colour = type)) + geom_line() + geom_point() +
    geom_errorbar(data = filter(annual_long, type != "Observed reported cases"), aes(ymin = lower, ymax = upper), width = .15) +
    labs(title = "E. Annual reported-case validation", x = NULL, y = "Cases") + theme_v4
  p6 <- ggplot(serology_windows, aes(site, state_immune_median)) +
    geom_errorbar(aes(ymin = state_immune_q025, ymax = state_immune_q975), width = .15, colour = "#0072B2") +
    geom_point(colour = "#0072B2") + geom_point(aes(y = observed_proportion), colour = "#D55E00", shape = 18, size = 2.5) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    labs(title = "F. External serology consistency", subtitle = "Blue: model; orange: survey; not fitted", x = NULL, y = "Immune proportion") + theme_v4

  subtitle <- sprintf(
    "node=%s | q_external=%.4f (fixed data, not estimated) | HMC gate: %s",
    node_id, bundle$config$q_external, if (hmc_pass) "PASS" else "FAIL"
  )
  figure <- (p1 | p2) / (p3 | p4) / (p5 | p6) + plot_annotation(
    title = "Ceara v4.1 modular-q renewal model",
    subtitle = subtitle, tag_levels = "A",
    theme = theme(plot.subtitle = element_text(face = "bold", colour = if (hmc_pass) "grey30" else "#D55E00"))
  )
  figure_path <- file.path(figure_dir, paste0("renewal_v4_1_modular_q_", node_id, "_six_panel", if (hmc_pass) "" else "_UNCONVERGED", ".pdf"))
  ggsave(figure_path, figure, width = 183, height = 225, units = "mm", device = cairo_pdf)
  message("[", node_id, "] figure saved: ", figure_path)
  invisible(figure_path)
}

if (sys.nframe() == 0L) run_plot_q_node_figures()
