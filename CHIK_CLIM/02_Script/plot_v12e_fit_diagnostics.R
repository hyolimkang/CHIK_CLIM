# ---------------------------------------------------------------------------
# plot_v12e_fit_diagnostics.R
#
# Posterior diagnostics for fit_chik_ceara_stan_weekly.R.
#
# Run after fit_chik_ceara_stan_weekly.R has completed in the same R session.
# Required objects:
#   - fit_v12e
#   - ce_fit
#   - years
#   - UF_CODE, UF_NAME
#
# By default this script only prints plots in the active R plotting device.
# Set SAVE_PLOTS <- TRUE before sourcing if you also want PNG files.
# ---------------------------------------------------------------------------

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "posterior", "scales")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
  library(scales)
})

utils::globalVariables(c(
  "annual_infections", "cases_rep", "cases", "lwr50", "lwr95", "median",
  "parameter", "population", "t", "upr50", "upr95", "week_start",
  "year_index"
))

required_objects <- c("fit_v12e", "ce_fit", "years", "UF_CODE", "UF_NAME")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    "\nRun source('02_Script/fit_chik_ceara_stan_weekly.R') first."
  )
}

if (!exists("SAVE_PLOTS")) SAVE_PLOTS <- FALSE
fig_dir <- "03_Output/figures"
if (SAVE_PLOTS) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

extract_fit_draws <- function(fit) {
  if ("CmdStanMCMC" %in% class(fit)) {
    as.data.frame(fit$draws(format = "draws_df"))
  } else {
    as.data.frame(fit)
  }
}

extract_index <- function(x) {
  as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", x))
}

summarise_indexed_draws <- function(draws_df, var, index_df, value_name = "value") {
  cols <- grep(paste0("^", var, "\\["), names(draws_df), value = TRUE)
  if (length(cols) == 0) {
    stop("No posterior columns found for variable: ", var)
  }

  cols <- cols[order(extract_index(cols))]

  draws_df |>
    dplyr::select(dplyr::all_of(cols)) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(".draw"),
      names_to = "parameter",
      values_to = value_name
    ) |>
    dplyr::mutate(t = extract_index(.data$parameter)) |>
    dplyr::group_by(t) |>
    dplyr::summarise(
      median = stats::median(.data[[value_name]], na.rm = TRUE),
      lwr95 = stats::quantile(.data[[value_name]], 0.025, na.rm = TRUE),
      upr95 = stats::quantile(.data[[value_name]], 0.975, na.rm = TRUE),
      lwr50 = stats::quantile(.data[[value_name]], 0.25, na.rm = TRUE),
      upr50 = stats::quantile(.data[[value_name]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(index_df, by = "t")
}

check_time_index <- function(index_df) {
  if (any(is.na(index_df$week_start))) {
    stop("`week_start` contains missing values; cannot plot diagnostics.")
  }

  duplicate_dates <- index_df$week_start[duplicated(index_df$week_start)]
  if (length(duplicate_dates) > 0) {
    warning(
      "Duplicate `week_start` values found in ce_fit: ",
      paste(unique(duplicate_dates), collapse = ", "),
      call. = FALSE
    )
  }

  bad_order <- which(diff(index_df$week_start) <= 0)
  if (length(bad_order) > 0) {
    warning(
      "`ce_fit` is not strictly ordered by `week_start`. ",
      "The posterior draws are still joined by Stan time index `t`, ",
      "but plots will be arranged by date for display. ",
      "Refit with ce_fit arranged by `week_start` if the model itself was fit in the wrong order.",
      call. = FALSE
    )
  }

  invisible(index_df)
}

plot_posterior_ribbon <- function(df, y_label, title, subtitle = NULL,
                                  observed_df = NULL, observed_y = NULL,
                                  use_pseudo_log = FALSE) {
  df <- df |>
    dplyr::arrange(.data$week_start)

  p <- ggplot(df, aes(x = .data$week_start)) +
    geom_ribbon(aes(ymin = .data$lwr95, ymax = .data$upr95), fill = "#9ecae1", alpha = 0.25) +
    geom_ribbon(aes(ymin = .data$lwr50, ymax = .data$upr50), fill = "#3182bd", alpha = 0.25) +
    geom_line(aes(y = .data$median), colour = "#08519c", linewidth = 0.8) +
    labs(
      x = NULL,
      y = y_label,
      title = title,
      subtitle = subtitle
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  if (!is.null(observed_df) && !is.null(observed_y)) {
    observed_df <- observed_df |>
      dplyr::arrange(.data$week_start)

    p <- p +
      geom_point(
        data = observed_df,
        aes(x = .data$week_start, y = .data[[observed_y]]),
        inherit.aes = FALSE,
        size = 0.7,
        alpha = 0.65,
        colour = "black"
      )
  }

  if (use_pseudo_log) {
    p <- p +
      scale_y_continuous(
       #trans = scales::pseudo_log_trans(base = 10),
        labels = scales::label_number()
      )
  } else {
    p <- p + scale_y_continuous(labels = scales::label_number())
  }

  p
}

draws_df <- extract_fit_draws(fit_v12e)

time_index <- ce_fit |>
  dplyr::mutate(week_start = as.Date(.data$week_start)) |>
  dplyr::mutate(t = dplyr::row_number()) |>
  dplyr::select(t, week_start, year, week_of_year, cases, population)
check_time_index(time_index)

s_out <- summarise_indexed_draws(draws_df, "S_out", time_index)
s_frac <- s_out |>
  dplyr::mutate(
    median = median / population,
    lwr95 = lwr95 / population,
    upr95 = upr95 / population,
    lwr50 = lwr50 / population,
    upr50 = upr50 / population
  )

i_out <- summarise_indexed_draws(draws_df, "I_out", time_index)
new_inf <- summarise_indexed_draws(draws_df, "new_inf", time_index)
reff <- summarise_indexed_draws(draws_df, "Reff_t", time_index)
cases_rep <- summarise_indexed_draws(draws_df, "cases_rep", time_index, value_name = "cases_rep")

message("\nLargest posterior predictive upper intervals for reported cases:")
print(
  cases_rep |>
    dplyr::arrange(dplyr::desc(.data$upr95)) |>
    dplyr::select(
      week_start,
      cases,
      median,
      lwr95,
      upr95,
      lwr50,
      upr50
    ) |>
    utils::head(10)
)

p_s_count <- plot_posterior_ribbon(
  s_out,
  y_label = "Susceptible people",
  title = sprintf("Posterior susceptible population S(t), UF %s", UF_CODE),
  subtitle = "Line = posterior median; ribbons = 50% and 95% credible intervals"
)

p_s_frac <- plot_posterior_ribbon(
  s_frac,
  y_label = "Susceptible fraction",
  title = sprintf("Posterior susceptible fraction S(t) / N(t), UF %s", UF_CODE),
  subtitle = "Useful for checking depletion after large outbreaks"
) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA))

p_i <- plot_posterior_ribbon(
  i_out,
  y_label = "Current infectious people",
  title = sprintf("Posterior infectious state I(t), UF %s", UF_CODE),
  subtitle = "Latent state implied by fitted transmission dynamics",
  use_pseudo_log = FALSE
)

p_new_inf <- plot_posterior_ribbon(
  new_inf,
  y_label = "Weekly infections",
  title = sprintf("Posterior weekly infections, UF %s", UF_CODE),
  subtitle = "These are inferred infections, not reported cases",
  use_pseudo_log = FALSE
)

p_reff <- plot_posterior_ribbon(
  reff,
  y_label = "Effective reproduction number",
  title = sprintf("Posterior R_eff(t), UF %s", UF_CODE),
  subtitle = "Dashed line marks R_eff = 1"
) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey35")

p_cases_ppc <- plot_posterior_ribbon(
  cases_rep,
  y_label = "Weekly reported cases",
  title = sprintf("Posterior predictive check: reported cases, UF %s", UF_CODE),
  subtitle = "Black dots = observed cases; ribbons = posterior predictive intervals",
  observed_df = time_index,
  observed_y = "cases",
  use_pseudo_log = FALSE
)

annual_cols <- grep("^annual_infections\\[", names(draws_df), value = TRUE)
annual_plot <- NULL
if (length(annual_cols) > 0) {
  annual_cols <- annual_cols[order(extract_index(annual_cols))]
  annual_summary <- draws_df |>
    dplyr::select(dplyr::all_of(annual_cols)) |>
    dplyr::mutate(.draw = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(".draw"),
      names_to = "parameter",
      values_to = "annual_infections"
    ) |>
    dplyr::mutate(year_index = extract_index(parameter)) |>
    dplyr::group_by(year_index) |>
    dplyr::summarise(
      median = stats::median(annual_infections, na.rm = TRUE),
      lwr95 = stats::quantile(annual_infections, 0.025, na.rm = TRUE),
      upr95 = stats::quantile(annual_infections, 0.975, na.rm = TRUE),
      lwr50 = stats::quantile(annual_infections, 0.25, na.rm = TRUE),
      upr50 = stats::quantile(annual_infections, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(year = years[year_index])

  annual_plot <- ggplot(annual_summary, aes(x = year, y = median)) +
    geom_linerange(aes(ymin = lwr95, ymax = upr95), linewidth = 1.1, colour = "#9ecae1") +
    geom_linerange(aes(ymin = lwr50, ymax = upr50), linewidth = 3.0, colour = "#3182bd") +
    geom_point(size = 2, colour = "#08519c") +
    labs(
      x = NULL,
      y = "Annual infections",
      title = sprintf("Posterior annual infections, UF %s", UF_CODE),
      subtitle = "Thick intervals = 50%; thin intervals = 95%"
    ) +
    scale_y_continuous(labels = scales::label_number()) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

plots <- list(
  susceptible_count = p_s_count,
  susceptible_fraction = p_s_frac,
  infectious_state = p_i,
  weekly_infections = p_new_inf,
  reff = p_reff,
  posterior_predictive_cases = p_cases_ppc
)

if (!is.null(annual_plot)) {
  plots$annual_infections <- annual_plot
}

if (SAVE_PLOTS) {
  for (nm in names(plots)) {
    out_file <- file.path(fig_dir, sprintf("%s_v12e_%s.png", UF_NAME, nm))
    ggsave(out_file, plots[[nm]], width = 8.5, height = 4.8, dpi = 300)
    message("Saved: ", out_file)
  }
}

print(p_s_frac)
print(p_cases_ppc)
print(p_reff)

