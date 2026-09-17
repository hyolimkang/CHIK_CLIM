# =============================================================================
# Ceara major-outbreak overlay: observed incidence, annual-model S/N,
# short-window Re, and lagged climate summaries.
#
# This is a descriptive linkage of separately estimated components.  It does
# not add climate to, or refit, either the annual FOI model or renewal model.
# S_onset is necessarily the annual posterior S/N at the START of the onset
# calendar year: the annual FOI model does not identify a within-year S(t).
# =============================================================================

required_packages <- c("here", "rstan", "dplyr", "tidyr", "readr", "ggplot2", "scales", "patchwork")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

project_root_v2 <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not locate the inner CHIK_CLIM R project.")
}

qfun <- function(x, probability) unname(stats::quantile(x, probability, na.rm = TRUE))

read_major_episode_table <- function(root) {
  window_path <- file.path(
    root, "03_Output", "tables", "renewal_major_episode_re",
    "ceara_major_outbreaks_6wk", "major_episode_six_week_windows.csv"
  )
  re6_path <- file.path(
    root, "03_Output", "tables", "renewal_major_episode_re",
    "ceara_major_outbreaks_6wk", "major_episode_re6_summary.csv"
  )
  re_all_path <- file.path(
    root, "03_Output", "tables", "renewal_episode_re", "ceara_early_phase_re",
    "episode_re_summary.csv"
  )
  if (!all(file.exists(c(window_path, re6_path, re_all_path)))) {
    stop("Required major-episode renewal outputs are missing.")
  }
  episode_meta <- read_csv(window_path, show_col_types = FALSE) |>
    distinct(episode_id, wave_id, onset, peak_week, peak_cases,
             total_episode_cases, duration_weeks, cases_first_4wk, cases_first_8wk) |>
    mutate(
      onset = as.Date(onset), peak_week = as.Date(peak_week),
      short_episode_id = sub("^CE_major_", "CE_", episode_id)
    )
  re6 <- read_csv(re6_path, show_col_types = FALSE) |>
    select(episode_id, Re6_median, Re6_q975) |>
    rename(Re_early_6 = Re6_median, Re_early_6_hi = Re6_q975)
  re_sensitivity <- read_csv(re_all_path, show_col_types = FALSE) |>
    dplyr::filter(window_weeks %in% c(4L, 8L)) |>
    transmute(
      short_episode_id = episode_id,
      window_weeks,
      Re_median,
      Re_q975
    ) |>
    pivot_wider(
      names_from = window_weeks, values_from = c(Re_median, Re_q975),
      names_glue = "{.value}_early_{window_weeks}"
    ) |>
    rename(
      Re_early_4 = Re_median_early_4,
      Re_early_8 = Re_median_early_8,
      Re_early_4_hi = Re_q975_early_4,
      Re_early_8_hi = Re_q975_early_8
    )
  episode_meta |>
    left_join(re6, by = "episode_id") |>
    left_join(re_sensitivity, by = "short_episode_id") |>
    mutate(Re_early_primary = Re_early_6) |>
    arrange(onset)
}

annual_s_onset <- function(root, episodes, primary_fit) {
  fit_path <- file.path(root, "02_Script", "stan", paste0("annual_foi_shape_v2_", primary_fit, ".rds"))
  if (!file.exists(fit_path)) stop("Primary annual posterior fit is missing: ", fit_path)
  bundle <- readRDS(fit_path)
  draws <- rstan::extract(bundle$fit, pars = "S_start", permuted = TRUE)$S_start
  annual <- bundle$prepared$annual
  episode_year <- as.integer(format(episodes$onset, "%Y"))
  annual_index <- match(episode_year, annual$year)
  if (anyNA(annual_index)) stop("An episode onset is outside the annual posterior time range.")
  s_prop <- vapply(seq_along(annual_index), function(i) {
    draws[, annual_index[i]] / annual$N_start[annual_index[i]]
  }, numeric(nrow(draws)))
  # vapply produces draws x episodes (including the one-episode edge case).
  if (is.null(dim(s_prop))) s_prop <- matrix(s_prop, ncol = 1L)
  tibble(
    episode_id = episodes$episode_id,
    S_onset_median = apply(s_prop, 2, qfun, probability = .5),
    S_onset_lo = apply(s_prop, 2, qfun, probability = .025),
    S_onset_hi = apply(s_prop, 2, qfun, probability = .975),
    S_onset_basis = "Annual FOI posterior S/N at start of onset calendar year"
  )
}

build_ceara_climate <- function(root) {
  climate_path <- file.path(root, "01_Data", "chik_dlnm_panel_muni_week_2015_2025.rds")
  if (!file.exists(climate_path)) stop("Climate panel is missing: ", climate_path)
  readRDS(climate_path) |>
    mutate(
      week_start = as.Date(week_start),
      muni6 = as.character(muni6),
      uf_code = substr(muni6, 1L, 2L)
    ) |>
    dplyr::filter(uf_code == "23") |>
    group_by(week_start) |>
    summarise(
      temp = stats::weighted.mean(Tmean, population),
      rain = stats::weighted.mean(PRCP, population),
      .groups = "drop"
    ) |>
    arrange(week_start)
}

episode_climate_summaries <- function(episodes, climate) {
  bind_rows(lapply(seq_len(nrow(episodes)), function(i) {
    onset <- episodes$onset[i]
    # Inclusive lag-2 through lag-6 weekly climate: five preceding weeks.
    window <- climate |>
      dplyr::filter(week_start >= onset - 42, week_start <= onset - 14)
    if (nrow(window) != 5L) stop("Incomplete 2--6 week climate lag window for ", episodes$episode_id[i])
    tibble(
      episode_id = episodes$episode_id[i],
      temp_lag_2_6w = mean(window$temp),
      rain_lag_2_6w = mean(window$rain)
    )
  }))
}

episode_recurrence_covariates <- function(episodes) {
  episodes |>
    arrange(onset) |>
    mutate(
      recurrence_order = row_number(),
      time_since_previous_epidemic = c(
        NA_real_,
        # Audit duration includes the onset week, hence duration - 1 intervals.
        as.numeric(onset[-1] - (onset[-n()] + 7L * (duration_weeks[-n()] - 1L)),
                   units = "days") / 7
      ),
      previous_epidemic_size = c(NA_real_, head(total_episode_cases, -1L))
    )
}

write_q_implied_draws <- function(root, table_dir) {
  paths <- list.files(
    file.path(root, "02_Script", "stan"),
    pattern = "^annual_foi_shape_v2_.*\\.rds$", full.names = TRUE
  )
  paths <- paths[!grepl("chik_dynamic", paths)]
  q_draws <- bind_rows(lapply(paths, function(path) {
    bundle <- readRDS(path)
    fit_id <- sub("^annual_foi_shape_v2_(.*)\\.rds$", "\\1", basename(path))
    draws <- rstan::extract(bundle$fit, permuted = TRUE)
    observed_total <- sum(bundle$prepared$annual$observed_cases)
    latent_total <- rowSums(draws$X)
    tibble(
      fit_id = fit_id,
      observation_model = bundle$config$model,
      posterior_draw = seq_along(latent_total),
      observed_cases_total = observed_total,
      latent_infections_total = latent_total,
      q_implied_cases_per_infection = observed_total / latent_total,
      infections_per_reported_case = latent_total / observed_total,
      q_implied_gt_one = observed_total / latent_total > 1,
      # Use exact list indexing: `$q` would partially match `q_implied` in
      # SHAPE fits and can silently turn this into a matrix/list column.
      q_parameter = if ("q" %in% names(draws)) as.numeric(draws[["q"]]) else NA_real_
    )
  }))
  q_summary <- q_draws |>
    group_by(fit_id, observation_model) |>
    summarise(
      q_implied_q025 = qfun(q_implied_cases_per_infection, .025),
      q_implied_median = qfun(q_implied_cases_per_infection, .5),
      q_implied_q975 = qfun(q_implied_cases_per_infection, .975),
      infections_per_case_median = qfun(infections_per_reported_case, .5),
      Pr_q_implied_gt_one = mean(q_implied_gt_one),
      .groups = "drop"
    )
  write_csv(q_draws, file.path(table_dir, "q_implied_by_posterior_draw.csv"))
  write_csv(q_summary, file.path(table_dir, "q_implied_posterior_summary.csv"))
  q_summary
}

make_overlay_plot <- function(root, episodes, climate, episode_table, figure_dir) {
  cases <- readRDS(file.path(root, "01_Data", "ce_weekly_2014_2025.rds")) |>
    transmute(week_start = as.Date(week_start), cases = as.numeric(cases))
  annual <- read_csv(
    file.path(root, "03_Output", "tables", "annual_foi_shape_v2",
              "annual_foi_shape_v2_annual_posterior_summary.csv"),
    show_col_types = FALSE
  ) |>
    dplyr::filter(fit_id == "M1v2_SHAPE_s0_75") |>
    transmute(
      date = as.Date(paste0(year, "-01-01")),
      S_onset_lo = S_start_prop_q025,
      S_onset_median = S_start_prop_median,
      S_onset_hi = S_start_prop_q975
    )
  episode_band <- episodes |>
    transmute(onset, end = peak_week + 7, label = paste0(format(peak_week, "%Y"), " peak"))
  re_points <- episode_table |>
    left_join(episodes |> select(episode_id, peak_week), by = "episode_id") |>
    mutate(
      # Label by audit-defined onset, not peak year: this is an early-onset Re.
      label = sprintf("onset %s\\nRe = %.2f", format(onset_week, "%Y-%m"), Re_early_primary)
    )

  theme_overlay <- theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 10), plot.margin = margin(4, 8, 4, 8))
  p_cases <- ggplot(cases, aes(week_start, cases)) +
    geom_rect(data = episode_band, aes(xmin = onset, xmax = end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "#E69F00", alpha = .14) +
    geom_line(colour = "#D55E00", linewidth = .35) +
    scale_y_continuous(trans = scales::pseudo_log_trans(10),
                       labels = scales::label_number(big.mark = ",")) +
    labs(title = "Observed reported cases (gold: onset to peak of major episode)", x = NULL, y = "Weekly cases") + theme_overlay
  p_s <- ggplot(annual, aes(date, S_onset_median)) +
    geom_ribbon(aes(ymin = S_onset_lo, ymax = S_onset_hi), fill = "#56B4E9", alpha = .3) +
    geom_step(colour = "#0072B2", linewidth = .75) +
    geom_point(data = episode_table, aes(x = onset_week, y = S_onset_median), inherit.aes = FALSE,
               colour = "#0072B2", size = 1.8) +
    scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, 1)) +
    labs(title = "Inferred susceptible fraction", subtitle = "Annual-model S/N at start of each onset calendar year", x = NULL, y = "S/N") + theme_overlay
  p_re <- ggplot(re_points, aes(onset_week, Re_early_primary)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    geom_point(colour = "#CC79A7", size = 2.4) +
    geom_text(aes(label = label), vjust = -1, size = 2.7) +
    # The upper axis limit is derived from the maximum 95% interval endpoint,
    # not posterior medians. This guarantees no displayed Re interval clips.
    coord_cartesian(ylim = c(0, max(c(episodes$Re_early_4_hi,
                                     episodes$Re_early_6_hi,
                                     episodes$Re_early_8_hi), na.rm = TRUE) * 1.10)) +
    labs(title = "Early effective reproduction number", subtitle = "Reported-incidence renewal estimate; fixed first 6 weeks", x = NULL, y = expression(R[e])) + theme_overlay
  p_temp <- ggplot(climate, aes(week_start, temp)) +
    geom_line(colour = "#E69F00", linewidth = .35) +
    geom_point(data = episode_table, aes(x = onset_week, y = temp_lag_2_6w), inherit.aes = FALSE,
               colour = "#E69F00", size = 1.7) +
    labs(title = "Temperature", subtitle = "Points: population-weighted mean at lag 2--6 weeks", x = NULL, y = "Mean temperature") + theme_overlay
  p_rain <- ggplot(climate, aes(week_start, rain)) +
    geom_line(colour = "#009E73", linewidth = .35) +
    geom_point(data = episode_table, aes(x = onset_week, y = rain_lag_2_6w), inherit.aes = FALSE,
               colour = "#009E73", size = 1.7) +
    labs(title = "Rainfall", subtitle = "Points: population-weighted mean at lag 2--6 weeks", x = NULL, y = "Weekly precipitation") + theme_overlay
  combined <- p_cases / p_s / p_re / p_temp / p_rain
  ggsave(file.path(figure_dir, "ceara_major_episode_observed_overlay_inferred_drivers.pdf"),
         combined, width = 220, height = 270, units = "mm", device = cairo_pdf)
}

build_ceara_episode_driver_overlay <- function(primary_fit = "M1v2_SHAPE_s0_75") {
  root <- project_root_v2()
  table_dir <- file.path(root, "03_Output", "tables", "annual_foi_shape_v2", "ceara_major_episode_drivers")
  figure_dir <- file.path(root, "03_Output", "figures", "annual_foi_shape_v2", "ceara_major_episode_drivers")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  episodes <- read_major_episode_table(root) |> episode_recurrence_covariates()
  climate <- build_ceara_climate(root)
  episode_table <- episodes |>
    left_join(annual_s_onset(root, episodes, primary_fit), by = "episode_id") |>
    left_join(episode_climate_summaries(episodes, climate), by = "episode_id") |>
    transmute(
      episode_id,
      onset_week = onset,
      Re_early_4,
      Re_early_6,
      Re_early_8,
      Re_early_primary,
      S_onset_median,
      S_onset_lo,
      S_onset_hi,
      temp_lag_2_6w,
      rain_lag_2_6w,
      time_since_previous_epidemic,
      previous_epidemic_size,
      recurrence_order
    )
  write_csv(episode_table, file.path(table_dir, "ceara_major_episode_driver_table.csv"))
  write_csv(episodes, file.path(table_dir, "ceara_major_episode_metadata.csv"))
  write_csv(
    tibble(
      field = c("episode inclusion", "Re_early_primary", "S_onset", "climate lag", "time_since_previous_epidemic"),
      definition = c(
        "Existing audit-defined Ceara episodes with total reported cases >= 4,000.",
        "Median effective reproduction number from the fixed first 6 reported weeks.",
        "Annual FOI posterior S/N at the start of the onset calendar year; it is not a weekly S(t) estimate.",
        "Population-weighted state mean of weeks 2 through 6 before the audit-defined onset.",
        "Weeks from the prior included episode's audit end to this episode's audit onset."
      )
    ),
    file.path(table_dir, "ceara_major_episode_driver_table_metadata.csv")
  )
  write_csv(climate, file.path(table_dir, "ceara_population_weighted_weekly_climate.csv"))
  q_summary <- write_q_implied_draws(root, table_dir)
  make_overlay_plot(root, episodes, climate, episode_table, figure_dir)
  message("[save] ", file.path(table_dir, "ceara_major_episode_driver_table.csv"))
  message("[save] ", file.path(figure_dir, "ceara_major_episode_observed_overlay_inferred_drivers.pdf"))
  invisible(list(episode_table = episode_table, q_summary = q_summary))
}

if (sys.nframe() == 0L) build_ceara_episode_driver_overlay()
