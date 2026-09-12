# Master per-wave analysis table for the provisional Re x climate x
# susceptibility checkpoint: joins the already-frozen wave census + early-Re
# outputs (brazil_chik_wave_analysis_master.csv) with pre-onset climate
# windows and the wave-onset susceptibility (from script 10, tagged by
# FOI-scale scenario).
#
# Analysis population: major_epidemic_primary == TRUE & re6_analysis_usable
# == TRUE (Section 1). Susceptibility source is ALWAYS the M0 (mean-
# preserving) annual-SHAPE-based weekly reconstruction, at weekly resolution
# for every wave (Section 2) -- no annual-resolution fallback is needed since
# brazil_chik_wave_onset_susceptibility_summary_scale*.csv covers all major
# waves at weekly resolution.

required_episode_packages <- c("here", "dplyr", "readr", "lubridate", "tibble")
missing_episode_packages <- required_episode_packages[
  !vapply(required_episode_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_episode_packages)) stop("Missing package(s): ", paste(missing_episode_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(lubridate); library(tibble) })

episode_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(episode_root(), "02_Script", "40_renewal_model", "15_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

episode_settings <- list(
  foi_scale_tag = Sys.getenv("FOI_SCALE_TAG", "1_0"),
  windows = list(
    pre6_2 = -6:-2,  # PRIMARY (Section 3)
    pre4_1 = -4:-1,  # sensitivity
    pre8_4 = -8:-4   # sensitivity
  )
)

# For one wave onset, compute mean Tmean and summed PRCP over
# [onset_week + window[1]*7 days, onset_week + window[length]*7 days],
# i.e. whole-week offsets strictly before onset (all offsets are negative).
climate_window_summary <- function(climate_state, onset_week, window_offsets) {
  target_weeks <- onset_week + window_offsets * 7L
  rows <- climate_state |> filter(week_start %in% target_weeks)
  if (nrow(rows) != length(window_offsets)) {
    return(list(temp = NA_real_, precip = NA_real_, n_weeks_found = nrow(rows)))
  }
  list(temp = mean(rows$Tmean, na.rm = FALSE), precip = sum(rows$PRCP, na.rm = FALSE), n_weeks_found = nrow(rows))
}

run_build_episode_climate_dataset <- function() {
  paths <- ensure_national_output_dirs()
  master_path <- file.path(paths$table, "brazil_chik_wave_analysis_master.csv")
  climate_path <- file.path(paths$table, "brazil_chik_uf_weekly_climate.csv")
  s_summary_path <- file.path(paths$table, sprintf("brazil_chik_wave_onset_susceptibility_summary_scale%s.csv", episode_settings$foi_scale_tag))
  for (p in c(master_path, climate_path, s_summary_path)) if (!file.exists(p)) stop("Required input missing: ", p)

  master <- read_csv(master_path, show_col_types = FALSE) |> mutate(onset_week = as.Date(onset_week))
  climate <- read_csv(climate_path, show_col_types = FALSE) |> mutate(week_start = as.Date(week_start))
  s_draws_summary <- read_csv(s_summary_path, show_col_types = FALSE)

  # Section 1: analysis population. brazil_chik_wave_analysis_master.csv is
  # already pre-filtered to major_epidemic_primary == TRUE internally by
  # 04_build_brazil_wave_analysis_master.R (that now-constant column is
  # dropped from its output) -- only re6_analysis_usable needs filtering here.
  if ("major_epidemic_primary" %in% names(master)) {
    population <- master |> filter(major_epidemic_primary, re6_analysis_usable)
  } else {
    population <- master |> filter(re6_analysis_usable)
  }
  n_before_climate <- nrow(population)
  message("[episode-dataset] analysis population (major_epidemic_primary & re6_analysis_usable): ", n_before_climate, " waves")

  exclusions <- list()
  rows <- vector("list", nrow(population))
  for (i in seq_len(nrow(population))) {
    wave <- population[i, ]
    climate_state <- climate |> filter(state == wave$state)
    w1 <- climate_window_summary(climate_state, wave$onset_week, episode_settings$windows$pre6_2)
    w2 <- climate_window_summary(climate_state, wave$onset_week, episode_settings$windows$pre4_1)
    w3 <- climate_window_summary(climate_state, wave$onset_week, episode_settings$windows$pre8_4)

    missing_reason <- NULL
    if (anyNA(c(w1$temp, w1$precip))) missing_reason <- "missing PRIMARY pre-onset climate window (-6:-2)"

    s_row <- s_draws_summary |> filter(state == wave$state, wave_id == wave$wave_id)
    if (!nrow(s_row)) {
      missing_reason <- c(missing_reason, "no wave-onset susceptibility draw summary (script 10)")
    }

    if (!is.null(missing_reason)) {
      exclusions[[length(exclusions) + 1L]] <- tibble(state = wave$state, wave_id = wave$wave_id, reason = paste(missing_reason, collapse = "; "))
      rows[[i]] <- NULL
      next
    }

    rows[[i]] <- tibble(
      state = wave$state, wave_id = wave$wave_id, wave_order = wave$wave_order, is_recurrence = wave$is_recurrence,
      onset_week = wave$onset_week, onset_year = wave$calendar_year,
      Re_early_6_median = wave$Re_early_6_median, Re_early_6_q025 = wave$Re_early_6_q025, Re_early_6_q975 = wave$Re_early_6_q975,
      Re_early_4_median = wave$Re_early_4_median, Re_early_8_median = wave$Re_early_8_median,
      S_prop = s_row$S_onset_median, S_prop_q025 = s_row$S_onset_q025, S_prop_q975 = s_row$S_onset_q975,
      susceptibility_source = "M0_mean_preserving_annual_shape",
      susceptibility_time_resolution = "weekly_at_onset",
      temp_pre6_2 = w1$temp, precip_pre6_2 = w1$precip,
      temp_pre4_1 = w2$temp, precip_pre4_1 = w2$precip,
      temp_pre8_4 = w3$temp, precip_pre8_4 = w3$precip,
      weeks_since_previous_wave = wave$weeks_since_previous_wave, previous_wave_total_cases = wave$previous_wave_total_cases
    )
  }
  episode <- bind_rows(rows)
  exclusion_log <- bind_rows(exclusions)

  n_excluded <- if (nrow(exclusion_log)) nrow(exclusion_log) else 0L
  message("[episode-dataset] retained ", nrow(episode), " / ", n_before_climate, " waves (", n_excluded, " excluded for missing predictor data)")
  if (n_excluded) {
    message("[episode-dataset] exclusion reasons:")
    print(as.data.frame(exclusion_log))
  }
  if (any(episode$Re_early_6_median <= 0, na.rm = TRUE)) stop("Re_early_6_median <= 0 present -- required > 0 per Section 5.")

  out_path <- file.path(paths$table, sprintf("brazil_chik_Re_S_climate_episode_analysis_scale%s.csv", episode_settings$foi_scale_tag))
  write_csv(episode, out_path)
  exclusion_path <- file.path(paths$table, sprintf("brazil_chik_Re_S_climate_episode_exclusions_scale%s.csv", episode_settings$foi_scale_tag))
  write_csv(exclusion_log, exclusion_path)
  message("[episode-dataset] saved: ", out_path)
  message("[episode-dataset] saved: ", exclusion_path)
  invisible(list(episode = episode, exclusions = exclusion_log))
}

if (sys.nframe() == 0L) run_build_episode_climate_dataset()
