# Build a weekly post-first-epidemic climate x susceptibility opportunity panel
# for the primary recurrent-analysis states (BA, RJ, MT).
#
# This is additive post-processing only. It reuses the already fitted
# first-epidemic GAM and each state's already fitted long-term reconstruction;
# it does not refit a climate or susceptibility model. The weekly panel is
# descriptive context for recurrent-wave positions, not an outbreak forecast.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "lubridate", "mgcv")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
  library(readr)
  library(tibble)
  library(lubridate)
  library(mgcv)
})

root <- here::here()
if (!file.exists(file.path(root, "CHIK_CLIM.Rproj"))) root <- file.path(root, "CHIK_CLIM")
if (!file.exists(file.path(root, "CHIK_CLIM.Rproj"))) stop("Could not identify the inner CHIK_CLIM R-project root.")
source(file.path(root, "02_Script/07_national_pipeline/06_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

PRIMARY_STATES <- c("BA", "RJ", "MT")
PRIMARY_WINDOW_WEEKS <- 6L
state_fit_paths <- c(
  BA = file.path(root, "03_Output/03_bahia_pipeline/model_fits/v4_9_replication/outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds"),
  RJ = file.path(root, "03_Output/05_rio_de_janeiro_pipeline/model_fits/v4_9_replication/outputs/caseonly/rj_global_q_case_only.rds"),
  MT = file.path(root, "03_Output/06_mato_grosso_pipeline/model_fits/v4_9_replication/outputs/caseonly/mt_global_q_case_only.rds")
)
if (!all(file.exists(state_fit_paths))) {
  stop("One or more primary-state fitted reconstructions are missing. Run the existing state analyses first.")
}

state_labels <- c(BA = "Bahia", RJ = "Rio de Janeiro", MT = "Mato Grosso")

# The GAM uses mean temperature at anchor weeks t-2:t and cumulative
# precipitation over t-6:t-2. The corresponding pre-growth susceptibility is
# six weeks before the anchor: the week immediately before a 6-week window.
weekly_climate_summary <- function(climate_state) {
  climate_state <- climate_state |> arrange(week_start)
  n <- nrow(climate_state)
  tibble(
    state = climate_state$state,
    anchor_week = climate_state$week_start,
    temperature_summary = vapply(seq_len(n), function(i) {
      idx <- i + c(-2L, -1L, 0L)
      if (any(idx < 1L)) return(NA_real_)
      mean(climate_state$Tmean[idx])
    }, numeric(1)),
    precipitation_summary = vapply(seq_len(n), function(i) {
      idx <- i + (-6L:-2L)
      if (any(idx < 1L)) return(NA_real_)
      sum(climate_state$PRCP[idx])
    }, numeric(1))
  )
}

summarise_weekly_susceptibility <- function(state, fit_path) {
  message("[weekly opportunities] extracting S_prop for ", state)
  fit_bundle <- readRDS(fit_path)
  weekly_dates <- as.Date(fit_bundle$weekly_data$week_start)
  s_prop <- rstan::extract(fit_bundle$fit, pars = "S_prop")$S_prop
  if (ncol(s_prop) != length(weekly_dates)) {
    stop("S_prop dimensions do not match weekly dates for ", state)
  }
  tibble(
    state = state,
    S_week = weekly_dates,
    S_pre_median = apply(s_prop, 2L, median),
    S_pre_lower = apply(s_prop, 2L, quantile, probs = .025),
    S_pre_upper = apply(s_prop, 2L, quantile, probs = .975)
  )
}

climate <- read_csv(file.path(paths$table, "brazil_chik_uf_weekly_climate.csv"), show_col_types = FALSE) |>
  mutate(week_start = as.Date(week_start)) |>
  dplyr::filter(state %in% PRIMARY_STATES)
training <- read_csv(file.path(paths$table, "first_epidemic_climate_training_data.csv"), show_col_types = FALSE)
gam_fit <- readRDS(file.path(paths$table, "first_epidemic_climate_gam.rds"))
temp_range <- range(training$temperature)
precip_range <- range(training$precipitation)

climate_panel <- bind_rows(lapply(PRIMARY_STATES, function(state) {
  weekly_climate_summary(climate |> dplyr::filter(state == !!state))
})) |>
  mutate(
    # This row identifier must be created before incomplete lag rows are
    # removed for prediction.  Re-numbering only the prediction subset would
    # shift every fitted value when it is joined back to the full weekly panel.
    climate_row = row_number(),
    S_week = anchor_week - PRIMARY_WINDOW_WEEKS * 7L,
    climate_in_training_support =
      !is.na(temperature_summary) & !is.na(precipitation_summary) &
      temperature_summary >= temp_range[1] & temperature_summary <= temp_range[2] &
      precipitation_summary >= precip_range[1] & precipitation_summary <= precip_range[2]
  )

climate_prediction_input <- climate_panel |>
  dplyr::filter(!is.na(temperature_summary), !is.na(precipitation_summary))
prediction <- predict(
  gam_fit,
  newdata = climate_prediction_input |> select(temperature = temperature_summary, precipitation = precipitation_summary),
  se.fit = TRUE, exclude = "s(UF)", newdata.guaranteed = TRUE
)
climate_panel <- climate_panel |>
  left_join(
    climate_prediction_input |>
      transmute(
        climate_row,
        R_climate_hat = exp(as.numeric(prediction$fit)),
        R_climate_lower = exp(as.numeric(prediction$fit) - 1.96 * as.numeric(prediction$se.fit)),
        R_climate_upper = exp(as.numeric(prediction$fit) + 1.96 * as.numeric(prediction$se.fit))
      ),
    by = "climate_row"
  ) |>
  select(-climate_row)

weekly_s <- bind_rows(lapply(PRIMARY_STATES, function(state) {
  summarise_weekly_susceptibility(state, state_fit_paths[[state]])
}))

master <- read_csv(file.path(paths$table, "brazil_chik_wave_analysis_master.csv"), show_col_types = FALSE) |>
  mutate(across(c(onset_week, end_week), as.Date))
first_wave_end <- master |>
  dplyr::filter(state %in% PRIMARY_STATES, wave_order == 1L) |>
  group_by(state) |>
  summarise(first_wave_end = max(end_week), .groups = "drop")
if (nrow(first_wave_end) != length(PRIMARY_STATES)) stop("Could not identify the first qualifying wave in every primary state.")

# Keep the existing, frozen definition of a recurrent outbreak.  This avoids
# treating short or isolated rises in reported cases as a new epidemic merely
# to create a negative-control window.  A future outcome is attached below
# only when its full 6- or 8-week follow-up is observable.
recurrent_major_waves <- master |>
  dplyr::filter(state %in% PRIMARY_STATES, wave_order > 1L) |>
  transmute(
    state,
    outcome_episode_id = wave_id,
    outcome_episode_onset = onset_week,
    outcome_episode_end = end_week,
    outcome_Re_early_6_median = Re_early_6_median,
    outcome_Re_early_6_lower = Re_early_6_q025,
    outcome_Re_early_6_upper = Re_early_6_q975,
    outcome_reported_cumulative_incidence_per_100k = cumulative_incidence_per_100k
  ) |>
  arrange(state, outcome_episode_onset)

# The confirmed-case panel used by the frozen census ends on 2025-12-28.
# Climate can be available beyond the surveillance series in future updates;
# anchoring this explicitly prevents right-censored windows from being
# classified as non-outbreak opportunities.
OBSERVATION_END <- as.Date("2025-12-28")

opportunities <- climate_panel |>
  left_join(weekly_s, by = c("state", "S_week")) |>
  left_join(first_wave_end, by = "state") |>
  dplyr::filter(
    climate_in_training_support,
    !is.na(S_pre_median),
    anchor_week >= first_wave_end + PRIMARY_WINDOW_WEEKS * 7L
  ) |>
  mutate(
    R_pred_median = R_climate_hat * S_pre_median,
    opportunity_id = paste(state, format(anchor_week, "%Y%m%d"), sep = "_"),
    state_name = factor(recode(state, !!!state_labels), levels = unname(state_labels))
  ) |>
  arrange(state, anchor_week)
if (!nrow(opportunities)) stop("No post-first weekly opportunities remain after climate-support filtering.")

# Outcome definition for the all-opportunities phase plane.  Every retained
# row is a *sliding six-week recurrent-growth opportunity* whose GAM climate
# covariates end at `anchor_week`; S_pre is the fitted susceptible fraction in
# the week immediately before that six-week window.  The primary outcome is a
# frozen major-wave onset in the subsequent 8 weeks.  The corresponding 6-week
# indicator is retained as a prespecified shorter-horizon sensitivity.  For an
# outcome event, its full reported cumulative incidence and independently
# estimated early Re are attached; non-outbreak windows receive zero future
# major-wave incidence and NA early Re rather than invented values.
opportunity_outcomes <- bind_rows(lapply(seq_len(nrow(opportunities)), function(i) {
  opportunity <- opportunities[i, ]
  followup_start <- opportunity$anchor_week + 7L
  followup_end_6w <- opportunity$anchor_week + 6L * 7L
  followup_end_8w <- opportunity$anchor_week + 8L * 7L
  followup_complete <- followup_end_8w <= OBSERVATION_END
  candidate_waves <- recurrent_major_waves |>
    dplyr::filter(
      state == opportunity$state,
      outcome_episode_onset >= followup_start,
      outcome_episode_onset <= followup_end_8w
    ) |>
    arrange(outcome_episode_onset)
  selected_wave <- candidate_waves |> slice_head(n = 1L)

  tibble(
    opportunity_id = opportunity$opportunity_id,
    outcome_followup_start = followup_start,
    outcome_followup_end_6w = followup_end_6w,
    outcome_followup_end_8w = followup_end_8w,
    outcome_followup_complete_8w = followup_complete,
    major_wave_onset_next_6w = if (followup_complete) any(candidate_waves$outcome_episode_onset <= followup_end_6w) else NA,
    major_wave_onset_next_8w = if (followup_complete) nrow(candidate_waves) > 0L else NA,
    n_major_wave_onsets_next_8w = if (followup_complete) nrow(candidate_waves) else NA_integer_,
    outcome_episode_id = if (nrow(selected_wave)) selected_wave$outcome_episode_id else NA_character_,
    outcome_episode_onset = if (nrow(selected_wave)) selected_wave$outcome_episode_onset else as.Date(NA),
    outcome_episode_end = if (nrow(selected_wave)) selected_wave$outcome_episode_end else as.Date(NA),
    outcome_Re_early_6_median = if (nrow(selected_wave)) selected_wave$outcome_Re_early_6_median else NA_real_,
    outcome_Re_early_6_lower = if (nrow(selected_wave)) selected_wave$outcome_Re_early_6_lower else NA_real_,
    outcome_Re_early_6_upper = if (nrow(selected_wave)) selected_wave$outcome_Re_early_6_upper else NA_real_,
    outcome_reported_cumulative_incidence_per_100k = if (nrow(selected_wave)) selected_wave$outcome_reported_cumulative_incidence_per_100k else 0
  )
}))

opportunities <- opportunities |>
  left_join(opportunity_outcomes, by = "opportunity_id") |>
  mutate(
    outcome_status_8w = case_when(
      !outcome_followup_complete_8w ~ "Right-censored follow-up",
      major_wave_onset_next_8w ~ "Major wave within 8 weeks",
      TRUE ~ "No major wave within 8 weeks"
    ),
    outcome_status_8w = factor(
      outcome_status_8w,
      levels = c("No major wave within 8 weeks", "Major wave within 8 weeks", "Right-censored follow-up")
    )
  )

# Existing recurrent early-growth estimates are the outcome-coded overlays.
episodes <- read_csv(file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_EPISODES.csv"), show_col_types = FALSE) |>
  dplyr::filter(UF %in% PRIMARY_STATES) |>
  mutate(episode_start = as.Date(episode_start), episode_end = as.Date(episode_end))
episode_windows <- read_csv(file.path(paths$table, "RECURRENT_EPISODES.csv"), show_col_types = FALSE) |>
  mutate(early_growth_end = as.Date(early_growth_end)) |>
  dplyr::filter(UF %in% PRIMARY_STATES)
episode_size <- master |>
  transmute(UF = state, episode_id = wave_id, reported_incidence_per_100k = cumulative_incidence_per_100k)

events <- episodes |>
  left_join(episode_windows |> select(UF, episode_id, early_growth_end), by = c("UF", "episode_id")) |>
  left_join(episode_size, by = c("UF", "episode_id")) |>
  transmute(
    state = UF,
    episode_id,
    episode_start,
    episode_end,
    anchor_week = early_growth_end,
    temperature_summary_event = temperature_summary,
    precipitation_summary_event = precipitation_summary,
    R_climate_hat_event = R_climate_hat,
    S_pre_median_event = S_pre_median,
    Re_pred_median,
    Re_obs_median,
    Re_obs_lower,
    Re_obs_upper,
    reported_incidence_per_100k,
    state_name = factor(recode(UF, !!!state_labels), levels = unname(state_labels))
  ) |>
  arrange(state, anchor_week)

# The central climate and susceptibility estimates must reproduce the existing
# recurrent episode table at the same anchor dates. Re_pred differs slightly
# because its stored value propagates Monte Carlo uncertainty.
event_alignment <- events |>
  left_join(
    opportunities |> select(state, anchor_week,
                            temperature_summary_weekly = temperature_summary,
                            precipitation_summary_weekly = precipitation_summary,
                            R_climate_hat_weekly = R_climate_hat, S_pre_median_weekly = S_pre_median),
    by = c("state", "anchor_week")
  ) |>
  mutate(
    climate_abs_diff = abs(R_climate_hat_event - R_climate_hat_weekly),
    susceptibility_abs_diff = abs(S_pre_median_event - S_pre_median_weekly)
  )
if (anyNA(event_alignment$climate_abs_diff) || anyNA(event_alignment$susceptibility_abs_diff) ||
    max(event_alignment$climate_abs_diff) > 1e-8 || max(event_alignment$susceptibility_abs_diff) > 1e-8) {
  message("=== Opportunity/event alignment diagnostic ===")
  print(as.data.frame(event_alignment |>
    select(state, episode_id, anchor_week,
           temperature_summary_event, temperature_summary_weekly,
           precipitation_summary_event, precipitation_summary_weekly,
           R_climate_hat_event, R_climate_hat_weekly, climate_abs_diff,
           S_pre_median_event, S_pre_median_weekly, susceptibility_abs_diff)), digits = 12)
  stop("Weekly opportunity panel does not reproduce the existing recurrent episode climate/S estimates.")
}

alignment_summary <- tibble(
  n_weekly_opportunities = nrow(opportunities),
  n_complete_8w_opportunities = sum(opportunities$outcome_followup_complete_8w),
  n_major_wave_outcomes_next_6w = sum(opportunities$major_wave_onset_next_6w, na.rm = TRUE),
  n_major_wave_outcomes_next_8w = sum(opportunities$major_wave_onset_next_8w, na.rm = TRUE),
  surveillance_observation_end = OBSERVATION_END,
  n_event_overlays = nrow(events),
  max_abs_climate_difference = max(event_alignment$climate_abs_diff),
  max_abs_susceptibility_difference = max(event_alignment$susceptibility_abs_diff),
  climate_support_temperature_min = temp_range[1],
  climate_support_temperature_max = temp_range[2],
  climate_support_precipitation_min = precip_range[1],
  climate_support_precipitation_max = precip_range[2]
)

write_csv(opportunities, file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITIES_WEEKLY.csv"))
write_csv(events, file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_PRIMARY_EVENTS.csv"))
write_csv(opportunity_outcomes, file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITY_OUTCOMES.csv"))
write_csv(alignment_summary, file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITY_ALIGNMENT.csv"))
message("[saved] RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITIES_WEEKLY.csv")
message("[saved] RECURRENT_CLIMATE_SUSCEPTIBILITY_PRIMARY_EVENTS.csv")
message("[saved] RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITY_OUTCOMES.csv")
message("[saved] RECURRENT_CLIMATE_SUSCEPTIBILITY_OPPORTUNITY_ALIGNMENT.csv")
message("[weekly opportunities] ", nrow(opportunities), " climate-supported post-first weekly opportunities; ",
        sum(opportunities$outcome_followup_complete_8w), " complete 8-week outcomes; ",
        sum(opportunities$major_wave_onset_next_8w, na.rm = TRUE), " major-wave outcomes; ",
        nrow(events), " event overlays")
