# =============================================================================
# Nationwide Brazilian chikungunya epidemic-episode census (detection only).
#
# This is deliberately upstream of renewal, Re, susceptibility, FOI, climate,
# and recurrence modelling.  Episode definitions are fixed solely from weekly
# reported incidence and are assessed under one prespecified sensitivity set.
# =============================================================================

required_packages <- c("here", "dplyr", "tidyr", "readr", "ggplot2", "scales", "patchwork")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

project_root_census <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not locate the inner CHIK_CLIM R project.")
}

# All thresholds below are prespecified working definitions, not fitted values.
definitions <- tibble::tribble(
  ~definition,   ~onset_threshold, ~quiet_weeks, ~minimum_cumulative_incidence,
  "permissive",              0.1,            6L,                           2,
  "primary",                 0.2,            8L,                           5,
  "strict",                  0.5,           12L,                          10
)
quiet_threshold <- 0.1
generation_history_weeks <- 8L # existing renewal generation-interval support
backlog_share_flag <- 0.80     # transparent pre-fit screen; not used for episodes

uf_lookup <- tibble::tribble(
  ~uf_code, ~state,
  "11", "RO", "12", "AC", "13", "AM", "14", "RR", "15", "PA", "16", "AP", "17", "TO",
  "21", "MA", "22", "PI", "23", "CE", "24", "RN", "25", "PB", "26", "PE", "27", "AL", "28", "SE", "29", "BA",
  "31", "MG", "32", "ES", "33", "RJ", "35", "SP", "41", "PR", "42", "SC", "43", "RS",
  "50", "MS", "51", "MT", "52", "GO", "53", "DF"
)

consecutive_quiet_before <- function(incidence, onset_index) {
  if (onset_index <= 1L) return(0L)
  n_quiet <- 0L
  index <- onset_index - 1L
  while (index >= 1L && is.finite(incidence[index]) && incidence[index] < quiet_threshold) {
    n_quiet <- n_quiet + 1L
    index <- index - 1L
  }
  n_quiet
}

completed_quiet_run_length <- function(incidence, quiet_start, minimum_length) {
  if (is.na(quiet_start)) return(NA_integer_)
  index <- quiet_start
  n_quiet <- 0L
  while (index <= length(incidence) && is.finite(incidence[index]) && incidence[index] < quiet_threshold) {
    n_quiet <- n_quiet + 1L
    index <- index + 1L
  }
  max(n_quiet, minimum_length)
}

find_first_quiet_run <- function(incidence, first_index, quiet_length) {
  last_start <- length(incidence) - quiet_length + 1L
  if (first_index > last_start) return(NA_integer_)
  for (index in seq.int(first_index, last_start)) {
    values <- incidence[index:(index + quiet_length - 1L)]
    if (all(is.finite(values) & values < quiet_threshold)) return(index)
  }
  NA_integer_
}

first_window_summary <- function(cases, onset_index, window_length, n_total) {
  end_index <- onset_index + window_length - 1L
  if (end_index > n_total) {
    return(list(complete = FALSE, total = NA_real_, nonzero = NA_integer_, max_share = NA_real_))
  }
  values <- cases[onset_index:end_index]
  total <- sum(values)
  list(
    complete = TRUE,
    total = total,
    nonzero = sum(values > 0),
    max_share = if (total > 0) max(values) / total else NA_real_
  )
}

eligibility_reason <- function(six_week, has_history, threshold_cases) {
  reasons <- character()
  if (!six_week$complete) reasons <- c(reasons, "fewer than 6 complete surveillance weeks after onset")
  if (six_week$complete && six_week$nonzero < 4L) reasons <- c(reasons, "fewer than 4 nonzero weeks in first 6")
  if (six_week$complete && six_week$total < threshold_cases) {
    reasons <- c(reasons, paste0("fewer than ", threshold_cases, " reported cases in first 6 weeks"))
  }
  if (!has_history) reasons <- c(reasons, "fewer than 8 pre-onset weeks for generation-interval history")
  if (six_week$complete && is.finite(six_week$max_share) && six_week$max_share >= backlog_share_flag) {
    reasons <- c(reasons, "first-six cases concentrated in one week (>=80%; backlog review flag)")
  }
  paste(reasons, collapse = "; ")
}

detect_episodes_one_state <- function(df, definition) {
  df <- arrange(df, week_start)
  n_total <- nrow(df)
  onset_active <- is.finite(df$weekly_incidence_per_100k) &
    df$weekly_incidence_per_100k >= definition$onset_threshold
  rows <- list()
  cursor <- 1L
  episode_number <- 0L
  # Index immediately after the quiet run that most recently placed the
  # detector in an inter-episode state. NA means no qualifying quiet run has
  # yet been observed since the beginning of the data / a left-censored event.
  ready_onset_index <- NA_integer_
  qualifying_quiet_length <- NA_integer_

  while (cursor < n_total) {
    # A completed quiet period changes the detector into an inter-episode
    # state. It need not be immediately adjacent to onset: isolated low-level
    # reports after that quiet period do not themselves constitute an episode.
    # This is the intended hysteresis, rather than a requirement that the last
    # eight observations before a two-week onset all equal "quiet".
    left_censored_start <- cursor == 1L && onset_active[1L] && onset_active[2L]
    if (left_censored_start) {
      onset_index <- 1L
    } else {
      if (is.na(ready_onset_index)) {
        pre_onset_quiet_start <- find_first_quiet_run(
          df$weekly_incidence_per_100k, cursor, definition$quiet_weeks
        )
        if (is.na(pre_onset_quiet_start)) break
        ready_onset_index <- pre_onset_quiet_start + definition$quiet_weeks
        qualifying_quiet_length <- completed_quiet_run_length(
          df$weekly_incidence_per_100k, pre_onset_quiet_start, definition$quiet_weeks
        )
      }
      if (ready_onset_index >= n_total) break
      candidate_index <- seq.int(ready_onset_index, n_total - 1L)
      candidates <- candidate_index[onset_active[candidate_index] & onset_active[candidate_index + 1L]]
      if (!length(candidates)) break
      onset_index <- candidates[1L]
    }

    quiet_start <- find_first_quiet_run(
      df$weekly_incidence_per_100k, onset_index + 2L, definition$quiet_weeks
    )
    right_censored_end <- is.na(quiet_start)
    end_index <- if (right_censored_end) n_total else quiet_start - 1L
    next_cursor <- if (right_censored_end) n_total + 1L else quiet_start + definition$quiet_weeks
    episode_incidence <- sum(df$weekly_incidence_per_100k[onset_index:end_index])

    if (episode_incidence >= definition$minimum_cumulative_incidence) {
      episode_number <- episode_number + 1L
      window4 <- first_window_summary(df$reported_cases, onset_index, 4L, n_total)
      window6 <- first_window_summary(df$reported_cases, onset_index, 6L, n_total)
      window8 <- first_window_summary(df$reported_cases, onset_index, 8L, n_total)
      peak_index <- onset_index - 1L + which.max(df$reported_cases[onset_index:end_index])
      has_history <- onset_index > generation_history_weeks
      reason30 <- eligibility_reason(window6, has_history, threshold_cases = 30L)
      rows[[length(rows) + 1L]] <- tibble(
        state = df$state[1],
        onset_week = df$week_start[onset_index],
        peak_week = df$week_start[peak_index],
        end_week = df$week_start[end_index],
        duration_weeks = end_index - onset_index + 1L,
        weeks_onset_to_peak = peak_index - onset_index,
        episode_case_count = sum(df$reported_cases[onset_index:end_index]),
        episode_cumulative_incidence_per_100k = episode_incidence,
        peak_weekly_cases = df$reported_cases[peak_index],
        peak_weekly_incidence_per_100k = df$weekly_incidence_per_100k[peak_index],
        # Length of the completed quiet run that enabled this onset. It can
        # end before onset if isolated sub-threshold reports occur later.
        quiet_weeks_before_episode = if (left_censored_start) NA_integer_ else qualifying_quiet_length,
        first4_case_count = window4$total,
        first6_case_count = window6$total,
        first8_case_count = window8$total,
        first6_nonzero_weeks = window6$nonzero,
        re_prefit_eligible = !nzchar(reason30),
        re_prefit_eligible_n20 = !nzchar(eligibility_reason(window6, has_history, threshold_cases = 20L)),
        re_prefit_eligible_n50 = !nzchar(eligibility_reason(window6, has_history, threshold_cases = 50L)),
        reason_if_not_eligible = ifelse(nzchar(reason30), reason30, NA_character_),
        left_censored_start = left_censored_start,
        right_censored_end = right_censored_end,
        first6_max_week_share = window6$max_share
      )
    }
    # The quiet run that ended this episode is immediately the qualifying
    # inter-episode state for the next onset. Do not require a second run.
    if (!right_censored_end) {
      ready_onset_index <- quiet_start + definition$quiet_weeks
      qualifying_quiet_length <- completed_quiet_run_length(
        df$weekly_incidence_per_100k, quiet_start, definition$quiet_weeks
      )
    }
    cursor <- next_cursor
  }
  if (!length(rows)) return(tibble())
  bind_rows(rows) |>
    arrange(onset_week) |>
    mutate(
      episode_id = paste0(state, "_episode_", sprintf("%02d", row_number())),
      recurrence_order = row_number(),
      is_recurrence = recurrence_order >= 2L
    ) |>
    select(state, episode_id, everything())
}

detect_national_episodes <- function(state_week, definition) {
  by_state <- split(state_week, state_week$state)
  accepted <- bind_rows(lapply(by_state, detect_episodes_one_state, definition = definition))
  if (!nrow(accepted)) return(accepted)
  accepted |>
    mutate(definition = definition$definition, .before = 1) |>
    arrange(state, onset_week)
}

match_primary_stability <- function(primary, sensitivity) {
  if (!nrow(primary)) return(primary)
  comparison_names <- c("permissive", "strict")
  match_one <- function(row, def_name) {
    candidate <- filter(sensitivity, definition == def_name, state == row$state)
    if (!nrow(candidate)) return(FALSE)
    any(candidate$onset_week <= row$end_week & candidate$end_week >= row$onset_week)
  }
  stable <- vapply(seq_len(nrow(primary)), function(i) {
    all(vapply(comparison_names, function(def_name) match_one(primary[i, ], def_name), logical(1)))
  }, logical(1))
  notes <- vapply(seq_len(nrow(primary)), function(i) {
    found <- vapply(comparison_names, function(def_name) match_one(primary[i, ], def_name), logical(1))
    paste(paste0(comparison_names, ": ", ifelse(found, "overlapping accepted episode", "not accepted/matched")), collapse = "; ")
  }, character(1))
  primary |>
    mutate(
      episode_definition_stable = stable,
      sensitivity_notes = notes
    )
}

summarise_definition <- function(census, definition_name, all_states) {
  n_by_state <- census |>
    count(state, name = "n_episodes")
  affected <- n_by_state$state
  recurrence <- sum(census$is_recurrence)
  re_eligible <- sum(census$re_prefit_eligible)
  states_two <- n_by_state |> filter(n_episodes >= 2L) |> pull(state)
  re_by_state <- census |> filter(re_prefit_eligible) |> count(state, name = "n_re")
  states_two_re <- re_by_state |> filter(n_re >= 2L) |> pull(state)
  tibble(
    definition = definition_name,
    metric = c(
      "total_ufs", "ufs_with_at_least_one_accepted_episode", "ufs_with_no_accepted_episode",
      "total_accepted_episodes", "first_episodes", "recurrent_episodes",
      "episodes_primary_Re_prefit_eligible", "recurrent_episodes_primary_Re_prefit_eligible",
      "states_contributing_at_least_2_accepted_episodes", "states_contributing_at_least_2_Re_eligible_episodes",
      "median_episodes_per_affected_state", "maximum_episodes_in_one_state"
    ),
    value = c(
      length(all_states), length(affected), length(setdiff(all_states, affected)), nrow(census),
      sum(census$recurrence_order == 1L), recurrence, re_eligible,
      sum(census$is_recurrence & census$re_prefit_eligible), length(states_two), length(states_two_re),
      if (length(affected)) median(n_by_state$n_episodes) else NA_real_,
      if (length(affected)) max(n_by_state$n_episodes) else NA_real_
    ),
    states = c(rep(NA_character_, 10L), NA_character_, NA_character_)
  ) |>
    mutate(states = case_when(
      metric == "ufs_with_no_accepted_episode" ~ paste(sort(setdiff(all_states, affected)), collapse = "; "),
      metric == "states_contributing_at_least_2_accepted_episodes" ~ paste(sort(states_two), collapse = "; "),
      metric == "states_contributing_at_least_2_Re_eligible_episodes" ~ paste(sort(states_two_re), collapse = "; "),
      TRUE ~ states
    ))
}

ceara_quiet_audit <- function(state_week, primary_census, old_audit_path) {
  ce <- filter(state_week, state == "CE", week_start >= as.Date("2016-01-01"), week_start <= as.Date("2017-12-31"))
  quiet_run_start <- find_first_quiet_run(ce$weekly_incidence_per_100k, 1L, 8L)
  # Longest contiguous sub-threshold run is reported even if it is not 8 weeks.
  below <- ce$weekly_incidence_per_100k < quiet_threshold
  rr <- rle(ifelse(is.na(below), FALSE, below))
  longest_quiet <- if (any(rr$values)) max(rr$lengths[rr$values]) else 0L
  classification <- if (!is.na(quiet_run_start)) {
    "split: qualifying 8-week quiet interval exists between the 2016--2017 calendar years"
  } else {
    "retain: no qualifying 8-week quiet interval between the 2016--2017 calendar years"
  }
  old_ce <- if (file.exists(old_audit_path)) {
    read_csv(old_audit_path, show_col_types = FALSE) |>
      filter(uf == "CE") |>
      transmute(previous_wave_id = wave_id, previous_start_week = as.Date(start_week),
                previous_peak_week = as.Date(peak_week), previous_end_week = as.Date(end_week))
  } else tibble()
  result <- tibble(
    audit_period_start = min(ce$week_start), audit_period_end = max(ce$week_start),
    quiet_threshold_per_100k = quiet_threshold, quiet_duration_weeks = 8L,
    longest_consecutive_quiet_weeks = longest_quiet,
    qualifying_quiet_interval_start = ifelse(is.na(quiet_run_start), NA_character_, as.character(ce$week_start[quiet_run_start])),
    qualifying_quiet_interval_end = ifelse(is.na(quiet_run_start), NA_character_, as.character(ce$week_start[quiet_run_start + 7L])),
    classification = classification,
    primary_episodes_overlapping_2016_2017 = sum(primary_census$state == "CE" & primary_census$onset_week <= max(ce$week_start) & primary_census$end_week >= min(ce$week_start))
  )
  list(result = result, old_ce = old_ce)
}

episode_page <- function(state_data, state_episodes, state_label) {
  bands <- state_episodes |>
    mutate(recurrence_class = ifelse(is_recurrence, "Recurrence", "First epidemic"))
  markers <- bands |>
    mutate(re_eligibility = ifelse(re_prefit_eligible, "Re-eligible", "Re-ineligible"))
  common_theme <- theme_classic(base_size = 8) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10))
  p_cases <- ggplot(state_data, aes(week_start, reported_cases)) +
    geom_rect(data = bands, aes(xmin = onset_week, xmax = end_week, ymin = -Inf, ymax = Inf, fill = recurrence_class), inherit.aes = FALSE, alpha = .18) +
    geom_line(colour = "#333333", linewidth = .28) +
    scale_y_continuous(trans = scales::pseudo_log_trans(10), labels = scales::label_number(big.mark = ",")) +
    scale_fill_manual(values = c("First epidemic" = "#56B4E9", "Recurrence" = "#E69F00"),
                      limits = c("First epidemic", "Recurrence"), drop = FALSE, name = NULL) +
    labs(title = paste0(state_label, ": raw weekly reported cases"), x = NULL, y = "Cases") + common_theme
  p_incidence <- ggplot(state_data, aes(week_start, weekly_incidence_per_100k)) +
    geom_rect(data = bands, aes(xmin = onset_week, xmax = end_week, ymin = -Inf, ymax = Inf, fill = recurrence_class), inherit.aes = FALSE, alpha = .18) +
    geom_line(colour = "#333333", linewidth = .28) +
    geom_hline(yintercept = c(quiet_threshold, 0.2), linetype = c(2, 3), colour = c("#666666", "#D55E00"), linewidth = .35) +
    geom_point(data = markers, aes(x = onset_week, y = 0, colour = recurrence_class, shape = re_eligibility), inherit.aes = FALSE, size = 2.1) +
    geom_point(data = markers, aes(x = peak_week, y = peak_weekly_incidence_per_100k), inherit.aes = FALSE, shape = 4, size = 1.9) +
    geom_point(data = markers, aes(x = end_week, y = 0), inherit.aes = FALSE, shape = 6, size = 1.9) +
    scale_y_continuous(trans = scales::pseudo_log_trans(0.05), labels = scales::label_number(accuracy = .1)) +
    scale_fill_manual(values = c("First epidemic" = "#56B4E9", "Recurrence" = "#E69F00"),
                      limits = c("First epidemic", "Recurrence"), drop = FALSE, name = NULL) +
    scale_colour_manual(values = c("First epidemic" = "#0072B2", "Recurrence" = "#D55E00"),
                        limits = c("First epidemic", "Recurrence"), drop = FALSE, name = NULL) +
    scale_shape_manual(values = c("Re-eligible" = 16, "Re-ineligible" = 1),
                       limits = c("Re-eligible", "Re-ineligible"), drop = FALSE, name = NULL) +
    labs(title = "Weekly incidence and primary-definition episode audit", subtitle = "Onset = circle; peak = cross; end = triangle; colour = first/recurrence; fill = episode; shape = Re pre-fit eligibility", x = NULL, y = "Cases per 100,000") + common_theme
  p_cases / p_incidence
}

write_visual_audit <- function(state_week, primary_census, ceara_result, output_path) {
  grDevices::pdf(output_path, width = 10, height = 7.4, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (state_name in sort(unique(state_week$state))) {
    state_data <- filter(state_week, state == state_name)
    state_episodes <- filter(primary_census, state == state_name)
    print(episode_page(state_data, state_episodes, state_name))
  }
  ce_data <- filter(state_week, state == "CE")
  ce_episodes <- filter(primary_census, state == "CE")
  print(
    episode_page(ce_data, ce_episodes, "CEARA: dedicated 2016--2017 audit") +
      plot_annotation(
        title = ceara_result$classification,
        subtitle = paste0("Longest sub-threshold run in 2016--2017: ", ceara_result$longest_consecutive_quiet_weeks,
                          " weeks; primary detector finds ", ceara_result$primary_episodes_overlapping_2016_2017,
                          " accepted episode(s) overlapping this interval.")
      )
  )
}

build_brazil_chik_episode_census <- function() {
  root <- project_root_census()
  table_dir <- file.path(root, "03_Output", "tables", "national_episode_census")
  figure_dir <- file.path(root, "03_Output", "figures", "national_episode_census")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  panel_path <- file.path(root, "01_Data", "chik_brazil_muni_week_2015_2025.rds")
  population_path <- file.path(root, "01_Data", "ibge_population_projection_uf_2024revision.rds")
  if (!file.exists(panel_path) || !file.exists(population_path)) stop("Required national panel or UF population input is missing.")

  # `cases_confirmed` is used consistently for every UF. It exactly matches
  # the existing 2015--2025 Ceara early-Re weekly input, whereas notified
  # counts retain a high non-epidemic background that defeats the specified
  # <0.1/100k quiet-period rule. The 2014-12-28 partial week is excluded.
  panel <- readRDS(panel_path) |>
    transmute(muni6 = sprintf("%06d", as.integer(muni6)), week_start = as.Date(week_start), reported_cases = as.numeric(cases_confirmed)) |>
    filter(week_start >= as.Date("2015-01-04"), week_start <= as.Date("2025-12-28")) |>
    mutate(uf_code = substr(muni6, 1L, 2L)) |>
    inner_join(uf_lookup, by = "uf_code")
  population <- readRDS(population_path) |>
    transmute(uf_code = sprintf("%02d", as.integer(uf_code)), year = as.integer(year), population = pop_total)
  dates <- seq(as.Date("2015-01-04"), as.Date("2025-12-28"), by = "week")
  state_week <- tidyr::expand_grid(state = uf_lookup$state, week_start = dates) |>
    left_join(uf_lookup, by = "state") |>
    left_join(panel |> group_by(state, week_start) |> summarise(reported_cases = sum(reported_cases), .groups = "drop"), by = c("state", "week_start")) |>
    mutate(reported_cases = coalesce(reported_cases, 0), year = as.integer(format(week_start, "%Y"))) |>
    left_join(population, by = c("uf_code", "year")) |>
    transmute(state, week_start, reported_cases, population, weekly_incidence_per_100k = reported_cases / population * 1e5) |>
    arrange(state, week_start)
  if (anyNA(state_week$population) || any(!is.finite(state_week$weekly_incidence_per_100k))) stop("Population join failed for state-week census input.")
  if (n_distinct(state_week$state) != 27L || any(count(state_week, state)$n != length(dates))) stop("State-week input is incomplete.")
  write_csv(state_week, file.path(table_dir, "brazil_chik_episode_census_state_week_input.csv"))
  write_csv(definitions, file.path(table_dir, "brazil_chik_episode_census_definitions.csv"))
  write_csv(
    tibble(
      field = c("analysis_period", "case_definition", "weekly_case_input", "population_input",
                "quiet_period_interpretation", "quiet_weeks_before_episode", "downstream_model_fitting"),
      value = c(
        "2015-01-04 through 2025-12-28; the partial 2014-12-28 week is excluded",
        "cases_confirmed, consistently for all 27 UFs; matches the existing Ceara early-Re weekly input",
        "01_Data/chik_brazil_muni_week_2015_2025.rds",
        "01_Data/ibge_population_projection_uf_2024revision.rds",
        "A completed quiet period enables the next onset search; isolated low-level reports after it do not constitute a new episode.",
        "Length of the most recent completed qualifying quiet run that enabled onset; it is not necessarily contiguous with onset.",
        "None: this task performs episode detection and Re pre-fit eligibility assessment only."
      )
    ),
    file.path(table_dir, "brazil_chik_episode_census_metadata.csv")
  )

  all_census <- bind_rows(lapply(seq_len(nrow(definitions)), function(i) detect_national_episodes(state_week, definitions[i, ])))
  primary <- filter(all_census, definition == "primary") |>
    select(-definition) |>
    match_primary_stability(all_census)
  required_columns <- c(
    "state", "episode_id", "onset_week", "peak_week", "end_week", "recurrence_order", "is_recurrence",
    "duration_weeks", "weeks_onset_to_peak", "episode_case_count", "episode_cumulative_incidence_per_100k",
    "peak_weekly_cases", "peak_weekly_incidence_per_100k", "quiet_weeks_before_episode", "first4_case_count",
    "first6_case_count", "first8_case_count", "first6_nonzero_weeks", "re_prefit_eligible", "re_prefit_eligible_n20",
    "re_prefit_eligible_n50", "reason_if_not_eligible", "episode_definition_stable", "sensitivity_notes"
  )
  census_output <- primary |> select(all_of(required_columns))
  write_csv(census_output, file.path(table_dir, "brazil_chik_episode_census.csv"))
  write_csv(all_census, file.path(table_dir, "brazil_chik_episode_census_sensitivity_episodes.csv"))
  write_csv(primary |> select(state, episode_id, left_censored_start, right_censored_end, first6_max_week_share),
            file.path(table_dir, "brazil_chik_episode_census_prefit_audit_details.csv"))

  summary_table <- bind_rows(lapply(definitions$definition, function(name) {
    summarise_definition(filter(all_census, definition == name), name, uf_lookup$state)
  }))
  write_csv(summary_table, file.path(table_dir, "brazil_chik_episode_census_summary.csv"))
  old_audit_path <- file.path(root, "03_Output", "tables", "chik_state_epidemic_wave_audit.csv")
  ce_audit <- ceara_quiet_audit(state_week, primary, old_audit_path)
  write_csv(ce_audit$result, file.path(table_dir, "ceara_2016_2017_quiet_interval_audit.csv"))
  write_csv(ce_audit$old_ce, file.path(table_dir, "ceara_previous_wave_audit_comparison.csv"))
  write_visual_audit(state_week, primary, ce_audit$result,
                     file.path(figure_dir, "brazil_chik_episode_census_visual_audit.pdf"))
  message("[save] primary census: ", file.path(table_dir, "brazil_chik_episode_census.csv"))
  message("[save] visual audit: ", file.path(figure_dir, "brazil_chik_episode_census_visual_audit.pdf"))
  invisible(list(census = census_output, summary = summary_table, ceara = ce_audit$result))
}

if (sys.nframe() == 0L) build_brazil_chik_episode_census()
