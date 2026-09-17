# =============================================================================
# Nationwide Brazilian chikungunya epidemic-wave census v2.
#
# Detection uses a 3-week smoothed incidence curve only. All burden and
# early-Re pre-fit quantities use the preserved raw weekly confirmed cases.
# No Re, renewal, susceptibility, climate, or recurrence-hazard model is fit.
# =============================================================================

required_packages <- c("here", "dplyr", "tidyr", "readr", "ggplot2", "scales", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(scales); library(patchwork)
})

project_root_wave_v2 <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not locate the inner CHIK_CLIM R project.")
}

uf_lookup <- tibble::tribble(
  ~uf_code, ~state,
  "11", "RO", "12", "AC", "13", "AM", "14", "RR", "15", "PA", "16", "AP", "17", "TO",
  "21", "MA", "22", "PI", "23", "CE", "24", "RN", "25", "PB", "26", "PE", "27", "AL", "28", "SE", "29", "BA",
  "31", "MG", "32", "ES", "33", "RJ", "35", "SP", "41", "PR", "42", "SC", "43", "RS",
  "50", "MS", "51", "MT", "52", "GO", "53", "DF"
)

wave_settings <- list(
  primary = list(label = "primary", smooth_weeks = 3L, trough_fraction = .30, trough_persistence = 3L, onset_fraction = .05),
  permissive = list(label = "permissive", smooth_weeks = 3L, trough_fraction = .40, trough_persistence = 2L, onset_fraction = .025),
  strict = list(label = "strict", smooth_weeks = 3L, trough_fraction = .20, trough_persistence = 4L, onset_fraction = .10),
  generation_weights = { w <- diff(pgamma(0:8, shape = 4, rate = 2)); w / sum(w) },
  backlog_share_flag = .80,
  peak_neighborhood_weeks = 3L
)

centered_mean <- function(x, width = 3L) {
  half <- width %/% 2L
  vapply(seq_along(x), function(i) mean(x[max(1L, i - half):min(length(x), i + half)]), numeric(1))
}

run_starts <- function(condition, min_length) {
  rr <- rle(ifelse(is.na(condition), FALSE, condition))
  ends <- cumsum(rr$lengths)
  starts <- ends - rr$lengths + 1L
  starts[rr$values & rr$lengths >= min_length]
}

local_peaks <- function(x, neighborhood = 3L) {
  # Plateau-aware meaningful local maxima. A candidate must be the maximum
  # over a +/- 3-week neighbourhood; this does not further smooth incidence,
  # but prevents every short post-smoothing wiggle becoming a new peak.
  rr <- rle(x)
  ends <- cumsum(rr$lengths)
  starts <- ends - rr$lengths + 1L
  out <- integer()
  for (i in seq_along(rr$values)) {
    if (rr$values[i] <= 0 || i == 1L || i == length(rr$values)) next
    left <- rr$values[i - 1L]; right <- rr$values[i + 1L]
    lo <- max(1L, starts[i] - neighborhood); hi <- min(length(x), ends[i] + neighborhood)
    locally_highest <- rr$values[i] >= max(x[lo:hi])
    if (locally_highest && rr$values[i] >= left && rr$values[i] >= right && (rr$values[i] > left || rr$values[i] > right)) {
      out <- c(out, floor((starts[i] + ends[i]) / 2))
    }
  }
  out
}

qualifying_trough <- function(smooth, left_peak, right_peak, fraction, persistence) {
  if (right_peak - left_peak <= 1L) return(NULL)
  peak_small <- min(smooth[left_peak], smooth[right_peak])
  eligible <- smooth[(left_peak + 1L):(right_peak - 1L)] <= fraction * peak_small
  starts <- run_starts(eligible, persistence)
  if (!length(starts)) return(NULL)
  indices <- (left_peak + 1L):(right_peak - 1L)
  start_index <- indices[starts[1L]]
  # Use the whole qualifying low run around its minimum as the trough interval.
  run_end <- start_index
  while (run_end < right_peak && smooth[run_end + 1L] <= fraction * peak_small) run_end <- run_end + 1L
  min_index <- start_index - 1L + which.min(smooth[start_index:run_end])
  list(start = start_index, end = run_end, minimum = min_index,
       value = smooth[min_index], ratio = smooth[min_index] / peak_small,
       duration = run_end - start_index + 1L)
}

segment_peaks <- function(smooth, settings) {
  peaks <- local_peaks(smooth, wave_settings$peak_neighborhood_weeks)
  if (!length(peaks)) return(list(clusters = list(), boundaries = list(), peak_indices = integer()))
  if (length(peaks) == 1L) return(list(clusters = list(peaks), boundaries = list(), peak_indices = peaks))
  # Hierarchical peak merging. A small local maximum must not permanently
  # bridge two large waves merely because its own tiny height makes the
  # pairwise trough criterion unrealistically close to zero. After merging a
  # non-separated small peak, compare the dominant peaks of the new clusters.
  clusters <- lapply(peaks, function(x) x)
  changed <- TRUE
  while (changed && length(clusters) > 1L) {
    changed <- FALSE
    merged <- list()
    index <- 1L
    while (index <= length(clusters)) {
      current <- clusters[[index]]
      while (index < length(clusters)) {
        next_cluster <- clusters[[index + 1L]]
        current_peak <- current[which.max(smooth[current])]
        next_peak <- next_cluster[which.max(smooth[next_cluster])]
        trough <- qualifying_trough(smooth, current_peak, next_peak,
                                    settings$trough_fraction, settings$trough_persistence)
        if (!is.null(trough)) break
        current <- c(current, next_cluster)
        index <- index + 1L
        changed <- TRUE
      }
      merged[[length(merged) + 1L]] <- current
      index <- index + 1L
    }
    clusters <- merged
  }
  boundaries <- lapply(seq_len(max(0L, length(clusters) - 1L)), function(i) {
    left_peak <- clusters[[i]][which.max(smooth[clusters[[i]]])]
    right_peak <- clusters[[i + 1L]][which.max(smooth[clusters[[i + 1L]]])]
    qualifying_trough(smooth, left_peak, right_peak,
                      settings$trough_fraction, settings$trough_persistence)
  })
  list(clusters = clusters, boundaries = boundaries, peak_indices = peaks)
}

find_onset <- function(smooth, trough_index, peak_index, onset_fraction) {
  baseline <- smooth[trough_index]
  threshold <- baseline + onset_fraction * (smooth[peak_index] - baseline)
  if (peak_index <= trough_index + 1L) return(peak_index)
  candidates <- seq.int(trough_index + 1L, peak_index - 1L)
  sustained <- candidates[candidates + 1L <= peak_index & smooth[candidates] >= threshold & smooth[candidates + 1L] >= threshold]
  if (length(sustained)) sustained[1L] else peak_index
}

last_low_run_after_peak <- function(smooth, peak_index, fraction, persistence) {
  if (peak_index >= length(smooth)) return(NULL)
  eligible <- smooth[(peak_index + 1L):length(smooth)] <= fraction * smooth[peak_index]
  starts <- run_starts(eligible, persistence)
  if (!length(starts)) return(NULL)
  start <- peak_index + starts[1L]
  end <- start
  while (end < length(smooth) && smooth[end + 1L] <= fraction * smooth[peak_index]) end <- end + 1L
  list(start = start, end = end)
}

window_metrics <- function(raw_cases, onset_index, width, generation_weights, backlog_flag) {
  last <- onset_index + width - 1L
  if (last > length(raw_cases)) return(list(complete = FALSE, cases = NA_real_, nonzero = NA_integer_, backlog = NA, lambda_nonzero = NA_integer_))
  values <- raw_cases[onset_index:last]
  lambda_nonzero <- sum(vapply(onset_index:last, function(t) {
    history <- t - seq_along(generation_weights)
    all(history >= 1L) && sum(generation_weights * raw_cases[history]) > 0
  }, logical(1)))
  list(complete = TRUE, cases = sum(values), nonzero = sum(values > 0),
       backlog = if (sum(values) > 0) max(values) / sum(values) >= backlog_flag else FALSE,
       lambda_nonzero = lambda_nonzero)
}

eligibility_reason <- function(metrics, width, case_minimum, nonzero_minimum) {
  reasons <- character()
  if (!metrics$complete) reasons <- c(reasons, paste0("fewer than ", width, " complete surveillance weeks"))
  if (metrics$complete && metrics$nonzero < nonzero_minimum) reasons <- c(reasons, "too few nonzero raw incidence weeks")
  if (metrics$complete && metrics$cases < case_minimum) reasons <- c(reasons, paste0("fewer than ", case_minimum, " raw cases"))
  if (metrics$complete && metrics$lambda_nonzero < nonzero_minimum) reasons <- c(reasons, "too few nonzero renewal infectiousness contributions")
  if (metrics$complete && isTRUE(metrics$backlog)) reasons <- c(reasons, "one week contains >=80% of early-window cases (backlog screen)")
  paste(reasons, collapse = "; ")
}

detect_waves_one_state <- function(df, settings, primary = FALSE) {
  smooth <- df$incidence_smooth
  segments <- segment_peaks(smooth, settings)
  if (!length(segments$clusters)) return(tibble())
  n_wave <- length(segments$clusters)
  rows <- vector("list", n_wave)
  for (i in seq_len(n_wave)) {
    cluster <- segments$clusters[[i]]
    dominant_peak <- cluster[which.max(smooth[cluster])]
    trough_before <- if (i == 1L) {
      # Last baseline minimum before the first dominant peak.
      candidates <- seq_len(dominant_peak)
      candidates[max(which(smooth[candidates] == min(smooth[candidates])))]
    } else segments$boundaries[[i - 1L]]$minimum
    onset_index <- find_onset(smooth, trough_before, dominant_peak, settings$onset_fraction)
    if (i < n_wave) {
      next_boundary <- segments$boundaries[[i]]
      # The trough minimum is the non-overlapping boundary: assigning the
      # entire low interval to the preceding wave can otherwise make its
      # shading extend beyond the next wave's rising-limb onset.
      end_index <- next_boundary$minimum
      censored <- FALSE
    } else {
      terminal <- last_low_run_after_peak(smooth, dominant_peak, settings$trough_fraction, settings$trough_persistence)
      censored <- is.null(terminal)
      end_index <- if (censored) length(smooth) else terminal$end
      next_boundary <- NULL
    }
    first4 <- window_metrics(df$reported_cases, onset_index, 4L, wave_settings$generation_weights, wave_settings$backlog_share_flag)
    first6 <- window_metrics(df$reported_cases, onset_index, 6L, wave_settings$generation_weights, wave_settings$backlog_share_flag)
    first8 <- window_metrics(df$reported_cases, onset_index, 8L, wave_settings$generation_weights, wave_settings$backlog_share_flag)
    reason6 <- eligibility_reason(first6, 6L, 30L, 4L)
    rows[[i]] <- tibble(
      state = df$state[1], onset_week = df$week_start[onset_index], peak_week = df$week_start[dominant_peak],
      end_week = if (censored) as.Date(NA) else df$week_start[end_index], wave_end_censored = censored,
      peak_incidence_raw = df$incidence_raw[dominant_peak], peak_incidence_smooth = smooth[dominant_peak], peak_cases = df$reported_cases[dominant_peak],
      preceding_trough_week = df$week_start[trough_before], preceding_trough_incidence = smooth[trough_before],
      trough_ratio = if (i == 1L) NA_real_ else segments$boundaries[[i - 1L]]$ratio,
      trough_duration_weeks = if (i == 1L) NA_integer_ else segments$boundaries[[i - 1L]]$duration,
      duration_weeks = end_index - onset_index + 1L, weeks_onset_to_peak = dominant_peak - onset_index,
      total_cases = sum(df$reported_cases[onset_index:end_index]), cumulative_incidence_per_100k = sum(df$incidence_raw[onset_index:end_index]),
      first4_cases = first4$cases, first6_cases = first6$cases, first8_cases = first8$cases,
      first4_nonzero_weeks = first4$nonzero, first6_nonzero_weeks = first6$nonzero, first8_nonzero_weeks = first8$nonzero,
      re4_prefit_eligible = !nzchar(eligibility_reason(first4, 4L, 20L, 3L)),
      re6_prefit_eligible = !nzchar(reason6), re8_prefit_eligible = !nzchar(eligibility_reason(first8, 8L, 40L, 6L)),
      early_peak_flag = dominant_peak - onset_index <= 3L, reporting_anomaly_flag = isTRUE(first6$backlog),
      n_candidate_peaks_in_segment = length(cluster),
      qc_duration_gt104 = (end_index - onset_index + 1L) > 104L,
      qc_peak_late_gt52 = (dominant_peak - onset_index) > 52L,
      qc_multiple_peaks = length(cluster) > 1L,
      qc_no_meaningful_fall = censored,
      reason_re6_not_eligible = ifelse(nzchar(reason6), reason6, NA_character_)
    )
  }
  bind_rows(rows) |>
    arrange(onset_week) |>
    mutate(wave_id = paste0(state, "_wave_", sprintf("%02d", row_number())), .after = state)
}

apply_major_classification <- function(waves) {
  if (!nrow(waves)) return(waves)
  waves |>
    mutate(
      candidate_wave = TRUE,
      major_epidemic_abs50 = total_cases >= 50L,
      major_epidemic_abs100 = total_cases >= 100L,
      major_epidemic_abs200 = total_cases >= 200L,
      major_epidemic_inc2 = cumulative_incidence_per_100k >= 2,
      major_epidemic_inc5 = cumulative_incidence_per_100k >= 5,
      major_epidemic_inc10 = cumulative_incidence_per_100k >= 10,
      major_epidemic_primary = major_epidemic_abs100 & major_epidemic_inc5,
      reason_not_major = case_when(
        major_epidemic_primary ~ NA_character_,
        !major_epidemic_abs100 & !major_epidemic_inc5 ~ "fewer than 100 total cases and cumulative incidence below 5/100k",
        !major_epidemic_abs100 ~ "fewer than 100 total cases",
        TRUE ~ "cumulative incidence below 5/100k"
      )
    ) |>
    group_by(state) |>
    mutate(
      wave_order = ifelse(major_epidemic_primary, cumsum(major_epidemic_primary), NA_integer_),
      is_recurrence = ifelse(major_epidemic_primary, wave_order >= 2L, NA)
    ) |>
    ungroup()
}

detect_wave_census <- function(state_week, setting) {
  per_state <- split(state_week, state_week$state)
  bind_rows(lapply(per_state, detect_waves_one_state, settings = setting)) |>
    apply_major_classification() |>
    mutate(segmentation_definition = setting$label, .before = 1)
}

overlaps <- function(a, b) {
  # Censored final waves retain end_week = NA in the published table; for
  # sensitivity matching only, their observed interval extends to data end.
  a_end <- dplyr::coalesce(a$end_week, as.Date("2025-12-28"))
  b_end <- dplyr::coalesce(b$end_week, as.Date("2025-12-28"))
  a$onset_week <= b_end & b$onset_week <= a_end
}

attach_stability <- function(primary, all_waves) {
  if (!nrow(primary)) return(primary)
  has_match <- function(row, definition, major_only = FALSE) {
    candidates <- dplyr::filter(all_waves, segmentation_definition == definition, state == row$state)
    if (major_only) candidates <- dplyr::filter(candidates, major_epidemic_primary)
    if (!nrow(candidates)) return(FALSE)
    any(vapply(seq_len(nrow(candidates)), function(i) overlaps(row, candidates[i, ]), logical(1)))
  }
  stable_segment <- vapply(seq_len(nrow(primary)), function(i) {
    has_match(primary[i, ], "permissive") && has_match(primary[i, ], "strict")
  }, logical(1))
  stable_major <- vapply(seq_len(nrow(primary)), function(i) {
    primary$major_epidemic_primary[i] && primary$major_epidemic_abs200[i] && primary$major_epidemic_inc10[i] &&
      has_match(primary[i, ], "permissive", TRUE) && has_match(primary[i, ], "strict", TRUE)
  }, logical(1))
  primary |>
    mutate(wave_definition_stable = stable_segment, major_wave_definition_stable = stable_major)
}

summarise_census <- function(waves, definition, all_states) {
  major <- dplyr::filter(waves, major_epidemic_primary)
  n_by_state <- major |> count(state, name = "n_major")
  re_by_state <- major |> dplyr::filter(re6_prefit_eligible) |> count(state, name = "n_re")
  tibble(
    segmentation_definition = definition,
    metric = c("total_ufs", "states_with_at_least_1_major_epidemic", "states_with_no_major_epidemic",
               "total_candidate_waves", "total_primary_major_epidemic_waves", "first_major_epidemics",
               "recurrent_major_epidemic_waves", "primary_re6_eligible_major_waves",
               "recurrent_re6_eligible_major_waves", "major_waves_without_primary_re6_eligibility",
               "states_with_at_least_2_major_waves", "states_with_at_least_2_re6_eligible_waves",
               "median_major_wave_duration_weeks", "maximum_major_wave_duration_weeks"),
    value = c(length(all_states), nrow(n_by_state), length(setdiff(all_states, n_by_state$state)), nrow(waves), nrow(major),
              sum(major$wave_order == 1L), sum(major$is_recurrence), sum(major$re6_prefit_eligible),
              sum(major$is_recurrence & major$re6_prefit_eligible), sum(!major$re6_prefit_eligible),
              sum(n_by_state$n_major >= 2L), sum(re_by_state$n_re >= 2L),
              if (nrow(major)) median(major$duration_weeks) else NA_real_, if (nrow(major)) max(major$duration_weeks) else NA_real_),
    states = NA_character_
  ) |>
    mutate(states = case_when(
      metric == "states_with_no_major_epidemic" ~ paste(sort(setdiff(all_states, n_by_state$state)), collapse = "; "),
      metric == "states_with_at_least_2_major_waves" ~ paste(sort(n_by_state$state[n_by_state$n_major >= 2L]), collapse = "; "),
      metric == "states_with_at_least_2_re6_eligible_waves" ~ paste(sort(re_by_state$state[re_by_state$n_re >= 2L]), collapse = "; "),
      TRUE ~ states
    ))
}

summarise_major_size_sensitivity <- function(primary_waves) {
  grid <- tidyr::expand_grid(absolute_case_minimum = c(50L, 100L, 200L),
                             cumulative_incidence_minimum = c(2, 5, 10))
  bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    accepted <- primary_waves |>
      dplyr::filter(total_cases >= grid$absolute_case_minimum[i],
             cumulative_incidence_per_100k >= grid$cumulative_incidence_minimum[i]) |>
      arrange(state, onset_week) |>
      group_by(state) |>
      mutate(size_sensitivity_order = row_number(), size_sensitivity_recurrence = size_sensitivity_order >= 2L) |>
      ungroup()
    tibble(
      absolute_case_minimum = grid$absolute_case_minimum[i],
      cumulative_incidence_minimum = grid$cumulative_incidence_minimum[i],
      accepted_major_waves = nrow(accepted),
      recurrent_major_waves = sum(accepted$size_sensitivity_recurrence),
      re6_eligible_major_waves = sum(accepted$re6_prefit_eligible),
      recurrent_re6_eligible_major_waves = sum(accepted$size_sensitivity_recurrence & accepted$re6_prefit_eligible),
      states_with_at_least_one_major_wave = n_distinct(accepted$state)
    )
  }))
}

make_qc_table <- function(primary) {
  primary |>
    mutate(qc_any_flag = qc_duration_gt104 | qc_peak_late_gt52 | qc_multiple_peaks | qc_no_meaningful_fall | reporting_anomaly_flag | early_peak_flag) |>
    dplyr::filter(qc_any_flag) |>
    transmute(state, wave_id, onset_week, peak_week, end_week, major_epidemic_primary,
              qc_duration_gt104, qc_peak_late_gt52, qc_multiple_peaks, qc_no_meaningful_fall,
              reporting_anomaly_flag, early_peak_flag, n_candidate_peaks_in_segment)
}

make_v1_comparison <- function(root, primary) {
  v1_path <- file.path(root, "03_Output", "tables", "national_episode_census", "brazil_chik_episode_census.csv")
  if (!file.exists(v1_path)) return(tibble())
  v1 <- read_csv(v1_path, show_col_types = FALSE) |>
    mutate(onset_week = as.Date(onset_week), end_week = as.Date(end_week))
  v2_major <- dplyr::filter(primary, major_epidemic_primary)
  long_v1 <- dplyr::filter(v1, duration_weeks > 104L)
  split_count <- sum(vapply(seq_len(nrow(long_v1)), function(i) {
    candidates <- dplyr::filter(v2_major, state == long_v1$state[i], onset_week <= long_v1$end_week[i], end_week >= long_v1$onset_week[i])
    nrow(candidates) >= 2L
  }, logical(1)))
  bind_rows(
    tibble(version = c("V1 activity-period", "V2 peak-trough major-wave"),
           metric = "number_of_events", value = c(nrow(v1), nrow(v2_major))),
    tibble(version = c("V1 activity-period", "V2 peak-trough major-wave"), metric = "median_duration_weeks", value = c(median(v1$duration_weeks), median(v2_major$duration_weeks))),
    tibble(version = c("V1 activity-period", "V2 peak-trough major-wave"), metric = "maximum_duration_weeks", value = c(max(v1$duration_weeks), max(v2_major$duration_weeks))),
    tibble(version = c("V1 activity-period", "V2 peak-trough major-wave"), metric = "recurrent_event_count", value = c(sum(v1$is_recurrence), sum(v2_major$is_recurrence))),
    tibble(version = c("V1 activity-period", "V2 peak-trough major-wave"), metric = "re6_eligible_count", value = c(sum(v1$re_prefit_eligible), sum(v2_major$re6_prefit_eligible))),
    tibble(version = "V1 activity-period", metric = "V1_multiyear_episodes_split_into_at_least_2_V2_major_waves", value = split_count)
  )
}

state_wave_page <- function(df, waves) {
  major <- dplyr::filter(waves, major_epidemic_primary) |>
    mutate(recurrence_class = ifelse(is_recurrence, "Recurrence", "First major wave"),
           re_class = ifelse(re6_prefit_eligible, "Re6 eligible", "Re6 ineligible"),
           re_window_end = onset_week + 35,
           plot_end = coalesce(end_week, max(df$week_start)))
  minor <- dplyr::filter(waves, !major_epidemic_primary)
  theme_wave <- theme_classic(base_size = 8) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10))
  p_cases <- ggplot(df, aes(week_start, reported_cases)) +
    geom_rect(data = major, aes(xmin = onset_week, xmax = plot_end, ymin = -Inf, ymax = Inf, fill = recurrence_class), inherit.aes = FALSE, alpha = .16) +
    geom_rect(data = major, aes(xmin = onset_week, xmax = re_window_end, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = "#333333", alpha = .22) +
    geom_line(colour = "#222222", linewidth = .28) +
    scale_y_continuous(trans = scales::pseudo_log_trans(10), labels = scales::label_number(big.mark = ",")) +
    scale_fill_manual(values = c("First major wave" = "#56B4E9", "Recurrence" = "#E69F00"), limits = c("First major wave", "Recurrence"), drop = FALSE, name = NULL) +
    labs(title = paste0(df$state[1], ": raw weekly confirmed cases"), subtitle = "Light shading: full major wave; dark shading: exact first 6-week Re window", x = NULL, y = "Cases") + theme_wave
  p_inc <- ggplot(df, aes(week_start)) +
    geom_line(aes(y = incidence_raw), colour = "grey55", linewidth = .25) +
    geom_line(aes(y = incidence_smooth), colour = "#0072B2", linewidth = .45) +
    geom_point(data = waves, aes(x = peak_week, y = peak_incidence_smooth), inherit.aes = FALSE, shape = 4, colour = "grey20", size = 1.8) +
    geom_point(data = minor, aes(x = peak_week, y = peak_incidence_smooth), inherit.aes = FALSE, shape = 1, colour = "grey35", size = 1.7) +
    geom_point(data = major, aes(x = onset_week, y = preceding_trough_incidence, colour = recurrence_class, shape = re_class), inherit.aes = FALSE, size = 2.1) +
    geom_point(data = major |> dplyr::filter(!wave_end_censored), aes(x = end_week, y = 0), inherit.aes = FALSE, shape = 6, size = 1.8) +
    geom_point(data = major |> dplyr::filter(!is.na(trough_ratio)), aes(x = preceding_trough_week, y = preceding_trough_incidence), inherit.aes = FALSE, shape = 25, fill = "white", size = 1.8) +
    scale_y_continuous(trans = scales::pseudo_log_trans(.05), labels = scales::label_number(accuracy = .1)) +
    scale_colour_manual(values = c("First major wave" = "#0072B2", "Recurrence" = "#D55E00"), limits = c("First major wave", "Recurrence"), drop = FALSE, name = NULL) +
    scale_shape_manual(values = c("Re6 eligible" = 16, "Re6 ineligible" = 1), limits = c("Re6 eligible", "Re6 ineligible"), drop = FALSE, name = NULL) +
    labs(title = "Raw and 3-week smoothed incidence", subtitle = "Cross: candidate peak; open circle: minor candidate; coloured onset: major wave; triangle: end; diamond: qualifying trough", x = NULL, y = "Incidence per 100,000") + theme_wave
  p_cases / p_inc
}

write_visual_audit <- function(state_week, primary, path) {
  grDevices::pdf(path, width = 10, height = 7.4, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (state_name in sort(unique(state_week$state))) {
    print(state_wave_page(dplyr::filter(state_week, state == state_name), dplyr::filter(primary, state == state_name)))
  }
}

build_brazil_chik_wave_census_v2 <- function() {
  root <- project_root_wave_v2()
  table_dir <- file.path(root, "03_Output", "tables", "national_wave_census_v2")
  figure_dir <- file.path(root, "03_Output", "figures", "national_wave_census_v2")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE); dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  panel_path <- file.path(root, "01_Data", "chik_brazil_muni_week_2015_2025.rds")
  population_path <- file.path(root, "01_Data", "ibge_population_projection_uf_2024revision.rds")
  panel <- readRDS(panel_path) |>
    transmute(muni6 = sprintf("%06d", as.integer(muni6)), week_start = as.Date(week_start), reported_cases = as.numeric(cases_confirmed)) |>
    dplyr::filter(week_start >= as.Date("2015-01-04"), week_start <= as.Date("2025-12-28")) |>
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
    transmute(state, week_start, reported_cases, population, incidence_raw = reported_cases / population * 1e5) |>
    group_by(state) |>
    mutate(incidence_smooth = centered_mean(incidence_raw, 3L)) |>
    ungroup() |>
    arrange(state, week_start)
  if (anyNA(state_week$population) || n_distinct(state_week$state) != 27L) stop("Incomplete state-week input.")
  write_csv(state_week, file.path(table_dir, "brazil_chik_wave_census_v2_state_week_diagnostic.csv"))
  write_csv(tibble(
    setting = c("permissive", "primary", "strict"), smoothing_weeks = 3L, peak_neighborhood_weeks = wave_settings$peak_neighborhood_weeks,
    trough_fraction = c(.40, .30, .20), trough_persistence_weeks = c(2L, 3L, 4L), onset_amplitude_fraction = c(.025, .05, .10),
    absolute_case_minimum = c(50L, 100L, 200L), cumulative_incidence_minimum = c(2, 5, 10)
  ), file.path(table_dir, "brazil_chik_wave_census_v2_definitions.csv"))

  all_waves <- bind_rows(lapply(wave_settings[c("permissive", "primary", "strict")], detect_wave_census, state_week = state_week))
  primary <- dplyr::filter(all_waves, segmentation_definition == "primary") |> attach_stability(all_waves)
  required <- c("state", "wave_id", "onset_week", "peak_week", "end_week", "wave_end_censored", "wave_order", "is_recurrence",
                "peak_incidence_raw", "peak_incidence_smooth", "peak_cases", "preceding_trough_week", "preceding_trough_incidence", "trough_ratio", "trough_duration_weeks",
                "duration_weeks", "weeks_onset_to_peak", "total_cases", "cumulative_incidence_per_100k", "first4_cases", "first6_cases", "first8_cases",
                "first4_nonzero_weeks", "first6_nonzero_weeks", "first8_nonzero_weeks", "candidate_wave", "major_epidemic_primary",
                "major_epidemic_abs50", "major_epidemic_abs100", "major_epidemic_abs200", "major_epidemic_inc2", "major_epidemic_inc5", "major_epidemic_inc10",
                "re4_prefit_eligible", "re6_prefit_eligible", "re8_prefit_eligible", "early_peak_flag", "reporting_anomaly_flag", "wave_definition_stable",
                "major_wave_definition_stable", "reason_not_major", "reason_re6_not_eligible")
  output <- primary |> select(all_of(required))
  write_csv(output, file.path(table_dir, "brazil_chik_wave_census_v2.csv"))
  write_csv(all_waves, file.path(table_dir, "brazil_chik_wave_census_v2_sensitivity_waves.csv"))
  summary <- bind_rows(lapply(c("permissive", "primary", "strict"), function(label) summarise_census(dplyr::filter(all_waves, segmentation_definition == label), label, uf_lookup$state)))
  write_csv(summary, file.path(table_dir, "brazil_chik_wave_census_v2_summary.csv"))
  size_sensitivity <- summarise_major_size_sensitivity(primary)
  write_csv(size_sensitivity, file.path(table_dir, "brazil_chik_wave_census_v2_major_size_sensitivity.csv"))
  qc <- make_qc_table(primary); write_csv(qc, file.path(table_dir, "brazil_chik_wave_census_v2_qc_flags.csv"))
  comparison <- make_v1_comparison(root, primary); write_csv(comparison, file.path(table_dir, "brazil_chik_wave_census_v1_v2_comparison.csv"))
  write_csv(tibble(
    field = c("case_definition", "smoothing", "raw_incidence_use", "reporting_anomaly_screen", "downstream_fitting"),
    value = c("cases_confirmed for all UFs; same definition as Ceara early-Re input", "Centered 3-week mean used exclusively for detection and boundaries", "All cases, burden, and early-Re pre-fit metrics use raw weekly cases", "Automated flag when one raw week accounts for >=80% of first-six-week cases; no external anomaly register was supplied", "None: no Re or renewal fit is run by this script")
  ), file.path(table_dir, "brazil_chik_wave_census_v2_metadata.csv"))
  write_visual_audit(state_week, primary, file.path(figure_dir, "brazil_chik_wave_census_v2_visual_audit.pdf"))
  message("[save] ", file.path(table_dir, "brazil_chik_wave_census_v2.csv"))
  message("[save] ", file.path(figure_dir, "brazil_chik_wave_census_v2_visual_audit.pdf"))
  invisible(list(waves = output, summary = summary, qc = qc, comparison = comparison))
}

if (sys.nframe() == 0L) build_brazil_chik_wave_census_v2()
