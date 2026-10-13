# PHASE 3 safeguard #4 -- a DEDICATED seed-free training-data builder.
#
# The accepted make_v4_8_data()/make_v4_9_data() chain REQUIRES the 2022
# seed date to exist in the supplied weekly series (it calls
# `which(dates == seed_start_date)` and `stop()`s if not found). Both H1
# (train ends 2019-12-29) and H2 (train ends 2021-12-26) end BEFORE
# 2022-01-02, so that date never occurs in training data. Rather than
# silently coercing the existing builder to accept this, this is a
# SEPARATE function, built on make_v4_3_data() (which has no seed logic
# at all) plus the identical second-harmonic addition make_v4_9_data()
# applies -- but with is_seed/X_seed explicitly zero and fit_index
# covering every training week (no week is ever excluded from the
# training likelihood).
#
# Source AFTER 00_config_year_effect_validation.R and the v4.0/v4.3
# legacy_model_function chain (for generation_weights(), make_v4_3_data(),
# find_sero_week_index(), build_state_weekly()).

build_training_only_tau_data <- function(train_start, train_end, tau_prior_scale, serology_cutoff = NULL) {
  panel <- readRDS(file.path(ROOT, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
  weekly_cases <- build_state_weekly(panel, "23") |> dplyr::filter(week_start >= train_start, week_start <= train_end)
  weekly_demography <- readRDS(file.path(ROOT, "01_Data/ceara_weekly_demography_2015_2025.rds")) |>
    dplyr::filter(week_start >= train_start, week_start <= train_end)
  weekly <- weekly_cases |> dplyr::select(-population) |> dplyr::inner_join(weekly_demography, by = "week_start") |> dplyr::arrange(week_start)
  if (nrow(weekly) != nrow(weekly_cases) || nrow(weekly) != nrow(weekly_demography)) stop("Panel mismatch in training window.")

  # Safeguard #5: serology (2018-09-16 midpoint) only included if it
  # precedes the training cutoff -- true for both H1 and H2, checked explicitly here.
  sero_week <- find_sero_week_index(as.Date(weekly$week_start))
  sero_cutoff <- if (is.null(serology_cutoff)) train_end else serology_cutoff
  serology_included <- sero_week$matched_week_start <= sero_cutoff
  if (!serology_included) stop("Serology observation date postdates the training cutoff -- would leak holdout information. Refusing to build training data with serology included.")

  prepared <- make_v4_3_data(weekly, FIXED_FULL$imports_per_week, tau_prior_scale,
                              sero_week$index, FIXED_FULL$sero_pos, FIXED_FULL$sero_n, FIXED_FULL$kappa_sero)
  n <- nrow(weekly)
  prepared$stan_data$seasonal_sin2 <- sin(4 * pi * seq_len(n) / 52.1775) # identical to make_v4_9_data()'s own addition
  prepared$stan_data$seasonal_cos2 <- cos(4 * pi * seq_len(n) / 52.1775)
  prepared$stan_data$tau_prior_scale <- prepared$stan_data$year_effect_prior_sd
  prepared$stan_data$year_effect_prior_sd <- NULL

  # NO seed, by construction -- every training week is used in the likelihood.
  prepared$stan_data$is_seed <- integer(n)
  prepared$stan_data$X_seed <- numeric(n)
  prepared$stan_data$N_fit <- n
  prepared$stan_data$fit_index <- seq_len(n)
  prepared$stan_data$q <- QVAL

  anom <- read_csv(ACCEPTED_ANOMALY_CSV, show_col_types = FALSE) |>
    dplyr::filter(state == "CE", lag_spec == "primary") |> dplyr::mutate(week_start = as.Date(week_start))
  joined <- tibble::tibble(week_start = as.Date(weekly$week_start)) |> dplyr::left_join(anom |> dplyr::select(week_start, z_T_anom, z_P_anom), by = "week_start")
  if (nrow(joined) != nrow(weekly)) stop("Climate anomaly join produced wrong row count (training window).")
  joined$z_T_anom[is.na(joined$z_T_anom)] <- 0
  joined$z_P_anom[is.na(joined$z_P_anom)] <- 0
  prepared$stan_data$z_T_anom <- joined$z_T_anom
  prepared$stan_data$z_P_anom <- joined$z_P_anom

  list(stan_data = prepared$stan_data, weekly = weekly, sero_week = sero_week, years = prepared$years,
       serology_included = serology_included, n_weeks = n, train_start = train_start, train_end = train_end)
}
