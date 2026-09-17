# Recurrent-epidemic climate x susceptibility phase analysis -- Sections
# 1-5: episode selection, pre-wave susceptibility (from each state's
# ALREADY FITTED, existing long-term v4.9 reconstruction -- no new long-
# term model fit), climate-GAM prediction (population-level fixed effect
# only, reusing the ALREADY FITTED first-epidemic climate GAM -- no
# refit), mechanistic Re_pred = R_climate_hat * S_pre_prop with Monte
# Carlo uncertainty propagation, and independent observed early Re (from
# the existing national early-Re pipeline, NOT derived from long-term S).
#
# Primary states: BA (serology anchored), RJ (partially identified), MT
# (case-based, identification diagnostics completed -- fixed-q profile
# classified PARTIAL IDENTIFICATION). CE is CONDITIONAL/SENSITIVITY only,
# reported across its 4 predefined fixed-q scenarios (0.05/0.10/0.15/0.20)
# -- never collapsed to one point. PE is excluded (absolute S not
# robustly identified).

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "mgcv")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(mgcv) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/07_national_pipeline/06_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

PRIMARY_WINDOW_WEEKS <- 6L # same "primary" window used throughout the first-epidemic climate GAM

# ================================================================
# Section 1: recurrent episode selection
# ================================================================
master <- read_csv(file.path(paths$table, "brazil_chik_wave_analysis_master.csv"), show_col_types = FALSE) |>
  mutate(onset_week = as.Date(onset_week), end_week = as.Date(end_week))

recurrent <- master |> dplyr::filter(state %in% c("BA", "RJ", "MT", "CE"), wave_order > 1, re6_analysis_usable) |>
  transmute(UF = state, episode_id = wave_id,
            episode_start = onset_week, episode_end = end_week,
            early_growth_start = onset_week,
            early_growth_end = onset_week + (PRIMARY_WINDOW_WEEKS - 1L) * 7L,
            total_cases, Re_early_6_median, Re_early_6_q025, Re_early_6_q975) |>
  arrange(UF, episode_start)

message("=== Section 1: recurrent episodes selected ===")
print(as.data.frame(recurrent |> select(UF, episode_id, episode_start, episode_end, early_growth_start, early_growth_end)))
write_csv(recurrent |> select(UF, episode_id, episode_start, episode_end, early_growth_start, early_growth_end),
          file.path(paths$table, "RECURRENT_EPISODES.csv"))
message("[saved] ", file.path(paths$table, "RECURRENT_EPISODES.csv"))

# ================================================================
# Section 2: pre-wave susceptibility from each state's existing
# long-term reconstruction. IDENTIFICATION STATUS is carried explicitly.
# ================================================================
state_fits <- list(
  BA = list(path = file.path(root, "03_Output/model_fits/bahia/v4_9_replication/outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds"),
            identification_status = "serology_anchored", q_scenario = NA_character_),
  RJ = list(path = file.path(root, "03_Output/model_fits/rio_de_janeiro/v4_9_replication/outputs/caseonly/rj_global_q_case_only.rds"),
            identification_status = "partially_identified_case_only", q_scenario = NA_character_),
  MT = list(path = file.path(root, "03_Output/model_fits/mato_grosso/v4_9_replication/outputs/caseonly/mt_global_q_case_only.rds"),
            identification_status = "case_based_partial_identification", q_scenario = NA_character_)
)
CE_Q_SCENARIOS <- c("0.05", "0.10", "0.15", "0.20")
for (q in CE_Q_SCENARIOS) {
  state_fits[[paste0("CE_q", q)]] <- list(
    path = file.path(root, sprintf("03_Output/model_fits/ceara/v4_9_hierarchical_seasonality/outputs/q%s/renewal_ceara_v4_9_fit_q%s.rds", q, q)),
    identification_status = "conditional_fixed_q_sensitivity", q_scenario = q
  )
}

get_S_pre <- function(fit_bundle, t_pre_date) {
  dates <- as.Date(fit_bundle$weekly_data$week_start)
  idx <- which(dates == t_pre_date)
  if (length(idx) == 0) idx <- which.min(abs(dates - t_pre_date))
  S_draws <- rstan::extract(fit_bundle$fit, pars = "S_prop")$S_prop[, idx]
  list(median = median(S_draws), lower = quantile(S_draws, .025), upper = quantile(S_draws, .975), draws = S_draws, matched_date = dates[idx])
}

message("\n=== Section 2: loading state fits ===")
loaded_fits <- lapply(state_fits, function(x) { message("  loading ", x$path); readRDS(x$path) })
names(loaded_fits) <- names(state_fits)

s_pre_rows <- list()
for (i in seq_len(nrow(recurrent))) {
  ep <- recurrent[i, ]
  t_pre <- ep$early_growth_start - 7L # week immediately BEFORE the recurrent episode's early-growth window
  if (ep$UF == "CE") {
    for (q in CE_Q_SCENARIOS) {
      key <- paste0("CE_q", q)
      s <- get_S_pre(loaded_fits[[key]], t_pre)
      s_pre_rows[[length(s_pre_rows) + 1L]] <- tibble(
        UF = ep$UF, episode_id = ep$episode_id, q_scenario = q, identification_status = state_fits[[key]]$identification_status,
        t_pre = t_pre, t_pre_matched = s$matched_date,
        S_pre_median = s$median, S_pre_lower = s$lower, S_pre_upper = s$upper
      )
    }
  } else {
    s <- get_S_pre(loaded_fits[[ep$UF]], t_pre)
    s_pre_rows[[length(s_pre_rows) + 1L]] <- tibble(
      UF = ep$UF, episode_id = ep$episode_id, q_scenario = NA_character_, identification_status = state_fits[[ep$UF]]$identification_status,
      t_pre = t_pre, t_pre_matched = s$matched_date,
      S_pre_median = s$median, S_pre_lower = s$lower, S_pre_upper = s$upper
    )
  }
}
s_pre_tbl <- bind_rows(s_pre_rows)
message("=== Section 2: pre-wave susceptibility ===")
print(as.data.frame(s_pre_tbl), digits = 3)

# ================================================================
# Section 3: climate-GAM prediction (population-level fixed effect only)
# ================================================================
climate <- read_csv(file.path(paths$table, "brazil_chik_uf_weekly_climate.csv"), show_col_types = FALSE) |> mutate(week_start = as.Date(week_start))
gam_fit <- readRDS(file.path(paths$table, "first_epidemic_climate_gam.rds"))
training <- read_csv(file.path(paths$table, "first_epidemic_climate_training_data.csv"), show_col_types = FALSE)
TEMP_RANGE <- range(training$temperature); PRECIP_RANGE <- range(training$precipitation)

climate_window_summary <- function(climate_uf, anchor_week) {
  temp_weeks <- anchor_week + (-2:0) * 7L
  precip_weeks <- anchor_week + (-6:-2) * 7L
  temp_rows <- climate_uf |> dplyr::filter(week_start %in% temp_weeks)
  precip_rows <- climate_uf |> dplyr::filter(week_start %in% precip_weeks)
  list(temp = if (nrow(temp_rows) == 3) mean(temp_rows$Tmean) else NA_real_,
       precip = if (nrow(precip_rows) == 5) sum(precip_rows$PRCP) else NA_real_)
}

climate_pred_rows <- bind_rows(lapply(seq_len(nrow(recurrent)), function(i) {
  ep <- recurrent[i, ]
  climate_uf <- climate |> dplyr::filter(state == ep$UF)
  w <- climate_window_summary(climate_uf, ep$early_growth_end)
  newdat <- data.frame(temperature = w$temp, precipitation = w$precip)
  pred <- predict(gam_fit, newdata = newdat, se.fit = TRUE, exclude = "s(UF)", newdata.guaranteed = TRUE)
  log_fit <- as.numeric(pred$fit); log_se <- as.numeric(pred$se.fit)
  in_support <- !is.na(w$temp) && !is.na(w$precip) &&
    w$temp >= TEMP_RANGE[1] && w$temp <= TEMP_RANGE[2] && w$precip >= PRECIP_RANGE[1] && w$precip <= PRECIP_RANGE[2]
  tibble(UF = ep$UF, episode_id = ep$episode_id, temperature_summary = w$temp, precipitation_summary = w$precip,
         R_climate_hat = exp(log_fit), R_climate_lower = exp(log_fit - 1.96 * log_se), R_climate_upper = exp(log_fit + 1.96 * log_se),
         R_climate_log_fit = log_fit, R_climate_log_se = log_se,
         climate_in_training_support = if (in_support) "YES" else "NO")
}))
message("\n=== Section 3: climate-GAM predictions (fixed effect only, UF random effect excluded) ===")
print(as.data.frame(climate_pred_rows), digits = 3)

# ================================================================
# Section 4: Re_pred = R_climate_hat * S_pre_prop, Monte Carlo uncertainty
# ================================================================
N_MC <- 2000L
set.seed(20260916L)
re_pred_rows <- bind_rows(lapply(seq_len(nrow(s_pre_tbl)), function(i) {
  s_row <- s_pre_tbl[i, ]
  c_row <- climate_pred_rows |> dplyr::filter(UF == s_row$UF, episode_id == s_row$episode_id)
  if (s_row$UF == "CE") {
    S_draws <- get_S_pre(loaded_fits[[paste0("CE_q", s_row$q_scenario)]], s_row$t_pre)$draws
  } else {
    S_draws <- get_S_pre(loaded_fits[[s_row$UF]], s_row$t_pre)$draws
  }
  S_mc <- sample(S_draws, N_MC, replace = TRUE)
  R_climate_mc <- exp(rnorm(N_MC, c_row$R_climate_log_fit, c_row$R_climate_log_se))
  Re_pred_mc <- R_climate_mc * S_mc
  tibble(UF = s_row$UF, episode_id = s_row$episode_id, q_scenario = s_row$q_scenario,
         Re_pred_median = median(Re_pred_mc), Re_pred_lower = quantile(Re_pred_mc, .025), Re_pred_upper = quantile(Re_pred_mc, .975))
}))
message("\n=== Section 4: Re_pred = R_climate_hat x S_pre_prop (Monte Carlo, N=", N_MC, ") ===")
print(as.data.frame(re_pred_rows), digits = 3)

# ================================================================
# Combine everything + Section 5 (Re_obs already in `recurrent`, from the
# independent national early-Re pipeline -- NOT derived from long-term S)
# ================================================================
combined <- s_pre_tbl |>
  left_join(climate_pred_rows, by = c("UF", "episode_id")) |>
  left_join(re_pred_rows, by = c("UF", "episode_id", "q_scenario")) |>
  left_join(recurrent |> select(UF, episode_id, episode_start, episode_end, total_cases, Re_obs_median = Re_early_6_median, Re_obs_lower = Re_early_6_q025, Re_obs_upper = Re_early_6_q975),
            by = c("UF", "episode_id"))

detail <- read_csv(file.path(paths$table, "brazil_chik_major_wave_early_re.csv"), show_col_types = FALSE) |> select(state, wave_id, peak_cases)
combined <- combined |> left_join(detail, by = c("UF" = "state", "episode_id" = "wave_id"))

message("\n=== Combined recurrent climate x susceptibility table ===")
print(as.data.frame(combined), digits = 3)

saveRDS(list(combined = combined, loaded_fits_paths = lapply(state_fits, `[[`, "path")), file.path(paths$table, "recurrent_climate_susceptibility_intermediate.rds"))
write_csv(combined, file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_EPISODES.csv"))
message("\n[saved] ", file.path(paths$table, "RECURRENT_CLIMATE_SUSCEPTIBILITY_EPISODES.csv"))
message("[saved] ", file.path(paths$table, "recurrent_climate_susceptibility_intermediate.rds"))
