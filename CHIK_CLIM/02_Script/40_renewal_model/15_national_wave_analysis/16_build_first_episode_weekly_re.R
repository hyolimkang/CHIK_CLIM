# Brazil first-epidemic climate-transmission pilot -- Section 2: weekly
# early-wave Re(t) for eligible first epidemics.
#
# IMPLEMENTATION NOTE (documented per instruction to reuse the frozen
# method exactly, not invent a new one): the existing frozen early-Re
# Stan model (renewal_ceara_episode_re.stan, run once nationally in
# 01_fit_brazil_major_wave_early_re.R) estimates ONE constant Re per
# (wave, window-length) using windows of 4, 6 and 8 weeks from the
# episode onset -- it does NOT produce a continuously time-varying
# weekly Re(t) trajectory, and building one would require a new
# estimation method, which is explicitly out of scope here. To satisfy
# "one row per usable week" while reusing the frozen fits with ZERO
# re-estimation, each of the three already-fitted window-length Re
# values is treated as ONE early-growth observation anchored to the
# LAST week of its window (the most recent, most informative week of
# that estimate): week = onset_week + window_weeks - 1. This yields up
# to 3 weekly-anchored Re observations per eligible UF, all drawn from
# the exact frozen Re_early_{4,6,8} posteriors (median/q025/q975 from
# the master table; full posterior draws from
# brazil_chik_major_wave_early_re_posterior_draws.rds for Section 6).

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/40_renewal_model/15_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

eligibility <- read_csv(file.path(paths$table, "FIRST_EPISODE_CLIMATE_ELIGIBILITY.csv"), show_col_types = FALSE) |>
  filter(eligible_for_climate_analysis == "YES")
master <- read_csv(file.path(paths$table, "brazil_chik_wave_analysis_master.csv"), show_col_types = FALSE) |>
  mutate(onset_week = as.Date(onset_week))

first_ep <- master |> filter(wave_order == 1, state %in% eligibility$state) |>
  select(state, wave_id, onset_week,
         Re_early_4_median, Re_early_4_q025, Re_early_4_q975,
         Re_early_6_median, Re_early_6_q025, Re_early_6_q975,
         Re_early_8_median, Re_early_8_q025, Re_early_8_q975)

weekly_re <- bind_rows(lapply(c(4L, 6L, 8L), function(w) {
  prefix <- paste0("Re_early_", w)
  first_ep |>
    transmute(UF = state, episode_id = wave_id, window_weeks = w,
              week = onset_week + (w - 1L) * 7L,
              Re_median = .data[[paste0(prefix, "_median")]],
              Re_lower = .data[[paste0(prefix, "_q025")]],
              Re_upper = .data[[paste0(prefix, "_q975")]])
})) |>
  filter(!is.na(Re_median)) |>
  arrange(UF, week)

message("=== FIRST_EPISODE_WEEKLY_RE: rows per UF ===")
print(as.data.frame(weekly_re |> count(UF, name = "n_weekly_re_obs")))
message(sprintf("\nTotal usable weekly-Re observations: %d (from %d eligible UFs)", nrow(weekly_re), n_distinct(weekly_re$UF)))

write_csv(weekly_re, file.path(paths$table, "FIRST_EPISODE_WEEKLY_RE.csv"))
message("[saved] ", file.path(paths$table, "FIRST_EPISODE_WEEKLY_RE.csv"))
