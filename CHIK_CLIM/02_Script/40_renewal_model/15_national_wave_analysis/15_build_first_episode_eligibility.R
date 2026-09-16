# Brazil first-epidemic climate-transmission pilot -- Section 1: identify
# the chronologically FIRST qualifying major epidemic in every one of the
# 27 UFs, using the ALREADY-FROZEN national wave census + early-Re outputs
# (wave_order == 1 among major_epidemic_primary waves). Eligibility for
# the climate analysis reuses the EXISTING re6_analysis_usable gate
# (re6_prefit_eligible & re6_fit_pass) -- the same objective rule already
# established for early-Re usability -- no new eligibility rule invented.

required_packages <- c("here", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/40_renewal_model/15_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

ALL_27_UFS <- c("AC","AL","AM","AP","BA","CE","DF","ES","GO","MA","MG","MS","MT","PA","PB",
                "PR","PE","PI","RJ","RN","RS","RO","RR","SC","SP","SE","TO")
stopifnot(length(ALL_27_UFS) == 27)

master <- read_csv(file.path(paths$table, "brazil_chik_wave_analysis_master.csv"), show_col_types = FALSE) |>
  mutate(onset_week = as.Date(onset_week), end_week = as.Date(end_week), peak_week = as.Date(peak_week))
detail <- read_csv(file.path(paths$table, "brazil_chik_major_wave_early_re.csv"), show_col_types = FALSE) |>
  select(state, wave_id, peak_cases, reason_not_major, reason_re6_not_eligible, re6_fit_status)

first_ep <- master |> filter(wave_order == 1) |>
  left_join(detail, by = c("state", "wave_id"))

eligibility <- tibble(state = ALL_27_UFS) |>
  left_join(
    first_ep |> transmute(
      state, first_episode_id = wave_id, episode_start = onset_week, episode_end = end_week,
      total_reported_cases = total_cases, peak_week, peak_cases,
      eligible_for_climate_analysis = if_else(re6_analysis_usable, "YES", "NO"),
      reason_for_exclusion = if_else(re6_analysis_usable, NA_character_, reason_re6_not_eligible)
    ),
    by = "state"
  ) |>
  mutate(
    eligible_for_climate_analysis = coalesce(eligible_for_climate_analysis, "NO"),
    reason_for_exclusion = if_else(
      is.na(first_episode_id), "no qualifying major epidemic wave in the national wave census (wave_order==1 does not exist for this UF)",
      reason_for_exclusion
    )
  ) |>
  arrange(state)

message("=== First-episode climate eligibility, all 27 UFs ===")
print(as.data.frame(eligibility), digits = 3)
message(sprintf("\nEligible: %d / 27 UFs", sum(eligibility$eligible_for_climate_analysis == "YES")))

write_csv(eligibility, file.path(paths$table, "FIRST_EPISODE_CLIMATE_ELIGIBILITY.csv"))
message("[saved] ", file.path(paths$table, "FIRST_EPISODE_CLIMATE_ELIGIBILITY.csv"))
