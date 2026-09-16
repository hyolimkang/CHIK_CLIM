# Brazil first-epidemic climate-transmission pilot -- Sections 3-4: join
# the ALREADY-BUILT UF-week climate table (brazil_chik_uf_weekly_climate.csv
# -- population-weighted mean Tmean [deg C] and population-weighted total
# weekly PRCP [mm], built by 09_build_brazil_uf_weekly_climate.R, reused
# unmodified) onto each weekly-Re observation using three PRE-SPECIFIED
# lag-window sets (primary + 2 sensitivity, fixed BEFORE inspecting
# results, per instruction -- no lag search).

required_packages <- c("here", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/40_renewal_model/15_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

weekly_re <- read_csv(file.path(paths$table, "FIRST_EPISODE_WEEKLY_RE.csv"), show_col_types = FALSE) |>
  mutate(week = as.Date(week))
climate <- read_csv(file.path(paths$table, "brazil_chik_uf_weekly_climate.csv"), show_col_types = FALSE) |>
  mutate(week_start = as.Date(week_start))

# Pre-specified lag-window sets: offsets relative to the Re(t) anchor week
# (0 = anchor week itself, negative = weeks before). All fixed before
# inspecting any result.
LAG_SPECS <- list(
  primary     = list(temp_offsets = -2:0, precip_offsets = -6:-2),
  sensitivity_A = list(temp_offsets = -1:0, precip_offsets = -4:-1),
  sensitivity_B = list(temp_offsets = -3:-1, precip_offsets = -8:-3)
)

climate_window_summary <- function(climate_uf, anchor_week, temp_offsets, precip_offsets) {
  temp_weeks <- anchor_week + temp_offsets * 7L
  precip_weeks <- anchor_week + precip_offsets * 7L
  temp_rows <- climate_uf |> filter(week_start %in% temp_weeks)
  precip_rows <- climate_uf |> filter(week_start %in% precip_weeks)
  list(
    temp = if (nrow(temp_rows) == length(temp_offsets)) mean(temp_rows$Tmean, na.rm = FALSE) else NA_real_,
    precip = if (nrow(precip_rows) == length(precip_offsets)) sum(precip_rows$PRCP, na.rm = FALSE) else NA_real_
  )
}

build_dataset <- function(spec_name, spec) {
  bind_rows(lapply(seq_len(nrow(weekly_re)), function(i) {
    row <- weekly_re[i, ]
    climate_uf <- climate |> filter(state == row$UF)
    w <- climate_window_summary(climate_uf, row$week, spec$temp_offsets, spec$precip_offsets)
    tibble(UF = row$UF, episode_id = row$episode_id, window_weeks = row$window_weeks, week = row$week,
           Re_median = row$Re_median, Re_lower = row$Re_lower, Re_upper = row$Re_upper,
           lag_spec = spec_name, temperature = w$temp, precipitation = w$precip)
  }))
}

all_specs <- bind_rows(lapply(names(LAG_SPECS), function(nm) build_dataset(nm, LAG_SPECS[[nm]])))

n_missing <- all_specs |> filter(is.na(temperature) | is.na(precipitation)) |> nrow()
message(sprintf("Rows with missing climate (window extends outside available climate coverage): %d / %d", n_missing, nrow(all_specs)))
if (n_missing > 0) {
  message("Excluded rows:")
  print(as.data.frame(all_specs |> filter(is.na(temperature) | is.na(precipitation)) |> select(UF, episode_id, window_weeks, week, lag_spec)))
}
all_specs_clean <- all_specs |> filter(!is.na(temperature), !is.na(precipitation))

# Primary-spec training dataset (Section 5 GAM uses this)
training <- all_specs_clean |> filter(lag_spec == "primary") |>
  select(UF, episode_id, window_weeks, week, Re_median, Re_lower, Re_upper, temperature, precipitation)

message("\n=== Primary-spec training dataset summary ===")
message(sprintf("N = %d rows, %d UFs", nrow(training), n_distinct(training$UF)))
message(sprintf("Temperature range: [%.2f, %.2f] deg C", min(training$temperature), max(training$temperature)))
message(sprintf("Precipitation range: [%.2f, %.2f] mm (5-week cumulative)", min(training$precipitation), max(training$precipitation)))
print(as.data.frame(training), digits = 3)

write_csv(all_specs_clean, file.path(paths$table, "FIRST_EPISODE_CLIMATE_ALL_SPECS.csv"))
write_csv(training, file.path(paths$table, "first_epidemic_climate_training_data.csv"))
message("\n[saved] ", file.path(paths$table, "FIRST_EPISODE_CLIMATE_ALL_SPECS.csv"))
message("[saved] ", file.path(paths$table, "first_epidemic_climate_training_data.csv"))
