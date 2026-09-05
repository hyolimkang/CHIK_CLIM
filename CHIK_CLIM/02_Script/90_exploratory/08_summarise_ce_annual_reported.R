# summarise_ce_annual_reported.R
# Annual reported chikungunya cases for Ceará from muni-week panel.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
})

YEAR_START <- 2015L
YEAR_END   <- 2025L
UF_CODE    <- "23"
CASE_VAR   <- "cases_confirmed"

panel_paths <- c(
  "01_Data/chik_brazil_muni_week_2015_2024.rds",
  "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_MORBID/CHIK_MORBID/01_Data/chik_brazil_muni_week_2015_2024.rds",
  "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_MORBID/CHIK_MORBID/01_Data/chik_brazil_muni_week_2015_2024.csv"
)
PANEL_PATH <- panel_paths[file.exists(panel_paths)][1]
if (is.na(PANEL_PATH)) stop("Muni-week panel not found.")

pop_paths <- c(
  "01_Data/ibge_pop_muni_year_2015_2024.rds",
  "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_MORBID/CHIK_MORBID/01_Data/ibge_pop_muni_year_2015_2024.csv"
)
POP_PATH <- pop_paths[file.exists(pop_paths)][1]

message("Panel: ", PANEL_PATH)
if (grepl("\\.rds$", PANEL_PATH, ignore.case = TRUE)) {
  panel <- readRDS(PANEL_PATH)
} else {
  panel <- read_csv(PANEL_PATH, show_col_types = FALSE)
}

if (!CASE_VAR %in% names(panel)) {
  stop("Column not found: ", CASE_VAR)
}

ce_panel <- panel |>
  mutate(
    muni6 = sprintf("%06d", as.integer(muni6)),
    week_start = as.Date(week_start),
    year = year(week_start)
  ) |>
  filter(substr(muni6, 1, 2) == UF_CODE)

data_max_year <- max(ce_panel$year, na.rm = TRUE)
data_max_date <- max(ce_panel$week_start, na.rm = TRUE)

ce_annual <- ce_panel |>
  filter(year >= YEAR_START, year <= YEAR_END) |>
  group_by(year) |>
  summarise(
    cases_confirmed = sum(.data[[CASE_VAR]], na.rm = TRUE),
    cases_notified  = if ("cases_notified" %in% names(ce_panel)) {
      sum(cases_notified, na.rm = TRUE)
    } else {
      NA_real_
    },
    n_weeks = n_distinct(week_start),
    .groups = "drop"
  ) |>
  complete(year = YEAR_START:YEAR_END, fill = list(
    cases_confirmed = 0,
    cases_notified = 0,
    n_weeks = 0
  )) |>
  arrange(year)

if (!is.na(POP_PATH)) {
  message("Population: ", POP_PATH)
  if (grepl("\\.rds$", POP_PATH, ignore.case = TRUE)) {
    pop <- readRDS(POP_PATH)
  } else {
    pop <- read_csv(POP_PATH, show_col_types = FALSE)
  }
  ce_pop <- pop |>
    mutate(muni6 = sprintf("%06d", as.integer(muni6))) |>
    filter(substr(muni6, 1, 2) == UF_CODE) |>
    group_by(year) |>
    summarise(population = sum(population, na.rm = TRUE), .groups = "drop")

  ce_annual <- ce_annual |>
    left_join(ce_pop, by = "year") |>
    mutate(
      attack_rate_pct    = round(100 * cases_confirmed / population, 3),
      attack_rate_per100k = round(1e5 * cases_confirmed / population, 1)
    )
}

ce_annual <- ce_annual |>
  mutate(
    data_note = dplyr::case_when(
      year > data_max_year ~ "no data in panel",
      year == data_max_year ~ paste0("partial through ", data_max_date),
      TRUE ~ "complete year"
    )
  )

out_dir <- "03_Output/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "ce_annual_reported_cases_2015_2025.csv")
write_csv(ce_annual, out_csv)

cat("\nCeará reported cases (", CASE_VAR, ") by calendar year\n", sep = "")
cat("Panel date range: ", as.character(min(ce_panel$week_start)), " to ",
    as.character(data_max_date), "\n\n", sep = "")
print(as.data.frame(ce_annual), row.names = FALSE)

cat("\nTotal ", YEAR_START, "-", min(YEAR_END, data_max_year), ": ",
    format(sum(ce_annual$cases_confirmed[ce_annual$year <= data_max_year],
               na.rm = TRUE), big.mark = ","), "\n", sep = "")
cat("Saved: ", out_csv, "\n", sep = "")
