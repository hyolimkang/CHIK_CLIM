# Rio de Janeiro v4.9 replication -- Sections 1-3: state-level weekly case
# data construction + audit, reusing the exact SAME pipeline already
# validated for CE/BA/PE (build_state_weekly(), national muni-week panel,
# national wave census v2). No new algorithm invented.

required_packages <- c("dplyr", "readr", "ggplot2", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/40_renewal_model/00_shared/01_build_ceara_state_weekly.R"))

base_dir <- file.path(root, "02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication")
figure_dir <- file.path(root, "03_Output/figures/rio_de_janeiro_v4_9_replication")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

UF_CODE <- "33" # Rio de Janeiro
DATE_START <- as.Date("2015-01-04")
DATE_END <- as.Date("2025-12-21") # same window convention as CE/BA/PE

panel <- readRDS(file.path(root, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds"))
rj_muni_n <- length(unique(panel$muni6[substr(panel$muni6, 1, 2) == UF_CODE]))
message("RJ municipalities in national panel: ", rj_muni_n, " (official count: 92)")

rj_weekly_full <- build_state_weekly(panel, UF_CODE)
message("build_state_weekly() returned ", nrow(rj_weekly_full), " weeks, range ",
        as.character(min(rj_weekly_full$week_start)), " to ", as.character(max(rj_weekly_full$week_start)))

# ---- Audit: missing weeks, duplicates, NA, negative values, week-53 handling ----
expected_weeks <- seq(min(rj_weekly_full$week_start), max(rj_weekly_full$week_start), by = "week")
missing_weeks <- setdiff(as.character(expected_weeks), as.character(rj_weekly_full$week_start))
dup_weeks <- rj_weekly_full$week_start[duplicated(rj_weekly_full$week_start)]
n_na_cases <- sum(is.na(rj_weekly_full$cases))
n_negative <- sum(rj_weekly_full$cases < 0, na.rm = TRUE)

rj_weekly <- rj_weekly_full |> filter(week_start >= DATE_START, week_start <= DATE_END)
n_weeks_model <- nrow(rj_weekly)
expected_model_weeks <- length(seq(DATE_START, DATE_END, by = "week"))

audit <- tibble(
  n_municipalities = rj_muni_n,
  first_observation_date = min(rj_weekly_full$week_start),
  last_observation_date = max(rj_weekly_full$week_start),
  n_weeks_full_series = nrow(rj_weekly_full),
  n_weeks_expected_full = length(expected_weeks),
  n_missing_weeks = length(missing_weeks),
  n_duplicate_weeks = length(dup_weeks),
  n_na_cases = n_na_cases,
  n_negative_cases = n_negative,
  n_weeks_model_window = n_weeks_model,
  n_weeks_model_window_expected = expected_model_weeks,
  total_cases_2015_2025 = sum(rj_weekly$cases, na.rm = TRUE),
  max_weekly_cases = max(rj_weekly$cases, na.rm = TRUE),
  max_weekly_cases_week = rj_weekly$week_start[which.max(rj_weekly$cases)],
  population_min = min(rj_weekly$population, na.rm = TRUE),
  population_max = max(rj_weekly$population, na.rm = TRUE)
)
message("\n=== RJ weekly case series audit ===")
print(as.data.frame(audit))
if (length(missing_weeks) > 0) message("Missing weeks: ", paste(missing_weeks, collapse = ", "))
if (length(dup_weeks) > 0) message("Duplicate weeks: ", paste(as.character(dup_weeks), collapse = ", "))

write_csv(audit, file.path(table_dir, "RJ_weekly_cases_audit.csv"))
write_csv(rj_weekly, file.path(table_dir, "RJ_weekly_cases_clean.csv"))
message("\n[saved] ", file.path(table_dir, "RJ_weekly_cases_audit.csv"))
message("[saved] ", file.path(table_dir, "RJ_weekly_cases_clean.csv"))

if (length(missing_weeks) > 0 || length(dup_weeks) > 0 || n_na_cases > 0 || n_negative > 0) {
  stop("STOP: unresolved data-processing anomalies detected (missing/duplicate weeks, NA, or negative cases) -- do not proceed to model fitting until resolved.")
} else {
  message("\nCONFIRMED: no missing weeks, no duplicates, no NA/negative case values. Safe to proceed.")
}

# ---- Annual summary ----
annual <- rj_weekly |>
  mutate(year = as.integer(format(week_start, "%Y"))) |>
  group_by(year) |>
  summarise(total_cases = sum(cases), mean_population = mean(population),
            incidence_per_100k = total_cases / mean(mean_population) * 1e5,
            peak_weekly_cases = max(cases), peak_week = week_start[which.max(cases)],
            n_weeks_nonzero = sum(cases > 0), .groups = "drop")
message("\n=== Annual summary ===")
print(as.data.frame(annual), digits = 4)
write_csv(annual, file.path(table_dir, "RJ_case_audit.csv"))
message("[saved] ", file.path(table_dir, "RJ_case_audit.csv"))

# ---- Figure ----
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- ggplot(rj_weekly, aes(week_start, cases)) +
  geom_line(colour = "#315A7D", linewidth = 0.4) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Rio de Janeiro: weekly reported chikungunya cases, 2015-2025",
       subtitle = sprintf("Total = %s cases across %d weeks (build_state_weekly, UF=33)", format(sum(rj_weekly$cases), big.mark=","), n_weeks_model),
       x = NULL, y = "Weekly reported cases") +
  theme_v4
ggsave(file.path(figure_dir, "RJ_weekly_cases_2015_2025.png"), p, width = 220, height = 110, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "RJ_weekly_cases_2015_2025.png"))
