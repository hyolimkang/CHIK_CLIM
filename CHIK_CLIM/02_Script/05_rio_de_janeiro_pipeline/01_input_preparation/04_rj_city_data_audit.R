# Rio de Janeiro CITY auxiliary analysis -- Sections 3-4: build and audit
# city weekly case series (muni6=330455, from the SAME validated national
# muni-week panel used everywhere else in this project) and city
# demographic inputs. Does NOT modify the RJ STATE model/data.

required_packages <- c("dplyr", "readr", "ggplot2", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/05_rio_de_janeiro_pipeline/tables/rio_de_janeiro_v4_9_replication/city")
figure_dir <- file.path(root, "03_Output/05_rio_de_janeiro_pipeline/figures/rio_de_janeiro_v4_9_replication/city")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

MUNI6 <- "330455"
DATE_START <- as.Date("2015-01-04")
DATE_END <- as.Date("2025-12-21") # same window as the state model

panel <- readRDS(file.path(root, "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds")) |> mutate(muni6 = as.character(muni6))
city_panel_full <- panel |> dplyr::filter(muni6 == MUNI6) |> arrange(week_start)
message("Rio de Janeiro CITY rows in national panel: ", nrow(city_panel_full),
        " (", as.character(min(city_panel_full$week_start)), " to ", as.character(max(city_panel_full$week_start)), ")")

city_cases <- city_panel_full |> dplyr::filter(week_start >= DATE_START, week_start <= DATE_END) |>
  select(week_start, week_of_year, t, year, cases_confirmed, population) |> rename(cases = cases_confirmed)

# ---- Audit ----
weeks_expected <- seq(DATE_START, DATE_END, by = "week")
missing_weeks <- setdiff(as.character(weeks_expected), as.character(city_cases$week_start))
dup_weeks <- city_cases$week_start[duplicated(city_cases$week_start)]

audit <- tibble(
  muni6 = MUNI6, n_weeks = nrow(city_cases), n_weeks_expected = length(weeks_expected),
  week_continuity_ok = identical(city_cases$week_start, weeks_expected),
  n_missing_weeks = length(missing_weeks), n_duplicate_weeks = length(dup_weeks),
  n_na_cases = sum(is.na(city_cases$cases)), n_negative_cases = sum(city_cases$cases < 0, na.rm = TRUE),
  total_cases = sum(city_cases$cases, na.rm = TRUE), max_weekly_cases = max(city_cases$cases, na.rm = TRUE),
  max_weekly_cases_week = city_cases$week_start[which.max(city_cases$cases)],
  population_min = min(city_cases$population, na.rm = TRUE), population_max = max(city_cases$population, na.rm = TRUE)
)
message("\n=== Rio de Janeiro CITY weekly case audit ===")
print(as.data.frame(audit))
write_csv(audit, file.path(table_dir, "RJ_CITY_case_audit.csv"))
write_csv(city_cases, file.path(table_dir, "RJ_CITY_weekly_cases_clean.csv"))

if (!audit$week_continuity_ok || audit$n_na_cases > 0 || audit$n_negative_cases > 0) {
  stop("STOP: unresolved RJ CITY data anomalies -- do not proceed.")
}
message("\nCONFIRMED: no missing weeks, no duplicates, no NA/negative cases.")

annual <- city_cases |> group_by(year) |>
  summarise(total_cases = sum(cases), mean_population = mean(population),
            incidence_per_100k = total_cases / mean(mean_population) * 1e5,
            peak_weekly_cases = max(cases), peak_week = week_start[which.max(cases)], .groups = "drop")
message("\n=== Annual summary (Rio de Janeiro CITY) ===")
print(as.data.frame(annual), digits = 4)

theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- ggplot(city_cases, aes(week_start, cases)) + geom_line(colour = "#315A7D", linewidth = 0.4) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Rio de Janeiro CITY: weekly reported chikungunya cases, 2015-2025",
       subtitle = sprintf("Total = %s cases (muni6=%s)", format(sum(city_cases$cases), big.mark = ","), MUNI6),
       x = NULL, y = "Weekly reported cases") + theme_v4
ggsave(file.path(figure_dir, "RJ_CITY_cases_2015_2025.png"), p, width = 220, height = 110, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "RJ_CITY_cases_2015_2025.png"))

# ---- Section 4: demographic inputs (city vs state comparison for the audit) ----
demography_path <- file.path(root, "01_Data/riodejaneiro_city_weekly_demography_2015_2025.rds")
if (!file.exists(demography_path)) stop("City demography not yet built -- run 09e/10e in 02_Script/00_data_prep first.")
city_demography <- readRDS(demography_path) |> dplyr::filter(week_start >= DATE_START, week_start <= DATE_END)
write_csv(city_demography, file.path(table_dir, "RJ_CITY_demographic_inputs.csv"))
message("\n[saved] ", file.path(table_dir, "RJ_CITY_demographic_inputs.csv"))

accounting_error <- max(abs(city_demography$N_end - city_demography$N_start - city_demography$births +
                               city_demography$all_cause_deaths - city_demography$net_population_reconciliation))
message(sprintf("Demographic accounting max abs error: %.2e", accounting_error))

weekly_full <- city_cases |> select(-population) |> inner_join(city_demography, by = "week_start") |> arrange(week_start)
if (nrow(weekly_full) != nrow(city_cases)) stop("STOP: city case series and demography weeks do not align.")
saveRDS(weekly_full, file.path(table_dir, "riodejaneiro_city_weekly_input.rds"))
message("\n[saved] ", file.path(table_dir, "riodejaneiro_city_weekly_input.rds"), " (audited weekly input for Stan)")
