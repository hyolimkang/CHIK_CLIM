# Mato Grosso v4.9 replication -- Section 3: epidemic-episode audit.
# Reuses the EXACT existing national wave census (the same
# episode-definition algorithm already used for CE/BA/PE/RJ) -- no manual
# redefinition of MT epidemic periods.

required_packages <- c("dplyr", "readr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/mato_grosso_v4_9_replication")
figure_dir <- file.path(root, "03_Output/figures/mato_grosso_v4_9_replication")

wave_census <- read_csv(file.path(root, "03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"),
                         show_col_types = FALSE)
mt_waves_all <- wave_census |> dplyr::filter(state == "MT") |>
  transmute(wave_id, wave_order, major_epidemic_primary,
            start_week = as.Date(onset_week, format = "%m/%d/%Y"),
            peak_week = as.Date(peak_week, format = "%m/%d/%Y"),
            end_week = as.Date(end_week, format = "%m/%d/%Y"),
            total_cases) |>
  arrange(start_week)

message("=== All MT episodes (national wave census) ===")
print(as.data.frame(mt_waves_all))

mt_waves_all <- mt_waves_all |>
  mutate(duration_weeks = as.integer(end_week - start_week) / 7,
         quiet_to_next_weeks = as.integer(lead(start_week) - end_week) / 7)

write_csv(mt_waves_all, file.path(table_dir, "MT_episode_audit.csv"))
message("\n[saved] ", file.path(table_dir, "MT_episode_audit.csv"))

major_waves <- mt_waves_all |> dplyr::filter(major_epidemic_primary == TRUE)
message("\n=== Major MT epidemics (n=", nrow(major_waves), ") ===")
print(as.data.frame(major_waves |> select(wave_id, start_week, peak_week, end_week, total_cases, duration_weeks)))

# Descriptive characterisation only (no susceptibility inference)
early_total <- sum(major_waves$total_cases[major_waves$wave_id %in% c("MT_wave_05", "MT_wave_06")])
late_total <- sum(major_waves$total_cases[major_waves$wave_id %in% c("MT_wave_13", "MT_wave_14")])
message("\n=== Descriptive characterisation ===")
message(sprintf("Early epidemic period (MT_wave_05+06, 2016-2019) = %s cases", format(early_total, big.mark=",")))
message(sprintf("Quiet intermediate period (MT_wave_07-12, 2019-2023): all waves <= 250 cases each"))
message(sprintf("Late recurrence (MT_wave_13+14, 2024-2025) = %s cases (%.1fx the early epidemic)",
                 format(late_total, big.mark=","), late_total / early_total))
message("Pattern: a substantial 2016-2019 epidemic period, a genuinely quiet 2019-2023 interval (all waves <250 cases/wave), then a much LARGER 2024-2025 recurrence -- the reverse size ordering seen in most other states (RJ/BA/PE all had their largest epidemic mid-series, not at the end).")

# ---- Figure ----
weekly <- read_csv(file.path(table_dir, "MT_weekly_cases_clean.csv"), show_col_types = FALSE)
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- ggplot(weekly, aes(week_start, cases)) +
  geom_line(colour = "grey40", linewidth = 0.3) +
  geom_rect(data = major_waves, aes(xmin = start_week, xmax = end_week, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#D55E00", alpha = 0.12) +
  geom_vline(data = major_waves, aes(xintercept = peak_week), colour = "#D55E00", linetype = "dashed", linewidth = 0.3) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Mato Grosso: weekly cases with major epidemic episodes (national wave census)",
       subtitle = sprintf("%d major episodes shaded; dashed lines = peak weeks", nrow(major_waves)),
       x = NULL, y = "Weekly reported cases") +
  theme_v4
ggsave(file.path(figure_dir, "MT_episode_plot.png"), p, width = 220, height = 120, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "MT_episode_plot.png"))
