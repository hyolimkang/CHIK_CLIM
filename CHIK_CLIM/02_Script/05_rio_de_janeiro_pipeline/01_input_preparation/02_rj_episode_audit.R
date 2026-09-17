# Rio de Janeiro v4.9 replication -- Section 3: epidemic-episode audit.
# Reuses the EXACT existing national wave census (the same
# episode-definition algorithm already used for CE/BA/PE) -- no manual
# redefinition of RJ epidemic periods.

required_packages <- c("dplyr", "readr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication")
figure_dir <- file.path(root, "03_Output/figures/rio_de_janeiro_v4_9_replication")

wave_census <- read_csv(file.path(root, "03_Output/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"),
                         show_col_types = FALSE)
rj_waves_all <- wave_census |> dplyr::filter(state == "RJ") |>
  transmute(wave_id, wave_order, major_epidemic_primary,
            start_week = as.Date(onset_week, format = "%m/%d/%Y"),
            peak_week = as.Date(peak_week, format = "%m/%d/%Y"),
            end_week = as.Date(end_week, format = "%m/%d/%Y"),
            total_cases) |>
  arrange(start_week)

message("=== All RJ episodes (national wave census) ===")
print(as.data.frame(rj_waves_all))

rj_waves_all <- rj_waves_all |>
  mutate(duration_weeks = as.integer(end_week - start_week) / 7,
         quiet_to_next_weeks = as.integer(lead(start_week) - end_week) / 7)

write_csv(rj_waves_all, file.path(table_dir, "RJ_episode_audit.csv"))
message("\n[saved] ", file.path(table_dir, "RJ_episode_audit.csv"))

major_waves <- rj_waves_all |> dplyr::filter(major_epidemic_primary == TRUE)
message("\n=== Major RJ epidemics (n=", nrow(major_waves), ") ===")
print(as.data.frame(major_waves |> select(wave_id, start_week, peak_week, end_week, total_cases, duration_weeks)))

# Descriptive characterisation only (no susceptibility inference)
message("\n=== Descriptive characterisation ===")
message(sprintf("RJ_wave_07 (2018-12-30 to 2020-12-20) alone = %s cases (%.1f%% of the 2015-2025 case total 126,428)",
                 format(major_waves$total_cases[major_waves$wave_id == "RJ_wave_07"], big.mark=","),
                 100 * major_waves$total_cases[major_waves$wave_id == "RJ_wave_07"] / 126428))
message("Pattern: one dominant, prolonged 'mega-epidemic' (2018-2020) preceded by two small early waves (2015-16, 2016-17) and followed by recurrent but much smaller waves (2023-24, 2024-25) -- descriptively resembles a 'delayed emergence into one massive epidemic, then recurrent minor waves' pattern, NOT a regularly recurring multi-wave pattern like Bahia/Pernambuco.")

# ---- Figure ----
weekly <- read_csv(file.path(table_dir, "RJ_weekly_cases_clean.csv"), show_col_types = FALSE)
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
p <- ggplot(weekly, aes(week_start, cases)) +
  geom_line(colour = "grey40", linewidth = 0.3) +
  geom_rect(data = major_waves, aes(xmin = start_week, xmax = end_week, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#D55E00", alpha = 0.12) +
  geom_vline(data = major_waves, aes(xintercept = peak_week), colour = "#D55E00", linetype = "dashed", linewidth = 0.3) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Rio de Janeiro: weekly cases with major epidemic episodes (national wave census)",
       subtitle = sprintf("%d major episodes shaded; dashed lines = peak weeks", nrow(major_waves)),
       x = NULL, y = "Weekly reported cases") +
  theme_v4
ggsave(file.path(figure_dir, "RJ_episode_plot.png"), p, width = 220, height = 120, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "RJ_episode_plot.png"))
