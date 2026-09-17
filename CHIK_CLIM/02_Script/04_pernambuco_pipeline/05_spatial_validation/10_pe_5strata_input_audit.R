# Pernambuco spatial model -- Section 3: pre-model spatial audit.
#
# Aggregate weekly reported cases and population for all 5 strata, report
# incidence, and quantify which stratum dominates each of the 9 previously
# defined PE major waves. Descriptive only -- no Stan fitting here.

required_packages <- c("dplyr", "readr", "tidyr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
spatial_table_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/tables/pernambuco_v4_9_replication/spatial_model")
spatial_figure_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/figures/pernambuco_v4_9_replication/spatial_model")
dir.create(spatial_figure_dir, recursive = TRUE, showWarnings = FALSE)
turnover_table_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")

stratum_lookup <- read_csv(file.path(spatial_table_dir, "PE_spatial_stratum_lookup.csv"), show_col_types = FALSE,
                            col_types = cols(muni6 = col_character()))
pe_panel <- readRDS(file.path(turnover_table_dir, "pe_municipality_week_panel.rds")) |>
  mutate(muni6 = as.character(muni6))
waves <- read_csv(file.path(turnover_table_dir, "pe_wave_definitions.csv"), show_col_types = FALSE) |> arrange(start_week)

panel_stratum <- pe_panel |>
  left_join(stratum_lookup |> select(muni6, final_model_stratum), by = "muni6")
if (any(is.na(panel_stratum$final_model_stratum))) {
  stop("STOP: some panel municipalities (case data present) have no stratum assignment -- investigate.")
}

# ---- Weekly cases / population / incidence by stratum ----------------------
weekly_stratum <- panel_stratum |>
  group_by(final_model_stratum, week_start) |>
  summarise(cases = sum(cases_confirmed, na.rm = TRUE), population = sum(population, na.rm = TRUE), .groups = "drop") |>
  mutate(incidence_100k = cases / population * 1e5)

write_csv(weekly_stratum, file.path(spatial_table_dir, "PE_5strata_weekly_cases.csv"))
message("[saved] ", file.path(spatial_table_dir, "PE_5strata_weekly_cases.csv"))

# Confirm summing the 5 strata reproduces the exact state series
state_check <- weekly_stratum |> group_by(week_start) |> summarise(cases_sum = sum(cases), .groups = "drop")
official <- pe_panel |> group_by(week_start) |> summarise(cases_official = sum(cases_confirmed, na.rm = TRUE), .groups = "drop")
recon <- inner_join(state_check, official, by = "week_start")
max_diff <- max(abs(recon$cases_sum - recon$cases_official))
message(sprintf("Sum of 5-strata weekly cases vs state series: max abs diff = %d across %d weeks", max_diff, nrow(recon)))
if (max_diff > 0) stop("STOP: 5-strata aggregation does not reproduce the state series exactly.")

stratum_pop_2025 <- weekly_stratum |> dplyr::filter(week_start == max(week_start)) |> select(final_model_stratum, population)
message("\n=== Stratum population (final week) ===")
print(as.data.frame(stratum_pop_2025))

# ---- Contribution of each stratum to every PE wave -------------------------
wave_contrib <- bind_rows(lapply(seq_len(nrow(waves)), function(i) {
  w <- waves[i, ]
  panel_stratum |>
    dplyr::filter(week_start >= w$start_week, week_start <= w$end_week) |>
    group_by(final_model_stratum) |>
    summarise(cases_in_wave = sum(cases_confirmed, na.rm = TRUE), .groups = "drop") |>
    mutate(wave_id = w$wave_id, wave_order = w$wave_order)
})) |>
  group_by(wave_id) |>
  mutate(total_wave_cases = sum(cases_in_wave), share_of_wave = cases_in_wave / total_wave_cases) |>
  ungroup() |>
  arrange(wave_order, desc(share_of_wave))

write_csv(wave_contrib, file.path(spatial_table_dir, "PE_5strata_wave_contributions.csv"))
message("\n[saved] ", file.path(spatial_table_dir, "PE_5strata_wave_contributions.csv"))

message("\n=== Dominant stratum per wave (share of that wave's cases) ===")
dominant <- wave_contrib |> group_by(wave_id) |> slice_max(share_of_wave, n = 1) |> ungroup() |>
  select(wave_id, wave_order, final_model_stratum, cases_in_wave, share_of_wave)
print(as.data.frame(dominant), digits = 3)

# Explicit answers for 2016 / 2021 / 2022 / later recurrence waves
message("\n=== Full stratum breakdown, key waves ===")
key_waves <- c("PE_wave_03" = "2016 index epidemic", "PE_wave_07" = "2021 epidemic",
               "PE_wave_08" = "2022 epidemic", "PE_wave_09" = "2022/23 recurrence",
               "PE_wave_10" = "2023/24 recurrence", "PE_wave_11" = "2024/25 recurrence")
for (wid in names(key_waves)) {
  message(sprintf("\n-- %s (%s) --", wid, key_waves[wid]))
  print(as.data.frame(wave_contrib |> dplyr::filter(wave_id == wid) |>
                         select(final_model_stratum, cases_in_wave, share_of_wave) |>
                         arrange(desc(share_of_wave))), digits = 3)
}

# ---- Figure: 5 weekly trajectories, common absolute-case and incidence scales ----
stratum_levels <- c("1_Recife", "2_Metropolitana_remainder", "3_Agreste", "4_Sertao", "5_Vale_Sao_Francisco_Araripe")
stratum_colours <- c("1_Recife" = "#D55E00", "2_Metropolitana_remainder" = "#E69F00",
                      "3_Agreste" = "#009E73", "4_Sertao" = "#0072B2", "5_Vale_Sao_Francisco_Araripe" = "#CC79A7")
weekly_stratum <- weekly_stratum |> mutate(final_model_stratum = factor(final_model_stratum, levels = stratum_levels))

theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

p_cases <- ggplot(weekly_stratum, aes(week_start, cases, colour = final_model_stratum)) +
  geom_line(linewidth = 0.4) +
  scale_colour_manual(values = stratum_colours, name = NULL) +
  labs(title = "Pernambuco 5-strata weekly reported cases (common absolute scale)", x = NULL, y = "Weekly cases") +
  theme_v4 + theme(legend.position = "top")

p_incidence <- ggplot(weekly_stratum, aes(week_start, incidence_100k, colour = final_model_stratum)) +
  geom_line(linewidth = 0.4) +
  scale_colour_manual(values = stratum_colours, name = NULL) +
  labs(title = "Pernambuco 5-strata weekly incidence per 100k (common incidence scale)", x = NULL, y = "Incidence per 100k") +
  theme_v4 + theme(legend.position = "none")

combined <- p_cases / p_incidence
ggsave(file.path(spatial_figure_dir, "PE_5strata_weekly_cases.png"), combined, width = 240, height = 220, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(spatial_figure_dir, "PE_5strata_weekly_cases.png"))

# ---- Write the required markdown summary -----------------------------------
md_lines <- c(
  "# Pernambuco 5-Strata Input Audit",
  "",
  "Descriptive pre-model audit only. No Stan fitting.",
  "",
  "## Stratum population (final observed week, 2025-12-21)",
  "",
  "| Stratum | Population |",
  "|---|---|",
  sprintf("| %s | %s |", stratum_pop_2025$final_model_stratum, format(stratum_pop_2025$population, big.mark = ",")),
  "",
  sprintf("Confirmed: summing the 5 strata's weekly cases reproduces the exact PE state-level series (max abs diff = %d across %d weeks).", max_diff, nrow(recon)),
  "",
  "## Dominant stratum per wave",
  "",
  "| Wave | Dominant stratum | Cases in wave | Share of wave |",
  "|---|---|---|---|",
  sprintf("| %s | %s | %d | %.1f%% |", dominant$wave_id, dominant$final_model_stratum, dominant$cases_in_wave, 100 * dominant$share_of_wave),
  "",
  "## Key-wave stratum breakdown"
)
for (wid in names(key_waves)) {
  sub <- wave_contrib |> dplyr::filter(wave_id == wid) |> arrange(desc(share_of_wave))
  md_lines <- c(md_lines, "", sprintf("### %s (%s)", wid, key_waves[wid]), "",
                "| Stratum | Cases | Share |", "|---|---|---|",
                sprintf("| %s | %d | %.1f%% |", sub$final_model_stratum, sub$cases_in_wave, 100 * sub$share_of_wave))
}
md_lines <- c(md_lines, "",
  "## Interpretation",
  "",
  "See PE_5strata_weekly_cases.png: the two-panel figure (common absolute-case",
  "scale on top, common incidence scale below) shows whether PE's recurrent",
  "state-level trajectory is visibly composed of asynchronous regional",
  "epidemics, or whether all 5 strata rise and fall together. See the",
  "dominant-stratum and key-wave tables above for the quantitative answer.")
pe_reports_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/reports/pernambuco_v4_9_replication")
dir.create(pe_reports_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(md_lines, file.path(pe_reports_dir, "PE_5strata_input_audit.md"))
message("\n[saved] ", file.path(pe_reports_dir, "PE_5strata_input_audit.md"))
