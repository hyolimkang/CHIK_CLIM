# Pernambuco spatial-turnover diagnostic -- explicit 2016 vs 2021/22 hotspot
# comparison (user-requested addition, beyond the Bahia-diagnostic template).
#
# PE's major-wave census gives: PE_wave_03 (onset 2015-11-29, peak
# 2016-02-28) as the index/"2016" epidemic, and PE_wave_07 (onset
# 2021-02-07) + PE_wave_08 (onset 2022-01-16) as the two subsequent
# "2021/22" epidemics. This script maps all three on a shared colour scale,
# classifies municipalities by top-20%-incidence "hotspot" status in the
# 2016 wave vs the pooled 2021+2022 period, and reports the direct
# muni-level correlation between the two eras.

required_packages <- c("dplyr", "readr", "sf", "ggplot2", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(sf); library(ggplot2); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/spatial_turnover_diagnostic")

pe_panel <- readRDS(file.path(table_dir, "pe_municipality_week_panel.rds"))
burden <- read_csv(file.path(table_dir, "pe_municipality_wave_burden.csv"), show_col_types = FALSE,
                    col_types = cols(muni6 = col_character()))
waves <- read_csv(file.path(table_dir, "pe_wave_definitions.csv"), show_col_types = FALSE) |> arrange(start_week)

poly <- readRDS(file.path(root, "01_Data/ibge_muni_polygons.rds")) |>
  dplyr::filter(substr(muni6, 1, 2) == "26")

WAVE_2016 <- "PE_wave_03"
WAVES_2021_22 <- c("PE_wave_07", "PE_wave_08")
stopifnot(all(c(WAVE_2016, WAVES_2021_22) %in% waves$wave_id))

message("=== Waves used for the 2016 vs 2021/22 comparison ===")
print(as.data.frame(waves |> dplyr::filter(wave_id %in% c(WAVE_2016, WAVES_2021_22)) |>
                       select(wave_id, wave_order, start_week, peak_week, end_week, total_cases)))

# ---- Pooled 2021/22 municipality burden (sum cases across both waves, ----
#      recompute incidence from the pooled case/population totals) ---------
pooled_2021_22 <- burden |>
  dplyr::filter(wave_id %in% WAVES_2021_22) |>
  group_by(muni6, name_muni) |>
  summarise(reported_cases_iw = sum(reported_cases_iw),
            population_mean_iw = mean(population_mean_iw), # same muni population base across the two adjacent waves
            .groups = "drop") |>
  mutate(wave_id = "PE_2021_22_pooled", incidence_iw = reported_cases_iw / population_mean_iw * 100000)

wave_2016 <- burden |> dplyr::filter(wave_id == WAVE_2016) |>
  select(muni6, name_muni, reported_cases_iw, population_mean_iw, incidence_iw, wave_id)

three_panel <- bind_rows(
  burden |> dplyr::filter(wave_id %in% WAVES_2021_22) |> select(muni6, wave_id, incidence_iw),
  wave_2016 |> select(muni6, wave_id, incidence_iw)
)

# ---- Panel maps: 2016, 2021, 2022, shared colour scale --------------------
theme_map <- theme_void(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold"), plot.subtitle = element_text(size = 7.5, colour = "grey40"),
        strip.text = element_text(size = 8, face = "bold"), legend.position = "bottom",
        legend.key.width = grid::unit(14, "pt"), legend.key.height = grid::unit(6, "pt"))

incidence_cap <- quantile(three_panel$incidence_iw, .98, na.rm = TRUE)
map_labels <- c(PE_wave_03 = "2016 epidemic (PE_wave_03)", PE_wave_07 = "2021 epidemic (PE_wave_07)", PE_wave_08 = "2022 epidemic (PE_wave_08)")

map_data <- poly |>
  left_join(three_panel, by = "muni6") |>
  mutate(wave_label = factor(map_labels[wave_id], levels = map_labels[c(WAVE_2016, WAVES_2021_22)]),
         incidence_capped = pmin(incidence_iw, incidence_cap))

p_three <- ggplot(map_data) +
  geom_sf(aes(fill = incidence_capped), colour = "white", linewidth = 0.05) +
  facet_wrap(~wave_label, ncol = 3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Incidence /100k\n(capped at 98th pct)", na.value = "grey90") +
  labs(title = "Pernambuco: 2016 index epidemic vs 2021 and 2022 epidemics (municipality level)",
       subtitle = sprintf("Shared colour scale across all three panels (capped at %.0f/100k, 98th pct)", incidence_cap)) +
  theme_map

ggsave(file.path(figure_dir, "pe_2016_vs_2021_2022_incidence_maps.png"), p_three, width = 260, height = 130, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pe_2016_vs_2021_2022_incidence_maps.png"))

# ---- Hotspot-status transition map: top-20% incidence in 2016 vs pooled 2021/22 ----
top20 <- function(df) {
  n <- ceiling(nrow(df) * 0.20)
  df |> arrange(desc(incidence_iw)) |> slice_head(n = n) |> pull(muni6)
}
hotspot_2016 <- top20(wave_2016)
hotspot_2021_22 <- top20(pooled_2021_22)

transition <- tibble(muni6 = poly$muni6) |>
  mutate(
    was_hotspot_2016 = muni6 %in% hotspot_2016,
    is_hotspot_2021_22 = muni6 %in% hotspot_2021_22,
    status = case_when(
      was_hotspot_2016 & is_hotspot_2021_22 ~ "Persistent hotspot (both eras)",
      was_hotspot_2016 & !is_hotspot_2021_22 ~ "2016-only hotspot (faded)",
      !was_hotspot_2016 & is_hotspot_2021_22 ~ "2021/22-only hotspot (new/emerging)",
      TRUE ~ "Not top-20% in either era"
    )
  )

write_csv(transition, file.path(table_dir, "pe_2016_vs_2021_22_hotspot_transition.csv"))
message("\n[saved] ", file.path(table_dir, "pe_2016_vs_2021_22_hotspot_transition.csv"))

message("\n=== Hotspot transition classification (top-20% incidence municipalities) ===")
print(as.data.frame(table(transition$status)))

hotspot_jaccard <- length(intersect(hotspot_2016, hotspot_2021_22)) / length(union(hotspot_2016, hotspot_2021_22))
message(sprintf("\nJaccard overlap of top-20%% hotspots, 2016 vs pooled 2021/22: %.3f", hotspot_jaccard))

map_transition <- poly |> left_join(transition, by = "muni6")

status_colours <- c("Persistent hotspot (both eras)" = "#D55E00",
                     "2016-only hotspot (faded)" = "#56B4E9",
                     "2021/22-only hotspot (new/emerging)" = "#009E73",
                     "Not top-20% in either era" = "grey85")

p_transition <- ggplot(map_transition) +
  geom_sf(aes(fill = status), colour = "white", linewidth = 0.08) +
  scale_fill_manual(values = status_colours, name = NULL) +
  labs(title = "Pernambuco: top-20% incidence hotspot status, 2016 index epidemic vs pooled 2021/22 epidemics",
       subtitle = sprintf("Jaccard overlap = %.3f. Clustering in the 'persistent' colour indicates hotspot stability; scattered 'new/emerging' colour indicates spatial turnover.", hotspot_jaccard)) +
  theme_map + theme(legend.position = "right")

ggsave(file.path(figure_dir, "pe_2016_vs_2021_2022_hotspot_transition_map.png"), p_transition, width = 200, height = 180, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pe_2016_vs_2021_2022_hotspot_transition_map.png"))

# ---- Direct muni-level scatter: 2016 incidence vs pooled 2021/22 incidence ----
scatter_df <- inner_join(
  wave_2016 |> select(muni6, name_muni, incidence_2016 = incidence_iw),
  pooled_2021_22 |> select(muni6, incidence_2021_22 = incidence_iw),
  by = "muni6"
)
rho_2016_vs_2021_22 <- cor(scatter_df$incidence_2016, scatter_df$incidence_2021_22, method = "spearman")
message(sprintf("\nSpearman correlation, municipality incidence: 2016 vs pooled 2021/22 = %.3f", rho_2016_vs_2021_22))

theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
p_scatter <- ggplot(scatter_df, aes(incidence_2016 + 1, incidence_2021_22 + 1)) +
  geom_point(alpha = 0.4, size = 1.2, colour = "#315A7D") +
  geom_smooth(method = "loess", colour = "#D55E00", se = TRUE, linewidth = 0.7) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Pernambuco municipality incidence: 2016 index epidemic vs pooled 2021/22 epidemics",
       subtitle = sprintf("Spearman rho = %.3f", rho_2016_vs_2021_22),
       x = "2016 (PE_wave_03) incidence per 100k (+1, log scale)",
       y = "Pooled 2021+2022 incidence per 100k (+1, log scale)") +
  theme_v4

ggsave(file.path(figure_dir, "pe_2016_vs_2021_2022_scatter.png"), p_scatter, width = 180, height = 150, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pe_2016_vs_2021_2022_scatter.png"))

summary_out <- tibble(
  jaccard_top20_2016_vs_2021_22 = hotspot_jaccard,
  spearman_2016_vs_2021_22 = rho_2016_vs_2021_22,
  n_persistent_hotspot = sum(transition$status == "Persistent hotspot (both eras)"),
  n_2016_only_hotspot = sum(transition$status == "2016-only hotspot (faded)"),
  n_2021_22_only_hotspot = sum(transition$status == "2021/22-only hotspot (new/emerging)"),
  n_neither = sum(transition$status == "Not top-20% in either era")
)
write_csv(summary_out, file.path(table_dir, "pe_2016_vs_2021_22_summary.csv"))
message("\n[saved] ", file.path(table_dir, "pe_2016_vs_2021_22_summary.csv"))
