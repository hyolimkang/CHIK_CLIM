# Bahia spatial-turnover diagnostic -- Section 6.
#
# Maps of municipality-level wave incidence (same colour scale across
# waves) using the existing official IBGE municipality geometry, plus a
# summary map of which wave produced the maximum incidence in each
# municipality. No allfoi used.

required_packages <- c("dplyr", "readr", "sf", "ggplot2", "viridis")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(sf); library(ggplot2); library(viridis) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
table_dir <- file.path(root, "03_Output/tables/bahia_spatial_turnover_diagnostic")
figure_dir <- file.path(root, "03_Output/figures/bahia_spatial_turnover_diagnostic")

burden <- read_csv(file.path(table_dir, "bahia_municipality_wave_burden.csv"), show_col_types = FALSE,
                    col_types = cols(muni6 = col_character()))
waves <- read_csv(file.path(table_dir, "bahia_wave_definitions.csv"), show_col_types = FALSE) |> arrange(start_week)

poly <- readRDS(file.path(root, "01_Data/ibge_muni_polygons.rds")) |>
  dplyr::filter(substr(muni6, 1, 2) == "29")

theme_map <- theme_void(base_size = 8) +
  theme(plot.title = element_text(size = 9, face = "bold"), plot.subtitle = element_text(size = 7.5, colour = "grey40"),
        strip.text = element_text(size = 7.5, face = "bold"), legend.position = "bottom",
        legend.key.width = grid::unit(14, "pt"), legend.key.height = grid::unit(6, "pt"))

# ---- Wave-specific incidence maps, SAME colour scale --------------------
map_data <- poly |> left_join(burden |> select(muni6, wave_id, incidence_iw), by = "muni6") |>
  mutate(wave_id = factor(wave_id, levels = waves$wave_id))

incidence_cap <- quantile(burden$incidence_iw, .98, na.rm = TRUE) # cap extreme outlier for a readable shared scale
map_data <- map_data |> mutate(incidence_capped = pmin(incidence_iw, incidence_cap))

p_maps <- ggplot(map_data) +
  geom_sf(aes(fill = incidence_capped), colour = "white", linewidth = 0.02) +
  facet_wrap(~wave_id, ncol = 3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Incidence /100k\n(capped at 98th pct)", na.value = "grey90") +
  labs(title = "Bahia municipality-level reported incidence by major epidemic wave",
       subtitle = sprintf("Same colour scale across waves (capped at %.0f/100k, 98th percentile) -- no allfoi used", incidence_cap)) +
  theme_map

ggsave(file.path(figure_dir, "bahia_wave_incidence_maps.png"), p_maps, width = 260, height = 300, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "bahia_wave_incidence_maps.png"))

# ---- Summary map: which wave produced the maximum incidence per muni ----
max_wave_df <- burden |>
  group_by(muni6) |>
  slice_max(incidence_iw, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(muni6, wave_id_max = wave_id, max_incidence = incidence_iw) |>
  mutate(wave_id_max = factor(wave_id_max, levels = waves$wave_id))

map_max <- poly |> left_join(max_wave_df, by = "muni6")

p_max <- ggplot(map_max) +
  geom_sf(aes(fill = wave_id_max), colour = "white", linewidth = 0.05) +
  scale_fill_viridis_d(option = "turbo", name = "Wave with max\nincidence", na.value = "grey90") +
  labs(title = "Bahia: which epidemic wave produced each municipality's peak incidence",
       subtitle = "Spatially clustered colours would indicate hotspot persistence; a scrambled pattern indicates hotspot migration") +
  theme_map + theme(legend.position = "right")

ggsave(file.path(figure_dir, "bahia_wave_max_incidence_summary_map.png"), p_max, width = 200, height = 180, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "bahia_wave_max_incidence_summary_map.png"))

message("\n=== Distribution of which wave produced max incidence (n municipalities) ===")
print(as.data.frame(table(max_wave_df$wave_id_max)))
