# Bahia global-q local-consistency audit -- Section 6.
# Maps of crude case-implied cumulative attack fraction A_i(2025 | q) for
# representative q values, using the same municipality geometry as the
# spatial-turnover diagnostic. Plausibility diagnostic only.

required_packages <- c("dplyr", "readr", "sf", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(sf); library(ggplot2) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/bahia_global_q_local_consistency_audit")
figure_dir <- file.path(root, "03_Output/figures/bahia_global_q_local_consistency_audit")

implied <- read_csv(file.path(table_dir, "bahia_municipality_implied_attack_by_q.csv"), show_col_types = FALSE,
                     col_types = cols(muni6 = col_character()))
poly <- readRDS(file.path(root, "01_Data/ibge_muni_polygons.rds")) |> filter(substr(muni6, 1, 2) == "29")

theme_map <- theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8, colour = "grey40"),
        legend.position = "bottom", legend.key.width = grid::unit(16, "pt"), legend.key.height = grid::unit(6, "pt"))

target_checkpoint <- as.Date("2025-12-21")
target_qs <- c(0.05, 0.136, 0.359)

for (q_val in target_qs) {
  df <- implied |> filter(checkpoint == target_checkpoint, q == q_val)
  map_df <- poly |> left_join(df |> select(muni6, A_i), by = "muni6")
  p <- ggplot(map_df) +
    geom_sf(aes(fill = pmin(A_i, 1)), colour = "white", linewidth = 0.02) +
    scale_fill_viridis_c(option = "inferno", direction = -1, limits = c(0, 1), name = "A_i (capped at 1.0)", na.value = "grey90") +
    labs(title = sprintf("Bahia: crude case-implied cumulative attack fraction A_i(2025 | q=%.3f)", q_val),
         subtitle = "CRUDE, assumes constant q and ignores transmission dynamics -- NOT a Stan estimate") +
    theme_map
  fname <- sprintf("bahia_municipality_implied_attack_q%s.png", gsub("\\.", "", sprintf("%.3f", q_val)))
  ggsave(file.path(figure_dir, fname), p, width = 180, height = 170, units = "mm", dpi = 300, bg = "white")
  message("[saved] ", file.path(figure_dir, fname))
}

# ---- A_i > 1 flag map, all three q values in one panel ---------------------
flag_df <- implied |> filter(checkpoint == target_checkpoint, q %in% target_qs) |>
  mutate(flag = A_i > 1, q_label = sprintf("q = %.3f", q))
map_flag <- poly |> left_join(flag_df |> select(muni6, q_label, flag), by = "muni6", relationship = "many-to-many")

p_flag <- ggplot(map_flag) +
  geom_sf(aes(fill = flag), colour = "white", linewidth = 0.02) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey92"), na.value = "grey92",
                     labels = c(`TRUE` = "A_i > 1 (flagged)", `FALSE` = "A_i <= 1"), name = NULL) +
  facet_wrap(~q_label) +
  labs(title = "Bahia: municipalities where crude case-implied attack fraction exceeds 1.0 (2025)",
       subtitle = "A falsification diagnostic under the lifelong-first-infection interpretation -- flagged, not auto-excluded") +
  theme_map + theme(strip.text = element_text(face = "bold"))
ggsave(file.path(figure_dir, "bahia_municipality_A_gt1_flag_map.png"), p_flag, width = 260, height = 140, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "bahia_municipality_A_gt1_flag_map.png"))
