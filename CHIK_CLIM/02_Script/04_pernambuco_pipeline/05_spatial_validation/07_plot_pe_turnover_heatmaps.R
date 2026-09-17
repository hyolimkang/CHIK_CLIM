# Pernambuco spatial-turnover diagnostic -- Section 7 (heatmaps) and the
# Section 4C/4D case-contributor-turnover figure.

required_packages <- c("dplyr", "readr", "ggplot2", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/spatial_turnover_diagnostic")

metrics <- read_csv(file.path(table_dir, "pe_wave_turnover_metrics.csv"), show_col_types = FALSE)
waves <- read_csv(file.path(table_dir, "pe_wave_definitions.csv"), show_col_types = FALSE) |> arrange(start_week)
wave_levels <- waves$wave_id

pairwise <- metrics |> dplyr::filter(comparison_type == "all_pairwise")
# Symmetrise for the heatmap (undirected pairwise comparison) and add the diagonal.
pairwise_sym <- bind_rows(
  pairwise,
  pairwise |> rename(wave_from2 = wave_to, wave_to2 = wave_from) |> rename(wave_from = wave_from2, wave_to = wave_to2)
)
diagonal <- tibble(wave_from = wave_levels, wave_to = wave_levels,
                    spearman_incidence_corr = 1, top10pct_jaccard = 1, top20pct_jaccard = 1, contrib80pct_overlap_jaccard = 1)
pairwise_sym <- bind_rows(pairwise_sym, diagonal) |>
  mutate(wave_from = factor(wave_from, levels = wave_levels), wave_to = factor(wave_to, levels = wave_levels))

theme_v4 <- theme_minimal(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

plot_heatmap <- function(df, value_col, title, low = "#F7F7F7", high = "#762A83", limits = c(0, 1)) {
  ggplot(df, aes(wave_from, wave_to, fill = .data[[value_col]])) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", .data[[value_col]])), size = 2.4) +
    scale_fill_gradient(low = low, high = high, limits = limits, na.value = "grey95") +
    labs(title = title, x = NULL, y = NULL) +
    coord_fixed() + theme_v4
}

p_corr <- plot_heatmap(pairwise_sym, "spearman_incidence_corr", "Pairwise Spearman correlation of municipality incidence between waves", limits = c(-0.2, 1))
p_top10 <- plot_heatmap(pairwise_sym, "top10pct_jaccard", "Top-10% municipality Jaccard overlap between waves")
p_top20 <- plot_heatmap(pairwise_sym, "top20pct_jaccard", "Top-20% municipality Jaccard overlap between waves")
p_contrib80 <- plot_heatmap(pairwise_sym, "contrib80pct_overlap_jaccard", "80%-case-contributor municipality Jaccard overlap between waves")

library(patchwork)
ggsave(file.path(figure_dir, "pe_wave_pairwise_correlation_heatmap.png"), p_corr, width = 180, height = 160, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pe_wave_pairwise_correlation_heatmap.png"))

overlap_panel <- (p_top10 | p_top20) / p_contrib80
ggsave(file.path(figure_dir, "pe_top_municipality_overlap_heatmap.png"), overlap_panel, width = 260, height = 260, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pe_top_municipality_overlap_heatmap.png"))

# ---- Section 4C/4D: case-contributor turnover (adjacent waves) -----------
adjacent <- metrics |> dplyr::filter(comparison_type == "adjacent") |>
  mutate(wave_to = factor(wave_to, levels = wave_levels))

turnover_long <- adjacent |>
  select(wave_to, share_w2_cases_from_bottom50_w1, share_w2_cases_from_bottom25_w1, share_w2_cases_from_zero_w1) |>
  pivot_longer(-wave_to, names_to = "source", values_to = "share") |>
  mutate(source = recode(source,
                          share_w2_cases_from_bottom50_w1 = "From bottom-50% incidence munis (prior wave)",
                          share_w2_cases_from_bottom25_w1 = "From bottom-25% incidence munis (prior wave)",
                          share_w2_cases_from_zero_w1 = "From zero-case munis (prior wave)"))

p_turnover <- ggplot(turnover_long, aes(wave_to, share, fill = source)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  scale_fill_manual(values = c("From bottom-50% incidence munis (prior wave)" = "#56B4E9",
                                "From bottom-25% incidence munis (prior wave)" = "#E69F00",
                                "From zero-case munis (prior wave)" = "#D55E00")) +
  labs(title = "Pernambuco: share of each wave's cases from previously low/zero-incidence municipalities",
       subtitle = "A key spatial-turnover indicator (Section 4D) -- higher bars indicate more of each wave's burden arising from places that were quiet last time",
       x = "Wave (vs immediately preceding wave)", y = "Share of state cases") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(figure_dir, "pe_wave_case_contributor_turnover.png"), p_turnover, width = 220, height = 140, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pe_wave_case_contributor_turnover.png"))
