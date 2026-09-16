# Pernambuco spatial-turnover diagnostic -- cumulative municipality attack
# fraction under candidate q (user-requested addition, adapted from
# 32_bahia_global_q_local_consistency_audit/scripts/03_municipality_plausibility_audit.R
# and 04_map_implied_attack.R).
#
# CRUDE case-implied cumulative infection fraction per municipality,
# A_i = cumulative_reported_incidence_i / q, assuming constant q and
# ignoring transmission dynamics entirely. This is NOT a Stan estimate --
# explicitly labelled throughout. Evaluated at PE's fixed-q sweep grid
# (0.05-0.30) plus PE's estimated global-q posterior (median and 95% CrI,
# from Model B / U14 serology, outputs/modelB_U14_full).

required_packages <- c("dplyr", "readr", "tibble", "tidyr", "sf", "ggplot2", "rstan")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(tidyr); library(sf); library(ggplot2); library(rstan) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
pe_root <- file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/spatial_turnover_diagnostic")

pe_panel <- readRDS(file.path(table_dir, "pe_municipality_week_panel.rds"))

# ---- Pull PE's global-q posterior (Model B, U14 serology, full run) --------
b <- readRDS(file.path(pe_root, "outputs/modelB_U14_full/renewal_pe_global_q_fit_modelB_U14_full.rds"))
fit <- b$fit
q_draws <- plogis(as.vector(rstan::extract(fit, pars = "logit_q", permuted = FALSE)[, , 1]))
q_median <- median(q_draws); q_lo <- quantile(q_draws, .025); q_hi <- quantile(q_draws, .975)
message(sprintf("PE global-q (Model B, U14 serology): median=%.4f, 95%% CrI [%.4f, %.4f]", q_median, q_lo, q_hi))

# Fixed-q sweep grid (0.05-0.30) used elsewhere in this PE replication, plus
# the global-q posterior median and CrI bounds.
q_grid <- tibble(
  q = c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, q_median, q_lo, q_hi),
  q_label = c("fixed q=0.05", "fixed q=0.10", "fixed q=0.15", "fixed q=0.20", "fixed q=0.25", "fixed q=0.30",
              "global-q median", "global-q 2.5%", "global-q 97.5%")
)

checkpoints <- as.Date(c("2016-12-31", "2018-12-31", "2021-12-31", "2022-12-31", "2025-12-21"))

# Cumulative reported incidence R_i(t) per municipality at each checkpoint,
# using each municipality's population AT that checkpoint week.
cum_incidence <- bind_rows(lapply(checkpoints, function(cp) {
  up_to <- pe_panel |> filter(week_start <= cp)
  pop_at_cp <- pe_panel |> filter(week_start == max(week_start[week_start <= cp])) |>
    select(muni6, population_at_checkpoint = population)
  up_to |>
    group_by(muni6, name_muni) |>
    summarise(cumulative_cases = sum(cases_confirmed, na.rm = TRUE), .groups = "drop") |>
    left_join(pop_at_cp, by = "muni6") |>
    mutate(checkpoint = cp, cumulative_reported_incidence = cumulative_cases / population_at_checkpoint)
}))

# Crude case-implied cumulative infection fraction under each candidate q.
implied <- cum_incidence |>
  crossing(q_grid) |>
  mutate(A_i = cumulative_reported_incidence / q)

write_csv(implied, file.path(table_dir, "pe_municipality_implied_attack_by_q.csv"))
message("[saved] ", file.path(table_dir, "pe_municipality_implied_attack_by_q.csv"), " (", nrow(implied), " rows)")

# ---- Summary distribution + flagging, per checkpoint x q -------------------
summary_tbl <- implied |>
  group_by(checkpoint, q, q_label) |>
  summarise(
    n_munis = n(),
    median_A = median(A_i), pop_weighted_mean_A = weighted.mean(A_i, w = population_at_checkpoint),
    p90_A = quantile(A_i, .90), p95_A = quantile(A_i, .95), max_A = max(A_i),
    n_gt_025 = sum(A_i > 0.25), pct_gt_025 = 100 * mean(A_i > 0.25),
    n_gt_050 = sum(A_i > 0.50), pct_gt_050 = 100 * mean(A_i > 0.50),
    n_gt_075 = sum(A_i > 0.75), pct_gt_075 = 100 * mean(A_i > 0.75),
    n_gt_100 = sum(A_i > 1.00), pct_gt_100 = 100 * mean(A_i > 1.00),
    .groups = "drop"
  )

message("\n=== 2025 checkpoint: municipality-level crude implied attack fraction summary ===")
print(as.data.frame(summary_tbl |> filter(checkpoint == max(checkpoints)) |>
                       select(q_label, median_A, pop_weighted_mean_A, max_A, pct_gt_050, pct_gt_100)), digits = 3)

write_csv(summary_tbl, file.path(table_dir, "pe_municipality_implied_attack_summary.csv"))
message("\n[saved] ", file.path(table_dir, "pe_municipality_implied_attack_summary.csv"))

# ---- Flag A_i > 1 -----------------------------------------------------------
flagged <- implied |> filter(A_i > 1) |>
  left_join(pe_panel |> group_by(muni6) |> summarise(muni_max_pop = max(population), .groups = "drop"), by = "muni6")
message(sprintf("\nTotal municipality x checkpoint x q rows with A_i > 1: %d (out of %d)", nrow(flagged), nrow(implied)))
if (nrow(flagged) > 0) {
  flagged_at_median_q <- flagged |> filter(q_label == "global-q median") |> arrange(desc(A_i)) |>
    select(muni6, name_muni, checkpoint, cumulative_cases, population_at_checkpoint, cumulative_reported_incidence, A_i)
  message(sprintf("\nAt global-q median (q=%.4f), municipality x checkpoint rows with A_i > 1 (largest first):", q_median))
  print(as.data.frame(head(flagged_at_median_q, 20)), digits = 3)
  write_csv(flagged, file.path(table_dir, "pe_municipality_flagged_A_gt1.csv"))
  message("[saved] ", file.path(table_dir, "pe_municipality_flagged_A_gt1.csv"))
} else {
  message("No municipality x checkpoint x q combination exceeds A_i = 1.")
}

# ---- Maps: A_i(2025 | q) for the fixed-q grid AND global-q median ----------
poly <- readRDS(file.path(root, "01_Data/ibge_muni_polygons.rds")) |> filter(substr(muni6, 1, 2) == "26")
theme_map <- theme_void(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8, colour = "grey40"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.key.width = grid::unit(16, "pt"), legend.key.height = grid::unit(6, "pt"))

target_checkpoint <- as.Date("2025-12-21")
map_labels <- c("fixed q=0.05", "fixed q=0.10", "fixed q=0.20", "fixed q=0.30", "global-q median")

map_df <- implied |> filter(checkpoint == target_checkpoint, q_label %in% map_labels) |>
  mutate(q_label = factor(q_label, levels = map_labels))
map_data <- poly |> left_join(map_df |> select(muni6, q_label, A_i), by = "muni6", relationship = "many-to-many")

p_grid <- ggplot(map_data) +
  geom_sf(aes(fill = pmin(A_i, 1)), colour = "white", linewidth = 0.05) +
  facet_wrap(~q_label, ncol = 3) +
  scale_fill_viridis_c(option = "inferno", direction = -1, limits = c(0, 1), name = "A_i (capped at 1.0)", na.value = "grey90") +
  labs(title = "Pernambuco: crude case-implied cumulative attack fraction A_i(2025 | q)",
       subtitle = "CRUDE, assumes constant q and ignores transmission dynamics -- NOT a Stan estimate") +
  theme_map

ggsave(file.path(figure_dir, "pe_municipality_implied_attack_by_q_maps.png"), p_grid, width = 260, height = 200, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pe_municipality_implied_attack_by_q_maps.png"))

# ---- A_i > 1 flag map across the same q values -----------------------------
flag_df <- implied |> filter(checkpoint == target_checkpoint, q_label %in% map_labels) |>
  mutate(flag = A_i > 1, q_label = factor(q_label, levels = map_labels))
map_flag <- poly |> left_join(flag_df |> select(muni6, q_label, flag), by = "muni6", relationship = "many-to-many")

p_flag <- ggplot(map_flag) +
  geom_sf(aes(fill = flag), colour = "white", linewidth = 0.05) +
  scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey92"), na.value = "grey92",
                     labels = c(`TRUE` = "A_i > 1 (flagged)", `FALSE` = "A_i <= 1"), name = NULL) +
  facet_wrap(~q_label, ncol = 3) +
  labs(title = "Pernambuco: municipalities where crude case-implied attack fraction exceeds 1.0 (2025)",
       subtitle = "A falsification diagnostic under the lifelong-first-infection interpretation -- flagged, not auto-excluded") +
  theme_map

ggsave(file.path(figure_dir, "pe_municipality_A_gt1_flag_map.png"), p_flag, width = 260, height = 200, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pe_municipality_A_gt1_flag_map.png"))
