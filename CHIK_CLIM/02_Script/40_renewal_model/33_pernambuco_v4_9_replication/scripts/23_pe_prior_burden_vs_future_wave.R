# Pernambuco spatial-turnover diagnostic -- Section 5.
#
# Descriptive test: do municipalities with greater cumulative PRIOR burden
# (at wave onset) contribute less to the CURRENT wave? Quintiles of
# cumulative_previous_incidence_iw vs subsequent-wave incidence_iw. This is
# NOT a causal susceptible-depletion estimate (reporting heterogeneity,
# climate, importation etc. may confound it) -- purely descriptive.

required_packages <- c("dplyr", "readr", "tibble", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(ggplot2) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/spatial_turnover_diagnostic")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

burden <- read_csv(file.path(table_dir, "pe_municipality_wave_burden.csv"), show_col_types = FALSE)
waves <- read_csv(file.path(table_dir, "pe_wave_definitions.csv"), show_col_types = FALSE) |> arrange(start_week)

# Exclude the first wave (PE_wave_03) since it has no meaningful prior
# cumulative incidence within the 2015-2025 window (near-zero for all munis).
analysis_waves <- waves$wave_id[-1]

quintile_df <- burden |>
  filter(wave_id %in% analysis_waves) |>
  filter(!is.na(cumulative_previous_incidence_iw)) |>
  group_by(wave_id) |>
  mutate(prior_quintile = ntile(cumulative_previous_incidence_iw, 5)) |>
  ungroup()

quintile_summary <- quintile_df |>
  group_by(wave_id, prior_quintile) |>
  summarise(
    n_munis = n(),
    median_prior_incidence = median(cumulative_previous_incidence_iw),
    median_subsequent_incidence = median(incidence_iw),
    q25_subsequent = quantile(incidence_iw, .25),
    q75_subsequent = quantile(incidence_iw, .75),
    pop_weighted_mean_subsequent = weighted.mean(incidence_iw, w = population_mean_iw),
    share_of_state_cases = sum(reported_cases_iw) / sum(quintile_df$reported_cases_iw[quintile_df$wave_id == wave_id[1]]),
    .groups = "drop"
  )

message("=== Quintile summary (by wave) ===")
print(as.data.frame(quintile_summary), digits = 3)

write_csv(quintile_summary, file.path(table_dir, "pe_prior_burden_future_wave.csv"))
message("\n[saved] ", file.path(table_dir, "pe_prior_burden_future_wave.csv"))

# Overall (pooled across waves) Spearman correlation
overall_rho <- cor(quintile_df$cumulative_previous_incidence_iw, quintile_df$incidence_iw, method = "spearman")
message(sprintf("\nPooled (all waves combined) Spearman correlation, prior cumulative incidence vs subsequent-wave incidence: %.3f", overall_rho))

per_wave_rho <- quintile_df |>
  group_by(wave_id) |>
  summarise(rho = cor(cumulative_previous_incidence_iw, incidence_iw, method = "spearman"), .groups = "drop")
message("\nPer-wave Spearman correlation:")
print(as.data.frame(per_wave_rho))
write_csv(per_wave_rho, file.path(table_dir, "pe_prior_burden_per_wave_spearman.csv"))

# ---- Figure: scatter + binned (quintile) summary ---------------------------
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

p_scatter <- ggplot(quintile_df, aes(cumulative_previous_incidence_iw + 1, incidence_iw + 1)) +
  geom_point(alpha = 0.15, size = 0.6, colour = "#315A7D") +
  geom_smooth(method = "loess", colour = "#D55E00", se = TRUE, linewidth = 0.7) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~wave_id, scales = "free") +
  labs(title = "Pernambuco: cumulative prior incidence vs subsequent-wave incidence (municipality level)",
       subtitle = sprintf("Pooled Spearman rho = %.3f. Descriptive only -- not a causal depletion estimate.", overall_rho),
       x = "Cumulative incidence before wave onset (+1, log scale)", y = "Wave incidence per 100k (+1, log scale)") +
  theme_v4 + theme(strip.background = element_blank())

ggsave(file.path(figure_dir, "pe_prior_burden_vs_future_incidence.png"), p_scatter, width = 260, height = 200, units = "mm", dpi = 300)
message("[saved] ", file.path(figure_dir, "pe_prior_burden_vs_future_incidence.png"))

p_quintile <- ggplot(quintile_summary, aes(factor(prior_quintile), median_subsequent_incidence)) +
  geom_col(fill = "#76558F", alpha = 0.8) +
  geom_errorbar(aes(ymin = q25_subsequent, ymax = q75_subsequent), width = 0.2, colour = "grey30") +
  facet_wrap(~wave_id, scales = "free_y") +
  labs(title = "Median subsequent-wave incidence by prior-cumulative-incidence quintile",
       subtitle = "Quintile 1 = lowest prior burden, Quintile 5 = highest prior burden. Error bars = IQR.",
       x = "Prior cumulative incidence quintile", y = "Median wave incidence per 100k") +
  theme_v4 + theme(strip.background = element_blank())

ggsave(file.path(figure_dir, "pe_prior_burden_quintile_summary.png"), p_quintile, width = 260, height = 200, units = "mm", dpi = 300)
message("[saved] ", file.path(figure_dir, "pe_prior_burden_quintile_summary.png"))
