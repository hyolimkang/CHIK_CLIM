# Main Figures 1-3 (Section 8), journal-style theme, exported as PDF/SVG/PNG(>=600dpi).

required_packages <- c("here", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales", "svglite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork); library(scales) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
script_dir <- file.path(root, "02_Script/02_ceara_pipeline/07_counterfactuals")
source(file.path(script_dir, "01_config.R"))

res <- readRDS(file.path(DIR_RESULTS, "three_arm_results.rds"))
totals <- res$totals; vax <- res$vax_bookkeeping; dates <- res$dates
agg_Inew <- res$agg_Inew; agg_S <- res$agg_S; agg_N <- res$agg_N

COL_NOVACC <- "#4A4A4A"; COL_VACC <- "#0072B2"; COL_DIRECT <- "#D55E00"; COL_INDIRECT <- "#009E73"; COL_TOTAL <- "#000000"
AGE_PALETTE <- c("0-11" = "#5B8FA8", "12" = "#D55E00", "13-17" = "#E8A951", "18-64" = "#7A7A7A", "65+" = "#3D3D3D")

theme_journal <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", colour = "black"),
      plot.title = element_text(size = base_size, face = "plain", hjust = 0),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.title = element_text(size = base_size, colour = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_blank(), legend.key = element_blank(),
      legend.text = element_text(size = base_size - 1), legend.title = element_blank(),
      legend.position = "top", legend.justification = "left",
      strip.background = element_blank(), strip.text = element_text(size = base_size, face = "plain"),
      plot.margin = margin(4, 6, 4, 4)
    )
}
export_figure <- function(plot, stem, w = 180, h = 150) {
  ggsave(file.path(DIR_FIG_MAIN, paste0(stem, ".pdf")), plot, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(DIR_FIG_MAIN, paste0(stem, ".svg")), plot, width = w, height = h, units = "mm")
  ggsave(file.path(DIR_FIG_MAIN, paste0(stem, ".png")), plot, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  message("[saved] ", stem, " (.pdf/.svg/.png)")
}
date_scale <- scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2022-01-01"), "1 year"), date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01)))
summarise_weekly <- function(vec_by_draw_matrix) tibble(median = apply(vec_by_draw_matrix, 2, median), lo95 = apply(vec_by_draw_matrix, 2, quantile, .025), hi95 = apply(vec_by_draw_matrix, 2, quantile, .975))

# ================================================================
# FIGURE 1 -- historical vaccine counterfactual
# ================================================================
X_A_mat <- matrix(totals$X_A, ncol = length(dates), byrow = TRUE); X_B_mat <- matrix(totals$X_B, ncol = length(dates), byrow = TRUE)
wk_A <- summarise_weekly(X_A_mat); wk_A$week_start <- dates; wk_A$arm <- "No vaccine"
wk_B <- summarise_weekly(X_B_mat); wk_B$week_start <- dates; wk_B$arm <- "Routine age-12 (stress test)"
wk_df <- bind_rows(wk_A, wk_B)

p1A <- ggplot(wk_df, aes(week_start, median, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12 (stress test)" = COL_VACC)) +
  scale_fill_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12 (stress test)" = COL_VACC)) +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) + date_scale +
  labs(title = "Weekly latent infections", x = NULL, y = "Infections/week") + theme_journal()

cum_A_mat <- t(apply(X_A_mat, 1, cumsum)); cum_B_mat <- t(apply(X_B_mat, 1, cumsum))
cwk_A <- summarise_weekly(cum_A_mat); cwk_A$week_start <- dates; cwk_A$arm <- "No vaccine"
cwk_B <- summarise_weekly(cum_B_mat); cwk_B$week_start <- dates; cwk_B$arm <- "Routine age-12 (stress test)"
p1B <- ggplot(bind_rows(cwk_A, cwk_B), aes(week_start, median, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12 (stress test)" = COL_VACC)) +
  scale_fill_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12 (stress test)" = COL_VACC)) +
  scale_y_continuous(labels = label_number(scale = 1e-6, suffix = "M")) + date_scale +
  labs(title = "Cumulative infections", x = NULL, y = "Cumulative infections") + theme_journal()

age_averted <- bind_rows(lapply(REPORTING_AGE_LABELS, function(g) {
  averted_by_draw <- rowSums(agg_Inew[, , g, "A"]) - rowSums(agg_Inew[, , g, "B"])
  tibble(age_group = g, median = median(averted_by_draw), lo95 = quantile(averted_by_draw, .025), hi95 = quantile(averted_by_draw, .975))
})) |> mutate(age_group = factor(age_group, levels = REPORTING_AGE_LABELS))
p1C <- ggplot(age_averted, aes(age_group, median, fill = age_group)) +
  geom_col(width = 0.65) + geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.15, linewidth = 0.35, colour = "black") +
  scale_fill_manual(values = AGE_PALETTE, guide = "none") +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  labs(title = "Cumulative infections averted by age", x = NULL, y = "Infections averted") + theme_journal()

decomp <- read_csv(file.path(DIR_TABLES, "TABLE3_effect_decomposition.csv"), show_col_types = FALSE) |> dplyr::filter(quantity != "Indirect fraction of total effect")
decomp$quantity <- factor(c("Total effect", "Direct-only", "Transmission-mediated\nindirect"), levels = c("Total effect", "Direct-only", "Transmission-mediated\nindirect"))
p1D <- ggplot(decomp, aes(quantity, median, fill = quantity)) +
  geom_col(width = 0.6) + geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0.15, linewidth = 0.35, colour = "black") +
  scale_fill_manual(values = c("Total effect" = "black", "Direct-only" = COL_DIRECT, "Transmission-mediated\nindirect" = COL_INDIRECT), guide = "none") +
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  labs(title = "Effect decomposition", x = NULL, y = "Infections averted") + theme_journal()

fig1 <- (p1A | p1B) / (p1C | p1D) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 9, face = "bold"))
export_figure(fig1, "fig1_historical_age12_vaccine_impact", w = 180, h = 150)

# ================================================================
# FIGURE 2 -- transmission mechanism
# ================================================================
R0_mat <- matrix(totals$R0, ncol = length(dates), byrow = TRUE)
p2A <- ggplot(tibble(week_start = dates, R0 = apply(R0_mat, 2, median)), aes(week_start, R0)) +
  geom_line(colour = "black", linewidth = 0.5) + date_scale +
  labs(title = "Historical R0(t) -- identical across all three arms", x = NULL, y = expression(R[0](t))) + theme_journal()

Reff_A <- summarise_weekly(matrix(totals$Reff_A, ncol = length(dates), byrow = TRUE)); Reff_A$week_start <- dates; Reff_A$arm <- "No vaccine"
Reff_B <- summarise_weekly(matrix(totals$Reff_B, ncol = length(dates), byrow = TRUE)); Reff_B$week_start <- dates; Reff_B$arm <- "Routine age-12"
p2B <- ggplot(bind_rows(Reff_A, Reff_B), aes(week_start, median, colour = arm, fill = arm)) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12" = COL_VACC)) + scale_fill_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12" = COL_VACC)) +
  date_scale + labs(title = expression(R[eff](t)), x = NULL, y = expression(R[eff])) + theme_journal()

lambda_A <- summarise_weekly(matrix(totals$lambda_A, ncol = length(dates), byrow = TRUE)); lambda_A$week_start <- dates; lambda_A$arm <- "No vaccine"
lambda_B <- summarise_weekly(matrix(totals$lambda_B, ncol = length(dates), byrow = TRUE)); lambda_B$week_start <- dates; lambda_B$arm <- "Routine age-12"
p2C <- ggplot(bind_rows(lambda_A, lambda_B), aes(week_start, pmax(median, 1e-8), colour = arm)) +
  geom_line(linewidth = 0.5) + scale_y_log10() +
  scale_colour_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12" = COL_VACC)) +
  date_scale + labs(title = "Force of infection (log scale)", x = NULL, y = expression(lambda(t))) + theme_journal()

infness_A <- summarise_weekly(matrix(totals$infectiousness_A, ncol = length(dates), byrow = TRUE)); infness_A$week_start <- dates; infness_A$arm <- "No vaccine"
infness_B <- summarise_weekly(matrix(totals$infectiousness_B, ncol = length(dates), byrow = TRUE)); infness_B$week_start <- dates; infness_B$arm <- "Routine age-12"
p2D <- ggplot(bind_rows(infness_A, infness_B), aes(week_start, median, colour = arm)) +
  geom_line(linewidth = 0.5) + scale_colour_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12" = COL_VACC)) +
  date_scale + labs(title = "Infectiousness(t)", x = NULL, y = "Infectiousness") + theme_journal()

fig2 <- (p2A | p2B) / (p2C | p2D) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 9, face = "bold"))
export_figure(fig2, "fig2_transmission_feedback_mechanism", w = 180, h = 150)

# ================================================================
# FIGURE 3 -- age/cohort vaccination
# ================================================================
entrants_wk <- vax |> group_by(t) |> summarise(median = median(entrants_total_B), .groups = "drop") |> mutate(week_start = dates[t])
p3A <- ggplot(entrants_wk, aes(week_start, median)) + geom_col(fill = "grey60", width = 5) + date_scale +
  labs(title = "Weekly entrants to age 12", x = NULL, y = "Persons/week") + theme_journal()

doses_wk <- vax |> group_by(t) |> summarise(median = median(doses_B), .groups = "drop") |> mutate(week_start = dates[t])
p3B <- ggplot(doses_wk, aes(week_start, median)) + geom_col(fill = COL_VACC, width = 5) + date_scale +
  labs(title = "Doses administered", x = NULL, y = "Doses/week") + theme_journal()

Uvac_total_wk <- summarise_weekly(matrix(totals$Uvac_B, ncol = length(dates), byrow = TRUE)); Uvac_total_wk$week_start <- dates
p3C <- ggplot(Uvac_total_wk, aes(week_start, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = COL_VACC, alpha = 0.15) + geom_line(colour = COL_VACC, linewidth = 0.5) +
  date_scale + labs(title = "Vaccine-derived immune population (all ages)", x = NULL, y = "Persons") + theme_journal()

# Panel D: 3 illustrative birth cohorts' doses-received timing (median draw)
cohort_example <- vax |> dplyr::filter(doses_B > 1e-6) |> mutate(year = as.integer(format(week_start, "%Y")), birth_cohort = year - TARGET_AGE) |>
  group_by(draw, birth_cohort) |> summarise(week_vaccinated_start = min(week_start), week_vaccinated_end = max(week_start), total_doses = sum(doses_B), .groups = "drop") |>
  group_by(birth_cohort) |> summarise(week_vaccinated_start = median(week_vaccinated_start), week_vaccinated_end = median(week_vaccinated_end), total_doses = median(total_doses), .groups = "drop")
p3D <- ggplot(cohort_example, aes(y = factor(birth_cohort))) +
  geom_segment(aes(x = week_vaccinated_start, xend = week_vaccinated_end, yend = factor(birth_cohort)), colour = COL_VACC, linewidth = 3, lineend = "butt") +
  date_scale + labs(title = "Vaccination window per birth cohort (exactly once, at age 12)", x = NULL, y = "Birth cohort") + theme_journal()

fig3 <- (p3A | p3B) / (p3C | p3D) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 9, face = "bold"))
export_figure(fig3, "fig3_age12_cohort_vaccination", w = 180, h = 150)

message("\n[09_make_main_figures] All main figures exported (PDF/SVG/PNG>=600dpi) to ", DIR_FIG_MAIN)
