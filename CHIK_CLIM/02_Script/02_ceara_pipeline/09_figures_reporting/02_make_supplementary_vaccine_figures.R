# Supplementary Figures S1-S8 (Section 9), same journal-style theme.

required_packages <- c("here", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork); library(scales) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
script_dir <- file.path(root, "02_Script/02_ceara_pipeline/07_counterfactuals")
source(file.path(script_dir, "01_config.R"))

res <- readRDS(file.path(DIR_RESULTS, "three_arm_results.rds"))
totals <- res$totals; vax <- res$vax_bookkeeping; dates <- res$dates
agg_Inew <- res$agg_Inew; agg_S <- res$agg_S; agg_N <- res$agg_N
N <- length(dates)

COL_NOVACC <- "#4A4A4A"; COL_VACC <- "#0072B2"; COL_DIRECT <- "#D55E00"
AGE_PALETTE <- c("0-11" = "#5B8FA8", "12" = "#D55E00", "13-17" = "#E8A951", "18-64" = "#7A7A7A", "65+" = "#3D3D3D")
theme_journal <- function(base_size = 8) {
  theme_classic(base_size = base_size) + theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), axis.ticks = element_line(colour = "black", linewidth = 0.3),
    panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "top", legend.justification = "left", legend.title = element_blank(), strip.background = element_blank())
}
date_scale <- scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2022-01-01"), "1 year"), date_labels = "%Y")
export_supp <- function(plot, stem, w = 180, h = 130) { ggsave(file.path(DIR_FIG_SUPP, paste0(stem, ".png")), plot, width = w, height = h, units = "mm", dpi = 600, bg = "white"); message("[saved] ", stem, ".png") }

# ---- S1: no-vaccine replay QA ----
climate <- readRDS(ACCEPTED_CLIMATE_FIT)
draws_v49 <- rstan::extract(climate$fit, pars = c("S", "U", "X"), permuted = TRUE)
rep_draw <- res$draw_idx[1]
rep_sub <- totals |> dplyr::filter(draw == rep_draw)
sS1 <- (p1 <- ggplot(rep_sub, aes(week_start)) + geom_line(aes(y = X_A, colour = "Age-model closed-loop"), linewidth = 0.4, linetype = "dashed") +
  geom_line(aes(y = draws_v49$X[rep_draw, seq_len(N)], colour = "Accepted climate-v4.9"), linewidth = 0.7, alpha = 0.6) +
  scale_colour_manual(values = c("Accepted climate-v4.9" = COL_NOVACC, "Age-model closed-loop" = COL_VACC)) +
  date_scale + labs(title = "S1: No-vaccine replay QA -- X(t)", x = NULL, y = "Infections/week") + theme_journal())
export_supp(sS1, "S1_no_vaccine_replay_QA")

# ---- S2: population/compartment conservation ----
conservation_A <- totals |> mutate(err_check = S_A + U_A - N_A, arm = "Arm A (S+Uinf-N)")
conservation_B <- totals |> mutate(err_check = S_B + U_B + Uvac_B - N_B, arm = "Arm B (S+Uinf+Uvac-N)")
conservation_C <- totals |> mutate(err_check = S_C + U_C + Uvac_C - N_C, arm = "Arm C (S+Uinf+Uvac-N)")
conservation_all <- bind_rows(conservation_A, conservation_B, conservation_C)
sS2 <- ggplot(conservation_all, aes(week_start, err_check, group = draw)) + geom_line(alpha = 0.1, colour = COL_VACC) +
  facet_wrap(~arm, ncol = 1) + date_scale + labs(title = "S2: population/compartment conservation, all arms", x = NULL, y = "Residual (should be ~0)") + theme_journal()
export_supp(sS2, "S2_compartment_conservation", h = 180)

# ---- S3: Arm C direct-only diagnostic -- non-target groups overlap Arm A ----
s3_data <- bind_rows(lapply(c("0-11", "65+"), function(g) {
  tibble(week_start = dates, age_group = g, A = apply(agg_Inew[, , g, "A"], 2, median), C = apply(agg_Inew[, , g, "C"], 2, median))
}))
sS3 <- ggplot(s3_data, aes(week_start)) + geom_line(aes(y = A, colour = "Arm A (no vaccine)"), linewidth = 0.6) +
  geom_line(aes(y = C, colour = "Arm C (direct-only)"), linewidth = 0.4, linetype = "dashed") +
  facet_wrap(~age_group, ncol = 1, scales = "free_y") + scale_colour_manual(values = c("Arm A (no vaccine)" = COL_NOVACC, "Arm C (direct-only)" = COL_DIRECT)) +
  date_scale + labs(title = "S3: Arm C non-target groups overlap Arm A (no legitimate pathway to differ)", x = NULL, y = "Infections/week") + theme_journal()
export_supp(sS3, "S3_arm_C_direct_only_diagnostic", h = 160)

# ---- S4: age-specific weekly infection curves (all 5 groups, Arm A) ----
s4_data <- bind_rows(lapply(REPORTING_AGE_LABELS, function(g) tibble(week_start = dates, age_group = g, median = apply(agg_Inew[, , g, "A"], 2, median))))
s4_data$age_group <- factor(s4_data$age_group, levels = REPORTING_AGE_LABELS)
sS4 <- ggplot(s4_data, aes(week_start, median, colour = age_group)) + geom_line(linewidth = 0.5) +
  scale_colour_manual(values = AGE_PALETTE) + date_scale + labs(title = "S4: age-specific weekly infections, Arm A (no vaccine)", x = NULL, y = "Infections/week") + theme_journal()
export_supp(sS4, "S4_age_specific_weekly_curves")

# ---- S5: posterior distribution of total infections averted ----
decomp_by_draw <- read_csv(file.path(DIR_TABLES, "TABLE3_decomposition_by_draw.csv"), show_col_types = FALSE)
sS5 <- ggplot(decomp_by_draw, aes(total_effect)) + geom_histogram(bins = 30, fill = COL_VACC, colour = "white", linewidth = 0.2) +
  labs(title = "S5: posterior distribution, total infections averted", x = "Infections averted", y = "Posterior draws") + theme_journal()
export_supp(sS5, "S5_posterior_total_averted", h = 100)

# ---- S6: posterior distribution of indirect-effect fraction ----
sS6 <- ggplot(decomp_by_draw, aes(indirect_fraction)) + geom_histogram(bins = 30, fill = "#009E73", colour = "white", linewidth = 0.2) +
  scale_x_continuous(labels = label_percent()) + labs(title = "S6: posterior distribution, indirect-effect fraction", x = "Indirect fraction of total effect", y = "Posterior draws") + theme_journal()
export_supp(sS6, "S6_posterior_indirect_fraction", h = 100)

# ---- S7: climate_multiplier(t) alongside epidemic trajectory ----
climate_mult <- rstan::extract(climate$fit, pars = "climate_multiplier")$climate_multiplier
mult_med <- apply(climate_mult[, seq_len(N)], 2, median)
s7_top <- ggplot(tibble(week_start = dates, mult = mult_med), aes(week_start, mult)) + geom_line(colour = "#CC79A7", linewidth = 0.5) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey60") + date_scale + labs(title = "S7: climate_multiplier(t)", x = NULL, y = "Multiplier") + theme_journal()
s7_bottom <- ggplot(rep_sub, aes(week_start, X_A)) + geom_line(colour = COL_NOVACC, linewidth = 0.5) + date_scale + labs(x = NULL, y = "Infections/week (Arm A)") + theme_journal()
sS7 <- s7_top / s7_bottom
export_supp(sS7, "S7_climate_multiplier_vs_epidemic", h = 140)

# ---- S8: Reff threshold crossing, no-vaccine vs vaccine ----
reff_cross <- totals |> group_by(t) |> summarise(Reff_A = median(Reff_A), Reff_B = median(Reff_B), .groups = "drop") |> mutate(week_start = dates[t])
n_weeks_above1_A <- sum(reff_cross$Reff_A > 1); n_weeks_above1_B <- sum(reff_cross$Reff_B > 1)
sS8 <- ggplot(reff_cross, aes(week_start)) + geom_hline(yintercept = 1, linetype = "22", colour = "grey60") +
  geom_line(aes(y = Reff_A, colour = "No vaccine"), linewidth = 0.5) + geom_line(aes(y = Reff_B, colour = "Routine age-12"), linewidth = 0.5) +
  scale_colour_manual(values = c("No vaccine" = COL_NOVACC, "Routine age-12" = COL_VACC)) + date_scale +
  labs(title = sprintf("S8: R_eff threshold crossing (weeks Reff>1: A=%d, B=%d)", n_weeks_above1_A, n_weeks_above1_B), x = NULL, y = expression(R[eff])) + theme_journal()
export_supp(sS8, "S8_Reff_threshold_crossing")

message("\n[10_make_supplementary_figures] All supplementary figures saved to ", DIR_FIG_SUPP)
