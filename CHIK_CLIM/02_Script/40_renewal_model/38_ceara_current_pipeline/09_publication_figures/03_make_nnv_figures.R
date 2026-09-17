# Publication-ready NNV figures (NNV-1 .. NNV-5), same journal style as the
# main/supplementary figures already produced for this analysis.
#
# REVISION (this version): terminology and visualisation only -- reuses the
# already-saved posterior draw-level results from 11_calculate_nnv.R.
# No transmission/vaccine simulation is rerun or altered.
#   - fig_nnv2: continuous line plot -> discrete YEAR-END CUMULATIVE NNV bar
#     chart (draw-level ratio computed first, then summarised).
#   - fig_nnv4: violin plots -> point + 50%/95% CrI interval estimates for
#     the main figure; the violin version is kept ONLY as a supplementary
#     diagnostic (figures/supplementary/fig_s_nnv_posterior_distributions_violin).
#   - fig_nnv5: age-specific dose/effect RATIO ("efficiency") removed from
#     all main outputs; replaced by an age-specific BENEFIT figure
#     (infections averted + % reduction only). No age group is labelled
#     with its own NNV anywhere in this script.
#   - All axis/label text uses scales::label_comma()/label_number() --
#     no scientific notation (1e+05 style) anywhere.

required_packages <- c("here", "dplyr", "readr", "tibble", "tidyr", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(tibble); library(tidyr); library(ggplot2); library(patchwork); library(scales) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
script_dir <- file.path(root, "02_Script/40_renewal_model/38_ceara_current_pipeline/06_vaccine_counterfactual")
source(file.path(script_dir, "01_config.R"))

DIR_RESULTS_NNV <- file.path(DIR_RESULTS, "nnv") # DIR_FIG_NNV, DIR_FIG_SUPP come from 00_config.R (03_Output/figures/...)

nnv_draws <- readRDS(file.path(DIR_RESULTS_NNV, "nnv_posterior_draws.rds"))
nnv_by_horizon <- readRDS(file.path(DIR_RESULTS_NNV, "nnv_by_horizon.rds"))
nnv_year_end_by_draw <- readRDS(file.path(DIR_RESULTS_NNV, "nnv_year_end_cumulative_by_draw.rds"))
age_effect_by_draw <- readRDS(file.path(DIR_RESULTS_NNV, "nnv_age_specific_by_draw.rds"))
nnv_susceptibility <- readRDS(file.path(DIR_RESULTS_NNV, "nnv_susceptibility_link.rds"))
dates <- readRDS(file.path(DIR_RESULTS, "three_arm_results.rds"))$dates

COL_INF <- "#0072B2"; COL_REP <- "#D55E00"; COL_DOSES <- "#4A4A4A"; COL_TOTAL <- "#0072B2"; COL_DIRECT <- "#D55E00"; COL_INDIRECT <- "#009E73"; COL_TARGET <- "#D55E00"; COL_NONTARGET <- "#0072B2"
theme_journal <- function(base_size = 8) {
  theme_classic(base_size = base_size) + theme(
    axis.line = element_line(colour = "black", linewidth = 0.3), axis.ticks = element_line(colour = "black", linewidth = 0.3),
    panel.grid = element_blank(), plot.background = element_rect(fill = "white", colour = NA),
    legend.position = "top", legend.justification = "left", legend.title = element_blank(), strip.background = element_blank())
}
export_figure <- function(plot, stem, dir = DIR_FIG_NNV, w = 180, h = 150) {
  ggsave(file.path(dir, paste0(stem, ".pdf")), plot, width = w, height = h, units = "mm", device = grDevices::cairo_pdf, bg = "white")
  ggsave(file.path(dir, paste0(stem, ".svg")), plot, width = w, height = h, units = "mm", bg = "white")
  ggsave(file.path(dir, paste0(stem, ".png")), plot, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  message("[saved] ", stem, ".{pdf,svg,png}")
}
qCrI <- function(x, probs, na.rm = TRUE) quantile(x, probs, na.rm = na.rm, names = FALSE)
comma_num <- scales::label_comma(accuracy = 1) # human-readable integers everywhere -- never 1e+05 style
comma_dec <- scales::label_comma(accuracy = 0.01)

# Remove any old figure files being superseded/renamed by this revision, so
# stale outputs don't linger under the new naming scheme.
old_files <- c("fig_nnv2_nnv_over_time", "fig_nnv4_direct_indirect_nnv", "fig_nnv5_age_specific_efficiency")
invisible(lapply(old_files, function(stem) file.remove(Sys.glob(file.path(DIR_FIG_NNV, paste0(stem, ".*"))))))

# ================================================================
# FIGURE NNV-1: primary NNV, forest-style across evaluation horizons
# (unchanged design -- only number formatting fixed, no scientific notation)
# ================================================================
horizon_labels <- c("1yr" = "1 year", "2yr" = "2 years", "3yr" = "3 years", "5yr" = "5 years", "end_of_analysis_2021" = "Full window\n(2015-2021)")
forest_data <- nnv_by_horizon |> group_by(horizon) |>
  summarise(inf_med = median(NNV_infection_h, na.rm = TRUE), inf_lo50 = qCrI(NNV_infection_h, .25), inf_hi50 = qCrI(NNV_infection_h, .75),
            inf_lo95 = qCrI(NNV_infection_h, .025), inf_hi95 = qCrI(NNV_infection_h, .975),
            rep_med = median(NNV_reported_case_h, na.rm = TRUE), rep_lo50 = qCrI(NNV_reported_case_h, .25), rep_hi50 = qCrI(NNV_reported_case_h, .75),
            rep_lo95 = qCrI(NNV_reported_case_h, .025), rep_hi95 = qCrI(NNV_reported_case_h, .975), .groups = "drop") |>
  mutate(horizon_lab = factor(horizon_labels[as.character(horizon)], levels = rev(horizon_labels)))

panelA <- ggplot(forest_data, aes(x = inf_med, y = horizon_lab)) +
  geom_linerange(aes(xmin = inf_lo95, xmax = inf_hi95), linewidth = 0.4, colour = COL_INF) +
  geom_linerange(aes(xmin = inf_lo50, xmax = inf_hi50), linewidth = 1.1, colour = COL_INF) +
  geom_point(size = 1.6, colour = COL_INF) + scale_x_log10(labels = comma_num) +
  labs(title = "A", x = "Number needed to vaccinate\nper infection averted (log scale)", y = NULL) + theme_journal()
panelB <- ggplot(forest_data, aes(x = rep_med, y = horizon_lab)) +
  geom_linerange(aes(xmin = rep_lo95, xmax = rep_hi95), linewidth = 0.4, colour = COL_REP) +
  geom_linerange(aes(xmin = rep_lo50, xmax = rep_hi50), linewidth = 1.1, colour = COL_REP) +
  geom_point(size = 1.6, colour = COL_REP) + scale_x_log10(labels = comma_num) +
  labs(title = "B", x = "Number needed to vaccinate\nper reported case averted (log scale)", y = NULL) + theme_journal()
fig_nnv1 <- panelA / panelB
export_figure(fig_nnv1, "fig_nnv1_primary_nnv", h = 115)

# ================================================================
# FIGURE NNV-2 (REVISED): four panels giving both the ANNUAL (year-only)
# picture and the CUMULATIVE picture side by side, so the reader can see
# WHY cumulative NNV moves the way it does:
#   A. Annual infections, no-vaccine vs vaccine (counts, not a ratio)
#   B. Annual infections averted, I^NV_y - I^V_y, with a zero baseline --
#      makes negative (rebound) years immediately visible
#   C. Cumulative infections averted (entire population)
#   D. Cumulative NNV per infection averted (draw-level ratio computed
#      first in 11_calculate_nnv.R, then summarised; 2015 flagged as
#      unstable rather than plotted as a falsely precise bar -- Section 7
#      Option B, linear axis retained since 2016-2021 span <4x)
# Panels A-B use annual_by_draw (year-only sums); panels C-D use
# nnv_year_end_by_draw (cumulative-since-programme-start through each
# year-end) -- these are deliberately different data sources, not the same
# numbers plotted twice.
# ================================================================
nnv_annual_by_draw <- readRDS(file.path(DIR_RESULTS_NNV, "nnv_annual_by_draw.rds"))
annual_summary <- nnv_annual_by_draw |> group_by(year) |>
  summarise(novacc_med = median(infections_no_vaccine_annual), novacc_lo = qCrI(infections_no_vaccine_annual, .025), novacc_hi = qCrI(infections_no_vaccine_annual, .975),
            vacc_med = median(infections_vaccine_annual), vacc_lo = qCrI(infections_vaccine_annual, .025), vacc_hi = qCrI(infections_vaccine_annual, .975),
            avert_med = median(infections_averted_annual), avert_lo = qCrI(infections_averted_annual, .025), avert_hi = qCrI(infections_averted_annual, .975), .groups = "drop") |>
  mutate(year = factor(year, levels = sort(unique(year))), avert_sign = ifelse(avert_med >= 0, "Positive (benefit)", "Negative (rebound)"))

year_summary <- nnv_year_end_by_draw |> group_by(year) |>
  summarise(inf_avert_med = median(cumulative_infections_averted_y), inf_avert_lo = qCrI(cumulative_infections_averted_y, .025), inf_avert_hi = qCrI(cumulative_infections_averted_y, .975),
            nnv_inf_med = median(NNV_infection_y, na.rm = TRUE), nnv_inf_lo = qCrI(NNV_infection_y, .025), nnv_inf_hi = qCrI(NNV_infection_y, .975),
            .groups = "drop") |>
  mutate(year = factor(year, levels = sort(unique(year))), stable = inf_avert_med >= 1000) # matches 11_calculate_nnv.R's stability rule
year_levels_chr <- levels(year_summary$year) # pin chronological order explicitly on every panel below -- a factor's
# own level order is not always enough once a layer is built from a row-filtered subset (see prior fix in this file).

panelA2 <- ggplot(annual_summary, aes(x = year)) +
  geom_col(aes(y = novacc_med, fill = "No vaccine"), position = position_nudge(x = -0.2), width = 0.35) +
  geom_errorbar(aes(ymin = novacc_lo, ymax = novacc_hi), position = position_nudge(x = -0.2), width = 0.12, linewidth = 0.35) +
  geom_col(aes(y = vacc_med, fill = "Vaccine (age-12 routine)"), position = position_nudge(x = 0.2), width = 0.35) +
  geom_errorbar(aes(ymin = vacc_lo, ymax = vacc_hi), position = position_nudge(x = 0.2), width = 0.12, linewidth = 0.35) +
  scale_fill_manual(values = c("No vaccine" = COL_DOSES, "Vaccine (age-12 routine)" = COL_INF), name = NULL) +
  scale_x_discrete(limits = year_levels_chr) + scale_y_continuous(labels = comma_num) +
  labs(title = "A", x = NULL, y = "Annual infections") + theme_journal()

panelB2 <- ggplot(annual_summary, aes(x = year, y = avert_med, fill = avert_sign)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.35) +
  geom_col(width = 0.6) + geom_errorbar(aes(ymin = avert_lo, ymax = avert_hi), width = 0.15, linewidth = 0.4) +
  scale_fill_manual(values = c("Positive (benefit)" = COL_INF, "Negative (rebound)" = COL_REP), name = NULL) +
  scale_x_discrete(limits = year_levels_chr) + scale_y_continuous(labels = comma_num) +
  labs(title = "B", x = NULL, y = expression(paste("Annual infections averted, ", I^{NV}[y] - I^{V}[y]))) + theme_journal()

panelC2 <- ggplot(year_summary, aes(year, inf_avert_med)) + geom_col(fill = "#CC79A7", width = 0.6) +
  geom_errorbar(aes(ymin = inf_avert_lo, ymax = inf_avert_hi), width = 0.15, linewidth = 0.4) +
  scale_x_discrete(limits = year_levels_chr) + scale_y_continuous(labels = comma_num) +
  labs(title = "C", x = NULL, y = "Cumulative infections\naverted (entire population)") + theme_journal()

y_max_inf <- max(year_summary$nnv_inf_hi[year_summary$stable])
panelD2 <- ggplot(year_summary, aes(x = year)) +
  geom_col(data = function(d) filter(d, stable), aes(y = nnv_inf_med), fill = COL_INF, width = 0.6) +
  geom_errorbar(data = function(d) filter(d, stable), aes(ymin = nnv_inf_lo, ymax = nnv_inf_hi), width = 0.15, linewidth = 0.4) +
  geom_text(data = function(d) filter(d, !stable), aes(y = 0.06 * y_max_inf), label = "Insufficient accrued\nbenefit for stable NNV", size = 2.1, angle = 90, fontface = "italic", colour = "grey40", hjust = 0) +
  scale_x_discrete(limits = year_levels_chr) + scale_y_continuous(labels = comma_num) +
  labs(title = "D", x = "Calendar year", y = "Cumulative NNV\nper infection averted") + theme_journal()

fig_nnv2 <- (panelA2 | panelB2) / (panelC2 | panelD2)
export_figure(fig_nnv2, "fig_nnv2_year_end_cumulative_nnv", h = 170)

# ================================================================
# FIGURE NNV-4 (REVISED): effect decomposition, MAIN = point + interval
# (median, 50% CrI thick, 95% CrI thin) -- replaces violin plots in the
# main figure. Linear axes throughout (values span <10x, a log axis is not
# needed for readability here per the "avoid log unless required" rule).
# ================================================================
decomp_long_counts <- nnv_draws |> select(draw, `Total\n(A-B)` = infections_averted, `Direct-only\n(A-C)` = direct_only_effect, `Indirect\n(C-B)` = indirect_effect) |>
  pivot_longer(-draw, names_to = "component", values_to = "value") |> mutate(component = factor(component, levels = c("Total\n(A-B)", "Direct-only\n(A-C)", "Indirect\n(C-B)")))
decomp_long_ratio <- nnv_draws |> select(draw, `Total\n(A-B)` = NNV_total_infection, `Direct-only\n(A-C)` = NNV_direct_only_infection, `Indirect\n(C-B)` = NNV_indirect_infection) |>
  pivot_longer(-draw, names_to = "component", values_to = "value") |> mutate(component = factor(component, levels = c("Total\n(A-B)", "Direct-only\n(A-C)", "Indirect\n(C-B)")))
comp_pal <- c("Total\n(A-B)" = COL_TOTAL, "Direct-only\n(A-C)" = COL_DIRECT, "Indirect\n(C-B)" = COL_INDIRECT)

interval_summary <- function(df) df |> group_by(component) |> summarise(med = median(value, na.rm = TRUE), lo50 = qCrI(value, .25), hi50 = qCrI(value, .75), lo95 = qCrI(value, .025), hi95 = qCrI(value, .975), .groups = "drop")
counts_summary <- interval_summary(decomp_long_counts); ratio_summary <- interval_summary(decomp_long_ratio)

panelA4 <- ggplot(counts_summary, aes(x = med, y = component, colour = component)) +
  geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.4) + geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 1.3) +
  geom_point(size = 2) + scale_colour_manual(values = comp_pal) + scale_x_continuous(labels = comma_num) +
  labs(title = "A", x = "Infections averted", y = NULL) + theme_journal() + theme(legend.position = "none")
panelB4 <- ggplot(ratio_summary, aes(x = med, y = component, colour = component)) +
  geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.4) + geom_linerange(aes(xmin = lo50, xmax = hi50), linewidth = 1.3) +
  geom_point(size = 2) + scale_colour_manual(values = comp_pal) + scale_x_continuous(labels = comma_num) +
  labs(title = "B", x = "Doses / effect (diagnostic ratio -- NOT age-specific NNV)", y = NULL) + theme_journal() + theme(legend.position = "none")
fig_nnv4 <- panelA4 / panelB4
export_figure(fig_nnv4, "fig_nnv4_effect_decomposition_interval", h = 110)

# ---- Supplementary: violin version retained as a distributional diagnostic ----
panelA4v <- ggplot(decomp_long_counts, aes(component, value, fill = component)) + geom_violin(colour = NA, alpha = 0.8) +
  stat_summary(fun = median, geom = "point", size = 1.6, colour = "black") + scale_fill_manual(values = comp_pal) + scale_y_continuous(labels = comma_num) +
  labs(title = "A", x = NULL, y = "Infections averted") + theme_journal() + theme(legend.position = "none")
panelB4v <- ggplot(decomp_long_ratio, aes(component, value, fill = component)) + geom_violin(colour = NA, alpha = 0.8) +
  stat_summary(fun = median, geom = "point", size = 1.6, colour = "black") + scale_fill_manual(values = comp_pal) + scale_y_continuous(labels = comma_num) +
  labs(title = "B", x = NULL, y = "Doses / effect (diagnostic ratio)") + theme_journal() + theme(legend.position = "none")
fig_s_violin <- panelA4v / panelB4v
export_figure(fig_s_violin, "fig_s_nnv_posterior_distributions_violin", dir = DIR_FIG_SUPP, h = 150)

# ================================================================
# FIGURE NNV-5 (REVISED): age-specific BENEFIT ONLY -- infections averted
# and percent reduction. No age-group "NNV" anywhere (doses were
# administered only to the age-12 entrant flow; other groups' benefit is
# transmission-mediated, not their own dose-denominated quantity).
# ================================================================
age_summary <- age_effect_by_draw |> mutate(age_group = factor(age_group, levels = REPORTING_AGE_LABELS)) |> group_by(age_group) |>
  summarise(inf_avert_med = median(infections_averted_g), inf_avert_lo = qCrI(infections_averted_g, .025), inf_avert_hi = qCrI(infections_averted_g, .975),
            pct_red_med = median(percent_reduction_g), pct_red_lo = qCrI(percent_reduction_g, .025), pct_red_hi = qCrI(percent_reduction_g, .975), .groups = "drop") |>
  mutate(is_target = age_group == "12")

panelA5 <- ggplot(age_summary, aes(age_group, inf_avert_med, fill = is_target)) + geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = inf_avert_lo, ymax = inf_avert_hi), width = 0.15, linewidth = 0.4) +
  scale_fill_manual(values = c(`TRUE` = COL_TARGET, `FALSE` = COL_NONTARGET), labels = c(`TRUE` = "Vaccine-targeted (age 12)", `FALSE` = "Not vaccine-targeted"), name = NULL) +
  scale_y_continuous(labels = comma_num) + labs(title = "A", x = "Age group", y = "Infections averted") + theme_journal()
panelB5 <- ggplot(age_summary, aes(age_group, pct_red_med, fill = is_target)) + geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = pct_red_lo, ymax = pct_red_hi), width = 0.15, linewidth = 0.4) +
  scale_fill_manual(values = c(`TRUE` = COL_TARGET, `FALSE` = COL_NONTARGET)) +
  scale_y_continuous(labels = scales::label_number(suffix = "%", accuracy = 1)) + labs(title = "B", x = "Age group", y = "Percentage reduction\nin infections") + theme_journal()
fig_nnv5 <- panelA5 / panelB5 +
  plot_annotation(caption = "Age 12 = directly vaccine-targeted. Reductions in all other groups reflect transmission-mediated population benefit, not a group-specific NNV.",
                   theme = theme(plot.caption = element_text(size = 6, hjust = 0)))
export_figure(fig_nnv5, "fig_nnv5_age_specific_benefit", h = 150)

# ================================================================
# FIGURE NNV-3: baseline susceptibility vs NNV (descriptive) -- unchanged
# design, number formatting fixed defensively.
# ================================================================
susc_note <- if (all(nnv_susceptibility$S_total_N_baseline == 1, na.rm = TRUE)) {
  "Baseline S/N = 1.0 for every draw (2015 predates CE's first major epidemic) -- no cross-draw variation to relate to NNV in this window"
} else "Descriptive only -- not a causal estimate"
panelA3 <- ggplot(nnv_susceptibility, aes(S_total_N_baseline, NNV_infection)) + geom_point(alpha = 0.4, size = 0.8, colour = COL_INF) +
  scale_x_continuous(labels = comma_dec) + scale_y_continuous(labels = comma_dec) +
  labs(title = "A", subtitle = "Overall S/N at programme start", x = "Baseline susceptible fraction (S/N)", y = "NNV per infection averted") + theme_journal()
panelB3 <- ggplot(nnv_susceptibility, aes(S_age12_N_baseline, NNV_infection)) + geom_point(alpha = 0.4, size = 0.8, colour = COL_REP) +
  scale_x_continuous(labels = comma_dec) + scale_y_continuous(labels = comma_dec) +
  labs(title = "B", subtitle = "Age-12 S/N at programme start", x = "Baseline susceptible fraction (S/N)", y = "NNV per infection averted") + theme_journal()
fig_nnv3 <- (panelA3 + panelB3) + plot_annotation(caption = susc_note, theme = theme(plot.caption = element_text(size = 6, hjust = 0)))
export_figure(fig_nnv3, "fig_nnv3_susceptibility_vs_nnv", h = 100)

message("\n[12_plot_nnv] All revised NNV figures saved to ", DIR_FIG_NNV, " (violin diagnostic to ", DIR_FIG_SUPP, ")")
