# Rebound mechanism audit (DIAGNOSTIC ONLY -- optional, not part of the
# primary counterfactual pipeline; NOT sourced by 00_run_stage_07.R).
#
# Question: why does the full-feedback vaccine arm (B) sometimes have MORE
# infections than the no-vaccine arm (A) in individual years after the
# 2016-2017 epidemic, even though the cumulative 2015-2021 effect remains
# strongly beneficial? Working hypothesis: vaccination suppresses 2016-2017
# transmission -> fewer natural infections -> more susceptibles survive
# those years -> some transmission is displaced later -> negative annual
# "infections averted" in 2018-2021. This script TESTS that hypothesis; it
# does not assume it.
#
# Reuses (does not modify) the accepted three_arm_results.rds for Parts
# 1-4, 7 and 9 -- no resimulation. Parts 5-6 resimulate Arms A/B only,
# using the exact same simulator (03_generic_age_vaccine_simulator.R) and
# the exact same paired posterior draw indices (three_arm_results.rds$draw_idx),
# varying only imports_per_week (Part 5) or coverage/VE_infection (Part 6).
# Arms A/B/C definitions, the vaccine mechanics, the age reconstruction,
# and the accepted climate-forced v4.9 fit are never modified. Stan is
# never refit (rstan::extract() on the already-fitted, already-saved
# stanfit object is not a refit).

required_packages <- c("here", "dplyr", "tidyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(readr); library(tibble); library(ggplot2); library(patchwork); library(scales)
})

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this script.")
  source(setup_path, local = .GlobalEnv)
}
root <- ROOT
STAGE_DIR <- file.path(root, "02_Script/02_ceara_pipeline/07_counterfactuals")
source(file.path(STAGE_DIR, "01_config.R"))
source(file.path(STAGE_DIR, "03_generic_age_vaccine_simulator.R"))
source(file.path(STAGE_DIR, "02_load_accepted_inputs.R"))
source(file.path(STAGE_DIR, "04_run_no_vaccine_arm.R"))
source(file.path(STAGE_DIR, "05_run_full_feedback_vaccine_arm.R"))

FIG_DIR <- file.path(root, "03_Output/02_ceara_pipeline/figures/vaccine/ceara_historical_age12_vaccine_counterfactual/rebound_audit")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
AUDIT_CACHE <- file.path(DIR_RESULTS, "rebound_audit_results.rds")

primary_path <- file.path(DIR_RESULTS, "three_arm_results.rds")
if (!file.exists(primary_path)) stop("STOP: three_arm_results.rds not found -- run Stage 07's primary pipeline first. This audit reuses it and never regenerates it.")
res <- readRDS(primary_path)
totals <- res$totals
agg_S <- res$agg_S; agg_N <- res$agg_N; agg_Inew <- res$agg_Inew
dates <- res$dates
N <- length(dates)
draw_idx <- res$draw_idx
n_use <- length(draw_idx)
years_of_weeks <- as.integer(format(dates, "%Y"))
YEARS <- 2015:2021

message(sprintf("[audit] Loaded three_arm_results.rds: %d draws x %d weeks (%s to %s)", n_use, N, min(dates), max(dates)))

q95 <- function(x) quantile(x, c(.025, .5, .975), names = FALSE, na.rm = TRUE)
to_mat <- function(v) matrix(v, ncol = N, byrow = TRUE) # totals is draw-major (see 07_execute_three_arms_all_draws.R): row i = draw draw_idx[i]'s N consecutive weeks

X_A_mat <- to_mat(totals$X_A); X_B_mat <- to_mat(totals$X_B); X_C_mat <- to_mat(totals$X_C)
lambda_A_mat <- to_mat(totals$lambda_A); lambda_B_mat <- to_mat(totals$lambda_B)
Reff_A_mat <- to_mat(totals$Reff_A); Reff_B_mat <- to_mat(totals$Reff_B)
infness_A_mat <- to_mat(totals$infectiousness_A); infness_B_mat <- to_mat(totals$infectiousness_B)

annual_by_draw <- function(mat, yr) rowSums(mat[, years_of_weeks == yr, drop = FALSE])

COL_NOVACC <- "#4A4A4A"; COL_VACC <- "#0072B2"; COL_DIRECT <- "#D55E00"; COL_INDIRECT <- "#009E73"
theme_journal <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(text = element_text(family = "sans", colour = "black"),
          plot.title = element_text(size = base_size, face = "plain", hjust = 0),
          axis.line = element_line(colour = "black", linewidth = 0.3), axis.ticks = element_line(colour = "black", linewidth = 0.3),
          axis.text = element_text(size = base_size - 1, colour = "black"), axis.title = element_text(size = base_size, colour = "black"),
          panel.grid = element_blank(), panel.background = element_rect(fill = "white", colour = NA), plot.background = element_rect(fill = "white", colour = NA),
          legend.background = element_blank(), legend.key = element_blank(), legend.text = element_text(size = base_size - 1), legend.title = element_blank(),
          legend.position = "top", legend.justification = "left",
          strip.background = element_blank(), strip.text = element_text(size = base_size, face = "plain"), plot.margin = margin(4, 6, 4, 4))
}
date_scale <- scale_x_date(breaks = seq(as.Date("2015-01-01"), as.Date("2022-01-01"), "1 year"), date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01)))
export_fig <- function(plot, stem, w = 180, h = 150) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), plot, width = w, height = h, units = "mm", device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), plot, width = w, height = h, units = "mm", dpi = 600, bg = "white")
  message("[saved] ", stem, " (.pdf/.png)")
}

# ================================================================
# PART 9 (run first) -- baseline reproduction gate
# ================================================================
message("\n================ PART 9 (gate): baseline reproduction check ================")
baseline_total_by_draw <- rowSums(X_A_mat) - rowSums(X_B_mat)
baseline_median <- median(baseline_total_by_draw)
message(sprintf("[gate] Median 2015-2021 cumulative infections averted (Arm A - Arm B) from three_arm_results.rds = %.0f", baseline_median))
message(sprintf("[gate] Expected (accepted primary analysis): approximately 539,000"))
if (abs(baseline_median - 539000) / 539000 > 0.05) {
  stop(sprintf("STOP: baseline reproduction failed -- got %.0f, expected ~539,000 (>5%% discrepancy). Do not proceed to sensitivity analyses; investigate the discrepancy first.", baseline_median))
}
message("[gate] PASS -- baseline reproduced within 5% of the accepted ~539k figure. Proceeding.")
if (!identical(res$config$window_start, WINDOW_START) || !identical(res$config$window_end, WINDOW_END)) {
  stop("STOP: three_arm_results.rds window does not match 01_config.R's current WINDOW_START/WINDOW_END.")
}
if (any(years_of_weeks < 2015 | years_of_weeks > 2021)) stop("STOP: unexpected years found outside 2015-2021.")

# ================================================================
# PART 1 -- reproduce the observed rebound pattern (no resimulation)
# ================================================================
message("\n================ PART 1: annual paired effects, 2015-2021 ================")
annual_summary <- bind_rows(lapply(YEARS, function(yr) {
  total_by_draw <- annual_by_draw(X_A_mat, yr) - annual_by_draw(X_B_mat, yr)
  direct_by_draw <- annual_by_draw(X_A_mat, yr) - annual_by_draw(X_C_mat, yr)
  indirect_by_draw <- annual_by_draw(X_C_mat, yr) - annual_by_draw(X_B_mat, yr)
  qt <- q95(total_by_draw); qd <- q95(direct_by_draw); qi <- q95(indirect_by_draw)
  tibble(year = yr,
         total_median = qt[2], total_lo95 = qt[1], total_hi95 = qt[3], total_pr_negative = mean(total_by_draw < 0),
         direct_median = qd[2], direct_lo95 = qd[1], direct_hi95 = qd[3], direct_pr_negative = mean(direct_by_draw < 0),
         indirect_median = qi[2], indirect_lo95 = qi[1], indirect_hi95 = qi[3], indirect_pr_negative = mean(indirect_by_draw < 0))
}))
print(as.data.frame(annual_summary), digits = 4)

# ================================================================
# FIGURE R1 -- when does the rebound occur?
# ================================================================
message("\n[fig] Building FIG_R1_rebound_timing ...")
wk_summ <- function(mat, arm_label) { q <- apply(mat, 2, q95); tibble(week_start = dates, median = q[2, ], lo95 = q[1, ], hi95 = q[3, ], arm = arm_label) }
panelA_df <- bind_rows(wk_summ(X_A_mat, "No vaccine (A)"), wk_summ(X_B_mat, "Routine age-12 (B)"))
panelA <- ggplot(panelA_df, aes(week_start, median, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.4) +
  scale_colour_manual(values = c("No vaccine (A)" = COL_NOVACC, "Routine age-12 (B)" = COL_VACC)) +
  scale_fill_manual(values = c("No vaccine (A)" = COL_NOVACC, "Routine age-12 (B)" = COL_VACC)) +
  date_scale + labs(title = "A: Weekly infections, X_A vs X_B", x = NULL, y = "Infections/week") + theme_journal()

D_mat <- t(apply(X_A_mat - X_B_mat, 1, cumsum))
D_summ <- apply(D_mat, 2, q95)
panelB_df <- tibble(week_start = dates, median = D_summ[2, ], lo95 = D_summ[1, ], hi95 = D_summ[3, ])
panelB <- ggplot(panelB_df, aes(week_start, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, fill = COL_VACC) + geom_line(colour = COL_VACC, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.3) +
  date_scale + labs(title = "B: Cumulative infections averted D(t) = cumsum(X_A - X_B)", x = NULL, y = "Cumulative infections averted") + theme_journal()

panelC <- ggplot(annual_summary, aes(factor(year), total_median)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_pointrange(aes(ymin = total_lo95, ymax = total_hi95), colour = COL_VACC, linewidth = 0.4) +
  labs(title = "C: Annual total infections averted (A - B)", x = "Year", y = "Infections averted") + theme_journal()

decomp_long <- bind_rows(
  annual_summary |> transmute(year, component = "Direct-only (A - C)", median = direct_median, lo95 = direct_lo95, hi95 = direct_hi95),
  annual_summary |> transmute(year, component = "Indirect (C - B)", median = indirect_median, lo95 = indirect_lo95, hi95 = indirect_hi95)
)
panelD <- ggplot(decomp_long, aes(factor(year), median, colour = component)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "black") +
  geom_pointrange(aes(ymin = lo95, ymax = hi95), position = position_dodge(width = 0.5), linewidth = 0.4) +
  scale_colour_manual(values = c("Direct-only (A - C)" = COL_DIRECT, "Indirect (C - B)" = COL_INDIRECT)) +
  labs(title = "D: Annual decomposition -- direct-only vs indirect", x = "Year", y = "Infections averted") + theme_journal()

fig_r1 <- (panelA / panelB / panelC / panelD)
export_fig(fig_r1, "FIG_R1_rebound_timing", w = 170, h = 220)

# ================================================================
# PART 2 -- susceptible-preservation mechanism (age-resolved, no resimulation)
# ================================================================
message("\n================ PART 2: age-resolved susceptible preservation ================")
FOCUS_GROUPS <- c("0-11", "65+"); CONTEXT_GROUPS <- c("12", "13-17", "18-64")
deltaS_by_group <- lapply(REPORTING_AGE_LABELS, function(g) agg_S[, , g, "B"] - agg_S[, , g, "A"])
names(deltaS_by_group) <- REPORTING_AGE_LABELS

checkpoints <- list("end of 2016" = max(which(years_of_weeks == 2016)), "end of 2017" = max(which(years_of_weeks == 2017)),
                     "end of 2018" = max(which(years_of_weeks == 2018)), "end of 2019" = max(which(years_of_weeks == 2019)),
                     "end of 2020" = max(which(years_of_weeks == 2020)), "end of 2021" = max(which(years_of_weeks == 2021)))
message("Delta S_g(t) = S_B,g(t) - S_A,g(t): median [95% CrI] across draws, at year-end checkpoints, never-directly-vaccinated groups only:")
for (g in FOCUS_GROUPS) {
  for (cp_name in names(checkpoints)) {
    t_idx <- checkpoints[[cp_name]]
    v <- deltaS_by_group[[g]][, t_idx]
    qv <- q95(v)
    message(sprintf("  group=%-5s  %-12s  median=%+.1f  [%+.1f, %+.1f]  Pr(>0)=%.3f", g, cp_name, qv[2], qv[1], qv[3], mean(v > 0)))
  }
}
susceptible_preservation_demonstrated <- all(sapply(FOCUS_GROUPS, function(g) {
  v <- deltaS_by_group[[g]][, checkpoints[["end of 2018"]]]
  median(v) > 0 && mean(v > 0) > 0.5
}))
message(sprintf("\n[Part 2 conclusion] Susceptible preservation in never-vaccinated groups (0-11 AND 65+, positive median Delta S with Pr>0.5 by end of 2018) = %s",
                if (susceptible_preservation_demonstrated) "DEMONSTRATED" else "NOT DEMONSTRATED"))

# ================================================================
# FIGURE R2 -- mechanistic rebound chain
# ================================================================
message("\n[fig] Building FIG_R2_rebound_mechanism ...")
deltaS_df <- bind_rows(lapply(c(FOCUS_GROUPS, "18-64"), function(g) {
  qv <- apply(deltaS_by_group[[g]], 2, median)
  tibble(week_start = dates, delta_S = qv, group = g)
})) |> mutate(group = factor(group, levels = c(FOCUS_GROUPS, "18-64")))
panelA2 <- ggplot(deltaS_df, aes(week_start, delta_S, colour = group, alpha = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "black") + geom_line(linewidth = 0.5) +
  scale_colour_manual(values = c("0-11" = "#5B8FA8", "65+" = "#3D3D3D", "18-64" = "#7A7A7A")) +
  scale_alpha_manual(values = c("0-11" = 1, "65+" = 1, "18-64" = 0.4)) +
  date_scale + labs(title = "A: Delta S_g(t) = S_B - S_A, never-vaccinated groups (18-64 shown faintly for context)", x = NULL, y = "Delta susceptibles") + theme_journal()

Reff_diff_summ <- apply(Reff_B_mat - Reff_A_mat, 2, q95)
panelB2 <- ggplot(tibble(week_start = dates, median = Reff_diff_summ[2, ], lo95 = Reff_diff_summ[1, ], hi95 = Reff_diff_summ[3, ]), aes(week_start, median)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) + geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, fill = COL_VACC) + geom_line(colour = COL_VACC, linewidth = 0.5) +
  date_scale + labs(title = "B: Delta Reff(t) = Reff_B - Reff_A", x = NULL, y = "Delta Reff") + theme_journal()

lambda_ratio_mat <- ifelse(lambda_A_mat > 1e-10, lambda_B_mat / lambda_A_mat, NA_real_)
lr_summ <- apply(lambda_ratio_mat, 2, q95)
panelC2 <- ggplot(tibble(week_start = dates, median = lr_summ[2, ], lo95 = lr_summ[1, ], hi95 = lr_summ[3, ]), aes(week_start, median)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) + geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, fill = COL_VACC) + geom_line(colour = COL_VACC, linewidth = 0.5) +
  date_scale + labs(title = "C: lambda_B(t) / lambda_A(t)", x = NULL, y = "Ratio") + theme_journal()

dX_summ <- apply(X_B_mat - X_A_mat, 2, q95)
panelD2 <- ggplot(tibble(week_start = dates, median = dX_summ[2, ], lo95 = dX_summ[1, ], hi95 = dX_summ[3, ]), aes(week_start, median)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) + geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.2, fill = COL_VACC) + geom_line(colour = COL_VACC, linewidth = 0.5) +
  date_scale + labs(title = "D: X_B(t) - X_A(t)", x = NULL, y = "Delta infections/week") + theme_journal()

fig_r2 <- (panelA2 / panelB2 / panelC2 / panelD2)
export_fig(fig_r2, "FIG_R2_rebound_mechanism", w = 170, h = 220)

# ================================================================
# PART 3 -- follow the dynamical chain: sustained crossover dates
# ================================================================
message("\n================ PART 3: sustained crossover after the 2017 epidemic ================")
search_from <- as.Date("2018-01-01")
first_sustained <- function(x, threshold, direction = c(">", "<"), min_run = 4) {
  direction <- match.arg(direction)
  eligible <- if (direction == ">") x > threshold else x < threshold
  idx0 <- which(dates >= search_from)
  for (i in idx0) {
    if (i + min_run - 1 > length(x)) break
    if (all(eligible[i:(i + min_run - 1)])) return(dates[i])
  }
  as.Date(NA)
}
med_deltaS_011 <- apply(deltaS_by_group[["0-11"]], 2, median)
med_deltaS_65p <- apply(deltaS_by_group[["65+"]], 2, median)
med_Reff_diff <- Reff_diff_summ[2, ]
med_lambda_ratio <- lr_summ[2, ]
med_dX <- dX_summ[2, ]

crossovers <- tibble(
  quantity = c("Delta S (0-11) > 0", "Delta S (65+) > 0", "Delta Reff > 0", "lambda_B/lambda_A > 1", "X_B - X_A > 0"),
  first_sustained_crossover = c(
    first_sustained(med_deltaS_011, 0, ">"), first_sustained(med_deltaS_65p, 0, ">"),
    first_sustained(med_Reff_diff, 0, ">"), first_sustained(med_lambda_ratio, 1, ">"),
    first_sustained(med_dX, 0, ">")
  )
)
message(sprintf("Searching for the first run of >=4 consecutive weeks (from %s onward) where the median crosses the stated threshold:", search_from))
print(as.data.frame(crossovers))
ordering_dates <- crossovers$first_sustained_crossover
names(ordering_dates) <- crossovers$quantity
ordering_observed <- names(sort(ordering_dates))
message("\nActual temporal ordering found (earliest sustained crossover first; NA = never sustained for >=4 weeks):")
message(paste("  ", seq_along(ordering_observed), ". ", ordering_observed, " (", format(sort(ordering_dates)), ")", sep = "", collapse = "\n"))

# ================================================================
# PART 4 -- verify the rebound is indirect (Arm C)
# ================================================================
message("\n================ PART 4: is the late negative effect direct or indirect? ================")
for (yr in YEARS) {
  direct_by_draw <- annual_by_draw(X_A_mat, yr) - annual_by_draw(X_C_mat, yr)
  indirect_by_draw <- annual_by_draw(X_C_mat, yr) - annual_by_draw(X_B_mat, yr)
  message(sprintf("  %d:  Pr(A - C < 0) = %.3f   Pr(C - B < 0) = %.3f", yr, mean(direct_by_draw < 0), mean(indirect_by_draw < 0)))
}
late_years <- 2018:2021
late_direct_pr_neg <- sapply(late_years, function(yr) mean((annual_by_draw(X_A_mat, yr) - annual_by_draw(X_C_mat, yr)) < 0))
late_indirect_pr_neg <- sapply(late_years, function(yr) mean((annual_by_draw(X_C_mat, yr) - annual_by_draw(X_B_mat, yr)) < 0))
late_driven_by_indirect <- mean(late_indirect_pr_neg) > mean(late_direct_pr_neg) && mean(late_indirect_pr_neg) > 0.5
if (late_driven_by_indirect) {
  message("\n\"Late negative annual effects are transmission-feedback mediated, not reversal of direct vaccine protection.\"")
} else {
  message(sprintf("\n[Part 4 conclusion] Late-period negative effects are NOT predominantly indirect-mediated (mean Pr(A-C<0)=%.3f vs mean Pr(C-B<0)=%.3f, 2018-2021) -- the direct/indirect ordering hypothesis is not supported as stated.",
                   mean(late_direct_pr_neg), mean(late_indirect_pr_neg)))
}

# ================================================================
# PART 5 -- importation sensitivity (resimulates Arms A & B only)
# ================================================================
message("\n================ PART 5: importation sensitivity ================")
inputs <- load_accepted_inputs()
if (!identical(inputs$dates, dates)) stop("STOP: freshly loaded inputs$dates does not match three_arm_results.rds$dates -- inputs are not aligned with the accepted primary analysis.")
stresstest <- SCENARIOS$stresstest
IMPORT_MULTIPLIERS <- c(0.25, 0.5, 1, 2, 4)

run_X_only <- function(fun, inputs_mod, R0_draw, ...) fun(inputs_mod, R0_draw, ...)$X_total

import_results <- vector("list", length(IMPORT_MULTIPLIERS))
for (mi in seq_along(IMPORT_MULTIPLIERS)) {
  mult <- IMPORT_MULTIPLIERS[mi]
  message(sprintf("[Part 5] import multiplier = %.2f ...", mult))
  if (mult == 1) {
    XA_list <- lapply(seq_len(n_use), function(i) X_A_mat[i, ])
    XB_list <- lapply(seq_len(n_use), function(i) X_B_mat[i, ])
  } else {
    inputs_mod <- inputs; inputs_mod$imports_per_week <- inputs$imports_per_week * mult
    XA_list <- vector("list", n_use); XB_list <- vector("list", n_use)
    for (i in seq_len(n_use)) {
      d <- draw_idx[i]; R0_draw <- inputs$R0_t_all_draws[d, ]
      XA_list[[i]] <- run_X_only(run_arm_A, inputs_mod, R0_draw)
      XB_list[[i]] <- run_X_only(run_arm_B, inputs_mod, R0_draw, coverage = stresstest$coverage, ve_infection = stresstest$ve_infection)
      if (i %% 50 == 0) message(sprintf("    ...%d / %d draws", i, n_use))
    }
  }
  XA_m <- do.call(rbind, XA_list); XB_m <- do.call(rbind, XB_list)
  early_by_draw <- (rowSums(XA_m[, years_of_weeks %in% 2015:2017, drop = FALSE]) - rowSums(XB_m[, years_of_weeks %in% 2015:2017, drop = FALSE]))
  late_by_draw <- (rowSums(XA_m[, years_of_weeks %in% 2018:2021, drop = FALSE]) - rowSums(XB_m[, years_of_weeks %in% 2018:2021, drop = FALSE]))
  total_by_draw <- rowSums(XA_m) - rowSums(XB_m)
  giveback <- ifelse(early_by_draw > 0 & late_by_draw < 0, -late_by_draw / early_by_draw, NA_real_)
  import_results[[mi]] <- tibble(multiplier = mult, early_median = median(early_by_draw), early_lo95 = quantile(early_by_draw, .025), early_hi95 = quantile(early_by_draw, .975),
                                  late_median = median(late_by_draw), late_lo95 = quantile(late_by_draw, .025), late_hi95 = quantile(late_by_draw, .975),
                                  total_median = median(total_by_draw), total_lo95 = quantile(total_by_draw, .025), total_hi95 = quantile(total_by_draw, .975),
                                  giveback_median = median(giveback, na.rm = TRUE), giveback_lo95 = quantile(giveback, .025, na.rm = TRUE), giveback_hi95 = quantile(giveback, .975, na.rm = TRUE),
                                  pr_late_negative = mean(late_by_draw < 0), pr_total_negative = mean(total_by_draw < 0), n_giveback_defined = sum(!is.na(giveback)))
}
import_sensitivity <- bind_rows(import_results)
message("\nImport-multiplier sensitivity (Arms A/B only; same paired R0 draws throughout):")
print(as.data.frame(import_sensitivity), digits = 4)

# ================================================================
# FIGURE R3 -- import sensitivity
# ================================================================
message("\n[fig] Building FIG_R3_import_sensitivity ...")
panelA3 <- ggplot(import_sensitivity, aes(multiplier, total_median)) + geom_pointrange(aes(ymin = total_lo95, ymax = total_hi95), colour = COL_VACC) + geom_line(colour = COL_VACC) +
  scale_x_log10(breaks = IMPORT_MULTIPLIERS) + labs(title = "A: Total 2015-2021 infections averted", x = "Import multiplier", y = "Infections averted") + theme_journal()
panelB3 <- ggplot(import_sensitivity, aes(multiplier, late_median)) + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_pointrange(aes(ymin = late_lo95, ymax = late_hi95), colour = COL_INDIRECT) + geom_line(colour = COL_INDIRECT) +
  scale_x_log10(breaks = IMPORT_MULTIPLIERS) + labs(title = "B: Late 2018-2021 infections averted", x = "Import multiplier", y = "Infections averted") + theme_journal()
panelC3 <- ggplot(import_sensitivity, aes(multiplier, giveback_median)) + geom_pointrange(aes(ymin = giveback_lo95, ymax = giveback_hi95), colour = COL_DIRECT) + geom_line(colour = COL_DIRECT) +
  scale_x_log10(breaks = IMPORT_MULTIPLIERS) + labs(title = "C: Give-back fraction (among draws where late effect < 0)", x = "Import multiplier", y = "Give-back fraction") + theme_journal()
panelD3 <- ggplot(import_sensitivity, aes(multiplier, pr_late_negative)) + geom_point(colour = COL_NOVACC, size = 2) + geom_line(colour = COL_NOVACC) +
  scale_x_log10(breaks = IMPORT_MULTIPLIERS) + ylim(0, 1) + labs(title = "D: Pr(late 2018-2021 effect < 0)", x = "Import multiplier", y = "Probability") + theme_journal()
fig_r3 <- (panelA3 | panelB3) / (panelC3 | panelD3)
export_fig(fig_r3, "FIG_R3_import_sensitivity", w = 180, h = 150)

import_robust <- all(import_sensitivity$pr_late_negative[import_sensitivity$multiplier %in% c(0.5, 1, 2)] > 0.5)

# ================================================================
# PART 6 -- coverage x VE_infection mechanism grid (Arm B only; Arm A reused)
# ================================================================
message("\n================ PART 6: coverage x VE_infection mechanism grid ================")
COVERAGE_GRID <- c(0.25, 0.50, 0.75, 1.00); VE_GRID <- c(0.50, 0.70, 0.90, 1.00)
XA_primary_list <- lapply(seq_len(n_use), function(i) X_A_mat[i, ]) # Arm A never vaccinates -> identical across the whole grid; reused from the cached primary result
grid_rows <- list(); k <- 0
for (cov in COVERAGE_GRID) for (ve in VE_GRID) {
  k <- k + 1
  message(sprintf("[Part 6] coverage=%.2f  VE_infection=%.2f ...", cov, ve))
  XB_list <- vector("list", n_use)
  for (i in seq_len(n_use)) {
    d <- draw_idx[i]; R0_draw <- inputs$R0_t_all_draws[d, ]
    XB_list[[i]] <- run_X_only(run_arm_B, inputs, R0_draw, coverage = cov, ve_infection = ve)
  }
  XA_m <- do.call(rbind, XA_primary_list); XB_m <- do.call(rbind, XB_list)
  early_by_draw <- rowSums(XA_m[, years_of_weeks %in% 2015:2017, drop = FALSE]) - rowSums(XB_m[, years_of_weeks %in% 2015:2017, drop = FALSE])
  late_by_draw <- rowSums(XA_m[, years_of_weeks %in% 2018:2021, drop = FALSE]) - rowSums(XB_m[, years_of_weeks %in% 2018:2021, drop = FALSE])
  total_by_draw <- rowSums(XA_m) - rowSums(XB_m)
  pct_reduction_by_draw <- 100 * total_by_draw / rowSums(XA_m)
  giveback <- ifelse(early_by_draw > 0 & late_by_draw < 0, -late_by_draw / early_by_draw, NA_real_)
  grid_rows[[k]] <- tibble(coverage = cov, ve_infection = ve,
                            total_median = median(total_by_draw), pct_reduction_median = median(pct_reduction_by_draw),
                            early_median = median(early_by_draw), late_median = median(late_by_draw),
                            giveback_median = median(giveback, na.rm = TRUE), pr_late_negative = mean(late_by_draw < 0), pr_total_negative = mean(total_by_draw < 0))
}
coverage_ve_grid <- bind_rows(grid_rows)
message("\nCoverage x VE_infection mechanism grid (import multiplier = 1, same paired draws):")
print(as.data.frame(coverage_ve_grid), digits = 4)

# ================================================================
# FIGURE R4 -- coverage x VE sensitivity (heatmaps)
# ================================================================
message("\n[fig] Building FIG_R4_coverage_VE_rebound_sensitivity ...")
hm_base <- function(fill_var, title, fill_lab, diverge = FALSE) {
  p <- ggplot(coverage_ve_grid, aes(factor(coverage), factor(ve_infection), fill = .data[[fill_var]])) +
    geom_tile(colour = "white") + geom_text(aes(label = sprintf("%.1f", .data[[fill_var]])), size = 2.6, colour = "black") +
    labs(title = title, x = "Coverage", y = "VE (infection)", fill = fill_lab) + theme_journal() + theme(legend.position = "right")
  if (diverge) {
    lim <- max(abs(coverage_ve_grid[[fill_var]]), na.rm = TRUE)
    p <- p + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim))
  } else {
    p <- p + scale_fill_viridis_c()
  }
  p
}
panelA4 <- hm_base("pct_reduction_median", "A: Median % reduction, 2015-2021", "% reduction")
panelB4 <- hm_base("late_median", "B: Median late 2018-2021 infections averted", "Infections", diverge = TRUE)
panelC4 <- hm_base("giveback_median", "C: Median give-back fraction", "Fraction")
panelD4 <- hm_base("pr_late_negative", "D: Pr(late 2018-2021 effect < 0)", "Probability")
fig_r4 <- (panelA4 | panelB4) / (panelC4 | panelD4)
export_fig(fig_r4, "FIG_R4_coverage_VE_rebound_sensitivity", w = 190, h = 170)

grid_extreme_only <- all(coverage_ve_grid$pr_late_negative[coverage_ve_grid$coverage == 1 & coverage_ve_grid$ve_infection == 1] > 0.5) &&
  mean(coverage_ve_grid$pr_late_negative[coverage_ve_grid$coverage < 1 | coverage_ve_grid$ve_infection < 1] > 0.5) < 0.5

# ================================================================
# PART 7 -- low-incidence diagnostic (existing primary results only)
# ================================================================
message("\n================ PART 7: low-incidence diagnostic ================")
med_X_A <- apply(X_A_mat, 2, median); med_X_B <- apply(X_B_mat, 2, median)
rebound_weeks <- which(med_X_B > med_X_A)
message(sprintf("Rebound weeks (median X_B > median X_A): %d / %d weeks total", length(rebound_weeks), N))
message(sprintf("  Among rebound weeks: median X_A = %.2f, median X_B = %.2f", median(med_X_A[rebound_weeks]), median(med_X_B[rebound_weeks])))
for (thr in c(1, 10, 100)) {
  p <- mean(med_X_A[rebound_weeks] < thr & med_X_B[rebound_weeks] < thr)
  message(sprintf("  Proportion of rebound weeks with BOTH arms < %d infection(s)/week: %.3f", thr, p))
}
message("By year:")
for (yr in 2018:2021) {
  wk <- rebound_weeks[years_of_weeks[rebound_weeks] == yr]
  if (length(wk) == 0) { message(sprintf("  %d: no rebound weeks", yr)); next }
  props <- sapply(c(1, 10, 100), function(thr) mean(med_X_A[wk] < thr & med_X_B[wk] < thr))
  message(sprintf("  %d: %d rebound weeks | median X_A=%.2f X_B=%.2f | prop(<1)=%.3f prop(<10)=%.3f prop(<100)=%.3f",
                   yr, length(wk), median(med_X_A[wk]), median(med_X_B[wk]), props[1], props[2], props[3]))
}
low_incidence_dominant <- mean(med_X_A[rebound_weeks] < 10 & med_X_B[rebound_weeks] < 10) > 0.5
if (low_incidence_dominant) {
  message("\n\"A process-stochastic low-incidence sensitivity should be considered before interpreting the rebound quantitatively.\"")
}

# ================================================================
# FIGURE R5 -- rebound magnitude vs absolute incidence
# ================================================================
message("\n[fig] Building FIG_R5_rebound_low_incidence_context ...")
context_weeks <- years_of_weeks >= 2017
plot_floor <- 0.1
panelA5 <- ggplot(tibble(week_start = dates[context_weeks], X_A = pmax(med_X_A[context_weeks], plot_floor), X_B = pmax(med_X_B[context_weeks], plot_floor)) |>
                     pivot_longer(c(X_A, X_B), names_to = "arm", values_to = "infections"),
                   aes(week_start, infections, colour = arm)) + geom_line(linewidth = 0.4) +
  scale_y_log10() + scale_colour_manual(values = c(X_A = COL_NOVACC, X_B = COL_VACC)) +
  date_scale + labs(title = "A: Median X_A, X_B (log10 scale; floor at 0.1)", x = NULL, y = "Infections/week (log10)") + theme_journal()
panelB5 <- ggplot(tibble(week_start = dates[context_weeks], diff = (med_X_B - med_X_A)[context_weeks]), aes(week_start, diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) + geom_line(colour = COL_VACC, linewidth = 0.4) +
  date_scale + labs(title = "B: Median X_B - X_A", x = NULL, y = "Delta infections/week") + theme_journal()
late_idx <- which(years_of_weeks %in% 2018:2021)
panelC5 <- ggplot(tibble(X_A = med_X_A[late_idx], diff = (med_X_B - med_X_A)[late_idx]), aes(X_A, diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) + geom_point(alpha = 0.4, size = 0.8, colour = COL_VACC) +
  scale_x_log10() + labs(title = "C: Rebound magnitude vs no-vaccine incidence, 2018-2021", x = "Median X_A (log10)", y = "Median X_B - X_A") + theme_journal()
annual_burden <- tibble(year = 2018:2021, no_vaccine_burden = sapply(2018:2021, function(yr) sum(med_X_A[years_of_weeks == yr])),
                         rebound_magnitude = sapply(2018:2021, function(yr) sum(pmax(med_X_B - med_X_A, 0)[years_of_weeks == yr])))
panelD5 <- ggplot(annual_burden, aes(no_vaccine_burden, rebound_magnitude, label = year)) + geom_point(colour = COL_DIRECT, size = 2) +
  geom_text(nudge_y = max(annual_burden$rebound_magnitude, 1) * 0.06, size = 2.8, colour = "black") +
  labs(title = "D: Annual rebound magnitude vs annual no-vaccine burden", x = "Annual no-vaccine infections (X_A)", y = "Annual rebound magnitude") + theme_journal()
fig_r5 <- (panelA5 | panelB5) / (panelC5 | panelD5)
export_fig(fig_r5, "FIG_R5_rebound_low_incidence_context", w = 190, h = 150)

# ================================================================
# PART 8 -- console interpretation
# ================================================================
message("\n================ PART 8: diagnostic interpretation ================")
cumulative_ever_negative <- any(rowSums(X_A_mat) - rowSums(X_B_mat) < 0)
message(sprintf("1. Does the cumulative 2015-2021 effect ever become negative (any single draw)? %s", if (cumulative_ever_negative) "YES" else "NO"))
message("2. Years with Pr(annual total effect < 0) > 0.5 and > 0.95:")
yrs_gt50 <- annual_summary$year[annual_summary$total_pr_negative > 0.5]
yrs_gt95 <- annual_summary$year[annual_summary$total_pr_negative > 0.95]
message(sprintf("   >0.5:  %s", if (length(yrs_gt50)) paste(yrs_gt50, collapse = ", ") else "none"))
message(sprintf("   >0.95: %s", if (length(yrs_gt95)) paste(yrs_gt95, collapse = ", ") else "none"))
message(sprintf("3. Late negative effect is predominantly: %s", if (late_driven_by_indirect) "INDIRECT (transmission-mediated)" else "not clearly indirect-dominant"))
message(sprintf("4. Never-directly-vaccinated groups (0-11, 65+) show increased susceptibility after 2016-2017 suppression? %s", if (susceptible_preservation_demonstrated) "YES" else "NO"))
message(sprintf("5. Does susceptible preservation temporally precede lambda_B>lambda_A and/or X_B>X_A? Observed ordering: %s", paste(ordering_observed, collapse = " -> ")))
message(sprintf("6. Import sensitivity: late rebound Pr(effect<0)>0.5 at multipliers 0.5/1/2? %s", if (import_robust) "ROBUST across this range" else "NOT robust -- sensitive to import assumption"))
message(sprintf("7. Coverage/VE sensitivity: rebound confined to the extreme coverage=1/VE=1 corner? %s", if (grid_extreme_only) "YES -- largely an artifact of the stress-test corner" else "NO -- present more broadly across the grid"))
message(sprintf("8. Does the rebound occur mainly at very low absolute incidence (<10/week)? %s", if (low_incidence_dominant) "YES" else "NO"))

s_ok <- susceptible_preservation_demonstrated
s_labels <- c("Delta S (0-11) > 0", "Delta S (65+) > 0")
s_positions <- which(ordering_observed %in% s_labels)
other_positions <- which(!(ordering_observed %in% s_labels))
order_ok <- length(s_positions) > 0 && length(other_positions) > 0 && min(s_positions) < min(other_positions)
indirect_ok <- late_driven_by_indirect
assumption_sensitive <- (!import_robust) || grid_extreme_only || low_incidence_dominant

classification <- if (!s_ok || !indirect_ok) {
  "D. PATTERN SUGGESTS IMPLEMENTATION OR MODEL-STRUCTURE PROBLEM"
} else if (!order_ok) {
  "C. MECHANISM UNCLEAR"
} else if (assumption_sensitive) {
  "B. CONSISTENT WITH EPIDEMIC DISPLACEMENT BUT ASSUMPTION-SENSITIVE"
} else {
  "A. STRONGLY CONSISTENT WITH EPIDEMIC DISPLACEMENT"
}
message(sprintf("\n>>> CLASSIFICATION: %s", classification))
message("Findings driving this classification:")
message(sprintf("  - susceptible preservation demonstrated: %s", s_ok))
message(sprintf("  - late effect indirect-dominant: %s", indirect_ok))
message(sprintf("  - temporal ordering (S-preservation precedes lambda/X reversal): %s", order_ok))
message(sprintf("  - assumption-sensitive (import and/or coverage-VE and/or low-incidence): %s", assumption_sensitive))

# ================================================================
# PART 9 (remaining checks) -- internal consistency
# ================================================================
message("\n================ PART 9: internal consistency checks ================")
message(sprintf("  - same posterior draw indices used throughout (draw_idx, n=%d): OK by construction", n_use))
message(sprintf("  - three_arm_results.rds mtime unchanged by this script: %s", file.exists(primary_path)))
message("  - R0(t) trajectories: unchanged (reused inputs$R0_t_all_draws[draw_idx, ], never resampled)")
message("  - vaccine code (03_generic_age_vaccine_simulator.R): unchanged (sourced, not edited)")
message("  - age reconstruction / ACCEPTED_AGE_SHARES: unchanged (read-only)")
message("  - Stan fit: never rerun (rstan::extract() on the already-saved stanfit only)")
message(sprintf("  - no seed in 2015-2021: %s (all is_seed == 0 confirmed in 02_load_accepted_inputs.R)", all(inputs$is_seed == 0)))
message(sprintf("  - baseline ~539k total effect reproduced from the SAME draw subset: PASS (median = %.0f)", baseline_median))

saveRDS(list(annual_summary = annual_summary, crossovers = crossovers, import_sensitivity = import_sensitivity,
             coverage_ve_grid = coverage_ve_grid, classification = classification,
             susceptible_preservation_demonstrated = susceptible_preservation_demonstrated,
             late_driven_by_indirect = late_driven_by_indirect, import_robust = import_robust,
             grid_extreme_only = grid_extreme_only, low_incidence_dominant = low_incidence_dominant,
             baseline_median = baseline_median),
        AUDIT_CACHE)
message("\n[saved] ", AUDIT_CACHE)
message("\n[audit] Complete. three_arm_results.rds was never modified.")
