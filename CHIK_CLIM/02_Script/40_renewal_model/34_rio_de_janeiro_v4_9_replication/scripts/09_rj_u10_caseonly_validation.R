# Rio de Janeiro -- U10 (Perisse 2020, Rio de Janeiro CITY) case-only
# EXTERNAL VALIDATION. NO refit. Uses the existing HMC-PASS global-q
# case-only posterior (outputs/caseonly/rj_global_q_case_only.rds, now the
# adapt_delta=0.98 full-PASS result) plus the 5 already-fitted targeted
# fixed-q posteriors. U10 is a CITY survey -- this is external geographic
# plausibility only, NOT a state==city assumption, and does NOT update q.

required_packages <- c("rstan", "dplyr", "readr", "ggplot2", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(ggplot2); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication")
figure_dir <- file.path(root, "03_Output/figures/rio_de_janeiro_v4_9_replication")

U10_N_TESTED <- 2120L
U10_N_POSITIVE <- 371L
U10_WEIGHTED_PREV <- 0.180
U10_CI_LO <- 0.148
U10_CI_HI <- 0.212
U10_SURVEY_START <- as.Date("2018-07-01")
U10_SURVEY_END <- as.Date("2018-10-31")

b <- readRDS(file.path(base_dir, "outputs/caseonly/rj_global_q_case_only.rds"))
fit <- b$fit
dates <- as.Date(b$weekly_data$week_start)
window_idx <- which(dates >= U10_SURVEY_START & dates <= U10_SURVEY_END)
stopifnot(length(window_idx) > 0, all(diff(window_idx) == 1L))
message(sprintf("U10 window mapped to %d model weeks: %s to %s", length(window_idx),
                 as.character(dates[min(window_idx)]), as.character(dates[max(window_idx)])))

immune_prop <- rstan::extract(fit, "immune_prop")$immune_prop # [iter, N]
p_state_u10 <- rowMeans(immune_prop[, window_idx, drop = FALSE])

message("\n=== Section 1: case-only statewide immune fraction over the U10 window ===")
p_median <- median(p_state_u10); p_lo50 <- quantile(p_state_u10, .25); p_hi50 <- quantile(p_state_u10, .75)
p_lo95 <- quantile(p_state_u10, .025); p_hi95 <- quantile(p_state_u10, .975)
message(sprintf("p_state_U10: median=%.4f, 50%% CrI=[%.4f,%.4f], 95%% CrI=[%.4f,%.4f]", p_median, p_lo50, p_hi50, p_lo95, p_hi95))
message(sprintf("Observed U10 (Rio CITY): weighted prevalence=%.3f, 95%% CI=[%.3f,%.3f]", U10_WEIGHTED_PREV, U10_CI_LO, U10_CI_HI))

# ---- Section 2: quantify compatibility ----
prob_in_ci <- mean(p_state_u10 >= U10_CI_LO & p_state_u10 <= U10_CI_HI)
diff_draws <- U10_WEIGHTED_PREV - p_state_u10
message(sprintf("\nPosterior probability that p_state_U10 falls inside the OBSERVED CI [%.3f,%.3f]: %.4f", U10_CI_LO, U10_CI_HI, prob_in_ci))
message(sprintf("Difference (observed weighted - posterior p_state_U10): median=%.4f, 95%% CrI=[%.4f,%.4f]",
                 median(diff_draws), quantile(diff_draws, .025), quantile(diff_draws, .975)))

# raw-count comparison too (beta-binomial style plausibility, informational only)
raw_prev <- U10_N_POSITIVE / U10_N_TESTED
message(sprintf("(For reference: U10 raw prevalence = %d/%d = %.4f)", U10_N_POSITIVE, U10_N_TESTED, raw_prev))

# ---- Fixed-q comparison ----
q_grid <- c(0.0100, 0.0125, 0.0150, 0.0175, 0.0200)
fixedq_dir <- file.path(base_dir, "outputs/fixedq")
fixedq_results <- bind_rows(lapply(q_grid, function(q) {
  fpath <- file.path(fixedq_dir, sprintf("rj_fixedq_q%.4f.rds", q))
  bb <- readRDS(fpath)
  imm <- rstan::extract(bb$fit, "immune_prop")$immune_prop
  p_draws <- rowMeans(imm[, window_idx, drop = FALSE])
  tibble(q_target = q, p_state_u10_median = median(p_draws), p_state_u10_lo95 = quantile(p_draws, .025),
         p_state_u10_hi95 = quantile(p_draws, .975), prob_in_observed_ci = mean(p_draws >= U10_CI_LO & p_draws <= U10_CI_HI))
}))
message("\n=== Fixed-q comparison at the U10 window ===")
print(as.data.frame(fixedq_results), digits = 4)

globalq_row <- tibble(q_target = NA_real_, p_state_u10_median = p_median, p_state_u10_lo95 = p_lo95,
                       p_state_u10_hi95 = p_hi95, prob_in_observed_ci = prob_in_ci, model = "global-q (free, HMC PASS)")
fixedq_results$model <- sprintf("fixed q=%.4f", fixedq_results$q_target)
combined <- bind_rows(globalq_row, fixedq_results)
write_csv(combined, file.path(table_dir, "RJ_U10_fixedq_comparison.csv"))
message("\n[saved] ", file.path(table_dir, "RJ_U10_fixedq_comparison.csv"))

# Which fixed-q value has p_state_u10 closest to the observed weighted prevalence?
closest_q <- fixedq_results$q_target[which.min(abs(fixedq_results$p_state_u10_median - U10_WEIGHTED_PREV))]
message(sprintf("\nFixed-q value with p_state_U10 closest to observed (0.180): q=%.4f", closest_q))

# ---- Section 5: classification ----
# A. COMPATIBLE: observed CI and posterior CI overlap substantially / prob_in_ci not tiny
# B. MODERATE: some gap but same order of magnitude, plausibly city>state (city usually >= state for a capital)
# C. MAJOR CONFLICT: negligible overlap, orders-of-magnitude gap
overlap <- !(p_hi95 < U10_CI_LO || p_lo95 > U10_CI_HI)
gap_ratio <- U10_WEIGHTED_PREV / p_median
classification <- if (overlap && prob_in_ci > 0.05) {
  "A. COMPATIBLE"
} else if (gap_ratio < 3 && gap_ratio > 0.33) {
  "B. MODERATE GEOGRAPHIC DISCREPANCY"
} else {
  "C. MAJOR EXTERNAL-VALIDATION CONFLICT"
}
message(sprintf("\n95%% CrI overlap between posterior and observed CI: %s", overlap))
message(sprintf("Ratio (observed/posterior median): %.2f", gap_ratio))
message(sprintf("\n=== CLASSIFICATION: %s ===", classification))

# ---- Figure: overlay ----
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
overlay_df <- tibble(
  source = c("Case-only statewide posterior\n(U10 window, HMC PASS)", "U10 observed\n(Rio de Janeiro CITY)"),
  median = c(p_median, U10_WEIGHTED_PREV), lo = c(p_lo95, U10_CI_LO), hi = c(p_hi95, U10_CI_HI)
)
p_overlay <- ggplot(overlay_df, aes(source, median, colour = source)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.8, fatten = 3) +
  scale_colour_manual(values = c("Case-only statewide posterior\n(U10 window, HMC PASS)" = "#315A7D", "U10 observed\n(Rio de Janeiro CITY)" = "#D55E00")) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA)) +
  labs(title = "RJ case-only STATEWIDE posterior vs U10 CITY-level observed seroprevalence",
       subtitle = sprintf("Window: %s to %s. NOT a state=city equivalence claim -- external plausibility check only. Classification: %s",
                           as.character(U10_SURVEY_START), as.character(U10_SURVEY_END), classification),
       x = NULL, y = "Seroprevalence / immune fraction") +
  theme_v4 + theme(legend.position = "none") + coord_flip()
ggsave(file.path(figure_dir, "RJ_U10_caseonly_overlay.png"), p_overlay, width = 190, height = 90, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "RJ_U10_caseonly_overlay.png"))
