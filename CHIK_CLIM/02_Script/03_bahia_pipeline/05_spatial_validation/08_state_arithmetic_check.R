# Bahia global-q local-consistency audit -- Sections 3 & 4.
#
# Determines how much of the low state-wide cumulative-infection estimate
# is simple arithmetic (reported cases / q / population) versus a Stan-
# specific artefact.

required_packages <- c("rstan", "dplyr", "readr", "tibble", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
table_dir <- file.path(root, "03_Output/03_bahia_pipeline/tables/bahia_global_q_local_consistency_audit")
figure_dir <- file.path(root, "03_Output/03_bahia_pipeline/figures/bahia_global_q_local_consistency_audit")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- readRDS(file.path(table_dir, "audit_inputs_bundle.rds"))
weekly_data <- inputs$fit_bundle$weekly_data
q_draws <- inputs$q_draws
X_draws <- inputs$X_draws
S_draws <- inputs$S_draws
U_draws <- inputs$U_draws

C_total <- sum(weekly_data$cases)
N_start <- weekly_data$N_start[1]
N_end <- weekly_data$N_end[nrow(weekly_data)]
N_mean <- mean((weekly_data$N_start + weekly_data$N_end) / 2)

message(sprintf("C_total (total reported Bahia cases, 2015-2025) = %d", C_total))
message(sprintf("N_start = %.0f | N_mean = %.0f | N_end = %.0f", N_start, N_mean, N_end))

# ---- Section 3: posterior-propagated crude arithmetic vs Stan quantities ----
X_crude_draw <- C_total / q_draws
A_crude_start <- X_crude_draw / N_start
A_crude_mean <- X_crude_draw / N_mean
A_crude_end <- X_crude_draw / N_end

sum_X_draw <- rowSums(X_draws)
N_T <- S_draws[, ncol(S_draws)] + U_draws[, ncol(U_draws)]
final_immune_fraction <- U_draws[, ncol(U_draws)] / N_T
sum_X_over_N <- sum_X_draw / N_T

summarise_q <- function(x, label) {
  tibble(quantity = label, median = median(x), lo50 = quantile(x, .25), hi50 = quantile(x, .75),
         lo95 = quantile(x, .025), hi95 = quantile(x, .975))
}

section3_summary <- bind_rows(
  summarise_q(A_crude_start, "C_total / q / N_start (crude)"),
  summarise_q(A_crude_mean, "C_total / q / N_mean (crude)"),
  summarise_q(A_crude_end, "C_total / q / N_end (crude)"),
  summarise_q(sum_X_over_N, "sum(X)/N_T (posterior, Stan)"),
  summarise_q(final_immune_fraction, "final immune fraction U[T]/(S[T]+U[T]) (posterior, Stan)")
)
message("\n=== Section 3: crude arithmetic vs posterior Stan quantities ===")
print(as.data.frame(section3_summary), digits = 4)

ratio_sumX_vs_crude <- sum_X_draw / X_crude_draw
ratio_immune_vs_crudeA <- final_immune_fraction / A_crude_mean

ratio_summary <- bind_rows(
  summarise_q(ratio_sumX_vs_crude, "sum(X) / (C_total/q)"),
  summarise_q(ratio_immune_vs_crudeA, "final_immune_fraction / (C_total/q/N_mean)")
)
message("\n=== Ratios (near 1.0 => low cumulative-infection result is mainly arithmetic, not a Stan artefact) ===")
print(as.data.frame(ratio_summary), digits = 4)

write_csv(bind_rows(
  section3_summary |> mutate(section = "3_state_arithmetic_vs_posterior"),
  ratio_summary |> mutate(section = "3_ratios")
), file.path(table_dir, "bahia_state_arithmetic_q_audit.csv"))

# ---- Section 4: fixed-q arithmetic comparison table ------------------------
q_grid <- c(0.05, 0.10, 0.136, 0.15, 0.20, 0.25, 0.30, 0.359)
fixed_q_table <- tibble(
  q = q_grid,
  total_implied_infections = C_total / q_grid,
  implied_statewide_attack_fraction_Nmean = (C_total / q_grid) / N_mean,
  implied_statewide_susceptible_fraction_Nmean = 1 - (C_total / q_grid) / N_mean
)
message("\n=== Section 4: fixed-q arithmetic scale diagnostic (NOT a fitted comparison) ===")
print(as.data.frame(fixed_q_table), digits = 4)
write_csv(fixed_q_table, file.path(table_dir, "bahia_state_arithmetic_q_audit_fixed_grid.csv"))

# ---- Figure: state attack fraction vs q ------------------------------------
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
q_seq <- seq(0.02, 0.6, length.out = 300)
curve_df <- tibble(q = q_seq, attack_fraction = (C_total / q_seq) / N_mean)

p <- ggplot(curve_df, aes(q, attack_fraction)) +
  geom_line(colour = "#315A7D", linewidth = 0.8) +
  geom_point(data = fixed_q_table, aes(q, implied_statewide_attack_fraction_Nmean), colour = "#D55E00", size = 2) +
  geom_text(data = fixed_q_table, aes(q, implied_statewide_attack_fraction_Nmean, label = sprintf("q=%.3f", q)),
            vjust = -0.8, size = 2.5, colour = "#D55E00") +
  geom_vline(xintercept = median(q_draws), linetype = 2, colour = "#76558F") +
  annotate("text", x = median(q_draws), y = max(curve_df$attack_fraction) * 0.9,
           label = sprintf("posterior median q=%.3f", median(q_draws)), colour = "#76558F", angle = 90, vjust = -0.5, size = 3) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "Bahia: crude statewide attack fraction vs assumed q",
       subtitle = sprintf("C_total = %d reported cases / q / N_mean = %.0f (arithmetic only, no transmission dynamics)", C_total, N_mean),
       x = "q (ascertainment fraction)", y = "Implied cumulative attack fraction") +
  theme_v4
ggsave(file.path(figure_dir, "bahia_state_attack_fraction_vs_q.png"), p, width = 200, height = 130, units = "mm", dpi = 300)
message("\n[saved] ", file.path(figure_dir, "bahia_state_attack_fraction_vs_q.png"))
