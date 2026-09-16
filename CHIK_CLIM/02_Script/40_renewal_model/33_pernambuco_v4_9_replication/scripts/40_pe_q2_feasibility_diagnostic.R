# Pernambuco Spatial-0-Q2 -- Section 1: cheap two-q feasibility diagnostic.
#
# NO Stan re-fit. A plain-R re-implementation of the EXACT Spatial-0
# deterministic S/U/X recursion (byte-equivalent to
# renewal_pernambuco_v4_9_spatial0.stan's transformed-parameters block),
# using the fitted model's representative (posterior-median) seasonal/year
# parameters, but GRIDDING alpha_R (the shared transmission-scale
# intercept) to ask: is there ANY common transmission scale at which
# Recife's own case series and observed U14 serology, and the other 4
# regions' case series, can ALL be reconciled via differential
# ascertainment (q_Recife, q_rest) alone -- without touching R0(t), S/U
# bookkeeping, or biological plausibility?
#
# Key structural fact this diagnostic exploits: in this model, q is a PURE
# observation-scaling parameter (expected_reported_cases = q*X) that does
# NOT feed back into the S/U/X recursion at all. So for a FIXED alpha_R
# (hence fixed R0(t), fixed X/S/U/immune_prop trajectories for every
# stratum), the case fit for each region can be assessed independently by
# choosing whatever q best explains that region's OWN observed cases,
# while the U14 fit is checked directly against immune_prop -- q cannot
# move immune_prop at all. This lets us test coherence WITHOUT a full
# joint Bayesian re-fit.

required_packages <- c("dplyr", "readr", "ggplot2", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(tidyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
pe_root <- file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_model")
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/spatial0_q2")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

STRATUM_LEVELS <- c("1_Recife", "2_Metropolitana_remainder", "3_Agreste", "4_Sertao", "5_Vale_Sao_Francisco_Araripe")
RECIFE_INDEX <- 1L

built <- readRDS(file.path(table_dir, "PE_5strata_stan_data.rds"))
sd <- built$stan_data
rep_draw <- readRDS(file.path(table_dir, "PE_spatial0_representative_draw.rds"))

N <- sd$N; R <- sd$R; G <- sd$G

# ---- Deterministic S/U/X recursion, EXACT translation of the Stan block ----
run_recursion <- function(alpha_R) {
  log_R0 <- alpha_R + rep_draw$year_effect[sd$year_id] +
    rep_draw$A_year[sd$year_id] * (rep_draw$beta_sin1 * sd$seasonal_sin + rep_draw$beta_cos1 * sd$seasonal_cos) +
    rep_draw$beta_sin2 * sd$seasonal_sin2 + rep_draw$beta_cos2 * sd$seasonal_cos2
  R0_t <- exp(log_R0)

  S <- matrix(0, N, R); U <- matrix(0, N, R); X <- matrix(0, N, R)
  S[1, ] <- sd$N_start[1, ]; U[1, ] <- 0
  reconciliation <- sd$N_end - sd$N_start - sd$births + sd$deaths

  for (t in 1:N) {
    total <- S[t, ] + U[t, ]
    infectiousness <- rep(0, R)
    for (g in 1:G) if (t > g) infectiousness <- infectiousness + sd$w[g] * X[t - g, ]
    foi <- R0_t[t] * infectiousness / total + sd$imports_r[t, ] / total
    X[t, ] <- S[t, ] * (-expm1(-foi))
    if (t < N) {
      susc_after <- S[t, ] - X[t, ]
      imm_after <- U[t, ] + X[t, ]
      susc_frac <- susc_after / (susc_after + imm_after)
      S[t + 1, ] <- susc_after + sd$births[t, ] - sd$deaths[t, ] * susc_frac + reconciliation[t, ] * susc_frac
      U[t + 1, ] <- imm_after - sd$deaths[t, ] * (1 - susc_frac) + reconciliation[t, ] * (1 - susc_frac)
    }
  }
  immune_prop <- U / (S + U)
  S_prop <- S / (S + U)
  list(X = X, S = S, U = U, immune_prop = immune_prop, S_prop = S_prop, R0_t = R0_t)
}

sero_idx <- sd$sero_window_start_idx:(sd$sero_window_start_idx + sd$sero_window_n_weeks - 1)
observed_C <- sd$C
rest_idx <- setdiff(1:R, RECIFE_INDEX)

best_q <- function(observed, latent) sum(observed) / sum(latent) # OLS-style scale (Poisson/NB2 MLE for pure scale factor)

alpha_grid <- seq(rep_draw$alpha_R_fitted - 0.3, rep_draw$alpha_R_fitted + 1.5, by = 0.1)
message("Fitted (shared-q) alpha_R = ", round(rep_draw$alpha_R_fitted, 4), ". Grid: ",
        round(min(alpha_grid), 3), " to ", round(max(alpha_grid), 3), " (", length(alpha_grid), " points)")

results <- bind_rows(lapply(alpha_grid, function(a) {
  rec <- run_recursion(a)
  q_recife_opt <- best_q(observed_C[, RECIFE_INDEX], rec$X[, RECIFE_INDEX])
  q_rest_opt <- best_q(rowSums(observed_C[, rest_idx]), rowSums(rec$X[, rest_idx]))
  u14_prevalence <- mean(rec$immune_prop[sero_idx, RECIFE_INDEX])
  tibble(
    alpha_R = a, max_R0 = max(rec$R0_t),
    q_recife_optimal = q_recife_opt, q_rest_optimal = q_rest_opt,
    u14_implied_prevalence = u14_prevalence,
    min_S_prop_any_stratum = min(rec$S_prop),
    cumulative_infected_2025_recife = rec$immune_prop[N, RECIFE_INDEX],
    cumulative_infected_2025_max_stratum = max(rec$immune_prop[N, ])
  )
}))

write_csv(results, file.path(table_dir, "PE_Q2_feasibility_grid.csv"))
message("[saved] ", file.path(table_dir, "PE_Q2_feasibility_grid.csv"))
message("\n=== Feasibility grid ===")
print(as.data.frame(results), digits = 4)

# ---- Feasibility assessment ----------------------------------------------
in_range <- results |> mutate(
  q_recife_in_range = q_recife_optimal >= 0.005 & q_recife_optimal <= 0.03,
  q_rest_in_range = q_rest_optimal >= 0.03 & q_rest_optimal <= 0.12,
  u14_close = abs(u14_implied_prevalence - 0.372) <= 0.05, # within 5pp of observed
  u14_in_reported_ci = u14_implied_prevalence >= 0.340 & u14_implied_prevalence <= 0.404,
  biologically_plausible = max_R0 <= 10 & min_S_prop_any_stratum >= 0
)
feasible <- in_range |> filter(q_recife_in_range, q_rest_in_range, u14_close, biologically_plausible)

message("\n=== In-range flags ===")
print(as.data.frame(in_range |> select(alpha_R, q_recife_optimal, q_rest_optimal, u14_implied_prevalence,
                                         q_recife_in_range, q_rest_in_range, u14_close, biologically_plausible)), digits = 3)

write_csv(in_range, file.path(table_dir, "PE_Q2_feasibility_assessment.csv"))

# ---- Figure: does any alpha_R simultaneously satisfy all three? ----------
plot_df <- results |> select(alpha_R, q_recife_optimal, q_rest_optimal, u14_implied_prevalence) |>
  pivot_longer(-alpha_R, names_to = "quantity", values_to = "value")
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

p1 <- ggplot(results, aes(alpha_R)) +
  geom_ribbon(aes(ymin = 0.005, ymax = 0.03), fill = "#0072B2", alpha = 0.08) +
  geom_line(aes(y = q_recife_optimal), colour = "#0072B2", linewidth = 0.8) +
  geom_hline(yintercept = c(0.005, 0.03), linetype = "dashed", colour = "#0072B2", linewidth = 0.3) +
  labs(title = "Optimal q_Recife (best-fit to Recife's own cases) vs alpha_R",
       subtitle = "Shaded band = diagnostic grid bound [0.005, 0.03]", x = "alpha_R (shared transmission intercept)", y = "q_Recife (best-fit)") +
  theme_v4

p2 <- ggplot(results, aes(alpha_R)) +
  geom_ribbon(aes(ymin = 0.03, ymax = 0.12), fill = "#D55E00", alpha = 0.08) +
  geom_line(aes(y = q_rest_optimal), colour = "#D55E00", linewidth = 0.8) +
  geom_hline(yintercept = c(0.03, 0.12), linetype = "dashed", colour = "#D55E00", linewidth = 0.3) +
  labs(title = "Optimal q_rest (best-fit to other 4 regions' cases) vs alpha_R",
       subtitle = "Shaded band = diagnostic grid bound [0.03, 0.12]", x = "alpha_R", y = "q_rest (best-fit)") +
  theme_v4

p3 <- ggplot(results, aes(alpha_R, u14_implied_prevalence)) +
  geom_ribbon(aes(ymin = 0.340, ymax = 0.404), fill = "#009E73", alpha = 0.15) +
  geom_line(colour = "#009E73", linewidth = 0.8) +
  geom_hline(yintercept = 0.372, linetype = "dashed", colour = "#009E73") +
  labs(title = "Implied Recife U14-window prevalence vs alpha_R (independent of q by construction)",
       subtitle = "Shaded band = reported U14 95% CI [0.340, 0.404]", x = "alpha_R", y = "Implied Recife prevalence") +
  theme_v4

library(patchwork)
fig <- p1 / p2 / p3
ggsave(file.path(figure_dir, "PE_Q2_feasibility_diagnostic.png"), fig, width = 200, height = 240, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "PE_Q2_feasibility_diagnostic.png"))

message("\n=== VERDICT ===")
if (nrow(feasible) > 0) {
  message("FEASIBLE region(s) found -- a common alpha_R exists where q_Recife, q_rest, and U14 are jointly plausible:")
  print(as.data.frame(feasible), digits = 3)
} else {
  message("NO alpha_R in the grid simultaneously satisfies q_Recife in [0.005,0.03], q_rest in [0.03,0.12], AND U14 within 5pp of observed.")
  best_u14_row <- results[which.min(abs(results$u14_implied_prevalence - 0.372)), ]
  message("\nClosest U14 fit in the grid:")
  print(as.data.frame(best_u14_row), digits = 4)
}
