# Pernambuco -- FINAL structural feasibility test: fixed, zero-sum regional
# transmission intercepts (delta_region), ONE common q, common temporal
# shape (year effects, seasonality, A_year) fixed at the Spatial-0
# representative draw. NO Stan/HMC in this script -- a deterministic
# multi-start profile/optimisation, exactly mirroring the Q2 feasibility
# diagnostic's method (script 40) but now optimising over region-specific
# TRANSMISSION intercepts instead of ascertainment.
#
# log R0[r,t] = alpha_global + delta_region[r] + common_temporal[t],
# sum(delta_region) = 0 (5 values from 4 free contrasts + Recife as the
# dependent residual, per instruction: no sigma_region, no hierarchical
# geometry -- this script only PROFILES fixed values; the orthonormal
# zero-sum Stan parameterisation is deferred to the final candidate model,
# built only if this feasibility test passes).

required_packages <- c("dplyr", "readr", "ggplot2", "tidyr", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2); library(tidyr); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
pe_root <- file.path(root, "03_Output/04_pernambuco_pipeline/model_fits/v4_9_replication")
table_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/tables/pernambuco_v4_9_replication/spatial_model")
figure_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/figures/pernambuco_v4_9_replication/spatial_Rfixed_feasibility")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

STRATUM_LEVELS <- c("1_Recife", "2_Metropolitana_remainder", "3_Agreste", "4_Sertao", "5_Vale_Sao_Francisco_Araripe")
RECIFE_INDEX <- 1L

built <- readRDS(file.path(table_dir, "PE_5strata_stan_data.rds"))
sd <- built$stan_data
rep_draw <- readRDS(file.path(table_dir, "PE_spatial0_representative_draw.rds"))
N <- sd$N; R <- sd$R; G <- sd$G

common_temporal <- rep_draw$year_effect[sd$year_id] +
  rep_draw$A_year[sd$year_id] * (rep_draw$beta_sin1 * sd$seasonal_sin + rep_draw$beta_cos1 * sd$seasonal_cos) +
  rep_draw$beta_sin2 * sd$seasonal_sin2 + rep_draw$beta_cos2 * sd$seasonal_cos2

run_recursion_regional <- function(alpha_global, delta_region) {
  log_R0 <- outer(common_temporal, delta_region, "+") + alpha_global # [N,R]
  R0_t <- exp(log_R0)
  S <- matrix(0, N, R); U <- matrix(0, N, R); X <- matrix(0, N, R)
  S[1, ] <- sd$N_start[1, ]; U[1, ] <- 0
  reconciliation <- sd$N_end - sd$N_start - sd$births + sd$deaths
  for (t in 1:N) {
    total <- S[t, ] + U[t, ]
    infectiousness <- rep(0, R)
    for (g in 1:G) if (t > g) infectiousness <- infectiousness + sd$w[g] * X[t - g, ]
    foi <- R0_t[t, ] * infectiousness / total + sd$imports_r[t, ] / total
    X[t, ] <- S[t, ] * (-expm1(-foi))
    if (t < N) {
      susc_after <- S[t, ] - X[t, ]
      imm_after <- U[t, ] + X[t, ]
      susc_frac <- susc_after / (susc_after + imm_after)
      S[t + 1, ] <- susc_after + sd$births[t, ] - sd$deaths[t, ] * susc_frac + reconciliation[t, ] * susc_frac
      U[t + 1, ] <- imm_after - sd$deaths[t, ] * (1 - susc_frac) + reconciliation[t, ] * (1 - susc_frac)
    }
  }
  list(X = X, S = S, U = U, immune_prop = U / (S + U), S_prop = S / (S + U), R0_t = R0_t)
}

sero_idx <- sd$sero_window_start_idx:(sd$sero_window_start_idx + sd$sero_window_n_weeks - 1)
observed_C <- sd$C
kappa_sero <- sd$kappa_sero
dbetabinom_ll <- function(y, n, alpha, beta) lchoose(n, y) + lbeta(y + alpha, n - y + beta) - lbeta(alpha, beta)

poisson_deviance <- function(obs, mu) {
  mu <- pmax(mu, 1e-9)
  term <- ifelse(obs == 0, 0, obs * log(obs / mu))
  2 * sum(term - (obs - mu))
}

# ---- Objective: par = c(alpha_global, d2, d3, d4, d5); d1 (Recife) = -(d2+d3+d4+d5) ----
objective <- function(par, sero_weight = 8) {
  alpha_global <- par[1]
  d_rest <- par[2:5]
  delta_region <- c(-sum(d_rest), d_rest)
  rec <- run_recursion_regional(alpha_global, delta_region)
  q_hat <- sum(observed_C) / sum(rec$X)
  mu <- q_hat * rec$X
  dev_total <- poisson_deviance(observed_C, mu)
  p_recife <- mean(rec$immune_prop[sero_idx, RECIFE_INDEX])
  p_recife <- min(max(p_recife, 1e-6), 1 - 1e-6)
  sero_ll <- dbetabinom_ll(sd$sero_n_positive, sd$sero_n_tested, p_recife * kappa_sero, (1 - p_recife) * kappa_sero)
  dev_total / 100 + sero_weight * (-sero_ll)
}

message("=== Multi-start optimisation over (alpha_global, 4 free zero-sum regional contrasts) ===")
set.seed(20260915L)
n_starts <- 12L
starts <- lapply(seq_len(n_starts), function(i) {
  c(alpha_global = rep_draw$alpha_R_fitted + runif(1, -0.1, 0.3),
    d_rest = runif(4, -0.5, 0.5))
})
opt_results <- lapply(starts, function(s0) {
  tryCatch(optim(s0, objective, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-9)),
           error = function(e) list(par = s0, value = Inf, convergence = 99))
})
values <- sapply(opt_results, function(o) o$value)
best <- opt_results[[which.min(values)]]
message("Best objective value: ", round(best$value, 4), " (across ", n_starts, " starts, range ",
        round(min(values), 3), " to ", round(max(values), 3), ")")

alpha_global_best <- best$par[1]
delta_region_best <- c(-sum(best$par[2:5]), best$par[2:5])
names(delta_region_best) <- STRATUM_LEVELS
message("\nBest-fit alpha_global: ", round(alpha_global_best, 4))
message("Best-fit delta_region (log scale) / exp(delta_region) (multiplicative):")
print(round(rbind(delta_region = delta_region_best, multiplier = exp(delta_region_best)), 4))

rec_best <- run_recursion_regional(alpha_global_best, delta_region_best)
q_hat_best <- sum(observed_C) / sum(rec_best$X)
p_recife_best <- mean(rec_best$immune_prop[sero_idx, RECIFE_INDEX])
message(sprintf("\nCommon q (best-fit): %.4f", q_hat_best))
message(sprintf("Implied Recife U14 prevalence: %.4f (observed 0.372, 95%% CI [0.340,0.404])", p_recife_best))
message(sprintf("Max R0 (any region/week): %.3f", max(rec_best$R0_t)))
message(sprintf("Min S/N (any region/week): %.4f", min(rec_best$S_prop)))

# ---- Per-region case fit quality ----
regional_fit <- bind_rows(lapply(seq_len(R), function(r) {
  obs <- observed_C[, r]; mu <- q_hat_best * rec_best$X[, r]
  tibble(stratum = STRATUM_LEVELS[r], observed_total = sum(obs), predicted_total = sum(mu),
         ratio_pred_over_obs = sum(mu) / sum(obs),
         correlation = suppressWarnings(cor(obs, mu)),
         poisson_deviance = poisson_deviance(obs, mu),
         delta_region = delta_region_best[r], multiplier = exp(delta_region_best[r]),
         min_S_prop = min(rec_best$S_prop[, r]), cumulative_infected_2025 = rec_best$immune_prop[N, r])
}))
message("\n=== Per-region case fit at best-fit solution ===")
print(as.data.frame(regional_fit), digits = 3)
write_csv(regional_fit, file.path(table_dir, "PE_regional_R_profile.csv"))
message("[saved] ", file.path(table_dir, "PE_regional_R_profile.csv"))

# ---- Section 6: does the regional RANKING hold across epidemic years, or ----
#      does it require region-by-time heterogeneity? Check by comparing
#      each region's dominant-wave share (already computed empirically in
#      the spatial-turnover diagnostic) against what a FIXED multiplier
#      would predict (a fixed multiplier predicts the SAME region should
#      dominate every wave, scaled by common seasonality only).
wave_contrib <- read_csv(file.path(table_dir, "PE_5strata_wave_contributions.csv"), show_col_types = FALSE)
dominant_by_wave <- wave_contrib |> group_by(wave_id) |> slice_max(share_of_wave, n = 1) |> ungroup() |>
  select(wave_id, dominant_stratum = final_model_stratum, share_of_wave)
message("\n=== Empirical dominant stratum per wave (from the pre-model audit) ===")
print(as.data.frame(dominant_by_wave))
message("\nBest-fit static multiplier ranking (highest to lowest): ",
        paste(names(sort(exp(delta_region_best), decreasing = TRUE)), collapse = " > "))
n_distinct_dominant <- n_distinct(dominant_by_wave$dominant_stratum)
message(sprintf("Number of DISTINCT strata that dominate at least one wave: %d of %d strata", n_distinct_dominant, R))
if (n_distinct_dominant >= 3) {
  message("*** WARNING: the dominant stratum switches across >=3 different regions across waves. ***")
  message("*** A single FIXED regional ranking cannot represent 'who dominates when' -- this is evidence for region-BY-TIME heterogeneity, not fixed offsets. ***")
}

saveRDS(list(alpha_global_best = alpha_global_best, delta_region_best = delta_region_best,
             q_hat_best = q_hat_best, p_recife_best = p_recife_best, regional_fit = regional_fit,
             opt_results = opt_results, rec_best = rec_best),
        file.path(table_dir, "PE_regional_R_feasibility_result.rds"))
message("\n[saved] PE_regional_R_feasibility_result.rds")
