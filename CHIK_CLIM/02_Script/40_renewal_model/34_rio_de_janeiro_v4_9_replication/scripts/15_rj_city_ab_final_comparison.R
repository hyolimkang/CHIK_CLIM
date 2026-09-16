# RJ CITY -- final Model A vs Model B (ad099 rescue) comparison. NO new
# model is fit here -- this script only reads the two already-saved
# posteriors (Model A case-only, Model B city+U10 ad099 rescue) and
# computes the five comparison metrics the user asked for:
#   1. divergences == 0 for both
#   2. how much cor(logit_q, alpha_R) shrinks from A to B
#   3. how well the U10 posterior predictive covers the observed 18.0%
#   4. whether city case log-likelihood / PPC is preserved between A and B
#   5. how much the 2025 S/N (susceptible fraction) shifts from A to B

required_packages <- c("rstan", "dplyr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication/city")

a <- readRDS(file.path(base_dir, "outputs/city_caseonly/rj_city_global_q_case_only.rds"))
b <- readRDS(file.path(base_dir, "outputs/city_u10_ad099/rj_city_u10_ad099.rds"))

dates_a <- as.Date(a$weekly_data$week_start)
dates_b <- as.Date(b$weekly_data$week_start)
stopifnot(identical(dates_a, dates_b))
obs_cases <- a$weekly_data$cases

draws_a <- rstan::extract(a$fit, pars = c("logit_q", "alpha_R", "phi_obs", "expected_reported_cases", "C_pred", "S_prop"), permuted = TRUE)
draws_b <- rstan::extract(b$fit, pars = c("logit_q", "alpha_R", "phi_obs", "expected_reported_cases", "C_pred", "S_prop", "sero_pred", "p_city_window"), permuted = TRUE)

message("=== 1. Divergences ===")
message(sprintf("Model A: divergences = %d, hmc_pass = %s", a$hmc$divergences, a$hmc$hmc_pass))
message(sprintf("Model B (ad099): divergences = %d, hmc_pass = %s", b$hmc$divergences, b$hmc$hmc_pass))

message("\n=== 2. cor(logit_q, alpha_R): A -> B ===")
cor_a <- cor(draws_a$logit_q, draws_a$alpha_R)
cor_b <- cor(draws_b$logit_q, draws_b$alpha_R)
message(sprintf("Model A: %.4f", cor_a))
message(sprintf("Model B: %.4f", cor_b))
message(sprintf("Reduction: %.4f (%.1f%% relative shrinkage in |cor|)", abs(cor_a) - abs(cor_b), 100 * (abs(cor_a) - abs(cor_b)) / abs(cor_a)))

message("\n=== 3. U10 posterior predictive coverage of observed 18.0% ===")
sero_pred <- draws_b$sero_pred[, 1]
p_city_window <- draws_b$p_city_window[, 1]
obs_prev <- b$u10$weighted_prev; obs_lo <- b$u10$ci_lo; obs_hi <- b$u10$ci_hi
message(sprintf("Observed U10: %.3f, 95%% CI=[%.3f,%.3f]", obs_prev, obs_lo, obs_hi))
message(sprintf("Posterior p_city_window (model-implied true prevalence): median=%.4f, 95%% CrI=[%.4f,%.4f]",
                 median(p_city_window), quantile(p_city_window, .025), quantile(p_city_window, .975)))
message(sprintf("Posterior predictive sero_pred (incl. survey sampling noise): median=%.4f, 95%% PI=[%.4f,%.4f]",
                 median(sero_pred), quantile(sero_pred, .025), quantile(sero_pred, .975)))
prob_pp_covers_obs <- mean(sero_pred >= obs_lo & sero_pred <= obs_hi)
prob_pp_ge_obs_point <- mean(sero_pred >= obs_prev)
message(sprintf("P(sero_pred within observed 95%% CI [%.3f,%.3f]) = %.4f", obs_lo, obs_hi, prob_pp_covers_obs))
message(sprintf("P(sero_pred >= observed point estimate %.3f) = %.4f (two-sided extremity ~ %.4f)",
                 obs_prev, prob_pp_ge_obs_point, 2 * min(prob_pp_ge_obs_point, 1 - prob_pp_ge_obs_point)))

message("\n=== 4. City case log-likelihood / PPC: A vs B ===")
nb2_loglik_draw <- function(y, mu, phi) sum(dnbinom(y, mu = mu, size = phi, log = TRUE))
n_draws_a <- nrow(draws_a$expected_reported_cases)
n_draws_b <- nrow(draws_b$expected_reported_cases)
loglik_a <- vapply(seq_len(n_draws_a), function(i) nb2_loglik_draw(obs_cases, draws_a$expected_reported_cases[i, ], draws_a$phi_obs[i]), numeric(1))
loglik_b <- vapply(seq_len(n_draws_b), function(i) nb2_loglik_draw(obs_cases, draws_b$expected_reported_cases[i, ], draws_b$phi_obs[i]), numeric(1))
message(sprintf("Model A: total case log-lik median = %.1f (95%% CrI [%.1f,%.1f])", median(loglik_a), quantile(loglik_a,.025), quantile(loglik_a,.975)))
message(sprintf("Model B: total case log-lik median = %.1f (95%% CrI [%.1f,%.1f])", median(loglik_b), quantile(loglik_b,.025), quantile(loglik_b,.975)))
message(sprintf("Delta (B - A) median log-lik = %.1f", median(loglik_b) - median(loglik_a)))

ppc_lo_a <- apply(draws_a$C_pred, 2, quantile, .025); ppc_hi_a <- apply(draws_a$C_pred, 2, quantile, .975)
ppc_lo_b <- apply(draws_b$C_pred, 2, quantile, .025); ppc_hi_b <- apply(draws_b$C_pred, 2, quantile, .975)
coverage_a <- mean(obs_cases >= ppc_lo_a & obs_cases <= ppc_hi_a)
coverage_b <- mean(obs_cases >= ppc_lo_b & obs_cases <= ppc_hi_b)
rmse_a <- sqrt(mean((obs_cases - apply(draws_a$C_pred, 2, median))^2))
rmse_b <- sqrt(mean((obs_cases - apply(draws_b$C_pred, 2, median))^2))
message(sprintf("Model A: 95%% PPC weekly coverage = %.1f%%, RMSE = %.1f", 100 * coverage_a, rmse_a))
message(sprintf("Model B: 95%% PPC weekly coverage = %.1f%%, RMSE = %.1f", 100 * coverage_b, rmse_b))

message("\n=== 5. 2025 S/N (susceptible fraction) change: A -> B ===")
S_final_a <- draws_a$S_prop[, ncol(draws_a$S_prop)]
S_final_b <- draws_b$S_prop[, ncol(draws_b$S_prop)]
message(sprintf("Model A: S/N (final week, ~2025) median = %.4f, 95%% CrI=[%.4f,%.4f]", median(S_final_a), quantile(S_final_a,.025), quantile(S_final_a,.975)))
message(sprintf("Model B: S/N (final week, ~2025) median = %.4f, 95%% CrI=[%.4f,%.4f]", median(S_final_b), quantile(S_final_b,.025), quantile(S_final_b,.975)))
message(sprintf("Absolute shift (B - A) = %.4f (%.1f percentage points)", median(S_final_b) - median(S_final_a), 100 * (median(S_final_b) - median(S_final_a))))

comparison_tbl <- tibble(
  metric = c("divergences", "cor_logitq_alphaR", "case_loglik_median", "ppc_95_coverage_pct", "case_rmse", "S_final_2025_median"),
  model_A_caseonly = c(a$hmc$divergences, cor_a, median(loglik_a), 100 * coverage_a, rmse_a, median(S_final_a)),
  model_B_u10_ad099 = c(b$hmc$divergences, cor_b, median(loglik_b), 100 * coverage_b, rmse_b, median(S_final_b))
)
print(comparison_tbl, digits = 4)
readr::write_csv(comparison_tbl, file.path(table_dir, "RJ_CITY_AB_final_comparison.csv"))
message("\n[saved] ", file.path(table_dir, "RJ_CITY_AB_final_comparison.csv"))
