# Identification diagnostics for the Bahia global-q + multisite-serology
# pilot (Section "IDENTIFICATION DIAGNOSTICS" / "SEROLOGY INFORMATION FLOW"
# of the design spec). Reports q prior vs posterior, posterior correlations
# between q and key transmission/immunity quantities, and the per-survey
# eta_geo/p_state/p_site decomposition table. Diagnostic only -- does not
# alter the model.

required_packages <- c("rstan", "dplyr", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(readr) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/03_bahia_pipeline/model_fits/v4_9_replication")

b <- readRDS(file.path(base_dir, "outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds")) # "on_full2" matches the actual saved fit tag; the script's original "pilot" default never corresponded to any saved fit
fit <- b$fit
audit <- b$sero_audit
w <- b$weekly_data
dates <- as.Date(w$week_start)

message("=== q: prior vs posterior ===")
q_post <- as.vector(rstan::extract(fit, "q")$q)
q_prior_draw <- as.vector(rstan::extract(fit, "q_prior_draw")$q_prior_draw)
message(sprintf("Prior : median=%.4f | 50%% CrI=[%.4f, %.4f] | 95%% CrI=[%.4f, %.4f]",
                 median(q_prior_draw), quantile(q_prior_draw, .25), quantile(q_prior_draw, .75),
                 quantile(q_prior_draw, .025), quantile(q_prior_draw, .975)))
message(sprintf("Posterior: median=%.4f | 50%% CrI=[%.4f, %.4f] | 95%% CrI=[%.4f, %.4f]",
                 median(q_post), quantile(q_post, .25), quantile(q_post, .75),
                 quantile(q_post, .025), quantile(q_post, .975)))
prior_sd <- sd(q_prior_draw); post_sd <- sd(q_post)
contraction <- 1 - (post_sd / prior_sd)
message(sprintf("SD contraction: prior sd=%.4f -> posterior sd=%.4f (%.1f%% contraction)", prior_sd, post_sd, 100 * contraction))

message("\n=== Posterior correlations: q vs key quantities ===")
idx_2016 <- max(which(format(dates, "%Y") == "2016"))
idx_2018 <- max(which(format(dates, "%Y") == "2018"))
idx_2025 <- nrow(w)
immune_prop <- rstan::extract(fit, "immune_prop")$immune_prop
alpha_R <- as.vector(rstan::extract(fit, "alpha_R")$alpha_R)
A_year <- rstan::extract(fit, "A_year")$A_year
phi_obs <- as.vector(rstan::extract(fit, "phi_obs")$phi_obs)

cor_table <- tibble(
  quantity = c("immune_2016", "immune_2018", "immune_2025", "alpha_R", "phi_obs",
               paste0("A_year[", seq_len(ncol(A_year)), "] (", b$config$years, ")")),
  cor_with_q = c(
    cor(q_post, immune_prop[, idx_2016]), cor(q_post, immune_prop[, idx_2018]), cor(q_post, immune_prop[, idx_2025]),
    cor(q_post, alpha_R), cor(q_post, phi_obs),
    sapply(seq_len(ncol(A_year)), function(k) cor(q_post, A_year[, k]))
  )
)
print(as.data.frame(cor_table))

message("\n=== Serology information-flow decomposition (per survey) ===")
eta_geo <- rstan::extract(fit, "eta_geo")$eta_geo
p_state_window <- rstan::extract(fit, "p_state_window")$p_state_window
p_site_window <- rstan::extract(fit, "p_site_window")$p_site_window
decomp <- tibble(
  survey_id = audit$sero_id,
  observed_prevalence = audit$observed_prevalence,
  p_state_median = apply(p_state_window, 2, median),
  p_site_median = apply(p_site_window, 2, median),
  eta_geo_median = apply(eta_geo, 2, median),
  eta_geo_lo95 = apply(eta_geo, 2, quantile, .025),
  eta_geo_hi95 = apply(eta_geo, 2, quantile, .975)
)
print(as.data.frame(decomp))
message(sprintf("median |eta_geo| = %.3f | max |eta_geo| = %.3f", median(abs(decomp$eta_geo_median)), max(abs(decomp$eta_geo_median))))

message("\n=== Comparison vs fixed q=0.05 multisite pilot (same 6 surveys) ===")
fixed <- readRDS(file.path(base_dir, "outputs_multisite_serology/q0.05/renewal_bahia_multisite_fit_q0.05.rds"))
eta_geo_fixed <- rstan::extract(fixed$fit, "eta_geo")$eta_geo
fixed_decomp <- tibble(survey_id = fixed$sero_audit$sero_id, eta_geo_median_fixed_q0.05 = apply(eta_geo_fixed, 2, median))
comparison <- left_join(decomp |> select(survey_id, eta_geo_median_global_q = eta_geo_median), fixed_decomp, by = "survey_id")
print(as.data.frame(comparison))

table_dir <- file.path(root, "03_Output/03_bahia_pipeline/tables/renewal_bahia_v4_9_global_q_multisite_serology")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(cor_table, file.path(table_dir, "global_q_posterior_correlations.csv"))
write_csv(decomp, file.path(table_dir, "global_q_serology_decomposition.csv"))
write_csv(comparison, file.path(table_dir, "global_q_vs_fixed_q005_eta_geo_comparison.csv"))
message("\n[saved] tables under: ", table_dir)
