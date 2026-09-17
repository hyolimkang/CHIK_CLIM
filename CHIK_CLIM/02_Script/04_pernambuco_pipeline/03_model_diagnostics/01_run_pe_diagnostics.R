# Pernambuco external replication -- Sections 8, 9, 10 diagnostics.
# Model A (case-only) vs Model B (U14 Recife serology). Identification
# questions, wave-level case PPC (existing objective PE wave census, not
# redefined), and U14 serology diagnostics.

required_packages <- c("rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/model_fits/v4_9_replication")
table_dir <- file.path(root, "03_Output/04_pernambuco_pipeline/tables/pernambuco_v4_9_replication")

A <- readRDS(file.path(base_dir, "outputs/modelA_caseonly/renewal_pe_global_q_fit_modelA_caseonly.rds"))
B <- readRDS(file.path(base_dir, "outputs/modelB_U14/renewal_pe_global_q_fit_modelB_U14.rds"))

# ---- Section 8: identification questions -----------------------------------
q_A <- as.vector(rstan::extract(A$fit, "q")$q)
q_B <- as.vector(rstan::extract(B$fit, "q")$q)
q_prior_B <- as.vector(rstan::extract(B$fit, "q_prior_draw")$q_prior_draw)

summarise_q <- function(x) c(median = median(x), lo50 = quantile(x, .25), hi50 = quantile(x, .75), lo95 = quantile(x, .025), hi95 = quantile(x, .975))

message("=== Section 8: q summaries ===")
message("Model A (case-only): ", paste(sprintf("%s=%.4f", names(summarise_q(q_A)), summarise_q(q_A)), collapse = " | "))
message("Model B (U14 serology): ", paste(sprintf("%s=%.4f", names(summarise_q(q_B)), summarise_q(q_B)), collapse = " | "))
message("q prior (B): ", paste(sprintf("%s=%.4f", names(summarise_q(q_prior_B)), summarise_q(q_prior_B)), collapse = " | "))

dates_A <- as.Date(A$weekly_data$week_start); dates_B <- as.Date(B$weekly_data$week_start)
idx_2018_A <- max(which(format(dates_A, "%Y") == "2018")); idx_2025_A <- nrow(A$weekly_data)
idx_2018_B <- max(which(format(dates_B, "%Y") == "2018")); idx_2025_B <- nrow(B$weekly_data)

immune_A <- rstan::extract(A$fit, "immune_prop")$immune_prop
immune_B <- rstan::extract(B$fit, "immune_prop")$immune_prop
alpha_R_A <- as.vector(rstan::extract(A$fit, "alpha_R")$alpha_R); alpha_R_B <- as.vector(rstan::extract(B$fit, "alpha_R")$alpha_R)
phi_obs_A <- as.vector(rstan::extract(A$fit, "phi_obs")$phi_obs); phi_obs_B <- as.vector(rstan::extract(B$fit, "phi_obs")$phi_obs)
A_year_A <- rstan::extract(A$fit, "A_year")$A_year; A_year_B <- rstan::extract(B$fit, "A_year")$A_year

cor_table <- tibble(
  model = c(rep("A_case_only", 4), rep("B_U14_serology", 4)),
  quantity = rep(c("immune_2018", "immune_2025", "alpha_R", "phi_obs"), 2),
  cor_with_q = c(
    cor(q_A, immune_A[, idx_2018_A]), cor(q_A, immune_A[, idx_2025_A]), cor(q_A, alpha_R_A), cor(q_A, phi_obs_A),
    cor(q_B, immune_B[, idx_2018_B]), cor(q_B, immune_B[, idx_2025_B]), cor(q_B, alpha_R_B), cor(q_B, phi_obs_B)
  )
)
message("\n=== Posterior correlations with q ===")
print(as.data.frame(cor_table), digits = 3)

immune_checkpoints <- c("2016", "2017", "2018", "2019", "2022", "2025")
immune_summary <- bind_rows(lapply(immune_checkpoints, function(y) {
  idx_A <- if (y == "2025") nrow(A$weekly_data) else max(which(format(dates_A, "%Y") == y))
  idx_B <- if (y == "2025") nrow(B$weekly_data) else max(which(format(dates_B, "%Y") == y))
  tibble(year = y,
         A_median = median(immune_A[, idx_A]), A_lo95 = quantile(immune_A[, idx_A], .025), A_hi95 = quantile(immune_A[, idx_A], .975),
         B_median = median(immune_B[, idx_B]), B_lo95 = quantile(immune_B[, idx_B], .025), B_hi95 = quantile(immune_B[, idx_B], .975))
}))
message("\n=== Immune fraction by year, A vs B ===")
print(as.data.frame(immune_summary), digits = 3)

write_csv(cor_table, file.path(table_dir, "pernambuco_q_correlations_A_vs_B.csv"))
write_csv(immune_summary, file.path(table_dir, "pernambuco_immune_fraction_A_vs_B.csv"))
message("\n[saved] pernambuco_q_correlations_A_vs_B.csv, pernambuco_immune_fraction_A_vs_B.csv")

# ---- Section 9: wave-level case PPC (existing objective PE census) --------
wave_census <- read_csv(file.path(root, "03_Output/07_national_pipeline/tables/national_wave_census_v2/brazil_chik_wave_census_v2.csv"), show_col_types = FALSE)
pe_waves <- wave_census |> dplyr::filter(state == "PE", major_epidemic_primary == TRUE) |>
  transmute(wave_id, onset_week = as.Date(onset_week, format = "%m/%d/%Y"),
            peak_week_obs = as.Date(peak_week, format = "%m/%d/%Y"),
            end_week = as.Date(end_week, format = "%m/%d/%Y")) |>
  arrange(onset_week)

wave_ppc_one_model <- function(bundle, model_label) {
  dates <- as.Date(bundle$weekly_data$week_start)
  C_pred <- rstan::extract(bundle$fit, "C_pred")$C_pred
  bind_rows(lapply(seq_len(nrow(pe_waves)), function(i) {
    w <- pe_waves[i, ]
    idx <- which(dates >= w$onset_week & dates <= w$end_week)
    obs_series <- bundle$weekly_data$cases[idx]
    pred_series <- C_pred[, idx, drop = FALSE]
    obs_total <- sum(obs_series); pred_total <- rowSums(pred_series)
    obs_peak <- max(obs_series); pred_peak <- apply(pred_series, 1, max)
    obs_peak_week <- dates[idx][which.max(obs_series)]
    pred_peak_week_idx <- apply(pred_series, 1, which.max)
    tibble(model = model_label, wave_id = w$wave_id, onset_week = w$onset_week, end_week = w$end_week,
           epidemic_width_weeks = length(idx),
           obs_total = obs_total, pred_total_median = median(pred_total), pred_obs_ratio_total = median(pred_total) / obs_total,
           obs_peak = obs_peak, pred_peak_median = median(pred_peak), pred_obs_ratio_peak = median(pred_peak) / obs_peak,
           obs_peak_week = obs_peak_week, pred_peak_week_median = dates[idx][round(median(pred_peak_week_idx))])
  }))
}

wave_ppc <- bind_rows(wave_ppc_one_model(A, "A_case_only"), wave_ppc_one_model(B, "B_U14_serology"))
message("\n=== Section 9: wave-level case PPC ===")
print(as.data.frame(wave_ppc |> select(model, wave_id, obs_total, pred_total_median, pred_obs_ratio_total, obs_peak, pred_peak_median, pred_obs_ratio_peak)), digits = 3)
write_csv(wave_ppc, file.path(table_dir, "pernambuco_wave_ppc.csv"))
message("[saved] ", file.path(table_dir, "pernambuco_wave_ppc.csv"))

# ---- Section 10: U14 serology diagnostics (Model B only) -------------------
p_state_window <- rstan::extract(B$fit, "p_state_window")$p_state_window[, 1]
p_site_window <- rstan::extract(B$fit, "p_site_window")$p_site_window[, 1]
eta_geo <- rstan::extract(B$fit, "eta_geo")$eta_geo[, 1]
sero_pred <- rstan::extract(B$fit, "sero_pred")$sero_pred[, 1]

u14_summary <- tibble(
  observed_prevalence = 770 / 2070,
  observed_95CI_lo = 0.340, observed_95CI_hi = 0.404,
  pe_state_immune_window_median = median(p_state_window), pe_state_immune_window_lo95 = quantile(p_state_window, .025), pe_state_immune_window_hi95 = quantile(p_state_window, .975),
  geo_adjusted_recife_median = median(p_site_window), geo_adjusted_recife_lo95 = quantile(p_site_window, .025), geo_adjusted_recife_hi95 = quantile(p_site_window, .975),
  eta_geo_U14_median = median(eta_geo), eta_geo_U14_lo95 = quantile(eta_geo, .025), eta_geo_U14_hi95 = quantile(eta_geo, .975),
  sero_pred_ppc_median = median(sero_pred) / 2070, sero_pred_ppc_lo95 = quantile(sero_pred, .025) / 2070, sero_pred_ppc_hi95 = quantile(sero_pred, .975) / 2070
)
message("\n=== Section 10: U14 serology diagnostics (Model B) ===")
print(as.data.frame(u14_summary), digits = 4)
write_csv(u14_summary, file.path(table_dir, "pernambuco_U14_serology_diagnostics.csv"))
message("[saved] ", file.path(table_dir, "pernambuco_U14_serology_diagnostics.csv"))

message("\nNOTE: eta_geo_U14 is interpreted strictly as a local-to-state serology TRANSPORT discrepancy (Recife vs Pernambuco state), NOT specifically as reporting bias, mosquito exposure, or ascertainment heterogeneity (Section 10 instruction).")

# ---- HMC diagnostics table --------------------------------------------------
hmc_table <- bind_rows(
  tibble(model = "A_case_only", divergences = A$hmc$divergences, max_treedepth_hits = A$hmc$max_treedepth_hits,
         max_rhat = A$hmc$maximum_rhat, min_bulk_ess = A$hmc$minimum_bulk_ess, bfmi_min = min(A$hmc$bfmi), hmc_pass = A$hmc$hmc_pass),
  tibble(model = "B_U14_serology", divergences = B$hmc$divergences, max_treedepth_hits = B$hmc$max_treedepth_hits,
         max_rhat = B$hmc$maximum_rhat, min_bulk_ess = B$hmc$minimum_bulk_ess, bfmi_min = min(B$hmc$bfmi), hmc_pass = B$hmc$hmc_pass)
)
message("\n=== HMC diagnostics ===")
print(as.data.frame(hmc_table), digits = 4)
write_csv(hmc_table, file.path(table_dir, "pernambuco_hmc_diagnostics.csv"))
message("[saved] ", file.path(table_dir, "pernambuco_hmc_diagnostics.csv"))
