# Pernambuco Spatial Model 1 -- Sections 13-14: primary HMC gate and
# primary identification diagnostics.

required_packages <- c("rstan", "dplyr", "readr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
pe_root <- file.path(root, "03_Output/model_fits/pernambuco/v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_model")

b <- readRDS(file.path(pe_root, "outputs/spatial5/pe_5strata_pilot.rds"))
fit <- b$fit
stratum_levels <- b$stratum_levels
R <- b$stan_data$R

message("=== Section 13: primary HMC gate ===")
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
hmc_by_chain <- bind_rows(lapply(seq_along(sp), function(i) {
  tibble(chain = i, divergent = sum(sp[[i]][, "divergent__"]),
         max_treedepth_hits = sum(sp[[i]][, "treedepth__"] >= b$config$max_treedepth),
         mean_stepsize = mean(sp[[i]][, "stepsize__"]), mean_leapfrog = mean(sp[[i]][, "n_leapfrog__"]),
         mean_accept = mean(sp[[i]][, "accept_stat__"]))
}))
hmc_by_chain$bfmi <- rstan::get_bfmi(fit)
print(as.data.frame(hmc_by_chain), digits = 4)

# minimum tail ESS via monitor()
mon <- rstan::monitor(fit, print = FALSE)
min_tail_ess <- min(mon[, "Tail_ESS"], na.rm = TRUE)
message(sprintf("\nOverall: divergences=%d, max_treedepth_hits=%d, max Rhat=%.4f, min bulk ESS=%.1f, min tail ESS=%.1f, hmc_pass=%s",
                 b$hmc$divergences, b$hmc$max_treedepth_hits, b$hmc$maximum_rhat, b$hmc$minimum_bulk_ess,
                 min_tail_ess, b$hmc$hmc_pass))

write_csv(hmc_by_chain, file.path(table_dir, "PE_5strata_HMC.csv"))
message("[saved] ", file.path(table_dir, "PE_5strata_HMC.csv"))

hmc_pass <- b$hmc$hmc_pass
if (!hmc_pass) {
  message("\n*** HMC GATE FAIL. Per the pre-registered stop rule, substantive susceptibility estimates below are NOT to be interpreted as reliable -- reported for diagnosis only. ***")
}

message("\n=== Section 14: primary identification diagnostics ===")
draws <- rstan::extract(fit, permuted = TRUE)
q_draws <- plogis(draws$logit_q)
q_summary <- tibble(median = median(q_draws), q25 = quantile(q_draws, .25), q75 = quantile(q_draws, .75),
                     lo95 = quantile(q_draws, .025), hi95 = quantile(q_draws, .975))
message("q_PE posterior: median=", round(q_summary$median, 4), " 50% CrI=[", round(q_summary$q25, 4), ",", round(q_summary$q75, 4),
        "] 95% CrI=[", round(q_summary$lo95, 4), ",", round(q_summary$hi95, 4), "]")

q_prior_draws <- plogis(draws$logit_q_prior_draw)
message("q_PE prior (same draws' prior replicate): median=", round(median(q_prior_draws), 4),
        " 95% CrI=[", round(quantile(q_prior_draws, .025), 4), ",", round(quantile(q_prior_draws, .975), 4), "]")

cor_q_alpha <- cor(draws$logit_q, draws$alpha_global)
message(sprintf("cor(logit_q, alpha_global) = %.4f", cor_q_alpha))

immune_2025 <- draws$U_PE_prop[, dim(draws$U_PE_prop)[2]]
cor_q_immune <- cor(q_draws, immune_2025)
message(sprintf("cor(q_PE, statewide immune_2025) = %.4f", cor_q_immune))
message(sprintf("Statewide immune_2025: median=%.4f 95%% CrI=[%.4f,%.4f]", median(immune_2025), quantile(immune_2025,.025), quantile(immune_2025,.975)))

write_csv(q_summary, file.path(table_dir, "PE_5strata_q_summary.csv"))
message("[saved] ", file.path(table_dir, "PE_5strata_q_summary.csv"))

# ---- Regional immune fraction at key checkpoints ---------------------------
weeks <- readRDS(file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic/pe_municipality_week_panel.rds")) |>
  pull(week_start) |> unique() |> sort()
checkpoint_dates <- as.Date(c("2016-12-31", "2018-12-31", "2021-12-31", "2022-12-31", "2025-12-21"))
checkpoint_idx <- sapply(checkpoint_dates, function(d) which.min(abs(weeks - d)))

immune_prop_arr <- draws$immune_prop # [iter, N, R]
S_prop_arr <- draws$S_prop
R0_arr <- draws$R0_t

susceptibility_tbl <- bind_rows(lapply(seq_along(checkpoint_dates), function(i) {
  t_idx <- checkpoint_idx[i]
  bind_rows(lapply(seq_len(R), function(r) {
    tibble(checkpoint = checkpoint_dates[i], stratum = stratum_levels[r],
           immune_median = median(immune_prop_arr[, t_idx, r]),
           immune_lo95 = quantile(immune_prop_arr[, t_idx, r], .025),
           immune_hi95 = quantile(immune_prop_arr[, t_idx, r], .975),
           S_prop_median = median(S_prop_arr[, t_idx, r]))
  }))
}))
message("\n=== Regional immune fraction by checkpoint ===")
print(as.data.frame(susceptibility_tbl), digits = 3)
write_csv(susceptibility_tbl, file.path(table_dir, "PE_5strata_susceptibility.csv"))
message("[saved] ", file.path(table_dir, "PE_5strata_susceptibility.csv"))

# ---- Statewide S_PE(t), U_PE(t) summary at same checkpoints ----------------
statewide_tbl <- tibble(
  checkpoint = checkpoint_dates,
  S_PE_median = sapply(checkpoint_idx, function(t) median(draws$S_PE_prop[, t])),
  U_PE_median = sapply(checkpoint_idx, function(t) median(draws$U_PE_prop[, t]))
)
message("\n=== Statewide S_PE(t) / U_PE(t) at checkpoints ===")
print(as.data.frame(statewide_tbl), digits = 3)

# ---- Does the later large epidemic (2021/2022) occur in strata retaining ----
#      susceptibles, or does it require extreme R0 in already-depleted strata?
onset_2021 <- weeks[which.min(abs(weeks - as.Date("2021-02-07")))]
onset_2022 <- weeks[which.min(abs(weeks - as.Date("2022-01-16")))]
idx_2021 <- which(weeks == onset_2021); idx_2022 <- which(weeks == onset_2022)

message("\n=== S_prop and max R0 at 2021 / 2022 epidemic onset, by stratum ===")
onset_check <- bind_rows(lapply(seq_len(R), function(r) {
  tibble(stratum = stratum_levels[r],
         S_prop_at_2021_onset = median(S_prop_arr[, idx_2021, r]),
         R0_at_2021_onset = median(R0_arr[, idx_2021, r]),
         S_prop_at_2022_onset = median(S_prop_arr[, idx_2022, r]),
         R0_at_2022_onset = median(R0_arr[, idx_2022, r]))
}))
print(as.data.frame(onset_check), digits = 3)
write_csv(onset_check, file.path(table_dir, "PE_5strata_onset_S_R0_by_stratum.csv"))
message("[saved] ", file.path(table_dir, "PE_5strata_onset_S_R0_by_stratum.csv"))
