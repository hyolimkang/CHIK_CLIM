# Rio de Janeiro CITY auxiliary -- Model B diagnostics: HMC gate, U10 fit
# (Section 11 scientific gate), Model A vs Model B comparison (Section 10),
# and the city-vs-state comparison (Section 12).

required_packages <- c("rstan", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication/city")
figure_dir <- file.path(root, "03_Output/figures/rio_de_janeiro_v4_9_replication/city")

bA <- readRDS(file.path(base_dir, "outputs/city_caseonly/rj_city_global_q_case_only.rds"))
bB <- readRDS(file.path(base_dir, "outputs/city_u10/rj_city_u10.rds"))
fitB <- bB$fit
weekly <- bB$weekly_data
dates <- as.Date(weekly$week_start)

message("=== Model B (CITY + U10): HMC gate ===")
sp <- rstan::get_sampler_params(fitB, inc_warmup = FALSE)
hmc_by_chain <- bind_rows(lapply(seq_along(sp), function(i) {
  tibble(chain = i, divergent = sum(sp[[i]][, "divergent__"]),
         max_treedepth_hits = sum(sp[[i]][, "treedepth__"] >= bB$config$max_treedepth),
         mean_stepsize = mean(sp[[i]][, "stepsize__"]), mean_accept = mean(sp[[i]][, "accept_stat__"]))
}))
hmc_by_chain$bfmi <- rstan::get_bfmi(fitB)
print(as.data.frame(hmc_by_chain), digits = 4)
message(sprintf("\nOverall: divergences=%d, max_treedepth_hits=%d, max Rhat=%.4f, min bulk ESS=%.1f, hmc_pass=%s",
                 bB$hmc$divergences, bB$hmc$max_treedepth_hits, bB$hmc$maximum_rhat, bB$hmc$minimum_bulk_ess, bB$hmc$hmc_pass))
write_csv(hmc_by_chain, file.path(table_dir, "RJ_CITY_u10_HMC.csv"))

drawsB <- rstan::extract(fitB, pars = c("logit_q", "alpha_R", "S_prop", "immune_prop", "R0_t", "expected_reported_cases", "C_pred", "p_city_window", "sero_pred"), permuted = TRUE)
q_B <- plogis(drawsB$logit_q)
cor_qalpha_B <- cor(drawsB$logit_q, drawsB$alpha_R)
p_window_B <- drawsB$p_city_window[, 1]
immune_2025_B <- drawsB$immune_prop[, ncol(drawsB$immune_prop)]
min_S_B <- apply(drawsB$S_prop, 1, min); max_R0_B <- apply(drawsB$R0_t, 1, max)

message("\n=== Section 11: U10 scientific gate (Model B) ===")
message(sprintf("Posterior city prevalence, U10 window: median=%.4f, 50%% CrI=[%.4f,%.4f], 95%% CrI=[%.4f,%.4f]",
                 median(p_window_B), quantile(p_window_B,.25), quantile(p_window_B,.75), quantile(p_window_B,.025), quantile(p_window_B,.975)))
message(sprintf("Observed U10 weighted prevalence: %.3f, 95%% CI=[%.3f,%.3f]", bB$u10$weighted_prev, bB$u10$ci_lo, bB$u10$ci_hi))
prob_in_ci <- mean(p_window_B >= bB$u10$ci_lo & p_window_B <= bB$u10$ci_hi)
message(sprintf("Posterior probability p_city_window in observed CI: %.4f", prob_in_ci))
pp_sero <- drawsB$sero_pred[, 1]
message(sprintf("Posterior predictive (sero_pred) median=%.4f, 95%% PI=[%.4f,%.4f]", median(pp_sero), quantile(pp_sero,.025), quantile(pp_sero,.975)))

# ---- Model A vs Model B comparison (Section 10) ----
drawsA <- rstan::extract(bA$fit, pars = c("logit_q", "alpha_R", "S_prop", "immune_prop", "R0_t"), permuted = TRUE)
q_A <- plogis(drawsA$logit_q)
cor_qalpha_A <- cor(drawsA$logit_q, drawsA$alpha_R)
datesA <- as.Date(bA$weekly_data$week_start)
u10idx_A <- which(datesA >= as.Date("2018-07-01") & datesA <= as.Date("2018-10-31"))
p_window_A <- rowMeans(drawsA$immune_prop[, u10idx_A, drop = FALSE])
immune_2025_A <- drawsA$immune_prop[, ncol(drawsA$immune_prop)]
min_S_A <- apply(drawsA$S_prop, 1, min); max_R0_A <- apply(drawsA$R0_t, 1, max)

comparison <- tibble(
  metric = c("q_city (median [95% CrI])", "cor(logit_q, alpha_R)", "U10-window immunity (median [95% CrI])",
             "immune_2025 (median)", "min S/N (median)", "max R0 (median)", "HMC divergences", "HMC pass"),
  city_cases_only = c(sprintf("%.4f [%.4f,%.4f]", median(q_A), quantile(q_A,.025), quantile(q_A,.975)),
                       sprintf("%.3f", cor_qalpha_A), sprintf("%.4f [%.4f,%.4f]", median(p_window_A), quantile(p_window_A,.025), quantile(p_window_A,.975)),
                       sprintf("%.4f", median(immune_2025_A)), sprintf("%.4f", median(min_S_A)), sprintf("%.2f", median(max_R0_A)),
                       as.character(bA$hmc$divergences), as.character(bA$hmc$hmc_pass)),
  city_cases_plus_u10 = c(sprintf("%.4f [%.4f,%.4f]", median(q_B), quantile(q_B,.025), quantile(q_B,.975)),
                           sprintf("%.3f", cor_qalpha_B), sprintf("%.4f [%.4f,%.4f]", median(p_window_B), quantile(p_window_B,.025), quantile(p_window_B,.975)),
                           sprintf("%.4f", median(immune_2025_B)), sprintf("%.4f", median(min_S_B)), sprintf("%.2f", median(max_R0_B)),
                           as.character(bB$hmc$divergences), as.character(bB$hmc$hmc_pass))
)
message("\n=== Section 10: Model A vs Model B comparison ===")
print(as.data.frame(comparison))
write_csv(comparison, file.path(table_dir, "RJ_CITY_model_comparison.csv"))

u10_summary <- tibble(quantity = c("posterior_median", "posterior_lo50", "posterior_hi50", "posterior_lo95", "posterior_hi95",
                                    "observed_weighted", "observed_ci_lo", "observed_ci_hi", "prob_in_observed_ci"),
                       value = c(median(p_window_B), quantile(p_window_B,.25), quantile(p_window_B,.75), quantile(p_window_B,.025), quantile(p_window_B,.975),
                                 bB$u10$weighted_prev, bB$u10$ci_lo, bB$u10$ci_hi, prob_in_ci))
write_csv(u10_summary, file.path(table_dir, "RJ_CITY_U10_summary.csv"))

# ---- Section 12: city vs state q comparison ----
state_b <- readRDS(file.path(base_dir, "outputs/caseonly/rj_global_q_case_only.rds"))
draws_state <- rstan::extract(state_b$fit, pars = "logit_q")
q_state <- plogis(draws_state$logit_q)
message(sprintf("\n=== Section 12: city vs state q ===\nq_state (case-only) median=%.4f [%.4f,%.4f]\nq_city (case-only) median=%.4f [%.4f,%.4f]\nq_city (+U10) median=%.4f [%.4f,%.4f]",
                 median(q_state), quantile(q_state,.025), quantile(q_state,.975),
                 median(q_A), quantile(q_A,.025), quantile(q_A,.975),
                 median(q_B), quantile(q_B,.025), quantile(q_B,.975)))

# ---- Figures ----
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
summarise_traj <- function(m) tibble(median = apply(m, 2, median), lo95 = apply(m, 2, quantile, .025), hi95 = apply(m, 2, quantile, .975))

pred_summary <- summarise_traj(drawsB$C_pred); pred_summary$week_start <- dates
exp_summary <- summarise_traj(drawsB$expected_reported_cases); exp_summary$week_start <- dates
obs_df <- tibble(week_start = dates, observed = weekly$cases)
p_cases <- ggplot() +
  geom_ribbon(data = pred_summary, aes(week_start, ymin = lo95, ymax = hi95), fill = "#56B4E9", alpha = 0.2) +
  geom_line(data = exp_summary, aes(week_start, median), colour = "#315A7D", linewidth = 0.5) +
  geom_point(data = obs_df, aes(week_start, observed), colour = "black", size = 0.4, alpha = 0.6) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "RJ CITY (Model B, case + U10): weekly case PPC", x = NULL, y = "Weekly reported cases") + theme_v4

U_summary <- summarise_traj(drawsB$immune_prop); U_summary$week_start <- dates
u10_pt <- tibble(week_start = as.Date("2018-07-01") + as.numeric(as.Date("2018-10-31") - as.Date("2018-07-01")) / 2,
                  observed = bB$u10$weighted_prev, lo = bB$u10$ci_lo, hi = bB$u10$ci_hi)
p_immunity <- ggplot(U_summary, aes(week_start, median)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#76558F", alpha = 0.18) + geom_line(colour = "#76558F", linewidth = 0.6) +
  geom_errorbar(data = u10_pt, aes(week_start, y = observed, ymin = lo, ymax = hi), width = 60, colour = "black", inherit.aes = FALSE) +
  geom_point(data = u10_pt, aes(week_start, observed), shape = 18, size = 2.5, colour = "black", inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(title = "RJ CITY immune fraction (Model B) with U10 DIRECTLY fit (no eta_geo -- geography matched)",
       x = NULL, y = "Immune fraction") + theme_v4

figB <- p_cases / p_immunity + patchwork::plot_annotation(
  title = "Rio de Janeiro CITY -- Model B (case + U10, logit-normal, geography matched)",
  subtitle = sprintf("%s | q_city median=%.4f | U10 posterior probability in observed CI=%.3f",
                      if (bB$hmc$hmc_pass) "HMC PASS" else sprintf("HMC FAIL (div=%d)", bB$hmc$divergences), median(q_B), prob_in_ci))
ggsave(file.path(figure_dir, "RJ_CITY_U10_trajectories.png"), figB, width = 200, height = 160, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "RJ_CITY_U10_trajectories.png"))

# U10 PPC figure
p_ppc <- ggplot(tibble(pp = pp_sero), aes(pp)) +
  geom_histogram(bins = 60, fill = "#315A7D", alpha = 0.7) +
  geom_vline(xintercept = bB$u10$weighted_prev, colour = "#D55E00", linewidth = 1) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "RJ CITY Model B: U10 posterior predictive check",
       subtitle = sprintf("Observed weighted prevalence = %.3f (orange line)", bB$u10$weighted_prev),
       x = "Posterior predictive weighted prevalence", y = "Count") + theme_v4
ggsave(file.path(figure_dir, "RJ_CITY_U10_PPC.png"), p_ppc, width = 180, height = 100, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "RJ_CITY_U10_PPC.png"))

# City vs state q figure
q_compare_df <- bind_rows(tibble(model = "State (case-only)", q = q_state), tibble(model = "City (case-only)", q = q_A), tibble(model = "City (+U10)", q = q_B))
p_qcompare <- ggplot(q_compare_df, aes(q, fill = model)) + geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("State (case-only)" = "#315A7D", "City (case-only)" = "#E69F00", "City (+U10)" = "#009E73")) +
  scale_x_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "q posterior: RJ state vs RJ city (case-only vs +U10)", x = "q", y = "Density") + theme_v4 + theme(legend.position = "top", legend.title = element_blank())
ggsave(file.path(figure_dir, "RJ_CITY_vs_STATE_q.png"), p_qcompare, width = 180, height = 110, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "RJ_CITY_vs_STATE_q.png"))

# City vs state susceptibility figure
S_state <- summarise_traj(drawsA$S_prop) # reuse city Model A's own dates since identical window; also add state
S_state_true <- summarise_traj(rstan::extract(state_b$fit, "S_prop")$S_prop); S_state_true$week_start <- as.Date(state_b$weekly_data$week_start)
S_cityA <- summarise_traj(drawsA$S_prop); S_cityA$week_start <- datesA
S_cityB <- summarise_traj(drawsB$S_prop); S_cityB$week_start <- dates
susc_df <- bind_rows(S_state_true |> mutate(model = "State (case-only)"), S_cityA |> mutate(model = "City (case-only)"), S_cityB |> mutate(model = "City (+U10)"))
p_susc <- ggplot(susc_df, aes(week_start, median, colour = model, fill = model)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), alpha = 0.12, colour = NA) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("State (case-only)" = "#315A7D", "City (case-only)" = "#E69F00", "City (+U10)" = "#009E73")) +
  scale_fill_manual(values = c("State (case-only)" = "#315A7D", "City (case-only)" = "#E69F00", "City (+U10)" = "#009E73"), guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "S/N: RJ state vs RJ city (case-only vs +U10)", x = NULL, y = "S / N") + theme_v4 + theme(legend.position = "top", legend.title = element_blank())
ggsave(file.path(figure_dir, "RJ_CITY_vs_STATE_susceptibility.png"), p_susc, width = 220, height = 120, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "RJ_CITY_vs_STATE_susceptibility.png"))
