# Rio de Janeiro CITY auxiliary -- Model A diagnostics: HMC gate, q
# identification, case PPC, trajectory figure. Mirrors 05_rj_global_q_diagnostics.R.

required_packages <- c("rstan", "dplyr", "readr", "tibble", "ggplot2", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(readr); library(tibble); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/34_rio_de_janeiro_v4_9_replication")
table_dir <- file.path(root, "03_Output/tables/rio_de_janeiro_v4_9_replication/city")
figure_dir <- file.path(root, "03_Output/figures/rio_de_janeiro_v4_9_replication/city")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

b <- readRDS(file.path(base_dir, "outputs/city_caseonly/rj_city_global_q_case_only.rds"))
fit <- b$fit
weekly <- b$weekly_data
dates <- as.Date(weekly$week_start)

message("=== Model A (CITY case-only): HMC gate ===")
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
hmc_by_chain <- bind_rows(lapply(seq_along(sp), function(i) {
  tibble(chain = i, divergent = sum(sp[[i]][, "divergent__"]),
         max_treedepth_hits = sum(sp[[i]][, "treedepth__"] >= b$config$max_treedepth),
         mean_stepsize = mean(sp[[i]][, "stepsize__"]), mean_accept = mean(sp[[i]][, "accept_stat__"]))
}))
hmc_by_chain$bfmi <- rstan::get_bfmi(fit)
print(as.data.frame(hmc_by_chain), digits = 4)
message(sprintf("\nOverall: divergences=%d, max_treedepth_hits=%d, max Rhat=%.4f, min bulk ESS=%.1f, hmc_pass=%s",
                 b$hmc$divergences, b$hmc$max_treedepth_hits, b$hmc$maximum_rhat, b$hmc$minimum_bulk_ess, b$hmc$hmc_pass))
write_csv(hmc_by_chain, file.path(table_dir, "RJ_CITY_caseonly_HMC.csv"))

draws <- rstan::extract(fit, pars = c("logit_q", "alpha_R", "S_prop", "immune_prop", "R0_t", "expected_reported_cases", "C_pred"), permuted = TRUE)
q_draws <- plogis(draws$logit_q)
cor_q_alpha <- cor(draws$logit_q, draws$alpha_R)
u10_window_idx <- which(dates >= as.Date("2018-07-01") & dates <= as.Date("2018-10-31"))
p_city_u10window <- rowMeans(draws$immune_prop[, u10_window_idx, drop = FALSE])
immune_2025 <- draws$immune_prop[, ncol(draws$immune_prop)]
min_S <- apply(draws$S_prop, 1, min)
max_R0 <- apply(draws$R0_t, 1, max)

message("\n=== Model A identification diagnostics ===")
message(sprintf("q_city: median=%.4f, 95%% CrI=[%.4f,%.4f]", median(q_draws), quantile(q_draws,.025), quantile(q_draws,.975)))
message(sprintf("cor(logit_q, alpha_R) = %.4f", cor_q_alpha))
message(sprintf("immune fraction, U10 window: median=%.4f, 95%% CrI=[%.4f,%.4f]", median(p_city_u10window), quantile(p_city_u10window,.025), quantile(p_city_u10window,.975)))
message(sprintf("immune_2025: median=%.4f, 95%% CrI=[%.4f,%.4f]", median(immune_2025), quantile(immune_2025,.025), quantile(immune_2025,.975)))
message(sprintf("min S/N: median=%.4f | max R0: median=%.4f", median(min_S), median(max_R0)))

summary_tbl <- tibble(model = "CITY case-only", q_median = median(q_draws), q_lo95 = quantile(q_draws,.025), q_hi95 = quantile(q_draws,.975),
                       cor_logitq_alphaR = cor_q_alpha, u10window_immune_median = median(p_city_u10window),
                       u10window_immune_lo95 = quantile(p_city_u10window,.025), u10window_immune_hi95 = quantile(p_city_u10window,.975),
                       immune_2025_median = median(immune_2025), min_S_median = median(min_S), max_R0_median = median(max_R0),
                       hmc_pass = b$hmc$hmc_pass, divergences = b$hmc$divergences)
write_csv(summary_tbl, file.path(table_dir, "RJ_CITY_caseonly_summary.csv"))

# ---- Case PPC + trajectories figure ----
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
summarise_traj <- function(m) tibble(median = apply(m, 2, median), lo95 = apply(m, 2, quantile, .025), hi95 = apply(m, 2, quantile, .975))

pred_summary <- summarise_traj(draws$C_pred); pred_summary$week_start <- dates
exp_summary <- summarise_traj(draws$expected_reported_cases); exp_summary$week_start <- dates
obs_df <- tibble(week_start = dates, observed = weekly$cases)

p_cases <- ggplot() +
  geom_ribbon(data = pred_summary, aes(week_start, ymin = lo95, ymax = hi95), fill = "#56B4E9", alpha = 0.2) +
  geom_line(data = exp_summary, aes(week_start, median), colour = "#315A7D", linewidth = 0.5) +
  geom_point(data = obs_df, aes(week_start, observed), colour = "black", size = 0.4, alpha = 0.6) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "RJ CITY (Model A, case-only): weekly case PPC", x = NULL, y = "Weekly reported cases") + theme_v4

S_summary <- summarise_traj(draws$S_prop); S_summary$week_start <- dates
U_summary <- summarise_traj(draws$immune_prop); U_summary$week_start <- dates
p_immunity <- ggplot() +
  geom_ribbon(data = S_summary, aes(week_start, ymin = lo95, ymax = hi95), fill = "#009E73", alpha = 0.15) +
  geom_line(data = S_summary, aes(week_start, median, colour = "Susceptible"), linewidth = 0.6) +
  geom_ribbon(data = U_summary, aes(week_start, ymin = lo95, ymax = hi95), fill = "#76558F", alpha = 0.15) +
  geom_line(data = U_summary, aes(week_start, median, colour = "Immune"), linewidth = 0.6) +
  scale_colour_manual(values = c("Susceptible" = "#009E73", "Immune" = "#76558F"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1)) +
  labs(title = "RJ CITY: susceptibility and immunity (case-only)", x = NULL, y = "Population proportion") + theme_v4 + theme(legend.position = "top")

q_band <- tibble(week_start = dates, lo95 = quantile(q_draws, .025), hi95 = quantile(q_draws, .975), median = median(q_draws))
p_q <- ggplot(q_band, aes(week_start)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#E69F00", alpha = 0.25) + geom_line(aes(y = median), colour = "#E69F00", linewidth = 0.9) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA)) +
  labs(title = sprintf("q_city (case-only) posterior: median=%.4f, 95%% CrI=[%.4f,%.4f]", median(q_draws), quantile(q_draws,.025), quantile(q_draws,.975)),
       x = NULL, y = "q") + theme_v4

figure <- p_cases / p_immunity / p_q +
  patchwork::plot_annotation(title = "Rio de Janeiro CITY -- Model A (case-only, HMC gate: see title)",
                              subtitle = sprintf("%s | q_city median=%.4f | cor(logit_q,alpha_R)=%.3f",
                                                  if (b$hmc$hmc_pass) "HMC PASS" else sprintf("HMC FAIL (div=%d)", b$hmc$divergences),
                                                  median(q_draws), cor_q_alpha))
ggsave(file.path(figure_dir, "RJ_CITY_caseonly_trajectories.png"), figure, width = 200, height = 220, units = "mm", dpi = 300, bg = "white")
message("\n[saved] ", file.path(figure_dir, "RJ_CITY_caseonly_trajectories.png"))
