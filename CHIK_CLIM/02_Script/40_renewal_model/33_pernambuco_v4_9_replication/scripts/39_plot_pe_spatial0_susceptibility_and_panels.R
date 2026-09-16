# Pernambuco Spatial-0 -- regional susceptibility curves + the established
# 6/7-panel trajectory figure (same layout/theme as script 07), built at
# the STATE-AGGREGATE level (summed across the 3 chains that reached
# post-warmup sampling: chains 2, 4, 5 -- LOW/INTERMEDIATE/HIGH branches,
# which the branch-agreement diagnostic (script 38) confirmed converge to
# the same posterior). Chains 1, 3, 6 never reached sampling and are
# excluded (not "salvaging" a mode -- these 3 chains all reached the SAME
# region, per script 38).

required_packages <- c("rstan", "ggplot2", "dplyr", "patchwork", "scales", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(ggplot2); library(dplyr); library(patchwork); library(tibble); library(readr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
pe_root <- file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication")
chains_dir <- file.path(pe_root, "outputs/spatial0/chains")
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/spatial0")
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_model")
turnover_table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_replication/spatial_turnover_diagnostic")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

COL_INK <- "#202124"; COL_MUTED <- "#6B7280"; COL_GRID <- "#E5E7EB"
COL_NAVY <- "#315A7D"; COL_BLUE <- "#56B4E9"; COL_VERMILLION <- "#D55E00"
COL_ORANGE <- "#E69F00"; COL_GREEN <- "#009E73"; COL_PURPLE <- "#76558F"
STRATUM_LEVELS <- c("1_Recife", "2_Metropolitana_remainder", "3_Agreste", "4_Sertao", "5_Vale_Sao_Francisco_Araripe")
STRATUM_COLOURS <- c("1_Recife" = "#D55E00", "2_Metropolitana_remainder" = "#E69F00",
                      "3_Agreste" = "#009E73", "4_Sertao" = "#0072B2", "5_Vale_Sao_Francisco_Araripe" = "#CC79A7")

theme_publication <- function(base_size = 8.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, color = COL_MUTED, margin = margin(b = 7)),
      axis.title = element_text(size = base_size, color = COL_INK),
      axis.text = element_text(size = base_size - 0.5, color = COL_MUTED),
      axis.line = element_line(color = COL_INK, linewidth = 0.35),
      axis.ticks = element_line(color = COL_INK, linewidth = 0.35),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "top", legend.justification = "left", legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.width = grid::unit(13, "pt"), legend.key.height = grid::unit(7, "pt"),
      plot.margin = margin(7, 8, 7, 7)
    )
}
summarise_trajectory <- function(draw_matrix, prefix) {
  tibble(variable = prefix, t = seq_len(ncol(draw_matrix)),
         q025 = apply(draw_matrix, 2, quantile, .025), q25 = apply(draw_matrix, 2, quantile, .25),
         median = apply(draw_matrix, 2, median),
         q75 = apply(draw_matrix, 2, quantile, .75), q975 = apply(draw_matrix, 2, quantile, .975))
}

message("Reading chains 2, 4, 5 (the 3 chains that reached post-warmup sampling) via read_stan_csv()...")
csv_files <- file.path(chains_dir, sprintf("chain_%d.csv", c(2, 4, 5)))
fit <- rstan::read_stan_csv(csv_files)
message("Loaded stanfit with ", fit@sim$chains, " chains, ", fit@sim$n_save[1] - fit@sim$warmup2[1], " post-warmup draws/chain.")

weeks <- readRDS(file.path(turnover_table_dir, "pe_municipality_week_panel.rds")) |> pull(week_start) |> unique() |> sort()
stan_data <- readRDS(file.path(table_dir, "PE_5strata_stan_data.rds"))$stan_data # same 5-strata data as Spatial-1/0
R <- length(STRATUM_LEVELS)

draws <- rstan::extract(fit, pars = c("S_prop", "R0_t", "R_eff_t", "S_PE_prop", "U_PE_prop",
                                       "C_pred", "C_pred_state_sum", "expected_reported_cases", "p_site_window"),
                         permuted = TRUE)
q_draws <- plogis(as.vector(rstan::extract(fit, "logit_q")$logit_q))
q_prior_draws <- plogis(as.vector(rstan::extract(fit, "logit_q_prior_draw")$logit_q_prior_draw))

# ================================================================
# 1. Regional susceptibility curves (explicit user request)
# ================================================================
S_regional <- bind_rows(lapply(seq_len(R), function(r) {
  summarise_trajectory(draws$S_prop[, , r], STRATUM_LEVELS[r]) |> mutate(week_start = weeks[t])
}))
S_statewide <- summarise_trajectory(draws$S_PE_prop, "Statewide (S_PE)") |> mutate(week_start = weeks[t])

date_breaks <- seq(as.Date("2015-01-01"), as.Date("2026-01-01"), by = "1 year")
x_scale <- scale_x_date(breaks = date_breaks, date_labels = "%Y", expand = expansion(mult = c(0.01, 0.015)))

p_S_regional <- ggplot(S_regional, aes(week_start, median, colour = variable)) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = variable), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = STRATUM_COLOURS, name = NULL) + scale_fill_manual(values = STRATUM_COLOURS, guide = "none") +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
  x_scale +
  labs(title = "Pernambuco Spatial-0: regional susceptible proportion (S/N) by stratum",
       subtitle = "Posterior median, 95% CrI. Common R0(t) across strata -- differences arise ONLY from separate S/U histories.",
       x = NULL, y = "S / N") +
  theme_publication() + theme(legend.position = "top")

p_S_statewide <- ggplot(S_statewide, aes(week_start, median)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = "grey40", alpha = 0.18) +
  geom_line(colour = "grey15", linewidth = 0.7) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
  x_scale +
  labs(title = "Statewide susceptible proportion S_PE(t) (GENERATED from the 5 regional pools)", x = NULL, y = "S_PE / N_PE") +
  theme_publication()

fig_susc <- p_S_regional / p_S_statewide
ggsave(file.path(figure_dir, "PE_spatial0_regional_susceptibility.png"), fig_susc, width = 200, height = 200, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "PE_spatial0_regional_susceptibility.png"))

# ================================================================
# 2. Established 7-panel trajectory figure (state-aggregate level)
# ================================================================
observed_state <- rowSums(stan_data$C)
C_pred_state <- draws$C_pred_state_sum
expected_state <- apply(draws$expected_reported_cases, c(1, 2), sum) # sum over strata -> [iter, N]

cases_expected_summary <- summarise_trajectory(expected_state, "expected") |> mutate(week_start = weeks[t])
cases_predictive_summary <- summarise_trajectory(C_pred_state, "predictive") |> mutate(week_start = weeks[t])
cases_df <- tibble(week_start = weeks, observed = observed_state)

p_cases <- ggplot() +
  geom_ribbon(data = cases_predictive_summary, aes(week_start, ymin = q025, ymax = q975), fill = COL_BLUE, alpha = 0.15) +
  geom_ribbon(data = cases_expected_summary, aes(week_start, ymin = q25, ymax = q75), fill = COL_NAVY, alpha = 0.24) +
  geom_line(data = cases_expected_summary, aes(week_start, median, colour = "Model expectation"), linewidth = 0.65) +
  geom_point(data = cases_df, aes(week_start, observed, colour = "Observed cases"), size = 0.5, alpha = 0.65) +
  scale_colour_manual(values = c("Observed cases" = COL_INK, "Model expectation" = COL_NAVY)) +
  scale_y_continuous(labels = scales::label_number(big.mark = ","), expand = expansion(mult = c(0, 0.06))) +
  x_scale +
  labs(title = "Reported cases and posterior prediction (state total = SUM of 5 regional predictions)",
       subtitle = "Likelihood fit to the 5 regional series only (Section 10) -- this is a posterior predictive check, not a fit target",
       x = NULL, y = "Weekly reported cases") +
  theme_publication()

immunity <- bind_rows(
  summarise_trajectory(draws$S_PE_prop, "Susceptible") |> mutate(week_start = weeks[t]),
  summarise_trajectory(draws$U_PE_prop, "Infection-derived immune (statewide)") |> mutate(week_start = weeks[t])
)
p_immunity <- ggplot(immunity, aes(week_start, colour = variable, fill = variable)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), alpha = 0.10, colour = NA) +
  geom_line(aes(y = median), linewidth = 0.7) +
  scale_colour_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune (statewide)" = COL_PURPLE)) +
  scale_fill_manual(values = c("Susceptible" = COL_GREEN, "Infection-derived immune (statewide)" = COL_PURPLE)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent(accuracy = 1), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0.01))) +
  x_scale +
  labs(title = "Statewide susceptibility and immunity (GENERATED from 5 regional pools)", x = NULL, y = "Population proportion") +
  theme_publication()

reproduction <- bind_rows(
  summarise_trajectory(draws$R0_t, "R0(t)") |> mutate(week_start = weeks[t]), # common across strata in Spatial-0
  bind_rows(lapply(seq_len(R), function(r) summarise_trajectory(draws$R_eff_t[, , r], paste0("Reff: ", STRATUM_LEVELS[r])))) |> mutate(week_start = weeks[t])
)
reff_colours <- setNames(STRATUM_COLOURS, paste0("Reff: ", STRATUM_LEVELS))
repro_colours <- c("R0(t)" = COL_INK, reff_colours)
p_reproduction <- ggplot(reproduction, aes(week_start, colour = variable)) +
  geom_hline(yintercept = 1, colour = COL_MUTED, linewidth = 0.35, linetype = "22") +
  geom_line(aes(y = median), linewidth = 0.55) +
  scale_colour_manual(values = repro_colours, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.06))) +
  x_scale +
  labs(title = "Reproduction numbers: common R0(t) vs regional Reff(t)",
       subtitle = "R0(t) is IDENTICAL across strata by construction (Section F) -- only Reff(t)=R0(t)*S/N differs regionally",
       x = NULL, y = "Reproduction number") +
  theme_publication() + theme(legend.text = element_text(size = 6))

q_band <- tibble(week_start = weeks, post_lo95 = quantile(q_draws, .025), post_hi95 = quantile(q_draws, .975),
                  post_lo50 = quantile(q_draws, .25), post_hi50 = quantile(q_draws, .75), post_median = median(q_draws))
p_ascertainment <- ggplot(q_band, aes(week_start)) +
  geom_ribbon(aes(ymin = post_lo95, ymax = post_hi95), fill = COL_ORANGE, alpha = 0.18) +
  geom_ribbon(aes(ymin = post_lo50, ymax = post_hi50), fill = COL_ORANGE, alpha = 0.30) +
  geom_line(aes(y = post_median), colour = COL_ORANGE, linewidth = 0.9) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
  x_scale +
  labs(title = "Case ascertainment (q_PE), one common value across all 5 strata",
       subtitle = sprintf("Chains 2/4/5 (LOW/INTERMEDIATE/HIGH branches) agree: median=%.3f, 95%% CrI=[%.3f, %.3f]",
                           median(q_draws), quantile(q_draws, .025), quantile(q_draws, .975)),
       x = NULL, y = "Ascertainment fraction (q)") +
  theme_publication()

attack <- summarise_trajectory(draws$U_PE_prop, "immune") |> mutate(week_start = weeks[t])
p_attack <- ggplot(attack, aes(week_start)) +
  geom_ribbon(aes(ymin = q025, ymax = q975), fill = COL_GREEN, alpha = 0.14) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = COL_GREEN, alpha = 0.24) +
  geom_line(aes(y = median), colour = COL_GREEN, linewidth = 0.7) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA), expand = expansion(mult = c(0, 0.07))) +
  x_scale +
  labs(title = "Cumulative infection (statewide immune fraction)", x = NULL, y = "Cumulative proportion infected") +
  theme_publication()

p_site <- draws$p_site_window[, 1]
obs_prev <- stan_data$sero_n_positive[1] / stan_data$sero_n_tested[1]
obs_ci <- prop.test(stan_data$sero_n_positive[1], stan_data$sero_n_tested[1])$conf.int
sero_df <- tibble(quantity = c("Observed (U14)", "Model (Recife stratum, direct link)"),
                   median = c(obs_prev, median(p_site)), lo = c(obs_ci[1], quantile(p_site, .025)), hi = c(obs_ci[2], quantile(p_site, .975)))
p_sero <- ggplot(sero_df, aes(quantity, median, colour = quantity)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), fatten = 2) +
  scale_colour_manual(values = c("Observed (U14)" = COL_INK, "Model (Recife stratum, direct link)" = COL_PURPLE)) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(title = "U14 Recife serosurvey: observed vs Spatial-0 (direct Recife-stratum link, no eta_geo)", x = NULL, y = "seroprevalence") +
  theme_publication(base_size = 8) + theme(legend.position = "none")

figure <- (p_cases | p_reproduction) / (p_immunity | p_ascertainment) / (p_attack | p_sero) +
  patchwork::plot_annotation(title = "Pernambuco Spatial-0: 5-strata separate susceptible pools, one common q_PE (2015-2025)",
                              subtitle = sprintf("Chains 2, 4, 5 only (LOW/INTERMEDIATE/HIGH branches -- all converged to the same posterior, script 38). Chains 1,3,6 did not reach sampling."),
                              tag_levels = "A",
                              theme = theme(plot.title = element_text(size = 12, face = "bold", margin = margin(b = 3)),
                                            plot.subtitle = element_text(size = 8.5, colour = COL_MUTED, margin = margin(b = 8)),
                                            plot.tag = element_text(size = 11, face = "bold")))

ggsave(file.path(figure_dir, "PE_spatial0_trajectories_6panel.png"), figure, width = 200, height = 260, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "PE_spatial0_trajectories_6panel.png"))

message(sprintf("\nq_PE (chains 2/4/5): median=%.4f, 95%% CrI=[%.4f, %.4f]", median(q_draws), quantile(q_draws,.025), quantile(q_draws,.975)))
message(sprintf("Statewide immune_2025: median=%.4f", median(draws$U_PE_prop[, ncol(draws$U_PE_prop)])))
message(sprintf("Recife U14 fit: model=%.4f [%.4f,%.4f] vs observed=%.4f [%.4f,%.4f]",
                 median(p_site), quantile(p_site,.025), quantile(p_site,.975), obs_prev, obs_ci[1], obs_ci[2]))
