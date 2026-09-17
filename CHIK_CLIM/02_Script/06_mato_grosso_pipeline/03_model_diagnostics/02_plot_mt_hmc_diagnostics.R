# Mato Grosso -- HMC/chain diagnostic figures for the final (ad098) global-q
# case-only fit: per-chain traceplots (with divergence markers) for the key
# structural parameters, plus an Rhat/ESS summary plot. Same bayesplot
# convention as
# 33_pernambuco_v4_9_replication/scripts/14_plot_pe_ridge_reparam_traceplots.R.

required_packages <- c("rstan", "bayesplot", "ggplot2", "dplyr", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(bayesplot); library(ggplot2); library(dplyr); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "03_Output/06_mato_grosso_pipeline/model_fits/v4_9_replication")
figure_dir <- file.path(root, "03_Output/06_mato_grosso_pipeline/figures/mato_grosso_v4_9_replication")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

b <- readRDS(file.path(base_dir, "outputs/caseonly/mt_global_q_case_only.rds"))
fit <- b$fit

pars_to_trace <- c("logit_q", "alpha_R", "phi_obs", "beta_sin1", "beta_cos1", "sigma_season_year", "lp__")
posterior_array <- as.array(fit)
np <- bayesplot::nuts_params(fit)

p_trace <- bayesplot::mcmc_trace(posterior_array, pars = pars_to_trace, np = np,
                                  np_style = bayesplot::trace_style_np(div_color = "red", div_size = 0.4, div_alpha = 0.8)) +
  ggtitle("Mato Grosso global-q case-only (adapt_delta=0.98): chain traceplots",
          subtitle = sprintf("%d chains x %d post-warmup draws | divergences=%d (red ticks, if any) | HMC %s",
                              b$config$chains, b$config$sampling, b$hmc$divergences,
                              if (b$hmc$hmc_pass) "PASS" else "FAIL"))
ggsave(file.path(figure_dir, "MT_global_q_traceplots.png"), p_trace, width = 260, height = 200, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "MT_global_q_traceplots.png"))

# ---- Rhat / ESS summary across ALL model parameters (not just the traced subset) ----
mon <- rstan::monitor(fit, print = FALSE)
rhat_vals <- mon[, "Rhat"]
ess_bulk <- mon[, "Bulk_ESS"]
ess_tail <- mon[, "Tail_ESS"]

p_rhat <- bayesplot::mcmc_rhat_hist(rhat_vals) +
  ggtitle("Rhat distribution across all parameters", subtitle = sprintf("max Rhat = %.4f (threshold 1.01)", max(rhat_vals, na.rm = TRUE)))
p_ess_bulk <- bayesplot::mcmc_neff_hist(ess_bulk / (b$config$chains * b$config$sampling)) +
  ggtitle("Bulk ESS ratio (ESS / total draws)", subtitle = sprintf("min bulk ESS = %.0f", min(ess_bulk, na.rm = TRUE)))
p_ess_tail <- bayesplot::mcmc_neff_hist(ess_tail / (b$config$chains * b$config$sampling)) +
  ggtitle("Tail ESS ratio (ESS / total draws)", subtitle = sprintf("min tail ESS = %.0f", min(ess_tail, na.rm = TRUE)))

sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
accept_df <- dplyr::bind_rows(lapply(seq_along(sp), function(i) {
  data.frame(chain = factor(i), accept_stat = sp[[i]][, "accept_stat__"], stepsize = sp[[i]][, "stepsize__"][1])
}))
p_accept <- ggplot(accept_df, aes(chain, accept_stat, fill = chain)) +
  geom_violin(alpha = 0.6, show.legend = FALSE) +
  labs(title = "Per-chain acceptance statistic distribution", x = "Chain", y = "accept_stat__") +
  theme_classic(base_size = 10)

figure <- (p_rhat | p_ess_bulk) / (p_ess_tail | p_accept)
ggsave(file.path(figure_dir, "MT_global_q_HMC_summary.png"), figure, width = 220, height = 180, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "MT_global_q_HMC_summary.png"))
