# STEP 2 (deep-dive) -- P2 chain-level mixing diagnostics, and the
# posterior geometry between tau_year (year-effect SD) and
# sigma_season_year (seasonal-amplitude SD). Read-only against the
# ORIGINAL ce_tau_year_P2.rds (adapt_delta=0.95) -- no new Stan fitting.

required_packages <- c("here", "rstan", "dplyr", "readr", "tibble", "posterior", "bayesplot", "ggplot2", "tidyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(rstan); library(dplyr); library(readr); library(tibble); library(posterior); library(bayesplot); library(ggplot2); library(tidyr) })

source(file.path(ROOT, "02_Script/02_ceara_pipeline/14_hmc_reliability_audit/00_config_hmc_audit.R"))
OKABE <- c(blue = "#0072B2", orange = "#D55E00", yellow = "#E69F00", green = "#009E73", purple = "#CC79A7", grey = "#56B4E9")

bundle <- readRDS(MODELS_TO_AUDIT$P2$path)
fit <- bundle$fit
Y <- bundle$stan_data$Y
message(sprintf("[P2 deep-dive] Y=%d years, %d chains, %d post-warmup iter/chain", Y, fit@sim$chains, fit@sim$iter - fit@sim$warmup))

# ================================================================
# A. Chain-level mixing: traceplots, rank plots, ACF, per-chain summary,
#    explicit between/within variance decomposition (classic Gelman-Rubin).
# ================================================================
draws_mix <- as_draws_array(fit, pars = c("tau_year", "sigma_season_year", "alpha_R"))

p_trace <- mcmc_trace(draws_mix, pars = c("tau_year", "sigma_season_year"), facet_args = list(nrow = 2)) +
  labs(title = "P2 (primary): chain traceplots") + theme_classic(base_size = 9)
ggsave(file.path(DIR_FIGURES_HA, "P2_mixing_traceplot.png"), p_trace, width = 240, height = 150, units = "mm", dpi = 300, bg = "white")

p_rank <- mcmc_rank_overlay(draws_mix, pars = c("tau_year", "sigma_season_year")) +
  labs(title = "P2 (primary): rank plots (uniform = well-mixed chains)") + theme_classic(base_size = 9)
ggsave(file.path(DIR_FIGURES_HA, "P2_mixing_rankplot.png"), p_rank, width = 240, height = 120, units = "mm", dpi = 300, bg = "white")

p_acf <- mcmc_acf(draws_mix, pars = c("tau_year", "sigma_season_year"), lags = 30) +
  labs(title = "P2 (primary): within-chain autocorrelation") + theme_classic(base_size = 9)
ggsave(file.path(DIR_FIGURES_HA, "P2_mixing_acf.png"), p_acf, width = 240, height = 150, units = "mm", dpi = 300, bg = "white")

gelman_rubin_decompose <- function(draws_arr, par) {
  m <- dim(draws_arr)[2]; n <- dim(draws_arr)[1]
  chain_vals <- lapply(seq_len(m), function(c) as.numeric(draws_arr[, c, par]))
  chain_means <- vapply(chain_vals, mean, numeric(1)); chain_vars <- vapply(chain_vals, var, numeric(1))
  overall_mean <- mean(unlist(chain_vals))
  W <- mean(chain_vars)
  B_over_n <- sum((chain_means - overall_mean)^2) / (m - 1)
  var_hat <- (n - 1) / n * W + B_over_n
  tibble(parameter = par, chain = seq_len(m), chain_mean = chain_means, chain_sd = sqrt(chain_vars)) |>
    mutate(within_chain_var_W = W, between_chain_var_B_over_n = B_over_n, pooled_var_hat = var_hat,
           rhat_classic = sqrt(var_hat / W), frac_variance_between = B_over_n / var_hat)
}
gr_tau <- gelman_rubin_decompose(draws_mix, "tau_year")
gr_sigma <- gelman_rubin_decompose(draws_mix, "sigma_season_year")
gr_all <- bind_rows(gr_tau, gr_sigma)
write_csv(gr_all, file.path(DIR_TABLES_HA, "P2_mixing_gelman_rubin_decomposition.csv"))
message("\n=== P2: Gelman-Rubin within/between-chain variance decomposition ===")
print(as.data.frame(gr_all |> dplyr::select(parameter, chain, chain_mean, chain_sd, within_chain_var_W, between_chain_var_B_over_n, rhat_classic, frac_variance_between)), digits = 4)

ess_per_chain <- bind_rows(lapply(c("tau_year", "sigma_season_year"), function(par) {
  bind_rows(lapply(seq_len(dim(draws_mix)[2]), function(c) {
    x <- as_draws_array(draws_mix[, c, par, drop = FALSE])
    tibble(parameter = par, chain = c, ess_bulk_within_chain = ess_bulk(x), ess_tail_within_chain = ess_tail(x), autocorr_lag1 = cor(as.numeric(x)[-1], as.numeric(x)[-length(as.numeric(x))]))
  }))
}))
write_csv(ess_per_chain, file.path(DIR_TABLES_HA, "P2_mixing_ess_per_chain.csv"))
message("\n=== P2: within-chain ESS and lag-1 autocorrelation ==="); print(as.data.frame(ess_per_chain), digits = 4)

# ================================================================
# B. Posterior geometry: tau_year vs sigma_season_year joint shape.
# ================================================================
ext <- rstan::extract(fit, pars = c("tau_year", "sigma_season_year", "alpha_R"), permuted = FALSE)
n_chain <- dim(ext)[2]
geom_df <- bind_rows(lapply(seq_len(n_chain), function(c) tibble(chain = factor(c), tau_year = ext[, c, "tau_year"], sigma_season_year = ext[, c, "sigma_season_year"], alpha_R = ext[, c, "alpha_R"])))
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
geom_df$divergent <- unlist(lapply(sp, function(x) x[, "divergent__"]))

overall_cor <- cor(geom_df$tau_year, geom_df$sigma_season_year, method = "pearson")
overall_cor_spearman <- cor(geom_df$tau_year, geom_df$sigma_season_year, method = "spearman")
message(sprintf("\n[P2 geometry] tau_year vs sigma_season_year: Pearson r=%.4f, Spearman rho=%.4f (pooled across chains, post-warmup)", overall_cor, overall_cor_spearman))

p_joint <- ggplot(geom_df, aes(tau_year, sigma_season_year)) +
  geom_point(aes(colour = chain), alpha = 0.25, size = 0.7) +
  geom_point(data = geom_df |> dplyr::filter(divergent == 1), colour = OKABE["orange"], shape = 4, size = 2, stroke = 1.1) +
  geom_density_2d(colour = "grey30", linewidth = 0.3) +
  scale_colour_manual(values = unname(OKABE[c("blue","green","purple","yellow")])) +
  labs(title = "P2: joint posterior of tau_year and sigma_season_year",
       subtitle = sprintf("Pearson r=%.3f, Spearman rho=%.3f | orange X = divergent transitions | colour = chain", overall_cor, overall_cor_spearman),
       x = "tau_year (year-effect SD)", y = "sigma_season_year (seasonal-amplitude SD)") +
  theme_classic(base_size = 9)
ggsave(file.path(DIR_FIGURES_HA, "P2_geometry_tau_vs_sigma_season.png"), p_joint, width = 200, height = 160, units = "mm", dpi = 300, bg = "white")

# Heteroscedasticity / funnel check: does the SPREAD of sigma_season_year
# change across tau_year bins? (classic funnel signature if conditional
# variance shrinks/grows systematically with the other parameter)
geom_df <- geom_df |> mutate(tau_bin = cut(tau_year, breaks = quantile(tau_year, probs = seq(0, 1, 0.2)), include.lowest = TRUE, labels = paste0("Q", 1:5)))
funnel_check <- geom_df |> group_by(tau_bin) |> summarise(
  n = n(), mean_tau = mean(tau_year), median_sigma_season = median(sigma_season_year),
  sd_sigma_season = sd(sigma_season_year), iqr_sigma_season = IQR(sigma_season_year),
  mean_alpha_R = mean(alpha_R), sd_alpha_R = sd(alpha_R), n_divergent = sum(divergent), .groups = "drop")
write_csv(funnel_check, file.path(DIR_TABLES_HA, "P2_geometry_funnel_check_by_tau_quintile.csv"))
message("\n=== P2 geometry: sigma_season_year (and alpha_R) spread by tau_year quintile ===")
print(as.data.frame(funnel_check), digits = 4)

# Per-year year_effect vs A_year correlation (the quantities that
# ACTUALLY enter log_R0 for a given year, not just the hyperparameters).
ext_ye <- rstan::extract(fit, pars = c("z_year", "z_season_year", "tau_year", "sigma_season_year"), permuted = TRUE)
n_draws <- length(ext_ye$tau_year)
message(sprintf("[P2 geometry] dim(z_year)=%s, dim(z_season_year)=%s, length(tau_year)=%d, Y=%d",
                 paste(dim(ext_ye$z_year), collapse = "x"), paste(dim(ext_ye$z_season_year), collapse = "x"), n_draws, Y))
z_year_mat <- matrix(as.numeric(ext_ye$z_year), nrow = n_draws, ncol = Y)
z_season_mat <- matrix(as.numeric(ext_ye$z_season_year), nrow = n_draws, ncol = Y)
z_year_centered <- z_year_mat - matrix(rowMeans(z_year_mat), nrow = n_draws, ncol = Y)
year_effect_draws <- sweep(z_year_centered, 1, ext_ye$tau_year, "*")                       # [draws, Y]
A_year_draws <- exp(sweep(z_season_mat, 1, ext_ye$sigma_season_year, "*"))                 # [draws, Y]
per_year_cor <- tibble(year_index = seq_len(Y),
                        cor_year_effect_vs_A_year = vapply(seq_len(Y), function(y) cor(year_effect_draws[, y], A_year_draws[, y]), numeric(1)))
write_csv(per_year_cor, file.path(DIR_TABLES_HA, "P2_geometry_year_effect_vs_A_year_per_year_cor.csv"))
message("\n=== P2 geometry: per-year correlation between realised year_effect[y] and A_year[y] (the quantities entering log_R0) ===")
print(as.data.frame(per_year_cor), digits = 4)
message(sprintf("Mean per-year correlation: %.4f | range: [%.4f, %.4f]", mean(per_year_cor$cor_year_effect_vs_A_year), min(per_year_cor$cor_year_effect_vs_A_year), max(per_year_cor$cor_year_effect_vs_A_year)))

message("\n[P2 chain mixing + geometry] Complete -> ", DIR_TABLES_HA, " / ", DIR_FIGURES_HA)
