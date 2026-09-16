# Pernambuco Model C (ridge reparam) -- trace plots for key parameters,
# with divergent post-warmup draws marked (red ticks), to visualise where
# the 26 residual divergences (concentrated in chain 2, low-q region)
# occur relative to mixing.

required_packages <- c("rstan", "bayesplot", "ggplot2", "dplyr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(bayesplot); library(ggplot2); library(dplyr) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
FIT_TAG <- Sys.getenv("FIT_TAG", "pe_ridge_U14_full")
bundle <- readRDS(file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication/outputs/modelC", paste0(FIT_TAG, ".rds")))
fit <- bundle$fit

figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/modelC_ridge_reparam")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

pars_to_trace <- c("logit_q", "q", "alpha_R", "gamma_ridge", "phi_obs", "sigma_season_year")
posterior_array <- as.array(fit, pars = pars_to_trace)
np <- bayesplot::nuts_params(fit)

color_scheme_set("mix-blue-red")
p_trace <- bayesplot::mcmc_trace(posterior_array, pars = pars_to_trace, np = np, np_style = trace_style_np(div_color = "red", div_size = 0.4, div_alpha = 0.8)) +
  ggtitle("Pernambuco Model C (ridge reparam, FULL run): trace plots with divergent draws marked",
          subtitle = sprintf("%d post-warmup divergent transitions (mostly chain 2) -- concentrated near q's low end", rstan::get_num_divergent(fit)))
ggsave(file.path(figure_dir, "pernambuco_modelC_traceplots_full.png"), p_trace, width = 260, height = 200, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pernambuco_modelC_traceplots_full.png"))

# Pairs plot: q vs alpha_R vs gamma_ridge, divergences highlighted (the
# classic diagnostic for seeing WHERE in parameter space divergences occur).
p_pairs <- bayesplot::mcmc_pairs(posterior_array, pars = c("logit_q", "alpha_R", "gamma_ridge"),
                                  np = np, off_diag_args = list(size = 0.5, alpha = 0.3))
ggsave(file.path(figure_dir, "pernambuco_modelC_pairs_divergences_full.png"), p_pairs, width = 200, height = 200, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "pernambuco_modelC_pairs_divergences_full.png"))
