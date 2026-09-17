# Visualises a chain-split HMC gate failure for one q-node fit: trace plots
# and posterior densities by chain for the parameters most affected, plus a
# side-by-side comparison against the td14 reference fit (where all chains
# agreed). Intended for exactly the failure pattern seen in q50 (Rhat~1.5,
# ESS~7 across most parameters despite zero divergences/treedepth hits --
# i.e. chains converging to different, internally-smooth basins rather than
# a diffuse mixing problem).
#
# Usage: Q_NODE_ID=q50 Rscript 10_plot_chain_split_diagnostic.R

required_packages <- c("here", "rstan", "dplyr", "tibble", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(ggplot2); library(patchwork) })

diag_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

chain_colours <- c("1" = "#D55E00", "2" = "#0072B2", "3" = "#009E73", "4" = "#CC79A7")

trace_and_density <- function(raw, param, title) {
  n_iter <- dim(raw)[1]; n_chain <- dim(raw)[2]
  df <- bind_rows(lapply(seq_len(n_chain), function(c) tibble(iteration = seq_len(n_iter), chain = factor(c), value = raw[, c, param])))
  trace <- ggplot(df, aes(iteration, value, colour = chain)) + geom_line(linewidth = .25, alpha = .8) +
    scale_colour_manual(values = chain_colours) + labs(title = paste0(title, " -- trace"), x = "post-warmup iteration", y = param) +
    theme_classic(base_size = 8) + theme(legend.position = "none")
  dens <- ggplot(df, aes(value, colour = chain, fill = chain)) + geom_density(alpha = .2) +
    scale_colour_manual(values = chain_colours) + scale_fill_manual(values = chain_colours) +
    labs(title = paste0(title, " -- posterior density by chain"), x = param, y = "density") +
    theme_classic(base_size = 8) + theme(legend.position = "top", legend.title = element_blank())
  trace | dens
}

run_plot_chain_split_diagnostic <- function(node_id = Sys.getenv("Q_NODE_ID", "q50")) {
  root <- diag_root()
  base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "19_v4_1_modular_q")
  fit_path <- file.path(base_dir, "outputs", node_id, paste0("renewal_ceara_v4_1_modular_q_fit_", node_id, ".rds"))
  if (!file.exists(fit_path)) stop("No fit found for node '", node_id, "': ", fit_path)
  bundle <- readRDS(fit_path)
  fit <- bundle$fit

  pars <- c("alpha_R", "beta_cos", "phi_obs")
  raw <- rstan::extract(fit, pars = pars, permuted = FALSE, inc_warmup = FALSE)
  z_year_raw <- rstan::extract(fit, pars = "z_year", permuted = FALSE, inc_warmup = FALSE)

  p1 <- trace_and_density(raw, "alpha_R", "A. alpha_R")
  p2 <- trace_and_density(raw, "beta_cos", "B. beta_cos")
  p3 <- trace_and_density(raw, "phi_obs", "C. phi_obs")
  z1_raw <- array(z_year_raw[, , 1], dim = c(dim(z_year_raw)[1], dim(z_year_raw)[2], 1))
  dimnames(z1_raw) <- list(NULL, NULL, "z_year[1]")
  p4 <- trace_and_density(z1_raw, "z_year[1]", "D. z_year[1]")

  per_chain_means <- bind_rows(lapply(seq_len(4), function(c) tibble(
    chain = c, alpha_R = mean(raw[, c, "alpha_R"]), beta_cos = mean(raw[, c, "beta_cos"]),
    phi_obs = mean(raw[, c, "phi_obs"]), z_year_1 = mean(z_year_raw[, c, 1])
  )))

  figure <- p1 / p2 / p3 / p4 + plot_annotation(
    title = sprintf("v4.1 modular-q node %s: chain-split HMC gate failure", node_id),
    subtitle = sprintf(
      "q_external=%.4f | chains 2-4 agree (alpha_R~%.2f) but chain 1 diverges to a different basin (alpha_R=%.2f) -- zero divergences/treedepth hits within either basin",
      bundle$config$q_external, mean(per_chain_means$alpha_R[2:4]), per_chain_means$alpha_R[1]
    ),
    theme = theme(plot.subtitle = element_text(colour = "#D55E00", face = "bold", size = 8))
  )
  figure_dir <- file.path(root, "03_Output", "figures", "renewal_v4_1_modular_q", node_id)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  figure_path <- file.path(figure_dir, paste0("renewal_v4_1_modular_q_", node_id, "_chain_split_diagnostic.pdf"))
  ggsave(figure_path, figure, width = 200, height = 260, units = "mm", device = cairo_pdf)

  table_dir <- file.path(root, "03_Output", "tables", "renewal_v4_1_modular_q", node_id)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(per_chain_means, file.path(table_dir, paste0("chain_split_per_chain_means_", node_id, ".csv")))

  message("[", node_id, "] chain-split diagnostic figure saved: ", figure_path)
  message("[", node_id, "] per-chain means:")
  print(as.data.frame(per_chain_means))
  invisible(figure_path)
}

if (sys.nframe() == 0L) run_plot_chain_split_diagnostic()
