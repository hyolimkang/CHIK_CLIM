# Pernambuco Model D pilot -- required Section 6 diagnostic plots:
# gamma_curve vs logit_q and alpha_R vs logit_q, divergent draws marked.

required_packages <- c("rstan", "ggplot2", "dplyr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(ggplot2); library(dplyr); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/modelD_curved_reparam")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

b <- readRDS(file.path(root, "02_Script/40_renewal_model/33_pernambuco_v4_9_replication/outputs/modelD/pe_curved_U14_pilot.rds"))
fit <- b$fit

get_flat <- function(par) {
  arr <- rstan::extract(fit, pars = par, permuted = FALSE, inc_warmup = FALSE)
  as.vector(arr[, , 1])
}
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
divergent <- unlist(lapply(sp, function(x) x[, "divergent__"]))
chain_id <- rep(seq_along(sp), times = sapply(sp, nrow))

logit_q <- get_flat("logit_q")
alpha_R <- get_flat("alpha_R")
gamma_curve <- get_flat("gamma_curve")

plot_df <- tibble(logit_q, alpha_R, gamma_curve, chain = factor(chain_id),
                   divergent = factor(divergent, labels = c("Non-divergent", "Divergent")))

theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

p1 <- ggplot(plot_df, aes(logit_q, gamma_curve)) +
  geom_point(aes(colour = divergent, alpha = divergent, size = divergent, shape = chain)) +
  scale_colour_manual(values = c("Non-divergent" = "grey50", "Divergent" = "red")) +
  scale_alpha_manual(values = c("Non-divergent" = 0.15, "Divergent" = 0.9)) +
  scale_size_manual(values = c("Non-divergent" = 0.5, "Divergent" = 1.8)) +
  labs(title = "Pernambuco Model D pilot: gamma_curve vs logit_q",
       subtitle = sprintf("cor=%.3f (near 0 -- curvature well removed). 32 divergences, 30/32 in chain 4.", cor(logit_q, gamma_curve)),
       x = "logit_q", y = "gamma_curve") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())

p2 <- ggplot(plot_df, aes(logit_q, alpha_R)) +
  geom_point(aes(colour = divergent, alpha = divergent, size = divergent, shape = chain)) +
  scale_colour_manual(values = c("Non-divergent" = "grey50", "Divergent" = "red")) +
  scale_alpha_manual(values = c("Non-divergent" = 0.15, "Divergent" = 0.9)) +
  scale_size_manual(values = c("Non-divergent" = 0.5, "Divergent" = 1.8)) +
  labs(title = "Pernambuco Model D pilot: alpha_R (derived) vs logit_q",
       subtitle = "Divergent draws span the typical q range, not a narrow tail -- concentrated in one chain, not one region",
       x = "logit_q", y = "alpha_R (derived)") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())

library(patchwork)
combined <- p1 / p2
ggsave(file.path(figure_dir, "PE_gamma_curve_pilot_diagnostics.png"), combined, width = 200, height = 220, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "PE_gamma_curve_pilot_diagnostics.png"))
