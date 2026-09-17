# q posterior comparison: prior vs WITH serology vs WITHOUT serology
# (J_sero=0 identification ablation). Diagnostic only.

required_packages <- c("rstan", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(ggplot2) })

root <- ROOT # from 00_project_setup.R (sourced by the stage runner / .Rprofile) -- no scientific change
base_dir <- file.path(root, "03_Output/03_bahia_pipeline/model_fits/v4_9_replication")
figure_dir <- file.path(root, "03_Output/03_bahia_pipeline/figures/renewal_bahia_v4_9_global_q_multisite_serology/comparison")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

on <- readRDS(file.path(base_dir, "outputs_global_q_multisite_serology/on_full2/renewal_bahia_global_q_fit_on_full2.rds"))
off <- readRDS(file.path(base_dir, "outputs_global_q_multisite_serology/off/renewal_bahia_global_q_fit_off.rds"))

q_prior <- as.vector(rstan::extract(on$fit, "q_prior_draw")$q_prior_draw)
q_on <- as.vector(rstan::extract(on$fit, "q")$q)
q_off <- as.vector(rstan::extract(off$fit, "q")$q)

df <- bind_rows(
  tibble(q = q_prior, dist = "Prior"),
  tibble(q = q_off, dist = "Posterior, WITHOUT serology (J_sero=0)"),
  tibble(q = q_on, dist = "Posterior, WITH serology (6 surveys)")
)
fixed_q_grid <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

p <- ggplot(df, aes(q, fill = dist, colour = dist)) +
  geom_density(alpha = 0.30, linewidth = 0.6) +
  geom_vline(xintercept = fixed_q_grid, linetype = 3, colour = "grey60", linewidth = 0.3) +
  annotate("text", x = fixed_q_grid, y = Inf, label = sprintf("%.2f", fixed_q_grid),
           angle = 90, vjust = 1.2, hjust = 1.1, size = 2.3, colour = "grey40") +
  scale_x_continuous(limits = c(0, 0.7)) +
  scale_fill_manual(values = c("Prior" = "grey60",
                                "Posterior, WITHOUT serology (J_sero=0)" = "#D55E00",
                                "Posterior, WITH serology (6 surveys)" = "#76558F")) +
  scale_colour_manual(values = c("Prior" = "grey40",
                                  "Posterior, WITHOUT serology (J_sero=0)" = "#A34700",
                                  "Posterior, WITH serology (6 surveys)" = "#4B2E6B")) +
  labs(title = "Bahia global-q identification ablation: does serology inform q?",
       subtitle = "Cases alone barely move q off the prior and push it ABOVE the previously explored fixed-q grid;\nserology pulls q back down substantially -- meaningful (if incomplete) identification power",
       x = "q", y = "density") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())

ggsave(file.path(figure_dir, "fig5_sero_ablation_q_comparison.png"), p, width = 190, height = 125, units = "mm", dpi = 300)
message("[saved] ", file.path(figure_dir, "fig5_sero_ablation_q_comparison.png"))
