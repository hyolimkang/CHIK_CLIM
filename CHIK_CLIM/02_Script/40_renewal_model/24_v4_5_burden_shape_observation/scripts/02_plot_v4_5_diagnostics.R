# v4.5 diagnostic figure: (A) q posterior by chain, showing whether the
# low-q/high-q split persists under the burden+shape observation model;
# (B) yearly burden posterior-predictive vs observed, by chain, showing
# the burden PPC failure of the high-q chains in 2016/2017.

required_packages <- c("rstan", "dplyr", "tibble", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(tibble); library(ggplot2); library(patchwork) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
base_dir <- file.path(root, "02_Script/40_renewal_model/24_v4_5_burden_shape_observation")
bundle <- readRDS(file.path(base_dir, "outputs/v4_5/renewal_ceara_v4_5_fit_v4_5.rds"))
fit <- bundle$fit

q_by_chain <- rstan::extract(fit, "q", permuted = FALSE)
C_total_pred_by_chain <- rstan::extract(fit, "C_total_pred", permuted = FALSE)
years <- bundle$config$years
observed_year_totals <- as.numeric(tapply(bundle$weekly_data$cases, bundle$stan_data$year_id, sum))

q_df <- bind_rows(lapply(1:4, function(c) tibble(
  chain = factor(paste0("chain ", c, if (c <= 2) " (low-q init)" else " (high-q init)")),
  q = q_by_chain[, c, 1]
)))

theme_v4 <- theme_classic(base_size = 8.5) + theme(panel.grid.major.y = element_line(colour = "grey90"))

p1 <- ggplot(q_df, aes(q, fill = chain)) +
  geom_histogram(position = "identity", alpha = .55, bins = 60) +
  labs(title = "A. v4.5 q posterior by chain — persistent bimodality", x = "q", y = "count") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())

burden_df <- bind_rows(lapply(1:4, function(c) {
  pred <- C_total_pred_by_chain[, c, ]
  tibble(
    chain = factor(paste0("chain ", c, if (c <= 2) " (low-q)" else " (high-q)")),
    year = years,
    pred_median = apply(pred, 2, median),
    pred_lo = apply(pred, 2, quantile, .025),
    pred_hi = apply(pred, 2, quantile, .975),
    observed = observed_year_totals
  )
}))

p2 <- ggplot(burden_df, aes(x = factor(year))) +
  geom_pointrange(aes(y = pred_median, ymin = pred_lo, ymax = pred_hi, colour = chain),
                   position = position_dodge(width = .6), fatten = 1.5) +
  geom_point(aes(y = observed), shape = 4, size = 2.2, colour = "black") +
  scale_y_log10(labels = scales::comma) +
  labs(title = "B. Yearly burden PPC by chain (X = observed)", x = "year", y = "total cases (log scale)") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())

figure <- p1 / p2
figure_dir <- file.path(root, "03_Output/figures/renewal_v4_5_burden_shape_observation")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(figure_dir, "v4_5_q_bimodality_and_burden_ppc.png")
ggsave(out_path, figure, width = 180, height = 200, units = "mm", dpi = 300)
message("[v4.5] figure saved: ", out_path)
