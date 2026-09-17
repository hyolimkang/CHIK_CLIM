# Pernambuco Model D -- free quadratic ridge fit + GATE (Steps 1-6).
#
# The old affine slope (b1=-0.45) is a computational coordinate constant,
# not a scientific parameter -- no requirement that a fresh quadratic fit
# reproduce it. Re-centers at u0 = median(logit_q) (NOT the old
# logit_q_reference), fits alpha_R ~ c0 + b1*d + b2*d^2 FREELY on
# non-divergent draws, and gates whether to proceed to Stan based on
# whether the curved mean relationship is actually removed.

required_packages <- c("rstan", "dplyr", "ggplot2", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(rstan); library(dplyr); library(ggplot2); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
figure_dir <- file.path(root, "03_Output/figures/pernambuco_v4_9_replication/modelD_curved_reparam")
table_dir <- file.path(root, "03_Output/tables/pernambuco_v4_9_curved_reparam")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

b <- readRDS(file.path(root, "03_Output/model_fits/pernambuco/v4_9_replication/outputs/modelC/pe_ridge_U14_full.rds"))
fit <- b$fit

get_flat <- function(par) {
  arr <- rstan::extract(fit, pars = par, permuted = FALSE, inc_warmup = FALSE)
  as.vector(arr[, , 1])
}
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
divergent <- unlist(lapply(sp, function(x) x[, "divergent__"]))
nondiv <- divergent == 0

logit_q <- get_flat("logit_q")
alpha_R <- get_flat("alpha_R")

# ---- Step 2: re-center at posterior median (non-divergent draws) ---------
u0 <- median(logit_q[nondiv])
message(sprintf("u0 = median(logit_q, non-divergent) = %.6f  (old logit_q_reference = %.6f, DIFFERENT -- not reused)", u0, b$config$logit_q_reference))

d_all <- logit_q - u0
d_nd <- d_all[nondiv]
y_nd <- alpha_R[nondiv]

# ---- Step 3: free quadratic fit -------------------------------------------
fit_quad <- lm(y_nd ~ d_nd + I(d_nd^2))
c0 <- unname(coef(fit_quad)[1]); b1 <- unname(coef(fit_quad)[2]); b2 <- unname(coef(fit_quad)[3])
rmse_quad <- sqrt(mean(residuals(fit_quad)^2))
fit_lin_check <- lm(y_nd ~ d_nd) # for RMSE comparison only
rmse_lin <- sqrt(mean(residuals(fit_lin_check)^2))

message(sprintf("\nFree quadratic fit: alpha_R = %.4f + %.4f*d + %.4f*d^2", c0, b1, b2))
message(sprintf("RMSE: linear-only=%.5f -> quadratic=%.5f (%.1f%% reduction)", rmse_lin, rmse_quad, 100*(1-rmse_quad/rmse_lin)))

# ---- Step 4: gamma_curve_diagnostic for EVERY draw (incl. divergent) -----
gamma_curve_diagnostic <- alpha_R - c0 - b1 * d_all - b2 * d_all^2

# ---- Step 5: residual correlations + binned mean/SD ------------------------
cor_d <- cor(gamma_curve_diagnostic[nondiv], d_nd)
cor_d2 <- cor(gamma_curve_diagnostic[nondiv], d_nd^2)
message(sprintf("\nResidual (non-divergent) correlation with d: %.4f | with d^2: %.4f  (both should be ~0 if curvature removed)", cor_d, cor_d2))

n_bins <- 9
bin_breaks <- quantile(logit_q[nondiv], probs = seq(0, 1, length.out = n_bins + 1))
bin_id <- cut(logit_q, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE)
binned <- tibble(logit_q = logit_q, gamma_curve_diagnostic = gamma_curve_diagnostic, bin = bin_id, divergent = divergent) |>
  dplyr::filter(!is.na(bin)) |>
  group_by(bin) |>
  summarise(logit_q_mean = mean(logit_q), n = n(), n_divergent = sum(divergent),
            mean_gamma_curve = mean(gamma_curve_diagnostic), sd_gamma_curve = sd(gamma_curve_diagnostic), .groups = "drop")
message("\n=== Binned mean/SD of gamma_curve_diagnostic across logit_q (9 bins) ===")
print(as.data.frame(binned), digits = 4)

# GATE check: is the binned mean roughly flat (no residual trend across bins)?
bin_mean_range <- max(binned$mean_gamma_curve) - min(binned$mean_gamma_curve)
overall_sd <- sd(gamma_curve_diagnostic[nondiv])
gate_pass <- abs(cor_d) < 0.15 && abs(cor_d2) < 0.15 && bin_mean_range < overall_sd
message(sprintf("\nBinned mean range = %.4f | overall residual SD = %.4f | ratio = %.3f", bin_mean_range, overall_sd, bin_mean_range/overall_sd))
message(sprintf("\n=== GATE: %s ===", if (gate_pass) "PASS -- proceed to Stan implementation" else "FAIL -- curvature not adequately removed, STOP"))

results <- tibble(quantity = c("u0_median_logit_q", "c0", "b1", "b2", "rmse_linear", "rmse_quadratic",
                                "cor_resid_vs_d", "cor_resid_vs_d2", "bin_mean_range", "overall_resid_sd", "gate_pass"),
                   value = c(u0, c0, b1, b2, rmse_lin, rmse_quad, cor_d, cor_d2, bin_mean_range, overall_sd, gate_pass))
write.csv(results, file.path(table_dir, "PE_free_quadratic_ridge_fit.csv"), row.names = FALSE)
write.csv(binned, file.path(table_dir, "PE_free_quadratic_ridge_binned.csv"), row.names = FALSE)
saveRDS(list(u0 = u0, c0 = c0, b1 = b1, b2 = b2, gate_pass = gate_pass), file.path(table_dir, "PE_free_quadratic_ridge_constants.rds"))
message("\n[saved] PE_free_quadratic_ridge_fit.csv, PE_free_quadratic_ridge_binned.csv, PE_free_quadratic_ridge_constants.rds")

# ---- Figure: gamma_curve_diagnostic vs logit_q, divergences marked --------
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))
plot_df <- tibble(logit_q = logit_q, gamma_curve_diagnostic = gamma_curve_diagnostic,
                   divergent = factor(divergent, labels = c("Non-divergent", "Divergent")))
p <- ggplot(plot_df, aes(logit_q, gamma_curve_diagnostic)) +
  geom_point(aes(colour = divergent, alpha = divergent, size = divergent)) +
  geom_point(data = binned, aes(x = logit_q_mean, y = mean_gamma_curve), inherit.aes = FALSE, colour = "#315A7D", size = 2.5, shape = 18) +
  geom_errorbar(data = binned, aes(x = logit_q_mean, ymin = mean_gamma_curve - sd_gamma_curve, ymax = mean_gamma_curve + sd_gamma_curve),
                inherit.aes = FALSE, colour = "#315A7D", width = 0.03) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  scale_colour_manual(values = c("Non-divergent" = "grey60", "Divergent" = "red")) +
  scale_alpha_manual(values = c("Non-divergent" = 0.15, "Divergent" = 0.9)) +
  scale_size_manual(values = c("Non-divergent" = 0.5, "Divergent" = 1.8)) +
  labs(title = "Pernambuco: gamma_curve_diagnostic (residual after FREE quadratic fit) vs logit_q",
       subtitle = sprintf("cor(resid,d)=%.3f, cor(resid,d^2)=%.3f | blue diamonds = binned mean +/- SD (9 bins) | GATE: %s",
                           cor_d, cor_d2, if (gate_pass) "PASS" else "FAIL"),
       x = "logit_q", y = "gamma_curve_diagnostic (residual)") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figure_dir, "PE_curved_ridge_divergence_map.png"), p, width = 220, height = 150, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "PE_curved_ridge_divergence_map.png"))
