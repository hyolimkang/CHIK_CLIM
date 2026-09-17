# Pernambuco Model D (curved-ridge test) -- Sections 2 & 3 (revised).
#
# Uses the EXISTING affine-reparameterised FULL run (adapt_delta=0.95,
# 26 divergences). Retains the EXISTING affine ridge slope b1=-0.45 and
# reference u0=logit_q_reference (=qlogis(0.10)) exactly as already used
# to construct gamma_ridge in that Stan file:
#   gamma_ridge = alpha_R - b1*(logit_q - u0)     [already computed by Stan]
# Since gamma_ridge still shows a clear residual U-shape vs logit_q, this
# script estimates ONLY the residual quadratic curvature b2 from
#   gamma_ridge = c0 + b2*(u-u0)^2 + residual
# using NON-DIVERGENT post-warmup draws only. A diagnostic regression
# alpha_R ~ c0 + b1_check*(u-u0) + b2_check*(u-u0)^2 checks that
# b1_check is compatible with the existing fixed b1=-0.45. Geometry-fitting
# only -- b1/b2 are then passed to Stan as fixed DATA, not estimated there.

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

EXISTING_B1 <- b$config$ridge_slope           # -0.45, the fixed affine slope already used
u0 <- b$config$logit_q_reference              # qlogis(0.10), the fixed reference already used
message(sprintf("Retained existing affine constants: b1=%.6f | u0=%.6f (qlogis(0.10)=%.6f)", EXISTING_B1, u0, qlogis(0.10)))
if (!isTRUE(all.equal(u0, qlogis(0.10)))) stop("u0 does not match qlogis(0.10) -- investigate before proceeding.")

get_flat <- function(par) {
  arr <- rstan::extract(fit, pars = par, permuted = FALSE, inc_warmup = FALSE)
  as.vector(arr[, , 1])
}
sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
divergent <- unlist(lapply(sp, function(x) x[, "divergent__"]))
nondiv <- divergent == 0

logit_q <- get_flat("logit_q")
alpha_R <- get_flat("alpha_R")
gamma_ridge <- get_flat("gamma_ridge") # already = alpha_R - EXISTING_B1*(logit_q - u0), computed by Stan

d <- logit_q[nondiv] - u0
gamma_nd <- gamma_ridge[nondiv]
alpha_nd <- alpha_R[nondiv]

# ---- Primary fit: residual quadratic curvature on gamma_ridge only -------
fit_resid_quad <- lm(gamma_nd ~ I(d^2))
c0_resid <- unname(coef(fit_resid_quad)[1]); b2 <- unname(coef(fit_resid_quad)[2])
fit_resid_flat <- lm(gamma_nd ~ 1) # "linear" baseline here is just the constant (b1 already removed)

rmse_flat <- sqrt(mean(residuals(fit_resid_flat)^2))
rmse_quad <- sqrt(mean(residuals(fit_resid_quad)^2))

# ---- Diagnostic check: fit alpha_R ~ c0 + b1_check*d + b2_check*d^2 -------
fit_diag <- lm(alpha_nd ~ d + I(d^2))
b1_check <- unname(coef(fit_diag)[2]); b2_check <- unname(coef(fit_diag)[3])

message(sprintf("\nResidual quadratic fit (gamma_ridge ~ c0 + b2*d^2): c0=%.4f, b2=%.4f | RMSE flat=%.5f -> quad=%.5f (%.1f%% reduction)",
                 c0_resid, b2, rmse_flat, rmse_quad, 100*(1 - rmse_quad/rmse_flat)))
message(sprintf("Diagnostic check (alpha_R ~ c0 + b1_check*d + b2_check*d^2): b1_check=%.4f (existing fixed b1=%.4f) | b2_check=%.4f (vs residual-fit b2=%.4f)",
                 b1_check, EXISTING_B1, b2_check, b2))
message(sprintf("b1_check compatible with existing b1: %s (diff=%.4f)", abs(b1_check - EXISTING_B1) < 0.15, b1_check - EXISTING_B1))
message(sprintf("b2 consistency (two routes): diff=%.4f", b2_check - b2))

cor_resid_quad_vs_d2 <- cor(residuals(fit_resid_quad), d^2)
cor_resid_flat_vs_d2 <- cor(residuals(fit_resid_flat), d^2)
message(sprintf("\nResidual-vs-u^2 correlation: flat model=%.4f (should be large) -> quadratic model=%.4f (should be near 0)", cor_resid_flat_vs_d2, cor_resid_quad_vs_d2))

results <- tibble(
  quantity = c("u0", "existing_b1", "b2_residual_fit", "b1_check", "b2_check", "rmse_flat", "rmse_quad",
               "cor_resid_flat_vs_d2", "cor_resid_quad_vs_d2"),
  value = c(u0, EXISTING_B1, b2, b1_check, b2_check, rmse_flat, rmse_quad, cor_resid_flat_vs_d2, cor_resid_quad_vs_d2)
)
write.csv(results, file.path(table_dir, "PE_curved_ridge_regression_fit.csv"), row.names = FALSE)
message("\n[saved] ", file.path(table_dir, "PE_curved_ridge_regression_fit.csv"))
message(sprintf("\nFIXED CONSTANTS for Stan (data, NOT estimated): u0=%.6f, b1=%.6f (retained), b2=%.6f (new)", u0, EXISTING_B1, b2))
saveRDS(list(u0 = u0, b1 = EXISTING_B1, b2 = b2), file.path(table_dir, "PE_curved_ridge_constants.rds"))

# ---- Figures ---------------------------------------------------------------
theme_v4 <- theme_classic(base_size = 9) + theme(panel.grid.major.y = element_line(colour = "grey90"))

# (a) gamma_ridge vs logit_q, with flat vs quadratic residual fit
d_seq <- seq(min(d), max(d), length.out = 200)
curve_gamma <- tibble(logit_q = d_seq + u0, flat = coef(fit_resid_flat)[1], quad = c0_resid + b2 * d_seq^2)
plot_df_gamma <- tibble(logit_q = logit_q, gamma_ridge = gamma_ridge, divergent = factor(divergent, labels = c("Non-divergent", "Divergent")))
p_gamma <- ggplot(plot_df_gamma, aes(logit_q, gamma_ridge)) +
  geom_point(aes(colour = divergent, alpha = divergent, size = divergent)) +
  scale_colour_manual(values = c("Non-divergent" = "grey50", "Divergent" = "red")) +
  scale_alpha_manual(values = c("Non-divergent" = 0.15, "Divergent" = 0.9)) +
  scale_size_manual(values = c("Non-divergent" = 0.5, "Divergent" = 1.8)) +
  geom_line(data = curve_gamma, aes(x = logit_q, y = flat), inherit.aes = FALSE, colour = "#315A7D", linewidth = 0.8) +
  geom_line(data = curve_gamma, aes(x = logit_q, y = quad), inherit.aes = FALSE, colour = "#D55E00", linewidth = 0.8, linetype = 2) +
  labs(title = "Pernambuco: gamma_ridge vs logit_q -- residual U-shape after existing affine removal",
       subtitle = sprintf("Blue = flat (no curvature, RMSE=%.4f) | Orange dashed = quadratic residual fit (RMSE=%.4f, %.0f%% reduction). Red = divergent draws.",
                           rmse_flat, rmse_quad, 100*(1-rmse_quad/rmse_flat)),
       x = "logit_q", y = "gamma_ridge") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figure_dir, "PE_gamma_curve_logitq_geometry.png"), p_gamma, width = 220, height = 150, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "PE_gamma_curve_logitq_geometry.png"))

# (b) alpha_R vs logit_q, linear-only (existing) vs quadratic diagnostic fit
curve_alpha <- tibble(logit_q = d_seq + u0,
                       linear_existing = (coef(fit_resid_flat)[1]) + EXISTING_B1 * d_seq,
                       quadratic_diag = coef(fit_diag)[1] + b1_check * d_seq + b2_check * d_seq^2)
plot_df_alpha <- tibble(logit_q = logit_q, alpha_R = alpha_R, divergent = factor(divergent, labels = c("Non-divergent", "Divergent")))
p_alpha <- ggplot(plot_df_alpha, aes(logit_q, alpha_R)) +
  geom_point(aes(colour = divergent, alpha = divergent, size = divergent)) +
  scale_colour_manual(values = c("Non-divergent" = "grey50", "Divergent" = "red")) +
  scale_alpha_manual(values = c("Non-divergent" = 0.15, "Divergent" = 0.9)) +
  scale_size_manual(values = c("Non-divergent" = 0.5, "Divergent" = 1.8)) +
  geom_line(data = curve_alpha, aes(x = logit_q, y = linear_existing), inherit.aes = FALSE, colour = "#315A7D", linewidth = 0.8) +
  geom_line(data = curve_alpha, aes(x = logit_q, y = quadratic_diag), inherit.aes = FALSE, colour = "#D55E00", linewidth = 0.8, linetype = 2) +
  labs(title = "Pernambuco: alpha_R vs logit_q -- existing affine (linear) vs quadratic diagnostic fit",
       subtitle = sprintf("Blue = existing fixed affine (b1=%.3f) | Orange dashed = quadratic diagnostic (b1_check=%.3f, b2_check=%.3f). Red = divergent draws.",
                           EXISTING_B1, b1_check, b2_check),
       x = "logit_q", y = "alpha_R") +
  theme_v4 + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(figure_dir, "PE_alphaR_logitq_linear_vs_quadratic.png"), p_alpha, width = 220, height = 150, units = "mm", dpi = 300, bg = "white")
message("[saved] ", file.path(figure_dir, "PE_alphaR_logitq_linear_vs_quadratic.png"))
