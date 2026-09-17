# Brazil first-epidemic climate-transmission pilot -- Section 5: primary
# GAM. log(Re_median) ~ s(temperature) + s(precipitation) + s(UF, bs="re").
# Modest spline complexity (k conservative given ~18 independent UF
# clusters, each contributing 3 correlated within-episode observations --
# the UF random effect absorbs that clustering). This is a median-Re point
# estimate; Section 6 (script 19) repeats this over posterior Re draws.

required_packages <- c("here", "dplyr", "readr", "mgcv", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(mgcv); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/07_national_pipeline/06_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

training <- read_csv(file.path(paths$table, "first_epidemic_climate_training_data.csv"), show_col_types = FALSE) |>
  mutate(UF = factor(UF), log_Re = log(Re_median))

message(sprintf("Training data: N=%d rows, %d UFs (%.1f obs/UF)", nrow(training), n_distinct(training$UF), nrow(training) / n_distinct(training$UF)))
stopifnot(all(training$Re_median > 0))

K_SMOOTH <- 4L # conservative given ~18 independent UF clusters
gam_fit <- mgcv::gam(
  log_Re ~ s(temperature, k = K_SMOOTH) + s(precipitation, k = K_SMOOTH) + s(UF, bs = "re"),
  data = training, method = "REML"
)

message("\n=== Primary GAM summary ===")
print(summary(gam_fit))

gam_diag <- tibble(
  n_obs = nrow(training), n_ufs = n_distinct(training$UF),
  edf_temperature = summary(gam_fit)$s.table["s(temperature)", "edf"],
  edf_precipitation = summary(gam_fit)$s.table["s(precipitation)", "edf"],
  edf_UF_re = summary(gam_fit)$s.table["s(UF)", "edf"],
  p_temperature = summary(gam_fit)$s.table["s(temperature)", "p-value"],
  p_precipitation = summary(gam_fit)$s.table["s(precipitation)", "p-value"],
  deviance_explained = summary(gam_fit)$dev.expl,
  adj_r2 = summary(gam_fit)$r.sq,
  REML_score = gam_fit$gcv.ubre
)
message("\n=== GAM diagnostics ===")
print(as.data.frame(gam_diag), digits = 4)

message("\n=== Concurvity ===")
conc <- mgcv::concurvity(gam_fit, full = TRUE)
print(round(conc, 3))

message("\n=== Residual diagnostics ===")
resid_vals <- residuals(gam_fit)
message(sprintf("Residual SD=%.4f | Shapiro-Wilk p=%.4f (normality of residuals)", sd(resid_vals), shapiro.test(resid_vals)$p.value))

# ---- Prediction grid (Section 7) -- restricted to the observed climate range ----
temp_range <- range(training$temperature)
precip_range <- range(training$precipitation)
grid <- expand.grid(
  temperature = seq(temp_range[1], temp_range[2], length.out = 60),
  precipitation = seq(precip_range[1], precip_range[2], length.out = 60)
)
# Predict with UF averaged out (exclude the random-effect term) -- this is
# the population-level (average-UF) climate-response surface.
pred <- predict(gam_fit, newdata = grid, se.fit = TRUE, exclude = "s(UF)", newdata.guaranteed = TRUE)
fit_vec <- as.numeric(pred$fit); se_vec <- as.numeric(pred$se.fit)
grid$R0_climate_hat <- exp(fit_vec)
grid$R0_climate_hat_lo <- exp(fit_vec - 1.96 * se_vec)
grid$R0_climate_hat_hi <- exp(fit_vec + 1.96 * se_vec)

message(sprintf("\nPrediction grid: temperature [%.2f, %.2f] deg C x precipitation [%.1f, %.1f] mm (observed-data range only)",
                 temp_range[1], temp_range[2], precip_range[1], precip_range[2]))

# ---- Save objects ----
saveRDS(gam_fit, file.path(paths$table, "first_epidemic_climate_gam.rds"))
write_csv(grid, file.path(paths$table, "climate_R0_prediction_grid.csv"))
write_csv(gam_diag, file.path(paths$table, "first_epidemic_climate_model_diagnostics.csv"))
message("\n[saved] ", file.path(paths$table, "first_epidemic_climate_gam.rds"))
message("[saved] ", file.path(paths$table, "climate_R0_prediction_grid.csv"))
message("[saved] ", file.path(paths$table, "first_epidemic_climate_model_diagnostics.csv"))
