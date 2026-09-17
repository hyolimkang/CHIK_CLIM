# Brazil first-epidemic climate-transmission pilot -- Section 6: propagate
# Re estimation uncertainty by repeating the primary GAM over >=200
# posterior draws of early-wave Re (reusing the EXACT already-saved
# early-Re posterior draws -- no new Re model fitting). For each MC
# iteration, one Re draw is resampled per (state, wave_id, window_weeks)
# training row, the GAM is refit, and the climate-response grid/curves are
# predicted; results are combined into median + 95% interval bands.

required_packages <- c("here", "dplyr", "readr", "mgcv", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(here); library(dplyr); library(readr); library(mgcv); library(tibble) })

root <- "C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Documents/GitHub/CHIK_CLIM/CHIK_CLIM"
source(file.path(root, "02_Script/07_national_pipeline/06_national_wave_analysis/00_national_wave_analysis_helpers.R"))
paths <- ensure_national_output_dirs()

training <- read_csv(file.path(paths$table, "first_epidemic_climate_training_data.csv"), show_col_types = FALSE) |>
  mutate(UF = factor(UF))
re_draws <- readRDS(file.path(paths$table, "brazil_chik_major_wave_early_re_posterior_draws.rds")) |>
  rename(UF = state) |>
  semi_join(training |> distinct(UF, episode_id, window_weeks), by = c("UF", "wave_id" = "episode_id", "window_weeks"))

N_MC <- 200L
K_SMOOTH <- 4L
set.seed(20260916L)

temp_range <- range(training$temperature)
precip_range <- range(training$precipitation)
grid <- expand.grid(
  temperature = seq(temp_range[1], temp_range[2], length.out = 40),
  precipitation = seq(precip_range[1], precip_range[2], length.out = 40)
)

message(sprintf("Running %d Monte Carlo iterations over Re posterior draws...", N_MC))
mc_preds <- matrix(NA_real_, nrow = nrow(grid), ncol = N_MC)
mc_edf_temp <- numeric(N_MC); mc_edf_precip <- numeric(N_MC); mc_dev_expl <- numeric(N_MC)

max_draw <- max(re_draws$draw)
for (m in seq_len(N_MC)) {
  draw_id <- sample.int(max_draw, 1L)
  one_draw <- re_draws |> dplyr::filter(draw == draw_id) |>
    rename(episode_id = wave_id)
  td <- training |> select(-Re_median, -Re_lower, -Re_upper) |>
    inner_join(one_draw |> select(UF, episode_id, window_weeks, Re), by = c("UF", "episode_id", "window_weeks")) |>
    mutate(UF = factor(UF, levels = levels(training$UF)), log_Re = log(Re))
  fit_m <- tryCatch(
    mgcv::gam(log_Re ~ s(temperature, k = K_SMOOTH) + s(precipitation, k = K_SMOOTH) + s(UF, bs = "re"), data = td, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit_m)) next
  pred_m <- predict(fit_m, newdata = grid, exclude = "s(UF)", newdata.guaranteed = TRUE)
  mc_preds[, m] <- exp(as.numeric(pred_m))
  sm <- summary(fit_m)
  mc_edf_temp[m] <- sm$s.table["s(temperature)", "edf"]
  mc_edf_precip[m] <- sm$s.table["s(precipitation)", "edf"]
  mc_dev_expl[m] <- sm$dev.expl
}
n_success <- sum(!is.na(mc_preds[1, ]))
message(sprintf("Completed: %d / %d successful MC iterations", n_success, N_MC))

grid$R0_climate_hat_mc_median <- apply(mc_preds, 1, median, na.rm = TRUE)
grid$R0_climate_hat_mc_lo95 <- apply(mc_preds, 1, quantile, .025, na.rm = TRUE)
grid$R0_climate_hat_mc_hi95 <- apply(mc_preds, 1, quantile, .975, na.rm = TRUE)

mc_summary <- tibble(
  n_mc_iterations = N_MC, n_successful = n_success,
  edf_temperature_median = median(mc_edf_temp, na.rm = TRUE), edf_temperature_lo95 = quantile(mc_edf_temp, .025, na.rm = TRUE), edf_temperature_hi95 = quantile(mc_edf_temp, .975, na.rm = TRUE),
  edf_precipitation_median = median(mc_edf_precip, na.rm = TRUE), edf_precipitation_lo95 = quantile(mc_edf_precip, .025, na.rm = TRUE), edf_precipitation_hi95 = quantile(mc_edf_precip, .975, na.rm = TRUE),
  dev_explained_median = median(mc_dev_expl, na.rm = TRUE)
)
message("\n=== Monte Carlo summary ===")
print(as.data.frame(mc_summary), digits = 4)

write_csv(grid, file.path(paths$table, "climate_R0_prediction_grid_mc.csv"))
write_csv(mc_summary, file.path(paths$table, "first_epidemic_climate_gam_mc_summary.csv"))
message("\n[saved] ", file.path(paths$table, "climate_R0_prediction_grid_mc.csv"))
message("[saved] ", file.path(paths$table, "first_epidemic_climate_gam_mc_summary.csv"))
