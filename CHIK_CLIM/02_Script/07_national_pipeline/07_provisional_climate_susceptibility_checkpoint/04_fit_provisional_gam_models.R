# Provisional GAM comparison: does the current M0 (mean-preserving)
# FOI-informed susceptibility add explanatory/predictive value for early
# epidemic growth (Re_early_6) beyond climate alone?
#
# This is an EXPLORATORY CHECKPOINT (task Section 17), not a final causal
# model. Three model families (Sections 6-8), each refit across the 3
# pre-onset climate windows (Section 13) and the 3 FOI-scale episode
# datasets from script 11 (Section 14) -- 9 combinations total for the main
# comparison table. Monte Carlo uncertainty (Section 12) is handled
# separately by script 13, primary window/scale only.

required_gam_packages <- c("here", "dplyr", "readr", "tibble", "mgcv", "tidyr")
missing_gam_packages <- required_gam_packages[
  !vapply(required_gam_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_gam_packages)) stop("Missing package(s): ", paste(missing_gam_packages, collapse = ", "))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(tibble); library(mgcv); library(tidyr) })

gam_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(gam_root(), "02_Script", "07_national_pipeline", "06_national_wave_analysis", "00_national_wave_analysis_helpers.R"))

gam_settings <- list(
  climate_windows = list(
    pre6_2 = list(temp = "temp_pre6_2", precip = "precip_pre6_2", label = "pre6_2_PRIMARY"),
    pre4_1 = list(temp = "temp_pre4_1", precip = "precip_pre4_1", label = "pre4_1"),
    pre8_4 = list(temp = "temp_pre8_4", precip = "precip_pre8_4", label = "pre8_4")
  ),
  foi_scale_tags = c("0_5" = "low-FOI scenario", "1_0" = "baseline long-term FOI scenario (PRIMARY)", "2_0" = "high-FOI scenario")
)

prep_episode_data <- function(path, temp_col, precip_col) {
  d <- read_csv(path, show_col_types = FALSE)
  if (any(d$Re_early_6_median <= 0, na.rm = TRUE)) stop("Re_early_6_median <= 0 present in ", path)
  d |> transmute(
    state = factor(state), wave_id, is_recurrence,
    logRe = log(Re_early_6_median), Re_early_6_median,
    S_prop = .data[["S_prop"]],
    temp = .data[[temp_col]], log1p_precip = log1p(.data[[precip_col]])
  )
}

fit_model_set <- function(data) {
  list(
    M_climate = gam(logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + s(state, bs = "re"),
                     data = data, method = "REML"),
    M_climate_S = gam(logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + s(S_prop, k = 3) + s(state, bs = "re"),
                       data = data, method = "REML"),
    M_climate_S_linear = gam(logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + S_prop + s(state, bs = "re"),
                              data = data, method = "REML"),
    M_mechanistic_gate = gam(logRe ~ offset(log(S_prop)) + s(temp, k = 4) + s(log1p_precip, k = 4) + s(state, bs = "re"),
                              data = data, method = "REML")
  )
}

# Leave-one-STATE-out (not leave-one-wave-out): repeated waves from the same
# UF are not independent (Section 10/11), so the CV unit is the state.
# Prediction for the held-out state uses exclude = "s(state)" so the
# state-specific random effect (which by construction cannot be estimated
# for a truly unseen state) is zeroed out -- population-average prediction.
leave_one_state_out <- function(data, formula) {
  states <- levels(data$state)
  preds <- vector("list", length(states))
  for (i in seq_along(states)) {
    held_out <- states[i]
    train <- data |> dplyr::filter(state != held_out) |> droplevels()
    test <- data |> dplyr::filter(state == held_out)
    fit <- tryCatch(gam(formula, data = train, method = "REML"), error = function(e) NULL)
    if (is.null(fit)) next
    smooth_labels <- vapply(fit$smooth, function(s) s$label, character(1))
    state_smooth_label <- smooth_labels[grepl("^s\\(state\\)", smooth_labels)]
    pred <- tryCatch(
      predict(fit, newdata = test, exclude = state_smooth_label),
      error = function(e) rep(NA_real_, nrow(test))
    )
    preds[[i]] <- tibble(state = held_out, wave_id = test$wave_id, observed = test$logRe, predicted = as.numeric(pred))
  }
  bind_rows(preds)
}

model_summary_row <- function(model, model_name, data, cv_formula, climate_window, foi_scale_tag) {
  s <- summary(model)
  cv <- leave_one_state_out(data, cv_formula)
  cv <- cv |> dplyr::filter(is.finite(observed), is.finite(predicted))
  rmse <- sqrt(mean((cv$observed - cv$predicted)^2))
  mae <- mean(abs(cv$observed - cv$predicted))

  s_prop_row <- NA_character_; edf_s <- NA_real_; p_s <- NA_real_
  smooth_s_match <- if (!is.null(s$s.table)) grep("^s\\(S_prop\\)", rownames(s$s.table), value = TRUE) else character(0)
  if (length(smooth_s_match)) {
    edf_s <- s$s.table[smooth_s_match[1], "edf"]; p_s <- s$s.table[smooth_s_match[1], "p-value"]
    s_prop_row <- "smooth"
  } else if (!is.null(s$p.table) && "S_prop" %in% rownames(s$p.table)) {
    edf_s <- 1; p_s <- s$p.table["S_prop", "Pr(>|t|)"]
    s_prop_row <- if (s$p.table["S_prop", "Estimate"] > 0) "linear_positive" else "linear_negative"
  }

  tibble(
    model_name = model_name, n_waves = nrow(data), n_states = n_distinct(data$state),
    climate_window = climate_window, FOI_scale_scenario = foi_scale_tag,
    AIC = AIC(model), deviance_explained = s$dev.expl, adjusted_R2 = s$r.sq,
    leave_one_state_out_RMSE = rmse, leave_one_state_out_MAE = mae,
    direction_of_S_effect = s_prop_row, effective_df_S = edf_s,
    evidence_for_nonzero_S_effect = if (is.na(p_s)) NA_character_ else sprintf("p=%.4f", p_s)
  )
}

run_provisional_gam_comparison <- function() {
  paths <- ensure_national_output_dirs()
  all_rows <- list()
  is_recurrence_rows <- list()

  for (scale_tag in names(gam_settings$foi_scale_tags)) {
    episode_path <- file.path(paths$table, sprintf("brazil_chik_Re_S_climate_episode_analysis_scale%s.csv", scale_tag))
    if (!file.exists(episode_path)) stop("Missing episode dataset: ", episode_path, " -- run script 11 for this scale tag first.")

    for (window_name in names(gam_settings$climate_windows)) {
      w <- gam_settings$climate_windows[[window_name]]
      message("[gam] scale=", scale_tag, " window=", w$label)
      data <- prep_episode_data(episode_path, w$temp, w$precip)
      fits <- fit_model_set(data)

      all_rows[[length(all_rows) + 1L]] <- model_summary_row(
        fits$M_climate, "M_climate", data,
        logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + s(state, bs = "re"), w$label, scale_tag
      )
      all_rows[[length(all_rows) + 1L]] <- model_summary_row(
        fits$M_climate_S, "M_climate_S", data,
        logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + s(S_prop, k = 3) + s(state, bs = "re"), w$label, scale_tag
      )
      all_rows[[length(all_rows) + 1L]] <- model_summary_row(
        fits$M_climate_S_linear, "M_climate_S_linear", data,
        logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + S_prop + s(state, bs = "re"), w$label, scale_tag
      )
      all_rows[[length(all_rows) + 1L]] <- model_summary_row(
        fits$M_mechanistic_gate, "M_mechanistic_gate", data,
        logRe ~ offset(log(S_prop)) + s(temp, k = 4) + s(log1p_precip, k = 4) + s(state, bs = "re"), w$label, scale_tag
      )

      # Section 11: is_recurrence fixed-effect sensitivity, primary window/scale only.
      if (window_name == "pre6_2" && scale_tag == "1_0") {
        m_recur <- gam(logRe ~ s(temp, k = 4) + s(log1p_precip, k = 4) + s(S_prop, k = 3) + is_recurrence + s(state, bs = "re"),
                        data = data, method = "REML")
        recur_table <- summary(m_recur)$p.table
        recur_row <- recur_table[grepl("^is_recurrence", rownames(recur_table)), , drop = FALSE]
        is_recurrence_rows[[1]] <- as_tibble(recur_row, rownames = "term") |>
          rename(estimate = Estimate, std_error = `Std. Error`, statistic = `t value`, p_value = `Pr(>|t|)`) |>
          mutate(model = "M_climate_S_plus_is_recurrence")
      }
    }
  }

  comparison <- bind_rows(all_rows)
  out_path <- file.path(paths$table, "brazil_chik_provisional_gam_comparison.csv")
  write_csv(comparison, out_path)
  message("[gam] saved: ", out_path, " (", nrow(comparison), " rows)")

  if (length(is_recurrence_rows)) {
    is_recurrence_path <- file.path(paths$table, "brazil_chik_provisional_gam_is_recurrence_sensitivity.csv")
    write_csv(bind_rows(is_recurrence_rows), is_recurrence_path)
    message("[gam] saved: ", is_recurrence_path)
  }

  invisible(comparison)
}

if (sys.nframe() == 0L) run_provisional_gam_comparison()
