## Data

pop_data_state <- readRDS("01_Data/ibge_pop_uf_age_group_2015_2024.rds")

N_mg <- pop_data_state |>
  filter(uf_code == "31")

####
theme_lancet_like <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "grey10"),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 2,
        colour = "grey10",
        margin = ggplot2::margin(b = 4)
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.5,
        colour = "grey30",
        margin = ggplot2::margin(b = 8)
      ),
      plot.caption = ggplot2::element_text(
        size = base_size - 2.5,
        colour = "grey40",
        hjust = 0,
        margin = ggplot2::margin(t = 8)
      ),
      axis.title.y = ggplot2::element_text(
        size = base_size - 0.5,
        colour = "grey15",
        margin = ggplot2::margin(r = 8)
      ),
      axis.text = ggplot2::element_text(
        size = base_size - 1.5,
        colour = "grey25"
      ),
      axis.line = ggplot2::element_line(colour = "grey25", linewidth = 0.35),
      axis.ticks = ggplot2::element_line(colour = "grey25", linewidth = 0.3),
      axis.ticks.length = grid::unit(2.5, "pt"),
      panel.grid.major.y = ggplot2::element_line(colour = "grey88", linewidth = 0.3),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.text = ggplot2::element_text(size = base_size - 1.5),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}


### test prevacc simulator 
extract_index <- function(x) {
  as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", x))
}

get_draw_vector <- function(draws_df, var, draw_id) {
  cols <- grep(paste0("^", var, "\\["), names(draws_df), value = TRUE)
  if (length(cols) == 0) {
    stop("No posterior columns found for variable: ", var)
  }
  cols <- cols[order(extract_index(cols))]
  as.numeric(draws_df[draw_id, cols, drop = TRUE])
}

get_draw_scalar <- function(draws_df, var, draw_id) {
  if (!var %in% names(draws_df)) {
    stop("No posterior column found for variable: ", var)
  }
  as.numeric(draws_df[[var]][draw_id])
}

make_pars_draw <- function(draws_df, draw_id, gamma_fixed = NULL) {
  list(
    beta_t = get_draw_vector(draws_df, "beta_t", draw_id),
    import_count_year = get_draw_vector(draws_df, "import_count_year", draw_id),
    rho = get_draw_scalar(draws_df, "rho", draw_id),
    S0_frac = get_draw_scalar(draws_df, "S0_frac", draw_id),
    I0 = get_draw_scalar(draws_df, "I0", draw_id),
    gamma = if (!is.null(gamma_fixed)) gamma_fixed else {
      if (exists("GAMMA_WEEK")) GAMMA_WEEK else 1.0
    }
  )
}

simulate_age_sir <- function(pars, ce_fit, years, scenario = "baseline",
                             target_age = 12L,
                             routine_schedule = NULL,
                             ve_sus = 0,
                             mu_death_annual = 0.0066,
                             max_age = 100L,
                             vaccinate_full_cohort = TRUE,
                             add_births = TRUE) {
  
  ages <- 0:max_age
  A <- length(ages)
  N_time <- nrow(ce_fit)
  
  pmax0 <- function(x) pmax(x, 0)
  routine_schedule <- validate_routine_schedule(routine_schedule)
  
  # Weekly natural mortality probability (matches v12f).
  p_death <- 1 - exp(-mu_death_annual / 52)
  p_rec <- 1 - exp(-pars$gamma)
  
  pop_age <- make_initial_age_population(ce_fit$population[1], max_age)
  
  S_u <- pop_age * pars$S0_frac
  S_v <- rep(0, A)
  
  I_u <- pop_age / sum(pop_age) * pars$I0
  I_v <- rep(0, A)
  
  R_u <- pmax0(pop_age - S_u - I_u)
  R_v <- rep(0, A)
  
  out <- vector("list", N_time)
  
  for (tt in seq_len(N_time)) {
    
    year_idx <- match(ce_fit$year[tt], years)
    if (is.na(year_idx)) {
      stop("Year not found in years vector: ", ce_fit$year[tt])
    }
    
    beta_t <- pars$beta_t[tt]
    
    vaccinated_this_week <- 0
    vaccinated_s_this_week <- 0
    vaccinated_i_this_week <- 0
    vaccinated_r_this_week <- 0
    
    # ------------------------------------------------------------
    # 1. Vaccination event (first epi week of a scheduled year)
    # ------------------------------------------------------------
    first_week_start_this_year <- min(ce_fit$week_start[ce_fit$year == ce_fit$year[tt]])
    
    coverage_this_year <- if (scenario == "routine") {
      coverage_for_year(routine_schedule, ce_fit$year[tt])
    } else {
      0
    }
    
    vaccination_due <- ce_fit$week_start[tt] == first_week_start_this_year &&
      coverage_this_year > 0
    
    if (vaccination_due) {
      age_idx <- which(ages %in% target_age)
      if (length(age_idx) == 0) stop("target_age is outside the modelled age range.")
      
      if (vaccinate_full_cohort) {
        vacc_s <- coverage_this_year * S_u[age_idx]
        vacc_i <- coverage_this_year * I_u[age_idx]
        vacc_r <- coverage_this_year * R_u[age_idx]
        
        S_u[age_idx] <- S_u[age_idx] - vacc_s
        I_u[age_idx] <- I_u[age_idx] - vacc_i
        R_u[age_idx] <- R_u[age_idx] - vacc_r
        
        S_v[age_idx] <- S_v[age_idx] + vacc_s
        I_v[age_idx] <- I_v[age_idx] + vacc_i
        R_v[age_idx] <- R_v[age_idx] + vacc_r
        
        vaccinated_s_this_week <- sum(vacc_s)
        vaccinated_i_this_week <- sum(vacc_i)
        vaccinated_r_this_week <- sum(vacc_r)
        vaccinated_this_week <- vaccinated_s_this_week +
          vaccinated_i_this_week + vaccinated_r_this_week
      } else {
        vacc_s <- coverage_this_year * S_u[age_idx]
        S_u[age_idx] <- S_u[age_idx] - vacc_s
        S_v[age_idx] <- S_v[age_idx] + vacc_s
        vaccinated_s_this_week <- sum(vacc_s)
        vaccinated_this_week <- sum(vacc_s)
      }
    }
    
    # ------------------------------------------------------------
    # 2. Infection and recovery dynamics
    # ------------------------------------------------------------
    N_prev <- sum(S_u) + sum(S_v) + sum(I_u) + sum(I_v) + sum(R_u) + sum(R_v)
    
    if (tt == 1L) {
      lambda_t <- 0
      new_inf_u <- rep(0, A)
      new_inf_v <- rep(0, A)
    } else {
      import_count_week <- pars$import_count_year[year_idx] / 52
      
      infectious_pressure <- (sum(I_u) + sum(I_v)) / N_prev +
        import_count_week / N_prev
      
      lambda_t <- beta_t * infectious_pressure
      
      p_inf_u <- 1 - exp(-lambda_t)
      p_inf_v <- 1 - exp(-lambda_t * (1 - ve_sus))
      
      new_inf_u <- S_u * p_inf_u
      new_inf_v <- S_v * p_inf_v
      
      rec_u <- I_u * p_rec
      rec_v <- I_v * p_rec
      
      S_u <- S_u - new_inf_u
      S_v <- S_v - new_inf_v
      
      I_u <- I_u + new_inf_u - rec_u
      I_v <- I_v + new_inf_v - rec_v
      
      R_u <- R_u + rec_u
      R_v <- R_v + rec_v
      
      # Births enter susceptible age 0.
      if (add_births) {
        S_u[1] <- S_u[1] + ce_fit$births_weekly[tt]
      }
    }
    
    # ------------------------------------------------------------
    # 3. Natural mortality (constant per-capita, all ages/compartments)
    #    Matches v12f: total -> (total + births) * (1 - p_death)
    # ------------------------------------------------------------
    deaths_this_week <- 0
    if (tt > 1L) {
      surv <- 1 - p_death
      pre_death_total <- sum(S_u) + sum(S_v) + sum(I_u) + sum(I_v) + sum(R_u) + sum(R_v)
      
      S_u <- S_u * surv
      S_v <- S_v * surv
      I_u <- I_u * surv
      I_v <- I_v * surv
      R_u <- R_u * surv
      R_v <- R_v * surv
      
      deaths_this_week <- pre_death_total * p_death
    }
    
    S_u <- pmax0(S_u); S_v <- pmax0(S_v)
    I_u <- pmax0(I_u); I_v <- pmax0(I_v)
    R_u <- pmax0(R_u); R_v <- pmax0(R_v)
    
    # ------------------------------------------------------------
    # 4. Diagnostics and output
    # ------------------------------------------------------------
    total_model_pop <- sum(S_u) + sum(S_v) + sum(I_u) + sum(I_v) + sum(R_u) + sum(R_v)
    
    effective_S <- sum(S_u) + (1 - ve_sus) * sum(S_v)
    susceptible_total <- sum(S_u) + sum(S_v)
    vaccinated_total <- sum(S_v) + sum(I_v) + sum(R_v)
    
    out[[tt]] <- tibble::tibble(
      t = tt,
      week_start = ce_fit$week_start[tt],
      calendar_year = ce_fit$year[tt],
      week_of_year = ce_fit$week_of_year[tt],
      scenario = scenario,
      
      beta_t = beta_t,
      beta_eff = beta_t * effective_S / total_model_pop,
      
      S_unvaccinated = sum(S_u),
      S_vaccinated = sum(S_v),
      I_unvaccinated = sum(I_u),
      I_vaccinated = sum(I_v),
      R_unvaccinated = sum(R_u),
      R_vaccinated = sum(R_v),
      
      total_population_modelled = total_model_pop,
      total_population_data = ce_fit$population[tt],
      pop_balance_resid = total_model_pop - ce_fit$population[tt],
      
      new_inf_unvaccinated = sum(new_inf_u),
      new_inf_vaccinated = sum(new_inf_v),
      total_infections = sum(new_inf_u) + sum(new_inf_v),
      
      deaths_this_week = deaths_this_week,
      
      susceptible_frac = susceptible_total / total_model_pop,
      effective_susceptible_frac = effective_S / total_model_pop,
      
      vaccinated_susceptible_frac = sum(S_v) / total_model_pop,
      vaccinated_total_frac = vaccinated_total / total_model_pop,
      
      vaccinated_this_week = vaccinated_this_week,
      vaccinated_s_this_week = vaccinated_s_this_week,
      vaccinated_i_this_week = vaccinated_i_this_week,
      vaccinated_r_this_week = vaccinated_r_this_week,
      vaccination_coverage = coverage_this_year,
      
      lambda = lambda_t,
      Reff = beta_t / p_rec * effective_S / total_model_pop
    )
    
    # ------------------------------------------------------------
    # 5. Aging
    # ------------------------------------------------------------
    if (tt < N_time) {
      S_u <- age_one_week(S_u)
      S_v <- age_one_week(S_v)
      I_u <- age_one_week(I_u)
      I_v <- age_one_week(I_v)
      R_u <- age_one_week(R_u)
      R_v <- age_one_week(R_v)
    }
  }
  
  dplyr::bind_rows(out)
}


simulate_prevacc_infections <- function(draws_df,
                                        ce_fit,
                                        years,
                                        n_draws = 100,
                                        draw_ids = NULL,
                                        target_age = 12L,
                                        ve_sus = 0,
                                        mu_death_annual = 0.0066,
                                        max_age = 100L,
                                        gamma_fixed = NULL,
                                        seed = 123) {
  
  T <- nrow(ce_fit)
  
  if (is.null(draw_ids)) {
    set.seed(seed)
    available_draws <- seq_len(nrow(draws_df))
    draw_ids <- sample(available_draws, size = min(n_draws, length(available_draws)))
  } else {
    n_draws <- length(draw_ids)
  }
  
  # Store draw-level weekly quantities
  weekly_infections_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  weekly_reported_mean_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  weekly_lambda_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  weekly_local_I_frac_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  weekly_effective_S_frac_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  
  sim_results_list <- vector("list", n_draws)
  
  for (i in seq_len(n_draws)) {
    
    draw_id <- draw_ids[i]
    
    pars_draw <- make_pars_draw(
      draws_df = draws_df,
      draw_id = draw_id,
      gamma_fixed = gamma_fixed
    )
    
    sim_i <- simulate_age_sir(
      pars = pars_draw,
      ce_fit = ce_fit,
      years = years,
      scenario = "baseline",
      target_age = target_age,
      routine_schedule = NULL,
      ve_sus = ve_sus,
      mu_death_annual = mu_death_annual,
      max_age = max_age
    ) |>
      dplyr::mutate(
        draw_id = draw_id,
        draw_number = i,
        expected_reported = pars_draw$rho * total_infections
      )
    
    sim_results_list[[i]] <- sim_i
    
    weekly_infections_mat[, i] <- sim_i$total_infections
    weekly_reported_mean_mat[, i] <- sim_i$expected_reported
    weekly_lambda_mat[, i] <- sim_i$lambda
    weekly_local_I_frac_mat[, i] <- sim_i$local_infectious_frac
    weekly_effective_S_frac_mat[, i] <- sim_i$effective_susceptible_frac
  }
  
  # Weekly summaries
  df_weekly_infections <- tibble::tibble(
    week = seq_len(T),
    week_start = ce_fit$week_start,
    calendar_year = ce_fit$year,
    median = apply(weekly_infections_mat, 1, median, na.rm = TRUE),
    low95 = apply(weekly_infections_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    hi95 = apply(weekly_infections_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  df_weekly_reported_mean <- tibble::tibble(
    week = seq_len(T),
    week_start = ce_fit$week_start,
    calendar_year = ce_fit$year,
    median = apply(weekly_reported_mean_mat, 1, median, na.rm = TRUE),
    low95 = apply(weekly_reported_mean_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    hi95 = apply(weekly_reported_mean_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  df_weekly_lambda <- tibble::tibble(
    week = seq_len(T),
    week_start = ce_fit$week_start,
    calendar_year = ce_fit$year,
    median = apply(weekly_lambda_mat, 1, median, na.rm = TRUE),
    low95 = apply(weekly_lambda_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    hi95 = apply(weekly_lambda_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  df_weekly_local_I_frac <- tibble::tibble(
    week = seq_len(T),
    week_start = ce_fit$week_start,
    calendar_year = ce_fit$year,
    median = apply(weekly_local_I_frac_mat, 1, median, na.rm = TRUE),
    low95 = apply(weekly_local_I_frac_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    hi95 = apply(weekly_local_I_frac_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  df_weekly_effective_S_frac <- tibble::tibble(
    week = seq_len(T),
    week_start = ce_fit$week_start,
    calendar_year = ce_fit$year,
    median = apply(weekly_effective_S_frac_mat, 1, median, na.rm = TRUE),
    low95 = apply(weekly_effective_S_frac_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    hi95 = apply(weekly_effective_S_frac_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  # Annual draw matrix
  annual_infections_draws <- purrr::map_dfr(seq_len(n_draws), function(i) {
    sim_results_list[[i]] |>
      dplyr::group_by(draw_number, draw_id, calendar_year) |>
      dplyr::summarise(
        annual_infections = sum(total_infections, na.rm = TRUE),
        annual_expected_reported = sum(expected_reported, na.rm = TRUE),
        peak_weekly_infections = max(total_infections, na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  annual_summary <- annual_infections_draws |>
    dplyr::group_by(calendar_year) |>
    dplyr::summarise(
      annual_infections_median = median(annual_infections, na.rm = TRUE),
      annual_infections_low95 = quantile(annual_infections, 0.025, na.rm = TRUE),
      annual_infections_hi95 = quantile(annual_infections, 0.975, na.rm = TRUE),
      
      annual_expected_reported_median = median(annual_expected_reported, na.rm = TRUE),
      annual_expected_reported_low95 = quantile(annual_expected_reported, 0.025, na.rm = TRUE),
      annual_expected_reported_hi95 = quantile(annual_expected_reported, 0.975, na.rm = TRUE),
      
      peak_weekly_infections_median = median(peak_weekly_infections, na.rm = TRUE),
      peak_weekly_infections_low95 = quantile(peak_weekly_infections, 0.025, na.rm = TRUE),
      peak_weekly_infections_hi95 = quantile(peak_weekly_infections, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(
    draw_ids = draw_ids,
    n_draws = n_draws,
    
    sim_results_list = sim_results_list,
    
    weekly_infections_mat = weekly_infections_mat,
    weekly_reported_mean_mat = weekly_reported_mean_mat,
    weekly_lambda_mat = weekly_lambda_mat,
    weekly_local_I_frac_mat = weekly_local_I_frac_mat,
    weekly_effective_S_frac_mat = weekly_effective_S_frac_mat,
    
    df_weekly_infections = df_weekly_infections,
    df_weekly_reported_mean = df_weekly_reported_mean,
    df_weekly_lambda = df_weekly_lambda,
    df_weekly_local_I_frac = df_weekly_local_I_frac,
    df_weekly_effective_S_frac = df_weekly_effective_S_frac,
    
    annual_infections_draws = annual_infections_draws,
    annual_summary = annual_summary
  ))
}



prevacc_test <- simulate_prevacc_infections(
  draws_df = draws_df,
  ce_fit = ce_fit,
  years = years,
  n_draws = 50,
  target_age = target_age,
  ve_sus = ve_sus,
  mu_death_annual = mu_death_annual,
  max_age = max_age,
  gamma_fixed = if (exists("GAMMA_WEEK")) GAMMA_WEEK else NULL,
  seed = 123
)



# Lancet/medical-journal style palette
col_median <- "#0B2C4D"   # deep navy
col_ribbon <- "#6BAED6"   # muted blue

p_prevacc_weekly_pub <- ggplot(
  prevacc_test$df_weekly_infections,
  aes(x = week_start, y = median)
) +
  geom_ribbon(
    aes(ymin = low95, ymax = hi95),
    fill = col_ribbon,
    alpha = 0.28,
    linewidth = 0
  ) +
  geom_line(
    colour = col_median,
    linewidth = 0.85,
    lineend = "round"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Weekly infections",
    title = "Prevaccination infection trajectory",
    subtitle = "Posterior median with 95% uncertainty interval",
    caption = "Shaded area denotes the 2.5th–97.5th percentile interval across posterior draws."
  ) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(colour = "grey10"),
    plot.title = element_text(
      face = "bold",
      size = 13,
      colour = "grey10",
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 10.5,
      colour = "grey30",
      margin = margin(b = 8)
    ),
    plot.caption = element_text(
      size = 8.5,
      colour = "grey40",
      hjust = 0,
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      size = 10.5,
      colour = "grey15",
      margin = margin(r = 8)
    ),
    axis.text = element_text(
      size = 9.5,
      colour = "grey25"
    ),
    axis.line = element_line(colour = "grey25", linewidth = 0.35),
    axis.ticks = element_line(colour = "grey25", linewidth = 0.3),
    axis.ticks.length = unit(2.5, "pt"),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

p_prevacc_weekly_pub

col_model <- "#0B2C4D"      # deep navy
col_ribbon <- "#9ECAE1"     # pale blue
col_obs <- "#7A1E28"        # muted burgundy

p_reported_weekly_pub <- ggplot() +
  geom_ribbon(
    data = prevacc_test$df_weekly_reported_mean,
    aes(x = week_start, ymin = low95, ymax = hi95),
    fill = col_ribbon,
    alpha = 0.35,
    linewidth = 0
  ) +
  geom_line(
    data = prevacc_test$df_weekly_reported_mean,
    aes(x = week_start, y = median),
    colour = col_model,
    linewidth = 0.85,
    lineend = "round"
  ) +
  geom_point(
    data = ce_fit,
    aes(x = week_start, y = cases),
    colour = col_obs,
    fill = "white",
    shape = 21,
    stroke = 0.35,
    size = 1.3,
    alpha = 0.75
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Weekly reported cases",
    title = "Observed and model-expected weekly reported cases",
    subtitle = "Posterior median and 95% uncertainty interval; points show observed surveillance data",
    caption = "Expected reported cases are calculated as reporting probability × simulated infections."
  ) +
  theme_classic(base_size = 11) +
  theme(
    text = element_text(colour = "grey10"),
    plot.title = element_text(
      face = "bold",
      size = 13,
      colour = "grey10",
      margin = margin(b = 4)
    ),
    plot.subtitle = element_text(
      size = 10.5,
      colour = "grey30",
      margin = margin(b = 8)
    ),
    plot.caption = element_text(
      size = 8.5,
      colour = "grey40",
      hjust = 0,
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      size = 10.5,
      colour = "grey15",
      margin = margin(r = 8)
    ),
    axis.text = element_text(
      size = 9.5,
      colour = "grey25"
    ),
    axis.line = element_line(colour = "grey25", linewidth = 0.35),
    axis.ticks = element_line(colour = "grey25", linewidth = 0.3),
    axis.ticks.length = unit(2.5, "pt"),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

p_reported_weekly_pub

###############################################################################
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(scales)

extract_fit_draws <- function(fit) {
  if ("CmdStanMCMC" %in% class(fit)) {
    as.data.frame(fit$draws(format = "draws_df"))
  } else {
    as.data.frame(fit)
  }
}

extract_index <- function(x) {
  as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", x))
}

draws_df_v12g <- extract_fit_draws(fit_v12e)

case_rep_cols <- grep("^cases_rep\\[", names(draws_df_v12g), value = TRUE)
case_rep_cols <- case_rep_cols[order(extract_index(case_rep_cols))]

cases_rep_mat <- as.matrix(draws_df_v12g[, case_rep_cols, drop = FALSE])

annual_cases_rep_check_v12g <- purrr::map_dfr(seq_along(years), function(i) {
  idx <- which(ce_fit$year == years[i])
  annual_rep_draw <- rowSums(cases_rep_mat[, idx, drop = FALSE], na.rm = TRUE)
  
  tibble(
    year = years[i],
    observed = sum(ce_fit$cases[idx], na.rm = TRUE),
    pred_median = median(annual_rep_draw, na.rm = TRUE),
    pred_q025 = quantile(annual_rep_draw, 0.025, na.rm = TRUE),
    pred_q975 = quantile(annual_rep_draw, 0.975, na.rm = TRUE),
    observed_inside_95 = observed >= pred_q025 & observed <= pred_q975
  )
})

print(as.data.frame(annual_cases_rep_check_v12g), digits = 5)


prevacc_state_summary <- tibble::tibble(
  week_start = ce_fit$week_start,
  calendar_year = ce_fit$year,
  
  effective_S_median = apply(prevacc_test$weekly_effective_S_frac_mat, 1, median, na.rm = TRUE),
  effective_S_low95  = apply(prevacc_test$weekly_effective_S_frac_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  effective_S_hi95   = apply(prevacc_test$weekly_effective_S_frac_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
  
  local_I_median = apply(prevacc_test$weekly_local_I_frac_mat, 1, median, na.rm = TRUE),
  local_I_low95  = apply(prevacc_test$weekly_local_I_frac_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  local_I_hi95   = apply(prevacc_test$weekly_local_I_frac_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
  
  lambda_median = apply(prevacc_test$weekly_lambda_mat, 1, median, na.rm = TRUE),
  lambda_low95  = apply(prevacc_test$weekly_lambda_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  lambda_hi95   = apply(prevacc_test$weekly_lambda_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
)

beta_cols <- grep("^beta_t\\[", names(draws_df), value = TRUE)
beta_cols <- beta_cols[order(extract_index(beta_cols))]

beta_mat_all <- as.matrix(draws_df[, beta_cols, drop = FALSE])

# Use the same posterior draws used in prevacc_test
beta_mat <- beta_mat_all[prevacc_test$draw_ids, , drop = FALSE]

beta_summary <- tibble::tibble(
  week_start = ce_fit$week_start,
  calendar_year = ce_fit$year,
  beta_median = apply(beta_mat, 2, median, na.rm = TRUE),
  beta_low95  = apply(beta_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
  beta_hi95   = apply(beta_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
)


p_rec <- 1 - exp(-(if (exists("GAMMA_WEEK")) GAMMA_WEEK else 1.0))

# Reff draw matrix: beta[t, draw] / p_rec * S_eff[t, draw]
# beta_mat is [draw x week], weekly_effective_S_frac_mat is [week x draw]
reff_mat <- t(beta_mat) / p_rec * prevacc_test$weekly_effective_S_frac_mat

reff_summary <- tibble::tibble(
  week_start = ce_fit$week_start,
  calendar_year = ce_fit$year,
  Reff_median = apply(reff_mat, 1, median, na.rm = TRUE),
  Reff_low95  = apply(reff_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
  Reff_hi95   = apply(reff_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
)

p_susceptible <- ggplot(
  prevacc_state_summary,
  aes(x = week_start, y = effective_S_median)
) +
  geom_ribbon(
    aes(ymin = effective_S_low95, ymax = effective_S_hi95),
    fill = "#9ECAE1",
    alpha = 0.35,
    linewidth = 0
  ) +
  geom_line(
    colour = "#0B2C4D",
    linewidth = 0.85,
    lineend = "round"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Susceptible fraction",
    title = "Estimated susceptible trajectory",
    subtitle = "Posterior median and 95% interval across 50 draws",
    caption = "Susceptible fraction is calculated from the prevaccination replay simulation."
  ) +
  theme_lancet_like()

p_susceptible

p_beta <- ggplot(
  beta_summary,
  aes(x = week_start, y = beta_median)
) +
  geom_ribbon(
    aes(ymin = beta_low95, ymax = beta_hi95),
    fill = "#DDA0A6",
    alpha = 0.35,
    linewidth = 0
  ) +
  geom_line(
    colour = "#7A1E28",
    linewidth = 0.85,
    lineend = "round"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = NULL,
    y = expression(beta[t]),
    title = "Estimated transmission rate",
    subtitle = "Weekly posterior trajectory of beta[t]",
    caption = "beta[t] combines the fitted baseline, year effect, and seasonal component."
  ) +
  theme_lancet_like()

p_beta


p_reff <- ggplot(
  reff_summary,
  aes(x = week_start, y = Reff_median)
) +
  geom_ribbon(
    aes(ymin = Reff_low95, ymax = Reff_hi95),
    fill = "#A1D99B",
    alpha = 0.35,
    linewidth = 0
  ) +
  geom_hline(
    yintercept = 1,
    colour = "grey35",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  geom_line(
    colour = "#3B6E22",
    linewidth = 0.85,
    lineend = "round"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = NULL,
    y = expression(R[eff]),
    title = "Effective transmission index",
    subtitle = "Reff calculated from beta[t] and the simulated susceptible fraction",
    caption = "Dashed line denotes Reff = 1."
  ) +
  theme_lancet_like()

p_reff


p_local_I <- ggplot(
  prevacc_state_summary,
  aes(x = week_start, y = local_I_median)
) +
  geom_ribbon(
    aes(
      ymin = pmax(local_I_low95, 1e-12),
      ymax = pmax(local_I_hi95, 1e-12)
    ),
    fill = "#9ECAE1",
    alpha = 0.35,
    linewidth = 0
  ) +
  geom_line(
    colour = "#0B2C4D",
    linewidth = 0.85,
    lineend = "round"
  ) +
  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_log10(
   labels = scales::label_scientific(digits = 1),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Local infectious fraction",
    title = "Estimated local infectious reservoir",
    subtitle = "Posterior median and 95% interval, shown on log scale",
    caption = "Persistent low-level infection between epidemic peaks indicates fitted local transmission continuity."
  ) +
  theme_lancet_like()

p_local_I





