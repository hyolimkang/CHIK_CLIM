## Data

pop_data_state <- readRDS("01_Data/ibge_pop_uf_age_group_split12_2015_2024.rds")

N_mg <- pop_data_state |>
  filter(uf_code == "31")

age_gr_levels = c(
  "<1", "1-4", "5–9", "10-11", "12", "13-17", "18–19", "20–24", "25–29",
  "30–34", "35–39", "40–44", "45–49", "50–54", "55–59",
  "60–64", "65–69", "70–74", "75–79", "80–84", "85+"
)

ve_sus = 0.989 

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

## expand
make_age_pop_weekly <- function(N_mg, ce_fit, age_gr_levels) {
  
  normalise_age <- function(x) {
    x |>
      as.character() |>
      gsub("–", "-", x = _) |>
      gsub("—", "-", x = _) |>
      trimws()
  }
  
  age_levels_norm <- normalise_age(age_gr_levels)
  
  pop_df <- N_mg |>
    dplyr::mutate(age_group_norm = normalise_age(age_group))
  
  missing_ages <- setdiff(age_levels_norm, unique(pop_df$age_group_norm))
  
  if (length(missing_ages) > 0) {
    stop(
      "The following age groups are missing in N_mg: ",
      paste(missing_ages, collapse = ", ")
    )
  }
  
  extra_ages <- setdiff(unique(pop_df$age_group_norm), age_levels_norm)
  
  if (length(extra_ages) > 0) {
    warning(
      "The following age groups exist in N_mg but are not used in age_gr_levels: ",
      paste(extra_ages, collapse = ", ")
    )
  }
  
  missing_years <- setdiff(unique(ce_fit$year), unique(pop_df$year))
  
  if (length(missing_years) > 0) {
    stop(
      "N_mg is missing years required by ce_fit: ",
      paste(missing_years, collapse = ", ")
    )
  }
  
  T <- nrow(ce_fit)
  A <- length(age_gr_levels)
  
  age_pop_weekly <- matrix(NA_real_, nrow = A, ncol = T)
  
  for (tt in seq_len(T)) {
    
    yy <- ce_fit$year[tt]
    
    pop_year <- pop_df |>
      dplyr::filter(year == yy) |>
      dplyr::mutate(
        age_group_norm = factor(age_group_norm, levels = age_levels_norm)
      ) |>
      dplyr::arrange(age_group_norm)
    
    if (nrow(pop_year) != A) {
      stop(
        "Population data for year ", yy,
        " has ", nrow(pop_year),
        " rows, but expected ", A, " age groups."
      )
    }
    
    if (any(is.na(pop_year$age_group_norm))) {
      stop("Some age groups could not be matched after normalisation.")
    }
    
    age_pop_weekly[, tt] <- pop_year$population
  }
  
  rownames(age_pop_weekly) <- age_gr_levels
  colnames(age_pop_weekly) <- seq_len(T)
  
  age_pop_weekly
}

make_ageing_prob_annual <- function(N_mg, age_gr_levels) {
  
  normalise_age <- function(x) {
    x |>
      as.character() |>
      gsub("–", "-", x = _) |>
      gsub("—", "-", x = _)
  }
  
  age_meta <- N_mg |>
    dplyr::mutate(age_group_norm = normalise_age(age_group)) |>
    dplyr::distinct(age_group_norm, age_lower, age_upper) |>
    dplyr::mutate(
      age_group_norm = factor(
        age_group_norm,
        levels = normalise_age(age_gr_levels)
      )
    ) |>
    dplyr::arrange(age_group_norm)
  
  if (nrow(age_meta) != length(age_gr_levels)) {
    stop("age_meta does not match age_gr_levels.")
  }
  
  age_width <- age_meta$age_upper - age_meta$age_lower + 1
  age_width[length(age_width)] <- Inf
  
  p <- 1 / age_width
  p[!is.finite(p)] <- 0
  
  p
}


apply_ageing <- function(x, ageing_prob) {
  
  A <- length(x)
  
  out <- x * ageing_prob
  out[A] <- 0
  
  x_new <- x - out
  x_new[2:A] <- x_new[2:A] + out[1:(A - 1)]
  
  pmax(x_new, 0)
}


apply_annual_cohort_update <- function(S_prev,
                                       I_prev,
                                       R_prev,
                                       V_prev,
                                       target_N,
                                       ageing_prob_annual,
                                       anchor_V = FALSE) {
  
  pmax0 <- function(x) pmax(x, 0)
  
  S_age <- apply_ageing(S_prev, ageing_prob_annual)
  I_age <- apply_ageing(I_prev, ageing_prob_annual)
  R_age <- apply_ageing(R_prev, ageing_prob_annual)
  V_age <- apply_ageing(V_prev, ageing_prob_annual)
  
  prev_N_age <- S_age + I_age + R_age + V_age
  
  if (anchor_V) {
    
    scale_age <- ifelse(prev_N_age > 0, target_N / prev_N_age, 0)
    
    S_new <- S_age * scale_age
    I_new <- I_age * scale_age
    R_new <- R_age * scale_age
    V_new <- V_age * scale_age
    
  } else {
    
    V_keep <- pmin(V_age, target_N)
    
    target_nonV <- pmax0(target_N - V_keep)
    prev_nonV <- S_age + I_age + R_age
    
    scale_nonV <- ifelse(prev_nonV > 0, target_nonV / prev_nonV, 0)
    
    S_new <- S_age * scale_nonV
    I_new <- I_age * scale_nonV
    R_new <- R_age * scale_nonV
    V_new <- V_keep
    
    empty_nonV <- prev_nonV <= 0 & target_nonV > 0
    
    if (any(empty_nonV)) {
      S_new[empty_nonV] <- target_nonV[empty_nonV]
      I_new[empty_nonV] <- 0
      R_new[empty_nonV] <- 0
    }
  }
  
  list(
    S = pmax0(S_new),
    I = pmax0(I_new),
    R = pmax0(R_new),
    V = pmax0(V_new)
  )
}

simulate_age_sir_20groups_popdyn <- function(pars,
                                             ce_fit,
                                             years,
                                             age_pop_weekly,
                                             age_gr_levels = NULL) {
  
  pmax0 <- function(x) pmax(x, 0)
  
  T <- nrow(ce_fit)
  A <- nrow(age_pop_weekly)
  
  if (ncol(age_pop_weekly) != T) {
    stop("age_pop_weekly must have T columns matching nrow(ce_fit).")
  }
  
  p_rec <- 1 - exp(-pars$gamma)
  
  S <- matrix(0, nrow = A, ncol = T)
  I <- matrix(0, nrow = A, ncol = T)
  R <- matrix(0, nrow = A, ncol = T)
  
  age_new_inf <- matrix(0, nrow = A, ncol = T)
  age_recov <- matrix(0, nrow = A, ncol = T)
  
  lambda <- numeric(T)
  local_infectious_frac <- numeric(T)
  import_frac <- numeric(T)
  total_population_modelled <- numeric(T)
  
  # Initial age population
  N0 <- age_pop_weekly[, 1]
  
  S[, 1] <- N0 * pars$S0_frac
  I[, 1] <- N0 / sum(N0) * pars$I0
  R[, 1] <- pmax0(N0 - S[, 1] - I[, 1])
  
  total_population_modelled[1] <- sum(S[, 1]) + sum(I[, 1]) + sum(R[, 1])
  local_infectious_frac[1] <- sum(I[, 1]) / total_population_modelled[1]
  import_frac[1] <- 0
  lambda[1] <- 0
  
  for (tt in 2:T) {
    
    year_idx <- match(ce_fit$year[tt], years)
    if (is.na(year_idx)) {
      stop("Year not found in years vector: ", ce_fit$year[tt])
    }
    
    S_prev <- S[, tt - 1]
    I_prev <- I[, tt - 1]
    R_prev <- R[, tt - 1]
    
    # ------------------------------------------------------------------
    # Demographic update at first week of a new calendar year
    # Preserve compartment proportions within each age group.
    # ------------------------------------------------------------------
    if (ce_fit$year[tt] != ce_fit$year[tt - 1]) {
      target_N <- age_pop_weekly[, tt]
      prev_N_age <- S_prev + I_prev + R_prev
      
      scale_age <- ifelse(prev_N_age > 0, target_N / prev_N_age, 0)
      
      S_prev <- S_prev * scale_age
      I_prev <- I_prev * scale_age
      R_prev <- R_prev * scale_age
    }
    
    N_prev <- sum(S_prev) + sum(I_prev) + sum(R_prev)
    
    import_count_week <- pars$import_count_year[year_idx] / 52
    
    local_infectious_frac[tt] <- sum(I_prev) / N_prev
    import_frac[tt] <- import_count_week / N_prev
    
    infectious_pressure <- local_infectious_frac[tt] + import_frac[tt]
    lambda[tt] <- pars$beta_t[tt] * infectious_pressure
    
    p_inf <- 1 - exp(-lambda[tt])
    
    new_inf <- S_prev * p_inf
    recov <- I_prev * p_rec
    
    S_new <- S_prev - new_inf
    I_new <- I_prev + new_inf - recov
    R_new <- R_prev + recov
    
    S[, tt] <- pmax0(S_new)
    I[, tt] <- pmax0(I_new)
    R[, tt] <- pmax0(R_new)
    
    age_new_inf[, tt] <- new_inf
    age_recov[, tt] <- recov
    
    total_population_modelled[tt] <- sum(S[, tt]) + sum(I[, tt]) + sum(R[, tt])
  }
  
  if (!is.null(age_gr_levels)) {
    rownames(S) <- age_gr_levels
    rownames(I) <- age_gr_levels
    rownames(R) <- age_gr_levels
    rownames(age_new_inf) <- age_gr_levels
    rownames(age_recov) <- age_gr_levels
  }
  
  weekly_df <- tibble::tibble(
    week = seq_len(T),
    week_start = ce_fit$week_start,
    calendar_year = ce_fit$year,
    total_infections = colSums(age_new_inf),
    total_S = colSums(S),
    total_I = colSums(I),
    total_R = colSums(R),
    total_population_modelled = total_population_modelled,
    effective_susceptible_frac = colSums(S) / total_population_modelled,
    local_infectious_frac = local_infectious_frac,
    import_frac = import_frac,
    lambda = lambda
  )
  
  list(
    weekly_df = weekly_df,
    age_new_inf = age_new_inf,
    age_S = S,
    age_I = I,
    age_R = R,
    age_recov = age_recov
  )
}

simulate_prevacc_infections_age20 <- function(draws_df,
                                              ce_fit,
                                              years,
                                              N_mg,
                                              age_gr_levels,
                                              n_draws = 100,
                                              draw_ids = NULL,
                                              gamma_fixed = NULL,
                                              seed = 123) {
  
  T <- nrow(ce_fit)
  A <- length(age_gr_levels)
  
  age_pop_weekly <- make_age_pop_weekly(
    N_mg = N_mg,
    ce_fit = ce_fit,
    age_gr_levels = age_gr_levels
  )
  
  if (is.null(draw_ids)) {
    set.seed(seed)
    available_draws <- seq_len(nrow(draws_df))
    draw_ids <- sample(available_draws, size = min(n_draws, length(available_draws)))
  } else {
    n_draws <- length(draw_ids)
  }
  
  age_infections_array <- array(NA_real_, dim = c(A, T, n_draws))
  age_S_array <- array(NA_real_, dim = c(A, T, n_draws))
  age_I_array <- array(NA_real_, dim = c(A, T, n_draws))
  age_R_array <- array(NA_real_, dim = c(A, T, n_draws))
  
  weekly_infections_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
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
    
    sim_i <- simulate_age_sir_20groups_popdyn(
      pars = pars_draw,
      ce_fit = ce_fit,
      years = years,
      age_pop_weekly = age_pop_weekly,
      age_gr_levels = age_gr_levels
    )
    
    sim_results_list[[i]] <- sim_i
    
    age_infections_array[, , i] <- sim_i$age_new_inf
    age_S_array[, , i] <- sim_i$age_S
    age_I_array[, , i] <- sim_i$age_I
    age_R_array[, , i] <- sim_i$age_R
    
    weekly_infections_mat[, i] <- sim_i$weekly_df$total_infections
    weekly_lambda_mat[, i] <- sim_i$weekly_df$lambda
    weekly_local_I_frac_mat[, i] <- sim_i$weekly_df$local_infectious_frac
    weekly_effective_S_frac_mat[, i] <- sim_i$weekly_df$effective_susceptible_frac
  }
  
  dimnames(age_infections_array) <- list(
    age_group = age_gr_levels,
    week = seq_len(T),
    draw = seq_len(n_draws)
  )
  
  dimnames(age_S_array) <- dimnames(age_infections_array)
  dimnames(age_I_array) <- dimnames(age_infections_array)
  dimnames(age_R_array) <- dimnames(age_infections_array)
  
  weekly_summary <- tibble::tibble(
    week = seq_len(T),
    week_start = ce_fit$week_start,
    calendar_year = ce_fit$year,
    infections_median = apply(weekly_infections_mat, 1, median, na.rm = TRUE),
    infections_low95 = apply(weekly_infections_mat, 1, quantile, probs = 0.025, na.rm = TRUE),
    infections_hi95 = apply(weekly_infections_mat, 1, quantile, probs = 0.975, na.rm = TRUE),
    lambda_median = apply(weekly_lambda_mat, 1, median, na.rm = TRUE),
    local_I_median = apply(weekly_local_I_frac_mat, 1, median, na.rm = TRUE),
    effective_S_median = apply(weekly_effective_S_frac_mat, 1, median, na.rm = TRUE)
  )
  
  # Age × week summary
  age_week_summary <- purrr::map_dfr(seq_len(A), function(a) {
    purrr::map_dfr(seq_len(T), function(tt) {
      x <- age_infections_array[a, tt, ]
      tibble::tibble(
        age_group = age_gr_levels[a],
        age_index = a,
        week = tt,
        week_start = ce_fit$week_start[tt],
        calendar_year = ce_fit$year[tt],
        infections_median = median(x, na.rm = TRUE),
        infections_low95 = quantile(x, 0.025, na.rm = TRUE),
        infections_hi95 = quantile(x, 0.975, na.rm = TRUE)
      )
    })
  })
  
  annual_age_summary <- age_week_summary |>
    dplyr::group_by(age_group, age_index, calendar_year) |>
    dplyr::summarise(
      annual_infections_median = sum(infections_median, na.rm = TRUE),
      annual_infections_low95 = sum(infections_low95, na.rm = TRUE),
      annual_infections_hi95 = sum(infections_hi95, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    draw_ids = draw_ids,
    n_draws = n_draws,
    age_pop_weekly = age_pop_weekly,
    sim_results_list = sim_results_list,
    
    age_infections_array = age_infections_array,
    age_S_array = age_S_array,
    age_I_array = age_I_array,
    age_R_array = age_R_array,
    
    weekly_infections_mat = weekly_infections_mat,
    weekly_lambda_mat = weekly_lambda_mat,
    weekly_local_I_frac_mat = weekly_local_I_frac_mat,
    weekly_effective_S_frac_mat = weekly_effective_S_frac_mat,
    
    weekly_summary = weekly_summary,
    age_week_summary = age_week_summary,
    annual_age_summary = annual_age_summary
  )
}



simulate_age_sir_20groups_vacc <- function(pars,
                                           ce_fit,
                                           years,
                                           age_pop_weekly,
                                           age_gr_levels,
                                           target_age_index,
                                           routine_schedule,
                                           VE_inf = 0,
                                           ve_sus = NULL,
                                           scenario = "baseline",
                                           weekly_delivery_speed = 0.10,
                                           delay_weeks = 2L,
                                           ageing_prob_weekly = NULL,
                                           ageing_prob_annual = NULL,
                                           use_ageing = TRUE,
                                           ageing_mode = "annual",
                                           anchor_V = FALSE) {
  
  pmax0 <- function(x) pmax(x, 0)
  
  T <- nrow(ce_fit)
  A <- nrow(age_pop_weekly)
  
  ageing_mode <- match.arg(ageing_mode, c("annual", "weekly", "none"))
  
  if (!use_ageing) {
    ageing_mode <- "none"
  }
  
  if (!is.null(ve_sus)) {
    VE_inf <- ve_sus
  }
  
  if (length(target_age_index) < 1) {
    stop("target_age_index must contain at least one age-group index.")
  }
  
  if (any(target_age_index < 1 | target_age_index > A)) {
    stop("target_age_index contains invalid age-group index.")
  }
  
  if (is.null(ageing_prob_weekly)) {
    ageing_prob_weekly <- rep(0, A)
  }
  
  if (is.null(ageing_prob_annual)) {
    ageing_prob_annual <- rep(0, A)
  }
  
  if (length(ageing_prob_weekly) != A) {
    stop("length(ageing_prob_weekly) must equal number of age groups.")
  }
  
  if (length(ageing_prob_annual) != A) {
    stop("length(ageing_prob_annual) must equal number of age groups.")
  }
  
  if (weekly_delivery_speed <= 0 || weekly_delivery_speed > 1) {
    stop("weekly_delivery_speed must be > 0 and <= 1.")
  }
  
  delay_weeks <- as.integer(delay_weeks)
  if (delay_weeks < 0) {
    stop("delay_weeks must be >= 0.")
  }
  
  p_rec <- 1 - exp(-pars$gamma)
  
  S <- matrix(0, nrow = A, ncol = T)
  I <- matrix(0, nrow = A, ncol = T)
  R <- matrix(0, nrow = A, ncol = T)
  V <- matrix(0, nrow = A, ncol = T)
  
  age_new_inf <- matrix(0, nrow = A, ncol = T)
  age_recov <- matrix(0, nrow = A, ncol = T)
  
  raw_allocation_age <- matrix(0, nrow = A, ncol = T)
  vacc_delayed <- matrix(0, nrow = A, ncol = T)
  wasted_dose <- matrix(0, nrow = A, ncol = T)
  immunized_this_week <- matrix(0, nrow = A, ncol = T)
  
  age_vaccinated_this_week <- matrix(0, nrow = A, ncol = T)
  age_protected_this_week <- matrix(0, nrow = A, ncol = T)
  
  coverage_target_raw <- numeric(T)
  coverage_target_susceptible <- numeric(T)
  
  lambda <- numeric(T)
  local_infectious_frac <- numeric(T)
  import_frac <- numeric(T)
  total_population_modelled <- numeric(T)
  protected_frac <- numeric(T)
  effective_susceptible_frac <- numeric(T)
  
  annual_supply_age <- rep(0, A)
  annual_used_age <- rep(0, A)
  annual_used_susceptible_age <- rep(0, A)
  annual_target_pop <- 0
  annual_total_supply <- 0
  current_campaign_year <- NA_integer_
  
  vacc_start_week <- rep(NA_integer_, A)
  vacc_end_week <- rep(NA_integer_, A)
  
  campaign_log <- list()
  
  N0 <- age_pop_weekly[, 1]
  
  S[, 1] <- N0 * pars$S0_frac
  I[, 1] <- N0 / sum(N0) * pars$I0
  R[, 1] <- pmax0(N0 - S[, 1] - I[, 1])
  V[, 1] <- 0
  
  setup_annual_campaign <- function(tt) {
    
    yy <- ce_fit$year[tt]
    
    annual_supply_age <<- rep(0, A)
    annual_used_age <<- rep(0, A)
    annual_used_susceptible_age <<- rep(0, A)
    annual_target_pop <<- 0
    annual_total_supply <<- 0
    current_campaign_year <<- NA_integer_
    
    if (scenario != "routine") return(invisible(NULL))
    
    matched <- routine_schedule$calendar_year == yy
    if (!any(matched)) return(invisible(NULL))
    
    coverage_this_year <- routine_schedule$coverage[matched][1]
    if (is.na(coverage_this_year) || coverage_this_year <= 0) return(invisible(NULL))
    
    target_pop_year <- sum(age_pop_weekly[target_age_index, tt], na.rm = TRUE)
    if (target_pop_year <= 0) return(invisible(NULL))
    
    annual_total_supply <<- target_pop_year * coverage_this_year
    annual_target_pop <<- target_pop_year
    current_campaign_year <<- yy
    
    annual_supply_age[target_age_index] <<-
      annual_total_supply *
      age_pop_weekly[target_age_index, tt] /
      annual_target_pop
    
    campaign_log[[as.character(yy)]] <<- tibble::tibble(
      calendar_year = yy,
      coverage = coverage_this_year,
      target_population = target_pop_year,
      annual_total_supply = annual_total_supply
    )
    
    invisible(NULL)
  }
  
  allocate_weekly_doses <- function(tt, S_state, I_state, R_state, V_state) {
    
    if (scenario != "routine") return(invisible(NULL))
    if (is.na(current_campaign_year)) return(invisible(NULL))
    if (ce_fit$year[tt] != current_campaign_year) return(invisible(NULL))
    if (annual_total_supply <= 0 || annual_target_pop <= 0) return(invisible(NULL))
    
    if (sum(annual_used_age[target_age_index], na.rm = TRUE) >= annual_total_supply) {
      return(invisible(NULL))
    }
    
    weekly_dose_total <- annual_total_supply * weekly_delivery_speed
    rem <- weekly_dose_total
    
    for (a in target_age_index) {
      
      if (rem <= 0) next
      
      remaining_supply_a <- annual_supply_age[a] - annual_used_age[a]
      if (remaining_supply_a <= 0) next
      
      proposed_alloc_a <- ceiling(
        weekly_dose_total * age_pop_weekly[a, tt] / annual_target_pop
      )
      
      alloc <- min(
        proposed_alloc_a,
        rem,
        remaining_supply_a
      )
      
      if (alloc <= 0) next
      
      N_a_state <- S_state[a] + I_state[a] + R_state[a] + V_state[a]
      prop_S <- ifelse(N_a_state > 0, S_state[a] / N_a_state, 0)
      
      vacc_to_S <- alloc * prop_S
      wasted <- alloc - vacc_to_S
      
      raw_allocation_age[a, tt] <<- alloc
      vacc_delayed[a, tt] <<- vacc_to_S
      wasted_dose[a, tt] <<- wasted
      
      annual_used_age[a] <<- annual_used_age[a] + alloc
      annual_used_susceptible_age[a] <<- annual_used_susceptible_age[a] + vacc_to_S
      
      rem <- rem - alloc
      
      if (is.na(vacc_start_week[a])) vacc_start_week[a] <<- tt
      vacc_end_week[a] <<- tt
    }
    
    coverage_target_raw[tt] <<-
      sum(annual_used_age[target_age_index], na.rm = TRUE) / annual_target_pop
    
    coverage_target_susceptible[tt] <<-
      sum(annual_used_susceptible_age[target_age_index], na.rm = TRUE) / annual_target_pop
    
    invisible(NULL)
  }
  
  setup_annual_campaign(tt = 1)
  
  allocate_weekly_doses(
    tt = 1,
    S_state = S[, 1],
    I_state = I[, 1],
    R_state = R[, 1],
    V_state = V[, 1]
  )
  
  total_population_modelled[1] <- sum(S[, 1]) + sum(I[, 1]) + sum(R[, 1]) + sum(V[, 1])
  local_infectious_frac[1] <- sum(I[, 1]) / total_population_modelled[1]
  import_frac[1] <- 0
  lambda[1] <- 0
  protected_frac[1] <- sum(V[, 1]) / total_population_modelled[1]
  effective_susceptible_frac[1] <- sum(S[, 1]) / total_population_modelled[1]
  
  for (tt in 2:T) {
    
    year_idx <- match(ce_fit$year[tt], years)
    if (is.na(year_idx)) stop("Year not found in years vector: ", ce_fit$year[tt])
    
    S_prev <- S[, tt - 1]
    I_prev <- I[, tt - 1]
    R_prev <- R[, tt - 1]
    V_prev <- V[, tt - 1]
    
    first_week_this_year <- min(ce_fit$week_start[ce_fit$year == ce_fit$year[tt]])
    is_first_week_this_year <- ce_fit$week_start[tt] == first_week_this_year
    
    if (is_first_week_this_year && ce_fit$year[tt] != ce_fit$year[tt - 1]) {
      
      target_N <- age_pop_weekly[, tt]
      
      if (ageing_mode == "annual") {
        
        annual_updated <- apply_annual_cohort_update(
          S_prev = S_prev,
          I_prev = I_prev,
          R_prev = R_prev,
          V_prev = V_prev,
          target_N = target_N,
          ageing_prob_annual = ageing_prob_annual,
          anchor_V = anchor_V
        )
        
        S_prev <- annual_updated$S
        I_prev <- annual_updated$I
        R_prev <- annual_updated$R
        V_prev <- annual_updated$V
        
      } else {
        
        prev_N_age <- S_prev + I_prev + R_prev + V_prev
        
        if (anchor_V) {
          
          scale_age <- ifelse(prev_N_age > 0, target_N / prev_N_age, 0)
          
          S_prev <- S_prev * scale_age
          I_prev <- I_prev * scale_age
          R_prev <- R_prev * scale_age
          V_prev <- V_prev * scale_age
          
        } else {
          
          V_keep <- pmin(V_prev, target_N)
          
          target_nonV <- pmax0(target_N - V_keep)
          prev_nonV <- S_prev + I_prev + R_prev
          
          scale_nonV <- ifelse(prev_nonV > 0, target_nonV / prev_nonV, 0)
          
          S_prev <- S_prev * scale_nonV
          I_prev <- I_prev * scale_nonV
          R_prev <- R_prev * scale_nonV
          V_prev <- V_keep
          
          empty_nonV <- prev_nonV <= 0 & target_nonV > 0
          if (any(empty_nonV)) {
            S_prev[empty_nonV] <- target_nonV[empty_nonV]
            I_prev[empty_nonV] <- 0
            R_prev[empty_nonV] <- 0
          }
        }
      }
      
      setup_annual_campaign(tt = tt)
    }
    
    if (scenario == "routine") {
      
      delay_source_week <- tt - delay_weeks
      
      if (delay_source_week >= 1) {
        
        effective_dose <- vacc_delayed[, delay_source_week]
        immunized <- VE_inf * effective_dose
        
        immunized <- pmin(immunized, S_prev)
        
        S_prev <- S_prev - immunized
        V_prev <- V_prev + immunized
        
        immunized_this_week[, tt] <- immunized
        age_protected_this_week[, tt] <- immunized
      }
    }
    
    allocate_weekly_doses(
      tt = tt,
      S_state = S_prev,
      I_state = I_prev,
      R_state = R_prev,
      V_state = V_prev
    )
    
    if (coverage_target_raw[tt] == 0 && scenario == "routine") {
      if (!is.na(current_campaign_year) &&
          ce_fit$year[tt] == current_campaign_year &&
          annual_target_pop > 0) {
        
        coverage_target_raw[tt] <-
          sum(annual_used_age[target_age_index], na.rm = TRUE) / annual_target_pop
        
        coverage_target_susceptible[tt] <-
          sum(annual_used_susceptible_age[target_age_index], na.rm = TRUE) / annual_target_pop
        
      } else {
        
        coverage_target_raw[tt] <- coverage_target_raw[tt - 1]
        coverage_target_susceptible[tt] <- coverage_target_susceptible[tt - 1]
      }
    }
    
    N_prev <- sum(S_prev) + sum(I_prev) + sum(R_prev) + sum(V_prev)
    
    import_count_week <- pars$import_count_year[year_idx] / 52
    
    local_infectious_frac[tt] <- sum(I_prev) / N_prev
    import_frac[tt] <- import_count_week / N_prev
    
    lambda[tt] <- pars$beta_t[tt] * (local_infectious_frac[tt] + import_frac[tt])
    
    p_inf <- 1 - exp(-lambda[tt])
    
    new_inf <- S_prev * p_inf
    recov <- I_prev * p_rec
    
    S_new <- S_prev - new_inf
    I_new <- I_prev + new_inf - recov
    R_new <- R_prev + recov
    V_new <- V_prev
    
    if (ageing_mode == "weekly") {
      S_new <- apply_ageing(S_new, ageing_prob_weekly)
      I_new <- apply_ageing(I_new, ageing_prob_weekly)
      R_new <- apply_ageing(R_new, ageing_prob_weekly)
      V_new <- apply_ageing(V_new, ageing_prob_weekly)
    }
    
    S[, tt] <- pmax0(S_new)
    I[, tt] <- pmax0(I_new)
    R[, tt] <- pmax0(R_new)
    V[, tt] <- pmax0(V_new)
    
    age_new_inf[, tt] <- new_inf
    age_recov[, tt] <- recov
    
    age_vaccinated_this_week[, tt] <- raw_allocation_age[, tt]
    
    total_population_modelled[tt] <- sum(S[, tt]) + sum(I[, tt]) + sum(R[, tt]) + sum(V[, tt])
    protected_frac[tt] <- sum(V[, tt]) / total_population_modelled[tt]
    effective_susceptible_frac[tt] <- sum(S[, tt]) / total_population_modelled[tt]
  }
  
  rownames(S) <- age_gr_levels
  rownames(I) <- age_gr_levels
  rownames(R) <- age_gr_levels
  rownames(V) <- age_gr_levels
  rownames(age_new_inf) <- age_gr_levels
  rownames(age_recov) <- age_gr_levels
  rownames(raw_allocation_age) <- age_gr_levels
  rownames(vacc_delayed) <- age_gr_levels
  rownames(wasted_dose) <- age_gr_levels
  rownames(immunized_this_week) <- age_gr_levels
  rownames(age_vaccinated_this_week) <- age_gr_levels
  rownames(age_protected_this_week) <- age_gr_levels
  
  campaign_summary <- dplyr::bind_rows(campaign_log)
  
  if (nrow(campaign_summary) > 0) {
    
    campaign_summary <- campaign_summary |>
      dplyr::rowwise() |>
      dplyr::mutate(
        total_used = sum(
          raw_allocation_age[
            target_age_index,
            ce_fit$year == calendar_year,
            drop = FALSE
          ],
          na.rm = TRUE
        ),
        total_used_susceptible = sum(
          vacc_delayed[
            target_age_index,
            ce_fit$year == calendar_year,
            drop = FALSE
          ],
          na.rm = TRUE
        ),
        total_immunized = sum(
          immunized_this_week[
            ,
            ce_fit$year == calendar_year,
            drop = FALSE
          ],
          na.rm = TRUE
        ),
        achieved_raw_coverage = total_used / target_population,
        achieved_susceptible_coverage = total_used_susceptible / target_population,
        first_allocation_week = {
          ww <- which(
            colSums(
              raw_allocation_age[
                target_age_index,
                ce_fit$year == calendar_year,
                drop = FALSE
              ],
              na.rm = TRUE
            ) > 0
          )
          if (length(ww) == 0) NA_integer_ else which(ce_fit$year == calendar_year)[ww[1]]
        },
        last_allocation_week = {
          ww <- which(
            colSums(
              raw_allocation_age[
                target_age_index,
                ce_fit$year == calendar_year,
                drop = FALSE
              ],
              na.rm = TRUE
            ) > 0
          )
          if (length(ww) == 0) NA_integer_ else which(ce_fit$year == calendar_year)[tail(ww, 1)]
        },
        real_vaccine_duration_weeks =
          ifelse(
            is.na(first_allocation_week) | is.na(last_allocation_week),
            0,
            last_allocation_week - first_allocation_week + 1
          )
      ) |>
      dplyr::ungroup()
  }
  
  weekly_df <- tibble::tibble(
    week = seq_len(T),
    week_start = ce_fit$week_start,
    calendar_year = ce_fit$year,
    scenario = scenario,
    total_infections = colSums(age_new_inf),
    total_S = colSums(S),
    total_I = colSums(I),
    total_R = colSums(R),
    total_V = colSums(V),
    total_population_modelled = total_population_modelled,
    effective_susceptible_frac = effective_susceptible_frac,
    protected_frac = protected_frac,
    local_infectious_frac = local_infectious_frac,
    import_frac = import_frac,
    lambda = lambda,
    vaccinated_this_week = colSums(raw_allocation_age),
    susceptible_vaccinated_this_week = colSums(vacc_delayed),
    protected_this_week = colSums(immunized_this_week),
    wasted_dose_this_week = colSums(wasted_dose),
    coverage_target_raw = coverage_target_raw,
    coverage_target_susceptible = coverage_target_susceptible
  )
  
  list(
    weekly_df = weekly_df,
    campaign_summary = campaign_summary,
    
    age_new_inf = age_new_inf,
    age_S = S,
    age_I = I,
    age_R = R,
    age_V = V,
    age_recov = age_recov,
    
    raw_allocation_age = raw_allocation_age,
    vacc_delayed = vacc_delayed,
    wasted_dose = wasted_dose,
    immunized_this_week = immunized_this_week,
    
    age_vaccinated_this_week = age_vaccinated_this_week,
    age_protected_this_week = age_protected_this_week,
    
    vacc_start_week = vacc_start_week,
    vacc_end_week = vacc_end_week,
    real_vaccine_duration_age = ifelse(
      is.na(vacc_start_week) | is.na(vacc_end_week),
      0,
      vacc_end_week - vacc_start_week + 1
    )
  )
}
simulate_postvacc_infections_age20 <- function(draws_df,
                                               ce_fit,
                                               years,
                                               prevacc_age20_test,
                                               N_mg,
                                               age_gr_levels,
                                               target_age_index,
                                               routine_schedule,
                                               ve_sus = 0,
                                               VE_inf = ve_sus,
                                               weekly_delivery_speed = 0.10,
                                               delay_weeks = 2L,
                                               gamma_fixed = NULL,
                                               anchor_V = FALSE,
                                               use_ageing = TRUE,
                                               ageing_mode = "annual") {
  
  if (is.null(prevacc_age20_test$draw_ids)) {
    stop("prevacc_age20_test$draw_ids is missing.")
  }
  
  ageing_mode <- match.arg(ageing_mode, c("annual", "weekly", "none"))
  
  if (!use_ageing) {
    ageing_mode <- "none"
  }
  
  draw_ids <- prevacc_age20_test$draw_ids
  n_draws <- length(draw_ids)
  
  T <- nrow(ce_fit)
  A <- length(age_gr_levels)
  
  age_pop_weekly <- make_age_pop_weekly(
    N_mg = N_mg,
    ce_fit = ce_fit,
    age_gr_levels = age_gr_levels
  )
  
  ageing_prob_weekly <- make_ageing_prob_weekly(
    N_mg = N_mg,
    age_gr_levels = age_gr_levels
  )
  
  ageing_prob_annual <- make_ageing_prob_annual(
    N_mg = N_mg,
    age_gr_levels = age_gr_levels
  )
  
  # Arrays: [age, week, draw]
  baseline_age_inf_array <- array(NA_real_, dim = c(A, T, n_draws))
  routine_age_inf_array  <- array(NA_real_, dim = c(A, T, n_draws))
  
  baseline_age_S_array <- array(NA_real_, dim = c(A, T, n_draws))
  routine_age_S_array  <- array(NA_real_, dim = c(A, T, n_draws))
  
  baseline_age_I_array <- array(NA_real_, dim = c(A, T, n_draws))
  routine_age_I_array  <- array(NA_real_, dim = c(A, T, n_draws))
  
  baseline_age_R_array <- array(NA_real_, dim = c(A, T, n_draws))
  routine_age_R_array  <- array(NA_real_, dim = c(A, T, n_draws))
  
  routine_age_V_array <- array(NA_real_, dim = c(A, T, n_draws))
  routine_age_protected_array <- array(NA_real_, dim = c(A, T, n_draws))
  
  # Weekly total matrices: [week, draw]
  baseline_weekly_inf_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  routine_weekly_inf_mat  <- matrix(NA_real_, nrow = T, ncol = n_draws)
  
  baseline_effective_S_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  routine_effective_S_mat  <- matrix(NA_real_, nrow = T, ncol = n_draws)
  
  baseline_local_I_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  routine_local_I_mat  <- matrix(NA_real_, nrow = T, ncol = n_draws)
  
  routine_protected_frac_mat <- matrix(NA_real_, nrow = T, ncol = n_draws)
  
  baseline_results_list <- vector("list", n_draws)
  routine_results_list  <- vector("list", n_draws)
  
  for (i in seq_len(n_draws)) {
    
    draw_id <- draw_ids[i]
    
    pars_draw <- make_pars_draw(
      draws_df = draws_df,
      draw_id = draw_id,
      gamma_fixed = gamma_fixed
    )
    
    baseline_i <- simulate_age_sir_20groups_vacc(
      pars = pars_draw,
      ce_fit = ce_fit,
      years = years,
      age_pop_weekly = age_pop_weekly,
      age_gr_levels = age_gr_levels,
      target_age_index = target_age_index,
      routine_schedule = routine_schedule,
      VE_inf = VE_inf,
      ve_sus = NULL,
      scenario = "baseline",
      weekly_delivery_speed = weekly_delivery_speed,
      delay_weeks = delay_weeks,
      ageing_prob_weekly = ageing_prob_weekly,
      ageing_prob_annual = ageing_prob_annual,
      use_ageing = use_ageing,
      ageing_mode = ageing_mode,
      anchor_V = anchor_V
    )
    
    routine_i <- simulate_age_sir_20groups_vacc(
      pars = pars_draw,
      ce_fit = ce_fit,
      years = years,
      age_pop_weekly = age_pop_weekly,
      age_gr_levels = age_gr_levels,
      target_age_index = target_age_index,
      routine_schedule = routine_schedule,
      VE_inf = VE_inf,
      ve_sus = NULL,
      scenario = "routine",
      weekly_delivery_speed = weekly_delivery_speed,
      delay_weeks = delay_weeks,
      ageing_prob_weekly = ageing_prob_weekly,
      ageing_prob_annual = ageing_prob_annual,
      use_ageing = use_ageing,
      ageing_mode = ageing_mode,
      anchor_V = anchor_V
    )
    
    baseline_results_list[[i]] <- baseline_i
    routine_results_list[[i]]  <- routine_i
    
    baseline_age_inf_array[, , i] <- baseline_i$age_new_inf
    routine_age_inf_array[, , i]  <- routine_i$age_new_inf
    
    baseline_age_S_array[, , i] <- baseline_i$age_S
    routine_age_S_array[, , i]  <- routine_i$age_S
    
    baseline_age_I_array[, , i] <- baseline_i$age_I
    routine_age_I_array[, , i]  <- routine_i$age_I
    
    baseline_age_R_array[, , i] <- baseline_i$age_R
    routine_age_R_array[, , i]  <- routine_i$age_R
    
    routine_age_V_array[, , i] <- routine_i$age_V
    routine_age_protected_array[, , i] <- routine_i$age_protected_this_week
    
    baseline_weekly_inf_mat[, i] <- baseline_i$weekly_df$total_infections
    routine_weekly_inf_mat[, i]  <- routine_i$weekly_df$total_infections
    
    baseline_effective_S_mat[, i] <- baseline_i$weekly_df$effective_susceptible_frac
    routine_effective_S_mat[, i]  <- routine_i$weekly_df$effective_susceptible_frac
    
    baseline_local_I_mat[, i] <- baseline_i$weekly_df$local_infectious_frac
    routine_local_I_mat[, i]  <- routine_i$weekly_df$local_infectious_frac
    
    routine_protected_frac_mat[, i] <- routine_i$weekly_df$protected_frac
  }
  
  dimnames(baseline_age_inf_array) <- list(
    age_group = age_gr_levels,
    week = seq_len(T),
    draw = seq_len(n_draws)
  )
  dimnames(routine_age_inf_array) <- dimnames(baseline_age_inf_array)
  
  dimnames(baseline_age_S_array) <- dimnames(baseline_age_inf_array)
  dimnames(routine_age_S_array)  <- dimnames(baseline_age_inf_array)
  
  dimnames(baseline_age_I_array) <- dimnames(baseline_age_inf_array)
  dimnames(routine_age_I_array)  <- dimnames(baseline_age_inf_array)
  
  dimnames(baseline_age_R_array) <- dimnames(baseline_age_inf_array)
  dimnames(routine_age_R_array)  <- dimnames(baseline_age_inf_array)
  
  dimnames(routine_age_V_array) <- dimnames(baseline_age_inf_array)
  dimnames(routine_age_protected_array) <- dimnames(baseline_age_inf_array)
  
  # Weekly total summary
  weekly_summary <- dplyr::bind_rows(
    tibble::tibble(
      scenario = "baseline",
      week = seq_len(T),
      week_start = ce_fit$week_start,
      calendar_year = ce_fit$year,
      infections_median = apply(baseline_weekly_inf_mat, 1, median, na.rm = TRUE),
      infections_low95 = apply(
        baseline_weekly_inf_mat,
        1,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      infections_hi95 = apply(
        baseline_weekly_inf_mat,
        1,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      ),
      effective_S_median = apply(baseline_effective_S_mat, 1, median, na.rm = TRUE),
      local_I_median = apply(baseline_local_I_mat, 1, median, na.rm = TRUE),
      protected_frac_median = 0
    ),
    tibble::tibble(
      scenario = "routine",
      week = seq_len(T),
      week_start = ce_fit$week_start,
      calendar_year = ce_fit$year,
      infections_median = apply(routine_weekly_inf_mat, 1, median, na.rm = TRUE),
      infections_low95 = apply(
        routine_weekly_inf_mat,
        1,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),
      infections_hi95 = apply(
        routine_weekly_inf_mat,
        1,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      ),
      effective_S_median = apply(routine_effective_S_mat, 1, median, na.rm = TRUE),
      local_I_median = apply(routine_local_I_mat, 1, median, na.rm = TRUE),
      protected_frac_median = apply(routine_protected_frac_mat, 1, median, na.rm = TRUE)
    )
  )
  
  # Age × week summary
  age_week_summary <- purrr::map_dfr(seq_len(A), function(a) {
    purrr::map_dfr(seq_len(T), function(tt) {
      
      xb <- baseline_age_inf_array[a, tt, ]
      xr <- routine_age_inf_array[a, tt, ]
      
      tibble::tibble(
        age_group = age_gr_levels[a],
        age_index = a,
        week = tt,
        week_start = ce_fit$week_start[tt],
        calendar_year = ce_fit$year[tt],
        
        baseline_median = median(xb, na.rm = TRUE),
        baseline_low95 = quantile(xb, 0.025, na.rm = TRUE),
        baseline_hi95 = quantile(xb, 0.975, na.rm = TRUE),
        
        routine_median = median(xr, na.rm = TRUE),
        routine_low95 = quantile(xr, 0.025, na.rm = TRUE),
        routine_hi95 = quantile(xr, 0.975, na.rm = TRUE),
        
        infections_averted_median = median(xb - xr, na.rm = TRUE),
        percent_reduction_median = median((xb - xr) / pmax(xb, 1e-9), na.rm = TRUE)
      )
    })
  })
  
  # Annual age-specific impact, draw-level
  annual_age_draws <- purrr::map_dfr(seq_len(n_draws), function(i) {
    purrr::map_dfr(seq_len(A), function(a) {
      tibble::tibble(
        draw_number = i,
        draw_id = draw_ids[i],
        age_group = age_gr_levels[a],
        age_index = a,
        week = seq_len(T),
        calendar_year = ce_fit$year,
        baseline_infections = baseline_age_inf_array[a, , i],
        routine_infections = routine_age_inf_array[a, , i]
      ) |>
        dplyr::group_by(draw_number, draw_id, age_group, age_index, calendar_year) |>
        dplyr::summarise(
          annual_baseline = sum(baseline_infections, na.rm = TRUE),
          annual_routine = sum(routine_infections, na.rm = TRUE),
          annual_averted = annual_baseline - annual_routine,
          annual_percent_reduction = annual_averted / pmax(annual_baseline, 1e-9),
          .groups = "drop"
        )
    })
  })
  
  annual_age_summary <- annual_age_draws |>
    dplyr::group_by(age_group, age_index, calendar_year) |>
    dplyr::summarise(
      baseline_median = median(annual_baseline, na.rm = TRUE),
      routine_median = median(annual_routine, na.rm = TRUE),
      averted_median = median(annual_averted, na.rm = TRUE),
      percent_reduction_median = median(annual_percent_reduction, na.rm = TRUE),
      
      baseline_low95 = quantile(annual_baseline, 0.025, na.rm = TRUE),
      baseline_hi95 = quantile(annual_baseline, 0.975, na.rm = TRUE),
      routine_low95 = quantile(annual_routine, 0.025, na.rm = TRUE),
      routine_hi95 = quantile(annual_routine, 0.975, na.rm = TRUE),
      averted_low95 = quantile(annual_averted, 0.025, na.rm = TRUE),
      averted_hi95 = quantile(annual_averted, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  total_annual_summary <- annual_age_draws |>
    dplyr::group_by(draw_number, draw_id, calendar_year) |>
    dplyr::summarise(
      annual_baseline = sum(annual_baseline, na.rm = TRUE),
      annual_routine = sum(annual_routine, na.rm = TRUE),
      annual_averted = sum(annual_averted, na.rm = TRUE),
      annual_percent_reduction = annual_averted / pmax(annual_baseline, 1e-9),
      .groups = "drop"
    ) |>
    dplyr::group_by(calendar_year) |>
    dplyr::summarise(
      baseline_median = median(annual_baseline, na.rm = TRUE),
      routine_median = median(annual_routine, na.rm = TRUE),
      averted_median = median(annual_averted, na.rm = TRUE),
      percent_reduction_median = median(annual_percent_reduction, na.rm = TRUE),
      percent_reduction_low95 = quantile(annual_percent_reduction, 0.025, na.rm = TRUE),
      percent_reduction_hi95 = quantile(annual_percent_reduction, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    draw_ids = draw_ids,
    n_draws = n_draws,
    age_pop_weekly = age_pop_weekly,
    ageing_prob_weekly = ageing_prob_weekly,
    ageing_prob_annual = ageing_prob_annual,
    ageing_mode = ageing_mode,
    use_ageing = use_ageing,
    anchor_V = anchor_V,
    
    baseline_results_list = baseline_results_list,
    routine_results_list = routine_results_list,
    
    baseline_age_inf_array = baseline_age_inf_array,
    routine_age_inf_array = routine_age_inf_array,
    
    baseline_age_S_array = baseline_age_S_array,
    routine_age_S_array = routine_age_S_array,
    baseline_age_I_array = baseline_age_I_array,
    routine_age_I_array = routine_age_I_array,
    baseline_age_R_array = baseline_age_R_array,
    routine_age_R_array = routine_age_R_array,
    
    routine_age_V_array = routine_age_V_array,
    routine_age_protected_array = routine_age_protected_array,
    
    baseline_weekly_inf_mat = baseline_weekly_inf_mat,
    routine_weekly_inf_mat = routine_weekly_inf_mat,
    
    weekly_summary = weekly_summary,
    age_week_summary = age_week_summary,
    annual_age_draws = annual_age_draws,
    annual_age_summary = annual_age_summary,
    total_annual_summary = total_annual_summary
  )
}

target_age_group <- "12"

target_age_index <- which(normalise_age(age_gr_levels) == target_age_group)

if (length(target_age_index) != 1) {
  stop("target_age_index must identify exactly one age group.")
}

routine_schedule <- tibble::tibble(
  calendar_year = 2015:2019,
  coverage = rep(0.50, 5)
)



prevacc_age20_test <- simulate_prevacc_infections_age20(
  draws_df = draws_df,
  ce_fit = ce_fit,
  years = years,
  N_mg = N_mg,
  age_gr_levels = age_gr_levels,
  n_draws = 50,
  gamma_fixed = if (exists("GAMMA_WEEK")) GAMMA_WEEK else NULL,
  seed = 123
)


draw_id_test <- prevacc_age20_test$draw_ids[1]

pars_test <- make_pars_draw(
  draws_df = draws_df,
  draw_id = draw_id_test,
  gamma_fixed = if (exists("GAMMA_WEEK")) GAMMA_WEEK else NULL
)

age20_baseline_test <- simulate_age_sir_20groups_vacc(
  pars = pars_test,
  ce_fit = ce_fit,
  years = years,
  age_pop_weekly = age_pop_weekly,
  age_gr_levels = age_gr_levels,
  target_age_index = target_age_index,
  routine_schedule = routine_schedule,
  VE_inf = ve_sus,
  scenario = "baseline",
  ageing_prob_weekly = make_ageing_prob_weekly(N_mg, age_gr_levels),
  ageing_prob_annual = make_ageing_prob_annual(N_mg, age_gr_levels),
  use_ageing = TRUE,
  ageing_mode = "annual",
  anchor_V = FALSE
)

age20_routine_test <- simulate_age_sir_20groups_vacc(
  pars = pars_test,
  ce_fit = ce_fit,
  years = years,
  age_pop_weekly = age_pop_weekly,
  age_gr_levels = age_gr_levels,
  target_age_index = target_age_index,
  routine_schedule = routine_schedule,
  VE_inf = ve_sus,
  scenario = "routine",
  weekly_delivery_speed = 0.10,
  delay_weeks = 2L,
  ageing_prob_weekly = make_ageing_prob_weekly(N_mg, age_gr_levels),
  ageing_prob_annual = make_ageing_prob_annual(N_mg, age_gr_levels),
  use_ageing = TRUE,
  ageing_mode = "annual",
  anchor_V = FALSE
)



postvacc_age20_test <- simulate_postvacc_infections_age20(
  draws_df = draws_df,
  ce_fit = ce_fit,
  years = years,
  prevacc_age20_test = prevacc_age20_test,
  N_mg = N_mg,
  age_gr_levels = age_gr_levels,
  target_age_index = target_age_index,
  routine_schedule = routine_schedule,
  ve_sus = ve_sus,
  VE_inf = ve_sus,
  weekly_delivery_speed = 0.10,
  delay_weeks = 2L,
  gamma_fixed = if (exists("GAMMA_WEEK")) GAMMA_WEEK else NULL,
  anchor_V = FALSE,
  use_ageing = TRUE,
  ageing_mode = "annual"
)


#### plot
age_annual_curve_df <- postvacc_age20_test$annual_age_draws |>
  dplyr::select(
    draw_number,
    draw_id,
    age_group,
    age_index,
    calendar_year,
    annual_baseline,
    annual_routine
  ) |>
  tidyr::pivot_longer(
    cols = c(annual_baseline, annual_routine),
    names_to = "scenario",
    values_to = "annual_infections"
  ) |>
  dplyr::mutate(
    scenario = dplyr::recode(
      scenario,
      annual_baseline = "baseline",
      annual_routine = "routine"
    ),
    scenario = factor(scenario, levels = c("baseline", "routine")),
    age_group = factor(age_group, levels = age_gr_levels)
  ) |>
  dplyr::group_by(age_group, age_index, calendar_year, scenario) |>
  dplyr::summarise(
    median = median(annual_infections, na.rm = TRUE),
    low95 = quantile(annual_infections, 0.025, na.rm = TRUE),
    hi95 = quantile(annual_infections, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(
  age_annual_curve_df,
  aes(x = calendar_year, y = median, colour = scenario, fill = scenario)
) +
  geom_ribbon(
    aes(ymin = low95, ymax = hi95),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.2) +
  facet_wrap(~ age_group, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  scale_x_continuous(breaks = sort(unique(age_annual_curve_df$calendar_year))) +
  labs(
    x = "Calendar year",
    y = "Annual infections",
    colour = NULL,
    fill = NULL,
    title = "Age-specific annual infections",
    subtitle = "Baseline versus routine vaccination"
  ) +
  theme_lancet_like() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) + scale_y_log10()


total_S_weekly_df <- dplyr::bind_rows(
  purrr::imap_dfr(
    postvacc_age20_test$baseline_results_list,
    function(sim, i) {
      sim$weekly_df |>
        dplyr::transmute(
          draw_number = i,
          draw_id = postvacc_age20_test$draw_ids[i],
          scenario = "baseline",
          week,
          week_start,
          calendar_year,
          total_S,
          total_V,
          total_population_modelled,
          effective_susceptible_frac,
          protected_frac,
          coverage_target_raw = 0,
          coverage_target_susceptible = 0
        )
    }
  ),
  purrr::imap_dfr(
    postvacc_age20_test$routine_results_list,
    function(sim, i) {
      sim$weekly_df |>
        dplyr::transmute(
          draw_number = i,
          draw_id = postvacc_age20_test$draw_ids[i],
          scenario = "routine",
          week,
          week_start,
          calendar_year,
          total_S,
          total_V,
          total_population_modelled,
          effective_susceptible_frac,
          protected_frac,
          coverage_target_raw,
          coverage_target_susceptible
        )
    }
  )
)


total_S_summary <- total_S_weekly_df |>
  dplyr::group_by(scenario, week, week_start, calendar_year) |>
  dplyr::summarise(
    total_S_median = median(total_S, na.rm = TRUE),
    total_S_low95 = quantile(total_S, 0.025, na.rm = TRUE),
    total_S_hi95 = quantile(total_S, 0.975, na.rm = TRUE),
    total_V_median = median(total_V, na.rm = TRUE),
    effective_S_frac_median = median(effective_susceptible_frac, na.rm = TRUE),
    protected_frac_median = median(protected_frac, na.rm = TRUE),
    coverage_raw_median = median(coverage_target_raw, na.rm = TRUE),
    coverage_susceptible_median = median(coverage_target_susceptible, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    scenario = factor(scenario, levels = c("baseline", "routine"))
  )


ggplot(
  total_S_summary,
  aes(x = week_start, y = total_S_median, colour = scenario, fill = scenario)
) +
  geom_ribbon(
    aes(ymin = total_S_low95, ymax = total_S_hi95),
    alpha = 0.12,
    colour = NA
  ) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Week",
    y = "Total susceptible population",
    colour = NULL,
    fill = NULL,
    title = "Total susceptible trajectory",
    subtitle = "Baseline versus routine vaccination"
  ) +
  theme_lancet_like() +
  theme(legend.position = "top")

ggplot(
  total_S_summary,
  aes(x = week_start, y = effective_S_frac_median, colour = scenario)
) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    x = "Week",
    y = "Effective susceptible fraction",
    colour = NULL,
    title = "Effective susceptible fraction",
    subtitle = "S / total population, excluding vaccine-protected compartment"
  ) +
  theme_lancet_like() +
  theme(legend.position = "top")

routine_total_delivery_df <- total_S_summary |>
  dplyr::filter(scenario == "routine") |>
  dplyr::select(
    week_start,
    calendar_year,
    effective_S_frac_median,
    protected_frac_median,
    coverage_raw_median,
    coverage_susceptible_median
  ) |>
  tidyr::pivot_longer(
    cols = c(
      effective_S_frac_median,
      protected_frac_median,
      coverage_raw_median,
      coverage_susceptible_median
    ),
    names_to = "metric",
    values_to = "value"
  ) |>
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      effective_S_frac_median = "Effective susceptible fraction",
      protected_frac_median = "Protected fraction",
      coverage_raw_median = "Raw target coverage",
      coverage_susceptible_median = "Susceptible-dose coverage"
    )
  )

ggplot(
  routine_total_delivery_df,
  aes(x = week_start, y = value, colour = metric)
) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    x = "Week",
    y = "Fraction",
    colour = NULL,
    title = "Routine scenario: susceptible, coverage, and protection trajectories",
    subtitle = "Checking whether delivery variables are reflected in state trajectories"
  ) +
  theme_lancet_like() +
  theme(legend.position = "top")

age_S_weekly_df <- dplyr::bind_rows(
  purrr::imap_dfr(
    postvacc_age20_test$baseline_results_list,
    function(sim, i) {
      purrr::map_dfr(seq_along(age_gr_levels), function(a) {
        tibble::tibble(
          draw_number = i,
          draw_id = postvacc_age20_test$draw_ids[i],
          scenario = "baseline",
          age_group = age_gr_levels[a],
          age_index = a,
          week = seq_len(nrow(ce_fit)),
          week_start = ce_fit$week_start,
          calendar_year = ce_fit$year,
          S = sim$age_S[a, ],
          V = sim$age_V[a, ],
          N_age = postvacc_age20_test$age_pop_weekly[a, ]
        )
      })
    }
  ),
  purrr::imap_dfr(
    postvacc_age20_test$routine_results_list,
    function(sim, i) {
      purrr::map_dfr(seq_along(age_gr_levels), function(a) {
        tibble::tibble(
          draw_number = i,
          draw_id = postvacc_age20_test$draw_ids[i],
          scenario = "routine",
          age_group = age_gr_levels[a],
          age_index = a,
          week = seq_len(nrow(ce_fit)),
          week_start = ce_fit$week_start,
          calendar_year = ce_fit$year,
          S = sim$age_S[a, ],
          V = sim$age_V[a, ],
          N_age = postvacc_age20_test$age_pop_weekly[a, ]
        )
      })
    }
  )
) |>
  dplyr::mutate(
    age_group = factor(age_group, levels = age_gr_levels),
    scenario = factor(scenario, levels = c("baseline", "routine")),
    S_frac = S / N_age,
    V_frac = V / N_age
  )

age_S_summary <- age_S_weekly_df |>
  dplyr::group_by(scenario, age_group, age_index, week, week_start, calendar_year) |>
  dplyr::summarise(
    S_median = median(S, na.rm = TRUE),
    S_low95 = quantile(S, 0.025, na.rm = TRUE),
    S_hi95 = quantile(S, 0.975, na.rm = TRUE),
    S_frac_median = median(S_frac, na.rm = TRUE),
    V_median = median(V, na.rm = TRUE),
    V_frac_median = median(V_frac, na.rm = TRUE),
    .groups = "drop"
  )


ggplot(
  age_S_summary,
  aes(x = week_start, y = S_median, colour = scenario, fill = scenario)
) +
  geom_ribbon(
    aes(ymin = S_low95, ymax = S_hi95),
    alpha = 0.10,
    colour = NA
  ) +
  geom_line(linewidth = 0.65) +
  facet_wrap(~ age_group, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Week",
    y = "Susceptible population",
    colour = NULL,
    fill = NULL,
    title = "Age-specific susceptible trajectories",
    subtitle = "Baseline versus routine vaccination"
  ) +
  theme_lancet_like() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(
  age_S_summary,
  aes(x = week_start, y = S_frac_median, colour = scenario)
) +
  geom_line(linewidth = 0.65) +
  facet_wrap(~ age_group, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(
    x = "Week",
    y = "Susceptible fraction within age group",
    colour = NULL,
    title = "Age-specific susceptible fraction trajectories",
    subtitle = "S divided by age-specific population"
  ) +
  theme_lancet_like() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(
  age_S_summary |>
    dplyr::filter(scenario == "routine"),
  aes(x = week_start, y = V_frac_median)
) +
  geom_line(linewidth = 0.65, colour = "#0B2C4D") +
  facet_wrap(~ age_group, scales = "free_y", ncol = 4) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  labs(
    x = "Week",
    y = "V fraction within age group",
    title = "Age-specific vaccine-protected fraction trajectories",
    subtitle = "Routine scenario only; checks ageing and delayed protection"
  ) +
  theme_lancet_like() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
#########################################
### old 
weekly_plot_df <- postvacc_age20_test$weekly_summary |>
  dplyr::mutate(
    scenario = factor(scenario, levels = c("baseline", "routine"))
  )

ggplot(
  weekly_plot_df,
  aes(x = week_start, y = infections_median, colour = scenario, fill = scenario)
) +
  geom_ribbon(
    aes(ymin = infections_low95, ymax = infections_hi95),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Week",
    y = "Weekly infections",
    title = "Counterfactual weekly infections",
    subtitle = "Baseline versus routine vaccination of 12-year-old cohort"
  ) +
  theme_lancet_like()

### total aggregated
annual_total_long <- postvacc_age20_test$annual_age_draws |>
  dplyr::group_by(draw_number, draw_id, calendar_year) |>
  dplyr::summarise(
    baseline = sum(annual_baseline, na.rm = TRUE),
    routine = sum(annual_routine, na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = c(baseline, routine),
    names_to = "scenario",
    values_to = "annual_infections"
  ) |>
  dplyr::group_by(calendar_year, scenario) |>
  dplyr::summarise(
    median = median(annual_infections, na.rm = TRUE),
    low95 = quantile(annual_infections, 0.025, na.rm = TRUE),
    hi95 = quantile(annual_infections, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    scenario = factor(scenario, levels = c("baseline", "routine"))
  )

ggplot(
  annual_total_long,
  aes(x = factor(calendar_year), y = median, fill = scenario)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(
    aes(ymin = low95, ymax = hi95),
    position = position_dodge(width = 0.75),
    width = 0.2,
    linewidth = 0.35
  ) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Calendar year",
    y = "Annual infections",
    fill = NULL,
    title = "Annual infections under baseline and routine vaccination",
    subtitle = "Posterior median and 95% interval"
  ) +
  theme_lancet_like()

### annual % reduction
annual_reduction_df <- postvacc_age20_test$total_annual_summary |>
  dplyr::mutate(
    percent_reduction_median = 100 * percent_reduction_median,
    percent_reduction_low95 = 100 * percent_reduction_low95,
    percent_reduction_hi95 = 100 * percent_reduction_hi95
  )

ggplot(
  annual_reduction_df,
  aes(x = factor(calendar_year), y = percent_reduction_median)
) +
  geom_col(width = 0.65, fill = "#0B2C4D") +
  geom_errorbar(
    aes(ymin = percent_reduction_low95, ymax = percent_reduction_hi95),
    width = 0.2,
    linewidth = 0.35,
    colour = "grey25"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100)
  ) +
  labs(
    x = "Calendar year",
    y = "Reduction in annual infections",
    title = "Annual percentage reduction under routine vaccination",
    subtitle = "Routine vaccination of 12-year-old cohort"
  ) +
  theme_lancet_like()

#### 
age_cumulative_long <- postvacc_age20_test$annual_age_draws |>
  dplyr::group_by(draw_number, draw_id, age_group, age_index) |>
  dplyr::summarise(
    baseline = sum(annual_baseline, na.rm = TRUE),
    routine = sum(annual_routine, na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = c(baseline, routine),
    names_to = "scenario",
    values_to = "cumulative_infections"
  ) |>
  dplyr::group_by(age_group, age_index, scenario) |>
  dplyr::summarise(
    median = median(cumulative_infections, na.rm = TRUE),
    low95 = quantile(cumulative_infections, 0.025, na.rm = TRUE),
    hi95 = quantile(cumulative_infections, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    age_group = factor(age_group, levels = age_gr_levels),
    scenario = factor(scenario, levels = c("baseline", "routine"))
  )

ggplot(
  age_cumulative_long,
  aes(x = age_group, y = median, fill = scenario)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(
    aes(ymin = low95, ymax = hi95),
    position = position_dodge(width = 0.75),
    width = 0.2,
    linewidth = 0.3
  ) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Age group",
    y = "Cumulative infections",
    fill = NULL,
    title = "Cumulative infections by age group",
    subtitle = "Baseline versus routine vaccination"
  ) +
  theme_lancet_like() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


##
age_averted_df <- postvacc_age20_test$annual_age_draws |>
  dplyr::group_by(draw_number, draw_id, age_group, age_index) |>
  dplyr::summarise(
    cumulative_baseline = sum(annual_baseline, na.rm = TRUE),
    cumulative_routine = sum(annual_routine, na.rm = TRUE),
    cumulative_averted = cumulative_baseline - cumulative_routine,
    cumulative_percent_reduction =
      cumulative_averted / pmax(cumulative_baseline, 1e-9),
    .groups = "drop"
  ) |>
  dplyr::group_by(age_group, age_index) |>
  dplyr::summarise(
    averted_median = median(cumulative_averted, na.rm = TRUE),
    averted_low95 = quantile(cumulative_averted, 0.025, na.rm = TRUE),
    averted_hi95 = quantile(cumulative_averted, 0.975, na.rm = TRUE),
    percent_reduction_median = median(cumulative_percent_reduction, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    age_group = factor(age_group, levels = age_gr_levels),
    percent_reduction_median = 100 * percent_reduction_median
  )

ggplot(
  age_averted_df,
  aes(x = age_group, y = averted_median)
) +
  geom_col(fill = "#0B2C4D", width = 0.65) +
  geom_errorbar(
    aes(ymin = averted_low95, ymax = averted_hi95),
    width = 0.2,
    linewidth = 0.3,
    colour = "grey25"
  ) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(
    x = "Age group",
    y = "Cumulative infections averted",
    title = "Cumulative infections averted by age group",
    subtitle = "Routine vaccination of 12-year-old cohort"
  ) +
  theme_lancet_like() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(
  age_averted_df,
  aes(x = age_group, y = percent_reduction_median)
) +
  geom_col(fill = "#0B2C4D", width = 0.65) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100)
  ) +
  labs(
    x = "Age group",
    y = "Cumulative percent reduction",
    title = "Percent reduction in cumulative infections by age group",
    subtitle = "Routine vaccination of 12-year-old cohort"
  ) +
  theme_lancet_like() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

