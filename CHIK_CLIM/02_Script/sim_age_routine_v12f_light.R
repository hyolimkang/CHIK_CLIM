# ---------------------------------------------------------------------------
# sim_age_routine_v12f_light.R
#
# LIGHT version of sim_age_routine_v12f.R.
#
# Runs the age-structured v12f counterfactual with ONLY the posterior MEDIAN
# parameters (a single baseline run + a single routine run). Use this to quickly
# check that the simulation runs and that baseline reproduces the Stan median
# trajectory, before launching the full draw-level run in sim_age_routine_v12f.R.
#
# Run after fit_chik_ceara_stan_weekly.R (USE_DEMOGRAPHIC_FLOW <- TRUE).
#
# Objects created:
#   - baseline_sim_light, routine_sim_light, cf_weekly_light, cf_annual_light
#   - p_light_weekly, p_light_annual_averted, p_light_susceptible
# ---------------------------------------------------------------------------

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "posterior", "scales", "tibble")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
  library(scales)
  library(tibble)
})

required_objects <- c("fit_v12e", "ce_fit", "years", "UF_CODE", "UF_NAME")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    "\nRun source('02_Script/fit_chik_ceara_stan_weekly.R') first."
  )
}

# ---- Scenario settings -----------------------------------------------------
target_age <- 12L

routine_schedule <- tibble::tibble(
  calendar_year = c(2015L, 2016L, 2017, 2018, 2019),
  coverage = c(0.50, 0.50, 0.50, 0.5, 0.5)
)

ve_sus <- 0.989
max_age <- 100L

if (exists("stan_data_v12e") && !is.null(stan_data_v12e$mu_death_annual)) {
  mu_death_annual <- stan_data_v12e$mu_death_annual
} else if (exists("BRAZIL_ANNUAL_DEATH_RATE")) {
  mu_death_annual <- BRAZIL_ANNUAL_DEATH_RATE
} else {
  mu_death_annual <- 0.0066
}

utils::globalVariables(c(
  "calendar_year", "effective_susceptible_frac", "infections_averted",
  "scenario", "t", "total_infections", "coverage",
  "vaccinated_this_week", "vaccination_coverage", "week_start"
))

# ---- Helpers ---------------------------------------------------------------
extract_fit_draws <- function(fit) {
  if ("CmdStanMCMC" %in% class(fit)) {
    as.data.frame(fit$draws(format = "draws_df"))
  } else {
    as.data.frame(fit)
  }
}

extract_index <- function(x) as.integer(sub(".*\\[([0-9]+)\\].*", "\\1", x))

get_med_vector <- function(draws_df, var) {
  cols <- grep(paste0("^", var, "\\["), names(draws_df), value = TRUE)
  if (length(cols) == 0) stop("No posterior columns found for variable: ", var)
  cols <- cols[order(extract_index(cols))]
  apply(draws_df[, cols, drop = FALSE], 2, stats::median, na.rm = TRUE)
}

get_med_scalar <- function(draws_df, var) {
  if (!var %in% names(draws_df)) stop("No posterior column found for variable: ", var)
  stats::median(draws_df[[var]], na.rm = TRUE)
}

make_initial_age_population <- function(total_pop, max_age = 100L) {
  ages <- 0:max_age
  if (exists("age_pop_initial", inherits = TRUE)) {
    age_pop <- get("age_pop_initial", inherits = TRUE)
    if (!all(c("age", "population") %in% names(age_pop))) {
      stop("age_pop_initial must have columns: age, population")
    }
    out <- rep(0, length(ages))
    out[match(age_pop$age, ages, nomatch = 0)] <- age_pop$population[age_pop$age %in% ages]
    return(out / sum(out) * total_pop)
  }
  weights <- exp(-0.018 * ages)
  weights <- weights / sum(weights)
  weights * total_pop
}

age_one_week <- function(x) {
  out <- x
  move <- x / 52
  out <- out - move
  out[-1] <- out[-1] + move[-length(move)]
  out[length(out)] <- out[length(out)] + move[length(move)]
  out
}

validate_routine_schedule <- function(schedule) {
  if (is.null(schedule)) return(NULL)
  required_cols <- c("calendar_year", "coverage")
  if (length(setdiff(required_cols, names(schedule))) > 0) {
    stop("routine_schedule must contain columns: calendar_year, coverage", call. = FALSE)
  }
  schedule <- schedule |>
    dplyr::mutate(
      calendar_year = as.integer(.data$calendar_year),
      coverage = as.numeric(.data$coverage)
    ) |>
    dplyr::arrange(.data$calendar_year)
  if (any(is.na(schedule$calendar_year)) || any(is.na(schedule$coverage))) {
    stop("routine_schedule contains missing values.", call. = FALSE)
  }
  if (any(schedule$coverage < 0 | schedule$coverage > 1)) {
    stop("routine_schedule$coverage must be between 0 and 1.", call. = FALSE)
  }
  schedule
}

coverage_for_year <- function(schedule, calendar_year) {
  if (is.null(schedule)) return(0)
  year_match <- schedule$calendar_year == calendar_year
  if (!any(year_match)) return(0)
  schedule$coverage[which(year_match)[1]]
}

simulate_age_sir <- function(pars, ce_fit, years, scenario = "baseline",
                             target_age = 12L,
                             routine_schedule = NULL,
                             ve_sus = 0.80,
                             mu_death_annual = 0.0066,
                             max_age = 100L,
                             add_births = TRUE) {
  
  ages <- 0:max_age
  A <- length(ages)
  N_time <- nrow(ce_fit)
  
  pmax0 <- function(x) pmax(x, 0)
  routine_schedule <- validate_routine_schedule(routine_schedule)
  
  p_death <- 1 - exp(-mu_death_annual / 52)
  p_rec <- 1 - exp(-pars$gamma)
  
  # Initial age distribution
  pop_age <- make_initial_age_population(ce_fit$population[1], max_age)
  
  # All-or-nothing vaccine structure:
  # S = susceptible
  # I = infectious
  # R = recovered/immune after infection
  # V = vaccine-protected, not susceptible and not infectious
  S <- pop_age * pars$S0_frac
  I <- pop_age / sum(pop_age) * pars$I0
  R <- pmax0(pop_age - S - I)
  V <- rep(0, A)
  
  out <- vector("list", N_time)
  
  for (tt in seq_len(N_time)) {
    
    year_idx <- match(ce_fit$year[tt], years)
    if (is.na(year_idx)) {
      stop("Year not found in years vector: ", ce_fit$year[tt])
    }
    
    beta_t <- pars$beta_t[tt]
    
    vaccinated_this_week <- 0
    protected_this_week <- 0
    vaccine_failure_this_week <- 0
    
    # -----------------------------------------------------------------------
    # 1. Routine vaccination at the first modelled week of each calendar year
    # -----------------------------------------------------------------------
    first_week_start_this_year <- min(
      ce_fit$week_start[ce_fit$year == ce_fit$year[tt]]
    )
    
    coverage_this_year <- if (scenario == "routine") {
      coverage_for_year(routine_schedule, ce_fit$year[tt])
    } else {
      0
    }
    
    vaccination_due <- ce_fit$week_start[tt] == first_week_start_this_year &&
      coverage_this_year > 0
    
    if (vaccination_due) {
      
      age_idx <- which(ages %in% target_age)
      if (length(age_idx) == 0) {
        stop("target_age is outside the modelled age range.")
      }
      
      # Vaccination is applied to the susceptible part of the target-age cohort.
      # All-or-nothing interpretation:
      # - coverage_this_year = fraction offered/receiving vaccine among susceptible age-12
      # - ve_sus = probability of successful protection among vaccinated susceptibles
      # - successful vaccinees move to V
      # - vaccine failures remain in S
      vacc_s <- coverage_this_year * S[age_idx]
      protected_s <- ve_sus * vacc_s
      failed_s <- (1 - ve_sus) * vacc_s
      
      S[age_idx] <- S[age_idx] - protected_s
      V[age_idx] <- V[age_idx] + protected_s
      
      vaccinated_this_week <- sum(vacc_s)
      protected_this_week <- sum(protected_s)
      vaccine_failure_this_week <- sum(failed_s)
    }
    
    # -----------------------------------------------------------------------
    # 2. Infection and recovery
    # -----------------------------------------------------------------------
    N_prev <- sum(S) + sum(I) + sum(R) + sum(V)
    
    if (tt == 1L) {
      
      import_count_week <- 0
      local_infectious_frac <- 0
      import_frac <- 0
      infectious_pressure <- 0
      lambda_t <- 0
      
      new_inf <- rep(0, A)
      rec <- rep(0, A)
      
    } else {
      
      import_count_week <- pars$import_count_year[year_idx] / 52
      
      local_infectious_frac <- sum(I) / N_prev
      import_frac <- import_count_week / N_prev
      infectious_pressure <- local_infectious_frac + import_frac
      
      lambda_t <- beta_t * infectious_pressure
      
      p_inf <- 1 - exp(-lambda_t)
      
      # Only truly susceptible people can be infected.
      # V is fully protected under all-or-nothing successful vaccination.
      new_inf <- S * p_inf
      rec <- I * p_rec
      
      S <- S - new_inf
      I <- I + new_inf - rec
      R <- R + rec
      
      if (add_births) {
        S[1] <- S[1] + ce_fit$births_weekly[tt]
      }
    }
    
    # -----------------------------------------------------------------------
    # 3. Natural mortality
    # -----------------------------------------------------------------------
    deaths_this_week <- 0
    
    if (tt > 1L) {
      surv <- 1 - p_death
      
      pre_death_total <- sum(S) + sum(I) + sum(R) + sum(V)
      
      S <- S * surv
      I <- I * surv
      R <- R * surv
      V <- V * surv
      
      deaths_this_week <- pre_death_total * p_death
    }
    
    S <- pmax0(S)
    I <- pmax0(I)
    R <- pmax0(R)
    V <- pmax0(V)
    
    total_model_pop <- sum(S) + sum(I) + sum(R) + sum(V)
    
    # In all-or-nothing structure, effective susceptible pool is simply S.
    effective_S <- sum(S)
    
    out[[tt]] <- tibble::tibble(
      t = tt,
      week_start = ce_fit$week_start[tt],
      calendar_year = ce_fit$year[tt],
      scenario = scenario,
      
      total_infections = sum(new_inf),
      new_inf_unvaccinated = sum(new_inf),
      new_inf_vaccinated = 0,
      
      total_population_modelled = total_model_pop,
      pop_balance_resid = total_model_pop - ce_fit$population[tt],
      deaths_this_week = deaths_this_week,
      
      local_infectious_frac = local_infectious_frac,
      import_frac = import_frac,
      infectious_pressure = infectious_pressure,
      lambda = lambda_t,
      
      susceptible_count = sum(S),
      infectious_count = sum(I),
      recovered_count = sum(R),
      protected_count = sum(V),
      
      susceptible_frac = sum(S) / total_model_pop,
      infectious_frac = sum(I) / total_model_pop,
      recovered_frac = sum(R) / total_model_pop,
      protected_frac = sum(V) / total_model_pop,
      
      effective_susceptible_frac = effective_S / total_model_pop,
      
      vaccinated_this_week = vaccinated_this_week,
      protected_this_week = protected_this_week,
      vaccine_failure_this_week = vaccine_failure_this_week,
      vaccination_coverage = coverage_this_year
    )
    
    # -----------------------------------------------------------------------
    # 4. Ageing
    # -----------------------------------------------------------------------
    if (tt < N_time) {
      S <- age_one_week(S)
      I <- age_one_week(I)
      R <- age_one_week(R)
      V <- age_one_week(V)
    }
  }
  
  dplyr::bind_rows(out)
}
# ---- Median run ------------------------------------------------------------
draws_df <- extract_fit_draws(fit_v12e)

pars_med <- list(
  beta_t = get_med_vector(draws_df, "beta_t"),
  import_count_year = get_med_vector(draws_df, "import_count_year"),
  rho = get_med_scalar(draws_df, "rho"),
  S0_frac = get_med_scalar(draws_df, "S0_frac"),
  I0 = get_med_scalar(draws_df, "I0"),
  gamma = if (exists("GAMMA_WEEK")) GAMMA_WEEK else 1.0
)

routine_schedule <- validate_routine_schedule(routine_schedule)

baseline_sim_light <- simulate_age_sir(
  pars = pars_med, ce_fit = ce_fit, years = years, scenario = "baseline",
  target_age = target_age, ve_sus = ve_sus,
  mu_death_annual = mu_death_annual, max_age = max_age
)

routine_sim_light <- simulate_age_sir(
  pars = pars_med, ce_fit = ce_fit, years = years, scenario = "routine",
  target_age = target_age, routine_schedule = routine_schedule, ve_sus = ve_sus,
  mu_death_annual = mu_death_annual, max_age = max_age
)

# Sanity: baseline (median params) vs Stan median new_inf.
# NOTE: because the dynamics are nonlinear, median(params) propagated forward is
# not identical to the median of Stan trajectories, so a small gap is expected.
stan_new_inf_med <- get_med_vector(draws_df, "new_inf")
baseline_check_light <- tibble::tibble(
  max_abs_diff = max(abs(baseline_sim_light$total_infections - stan_new_inf_med), na.rm = TRUE),
  cor_baseline_stan = stats::cor(baseline_sim_light$total_infections, stan_new_inf_med),
  baseline_peak_week = ce_fit$week_start[which.max(baseline_sim_light$total_infections)],
  stan_peak_week = ce_fit$week_start[which.max(stan_new_inf_med)]
)

message("\n[LIGHT] Baseline (median params) vs Stan median new_inf:")
print(baseline_check_light)

message("\n[LIGHT] Demographic drift (baseline modelled vs observed population):")
print(
  baseline_sim_light |>
    dplyr::summarise(
      max_abs_resid = max(abs(.data$pop_balance_resid), na.rm = TRUE),
      final_resid = dplyr::last(.data$pop_balance_resid),
      total_deaths = sum(.data$deaths_this_week, na.rm = TRUE)
    )
)

# ---- Aggregate + compare ---------------------------------------------------
cf_weekly_light <- dplyr::bind_rows(baseline_sim_light, routine_sim_light)

cf_annual_light <- cf_weekly_light |>
  dplyr::group_by(scenario, calendar_year) |>
  dplyr::summarise(total_infections = sum(total_infections, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = scenario, values_from = total_infections) |>
  dplyr::mutate(
    infections_averted = baseline - routine,
    percent_reduction = infections_averted / pmax(baseline, 1e-9)
  )

message("\n[LIGHT] Annual infections averted:")
print(as.data.frame(cf_annual_light), digits = 6)

# ---- Diagnostics: is the averted effect plausible? -------------------------
# These quantify WHY routine can wipe out outbreaks from a modest vaccinated
# share: a near-threshold R0 (just above 1) makes the epidemic knife-edge
# sensitive, so a small susceptibility reduction crosses R_eff = 1.
p_rec_med <- 1 - exp(-pars_med$gamma)
R0_t <- pars_med$beta_t / p_rec_med   # basic reproduction number per week

n_last <- nrow(ce_fit)
base_S_last <- baseline_sim_light$effective_susceptible_frac[n_last] *
  baseline_sim_light$total_population_modelled[n_last]
rout_S_last <- routine_sim_light$effective_susceptible_frac[n_last] *
  routine_sim_light$total_population_modelled[n_last]

vacc_total <- sum(routine_sim_light$vaccinated_this_week, na.rm = TRUE)

diagnostics_light <- tibble::tibble(
  R0_median = stats::median(R0_t),
  R0_min = min(R0_t),
  R0_max = max(R0_t),
  weeks_R0_above_1 = sum(R0_t > 1),
  vaccinated_total = vacc_total,
  vaccinated_frac_of_mean_pop = vacc_total / mean(ce_fit$population),
  cumulative_infections_averted = sum(baseline_sim_light$total_infections) -
    sum(routine_sim_light$total_infections),
  final_effective_susceptible_gap = rout_S_last - base_S_last
)

message("\n[LIGHT] Plausibility diagnostics:")
message("  - If R0_median is only slightly > 1, the model is knife-edge: small")
message("    vaccination tips R_eff below 1 and removes outbreaks entirely.")
message("  - cumulative_infections_averted should be of the same order as")
message("    final_effective_susceptible_gap (averted = susceptibles not infected).")
print(as.data.frame(diagnostics_light), digits = 6)

# ---- Plots -----------------------------------------------------------------
routine_schedule_label <- paste(
  sprintf("%d: %.0f%%", routine_schedule$calendar_year, 100 * routine_schedule$coverage),
  collapse = ", "
)

p_light_weekly <- ggplot(
  cf_weekly_light,
  aes(x = week_start, y = total_infections, colour = scenario)
) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = scales::label_number()) +
  labs(
    x = NULL, y = "Weekly infections", colour = NULL,
    title = sprintf("[LIGHT] Median prevacc vs postvacc, UF %s", UF_CODE),
    subtitle = sprintf("Routine age %d (%s), VE_sus = %.1f%%",
                       target_age, routine_schedule_label, 100 * ve_sus)
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

p_light_annual_averted <- ggplot(cf_annual_light, aes(x = calendar_year, y = infections_averted)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(fill = "#3182bd") +
  scale_x_continuous(breaks = sort(unique(cf_annual_light$calendar_year))) +
  scale_y_continuous(labels = scales::label_number(), expand = expansion(mult = c(0.05, 0.10))) +
  labs(
    x = NULL, y = "Infections averted",
    title = sprintf("[LIGHT] Annual infections averted, age-%d routine", target_age),
    subtitle = "Median-parameter run (baseline minus routine)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p_light_susceptible <- ggplot(
  cf_weekly_light,
  aes(x = week_start, y = effective_susceptible_frac, colour = scenario)
) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0, NA)) +
  labs(
    x = NULL, y = "Effective susceptible fraction", colour = NULL,
    title = sprintf("[LIGHT] Effective susceptible fraction, UF %s", UF_CODE),
    subtitle = "Vaccinated susceptibles weighted by 1 - VE_sus"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

print(p_light_weekly)
print(p_light_annual_averted)
print(p_light_susceptible)
