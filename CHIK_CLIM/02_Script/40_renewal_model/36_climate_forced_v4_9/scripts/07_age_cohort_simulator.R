# Age-structured vaccine extension -- parsimonious age-cohort bookkeeping
# layered OVER the fitted climate-forced v4.9 process. NO age-specific
# transmission: every age experiences the SAME weekly force of infection
# lambda[t]. This is a CLOSED-LOOP simulator: infectiousness[t] is built
# from THIS function's OWN recursively-generated X_total history -- the
# original fit's X(t) is never read after week 1's initial condition
# (S[1]=N_start[1] fully susceptible, matching v4.9). Only R0_t[t]
# (already includes the fitted climate multiplier), the generation-
# interval weights, imports, demography, and the 2022 is_seed/X_seed
# conditioning are taken from the fit.
#
# Exact-identity design (see PHASE0_PHASE1_AGE_VACCINE_AUDIT.md): since
# (a) the hazard (1-exp(-lambda[t])) is applied identically to every age,
# (b) births/deaths/reconciliation are allocated across ages purely by
# population SHARE (using the same scalar susceptible-fraction sfrac as
# the aggregate v4.9 recursion), and (c) ageing is a zero-sum
# redistribution across ages -- summing the age-resolved recursion over
# age reproduces the aggregate v4.9 S[t]/U[t]/X[t] recursion term-by-term.
#
# Vaccination (routine age-12 mechanism): acts ONLY on the weekly
# age-11 -> age-12 AGEING FLOW (never on the standing age-12 stock), so
# each birth cohort is offered vaccination exactly once, the week it
# crosses into age 12. `doses_administered` counts ALL age-12 entrants at
# the given coverage (regardless of prior infection status, matching a
# real programme that does not screen for infection history);
# `effective_protected` counts only the susceptible-entrant subset who
# both receive the dose AND for whom it takes (coverage * VE_infection).
# Only the effective-protected mass moves S -> Uvac; vaccinated-but-not-
# protected individuals remain susceptible; already-immune entrants
# remain Uinf untouched.

AGEING_FRACTION <- 1 / 52.1775 # weekly fraction of each age cohort advancing to age+1

#' Simulate the age-cohort model for ONE posterior draw (closed-loop).
#'
#' @param R0_t_draw numeric[N] -- EXACT fitted R0(t) for this draw (already includes climate_multiplier)
#' @param births,deaths,reconciliation,N_start numeric[N] -- from stan_data/weekly_data (age-agnostic, exogenous)
#' @param w numeric[G] generation-interval weights (sums to 1)
#' @param imports_per_week scalar
#' @param is_seed integer[N], X_seed numeric[N] -- conditioned 2022 seed window, applied to X_total exactly as in v4.9
#' @param age_shares data.frame(year, age, share) -- self-normalised age shares per calendar year
#' @param years integer[N] calendar year per week (for age-share lookup)
#' @param max_age integer
#' @param vaccination NULL, or list(target_age, coverage_by_year [named by as.character(year)], ve_infection),
#'   or list(target_age, coverage_by_week [numeric[N], one-time-pulse-capable], ve_infection).
#'   Exactly one of coverage_by_year / coverage_by_week must be supplied when non-NULL.
#' @return list with age-resolved matrices, aggregate vectors, and closed-loop diagnostics
#'   (infectiousness, force_of_infection, R_eff -- all reconstructed from the simulator's OWN state,
#'   using the same definitions as the fitted v4.9 model), plus vaccination bookkeeping vectors.
simulate_age_cohort <- function(R0_t_draw, births, deaths, reconciliation, N_start, w, imports_per_week,
                                 is_seed, X_seed, age_shares, years, max_age = 100L, vaccination = NULL) {
  N <- length(R0_t_draw)
  A <- max_age + 1L # ages 0..max_age
  G <- length(w)
  if (!is.null(vaccination)) {
    stopifnot(xor(!is.null(vaccination$coverage_by_year), !is.null(vaccination$coverage_by_week)))
  }

  S <- matrix(0, N, A); Uinf <- matrix(0, N, A); Uvac <- matrix(0, N, A)
  Inew_age <- matrix(0, N, A)
  X_total <- numeric(N); hazard_vec <- numeric(N)
  infectiousness_vec <- numeric(N); lambda_vec <- numeric(N); Reff_vec <- numeric(N)
  doses_administered <- numeric(N); effective_protected <- numeric(N)
  entrants_to_target_total <- numeric(N); entrants_to_target_susceptible <- numeric(N)

  get_share_vec <- function(yr) {
    s <- age_shares$share[age_shares$year == yr]
    if (length(s) != A) stop("Age share vector length mismatch for year ", yr)
    s
  }
  coverage_at <- function(t_next, year_next) {
    if (is.null(vaccination)) return(0)
    cov <- if (!is.null(vaccination$coverage_by_year)) vaccination$coverage_by_year[as.character(year_next)] else vaccination$coverage_by_week[t_next]
    if (is.na(cov)) 0 else cov
  }

  # ---- initialise week 1: fully susceptible, distributed by age share (matches v4.9's S[1]=N_start[1], U[1]=0) ----
  share1 <- get_share_vec(years[1])
  S[1, ] <- share1 * N_start[1]
  Uinf[1, ] <- 0
  Uvac[1, ] <- 0

  for (t in seq_len(N)) {
    total_t <- S[t, ] + Uinf[t, ] + Uvac[t, ]
    N_t <- sum(total_t)
    stopifnot(abs(N_t - N_start[t]) < 1e-6 * max(1, N_start[t]))

    infectiousness <- 0
    if (t > 1) {
      for (g in seq_len(min(G, t - 1))) infectiousness <- infectiousness + w[g] * X_total[t - g]
    }
    lambda_t <- R0_t_draw[t] * infectiousness / N_t + imports_per_week / N_t
    hazard <- -expm1(-lambda_t) # 1 - exp(-lambda_t), same for every age (NO age-specific FOI)
    infectiousness_vec[t] <- infectiousness; lambda_vec[t] <- lambda_t; hazard_vec[t] <- hazard
    Reff_vec[t] <- R0_t_draw[t] * sum(S[t, ]) / N_t # SAME definition as v4.9's R_eff_t = R0_t * S/N

    if (is_seed[t] == 1) {
      # Conditioned seed week: distribute X_seed[t] across ages in proportion
      # to each age's CURRENT susceptible share (same principle as the
      # hazard-based step -- a uniform, non-age-specific allocation rule).
      # The ordinary lambda-based hazard does NOT determine total infections
      # this week; X_total[t] is forced to X_seed[t] exactly.
      S_share <- if (sum(S[t, ]) > 0) S[t, ] / sum(S[t, ]) else rep(1 / A, A)
      Inew <- pmin(S_share * X_seed[t], S[t, ])
    } else {
      Inew <- S[t, ] * hazard
    }
    X_total[t] <- sum(Inew)
    Inew_age[t, ] <- Inew

    S_after_inf <- S[t, ] - Inew
    Uinf_after_inf <- Uinf[t, ] + Inew
    Uvac_after_inf <- Uvac[t, ] # infection does not touch Uvac (lifelong protection assumed)

    if (t == N) break

    # ---- demographic step: allocate deaths/reconciliation by AGE SHARE, using
    # the SAME scalar susceptible-fraction split as the aggregate v4.9 model ----
    S_tot <- sum(S_after_inf); U_tot <- sum(Uinf_after_inf) + sum(Uvac_after_inf)
    sfrac <- S_tot / (S_tot + U_tot)
    share_t <- get_share_vec(years[t]) # allocate by THIS week's age composition

    deaths_to_S <- deaths[t] * share_t * sfrac
    deaths_to_Uinf <- deaths[t] * share_t * (1 - sfrac) * (sum(Uinf_after_inf) / (U_tot + 1e-12))
    deaths_to_Uvac <- deaths[t] * share_t * (1 - sfrac) * (sum(Uvac_after_inf) / (U_tot + 1e-12))
    recon_to_S <- reconciliation[t] * share_t * sfrac
    recon_to_Uinf <- reconciliation[t] * share_t * (1 - sfrac) * (sum(Uinf_after_inf) / (U_tot + 1e-12))
    recon_to_Uvac <- reconciliation[t] * share_t * (1 - sfrac) * (sum(Uvac_after_inf) / (U_tot + 1e-12))

    S_after_demo <- S_after_inf - deaths_to_S + recon_to_S
    Uinf_after_demo <- Uinf_after_inf - deaths_to_Uinf + recon_to_Uinf
    Uvac_after_demo <- Uvac_after_inf - deaths_to_Uvac + recon_to_Uvac

    # ---- ageing: zero-sum redistribution (fraction AGEING_FRACTION moves a -> a+1; top age retains inflow).
    # Vaccination intercepts ONLY the age(target-1) -> age(target) inflow slice, not the standing stock. ----
    inflow_S <- S_after_demo[-A] * AGEING_FRACTION   # inflow[k] = flow from age (k-1) into age k, k=1..A-1 (0-indexed age k)
    inflow_Uinf <- Uinf_after_demo[-A] * AGEING_FRACTION
    inflow_Uvac <- Uvac_after_demo[-A] * AGEING_FRACTION

    if (!is.null(vaccination)) {
      target_age <- vaccination$target_age
      cov <- coverage_at(t + 1, years[t + 1])
      if (cov > 0) {
        src_idx <- target_age # inflow_S[target_age] is the flow FROM age (target_age-1) INTO age target_age (1-indexed inflow vector: inflow_S[k] = flow into age k)
        entrant_S <- inflow_S[src_idx]; entrant_Uinf <- inflow_Uinf[src_idx]; entrant_Uvac <- inflow_Uvac[src_idx]
        entrants_total <- entrant_S + entrant_Uinf + entrant_Uvac
        doses <- cov * entrants_total
        effective <- cov * vaccination$ve_infection * entrant_S
        inflow_S[src_idx] <- inflow_S[src_idx] - effective
        inflow_Uvac[src_idx] <- inflow_Uvac[src_idx] + effective
        doses_administered[t + 1] <- doses; effective_protected[t + 1] <- effective
        entrants_to_target_total[t + 1] <- entrants_total; entrants_to_target_susceptible[t + 1] <- entrant_S
      }
    }

    age_step_apply <- function(v, inflow) {
      out <- v * (1 - AGEING_FRACTION)
      out[-1] <- out[-1] + inflow
      out[A] <- out[A] + v[A] * AGEING_FRACTION # open-ended top group: no outflow beyond max_age
      out
    }
    S_aged <- age_step_apply(S_after_demo, inflow_S)
    Uinf_aged <- age_step_apply(Uinf_after_demo, inflow_Uinf)
    Uvac_aged <- age_step_apply(Uvac_after_demo, inflow_Uvac)

    # ---- births enter age-0 susceptible AFTER ageing ----
    S_aged[1] <- S_aged[1] + births[t]

    S[t + 1, ] <- S_aged; Uinf[t + 1, ] <- Uinf_aged; Uvac[t + 1, ] <- Uvac_aged
  }

  list(S_age = S, Uinf_age = Uinf, Uvac_age = Uvac, Inew_age = Inew_age, X_total = X_total,
       S_total = rowSums(S), Uinf_total = rowSums(Uinf), Uvac_total = rowSums(Uvac),
       hazard = hazard_vec, infectiousness = infectiousness_vec, lambda = lambda_vec, R_eff = Reff_vec,
       doses_administered = doses_administered, effective_protected = effective_protected,
       entrants_to_target_total = entrants_to_target_total, entrants_to_target_susceptible = entrants_to_target_susceptible)
}
