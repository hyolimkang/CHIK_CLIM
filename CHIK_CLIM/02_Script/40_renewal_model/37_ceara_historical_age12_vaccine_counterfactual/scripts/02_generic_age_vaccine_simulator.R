# Generic age-cohort + vaccination simulator with a FEEDBACK SWITCH, used
# for all three arms (A/B/C) via scenario parameters -- the renewal engine
# itself is written ONCE and never copy-pasted between arms.
#
# This is a NEW file in the isolated analysis directory; it does NOT
# modify 36_climate_forced_v4_9/scripts/07_age_cohort_simulator.R. The
# transmission/demographic/ageing/vaccination mechanics are identical to
# that accepted, validated simulator (same exact-identity design); the
# only addition is `feedback_on`:
#
#   feedback_on = TRUE  (Arms A, B): infectiousness/lambda are generated
#     from THIS run's OWN recursively-built X_total history (closed loop).
#   feedback_on = FALSE (Arm C, diagnostic only): lambda[t] is taken
#     directly from a supplied EXTERNAL trajectory (Arm A's lambda), so
#     vaccination can reduce susceptibles and hence infections WITHOUT
#     those reduced infections feeding back into future transmission --
#     isolating the direct-only component of the vaccine effect. This is
#     NOT a physical scenario; it exists purely for the decomposition in
#     Section 5 of the task specification.

AGEING_FRACTION <- 1 / 52.1775

#' @param R0_t_draw numeric[N] fitted R0(t) (already includes climate multiplier)
#' @param feedback_on logical: TRUE = closed-loop (own X_total history drives lambda); FALSE = lambda taken from `lambda_external`
#' @param lambda_external numeric[N] or NULL: required when feedback_on=FALSE (e.g. Arm A's lambda, for Arm C)
#' @param vaccination NULL or list(target_age, coverage_by_year OR coverage_by_week, ve_infection)
simulate_age_vaccine <- function(R0_t_draw, births, deaths, reconciliation, N_start, w, imports_per_week,
                                  is_seed, X_seed, age_shares, years, max_age = 100L,
                                  vaccination = NULL, feedback_on = TRUE, lambda_external = NULL) {
  N <- length(R0_t_draw)
  A <- max_age + 1L
  G <- length(w)
  if (!feedback_on && is.null(lambda_external)) stop("feedback_on=FALSE requires lambda_external.")
  if (!is.null(vaccination)) stopifnot(xor(!is.null(vaccination$coverage_by_year), !is.null(vaccination$coverage_by_week)))

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

  share1 <- get_share_vec(years[1])
  S[1, ] <- share1 * N_start[1]; Uinf[1, ] <- 0; Uvac[1, ] <- 0

  for (t in seq_len(N)) {
    total_t <- S[t, ] + Uinf[t, ] + Uvac[t, ]
    N_t <- sum(total_t)
    stopifnot(abs(N_t - N_start[t]) < 1e-6 * max(1, N_start[t]))

    if (feedback_on) {
      infectiousness <- 0
      if (t > 1) for (g in seq_len(min(G, t - 1))) infectiousness <- infectiousness + w[g] * X_total[t - g]
      lambda_t <- R0_t_draw[t] * infectiousness / N_t + imports_per_week / N_t
    } else {
      lambda_t <- lambda_external[t] # Arm C: fixed to the no-vaccine trajectory, no feedback
      infectiousness <- NA_real_ # not used/meaningful in this mode; not part of this arm's own recursion
    }
    hazard <- -expm1(-lambda_t)
    infectiousness_vec[t] <- infectiousness; lambda_vec[t] <- lambda_t; hazard_vec[t] <- hazard
    Reff_vec[t] <- R0_t_draw[t] * sum(S[t, ]) / N_t

    if (is_seed[t] == 1) {
      S_share <- if (sum(S[t, ]) > 0) S[t, ] / sum(S[t, ]) else rep(1 / A, A)
      Inew <- pmin(S_share * X_seed[t], S[t, ])
    } else {
      Inew <- S[t, ] * hazard
    }
    X_total[t] <- sum(Inew)
    Inew_age[t, ] <- Inew

    S_after_inf <- S[t, ] - Inew
    Uinf_after_inf <- Uinf[t, ] + Inew
    Uvac_after_inf <- Uvac[t, ]

    if (t == N) break

    S_tot <- sum(S_after_inf); U_tot <- sum(Uinf_after_inf) + sum(Uvac_after_inf)
    sfrac <- S_tot / (S_tot + U_tot)
    share_t <- get_share_vec(years[t])

    deaths_to_S <- deaths[t] * share_t * sfrac
    deaths_to_Uinf <- deaths[t] * share_t * (1 - sfrac) * (sum(Uinf_after_inf) / (U_tot + 1e-12))
    deaths_to_Uvac <- deaths[t] * share_t * (1 - sfrac) * (sum(Uvac_after_inf) / (U_tot + 1e-12))
    recon_to_S <- reconciliation[t] * share_t * sfrac
    recon_to_Uinf <- reconciliation[t] * share_t * (1 - sfrac) * (sum(Uinf_after_inf) / (U_tot + 1e-12))
    recon_to_Uvac <- reconciliation[t] * share_t * (1 - sfrac) * (sum(Uvac_after_inf) / (U_tot + 1e-12))

    S_after_demo <- S_after_inf - deaths_to_S + recon_to_S
    Uinf_after_demo <- Uinf_after_inf - deaths_to_Uinf + recon_to_Uinf
    Uvac_after_demo <- Uvac_after_inf - deaths_to_Uvac + recon_to_Uvac

    inflow_S <- S_after_demo[-A] * AGEING_FRACTION
    inflow_Uinf <- Uinf_after_demo[-A] * AGEING_FRACTION
    inflow_Uvac <- Uvac_after_demo[-A] * AGEING_FRACTION

    if (!is.null(vaccination)) {
      target_age <- vaccination$target_age
      cov <- coverage_at(t + 1, years[t + 1])
      if (cov > 0) {
        src_idx <- target_age
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
      out[A] <- out[A] + v[A] * AGEING_FRACTION
      out
    }
    S_aged <- age_step_apply(S_after_demo, inflow_S)
    Uinf_aged <- age_step_apply(Uinf_after_demo, inflow_Uinf)
    Uvac_aged <- age_step_apply(Uvac_after_demo, inflow_Uvac)
    S_aged[1] <- S_aged[1] + births[t]

    S[t + 1, ] <- S_aged; Uinf[t + 1, ] <- Uinf_aged; Uvac[t + 1, ] <- Uvac_aged
  }

  list(S_age = S, Uinf_age = Uinf, Uvac_age = Uvac, Inew_age = Inew_age, X_total = X_total,
       S_total = rowSums(S), Uinf_total = rowSums(Uinf), Uvac_total = rowSums(Uvac),
       N_total = rowSums(S) + rowSums(Uinf) + rowSums(Uvac),
       hazard = hazard_vec, infectiousness = infectiousness_vec, lambda = lambda_vec, R_eff = Reff_vec,
       doses_administered = doses_administered, effective_protected = effective_protected,
       entrants_to_target_total = entrants_to_target_total, entrants_to_target_susceptible = entrants_to_target_susceptible,
       feedback_on = feedback_on)
}
