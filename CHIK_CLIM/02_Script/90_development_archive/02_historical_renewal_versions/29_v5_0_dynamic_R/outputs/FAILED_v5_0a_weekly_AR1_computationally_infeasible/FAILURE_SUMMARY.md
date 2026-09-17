# v5.0a (weekly AR(1) dynamic-R residual) — archived as computationally infeasible

**Classification: computationally infeasible parameterisation.**
**NOT evidence that short-timescale R0 variation is scientifically unsupported.**

## What was tried

`renewal_ceara_v5_0a_dynamic_R_direct_serology.stan`: v4.9's full harmonic +
A_year backbone (`mu_R[t]`), plus `log_R0[t] = mu_R[t] + delta_R[t]` where
`delta_R` is a non-centred stationary AR(1) process over all N=573 weeks
(one latent parameter per week, `z_delta[1..573]`), re-centred to zero mean
within each calendar year. Priors: `rho_R ~ beta(8,2)`, `sigma_R_dynamic ~
normal(0,0.15)` (half-normal), both pre-registered, not tuned. Fit at
q=0.05, same HMC configuration as v4.9 (dense_e, adapt_delta=0.95,
max_treedepth=14, 4 chains, 800 warmup / 400 sampling planned).

## What happened

After ~20 minutes, only ~50 of 1200 planned iterations had completed
across all 4 chains (partial sample files preserved in
`v5_0a_q0.05_partial_pilot/chains/`):

| chain | iterations reached | at max treedepth (14) | divergent |
|---|---|---|---|
| 1 | 52 | 44 (85%) | 4 |
| 2 | 51 | 43 (84%) | 5 |
| 3 | 52 | 43 (83%) | 7 |
| 4 | 50 | 43 (86%) | 5 |

~85% of transitions used the maximum 16,383 leapfrog steps, and real
divergences appeared even during early warmup (not merely slow — a
genuine geometric pathology). Projected full-run time: many hours.

## Why (diagnosis, not a fix attempted)

Most likely the anticipated identifiability tension between `A_year` and
`delta_R` (flagged pre-registration in the v5.0a spec): re-centering
`delta_R` to zero mean within each year removes the level but not the
*shape* of the AR(1) trajectory, and a highly persistent process
(`rho_R` prior mean 0.8) can still trace slow, smooth within-year curves
that overlap substantially with the harmonic+A_year term's own shape,
producing a funnel-like joint posterior across 573 sequentially-correlated
weekly latent states plus the existing hierarchical structure.

## Disposition

Per explicit instruction: do NOT increase max_treedepth, do NOT lower
adapt_delta to force speed, do NOT wait out the full run. The run was
stopped (`TaskStop`, verified no orphan `Rscript.exe` processes) and this
partial pilot preserved as evidence. Superseded by
`renewal_ceara_v5_0b_lowrank_dynamic_R.stan`, which tests the same
scientific hypothesis (limited short-timescale R0 deviation from the v4.9
backbone) with a low-dimensional (K=8 cyclic basis coefficients per year,
Y×K=88 parameters total) rather than a full weekly-AR(1) (N=573
parameters) representation.
