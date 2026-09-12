# Ceará renewal-model development log

This log records model versions and decisions.  A computationally successful
run is not, by itself, evidence that a model is scientifically adequate.

## v4.0 — minimal all-age no-vaccination benchmark

**Status:** the base sampler run (max_treedepth=12) did not meet the strict
HMC gate and remains preserved, unmodified, as an explicitly unconverged
diagnostic run in `renewal_v4_0_minimal_no_vaccine/`. The td14 control run
(identical model/data/priors/init, max_treedepth raised to 14 only) **PASSED**
the same gate in full: divergences=0, max-treedepth hits=0, max Rhat=1.00897,
min bulk ESS=480.8, min tail ESS=690.5, BFMI in [0.926, 1.083] across all 4
chains. This confirms the base run's 64 treedepth-12 hits (51/0/5/8 by chain)
were a numerical ceiling, not a geometry breakdown requiring reparameterisation
on their own -- see the audit below. Results are in
`renewal_v4_0_minimal_no_vaccine_td14/`.

### Why this version was introduced

The preceding v2.3 state model combined a weekly latent AR(1) transmission
process, seasonal terms, a fitted reporting fraction, fitted seed hazards,
and a tempered serology likelihood.  Its saved diagnostic file reported no
divergences but 4,000 maximum-treedepth transitions, a maximum R-hat of
1.0373, and 161 monitored parameters with R-hat above 1.01.  It is therefore
not an adequate posterior basis for a vaccine counterfactual.

The purpose of v4.0 is deliberately narrower: establish a stable all-age,
no-vaccination transmission baseline for Ceará, 2015-01-04 to 2025-12-21.
It is a retrospective benchmark, not a prospective outbreak forecast and
not yet a vaccination-impact model.

### Sampler-run record

- **base (2026-09-11):** 4 chains, 2,000 iterations (1,000 warm-up),
  `adapt_delta = 0.95`, maximum treedepth 12. It had zero divergences,
  maximum R-hat 1.00992, and minimum bulk ESS 424.96, but 64 maximum-
  treedepth hits. It is retained as an explicitly unconverged diagnostic run.
- **td14 (2026-09-11):** identical likelihood, priors, data, and
  initialisation; only the numerical NUTS cap raised from 12 to 14. **PASS**:
  divergences=0, max-treedepth hits=0, max Rhat=1.00897, min bulk ESS=480.8,
  min tail ESS=690.5, BFMI all >0.30. Saved under a distinct `*_td14` fit and
  output folder so it cannot overwrite the base run. This confirms v4.0 has a
  computationally acceptable numerical implementation, but does not by itself
  make v4.0 the final baseline (see the audit and v4.1 entry below).

### Audit (2026-09-11, before any v4.1 change)

Twelve prespecified checks were run against the base fit and its inputs
before touching any code: (1) generation-interval weights sum to 1 exactly;
(2) `year_id` correctly maps all 573 weeks to 11 calendar years (52/53-week
years correctly resolved); (3) `seasonal_sin`/`seasonal_cos` use a period of
52.1775 weeks = 365.2425 days = exactly one solar year, spanning 10.98
cycles over the full series; (4) `N_end = N_start + births - deaths +
reconciliation` holds with 0 error; (5) `S + U = N_start` holds to 5.2e-08
(floating-point) across all draws and weeks; (6)-(7) `q_fixed` and
`imports_per_week` are plain `data`, never declared in Stan's `parameters`
block; (8) `year_effect = year_effect_prior_sd * (z_year - mean(z_year))` is
reproduced exactly (max abs diff 1.1e-16) by re-deriving it in R from the
extracted `z_year` draws; (9) the Stan source contains no weekly AR process,
climate coefficient, serology likelihood, or fitted seed parameter --
confirmed by direct source read; (10) the posterior-predictive `C_pred` and
the likelihood both read from the same `expected_reported_cases` vector --
confirmed by source read; (11) the `_UNCONVERGED`-suffixed six-panel PDF's
timestamp (12:58) is consistent with the base fit bundle's save time (12:54),
not a stale cached object; (12) the 64 treedepth-12 hits are NOT confined to
extreme marginal values of `alpha_R`/`beta_sin`/`beta_cos`/`phi_obs` -- e.g.
chain 1's 51 hits have `alpha_R` in [0.112, 0.182] vs an overall chain range
of [0.098, 0.244] -- consistent with a mild funnel from the one
likelihood-redundant direction in `z_year - mean(z_year)` rather than a
single extreme parameter regime. All twelve checks passed; this motivated
fixing that redundant direction in v4.1 (below) even though td14 already
passes the gate.

### Formulation

For week `t`, the latent infections are

`X_t = S_t [1 - exp{- (R0_t * sum_g w_g X_(t-g) / N_t + m / N_t)}]`.

`R0_t` has only an intercept, a centred annual deviation (11 calendar-year
levels, fixed prior scale), and sine/cosine annual seasonality.  Susceptible
and infection-derived immune stocks evolve deterministically using weekly
births, deaths, and population reconciliation. Observed reported cases are
negative-binomial draws with mean `q_fixed * X_t`.

### Explicit restrictions that reduce non-identifiability

- No weekly AR(1) latent process or process-noise parameter.
- No fitted reporting fraction: `q_fixed = 0.10` is an observation-scale
  scenario input only. It is not used to calibrate the spatial susceptibility
  pipeline.
- No fitted seed hazard: `imports_per_week = 1` is fixed external infection
  pressure.
- No climate coefficient, no serology likelihood, and no local serology
  offsets. The available municipality surveys are shown only as an external
  posterior-predictive consistency check.
- No vaccination mechanism or age stratification.

These restrictions are modelling choices, not claims that reporting,
importation, climate, or age are unimportant. They are deferred until this
baseline satisfies convergence and recovery checks.

### Prespecified acceptance checks

The fit is labelled **PASS** only if: divergences = 0; maximum-treedepth
hits = 0; maximum R-hat <= 1.01; minimum bulk ESS >= 100; and each chain has
BFMI >= 0.30.  Outputs also include weekly posterior-predictive coverage,
annual case comparison, external serology consistency, and an independently
generated simulation-recovery exercise. Failed gates must remain labelled
unconverged and must not be used for vaccine-impact inference.

### Next version

See `18_v4_1_minimal_no_vaccine/MODEL_DEVELOPMENT.md` for v4.1, which changes
exactly two things relative to this HMC-passing td14 backbone: the annual-
effect parameterisation (explicit orthonormal sum-to-zero basis, learned
`sigma_year`) and the reporting fraction (estimated constant `q` with an
informative prior, replacing `q_fixed = 0.10`). Only after v4.1 (or a later
minimal revision) is confirmed a stabilized backbone should a separate
vaccination scenario layer be added: retaining all-age transmission and
making age-specific vaccination an externally specified removal from the
susceptible stock, with age allocation tracked for impact attribution rather
than attempting to identify fully age-specific transmission from these data.
