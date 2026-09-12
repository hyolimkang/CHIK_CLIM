# Ceará renewal-model development log — v4.1

See `../17_v4_0_minimal_no_vaccine/MODEL_DEVELOPMENT.md` for v4.0 (base run
unconverged at max_treedepth=12; td14 control run **PASSED** the strict HMC
gate, confirming the backbone is numerically sound).

## v4.1 — orthonormal annual effects + estimated constant reporting fraction

**Status:** IN PROGRESS (primary fit running as of 2026-09-11 14:00).

### Why this version was introduced

td14 passing HMC removed the computational-convergence problem, but v4.0 is
not scientifically adequate as a final no-vaccination baseline because: (1)
`q` is fixed exactly at 0.10 with zero uncertainty, so every downstream
susceptibility number is conditional on an arbitrary reporting-fraction
scenario; (2) the posterior-predictive fit remains visibly imperfect for some
epidemic years, particularly 2022; (3) the annual-effect parameterisation
(`z_year - mean(z_year)`) retains one likelihood-redundant direction (v4.0
audit item 12: treedepth-12 hits were not confined to extreme parameter
values, consistent with a mild funnel from this redundancy).

v4.1 makes exactly two substantive changes to the HMC-passing td14 backbone
-- nothing else:

**A. Annual-effect parameterisation.** Replaces `z_year - mean(z_year))` with
an explicit `Y x (Y-1)` orthonormal sum-to-zero basis `B_year` (normalised
Helmert contrasts, built in R and passed as `data`; numerically verified:
`max|t(B)B - I| = 2.2e-16`, `max|colSums(B)| = 1.1e-16` for Y=11). Only the
Y-1 identified degrees of freedom (`z_year_free`) are sampled, and a single
learned scale `sigma_year ~ half-normal(0, 0.35)` (via `real<lower=0>`
truncation) replaces v4.0's externally fixed `year_effect_prior_sd = 0.40`.
`sum(year_effect)` is checked to be ~0 for every posterior draw
(`04_check_orthonormal_year_basis.R`).

**B. Reporting fraction.** Replaces `q_fixed = 0.10` (data, zero uncertainty)
with one estimated constant `real<lower=0,upper=1> q` and an informative
prior. External epidemiological information is `p_symptomatic ~ Beta(30,28)`,
`p_detect_given_symptomatic ~ Beta(20,60)`; the case likelihood identifies
mainly their product, so `p_symptomatic` and `p_detect_given_symptomatic` are
**not** estimated separately (that would add a non-identifiable direction).
Instead `q ~ Beta(16.2, 108.7)` is used directly, after a Monte Carlo (n=2e6)
+ Kolmogorov-Smirnov check confirmed it adequately approximates the induced
product distribution:

| distribution | mean | sd | median | q025 | q975 |
|---|---|---|---|---|---|
| induced product (MC, n=2e6) | 0.1293 | 0.0299 | 0.1271 | 0.0773 | 0.1938 |
| Beta(16.2, 108.7) (previously suggested) | 0.1297 | 0.0300 | 0.1278 | 0.0770 | 0.1938 |
| Beta(16.18, 109.0) (moment-matched, for comparison) | 0.1293 | 0.0299 | 0.1273 | 0.0766 | 0.1932 |

KS distance of Beta(16.2,108.7) vs the induced product: **0.0096** (below the
prespecified 0.02 adequacy threshold) → the previously suggested prior is
used as-is; no re-derivation was needed. Full table:
`03_Output/tables/renewal_v4_1_minimal_no_vaccine/q_prior_monte_carlo_check.csv`
(script: `00_check_q_prior_approximation.R`).

`q` is not weekly, not yearly, not monotonic in time -- one constant for the
whole 2015-2025 series.

### Explicitly unchanged from v4.0/td14

Generation interval; susceptible/immune bookkeeping and the lifelong
infection-derived-immunity assumption; demographic accounting (independent of
q); one annual sine/cosine seasonal cycle (no added harmonics); fixed
`imports_per_week = 1`; negative-binomial observation model; all-age
structure; no vaccination; no climate coefficient; no serology likelihood;
no age structure; no weekly latent R process or AR(1).

### Sampler-run record

- **base/primary (2026-09-11):** 4 chains, 2,000 iterations (1,000
  warm-up), `adapt_delta = 0.95`, `max_treedepth = 14`, `q ~ Beta(16.2,
  108.7)`. *(diagnostics pending -- fit in progress)*
- **q_prior_narrow (planned, only if primary passes):** identical
  model/data/init; `q ~ Beta(32.4, 217.4)` (primary (a,b) × 2 -- half the
  prior variance, same prior mean 0.1293).
- **q_prior_wide (planned, only if primary passes):** identical
  model/data/init; `q ~ Beta(8.1, 54.35)` (primary (a,b) × 0.5 -- double the
  prior variance, same prior mean).

### Primary questions this fit must answer

1. Does v4.1 still pass the strict HMC gate after estimating q?
2. Does the posterior update q meaningfully from its epidemiological prior?
3. Does uncertainty in q substantially widen uncertainty in S(t)?
4. Does the model improve the 2022 posterior-predictive mismatch?
5. Is any improvement attributable mainly to q, to sigma_year, or both?
6. Are q, alpha_R, sigma_year, and S_2025 strongly confounded?

*(Answers pending -- see `03_Output/tables/renewal_v4_1_minimal_no_vaccine/`
once diagnostics have been run: `q_and_sigma_year_prior_vs_posterior.csv`,
`q_posterior_correlations.csv`, `annual_ppc_highlight_years.csv`,
`annual_mismatch_comparison_v4_0_vs_v4_1.csv`.)*

### Prespecified acceptance checks (unchanged from v4.0)

divergences = 0; max-treedepth hits = 0; max Rhat <= 1.01; min bulk ESS >=
100; BFMI >= 0.30 each chain. Also reported: min tail ESS, treedepth
quantiles, chain-specific diagnostics
(`02_diagnose_v4_1_minimal_no_vaccine.R`).

### Simulation recovery (mandatory before calling v4.1 stable)

`03_simulation_recovery_v4_1_minimal_no_vaccine.R` simulates a synthetic
Ceará trajectory from known `q` (= prior mean), `alpha_R`, `sigma_year`,
`beta_sin`, `beta_cos`, `phi_obs`, refits, and checks 95% CrI coverage of
those six parameters plus `S_2022`, `S_2025`, and cumulative infections
(2015-2025). *(Pending.)*

### Decision rule

`06_v4_1_decision_rule.R` reads the outputs above and labels v4.1 a
**stabilized historical no-vaccination backbone** only if all of: (1) HMC
passes; (2) recovery is acceptable; (3) q is not pathologically confounded
(|r| <= 0.9) with alpha_R, sigma_year, S_2025, R0_2017, or R0_2022; (4) 2017
and 2022 are reasonably reproduced (|relative residual| < 50%); (5) the
q-prior sensitivity comparison does not show S/cumulative-infections/R0
swinging by more than 25% across the narrow/primary/wide priors. If HMC
passes but q/S are essentially unmoved from their prior and sensitivity shows
S is prior-determined, the verdict is reported explicitly as *"computationally
stable but absolute susceptibility weakly identified"* rather than a
pass/fail label. If the 2022 mismatch persists despite learned sigma_year and
uncertain q, the prespecified response is to stop and review remaining
structural assumptions -- **not** to add time-varying q to chase visual fit.

*(Verdict pending all of the above.)*

### Next version (not implemented here)

v4.2+ extensions (climate, serology likelihood, age structure, vaccination,
time-varying q, weekly latent R) are deferred pending this decision-rule
verdict and separate scientific review.
