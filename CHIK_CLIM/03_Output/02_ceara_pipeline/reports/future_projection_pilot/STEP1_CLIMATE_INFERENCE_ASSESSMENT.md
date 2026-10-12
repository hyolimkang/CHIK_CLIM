# STEP 1 -- Climate-transmission inference audit (Ceara, v4.9 climate-forced canary, q=0.05)

Produced by `02_Script/02_ceara_pipeline/10_future_projection_pilot/01_step1_climate_posterior_reanalysis.R`
and `02_step1_figures.R`. Re-analyses the **already-fitted, accepted**
posteriors (`ce_climate_forced_canary_q0.05.rds`, `renewal_ceara_v4_9_fit_q0.05.rds`);
no Stan model was refit to produce this report. Extends, and does not
duplicate, the existing `03_Output/02_ceara_pipeline/reports/climate_forced_v4_9/CLIMATE_V49_CANARY_ASSESSMENT.md`.

## 1. Model reconstruction (Section 1.1)

The implemented expression, extracted directly from `renewal_ceara_v4_9_climate_forced_canary.stan`
and verified numerically (not just read) against the saved posterior draws:

```
log_R0(t) = alpha_R + year_effect[year_id(t)]
          + A_year[year_id(t)] * (beta_sin1*sin1(t) + beta_cos1*cos1(t))
          + beta_sin2*sin2(t) + beta_cos2*cos2(t)
          + beta_T * z_T_anom(t) + beta_P * z_P_anom(t)

year_effect[y] = 0.40 * (z_year[y] - mean(z_year))         [year_effect_prior_sd = 0.40, FIXED]
A_year[y]      = exp(sigma_season_year * z_season_year[y])  [year-varying 1st-harmonic amplitude]
sin1/cos1, sin2/cos2: harmonics of period 52.1775 weeks (1st and 2nd)
beta_T, beta_P ~ Normal(0, 0.15)                             [shared, time-invariant]
```

**Numerical verification**: reconstructing `log_R0(t)` from the extracted
`alpha_R`, `year_effect`, `A_year`, `beta_sin1/cos1/sin2/cos2`, `beta_T`,
`beta_P` and the raw regressor data, and comparing against the `log_R0`
generated quantity Stan itself saved, across 50 randomly sampled posterior
draws x 573 weeks: **max absolute error = 4.4e-16 (floating-point noise) -- PASS**
(`STEP1_equation_reconstruction_check.csv`). The equation above is exactly
what is executed, not an approximation of it.

**Exact specification** (`STEP1_lag_and_standardisation_specification.csv`):
- Generation interval: G=8 weeks, Gamma(shape=4, rate=2), discretised and renormalised.
- Temperature lag: mean of weeks (t-2, t-1, t). Precipitation lag: log1p(sum of weeks t-6..t-2).
- Climatology: mean lagged value by (state=CE, ISO week-of-year), averaged across all available years.
- Standardisation: pooled mean/SD across CE/BA/RJ/MT for the primary lag spec (T_anom_mean=-0.0038, T_anom_sd=1.097, P_anom_mean=0.0126, P_anom_sd=0.756).
- Boundary weeks (2015-01/02) with an incomplete lag window: z-anomaly set to 0 (disclosed, not imputed).
- q (ascertainment) = 0.05, FIXED, not estimated here or anywhere else in the current pipeline.
- Imports: 1/week, fixed, deterministic, added as a hazard term (never double-counted as cases).
- Importation and the 2022 conditioned seed (8 weeks from 2022-01-02, `X_seed = cases/q`, excluded from the likelihood) are unchanged from the non-climate baseline.

These definitions were not changed by this audit.

## 2. Posterior evidence (Section 1.2)

### A-C. beta_T, beta_P posteriors vs. prior N(0, 0.15) -- FIG01

| parameter | posterior median | 95% CrI | P(effect>0) | posterior SD / prior SD |
|---|---|---|---|---|
| beta_T (temperature) | -0.016 | [-0.071, 0.032] | 26% | 0.17 |
| beta_P (precipitation) | 0.024 | [-0.004, 0.052] | 96% | 0.09 |

Both posteriors are far narrower than the prior (SD ratios 0.09-0.17), so
the ~573 weeks of data substantially update the weakly-informative prior --
this is a real likelihood contribution, not prior-dominated shrinkage to
zero. **Temperature**: 95% CrI straddles zero on both sides -- no
detectable linear effect at this lag spec, consistent with the independent
first-epidemic climate GAM's non-significant temperature term. **Precipitation**:
96% of posterior mass is positive, directionally consistent with that GAM's
significant positive association, but the CrI's lower bound sits just
below zero -- suggestive, not decisive, at this single-state scale.

### D. Posterior correlations -- extended beyond the existing assessment

The existing `CLIMATE_V49_CANARY_ASSESSMENT.md` already reported small
correlations of beta_T/beta_P with alpha_R, beta_sin1/cos1,
sigma_season_year, phi_obs and terminal S/immune draws (all |cor| <= 0.21).
This audit extends that table to **each individual year_effect[y]**
(`STEP1_beta_vs_year_effect_per_year_correlations.csv`) and finds a
materially larger correlation that the aggregated hyperparameter view
misses entirely:

| | cor(beta_T, ·) | cor(beta_P, ·) |
|---|---|---|
| beta_T vs beta_P | 0.377 | -- |
| max over single year_effect[y] | **0.803** (year 2022) | -0.611 (year 2019) |
| mean over all 11 year_effect[y] | 0.396 | 0.184 |
| sigma_season_year, mean(A_year), phi_obs, terminal S/immune | all \|cor\| <= 0.21 | all \|cor\| <= 0.21 |

**This is a genuine identifiability finding not visible in the existing
canary assessment.** `beta_T` correlates 0.80 with `year_effect[2022]` --
specifically the one calendar year that also carries the externally
conditioned recurrence seed. Mechanistically this is explicable rather
than pathological: 2022's weekly climate anomalies and that year's
free year-level deviation are both competing to explain the same
unusually large recurrence, so the sampler trades off between them along
a partially shared direction. It does not show up as an HMC pathology
(0 divergences, clean treedepth, healthy ESS) and does not materially move
the aggregate reconstruction (Section 3 below), but it means **beta_T's
marginal posterior should not be read as a "pure", year-independent
temperature effect** -- part of its identification leans on one unusual
year. beta_P's largest single-year correlation (-0.61, with 2019, a quiet
inter-epidemic year) is smaller but non-trivial too. All other structural
correlations (seasonality, observation dispersion, terminal susceptibility)
remain small, confirming the existing assessment's "no material
confounding with the pre-existing harmonic/year structure" for everything
*except* this one specific year-level interaction.

### E-F. Historical climate multiplier -- FIG02

`climate_multiplier(t) = exp(beta_T*z_T_anom(t) + beta_P*z_P_anom(t))`
stays within **[0.930, 1.099]** (posterior median) across 2015-2025; the
widest 95% credible interval half-width at any single week is **0.106**
(i.e. the forcing is a modest, bounded multiplicative adjustment, not an
extreme one, and its uncertainty band never excludes 1 by a wide margin at
any single week).

### G. Sensitivity of inferred S(t) to including climate -- full 573-week series

Extending the existing 5-checkpoint comparison to **every week**: the
climate-forced and non-climate baseline reconstructions of S(t)/N(t) differ
by at most **0.42%** (relative) at any single week, mean **0.096%** across
the whole series (`STEP1_S_prop_full_series_comparison.csv`, FIG03 panel C).
Confirms, at full resolution, that adding climate forcing does not
materially change the inferred long-term susceptibility trajectory.

### H. HMC diagnostics

| | baseline v4.9 | climate-forced canary |
|---|---|---|
| divergences | 0 | 0 |
| max treedepth hits | 0 | 0 |
| max Rhat (classic / rank-normalised) | 1.0008 / 1.0029 | 1.0097 / 1.0098 |
| min ESS bulk / tail (rank-normalised) | 571 / 581 | 470 / 378 |

Both pass the repo's Rhat<=1.01 gate; the climate model sits at the
boundary rather than comfortably under it (reproducing the existing
assessment's finding). **Computational convergence is clean and is not
the limiting factor in this audit** -- the limiting factor is the modest,
partially year-entangled signal described above, which is a
statistical-identifiability and effect-size issue, not a sampler problem.

## 3. Climate vs. seasonality: variance decomposition (Section 1.3)

Per-draw decomposition of `log_R0(t) - alpha_R` into seasonal, year-effect
and climate components (400 posterior draws x 573 weeks;
`STEP1_variance_decomposition_per_draw.csv`, FIG04). Because the three
regressor sets are only approximately mutually orthogonal (Phase-1 check:
max \|cor\| = 0.019 for climate vs. harmonics/year), the variance shares
below do not sum to exactly 100% -- the small residual covariance term is
reported explicitly rather than hidden:

| component | posterior median share of Var(log R0 - alpha_R) | 95% CI |
|---|---|---|
| Seasonal (harmonics) | **80.4%** | [76.8%, 83.6%] |
| Year effect | **19.4%** | [15.7%, 23.0%] |
| Climate (beta_T, beta_P) | **1.1%** | [0.1%, 3.3%] |
| Covariance remainder | -0.7% | [-4.7%, 1.6%] |

This is a **descriptive mathematical decomposition of the fitted model's
own additive structure, not a causal attribution of what drives
transmission**. Seasonality and year-to-year heterogeneity dominate the
historical variation in log R0(t) by construction (they have the most
free parameters and the least regularisation); the climate anomaly terms,
deliberately regularised toward zero (`Normal(0, 0.15)`), contribute a
small but non-zero and non-degenerate share. This is consistent with, and
quantifies, the qualitative picture in Sections 2A-C above.

## 4. Existing comparator and disclosures (Section 1.4)

- **PPC**: unchanged from the existing assessment (weekly 95% coverage
  96.3%->96.9%, weekly RMSE 458.3->455.5); FIG03 panel A shows the annual
  climate-forced curve visually overlapping the baseline curve almost
  everywhere -- climate addition neither improves nor degrades fit
  materially.
- **q is fixed, not estimated** (`STEP1_disclosures_q_and_lag_spec.csv`).
  The existing q in {0.05,...,0.30} sweep is a **q-identifiability check
  for the non-climate v4.9 baseline only** -- it contains no climate
  forcing whatsoever and must not be read as a sensitivity analysis of
  beta_T/beta_P. No q-sensitivity sweep exists for the climate-forced
  model; this is a genuine gap, not something resolved here.
- **Lag-specification sensitivity is prepared but NOT fit.** Two alternate
  lag specs ("sensitivity_A": T mean(t-1,t), P sum(t-4..t-1);
  "sensitivity_B": T mean(t-3..t-1), P sum(t-8..t-3)) already have their
  standardised anomaly series computed (`climate_anomaly_covariates.csv`)
  but have never been fit to the climate-forced Stan model. Per
  instruction, this audit does not launch new Stan fits -- this is flagged
  as a **required, not-yet-performed** sensitivity check before the
  climate coefficients are trusted at a lag spec other than "primary".
- **No multi-state pooling exists.** beta_T/beta_P are CE-specific,
  single-state estimates; the existing assessment's own "Phase 4" (pooled
  BA/RJ/MT/CE climate coefficients) was explicitly deferred and remains
  undone.

## 5. Classification and decision

**Classification: B -- sufficient for uncertainty-driven sensitivity
scenarios, not for provisional posterior-based scenario exploration treated
as a settled climate-transmission relationship.**

Rationale:
- Computation is clean (0 divergences both models; Rhat at, not past, the
  gate) -- this audit does not raise a convergence concern.
- The climate signal is real but modest and only partially decisive:
  precipitation is suggestive (96% posterior mass positive) but not
  significant by a conventional 95% threshold; temperature shows no
  detectable effect at this lag spec. Both are the expected, intended
  outcome of a deliberately conservative, regularising prior on a
  single-state canary, not a disqualifying finding.
- Variance decomposition confirms climate explains a small (~1%, 95% CI up
  to ~3%) share of historical log R0(t) variation relative to seasonality
  (80%) and year effects (19%) -- consistent with "real but second-order."
- The newly surfaced beta_T/year_effect[2022] correlation (0.80) is a
  genuine, mechanistically explicable identifiability caveat: part of the
  temperature coefficient's posterior location is entangled with one
  externally-seeded recurrence year, not purely with climate variation
  itself. This specifically argues against treating beta_T/beta_P as an
  already-validated, portable climate-transmission relationship (class A).
- At the same time, nothing here shows the climate component is
  *non-identifiable* or *uninterpretable* (class C): the posterior is
  prior-updating, bounded, orthogonal to the pre-existing harmonic/year
  structure in aggregate, and the resulting reconstruction is stable.

**Implications for future climate projections (carried into Step 2):**
Step 2's structural pilot must **not** present beta_T/beta_P-driven future
R0(t) trajectories as validated climate-change projections. It should (and
does, per the task's own design) run climate-neutral (C0), resampled
historical-variability (C1), and a forcing-identity control (C2) as
parallel structural arms, explicitly labelled as such, rather than
presenting a single "the climate effect" trajectory. Before any actual
SSP/CMIP6 integration (Step 3), the lag-spec sensitivity fit and
multi-state pooling identified above should be completed, and the
beta_T/year_effect[2022] entanglement should be re-examined once a
longer or multi-state record is available to separate them.
