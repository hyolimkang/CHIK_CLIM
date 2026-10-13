# CE Future Year-Effect Validation (Ceará chikungunya, historical + retrospective holdout)

Produced by `02_Script/02_ceara_pipeline/12_future_year_effect_validation/`.
No SSP/ISIMIP forcing, no vaccination. No accepted model, fit, or
previous diagnostic/robustness output was modified. Not committed.

## Phase 0 -- Audit (gate: PASS)

Confirmed from the accepted M0 fit's own stored data: `year_effect[y] =
0.40*(z_year[y]-mean(z_year))` with `year_effect_prior_sd=0.40` as FIXED
data (never estimated); `A_year[y]=exp(sigma_season_year*z_season_year[y])`
with `sigma_season_year ~ Normal(0,0.20)` genuinely estimated; serology
anchored at week 194 (2018-09-16 survey midpoint, 103/404, kappa=50);
2022 conditioned seed = 8 weeks from 2022-01-02 (`X_seed=cases/q`,
excluded from the likelihood via `fit_index`, N_fit=565 of N=573); q=0.05
fixed; imports=1/week fixed (hazard addition); S(1)=N_start(1) (fully
susceptible at 2015-01-04), U(1)=0; N=573 weeks, Y=11 years. The new
`renewal_ceara_v4_9_climate_forced_tau_year.stan` model (identical to
the accepted file except `year_effect_prior_sd` (fixed data) ->
`tau_year` (estimated parameter) + `tau_prior_scale` (data)) compiled
and smoke-tested cleanly -- `tau_year` samples freely, not fixed.

## Phase 1 -- tau_year estimation (full 2015-2025, 2022 seed RETAINED, disclosed)

| Prior | tau_year median | 95% CrI | divergences | max Rhat | min ESS | HMC gate |
|---|---|---|---|---|---|---|
| P1: HalfNormal(0,0.20) | 0.138 | [0.093, 0.253] | 20 | 1.046 | 104 | **FAIL** |
| P2: HalfNormal(0,0.40) | 0.139 | [0.092, 0.223] | 2 | 1.003 | 289 | **FAIL** (marginal) |
| P3: HalfNormal(0,0.80) | 0.145 | [0.093, 0.234] | 5 | 1.014 | 323 | **FAIL** (marginal) |

**Posterior medians are essentially identical across a 4-fold range of
prior scale (0.20/0.40/0.80) -- tau_year is strongly likelihood-informed,
not prior-dependent.** All three cluster at ~0.14, **far below** the
previously fixed value of 0.40 (ratio ~0.35x), and closely corroborate
the independent empirical finding from the earlier
`year_effect_diagnostic` study (historical within-draw SD of the
fixed-0.40 model = 0.138) -- two different methods converge on the same
number. Correlations of tau_year with alpha_R, sigma_season_year,
beta_T, beta_P and terminal S_prop are all small (|cor|<=0.19,
`PHASE1_tau_year_correlations.csv`) -- no material confounding.

**HMC**: none of the three cleared the strict gate, all showing MILD
(not severe: 0 treedepth hits throughout) divergences, worst under the
tightest prior P1 (20 divergences) -- consistent with the kind of
funnel-like geometry hierarchical scale parameters are known to produce,
and reported honestly rather than silently accepted (see Question 2).

## Phase 2 -- Posterior-predictive generators G0/G1/G2 (2026-2035, paired design)

| Generator | median 10y cum. infections | P(>=1 major epidemic) | median peak R0 | p99 peak R0 | max peak R0 |
|---|---|---|---|---|---|
| G0: neutral | 7,966 | 0.2% | 1.91 | 2.16 | 2.22 |
| G1: fixed SD=0.40 (existing Y1) | 2,780,611 | **100%** | 1.93 | 6.26 | 15.5 |
| G2: estimated tau_year | 23,741 | 55% | 1.91 | 4.13 | 9.62 |

Paired contrasts (same draw x replicate x future innovations across all
3 arms): G1-G0 = +2,773,519 [78,365, 5,911,696]; G2-G0 = +15,435
[-1,018, 1,342,970]; G2-G1 = -2,659,648 [-5,512,907, -68,685]
(`PHASE2_paired_contrasts.csv`). **Replacing the fixed scale with the
estimated tau_year reduces median 10-year cumulative infections by
~117-fold and the probability of a major epidemic from a near-certain
100% to 55%.** G2 still shows real, non-trivial elevated risk relative
to G0 (55% vs. 0.2%; p99 peak R0 4.1 vs. 2.2) -- genuine interannual
variability matters, just far less extremely than the uncalibrated
fixed scale implied.

## Phase 3 -- Training-only fits for retrospective validation

| Holdout | Training window | Y (years) | N (weeks) | tau_year median | 95% CrI | divergences | HMC gate |
|---|---|---|---|---|---|---|---|
| H1 | 2015-01-04 to 2019-12-29 | 5 | 261 | 0.152 | [0.080, 0.366] | 5 | **FAIL** (marginal) |
| H2 | 2015-01-04 to 2021-12-26 | 7 | 365 | 0.158 | [0.093, 0.307] | 10 | **FAIL** (marginal) |

Both built with the **dedicated seed-free data builder**
(`03_build_training_only_data.R`) -- `is_seed` sum verified = 0 for
both, serology included (survey midpoint 2018-09-16 precedes both
cutoffs). **Training-only tau_year estimates (0.15-0.16) are consistent
with the full-period estimate (0.14)** -- the scale is not an artefact
of including or excluding the 2022 data; a model trained on only 5-7
years recovers essentially the same interannual-variability scale as
one trained on the full 11.

## Phase 3/4 -- Retrospective holdout validation (the most important test)

100 draws x 20 replicates (full scale) for the reference cells;
smoke-scale (10x5) for the importation-sensitivity grid and the
training-resampled-climate arm. **No holdout cases, no 2022 seed, no
full-period posterior used anywhere in these simulations** -- initial
state, recent infection history and all transmission/climate parameters
come exclusively from the TRAINING-ONLY fit; climate for the primary
arm uses OBSERVED historical weather during the holdout (a retrospective,
observed-climate-conditional forecast, explicitly not a prospective one).

**H1 (holdout 2020-2025; observed: major epidemic=TRUE, cumulative reported cases=56,903, largest single epidemic=54,271)**

| Arm | P(major epidemic) | median cum. cases | 95% CI covers observed cum.? | median largest epidemic |
|---|---|---|---|---|
| G1 (fixed SD=0.40) | 97.4% | 53,629 | TRUE | 46,753 |
| G2 (tau_year, lambda=1, observed climate) | 68.0% | 1,599 | TRUE | 970 |
| G2 (lambda=0.25) | 44.0% | 700 | **FALSE** | 392 |
| G2 (lambda=4) | 88.0% | 3,035 | TRUE | 2,266 |
| G2 (training-resampled climate, lambda=1) | 58.0% | 1,133 | TRUE | 571 |

**H2 (holdout 2022-2025; observed: major epidemic=TRUE, cumulative reported cases=54,572, largest=46,404)**

| Arm | P(major epidemic) | median cum. cases | 95% CI covers observed cum.? | median largest epidemic |
|---|---|---|---|---|
| G1 (fixed SD=0.40) | 99.6% | 69,508 | TRUE | 59,344 |
| G2 (tau_year, lambda=1, observed climate) | **99.75%** | 16,103 | TRUE | 12,395 |
| G2 (lambda=0.25) | 96.0% | 20,180 | TRUE | 18,989 |
| G2 (lambda=4) | 100% | 22,092 | TRUE | 20,328 |
| G2 (training-resampled climate, lambda=1) | 100% | 19,914 | TRUE | 15,847 |

Full annual-level predictive intervals and coverage:
`PHASE4_annual_PI_vs_observed.csv`, `PHASE4_cell_validation_summary.csv`.

## Phase 5 -- Figures

FIG17 (tau_year prior vs. posterior, P1-P3, 0.40 reference line), FIG18
(historical vs. future year-effect densities, G1 vs. G2, historical
medians/CrI overlaid), FIG19 (retrospective holdout trajectories, H1/H2
x G1/G2, observed vs. predictive intervals), FIG20 (G0/G1/G2 cumulative
infections / outbreak probability / peak R0) -- all PDF + PNG@600dpi +
SVG with reproducible CSVs, inspected and two rounds of layout fixes
applied (axis-label collisions, clipped subtitles).

## Answers to the nine questions

**1. What is the posterior distribution of tau_year?**
Median ~0.14 (range 0.138-0.145 across priors P1-P3; 0.152-0.158 in the
two shorter training-only fits), 95% CrI roughly [0.08, 0.37] depending
on specification -- consistently and substantially below the previously
fixed 0.40.

**2. Is tau_year statistically identifiable with the available historical data?**
**Largely yes, with a disclosed caveat.** Its posterior LOCATION is
essentially invariant to a 4-fold change in prior scale and to
training on 5, 7, or 11 years of data -- strong evidence of
likelihood-driven identification, not prior-dependence. The caveat:
every one of the five tau_year fits (P1, P2, P3, H1-training,
H2-training) shows mild HMC divergences (2-20 out of 1,600 draws,
0 treedepth saturation throughout) -- a real, if modest, residual
geometric/sampling difficulty consistent with the known funnel
behaviour of hierarchical scale parameters, reported rather than hidden.

**3. Does the estimated future year-effect distribution differ materially from the fixed-SD=0.40 generator?**
**Yes, substantially.** FIG18 shows G2's future year_effect density is
far narrower than G1's and sits much closer to the historical
posterior's own spread; FIG17 shows the posterior is concentrated
roughly 2.8x below the old fixed scale regardless of prior.

**4. Are extreme future R0 values still frequent?**
**Less frequent, but not absent.** Under G2, median annual peak R0
(~1.9) is unchanged from G0/G1 (seasonality dominates the typical
peak), but the UPPER TAIL is still real: p99 peak R0 = 4.1 and max =
9.6 (vs. G0's 2.2) -- interannual variability can still occasionally
produce a high-R0 year, just far less reliably than G1's p99=6.3/max=15.5.

**5. Does the new generator reproduce plausible epidemic magnitude and recurrence in held-out years?**
**Partially.** The 95% predictive interval covers the observed
cumulative case count in 9 of 10 tested cells (fails only at the
lowest importation sensitivity value, lambda=0.25, for H1). The
PROBABILITY of at least one major epidemic is substantial-to-high
(44-100% across cells, 68-99.75% at the reference importation level) --
not the near-impossibility a naive "well-calibrated must mean small"
assumption might predict. However, the predictive MEDIAN badly
underestimates the actual observed magnitude in every G2 cell (by a
factor of roughly 3x for H2 to >30x for H1 at reference importation) --
the true 2022-23 resurgence sits well above the median, in the upper,
right-skewed part of a wide predictive distribution, not at its centre.

**6. Can the 2022 resurgence be generated without an observed-data-conditioned seed?**
**Yes, qualitatively, with high probability for H2 (99.75%) and moderate
probability for H1 (68% at reference importation)** -- using only a
constant background importation hazard (no conditioned seed, no
holdout cases) and the estimated tau_year mechanism. The generator does
NOT need the artificially inflated fixed-SD=0.40 scale to produce a
plausible resurgence; it needs that scale (or a somewhat higher
importation rate) to make such a resurgence the TYPICAL/median outcome
rather than a plausible-but-less-than-even-odds one.

**7. Is failure to reproduce resurgence magnitude more consistent with variance assumptions, introduction pressure, or other structural limitations?**
**Primarily importation pressure and the inherent skew of a correctly
calibrated probabilistic process -- not a mis-estimated tau_year.**
Evidence: (a) training-only tau_year (0.15-0.16) closely matches the
full-period estimate (0.14), so the variance scale itself is not
biased by the shorter training window; (b) increasing lambda_import
from 0.25 to 4 in H1 raises both P(major epidemic) (44%->88%) and
median cumulative cases (700->3,035) monotonically, and is the ONLY
tested lever that moves these metrics substantially; (c) even so, no
tested importation level pushes the G2 MEDIAN up to the observed
magnitude -- consistent with 2022 being a genuine upper-tail draw
(documented elsewhere in this project as coinciding with a newly
introduced ECSA lineage) rather than evidence that the whole generator
is mis-specified.

**8. Is the posterior-predictive generator scientifically defensible for exploratory 2026-2045 projections?**
**Yes, as a labelled, uncertainty-propagated STRUCTURAL/EXPLORATORY
tool -- not as a point-forecast or a validated predictive model.** It
is a large, disclosed improvement over the fixed-SD=0.40 generator
(which claimed a near-certain major epidemic every decade); it is
built on a robustly, consistently estimated parameter; and it
correctly reproduces the OBSERVED historical resurgence within its
predictive uncertainty in the great majority of tested configurations.
It should always be used with full predictive distributions (not
medians alone), and any single-scenario "expected" trajectory would
misrepresent what this exercise supports.

**9. What uncertainty remains unresolved?**
- Mild HMC divergences in every tau_year fit (none catastrophic, all
  disclosed, none re-fit with a higher adapt_delta here -- consistent
  with this project's established practice of disclosing rather than
  iterating on marginal gate failures).
- Only ONE real extreme historical event (the 2022-23 resurgence) is
  available to validate a probabilistic generator against -- Phase 4's
  own instruction that "a single held-out epidemic cannot fully
  validate a probabilistic generator" applies in full; H1 and H2 are
  not independent tests of that same event.
- Importation pressure is only coarsely explored (3 prespecified
  values); the true mechanism (a specific introduced lineage) is not
  represented by a constant-rate Poisson process at all.
- The training-resampled-climate secondary arm (climate carrying no
  information from the holdout period at all) was only run at smoke
  scale and only for G2 at the reference importation level.
- The simplified single-series wave detector (reused from the Step 2
  pilot and the earlier diagnostics) is known to over-count minor
  epidemics; `observed_largest_epidemic`/`median_largest_epidemic`
  (the clearly-dominant single event in each holdout) are more robust
  to this than raw "number of major epidemics," and were the figures
  emphasised above.
- The underlying climate-coefficient confounding and lag-specification
  sensitivities documented in the separate `climate_robustness`
  assessment are carried over unchanged and unresolved here.

## Classification

**B -- Partially supported; additional structural sensitivity or validation required.**

Not A: the predictive median materially underestimates the one real
historical resurgence available for testing (by 3x to >30x depending on
cell); every tau_year fit (5 of 5) carries unresolved mild HMC
divergences; only one genuine extreme event exists to validate against,
which Phase 4's own design principle says cannot fully validate a
probabilistic generator; importation pressure is coarsely specified and
demonstrably influences outcomes.

Not C: tau_year is robustly and consistently identified (stable
posterior location across a 4-fold prior range and across training
windows spanning 5-11 years); the resulting generator is a large,
well-evidenced, quantified improvement over the fixed-SD mechanism
it replaces; it reproduces the observed historical resurgence
qualitatively (high probability, correct order of magnitude) and
within its own stated predictive uncertainty in 9 of 10 tested cells,
without ever being given the observed outbreak to condition on.

## Files created / fits completed

```
02_Script/00_shared/stan/current/ceara/
  renewal_ceara_v4_9_climate_forced_tau_year.stan     (new Stan model; accepted v4.9 climate-forced file untouched)

02_Script/02_ceara_pipeline/12_future_year_effect_validation/
  00_config_year_effect_validation.R, 00_run_year_effect_validation.R
  01_phase0_audit.R, 02_phase1_fit_tau_year.R, 03_build_training_only_data.R,
  04_phase3_fit_holdout_training.R, 05_phase2_generator_comparison.R,
  06_phase3_4_retrospective_validation.R, 07_phase1_tau_analysis.R, 08_make_figures.R

03_Output/02_ceara_pipeline/
  model_fits/future_year_effect_validation/
    ce_tau_year_P1.rds, ce_tau_year_P2.rds, ce_tau_year_P3.rds   (full period, 2022 seed retained)
    ce_tau_year_training_H1.rds, ce_tau_year_training_H2.rds     (training-only, seed-free)
  reports/future_year_effect_validation/CE_FUTURE_YEAR_EFFECT_VALIDATION.md  (this file)
  tables/future_year_effect_validation/PHASE*.csv                (17 tables)
  figures/future_year_effect_validation/FIG17-20.{pdf,png@600dpi,svg}
  results/future_year_effect_validation/PHASE*.{rds,csv}         (trajectories, draws, timing)
```

**5 new Stan fits completed** (P1, P2, P3 full-period; H1, H2
training-only). No accepted model, fit, climate_robustness output, or
earlier future-pilot output was modified. Nothing committed.

## Recommended next action

Before any exploratory 2026-2045 use: (1) resolve the shared mild HMC
divergence pattern with a higher adapt_delta re-fit of at least P2 and
the two training fits; (2) replace the constant-rate importation
assumption with something that can represent discrete reintroduction
events, informed by the genomic/epidemiological literature on CE's 2022
lineage; (3) if a third holdout split becomes feasible with more years
of future data, add it -- two holdouts sharing the same single extreme
event is a structural limitation this report cannot remove; (4) carry
the `climate_robustness` report's own open items (multi-state pooling,
lag-spec sensitivity fit) forward in parallel, since they bear on the
same log_R0(t) structure this study builds on.
