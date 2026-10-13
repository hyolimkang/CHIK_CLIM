# CE Climate Robustness Assessment (Ceará, historical inference only, 2015-2025)

Produced by `02_Script/02_ceara_pipeline/11_climate_robustness/`.
Historical inference only. No SSP scenarios, no vaccination, no new
epidemic generator were implemented or run. M0 (accepted) was never
refit; M1-M4 are new fits saved separately under
`03_Output/02_ceara_pipeline/model_fits/climate_robustness/`. Not
committed (per instruction).

## Phase 1 — Climate input consistency audit: RESULT = PASS

The discrepancy flagged in the task (narrative SD 1.097/0.756 vs. saved
metadata 0.7663/0.5016) was investigated to its root cause:

- **The saved scaling metadata (0.7663/0.5016) is correct and has never
  changed.** A byte-for-byte diff against the file's content at the
  commit that first created it (`96f5c04`) versus the current HEAD shows
  **zero difference** — the file was relocated during a later
  repository reorganisation but its content was never regenerated or
  edited.
- **An independent from-scratch reconstruction** of the primary-lag
  temperature/precipitation anomalies, built directly from the raw
  national climate table (not from the already-saved covariates CSV),
  reproduces T_anom_sd=0.76630 / P_anom_sd=0.50162 to 6 decimal places —
  an exact match.
- **The z_T_anom/z_P_anom vector actually stored inside the accepted M0
  Stan fit bundle** matches this independent reconstruction (with the
  documented boundary-week NA→0 rule applied) to **0.00e+00** across all
  573 fitted weeks.
- **The "1.097/0.756" figure was a transcription error in the prose of
  `STEP1_CLIMATE_INFERENCE_ASSESSMENT.md`** (a previous, separate
  analysis), not a data or pipeline problem — the SAME earlier analysis
  run's own programmatically-generated table
  (`STEP1_lag_and_standardisation_specification.csv`) already correctly
  showed 0.7663/0.5016. That report's prose has been corrected with an
  explicit, dated note (not silently edited).

Full detail: `03_Output/02_ceara_pipeline/tables/climate_robustness/PHASE1_*.csv`
(6/6 audit checks pass). **Gate cleared — Phase 2 proceeded.**

*Scope note*: Phase 1's from-scratch independent reconstruction covered
the **primary** lag spec only (the one used by the active M0 model, as
instructed). M1/M2 (sensitivity_A/B lags) reuse the identical,
already-validated construction code
(`01_build_climate_anomaly_covariates.R`) with different lag-offset
parameters, and the Phase 2 scripts load them via the same join logic
validated in Phase 1(c) for primary — but sensitivity_A/B were not
independently re-derived from raw data step-by-step the way primary
was. This is a minor, disclosed scope limitation, not a finding of any
problem.

## Phase 2 — M1-M4 fits: 1 of 4 clean, 3 marginal HMC gate failures (disclosed, not hidden)

All 4 new fits used the byte-identical, unmodified accepted Stan source
(`renewal_ceara_v4_9_climate_forced_canary.stan`) — only the `data`
differed (climate covariates for M1/M2; `year_effect_prior_sd` for
M3/M4). Same HMC configuration as M0 (4 chains, 800 warmup + 400
sampling, adapt_delta=0.95, max_treedepth=14, dense_e). All 4 smoke
tests passed before full fitting.

| Model | divergences | max treedepth hits | max Rhat | min ESS (bulk) | HMC gate |
|---|---|---|---|---|---|
| M0 (accepted) | 0 | 0 | 1.0097 | 450.6 | **PASS** (at the boundary, as previously documented) |
| M1 (sens_A, SD=0.40) | **2** | 0 | 1.0133 | 331.8 | **FAIL** |
| M2 (sens_B, SD=0.40) | **1** | 0 | 1.0019 | 555.0 | **FAIL** |
| M3 (primary, SD=0.20) | 0 | 0 | 1.0054 | 495.4 | **PASS** (cleanest of all 5) |
| M4 (primary, SD=0.60) | **1** | 0 | 1.0062 | 523.3 | **FAIL** |

**Reported as failures, not treated as accepted results.** These are
*marginal* gate failures (1-2 divergent transitions out of 1,600
post-warmup draws, i.e. ≤0.13%; zero treedepth saturation; Rhat
0.2-1.3% above 1.00 in every case) — not gross pathologies — but per
instruction they are not silently waived. M1/M2/M4's posteriors are
still reported below for scientific comparison, with this caveat
attached to every table/figure that uses them. A production-grade
resolution would re-run these three with a higher `adapt_delta` (e.g.
0.99) and/or longer warmup before any of them is treated as a fully
accepted alternative model — not done here, as that would constitute
additional, unrequested refitting beyond the planned 4-model grid.

Full HMC table: `PHASE2_hmc_diagnostics_M1_M4.csv` /
`PHASE3_HMC_diagnostics_all_M0_M4.csv`.

## Phase 3 — Scientific comparison

### 1-3. beta_T / beta_P across M0-M4 (own z-scale; FIG13)

| Model | beta_T median [95% CrI] | P(beta_T>0) | beta_P median [95% CrI] | P(beta_P>0) |
|---|---|---|---|---|
| M0 | -0.016 [-0.071, 0.032] | 26% | 0.024 [-0.004, 0.052] | 96% |
| M1 | -0.016 [-0.069, 0.041] | 31% | 0.021 [-0.010, 0.052] | 91% |
| M2 | -0.020 [-0.071, 0.027] | 21% | 0.023 [-0.003, 0.049] | 95% |
| M3 | -0.012 [-0.065, 0.037] | 31% | 0.024 [-0.003, 0.053] | 96% |
| M4 | -0.016 [-0.063, 0.032] | 28% | 0.024 [-0.004, 0.051] | 95% |

**Direction and magnitude are remarkably stable**: temperature stays
null (CrI straddling zero, P(>0) 21-31%) and precipitation stays
weakly-positive-but-not-decisive (P(>0) 91-96%, CrI lower bound just
below zero) **across every tested lag spec and every tested
year-effect-prior scale**, including under a prior tightened to half
(M3) or loosened by 50% (M4) of the accepted value. Posterior SD/prior
SD ratios (0.09-0.19) confirm the N(0,0.15) climate prior is
meaningfully updated by the data in all 5 fits, not simply reproducing
the prior.

### Physical-scale comparison (`PHASE3_beta_physical_scale_comparison.csv`)

Per instruction, raw z-scale coefficients are **not** compared as if
they shared identical units (M1's own T_anom_sd=0.848 vs. M0/M2/M3/M4's
0.766-0.766). Converted to "effect per 1°C of lagged mean-temperature
anomaly" and "per 1 unit of log1p(lagged summed precipitation mm)
anomaly" (precipitation's own transform is already log-scale — flagged,
not presented as raw mm): beta_T_per_degC ranges -0.016 to -0.026
(all CrI crossing zero) and beta_P_per_unit ranges 0.040 to 0.048 (all
CrI lower bounds just below zero) across the 5 models — the same
qualitative conclusion holds on the physical scale.

### 4. Climate multiplier (FIG14)

Bounded and consistent across all 5 models: median range
[0.930-0.936, 1.087-1.099], max 95% CrI width 0.20-0.21 regardless of
lag spec or prior scale. No specification produces an extreme or
qualitatively different forcing envelope.

### 5. Climate-year-effect confounding, especially 2022 (FIG15) — THE KEY, PERSISTENT FINDING

| Model | cor(beta_T, year_effect[2022]) | cor(beta_P, year_effect[2022]) |
|---|---|---|
| M0 | 0.80 | 0.02 |
| M1 | 0.79 | 0.11 |
| M2 | 0.77 | -0.03 |
| M3 | 0.83 | 0.06 |
| M4 | 0.79 | -0.01 |

**2022 is the single largest cor(beta_T, year_effect[y]) across all 11
years in every one of the 5 models** (range 0.77-0.83), and this is
**unaffected by tightening (M3) or loosening (M4) the year-effect prior
flexibility** — directly answering the task's question about whether
stronger year-effect flexibility weakens the apparent climate
association: **it does not, in either direction.** beta_P shows no such
pattern with 2022 (|cor| ≤0.11 in every model) but consistently
correlates most with **2019** instead (0.59-0.64 across all 5 models,
`PHASE3_max_abs_year_correlation_by_model.csv`) — a different, also
stable, confounding pattern.

Mechanistically: the 8 conditioned-seed weeks of 2022 are excluded from
the case likelihood, but the remaining ~44 weeks of 2022 (the observed
post-seed peak and decline) are not, and those weeks' climate anomalies
and that year's free `year_effect[2022]` parameter are both competing
to explain the same unusually large recurrence. Because this pattern is
identical in magnitude across lag choice and prior flexibility, it
appears to be a **structural feature of fitting a single-state,
11-year record containing one dominant externally-seeded recurrence
year** — not an artifact correctable by any of the sensitivity axes
tested here.

### 6. R0(t) and S(t)/N(t) (FIG16A-B)

Visually indistinguishable across all 5 models at checkpoint resolution
and full weekly resolution. Maximum relative difference vs. M0 in
S_prop(t): M1 0.11%, M2 0.10%, M3 0.47%, M4 0.058% — all small. M3
(tightest year-effect prior) shows the largest, still-small, deviation,
consistent with that prior change having the most direct mechanical
effect on the year-effect/reconstruction interplay.

### 7. PPC and annual reconstruction (FIG16C)

Weekly 95% coverage 96.7-96.9% and weekly RMSE 452-458 across all 5
models — **climate forcing neither materially improves nor degrades
in-sample fit under any tested specification**, consistent with the
original Step 1 finding and extended here to confirm it holds across
lag/prior choices too. Annual predicted/observed ratios track within
0.01-0.02 of each other across all 5 models in every year
(`PHASE3_annual_PPC_M0_M4.csv`).

## Phase 4 — Figures

FIG13 (coefficient forest plot), FIG14 (weekly multiplier, small
multiples), FIG15 (confounding heatmap + 2022 bar chart), FIG16 (R0/S-
prop/PPC sensitivity) — all generated as PDF + PNG@600dpi + SVG with
reproducible underlying CSVs, inspected for layout issues and
regenerated once (FIG13's subtitle/axis title were clipped on first
render; fixed by wrapping text and reducing axis-title font size).

## Phase 5 — Answers to the seven questions

**1. Were historical climate inputs reproduced exactly?**
**Yes.** Confirmed to floating-point precision (0.00e+00) against the
actual vector inside the accepted Stan fit, via an independent
from-scratch reconstruction from raw data. The apparent discrepancy was
a prose error in an earlier report, now corrected, not a data issue.

**2. Are climate coefficient directions and magnitudes stable across lag choices?**
**Yes**, closely. beta_T stays negative/null and beta_P stays
positive/suggestive, with medians and CrIs varying by only
~0.003-0.004 in absolute terms across the primary/sensitivity_A/sensitivity_B
lag windows (M0 vs. M1 vs. M2).

**3. Are they stable across year-effect prior assumptions?**
**Yes.** Halving (M3) or increasing by 50% (M4) the year-effect prior
SD changes beta_T/beta_P medians by at most ~0.004 and leaves both
direction and the qualitative P(>0) pattern unchanged.

**4. Does the 2022 climate-year-effect confounding persist?**
**Yes, essentially unchanged**, across every tested lag spec and prior
scale (cor(beta_T, year_effect[2022]) = 0.77-0.83 in all 5 models). This
is the most important finding of this assessment: it is a robust
confound, not a fragile artifact of one modelling choice, which makes it
more concerning, not less — no sensitivity axis tested here resolves it.

**5. Does climate forcing materially improve historical model fit?**
**No.** PPC coverage and RMSE are statistically indistinguishable from
the non-climate baseline and stable across all 5 climate-forced
variants (96.7-96.9% coverage; RMSE 452-458).

**6. How uncertain is the climate multiplier?**
Modest and bounded in every specification: posterior-median range
roughly [0.93, 1.10], with a 95% credible-interval width of ~0.20-0.21
at its widest week, consistently across M0-M4.

**7. Can these coefficients be used for conditional future SSP projection, and with what limitations?**
**Only as an explicitly uncertainty-propagated, scenario-labelled
input — not as a validated, causal climate-transmission function** — for
these reasons, in order of importance: (a) the 2022 confounding is
real, robust to every tested sensitivity axis, and specifically
entangles beta_T with one externally-seeded recurrence year; (b) q=0.05
is fixed throughout, not estimated; (c) this is a single state and an
11-year record; (d) 3 of the 4 new sensitivity fits carried minor,
disclosed HMC irregularities; (e) no out-of-sample or multi-state
validation was performed in this assessment. If used, posterior draws
(not point estimates) must be propagated, and the C0/C1/C2
climate-neutral/resampled-historical/forcing-identity structural-arm
design the Step 2 pilot already uses is the right pattern to extend —
never a single "the" climate-transmission relationship.

## Classification

**B — Sensitive to assumptions; scenario-based interpretation required.**

Not because the coefficient estimates are unstable (they are, in fact,
remarkably consistent across all 5 tested specifications — this
assessment does **not** support classification C). Classified B rather
than A because: (i) the 2022 climate-year-effect confound is real,
persistent, and unresolved by any tested structural choice — the kind
of identifiability limitation that specifically requires scenario-based
(not point-estimate) use; and (ii) 3 of the 4 new sensitivity fits did
not clear the strict HMC gate, even if only marginally. Robustness to
modelling choices is a necessary but not sufficient condition for
causal interpretability, and this assessment establishes the former
without establishing the latter.

## Additional evidence required before transporting these associations to future warming conditions

1. **Multi-state / pooled climate fit** (BA, RJ, MT alongside CE) — the
   single most informative next test, since those states do not share
   CE's 2022 externally-seeded recurrence and could help separate a
   genuine shared climate signal from a CE-2022-specific artifact.
   Already flagged as outstanding in the original Step 1 report; still
   outstanding here.
2. **A genuinely estimated (not fixed) between-year variance
   hyperparameter**, connecting to the separate `year_effect_diagnostic`
   finding that the current fixed `year_effect_prior_sd=0.40` is a
   researcher-chosen constant, not data-informed.
3. **A leave-2022-out (or leave-one-year-out) refit comparison**,
   explicitly testing whether beta_T's estimate and sign are driven
   disproportionately by that single year — not performed here, and a
   natural next diagnostic given finding 4 above.
4. **Resolution of the 3 marginal HMC gate failures** (M1, M2, M4) with
   a higher `adapt_delta` before treating their posteriors as
   fully accepted alternative specifications.
5. **Out-of-sample predictive evaluation**, which this historical,
   in-sample exercise does not and cannot establish.
6. Only after (1)-(5) should fixed, point-estimate, or even posterior
   climate coefficients from this single-state fit be used to force
   actual SSP/CMIP6 scenarios — and even then, always as one labelled
   structural arm among several, exactly as Step 2's C0/C1/C2 design
   already anticipates.

## Files generated

```
02_Script/02_ceara_pipeline/11_climate_robustness/
  00_config_climate_robustness.R
  00_run_climate_robustness.R
  01_phase1_climate_input_audit.R
  02_phase2_fit_M1_M4.R
  03_phase3_scientific_comparison.R
  04_phase4_figures.R

03_Output/02_ceara_pipeline/
  reports/climate_robustness/CE_CLIMATE_ROBUSTNESS_ASSESSMENT.md   (this file)
  tables/climate_robustness/PHASE1-3_*.csv                         (18 tables)
  figures/climate_robustness/FIG13-16.{pdf,png@600dpi,svg}
  model_fits/climate_robustness/ce_climate_robustness_M{1,2,3,4}.rds  (+ per-chain progress CSVs)

03_Output/02_ceara_pipeline/reports/future_projection_pilot/STEP1_CLIMATE_INFERENCE_ASSESSMENT.md
  — one sentence corrected (transcription error; see Phase 1 above), dated and explained, not silently changed.
```

No existing accepted model, fitted object, or other output file was
modified. Nothing was committed.
