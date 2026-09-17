# v4.8 Seeded-Recurrence Assessment

Scope: tests whether conditioning on an externally seeded 2022 CHIKV
reintroduction (motivated by independent genomic evidence of a newly
introduced ECSA lineage) resolves the systematic 2022 underprediction
found in v4.7's continuous fixed-q renewal model. No other structural
change was made relative to v4.7.

## 1. Exact v4.6 -> v4.7 model changes (verified by diff)

`git diff`-level comparison of `25_v4_6_fixed_q_sweep/scripts/01_fit_v4_6.R`
vs `26_v4_7_fixed_q_full_period/scripts/01_fit_v4_7.R` confirms the Stan
model file is **byte-identical and reused**
(`renewal_ceara_v4_6_fixed_q_sweep.stan` — N/Y are data, not hardcoded).
The only substantive change is on the R side:

- `date_end`: `2019-12-29` (v4.6) -> `2025-12-21` (v4.7), extending `Y`
  from 5 to 11 annual-effect years and `N` from 313 to 573 weeks.
- Folder/function/tag renames only (`v4_6_*` -> `v4_7_*`).
- Everything else identical: q grid (0.05-0.30), Juazeiro serology
  treatment (beta-binomial, `kappa_sero=50` fixed), sampler config
  (`dense_e`, `adapt_delta=0.95`, `max_treedepth=14`, 4 chains,
  400 warmup / 200 sampling), priors, GI kernel, init strategy.

## 2. Exact v4.7 -> v4.8 changes

New Stan file `renewal_ceara_v4_8_seeded_recurrence.stan`, derived from
`renewal_ceara_v4_6_fixed_q_sweep.stan` (i.e. the same model v4.7 runs)
with exactly these additions — no other line changed:

- **data**: `is_seed[N]`, `X_seed[N]`, `N_fit`, `fit_index[N_fit]`.
- **transformed parameters**: inside the weekly loop, `X[t]` is now
  `X_seed[t]` (with a `reject()` guard if it would exceed `S[t]`) when
  `is_seed[t]==1`, otherwise unchanged (`S[t] * (-expm1(-force_of_infection[t]))`).
  `force_of_infection[t]`, `R0_t[t]`, `R_eff_t[t]`, and the S/U update
  recursion are computed exactly as before, for every week, using
  whichever `X[t]` is available.
- **model**: the case likelihood becomes
  `C[fit_index] ~ neg_binomial_2(expected_reported_cases[fit_index], phi_obs)`
  instead of the unindexed `C ~ neg_binomial_2(...)` — this is the only
  likelihood change, and it only removes the 8 seed weeks.
- **generated quantities**: unchanged in form (still computed for all N
  weeks, for plotting continuity).

New R script `27_v4_8_seeded_recurrence/scripts/01_fit_v4_8.R`, derived
from `01_fit_v4_7.R`, adds `make_v4_8_data()` (builds `is_seed`/`X_seed`/
`fit_index` from a seed window specified as **R constants**,
`SEED_START_DATE`/`SEED_LENGTH`) and passes `WARMUP`/`ITER_SAMPLING`
through environment variables (used once, for a purely technical rerun of
q=0.20 at 2x iterations — see Section 7). Same q grid, same date range,
same sampler configuration, same priors as v4.7.

No file belonging to v1-v4.7 was modified.

## 3. Why q remains fixed

The mode-validity audit (v4.1-v4.5) already established that q is not
identifiable from weekly case counts plus one serosurvey point, regardless
of the observation-likelihood structure (weekly NB, 4-week block NB,
yearly-burden + Dirichlet-multinomial shape all showed persistent or
newly-collapsed-but-unresolved bimodality). v4.6/v4.7 deliberately reverted
to q as fixed data and instead read q's effect off a sensitivity sweep.
v4.8 tests a specific structural hypothesis (does re-seeding fix the 2022
underprediction) and must not simultaneously re-open the already-settled
q-identifiability question, per the explicit instruction not to make any
other substantive model extension at this stage. Fixing q also lets the
2022 comparison be attributed unambiguously to the seeding mechanism
rather than to a shifting q estimate.

## 4. Why waning immunity is not introduced

Waning immunity, a second serotype, immune escape, and time-varying q
are all flexible nuisance mechanisms that could independently "explain
away" the 2022 underprediction without testing any specific hypothesis.
The task specifies ONE structural test at a time: reintroduction-seeding
first, because there is independent genomic evidence for it (a newly
introduced ECSA lineage), whereas waning immunity, escape, or a second
serotype are not motivated by anything found in this dataset yet. Adding
any of them now would confound the seeding test's result. If v4.8 is
CASE 2 (fails), the pre-specified next candidate is spatial heterogeneity,
not waning immunity (Section 12).

## 5. How the 2022 seed is implemented

Seed window: **2022-01-02 to 2022-02-20 (8 consecutive weeks)**, chosen
and hard-coded only in R (`SEED_START_DATE`, `SEED_LENGTH` in
`01_fit_v4_8.R`), not in the `.stan` file. Rationale:

- Length = `G` = 8, the generation-interval kernel length used throughout
  v4.0-v4.7 — the minimum window needed to prime the renewal convolution
  (`sum_g w[g]*X[t-g]`) before handing control back to the ordinary
  renewal equation.
- Start date = the first week of calendar year 2022, matching the
  calendar-year definition of "the 2022 recurrence" used everywhere else
  in this analysis, and immediately preceding the period of clearly
  exponential growth in the raw data (44 -> 75 -> 112 -> 153 -> 309 ->
  388 -> 561 -> 647 cases/week across the 8 seed weeks).
- For each week `t` in the seed window, `X_seed[t] = C[t] / q`, using the
  exact same `expected_reported_cases = q * X` mapping v4.7 already
  assumes everywhere — no reporting delay or other observation mapping
  exists in v4.7 to respect beyond this.
- All 44 remaining weeks of 2022 (2022-02-27 onward — the entire
  observed peak of ~3,500/week in mid-May and the full decline back to
  baseline by December) are generated purely by the renewal equation and
  are the actual object of the test.
- No new estimated parameter was introduced: `X_seed` is a deterministic
  function of data (`C[t]/q`), not a fitted initial-state parameter,
  matching the "Option A" design (condition + exclude) rather than
  "Option B" (small anchored initial-state parameters), because v4.7's
  observation model is a direct `q * X` mapping with no delay
  distribution to invert, making Option A the cleaner fit to the existing
  architecture.
- Susceptibility never resets: `S[t]`/`U[t]` are updated by the same
  recursion for every week of the whole 2015-2025 series, consuming
  whichever `X[t]` is available (renewal-derived or seed-conditioned).
  `S` at the start of the seed window is exactly "susceptibility
  remaining after all earlier infections" — nothing before the seed
  window was touched.

## 6. How double use of seed observations is prevented

The 8 seed-week case counts (`C[t]` for `t` in the seed window) are used
**once**, to construct `X_seed[t] = C[t]/q` outside the likelihood, in R.
Inside the Stan model, `fit_index` excludes exactly these 8 indices from
`C[fit_index] ~ neg_binomial_2(...)`, so the same observations never also
contribute a likelihood term. `N_fit = N - 8 = 565` for every q (verified
in the fit logs: `N_fit=565 of N=573`). Generated quantities still compute
`C_pred`/`log_lik_cases` for the seed weeks for plotting continuity only —
these are clearly marked (grey band, Figure panel A) and excluded from
every PPC/error metric reported in this document and in
`v4_8_seeded_recurrence_comparison.csv` (all "2022" metrics there are
computed on the post-seed weeks only, `pred_2022_postseed_*`).

## 7. HMC results for every q

| q | divergences | max treedepth hits | max Rhat | min bulk ESS | BFMI (per chain) | HMC gate |
|---|---|---|---|---|---|---|
| 0.05 | 0 | 0 | 1.0074 | 502 | all >0.3 | PASS |
| 0.10 | 0 | 0 | 1.0060 | 451 | all >0.3 | PASS |
| 0.15 | 0 | 0 | 1.0024 | 601 | all >0.3 | PASS |
| 0.20\* | 0 | 0 | 1.0010 | 1441 | all >0.3 | PASS |
| 0.25 | 0 | 0 | 1.0091 | 711 | all >0.3 | PASS |
| 0.30 | 0 | 0 | 0.9976 | 720 | all >0.3 | PASS |

\*q=0.20's first run (same settings as the other five: 400 warmup / 200
sampling) had 0 divergences, 0 treedepth hits, ESS=314, BFMI healthy, but
max Rhat=1.0201 — a borderline miss on the strict Rhat threshold alone,
not a real geometric pathology. Per Stan's own diagnostic guidance
("running the chains for more iterations may help"), q=0.20 was rerun
with a purely technical change (800 warmup / 400 sampling, identical
sampler configuration otherwise) and passed cleanly (values shown above).
All other q values used the standard 400/200 settings unchanged from v4.7.

**All six q values pass the strict HMC gate.** Seeding did not introduce
any real sampling pathology — the geometry with q fixed is fundamentally
well-behaved, exactly as in v4.6/v4.7.

## 8. 2017 PPC comparison

| q | observed 2017 | v4.8 predicted median | 95% CrI |
|---|---|---|---|
| 0.05 | 107,421 | 103,433 | [73,989 , 142,749] |
| 0.10 | 107,421 | 106,379 | [75,089 , 148,833] |
| 0.15 | 107,421 | 107,478 | [77,180 , 155,897] |
| 0.20 | 107,421 | 112,217 | [80,001 , 157,768] |
| 0.25 | 107,421 | 115,650 | [82,039 , 164,100] |
| 0.30 | 107,421 | 111,155 | [80,338 , 153,897] |

Effectively unchanged from v4.7 (2017 is untouched by the 2022 seed
mechanism, as expected since the seed window only affects weeks in
2022). The 2017 fit is not degraded by adding the seed mechanism.

## 9. 2022 PPC comparison (post-seed weeks only, excludes the 8 conditioned weeks)

| q | observed (post-seed) | v4.7 (continuous, full 2022) | v4.8 (seeded, post-seed) |
|---|---|---|---|
| 0.05 | 48,697 | median 26,434 [19,369 , 36,155] | median 25,040 [19,137 , 32,908] |
| 0.10 | 48,697 | median 24,460 [17,861 , 32,681] | median 23,687 [17,419 , 31,315] |
| 0.15 | 48,697 | median 24,043 [17,418 , 32,742] | median 23,580 [17,613 , 32,082] |
| 0.20 | 48,697 | median 23,880 [17,246 , 32,575] | median 24,304 [18,347 , 32,685] |
| 0.25 | 48,697 | median 24,308 [17,839 , 33,035] | median 24,474 [18,276 , 33,479] |
| 0.30 | 48,697 | median 24,306 [17,790 , 32,757] | median 25,980 [18,897 , 35,056] |

**Seeding changes the 2022 prediction by essentially nothing** (differences
of a few hundred to ~1,500 cases out of a ~25,000 median, well within
Monte Carlo noise) — see Figure panel F, where the v4.7 and v4.8 curves
are visually superimposed and both sit roughly 2x below the observed
level (dashed line) at every q.

2022 weekly peak: observed peak = 3,521 cases (week of 2022-05-15). v4.8's
predicted-median peak across the q grid is 1,176-1,295 cases, occurring
2022-04-24 to 2022-05-08 — **peak timing is approximately correct
(within 1-3 weeks) but peak magnitude is underpredicted by roughly
2.7-3.0x**, consistently across the whole q grid.

Susceptible fraction immediately before the seed window (2021-12-26):
71.8% (q=0.05) rising to 94.7% (q=0.30) — i.e. even at the q value with
the most depletion, nearly three-quarters of the population is still
susceptible going into 2022. Susceptible exhaustion is not the limiting
factor at any q.

## 10. Does re-seeding resolve the structural underprediction?

**No.** Across the entire pre-registered q grid, conditioning on an
external 2022 reintroduction leaves the systematic underprediction
essentially unchanged: the seeded model's post-seed 2022 total remains
at roughly half the observed count (median ratio observed/predicted
~1.9-2.1x, matching the descriptive ascertainment-multiplier diagnostic
in Section 11), and the weekly peak is underpredicted by ~2.7-3x, with
correct timing but wrong amplitude. Since susceptibility is confirmed not
exhausted (71.8-94.7% remaining) at the point of reintroduction, the
model's inability to reproduce the wave's amplitude cannot be attributed
to a lack of susceptible individuals, nor is it fixed by simply
supplying the missing seed. Given a real 2022 introduction, the existing
(spatially homogeneous, state-level) transmission model and its existing
R0(t)/susceptibility bookkeeping are not sufficient on their own to
reproduce the observed amplitude of the recurrence.

## 11. Decision: PASS / FAIL

**FAIL (Case 2).** The reintroduction/seeding hypothesis, tested exactly
as specified (fixed q, no other structural change, minimum-length data-
conditioned seed, seed weeks excluded from the likelihood), does not
resolve the 2022 underprediction. Per the pre-registered decision rule,
we do NOT respond by fitting a 2022-specific ascertainment parameter or
introducing waning immunity. The descriptive-only diagnostic below
quantifies how large a pure-ascertainment explanation would have to be,
strictly as a benchmark:

| q | v4.7 predicted 2022 (median) | observed 2022 (full year) | required ascertainment multiplier |
|---|---|---|---|
| 0.05 | 26,434 | 50,986 | 1.93x |
| 0.10 | 24,460 | 50,986 | 2.08x |
| 0.15 | 24,043 | 50,986 | 2.12x |
| 0.20 | 23,880 | 50,986 | 2.14x |
| 0.25 | 24,308 | 50,986 | 2.10x |
| 0.30 | 24,306 | 50,986 | 2.10x |

(`v4_7_ascertainment_multiplier_diagnostic.csv`. This multiplier is
reported only as a descriptive benchmark for how large a surveillance-
practice change would need to be if invoked as an alternative
explanation — it was never used as a likelihood weight, prior, or fitted
quantity anywhere in v4.7 or v4.8.)

## 12. Next structural candidate (design note only, NOT implemented)

Per the pre-registered fallback: the 2022 Ceará epidemic is reported to
have concentrated in mesoregions that were less affected by the 2016-2017
epidemic. A state-level, spatially homogeneous susceptibility pool
(v4.0-v4.8's shared assumption) cannot represent this: it forces every
region's 2017 depletion to apply uniformly, which is exactly the
assumption that would suppress the 2022 amplitude if the true
2022-susceptible mesoregions were largely spared in 2017.

Proposed (not implemented) next model: replicate the v4.7/v4.8 renewal
process independently across Ceará's ~7 health mesoregions `r`:

- `I[r,t]`: latent weekly infections in mesoregion `r`.
- `S[r,t]`, `U[r,t]`: mesoregion-specific susceptible/immune bookkeeping,
  each with its own 2015-2019 depletion history (driven by each region's
  own case time series, same GI kernel and R0(t) form, possibly a shared
  seasonal term with mesoregion-specific annual effects).
- `C[r,t]`: observed weekly cases per mesoregion, same fixed-q
  observation mapping and same q sensitivity scenario applied uniformly
  across regions (do not re-open q per region without separate
  justification).
- The 2022 recurrence would then be tested against each mesoregion's
  own remaining susceptibility, which is the direct test of the
  "concentrated in previously-less-affected mesoregions" hypothesis.
- Juazeiro serology would anchor whichever mesoregion contains Juazeiro
  do Norte specifically, rather than the state-level aggregate.

This spatial model is a materially larger undertaking (mesoregion-level
demographic/case panels, 7x the latent state, likely per-mesoregion HMC
diagnostics) and should only be scoped in detail if explicitly requested.

## Files produced

- `02_Script/stan/renewal_ceara_v4_8_seeded_recurrence.stan`
- `02_Script/40_renewal_model/27_v4_8_seeded_recurrence/scripts/01_fit_v4_8.R`
- `02_Script/40_renewal_model/27_v4_8_seeded_recurrence/scripts/02_plot_v4_8_sweep.R`
- `02_Script/40_renewal_model/27_v4_8_seeded_recurrence/outputs/q{0.05..0.30}/renewal_ceara_v4_8_fit_q*.rds`
- `02_Script/40_renewal_model/27_v4_8_seeded_recurrence/v4_8_seeded_recurrence_comparison.csv`
- `02_Script/40_renewal_model/27_v4_8_seeded_recurrence/v4_7_ascertainment_multiplier_diagnostic.csv`
- `03_Output/figures/renewal_v4_8_seeded_recurrence/v4_8_seeded_recurrence_summary.png`
- This document.

All v1-v4.7 model versions, scripts, and outputs are untouched.
