# Bahia External Replication of the Frozen Ceará v4.9 Model

Scientific question: does the frozen Ceará long-horizon renewal-model
structure reproduce major chikungunya epidemic trajectories in Bahia
across fixed-ascertainment scenarios, and how sensitive is the implied
susceptibility history to ascertainment? **This is a replication
exercise, not a re-tuning exercise** — no Bahia-specific model change was
made based on how well Bahia fit.

## 1. The frozen Ceará model used

`renewal_ceara_v4_9_hierarchical_seasonality.stan`
(`28_v4_9_hierarchical_seasonality/`), confirmed via
`V5_0_DYNAMIC_R_ASSESSMENT.md`'s freeze decision as the final Ceará
structure (v5.0a/v5.0b were NOT used — both failed on computational/
geometric grounds before any scientific verdict on their dynamic-R idea).

| component | specification |
|---|---|
| R scripts | `17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R` (generation weights, `compute_hmc_gate`, `build_state_weekly`), `22_v4_3.../scripts/02_fit_v4_3.R` (`make_v4_3_data`, `find_sero_week_index`, `load_td14_scalar_inits`), `27_v4_8.../scripts/01_fit_v4_8.R` (seed-window machinery, unused for Bahia), `28_v4_9.../scripts/01_fit_v4_9.R` (`make_v4_9_data`, second-harmonic construction) |
| generation interval | 8-week discrete kernel, `diff(pgamma(0:8, shape=4, rate=2))`, normalised |
| R0(t) | `alpha_R + year_effect[year] + A_year[year]*(beta_sin1*sin1+beta_cos1*cos1) + beta_sin2*sin2 + beta_cos2*cos2`; `A_year = exp(sigma_season_year * z_season_year)` (non-centred hierarchical seasonal amplitude) |
| demographic S/U | exact reconciliation identity `N_end - N_start - births + deaths`; susceptible/immune split updated every week from `X[t]` |
| observation model | `C[fit_index] ~ neg_binomial_2(q*X + 1e-9, phi_obs)`; `phi_obs ~ gamma(2,0.1)` |
| imports | fixed `imports_per_week = 1` |
| seed/reintroduction | conditional mechanism (`is_seed`/`X_seed`/`fit_index`) exists in the shared Stan interface but is a NO-OP for Bahia (see Section 5) |
| q grid | 0.05, 0.10, 0.15, 0.20, 0.25, 0.30 — identical to the final Ceará sensitivity analysis, **no q privileged for Bahia** |
| priors | `alpha_R~N(log1.2,0.5)`, `z_year~N(0,1)`, `beta_sin1/cos1~N(0,0.25)`, `beta_sin2/cos2~N(0,0.10)`, `sigma_season_year~N+(0,0.20)`, `phi_obs~Gamma(2,0.1)` — all unchanged from Ceará |

## 2. Confirmation: no transmission priors or structure altered

The ONLY code change relative to frozen v4.9 is
`renewal_ceara_v4_9_hierarchical_seasonality_optional_serology.stan`,
which adds a single `use_serology` data flag gating the one serology
likelihood line (`if (use_serology) { sero_pos ~ beta_binomial(...); }`).
Verified: every other line — the renewal recursion, GI kernel, all
priors, the harmonic/A_year structure, the demographic accounting, the
seed mechanism — is byte-identical to
`renewal_ceara_v4_9_hierarchical_seasonality.stan`. Bahia runs with
`use_serology=0`; Ceará's own results are unaffected (never re-run with
this file). HMC configuration (dense_e, adapt_delta=0.95,
max_treedepth=14, 4 chains) is identical; only sampling iterations were
increased for q=0.05 and q=0.30 to resolve marginal Rhat (Section 6),
per the explicit exception allowing additional iterations.

## 3. Bahia data audit

`01_bahia_data_audit.R` (`03_Output/tables/renewal_bahia_v4_9_replication/bahia_data_audit.csv`).
Bahia demography was newly built
(`00_data_prep/09b_fetch_bahia_sinasc_weekly_births.R`,
`10b_build_bahia_weekly_demography.R`, reusing the SAME cached national
SINASC zip files as Ceará's own fetch, just filtered to UF prefix "29"
instead of "23", and the identical generic
`interpolate_population()`/`distribute_annual_flow_to_weekly()`
functions). All checks passed:

- 573 weeks, 2015-01-04 to 2025-12-21 — identical window and week count to Ceará v4.9.
- Dates continuous, no duplicates, case/demography joins matched exactly.
- Demographic accounting identity (`N_end-N_start-births+deaths-reconciliation`) max error = 0.
- Population 14.44M–14.88M (much larger than Ceará's ~9M).
- Total reported cases 2015-2025: 92,498 (mean 161/week, max 1,499 in the week of 2024-03-17).

## 4. Bahia epidemic-wave and seed/reintroduction definitions

**Waves**: taken as-is from the existing, already-computed national
wave-census pipeline
(`03_Output/tables/national_wave_analysis/brazil_chik_wave_analysis_master.csv`,
state="BA") — NOT re-derived. 9 major-epidemic waves, `BA_wave_02`
through `BA_wave_10` (onset 2016-01-24 to 2024-01-21), each with
objectively-defined onset/peak/end weeks from that pipeline's peak-trough
segmentation rule (trough_fraction=0.30, trough_persistence=3 weeks,
major_epidemic = ≥100 total cases AND ≥5/100k cumulative incidence).

**Seed/reintroduction**: **none used.** Applying the same objective
long-gap standard implied by Ceará's own 2022 event — the ONLY Ceará
transition ever requiring re-seeding had `weeks_since_previous_wave=399`
(~7.7 years), while every other Ceará transition (up to 56 weeks) was
handled by the continuous renewal process alone — to Bahia's 8 inter-wave
gaps (14, 21, 28, 35, 35, 42, 49, 49 weeks, all from
`weeks_since_previous_wave` in the same census file): **none approach
the ~7-year scale that triggered Ceará's exception.** Per the explicit
instruction not to force a reintroduction the data doesn't call for,
`is_seed` is all-zero and every week is included in the likelihood
(`fit_index = 1:N`). This is consistent with Bahia's own case series,
which shows near-annual recurring waves rather than one dominant epidemic
followed by a multi-year silence.

| wave_id | onset | gap since previous wave | seed used? |
|---|---|---|---|
| BA_wave_02 | 2016-01-24 | — (first) | no |
| BA_wave_03 | 2016-10-16 | 49 wk | no |
| BA_wave_04 | 2017-11-12 | 21 wk | no |
| BA_wave_05 | 2019-02-03 | 49 wk | no |
| BA_wave_06 | 2020-01-19 | 35 wk | no |
| BA_wave_07 | 2021-01-03 | 14 wk | no |
| BA_wave_08 | 2021-12-05 | 42 wk | no |
| BA_wave_09 | 2022-11-06 | 35 wk | no |
| BA_wave_10 | 2024-01-21 | 28 wk | no |

## 5. Computational pilots and HMC diagnostics by q

Pilots (150/400 warmup/sampling) at q=0.05 and q=0.10 both showed clean
geometry (0 divergences, 0 treedepth hits, healthy BFMI) with only a
marginal Rhat (~1.03) — diagnosed as short-pilot Monte Carlo noise (all 4
chains agreed closely on `alpha_R`/`phi_obs`; the only NA-Rhat parameters
found in a full scan were structurally-constant week-1 quantities
[`X[1]`,`U[1]`,`S_prop[1]`,`immune_prop[1..2]`] and unused serology
placeholders [`p_state_sero_at_anchor`,`log_lik_serology`,`sero_pred`] —
not evidence of genuine pathology), matching the same pattern resolved by
more iterations throughout Ceará's own v4.9 development. Proceeded
directly to the full q sweep (800/400) per protocol.

| q | divergences | max treedepth hits | max Rhat | min bulk ESS | min tail ESS | BFMI (min) | HMC gate | elapsed |
|---|---|---|---|---|---|---|---|---|
| 0.05 | 0 | 0 | 1.006 | 508 | 834 | 0.794 | **PASS** | 262 s |
| 0.10 | 0 | 0 | 1.001 | 668 | 863 | 0.798 | **PASS** | 278 s |
| 0.15 | 0 | 0 | 1.006 | 630 | 888 | 0.780 | **PASS** | 286 s |
| 0.20 | 0 | 0 | 1.008 | 580 | 835 | 0.829 | **PASS** | 340 s |
| 0.25 | 0 | 0 | 1.003 | 611 | 626 | 0.855 | **PASS** | 351 s |
| 0.30\* | 0 | 0 | 1.001 | 1377 | 1584 | 0.815 | **PASS** | 426 s |

\*q=0.30's first 800/400 run showed the identical marginal-Rhat-only
pattern (max Rhat 1.010, 0 divergences, 0 treedepth, ESS/BFMI healthy);
rerun at 1600/800 resolved cleanly. **All 6 q values pass the strict HMC
gate.**

## 6. Weekly case PPC by q

FIGURE 1 (`bahia_figure1_weekly_ppc_by_q.png`) and the individual
six-panel trajectory figures
(`bahia_v4_9_trajectories_q{0.05..0.30}_full.png`). Posterior-median
trajectories track the repeated near-annual wave pattern closely across
the whole 2015-2025 series at every q; predictive intervals widen
somewhat but the median shape is consistent. Unlike Ceará (one dominant
2017 epidemic dwarfing everything else), Bahia shows comparably-sized
waves nearly every year, especially from 2020 onward — a qualitatively
different epidemic regime.

## 7. Wave-total PPC by q

FIGURE 4 (`bahia_figure4_wave_total_ppc_by_q.png`);
full table in `bahia_wave_level_ppc_summary.csv`. Total-case ratios
(predicted median / observed) are **q-invariant and mostly well
calibrated**: 7 of 9 waves fall in [0.93, 1.09] across the ENTIRE q grid.
One systematic exception:

- **BA_wave_10 (2024-2025, the most recent/ongoing wave)**: ratio ≈0.85-0.86
  at every single q — a genuine, q-independent ~15% underprediction of
  total burden, not resolved anywhere in the grid.

## 8. Wave-peak / timing / width PPC by q

FIGURE 5 (`bahia_figure5_wave_peak_ppc_by_q.png`). Two systematic,
q-invariant patterns, in the OPPOSITE direction from Ceará's known 2017
peak-smoothing limitation:

- **Peak magnitude is generally OVERpredicted**: ratios range from ~1.10
  (BA_wave_02) up to ~1.95 (BA_wave_09, i.e. nearly double the observed
  peak) for 7 of 9 waves, at every q. Only BA_wave_10 is underpredicted
  (~0.75).
- **Peak timing is systematically delayed**: the predicted-median peak
  week is later than the observed peak week for essentially every wave
  (e.g. BA_wave_03: observed 2017-01-15 vs predicted ~2017-03-12/26, a
  ~8-9 week lag; BA_wave_09: observed 2023-02-05 vs predicted
  ~2023-03-26/04-09, a ~7-9 week lag), consistently across the whole q
  grid.
- **Epidemic width** (weeks containing the central 80% of cases) is
  generally wider in the model than observed (e.g. BA_wave_05: 24 obs vs
  35-37 predicted; BA_wave_10: 47 obs vs 54-55 predicted), consistent
  with a smoother, more drawn-out predicted epidemic shape even where the
  peak itself is overpredicted.

## 9. Susceptibility trajectories by q

FIGURE 2 (`bahia_figure2_susceptibility_by_q.png`); per-wave pre/post
values in `bahia_wave_level_ppc_summary.csv`. Because Bahia's population
(~14.5-14.9M) is much larger relative to its case counts than Ceará's,
depletion is modest even at the lowest q: susceptible fraction drops from
100% to only ~88% (q=0.05) or ~94-97% (q≥0.15) by the end of 2025.
Pre/post-wave susceptibility declines monotonically and coherently across
all 9 waves at every q — no discontinuities or implausible jumps. As
required, **no single q is identified as preferred**; results are
reported as `S_Bahia(t | q)` across the full grid.

## 10. Comparison of independent Ceará and Bahia posteriors

`ceara_vs_bahia_posterior_comparison.png`/`.csv`
(Ceará available at q=0.05-0.20 only, per the earlier decision to treat
q=0.05 as the Ceará reference and not complete its 0.25/0.30 grid — no
pooling performed, each state fit fully independently):

- **Heterogeneous**: `alpha_R` is much higher for Ceará at low q (0.294
  vs Bahia's 0.072 at q=0.05), converging toward similar small values by
  q=0.20 — reflecting Ceará's need for a much higher baseline transmission
  intercept to produce one enormous epidemic versus Bahia's steadier
  recurring waves. `sigma_season_year`/SD(A_year) is roughly **2x larger
  for Ceará** (0.63-0.68 vs Bahia's 0.24-0.39) — Ceará's transmission
  intensity varies far more dramatically year-to-year (one huge year,
  several quiet ones) than Bahia's more uniform annual pattern.
  `phi_obs` is **markedly higher (less overdispersed) for Ceará**
  (~7.6-7.8 vs Bahia's ~2.8-3.0) — Bahia's weekly counts are noisier
  relative to their mean.
- **Similar**: `A_year` median across years is close to 1 in both states
  at every q (correctly centred, as designed); first-harmonic amplitude
  coefficients are the same order of magnitude, though see below.
- **Possibly weakly identified / differing seasonal phase**: `beta_sin1`
  is smaller and `beta_cos1` larger for Bahia (0.15-0.18 / 0.28-0.29) than
  Ceará (0.24-0.26 / 0.16-0.17) — the two states' epidemics peak in
  different parts of the annual cycle, a real, interpretable difference
  rather than a fitting artefact.

## 11. Transferability classification: **PARTIALLY TRANSFERABLE**

Applying the pre-registered criteria:

- HMC is completely stable across the full q grid (Section 5) — rules
  out NOT TRANSFERABLE.
- Long-term/annual burden is reproduced reasonably (7/9 waves within
  ±10% at every q) — satisfies the burden criterion for either
  TRANSFERABLE or PARTIALLY TRANSFERABLE.
- However, epidemic **peak magnitude and timing are systematically
  biased** in the same direction across nearly every wave and across the
  entire q grid (general overprediction of peak height, a consistent
  multi-week delay in peak timing, and wider-than-observed epidemic
  width) — this is exactly the "systematically smoothed or biased"
  epidemic-shape failure mode the PARTIALLY TRANSFERABLE category
  describes, not random per-wave noise that would still count as
  TRANSFERABLE.
- One wave (BA_wave_10, the most recent) also shows a genuine,
  q-independent total-burden shortfall (~15%), on top of the shape bias.

**Conclusion: the frozen Ceará v4.9 structure is PARTIALLY TRANSFERABLE
to Bahia** — it captures the overall recurring-wave regime, annual
burden, and coherent susceptibility depletion well, but systematically
mis-shapes individual epidemic peaks (opposite in direction from Ceará's
own known 2017 under-shaped-peak limitation: Bahia's peaks are generally
over-predicted and delayed rather than under-predicted).

## 12. Explicit list of model-data discrepancies

1. Peak magnitude overpredicted for 7/9 major waves (ratio up to ~1.95),
   at every q — most extreme for BA_wave_09 (2022-2023).
2. Peak timing delayed by ~1-9 weeks relative to observed, for nearly
   every wave, at every q.
3. Epidemic width generally wider (smoother rise/decline) than observed.
4. BA_wave_10 (2024-2025, the most recent/still-partially-censored wave)
   is the sole exception in the opposite direction: total burden ~15% low
   and peak ~25% low, at every q — plausibly a right-censoring/data-
   completeness artefact given this wave's `end_week` extends to the
   final available data week (2025-12-28), rather than a genuine model
   failure; flagged for attention but not adjudicated here.
5. `phi_obs` (weekly NB2 dispersion) is markedly different between states
   (Ceará ~7.7, Bahia ~3.0) — the same fixed prior (`Gamma(2,0.1)`)
   accommodates both, but the posteriors do not resemble each other,
   consistent with genuinely different residual/overdispersion structure
   rather than a shared nuisance scale.
6. Seasonal phase (`beta_sin1`/`beta_cos1` ratio) differs between states —
   plausibly reflects real differing climate-driven seasonality, not
   treated as a discrepancy requiring correction.

None of these were used to modify the Bahia model. Per the external-
replication stop rule, **no Bahia-specific v4.10 was created.**

## 13. Recommendation on proceeding to the next independent state

**Proceed to at least one or two further independent states** before
considering any hierarchical pooling or structural change. Bahia's
peak-shape bias (over-prediction + delay) is the *opposite* pattern from
Ceará's under-prediction of its single sharp 2017 peak — with only two
states, it is not yet possible to tell whether this reflects (a) a
structural harmonic-R0 limitation that manifests differently depending on
each state's specific wave frequency/amplitude, (b) genuine
state-specific epidemiological differences (e.g. Bahia's larger, more
urbanized, more heterogeneous population smoothing local peaks
differently), or (c) chance. A third state would materially sharpen this
picture. Only after several states have been independently assessed
should hierarchical pooling (explicitly NOT attempted here, including for
q — Section 5's per-state fixed-q sensitivity approach was deliberately
retained) or a further structural change be considered.

## Files produced

- `02_Script/00_data_prep/09b_fetch_bahia_sinasc_weekly_births.R`, `10b_build_bahia_weekly_demography.R`
- `01_Data/bahia_sinasc_births_daily.rds`, `bahia_sinasc_births_weekly.rds`, `bahia_weekly_demography_2015_2025.rds`
- `02_Script/stan/renewal_ceara_v4_9_hierarchical_seasonality_optional_serology.stan`
- `30_bahia_v4_9_replication/scripts/01-05` (data audit, fit, six-panel plot, wave PPC + figures 1-6, Ceará-Bahia comparison)
- `30_bahia_v4_9_replication/outputs/q{0.05..0.30}/renewal_bahia_v4_9_fit_q*.rds`
- `03_Output/tables/renewal_bahia_v4_9_replication/` (data audit, HMC diagnostics, wave-level PPC summary, Ceará-Bahia comparison, per-q trajectory summaries)
- `03_Output/figures/renewal_bahia_v4_9_replication/` (6-panel trajectories per q, figures 1-6, diagnostics dashboard, Ceará-Bahia comparison)
- This document.

No hierarchical model was implemented. The frozen v4.9 structure was not
changed based on Bahia's fit quality.
