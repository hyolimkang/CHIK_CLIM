# Closed-Loop Replay Validation and Vaccine-Engine Unit Tests

## Clarification on scope

The age-cohort simulator (`07_age_cohort_simulator.R`) was **already
closed-loop** at the time of the earlier Phase-4 identity test — its
`infectiousness[t]` term is built from the function's own recursively
generated `X_total` history; the function signature never accepts the
fitted `X(t)` as an input at all. This document reports the **stricter**
validation battery now requested (infectiousness, force_of_infection,
R_eff comparisons; era-specific breakdowns; explicit divergence search),
which re-confirms and strengthens that earlier result, plus the full
vaccine-engine unit-test suite (Phase 3) and the routine age-12 mechanism
with cohort QA (Phases 4-5).

## 1. Closed-loop replay errors (Phases 1-2)

60 posterior draws, CE climate-forced canary fit, q=0.05. Vaccination
OFF. Every quantity compared against the ORIGINAL fitted transformed
parameters.

| quantity | max abs diff | max rel diff |
|---|---|---|
| S | 9.1e-08 | 1.6e-14 |
| U | 4.9e-08 | 1.4e-14 |
| X | 5.1e-10 | 1.1e-13 |
| infectiousness | 5.3e-10 | 1.0e-13 |
| force_of_infection (lambda) | 1.3e-16 | 9.0e-14 |
| R_eff | 2.9e-15 | 2.3e-15 |

**Era-specific breakdown** (max relative diff within era, all quantities): pre-2017
(≤4.1e-15), 2017 epidemic window (≤2.2e-14), inter-epidemic 2019-2021
(≤8.0e-14), 2022 seed window (≤1.6e-16 for X — seed weeks are exactly
conditioned, not hazard-driven), 2022 epidemic post-seed (≤6.3e-14),
post-2022 (≤1.1e-13). No era shows elevated error relative to any other.

Cumulative infections: max \|relative diff\| across draws = 5.97e-15.
Final susceptible fraction: max \|abs diff\| across draws = 9.88e-15.

**First-divergence-week search** (threshold 1e-6 relative error): **no
week, in any of the 60 draws, exceeds the threshold** for S, U, or X. All
observed error is floating-point accumulation noise, not systematic
drift — confirmed by the error staying at the same ~1e-13 to 1e-16 order
of magnitude across every era rather than growing over the 573-week
horizon.

**PHASE 1-2 VERDICT: PASS.**

## 2. Vaccine=0 (coverage=0) identity test — TEST A

`max|diff S_total| = 0.000e+00`, `max|diff X_total| = 0.000e+00` —
**exactly identical** to the no-vaccine simulation (not merely close).

**TEST A VERDICT: PASS.**

## 3. VE_infection=0 identity test — TEST B

Coverage = 0.9, VE_infection = 0. `max|diff S_total| = 0`,
`max|diff X_total| = 0` (epidemiologically identical to no-vaccine), while
`doses_administered` totals 1,251,252 over the historical window (doses
ARE counted) and `effective_protected` totals exactly 0 (no one is
actually moved to `Uvac`, since a 0%-effective vaccine protects no one).

**TEST B VERDICT: PASS** (identical epidemiology, doses still counted, as required).

## 4. Age-12 one-time pulse — direct/indirect effect diagnostic — TEST C

100% coverage, VE=0.8, applied for **one week only** (2018-06-10, a quiet
inter-epidemic week).

- Entrants to age 12 that week: 2,515.49 total, 1,835.60 susceptible.
- Doses = 2,515.49 (= all entrants, coverage=1.0); effective_protected =
  1,468.48 (= 0.8 x 1,835.60, exact VE arithmetic).
- **Same-week** `lambda(t)` difference vs. no-vaccine: **0** (force of
  infection is built from PAST infections only — vaccination cannot
  retroactively change the current week's transmission risk; causality
  confirmed).
- Age-12 susceptibles drop by exactly 1,468.48, starting the pulse week
  itself (direct effect).
- `X_total` first declines vs. baseline **exactly 1 week after** the
  pulse (2018-06-17) — the minimum possible causal lag, since this
  week's reduced age-12 infections feed into next week's infectiousness.
- Cumulative infections averted over the following ~3 years: age 12
  (direct) = 23.3; ages 0-11 (**never vaccinated**, indirect) = 196.8;
  ages 65+ (**never vaccinated**, indirect) = 124.9. The indirect effect
  in unvaccinated groups is far larger than the direct effect, because
  those groups are much larger populations sharing the same reduced
  transmission pressure — the expected herd-effect signature. See
  `PHASE3_testC_pulse_direct_indirect.png`: averted infections in ages
  0-11 and 65+ start at zero and grow with each subsequent transmission
  season (2019, 2020, 2021), confirming the causal chain `S_12 down ->
  X_total down -> infectiousness down -> lambda down -> infections down
  in every age group`.

**TEST C VERDICT: PASS** (causal sequence and indirect-effect mechanism confirmed).

## 5. Routine age-12 cohort eligibility table (Phases 4-5)

**Scope note**: the real 2026-2050 climate/transmission scenario has
**not** been chosen (per explicit instruction). To validate the routine-
vaccination *mechanism* over a realistic multi-decade horizon, this
section uses a clearly-flagged, **non-forecast placeholder scaffold**:
the fitted 2025 seasonal R0(t) pattern and 2025 demographic flow levels
are repeated forward unchanged through 2050. This is mechanism validation
only, not a projection.

Coverage=0.9, VE=0.8, target age=12, active from 2026.

| birth cohort | year turning 12 | eligible | doses | effective |
|---|---|---|---|---|
| 2014 | 2026 | 117,247 | 105,522 | 54,845 |
| 2020 | 2032 | 108,611 | 97,750 | 65,232 |
| 2030 | 2042 | 96,967 | 87,270 | 71,739 |
| 2038 | 2050 | 87,391 | 78,652 | 65,697 |

(full 25-row table: `PHASE5_cohort_qa_summary.csv`)

**Cohort QA**: 25 distinct birth cohorts (2014-2038), **25 rows** in the
summary table — **no birth cohort appears more than once**
(`n() > 1` check returns zero matches). Structurally, the vaccination
logic only ever reads/writes the single `target_age` column
(age 12 = column 13) of the age11->12 **inflow slice**, never the
standing age-12 stock and never age 13 — so by construction a cohort
cannot be vaccinated twice (once at 12, again the following year at 13).
Each cohort's ~117 vaccination weeks span almost exactly one calendar
year (52-53 weeks), consistent with the continuous/fractional ageing
flow spreading each birth cohort's arrival at age 12 across the whole
year it turns 12, rather than a single discrete date.

**PHASE 4-5 VERDICT: PASS** (mechanism validated; no duplicate vaccination).

## Outputs

Tables: `PHASE1_2_closed_loop_overall_summary.csv`,
`PHASE1_2_closed_loop_era_summary.csv`,
`PHASE1_2_closed_loop_cumulative_final_summary.csv`,
`PHASE1_2_divergence_weeks_above_threshold.csv`,
`PHASE3_unit_tests_A_B_summary.csv`, `PHASE3_testC_pulse_diagnostic.csv`,
`PHASE5_cohort_qa_weekly.csv`, `PHASE5_cohort_qa_summary.csv`.

Figures: `PHASE1_2_closed_loop_replay_validation.png`,
`PHASE3_testC_pulse_direct_indirect.png`, `PHASE5_cohort_qa_summary.png`,
`PHASE4_5_mechanism_scaffold_impact.png`.

## STOP

Per instruction: no real 2026-2050 climate/transmission scenario has been
chosen or invented. The vaccine-capable forward simulator and routine
age-12 mechanism are built and validated; the mechanism-test scaffold
used for cohort QA is explicitly NOT a forecast. Awaiting review and the
future-scenario decision before any policy-relevant 2026-2050 analysis.
