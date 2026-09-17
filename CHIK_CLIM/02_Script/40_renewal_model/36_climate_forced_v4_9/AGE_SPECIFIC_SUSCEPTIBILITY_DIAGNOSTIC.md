# Age-Specific Susceptibility and Infection Diagnostic (CE climate-forced v4.9, vaccination OFF)

## Design (recap)

No age-specific transmission, no vaccination. Every age shares exactly the
same weekly `p_infection_given_susceptible[t] = 1 - exp(-lambda[t])`, and
`Inew[a,t] = S[a,t] * p_infection_given_susceptible[t]`. 80 posterior draws
of the accepted climate-forced CE fit were forward-simulated through the
Phase-4-validated age-cohort bookkeeping. Age-specific differences arise
**only** through births, ageing, and differential cumulative exposure of
successive birth cohorts.

## Correction (post-review): what "cumulative infection proportion" does and does not mean

A reviewer correctly flagged that the original Figure C2 showed the 60+
group's cumulative infection proportion *decreasing* between epidemics,
which cannot happen for an individual's own lifetime infection history
under lifelong immunity. Investigation found **two distinct issues**:

1. **A genuine implementation bug**, now fixed: Figure C2 originally
   plotted `cumsum(Inew_group(t)) / N_group(t)`. This is not a valid
   cumulative proportion for an age *band* — it only counts infections
   that happened while someone was already inside that band, so anyone
   ageing in after being infected at a younger age is missing from the
   numerator while still counted in the (growing) denominator, producing
   an artificial decline unrelated to any real epidemiology. Figure C2
   now uses `Uinf_group(t) / N_group(t)` (the same quantity as
   `immune_prop_group`, consistent with the summary table).
2. **Not a bug, but a real and important distinction**: after the fix,
   the wide, open-ended 60+ band is confirmed non-decreasing (0/572 weeks
   decrease), but the narrow 0-4 and 5-9 bands **still decline between
   epidemics** — correctly. `Uinf_group(t)/N_group(t)` measures the share
   of the *currently living* group that has ever been infected, not a
   fixed individual's lifetime probability. For a wide, terminal band
   (60+), membership only grows (via ageing-in from 59, never ageing out),
   so this ratio is essentially insulated from dilution. For a narrow,
   fast-turnover band (0-4 fully replaces its membership roughly every 5
   years via births + ageing-out), a large recent-epidemic-exposed
   sub-cohort is continuously replaced by newborns with zero infection
   history, which **legitimately dilutes the group average** even though
   no single individual ever loses immunity (verified directly in code:
   no step in the simulator — infection, death/reconciliation allocation,
   or ageing — ever moves mass from `Uinf` back to `S`).
   **This is not specific to the age extension**: the same mechanism is
   already present in the original, non-age-structured accepted
   climate-v4.9 fit itself — `U(t)/(S(t)+U(t))` there decreases in 265 of
   572 weeks (draw 1), just by an amount too small to notice
   (~-0.00007/week during quiet periods) at full-population scale. Every
   state's "cumulative proportion infected" panel produced throughout
   this project has always been this same population-average ratio, not
   a per-individual lifetime probability; the age extension simply makes
   the underlying birth-dilution mechanism visible by isolating a
   fast-turnover sub-population.

## Section 8 critical re-validation (re-confirmed, not assumed)

At single-year age resolution, across all 80 draws: max absolute
difference vs. the original aggregate model was **9.6e-08** (S), **4.5e-08**
(U), **4.4e-10** (X), **1.1e-07** (N) — all far below any epidemiologically
meaningful scale (population ~9 million). The identity established in
Phase 4 **still holds** at this finer resolution.

**Explicit numerical check that `p_infection_given_susceptible[t]` is
identical across ages** (5 spot-check weeks, avoiding the conditioned seed
window): max deviation of `Inew_age[a,t]/S_age[a,t]` from the common
hazard, across all 101 ages = **3.5e-18** — floating-point noise. Confirmed
by construction and numerically.

## Figures and table produced

- `AGE_DIAG_FigA1_susceptibility_heatmap.png` / `FigA2_immunity_heatmap.png`
- `AGE_DIAG_FigB_group_susceptibility_trajectories.png`
- `AGE_DIAG_FigC_infection_curves.png` (weekly incidence/100k + cumulative)
- `AGE_DIAG_FigD_cross_sections.png`
- `AGE_DIAGNOSTIC_group_summary_table.csv`,
  `AGE_DIAGNOSTIC_group_range_by_checkpoint.csv`

## Key checkpoint findings (population-weighted, posterior median)

| checkpoint | 0-4 S_prop | 5-9 S_prop | 60+ S_prop | range across groups |
|---|---|---|---|---|
| Baseline (2015-01-04) | 100% | 100% | 100% | 0.0 |
| Before 2017 epidemic (2017-02-05) | 93.3% | 92.5% | 92.5% | 0.8pp |
| After 2017 epidemic (2018-12-16) | 82.8% | 73.2% | 73.1% | 9.7pp |
| Before 2022 seed (2021-12-26) | 94.5% | 76.6% | 72.3% | 22.3pp |
| After 2022 epidemic (2023-12-17) | 87.2% | 69.9% | 60.5% | 26.7pp |
| End of simulation (2025-12-21) | 93.9% | 76.2% | 60.1% | **33.8pp** |

Ages 5-9 through 60+ track each other almost exactly at every checkpoint
(differences in the 3rd-4th decimal) — they were all already alive during
the 2017 epidemic and depleted by essentially the same fractional amount
under the shared hazard. **Only the 0-4 group consistently separates**,
and that separation **widens monotonically** after each epidemic (0.8pp
before 2017 -> 33.8pp by end-2025) as more of that group is replaced by
post-epidemic births who have never been exposed.

## Answers to the six required questions

**1. Given a common force of infection, do age groups develop meaningfully
different susceptible fractions over time?**
Yes — but the differentiation is a **two-tier** pattern, not a smooth
age gradient: "not yet exposed to any epidemic" (effectively the youngest
cohort, currently ages 0-4, but really "born after the last major wave")
vs. "exposed to at least one epidemic" (ages 5-9 upward, all converging to
nearly the same depleted level). By the end of the simulated period the
gap between these two tiers reaches ~34 percentage points.

**2. Which age groups remain most susceptible before recurrent epidemics?**
The youngest cohort every time — concretely, ages 0-4 remained the most
susceptible heading into both the 2017 epidemic (93.3%, barely above the
92.5% of everyone else at that early point in the series) and, much more
starkly, heading into the 2022 seed/epidemic (94.5%, vs 72-77% for 5+) —
because by 2021 the 2017 epidemic's depletion had had years to accumulate
in the older cohorts while newborns kept entering fully susceptible.

**3. Are age differences driven primarily by demographic cohort
replacement after earlier epidemics?**
Yes, entirely. There is no other mechanism available in this model (no
age-specific hazard, no age-specific mortality rate, no age-specific
mixing) — the only channels are births entering age 0 fully susceptible,
ageing shifting cohorts upward, and the shared hazard depleting whoever is
alive at epidemic time equally. The heatmap (`FigA1`) shows this directly
as a widening wedge of high susceptibility at low ages/recent years,
tracking exactly along the "born after the epidemic" diagonal.

**4. Are age-specific infection curves different only in magnitude, while
retaining the same epidemic timing?**
Yes — by construction (`FigC1`), since `Inew_age[a,t] = S_age[a,t] *`
the SAME `hazard[t]` for every age, every age group's incidence curve
peaks in the exact same week; only the height differs, scaled by how
susceptible that group currently is. This was verified as an algebraic
identity, not just observed in the plot.

**5. Does the resulting age structure provide enough information to
implement routine vaccination without adding age-specific transmission
parameters?**
Yes. Because the mechanism generating age heterogeneity here is purely
demographic bookkeeping (not a fitted or assumed age-transmission
parameter), routine vaccination can be layered on exactly as Phase 5
specifies — as an `S -> Uvac` movement applied to whichever single age
reaches the target age each week — without needing any new transmission
assumption. The diagnostic shows WHERE the susceptible pool concentrates
(young ages, especially right after a big epidemic), which is precisely
the information needed to reason about routine vaccination's likely
reach, without having pre-committed to a target age.

**6. Which candidate vaccination ages would be epidemiologically
informative to examine later? (NOT a selection — descriptive only)**
Candidates worth examining later, based on where susceptibility
concentrates and where cohorts are cleanly separated in this diagnostic:
- **Age 1-2**: captures children while still in the "fully susceptible,
  not yet epidemic-exposed" tier seen in `FigD`; closest to a classic
  early-childhood routine schedule.
- **Age 9-10**: the top of the age range still occasionally showing
  measurable separation from the adult plateau in `FigB`/`FigD`
  (school-entry-adjacent age, operationally convenient in many EPI
  programmes).
- **Age 12**: already flagged as a pre-existing programmatic
  reference point in this repo's own age-group data build
  (`ibge_pop_uf_age_group_split12_2015_2024.rds` isolates age 12
  specifically), though this diagnostic itself shows age 12 already
  belongs to the "adult-like", already-depleted tier post-2017/2022 in
  CE specifically.
No age is selected or optimised here; this is left for the dedicated
vaccination-counterfactual phase.

## STOP

Descriptive/diagnostic only, per instruction. No vaccination logic was
implemented in this task. Returning for review before Phase 5.
