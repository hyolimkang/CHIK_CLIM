# Historical Age-12 Routine Vaccination Counterfactual: Ceará, 2015–2021

**Status: mechanism analysis, historical window only. Not a policy effectiveness
estimate. Not extended beyond 2021.**

## 0. Scope and what this document is not

This analysis asks a narrow, mechanistic question: *under the transmission
conditions Ceará actually experienced between 2015 and 2021, as reconstructed
by the accepted climate-forced v4.9 renewal model, how much infection burden
would a routine age-12 vaccination programme have prevented, and through what
pathway?*

It is **not**:
- a projection of future (2022+) burden or of any 2026–2050 scenario,
- an estimate of the impact of any specific, realistic coverage/efficacy
  combination — Scenario 1 uses coverage = 100% and infection-blocking
  VE = 100% with lifelong protection, deliberately chosen as an upper-bound
  **mechanism stress test**, not a plausible vaccine profile,
- a cost-effectiveness, DALY, or economic evaluation,
- a re-estimation of any parameter in the frozen v4.9 or climate-forced
  models. All transmission parameters, climate effects, and initial
  susceptibility/immunity conditions are taken as given from the accepted fit
  `36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds`,
  unmodified.

Everything below answers the 10 questions the task specified, using the
already-completed QA (`diagnostics/`), tables (`tables/TABLE1-3_*.csv`), and
figures (`figures/main/`, `figures/supplementary/`) in this directory.

## 1. QA pass/fail (Section 7 of the task spec)

All 10 required QA checks **PASS**. Full results in
`diagnostics/QA_overall_summary.csv` and the per-check CSVs in the same
folder.

| # | Check | Result |
|---|-------|--------|
| 1 | Arm A reproduces the accepted closed-loop no-vaccine replay | **PASS** (max abs diff S = 7.45e-08, X = 7.13e-10 vs the accepted fit, across sampled draws) |
| 2/3 | Population and compartment conservation (S+Uinf+Uvac=N) | **PASS** (max residual ≈ 7e-08, all arms — see Fig. S2) |
| 4 | No duplicate vaccination (each birth cohort vaccinated in exactly 1 calendar year) | **PASS** — confirmed both numerically and visually (Fig. 3D) |
| 5 | Arm C non-target groups match Arm A | **PASS for the two groups with no legitimate pathway to differ** (0–11, 65+: relative diff ≈ 0.03%, Fig. S3). **Not scored, but reported**: 13–17 and 18–64 diverge increasingly from Arm A in Arm C — this is expected direct-effect carryover through ageing, not a defect (see Section 5, Section 6 unexpected behaviour) |
| 6 | Arm B non-target reductions require propagation time (not instantaneous) | **PASS** — first divergence at week 9–14 depending on group, never week 1 |
| 7 | R0(t) identical across arms | **PASS by construction** (identical input vector passed to all three arm simulations) |
| 8 | R_eff(t) allowed to differ | Informational — differs as expected (see Section 5) |
| 9 | lambda(t) may diverge in Arm B after propagation | Informational — first divergence at week 8 |
| 10 | Numerical decomposition identity (total = direct-only + indirect) | **PASS exactly** — max relative error = 0.00e+00 across all 200 draws |

No QA failure occurred, so no halt was triggered.

## 2. Total effect (Table 1, Table 3)

Over 2015–2021, comparing Scenario 0 (no vaccine) against Scenario 1
(routine age-12 vaccination, 100% coverage / 100% infection-blocking VE,
lifelong protection — the **mechanism stress test**, not a policy value):

| Quantity | No vaccine | Vaccine (stress test) |
|---|---|---|
| Cumulative infections | 2,463,215 [2,271,193 – 2,745,044] | 1,927,811 [1,704,723 – 2,265,601] |
| Cumulative reported cases (q = 0.05, fixed) | 123,161 | 96,391 |
| Peak weekly infections | 112,013 | 72,156 |
| Final (end-2021) susceptible proportion | 0.7414 | 0.7140 |
| Doses administered | — | 913,032 |

**Total infections averted: 538,773 [473,374 – 567,351] (posterior median [95%
CrI]), a 21.8% relative reduction** in cumulative infections over the
2015–2021 window, at a cost of ~913,000 doses (i.e. ≈0.59 infections averted
per dose in this stress-test scenario — again, this ratio is a mechanical
property of the 100%/100% assumption, not a realistic dose-efficiency
estimate).

## 3. Formal indirect-effect evidence (Table 3, Section 5 of the task spec)

The three-arm decomposition (Arm A = no vaccine; Arm B = full-feedback
vaccine; Arm C = vaccine with transmission feedback disabled, using Arm A's
own lambda(t) fixed for every week) gives:

| Component | Definition | Median | 95% CrI |
|---|---|---|---|
| Total effect | A − B | 538,773 | [473,374, 567,351] |
| Direct-only | A − C | 82,295 | [75,823, 91,702] |
| Transmission-mediated indirect | C − B | 456,069 | [383,759, 490,061] |
| **Indirect fraction of total effect** | (C−B)/(A−B) | **84.7%** | **[80.8%, 86.6%]** |

The identity total = direct-only + indirect holds **exactly** (QA10, error =
0.00e+00), by construction of the fixed-lambda Arm C design.

**This is unambiguous, model-internal evidence of a transmission-mediated
indirect effect**: the large majority (≈85%) of the infections averted by
vaccinating 12-year-olds arise not from protecting those individuals
directly, but from reduced transmission pressure that also protects everyone
else. This is exactly what Arm C is built to isolate — it applies the same
vaccination to the same cohorts but never lets the resulting reduction in
infectiousness feed back into the force of infection, so any protection Arm C
fails to capture (relative to full-feedback Arm B) must be transmission-borne.

**Caveat (Section 8, "important interpretation rules")**: this 84.7% figure
is specific to the 100%-coverage/100%-VE stress test. Indirect (herd)
protection in renewal/SIR-type models scales non-linearly, and often
super-linearly at low-to-moderate coverage, with the fraction of the
population immunised. It should **not** be read as "84.7% of a realistic
vaccination programme's benefit will be indirect" — that requires re-running
Scenario 1 at realistic coverage/VE values, explicitly deferred (Section 8
below).

## 4. Age-specific effects (Table 2)

| Age group | Infections, no vaccine | Infections, vaccine | % reduction |
|---|---|---|---|
| 0–11 (never directly vaccinated in this window) | 428,670 | 348,500 | 18.85% |
| **12 (target age)** | 37,471 | 3,427 | **90.86%** |
| 13–17 (ages into target cohort's protection) | 201,317 | 114,426 | 43.18% |
| 18–64 | 1,563,869 | 1,269,587 | 18.94% |
| 65+ (never directly vaccinated) | 231,793 | 191,691 | 17.14% |

The target age group (12) shows the largest relative reduction (90.9%), as
expected for a 100%-VE stress test. But **every other age group also shows a
substantial reduction (17–43%)**, including 0–11 and 65+, neither of which
can ever receive the vaccine directly within a 2015–2021 window. Since Arm
A/C QA (Section 1, item 5) confirms 0–11 and 65+ have essentially zero
"direct-only" pathway to differ, their entire 17–19% reduction is indirect
protection — a second, independent confirmation of Section 3's formal
decomposition, this time by age group rather than by mechanism arm.

The elevated reduction in 13–17 (43.2%, compared to 18.9–18.9% in the
neighbouring 0–11 and 18–64 bands) reflects a mix of (a) genuine indirect
protection and (b) **direct protection carried forward through ageing**:
individuals vaccinated at 12 remain lifelong-protected in this stress test
and simply age into the 13–17 band over the following one to five years.
Table 2 and Figure 1C do not separate these two contributions within the
13–17 band — the arm-level decomposition in Section 3 is the only place they
are formally separated.

## 5. Does vaccination alter transmission dynamics? R0, R_eff, lambda, infectiousness

**R0(t): unchanged, exactly, by construction (QA7).** The vaccination
mechanism only reallocates susceptible mass between compartments (S → Uvac);
it never touches the harmonic-seasonal, climate-forced R0(t) baseline shared
by all three arms. This is the expected and required result — R0 is a
property of the transmission environment, not of who is immune within it.
Figure 2A shows a single line for all three arms.

**R_eff(t), lambda(t), infectiousness(t): altered, but not monotonically.**
During the two large epidemic years in this window (2016, 2017), all three
quantities are *lower* in the vaccine arm, as expected: fewer susceptibles
being effectively removed from transmission means less onward infection.
Concretely, median annual lambda drops from 6.9e-4 to 4.9e-4 in 2016 and from
2.13e-3 to 1.50e-3 in 2017, with matching reductions in R_eff and
infectiousness (Figure 2B–D).

**Unexpected, but mechanistically explainable finding:** in the four
subsequent, quieter years (2018–2021), the pattern *reverses* — median lambda,
R_eff, and infectiousness are all somewhat **higher** in the vaccine arm than
in the no-vaccine arm (e.g. lambda roughly 2–4× higher in 2019–2021, though
absolute levels are small in both arms during these low-transmission years).
This is a real epidemiological consequence of dampening the 2017 epidemic,
not a numerical artefact: by blunting the 2016–2017 wave, the vaccine arm
leaves a **larger pool of naturally (non-vaccine) susceptible individuals**
un-infected at the end of 2017 than the no-vaccine arm does (since fewer
people caught the disease during the big wave). That larger leftover
susceptible pool then fuels modestly larger secondary transmission in the
following quiet years. Annually, this shows up as *negative* "infections
averted" in 2018–2021 (roughly −25,000 to −62,000 infections per year, i.e.
the vaccine arm has slightly *more* infections than the no-vaccine arm in
those specific years) that partially offsets the very large positive effect
in 2016–2017. See Section 6 and Section 7 for the full annual breakdown —
the **net effect across all seven years remains strongly positive**
(538,773 infections averted); this is a re-timing/partial-rebound phenomenon
internal to a fixed, bounded analysis window, not evidence against a net
benefit.

## 6. Is each cohort vaccinated exactly once?

**Yes — QA4 passes exactly**, and this is also visually confirmed in Figure
3D: each of the birth cohorts shown (2003–2009) occupies exactly one,
non-overlapping ~1-calendar-year vaccination window at age 12, with no
cohort appearing more than once across draws or arms. This follows directly
from the simulator design (Section 4 of the task spec): vaccination acts only
on the weekly age(11)→age(12) ageing *inflow*, never on the standing age-12
stock, so an individual is only ever "new to age 12" in one week of their
life.

## 7. Is the effect concentrated around 2017, or spread across the whole window?

**Concentrated, and non-uniform in sign.** Breaking the total effect down by
calendar year (median across 200 draws):

| Year | Infections averted (median) | Direct-only | Indirect | % of total net effect |
|---|---|---|---|---|
| 2015 | ~30 | ~6 | ~24 | ~0.0% |
| 2016 | 151,419 | 13,387 | 137,998 | 28.0% |
| 2017 | 580,216 | 61,622 | 518,739 | 107.2% |
| 2018 | −25,211 | 1,376 | −26,554 | −4.7% |
| 2019 | −60,824 | 1,496 | −62,398 | −11.2% |
| 2020 | −62,303 | 1,547 | −63,840 | −11.5% |
| 2021 | −42,126 | 2,425 | −44,515 | −7.8% |

**Over 96% of the net effect is delivered in a single year, 2017** (Ceará's
largest epidemic in this window), with 2016 contributing a further 28%. The
2018–2021 "quiet years" collectively give back about 35 percentage points of
that total through the rebound mechanism described in Section 5, netting to
the reported 100% (538,773 infections). This is consistent with, and a
direct numerical restatement of, the lambda/R_eff pattern in Section 5: the
benefit of vaccination in this window is realised almost entirely by
blunting the one dominant epidemic wave, with a partial, smaller give-back
afterward. This pattern would very plausibly look different (larger net
benefit, smaller give-back) over a longer window or under recurring large
epidemics, since more of the "deferred" susceptibility would eventually be
consumed by future waves rather than sitting idle at the 2021 cutoff — but
testing that is explicitly out of scope here (no extension past 2021).

## 8. Are results driven by the 100%/100% stress-test assumption?

**Yes, substantially, and this must not be forgotten when reading any of the
above numbers.** Specifically:

- The **90.9% reduction in the target age-12 group** is close to what 100%
  coverage × 100% infection-blocking VE mechanically implies for direct
  protection; realistic vaccines (partial coverage, VE well under 100%, and
  typically not fully infection-blocking) would show a much smaller direct
  effect in this age band.
- The **84.7% indirect fraction** is a property of complete, permanent
  immunisation of an entire birth cohort every year. Herd-protection
  fractions in transmission models are usually non-linear in coverage — this
  number should not be linearly rescaled to, say, "50% coverage → ~42%
  indirect fraction." It requires rerunning the same pipeline (already built
  to support it — Section 4) with realistic coverage/VE.
- The **rebound/give-back pattern in 2018–2021** would likely be smaller
  (both in absolute infections and as a fraction of the total effect) under
  lower coverage or lower VE, since less susceptible mass would be "saved"
  from the 2017 wave in the first place.
- Doses administered (913,032) is coverage x eligible entrants by
  construction and scales linearly with any future coverage assumption.

None of the qualitative conclusions — that indirect effects exist, are
large relative to the direct effect, are delayed relative to direct effects,
and that R0 itself is unaffected — depend on the 100%/100% assumption. Only
the *magnitudes* do.

## 9. Is the simulator ready to explore evidence-based coverage and VE values?

**Yes, without further engineering work.** `04_run_full_feedback_vaccine_arm.R`
and `05_run_direct_only_diagnostic_arm.R` already take `coverage` and
`ve_infection` as parameters (currently hard-coded to the Scenario 1 stress
test values in `00_config.R$SCENARIOS$stresstest`); the generic simulator
(`02_generic_age_vaccine_simulator.R`) has no hard-coded dependence on those
values being 1.0. Running a realistic scenario requires only: (a) choosing
literature-based coverage/VE values (and optionally waning, not currently
modelled — protection is lifelong in the current code), (b) updating
`SCENARIOS` in `00_config.R`, and (c) re-running `06_execute_three_arms_all_draws.R`
onward. No changes to the transmission core, the age-cohort bookkeeping, or
the three-arm decomposition logic are needed.

## 10. Is the simulator ready for a 2026–2050 forward projection?

**Structurally, largely yes — but this question is answered here only on
readiness, per instruction, and no such projection has been run.** The prior
`36_climate_forced_v4_9` phase already built and unit-tested a mechanism-test
scaffold extending the demographic/vaccination bookkeeping to 2026–2050 (using
repeated 2025 R0/demography purely to exercise the code path, explicitly
flagged as non-forecast). Genuine forward projection would additionally
require: a real future climate forcing scenario (none exists yet — the
climate anomaly covariates are only defined against observed 2015–2021
climatology), a decision on how R0's harmonic/year-effect structure should be
extrapolated beyond the fitted period, and a chosen realistic vaccine
coverage/VE/waning profile (Section 9). None of these choices has been made,
and per the explicit stop rule for this task, none is made here.

## 11. Summary of unexpected behaviour

Two genuine, non-bug findings emerged during this analysis that a reader
should be aware of:

1. **Direct-protection carryover through ageing** (QA5, Section 4): under
   Arm C's diagnostic feedback-disabled design, the 13–17 and 18–64 age
   bands are *expected* to diverge from Arm A over time, because individuals
   vaccinated at 12 keep their (stress-test) lifelong protection as they age
   into these bands, independent of any transmission feedback. This is not
   part of the formal "indirect effect" (which requires feedback) but is
   part of the formal "direct-only" component (A−C) as defined by this
   task's own decomposition. Verified to start at exactly zero and grow
   smoothly on the ageing timescale (not an instantaneous jump), ruling out
   a bookkeeping bug.
2. **Post-epidemic transmission rebound** (Section 5, Section 7): in the four
   years after the dampened 2017 epidemic, lambda/R_eff/infectiousness are
   modestly *higher* in the vaccine arm than the no-vaccine arm, because
   blunting the 2017 wave leaves more naturally-susceptible (non-vaccinated)
   individuals available to sustain smaller subsequent waves. This produces
   small *negative* year-specific "infections averted" in 2018–2021 that
   partially offset the very large 2016–2017 benefit. The net effect over
   the full window remains strongly positive.

Both findings were investigated in detail (temporal-onset checks, magnitude
checks against known transmission dynamics) and are judged to be correct
properties of the model, not defects. Both are documented in code comments
(`07_validate_counterfactuals.R`, QA5 block) as well as here.

## 12. Assumptions disclosed

- All transmission parameters (harmonic seasonality, year effects, GI,
  climate multiplier) are frozen at their accepted climate-forced v4.9
  posterior values; no refitting was performed for this analysis.
- Force of infection is homogeneous across ages conditional on the
  aggregate susceptible fraction (no age-structured contact matrix); age
  only enters through population share and the vaccination-eligible ageing
  flow.
- Natural immunity is lifelong (S → Uinf is one-way, matching the accepted
  v4.9 model).
- Vaccine-derived protection in Scenario 1 is also modelled as lifelong and
  100% infection-blocking — a deliberate upper-bound stress test, not a
  claim about any real vaccine's durability or mechanism of action.
- No mortality differential by disease or vaccination status; deaths and
  births are allocated by population share only, exactly as in the accepted
  age-cohort reconstruction.
- The reporting window (2015–2021) is fixed and does not include the 2022+
  seed period; results describe only this bounded historical window and
  should not be extrapolated to it or beyond without further modelling
  (Section 10).

---
Generated as the final deliverable of the isolated analysis directory
`02_Script/40_renewal_model/37_ceara_historical_age12_vaccine_counterfactual/`.
No files outside this directory were modified.
