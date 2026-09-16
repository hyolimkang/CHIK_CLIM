# Bahia Global-q Local Consistency Audit

**Status: post-hoc plausibility/consistency audit only. No Stan model was fit or modified. q was not re-estimated. No spatial/municipality-specific q or susceptible compartment was added.**

## 0. Purpose

The current full Bahia model (frozen v4.9 + one constant global Bahia `q`
+ geo-adjusted multi-site serology) has clean HMC (Rhat=1.004, 0
divergences) and posterior `q` median = 0.136 (95% CrI 0.070-0.261), with
a correspondingly low state-wide cumulative infected fraction (~4.4% by
end-2025). This audit tests whether that result is (A) mostly arithmetic,
(B) reasonably consistent with municipality case burden and local
serology, or (C) difficult to reconcile locally -- using ONLY the existing
posterior, the validated municipality-week case panel, and the six
existing serosurveys. It is explicitly NOT a municipality-level
transmission model.

## 1-2. Validated inputs and reconstruction check

All inputs reused from already-validated sources (no rebuild): the
`chik_dlnm_panel_muni_week_2015_2025.rds` municipality-week panel (same
one validated in the spatial-turnover diagnostic), the current
`on_full2` global-q + serology posterior (`q`, `X[t]`, `S[t]`, `U[t]`
draws), and the six existing Bahia serosurveys with their already-audited
collection windows.

| Check | Result |
|---|---|
| Weeks | 573 |
| Total state cases (Stan input) | 92,498 |
| Total cases (municipality aggregation) | 92,498 |
| Absolute / relative discrepancy | 0 / 0 |
| Max weekly discrepancy | 0 |

**Confirmed exact reconstruction** -- proceeding.

U19-U22 are four distinct communities *within* Salvador municipality;
their municipality-level context is identical by construction and cannot
resolve within-Salvador heterogeneity (flagged per instruction).

## 3. State-level arithmetic check

| Quantity | Median | 50% CrI | 95% CrI |
|---|---|---|---|
| C_total / q / N_start (crude) | 0.0472 | 0.0374-0.0593 | 0.0246-0.0922 |
| C_total / q / N_mean (crude) | 0.0464 | 0.0367-0.0582 | 0.0241-0.0905 |
| C_total / q / N_end (crude) | 0.0458 | 0.0363-0.0575 | 0.0238-0.0894 |
| sum(X)/N_T (posterior, Stan) | 0.0460 | 0.0364-0.0579 | 0.0240-0.0904 |
| final immune fraction U[T]/(S[T]+U[T]) (posterior, Stan) | 0.0438 | 0.0347-0.0551 | 0.0228-0.0861 |

**Ratios** (near 1.0 = arithmetic, not a Stan artefact):

| Ratio | Median | 95% CrI |
|---|---|---|
| sum(X) / (C_total/q) | **1.003** | 0.933-1.078 |
| final immune fraction / (C_total/q/N_mean) | **0.945** | 0.879-1.015 |

**Both ratios are essentially 1.0.** The low state-wide cumulative
infection estimate is almost entirely explained by simple arithmetic
(reported cases / q / population), confirming this is not a hidden
transmission-dynamics or Stan-specific artefact.

## 4. Fixed-q arithmetic scale table (not a fitted comparison)

| q | Total implied infections | Implied attack fraction | Implied susceptible fraction |
|---|---|---|---|
| 0.05 | 1,849,960 | 12.58% | 87.42% |
| 0.10 | 924,980 | 6.29% | 93.71% |
| 0.136 | 680,132 | 4.62% | 95.38% |
| 0.15 | 616,653 | 4.19% | 95.81% |
| 0.20 | 462,490 | 3.14% | 96.86% |
| 0.25 | 369,992 | 2.52% | 97.48% |
| 0.30 | 308,327 | 2.10% | 97.90% |
| 0.359 | 257,655 | 1.75% | 98.25% |

## 5. Municipality-level plausibility audit (crude, not a Stan estimate)

At **q = 0.136**, across all 414 municipalities and all 6 checkpoints
(2016-2025): **maximum crude A_i = 0.968 -- zero municipality-checkpoint
combinations exceed 1.0.** At end-2025 specifically: median A_i = 0.0082,
population-weighted mean = 0.046, 95th pct = 0.178, 2.9% of municipalities
exceed A_i > 0.25, 0.97% exceed 0.50, and none exceed 1.0.

At q = 0.05, by contrast, 4 municipality-checkpoint combinations (all at
2020 or 2022, i.e. peak-epidemic small municipalities) exceed A_i = 1.0 --
a genuine falsification signal at that end of the grid, though a small
absolute count (4 of 2,484 municipality x checkpoint x q=0.05
combinations) consistent with a few small, hard-hit municipalities rather
than a systemic failure.

**Conclusion: q=0.136 produces no biologically impossible municipality-level
crude attack fractions anywhere in Bahia over 2015-2025.**

## 6. Maps

`bahia_municipality_implied_attack_q005.png` /
`..._q0136.png` / `..._q0359.png` show the crude 2025 attack-fraction map
under each q (capped at 1.0 for colour-scale readability).
`bahia_municipality_A_gt1_flag_map.png` shows that only a handful of
municipalities are ever flagged, and only at q=0.05.

## 7-10. Serology-window alignment and the central finding of this audit

Aligning each survey's municipality (U03 -> Riachão do Jacuípe; U04,
U19-U22 -> Salvador) to its exact collection window (same window weights
as the Stan model) gives a **municipality-level cumulative reported
incidence at survey time (`R_survey`)** that is dramatically smaller than
the observed local seroprevalence:

| Survey | Municipality | Window | R_survey (muni cumulative reported incidence) | Observed seroprevalence |
|---|---|---|---|---|
| U03 | Riachão do Jacuípe | 2016-04 | **0.0000** | 20.0% |
| U04 | Salvador | 2016-11 to 2017-02 | 0.0000256 (0.0026%) | 11.8% |
| U19 | Salvador | 2018-03 to 2018-10 | 0.0000568 (0.0057%) | 6.1% |
| U20 | Salvador | 2018-03 to 2018-10 | 0.0000568 (0.0057%) | 4.7% |
| U21 | Salvador | 2018-03 to 2018-10 | 0.0000568 (0.0057%) | 22.6% |
| U22 | Salvador | 2018-03 to 2018-10 | 0.0000568 (0.0057%) | 4.4% |

**Crude serology-implied q** (`R_survey / observed_prevalence`, Section 9,
Wilson-CI-propagated, untruncated):

| Survey | q_implied | Untruncated range |
|---|---|---|
| U03 | 0.000 | [0.000, 0.000] |
| U04 | 0.000217 | [0.000191, 0.000246] |
| U19 | 0.000929 | [0.000631, 0.001383] |
| U20 | 0.001197 | [0.000750, 0.001931] |
| U21 | 0.000251 | [0.000206, 0.000311] |
| U22 | 0.001303 | [0.000776, 0.002214] |

**Every single survey implies a q roughly 100-500x smaller than the
model's global-q posterior median (0.136), and even smaller relative to
the entire fixed-q grid (0.05-0.359).** Section 12 propagates the FULL
q-posterior (not just the median) through each survey's `R_survey`, and
the discrepancy remains 2-3 orders of magnitude even at the 95% CrI
(`bahia_global_q_posterior_vs_local_serology.png`; e.g. U19's
q-posterior-propagated case-implied fraction is 0.042% median [0.022%,
0.082%] vs an observed 6.1% -- roughly 100x too low even at the upper
credible bound).

**Section 10 -- U19-U22 within-Salvador heterogeneity:** all four surveys
share an *identical* `R_survey` (Salvador's municipality-wide trajectory)
by construction, yet observed community seroprevalence ranges 4.4%-22.6%
(a 5.2x spread) and q_implied consequently ranges 0.00025-0.0013 purely
from within-municipality community variation. **This variation must NOT be
read as evidence that q differs among these four communities** -- it is
direct, clean evidence that municipality (or state) aggregate case burden
cannot resolve, or substitute for, hyperlocal community-level exposure.

**Interpretation (cautious, per the audited-before-concluding principle):**
this gap is far too large and far too systematic (6 of 6 surveys, 2-3
orders of magnitude) to be a rounding or Stan artefact, but it should NOT
be read as "q=0.136 is impossible." The most plausible, non-exotic
explanations are:

1. **Severe sub-municipality reporting/access heterogeneity.** Pau da
   Lima, Alto do Cabrito, Marechal Rondon, Nova Constituinte, and Rio Sena
   are peripheral, lower-income Salvador communities; well-documented
   healthcare-access gradients within large Brazilian cities plausibly
   mean their contribution to Salvador's *aggregate* SINAN case count is
   disproportionately tiny relative to their true local exposure -- this
   is exactly the kind of heterogeneity `eta_geo` was designed to absorb,
   and is consistent with the large `eta_geo` values already estimated
   (median 0.9-3.0 log-odds) and with the sub-regional turnover already
   documented in `BAHIA_SPATIAL_TURNOVER_DIAGNOSTIC.md`.
2. **Pre-2015 exposure (Section 11).** No 2014 Bahia chikungunya
   surveillance data exist anywhere in this repository (earliest raw
   SINAN file is `CHIKBR15`). Chikungunya was already circulating in
   Brazil from 2014, so 2016-2018 seroprevalence may partly reflect
   pre-2015 infections that the model's 2015-start case series cannot
   capture at all. We do **not** estimate the magnitude of this bias --
   only note it is directionally consistent with q_implied being biased
   downward (true local q could be somewhat higher than computed here).
3. Possible case-attribution practices (e.g. notification tied to
   treating facility rather than residence) that could further dilute a
   peripheral community's true burden out of Salvador's aggregate count --
   plausible but not verifiable with the data used in this audit.

None of these fully closes a 100-500x gap on their own, and this audit
cannot apportion how much each contributes.

## 11. Pre-2015 exposure

Confirmed: no 2014 Bahia chikungunya file exists in this repository (raw
SINAN files start at `CHIKBR15.csv.zip`). Stated per instruction, without
estimating a magnitude: 2016-2018 seroprevalence may include pre-2015
infections, so q_implied values computed from 2015+ reported cases alone
could be downward-biased to an unknown degree.

## 12. Posterior propagation

`bahia_global_q_posterior_local_attack_summary.csv` (survey level, full
q-posterior) and `bahia_global_q_posterior_municipality_attack_2025.csv`
(all 414 municipalities, full q-posterior). Across all municipalities, 4
have a posterior-median A_i(2025) > 0.5, and 4 have a 95% CrI upper bound
exceeding 1.0 (i.e. cannot fully rule out implausible attack fractions in
a small number of small municipalities once q-uncertainty is
propagated) -- worth flagging, not alarming given the small count.

## 13. Three key questions

**Q1 -- Is the low Bahia-wide cumulative infection fraction largely
arithmetic (reported cases / q / population)?**
**Yes.** Section 3's ratios are 1.003 and 0.945 (both ~1.0) -- confirmed
directly, not inferred.

**Q2 -- Is q≈0.136 broadly compatible with municipality case burden and
the six local serosurveys?**
**Compatible with municipality-wide case burden (Section 5: no
municipality-checkpoint combination is biologically impossible at
q=0.136), but NOT compatible with the six local serosurveys at their
survey-specific municipality and time window** (Sections 7-10, 12: every
survey implies a q 100-500x smaller than 0.136, and this persists across
the full q-posterior, not just its median).

**Q3 -- Does the local evidence support ONE constant Bahia-wide q, or does
apparent compatibility vary too much by geography?**
**It varies enormously by geography and by spatial scale.** The
municipality/state-wide arithmetic check (Section 5) shows no
inconsistency; the hyperlocal survey-site check (Sections 7-10) shows a
2-3-orders-of-magnitude inconsistency, uniformly across all six surveys.
Combined with `BAHIA_SPATIAL_TURNOVER_DIAGNOSTIC.md`'s MIXED finding
(substantial sub-regional turnover in which municipalities drive each
statewide wave), the evidence points toward real, large, sub-municipality
heterogeneity that a single state-wide q cannot represent -- without
positively falsifying q≈0.136 as a *state-level* aggregate quantity.

## 14. Classification

**B. GLOBAL-q WEAKLY SUPPORTED**

q≈0.136 is not clearly impossible: it produces no biologically implausible
municipality-wide cumulative attack fractions anywhere in Bahia over
2015-2025 (Sections 3-6 all pass). However, local evidence is highly
variable and does not strongly validate one state-wide q: all six
serosurveys, individually and under full posterior propagation, imply a
locally-required q roughly 100-500x smaller than the model's estimate,
with the discrepancy driven by identifiable (though not fully
quantifiable) sub-municipality heterogeneity and possible pre-2015
exposure rather than an obvious data error. This is consistent with, and
strengthens, the MIXED classification from the companion spatial-turnover
diagnostic: geography matters more than a single constant q can
represent, even though the state-level aggregate is not itself falsified.

No new q was selected or estimated as part of this classification.

## 15. Output files

Tables (`03_Output/tables/bahia_global_q_local_consistency_audit/`):
`bahia_reconstruction_check.csv`, `bahia_state_arithmetic_q_audit.csv`,
`bahia_state_arithmetic_q_audit_fixed_grid.csv`,
`bahia_municipality_implied_attack_by_q.csv`,
`bahia_municipality_implied_attack_summary.csv`,
`bahia_municipality_flagged_A_gt1.csv`,
`bahia_serology_case_implied_attack_by_q.csv`,
`bahia_serology_implied_q_audit.csv`,
`bahia_global_q_posterior_local_attack_summary.csv`,
`bahia_global_q_posterior_municipality_attack_2025.csv`.

Figures (`03_Output/figures/bahia_global_q_local_consistency_audit/`):
`bahia_state_attack_fraction_vs_q.png`,
`bahia_municipality_implied_attack_q005.png`,
`bahia_municipality_implied_attack_q0136.png`,
`bahia_municipality_implied_attack_q0359.png`,
`bahia_municipality_A_gt1_flag_map.png`,
`bahia_serology_observed_vs_case_implied_by_q.png`,
`bahia_serology_crude_implied_q.png`,
`bahia_global_q_posterior_vs_local_serology.png`.

## 16. Stop rule (honoured)

No changes were made to Stan, q, priors, spatial compartments, or allfoi.
Returned for scientific review before any further model development.
