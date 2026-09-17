# Scientific observations noted during the structural refactor -- NOT changed

This refactor (creation of `38_ceara_current_pipeline/`) is a file-location
and execution-order change only. Per the task's own instruction, anything
scientifically questionable noticed along the way is recorded here rather
than modified.

## 1. Stage 01 (climate preparation) has an implicit dependency on part of Stage 02

`01_input_climate_preparation/02_check_climate_orthogonality_and_priors.R`
(originally `36_climate_forced_v4_9/scripts/02_prior_predictive_and_orthogonality_check.R`)
reads the **baseline (non-climate) v4.9 q=0.05 fit** — the output this
refactor now labels "Stage 02B" — to obtain the exact harmonic-regressor
values and week alignment used at fitting time, so the orthogonality check
compares against what was actually fit rather than a re-derivation.

This means Stage 01, as numbered in this pipeline, is not purely upstream
of Stage 02: script 2 of Stage 01 cannot run until the Stage 02B fit
already exists. This is not a new problem introduced by the refactor — the
original `36_climate_forced_v4_9/scripts/` execution order had exactly the
same dependency, because the baseline v4.9 fit already existed in the
project (from `28_v4_9_hierarchical_seasonality/`) before the climate-anomaly
work began. `01_input_climate_preparation/00_run_stage_01.R` checks for
this file explicitly and stops with a clear message pointing at Stage 02
if it is absent, rather than silently failing with a cryptic
`readRDS`/file-not-found error.

No change made: the dependency is documented, not restructured.

## 2. Previously documented, unrelated-to-this-refactor scientific findings

These were identified and already written up in earlier phases of this
project (not discovered during this refactor, and not touched by it):

- **Post-epidemic transmission rebound / later negative annual vaccine
  effects** (2018-2021 in the historical age-12 counterfactual): documented
  in `37_ceara_historical_age12_vaccine_counterfactual/reports/HISTORICAL_AGE12_VACCINE_COUNTERFACTUAL_ASSESSMENT.md`
  (Section 5/7) and `.../reports/NNV_ASSESSMENT.md` (Section 5). Confirmed
  as a genuine SIR-type dynamic (blunting one large epidemic leaves more
  natural susceptibles for smaller subsequent waves), not a bug.
- **Direct-protection carryover through ageing** in the Arm C diagnostic
  (13-17/18-64 bands diverging from Arm A even with feedback disabled):
  documented in `07_vaccine_validation/03_validate_historical_counterfactual.R`'s
  QA5 block and the same assessment report.
- **NNV's strong sensitivity to evaluation horizon**, including a year
  (2015) flagged as producing an unstable/near-meaningless cumulative NNV:
  documented in `NNV_ASSESSMENT.md` Sections 5 and 9.
- **q (ascertainment) is fixed, not estimated**, in every fit used by this
  pipeline (q=0.05 throughout) — a known, disclosed modelling choice
  carried through from the v4.0-v4.9 development history
  (`17_v4_0_minimal_no_vaccine/MODEL_DEVELOPMENT.md` and the v4.1-v4.7
  q-identifiability audit folders), not something this refactor evaluates
  or changes.
- **Deterministic weekly importation and fractional (continuous-rate)
  ageing** are structural modelling choices baked into the frozen v4.9
  Stan model and the age-cohort simulator respectively; both predate this
  refactor and are unchanged by it.

None of the above were modified, re-derived, or re-litigated as part of
this structural refactor.
