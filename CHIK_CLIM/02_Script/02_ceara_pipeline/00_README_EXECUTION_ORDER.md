# CÉARA PIPELINE

After opening:

    CHIK_CLIM.Rproj

Step 00 runs automatically (via `.Rprofile` -> `CHIK_CLIM/00_project_setup.R`):
project root, common packages, RStan options, and the `project_path()` /
`ROOT` / `DIR_*` helpers are all set up for you. Step 00 does not prepare
data, fit models, or run simulations.

Then the canonical order for the Ceará analysis is:

    01_input_preparation
    02_transmission_fitting
    03_model_diagnostics
    (04_sensitivity_identification -- not applicable to Ceara, skipped)
    (05_spatial_validation -- not applicable to Ceara, skipped)
    06_age_demographic_extension
    07_counterfactuals
    08_impact_summaries
    09_figures_reporting

Stage numbers stay fixed across every state pipeline (`02_ceara_pipeline`,
`03_bahia_pipeline`, ...) even where a given state has no content for a
slot — "stage 6 = age/demographic extension" means the same thing
everywhere. Ceará has no q-identification-sweep or spatial-model work of
its own (that lives in the frozen v4.1-v4.7 development history and was
never part of this pipeline), so stages 04 and 05 do not exist here.

Each stage has a `00_run_stage_XX.R` orchestrator that sources that
stage's scripts in the correct order and checks its prerequisites first —
if something upstream is missing, it stops with a message naming the exact
missing file and which stage to run to produce it, instead of a cryptic
`readRDS`/file-not-found error.

To run several stages in one call, use the top-level orchestrator:

```r
# full rebuild including Stan (only refits whatever is actually missing)
Sys.setenv(RUN_FROM_STAGE = "1", RUN_TO_STAGE = "9")
source("02_Script/02_ceara_pipeline/99_run_ceara_pipeline.R")

# reuse the accepted fit; rerun downstream analysis only (the safe default)
source("02_Script/02_ceara_pipeline/99_run_ceara_pipeline.R")

# figures only
Sys.setenv(RUN_FROM_STAGE = "9", RUN_TO_STAGE = "9")
source("02_Script/02_ceara_pipeline/99_run_ceara_pipeline.R")
```

The default (no environment variables set) runs **Stage 3 through Stage
9** — it never touches Stage 02, so it can never accidentally trigger an
expensive Stan refit.

## Stage summary

| Step | Folder | Main action | Refits Stan? | Main output |
|------|--------|-------------|:---:|-------------|
| 01 | `01_input_preparation` | Build & QC climate anomaly covariates (temperature/precipitation, 3 lag specs) | NO | `03_Output/tables/climate_forced_v4_9/climate_anomaly_covariates.csv` |
| 02 | `02_transmission_fitting` | Fit the climate-forced v4.9 model (reusing td14 init and the baseline v4.9 fit if present) | **YES** | `36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds` |
| 03 | `03_model_diagnostics` | Validate the accepted fit; plot the six-panel historical reconstruction | NO | Six-panel figure + fit diagnostics tables |
| 04 | *(not applicable)* | -- | -- | -- |
| 05 | *(not applicable)* | -- | -- | -- |
| 06 | `06_age_demographic_extension` | Build age population shares; reconstruct single-year age-specific S/U from the accepted fit; confirm the age-cohort simulator exactly reproduces the accepted fit with vaccination OFF | NO | `CE_age_population_shares.csv`, age-diagnostic arrays/figures, closed-loop replay validation |
| 07 | `07_counterfactuals` | Simulate the historical age-12 routine-vaccination counterfactual (3 arms x 200 posterior draws); vaccine engine unit tests, cohort eligibility QA, full A/B/C counterfactual QA (10 checks) | NO | `three_arm_results.rds`; QA diagnostic CSVs/markdown -- **fails loudly** on any substantive failure |
| 08 | `08_impact_summaries` | Summarise impacts (Tables 1-3) and calculate + validate NNV | NO | Impact/NNV tables |
| 09 | `09_figures_reporting` | Main, supplementary, and NNV publication figures | NO | `03_Output/figures/.../{main,supplementary,nnv}/` |

**Stage 02 is the ONLY stage that refits Stan.** Stages 03-09 all reuse the
accepted posterior. Stage 07 runs the vaccine counterfactual **simulation**
(the age-cohort engine driven by already-sampled posterior draws), which is
computationally nontrivial but is NOT Stan fitting. Stage 08 calculates
impacts and NNV purely from Stage 07's saved simulation output. Stage 09 is
plotting only.

## Scientific dependency graph

```
   raw epidemiological + demographic data
                   |
                   v
      01 input preparation (climate)
                   |
                   v
      02 historical Stan inference
                   |
                   | produces posterior R0(t), S(t), U(t), X(t)
                   v
         03 model diagnostics
                   |
                   v
      06 age/demographic extension
     (+ closed-loop replay validation)
                   |
                   v
         07 counterfactuals
     (+ vaccine engine validation)
                   |
                   v
         08 impact summaries
                   |
                   v
         09 figures/reporting
```

Age reconstruction and vaccination (Stage 06/07) are **not** separately
refitted Stan models — both are deterministic/simulation post-processing
of the single accepted Stage 02 posterior.

## Relationship to folders 36 and 37 (in `40_renewal_model/`)

- `40_renewal_model/36_climate_forced_v4_9/` and
  `40_renewal_model/37_ceara_historical_age12_vaccine_counterfactual/` keep
  their **outputs** (`outputs/`, `results/`, `tables/`, `diagnostics/`,
  `reports/`, and the corresponding `03_Output/figures/...` directories) —
  this refactor did not move or rename any existing result.
- Their `scripts/` subfolders are now **empty**: every active script that
  used to live there has moved here (see `REFACTOR_MIGRATION_MAP.csv` for
  the full old -> new mapping, including the later top-level relocation
  from `40_renewal_model/38_ceara_current_pipeline/` to this folder). A
  short `README_ACTIVE_SCRIPTS_MOVED.md` in each of those two folders
  points here.
- `36_climate_forced_v4_9/outputs/ce_canary/ce_climate_forced_canary_q0.05.rds`
  remains the single accepted upstream fit that Stages 03-09 all read (its
  eventual relocation to `03_Output/model_fits/ceara/current/` under the
  project's extension-based output rule is a separate, not-yet-completed
  step -- see the parent reorganisation notes).

## Relationship to the model-development history

Ceará's model-development lineage (`40_renewal_model/01_v1` through
`35_mato_grosso_v4_9_replication`) documents how the current model was
developed, replicated across other Brazilian states, or audited for
identifiability — it is frozen history and is **not** part of the current
Ceará execution path (destined for `90_development_archive/` in a later
step of the ongoing repository reorganisation). This pipeline sources a
handful of its functions directly where the current model still depends on
them (e.g. Stage 02 reuses the exact frozen v4.0/v4.9 fitting functions
unchanged), but you should not need to open or manually run those folders
yourself for the current analysis.

## Forcing a Stage 02 refit

Stage 02 checks each of its three required fits (td14 init, baseline
v4.9 q=0.05, climate-forced v4.9 q=0.05) for existence and reuses whatever
is already there. To force regenerating one that already exists, delete or
rename its `.rds` file, then either run Stage 02 directly or include Stage
02 in the top-level orchestrator's range — both require
`Sys.setenv(CEARA_ALLOW_REFIT = "TRUE")` (the top-level orchestrator sets
this automatically whenever Stage 02 falls inside `[RUN_FROM_STAGE,
RUN_TO_STAGE]`). A full refit is never triggered just because this file
was sourced.
