# PERNAMBUCO PIPELINE

Consolidated from `40_renewal_model/33_pernambuco_v4_9_replication/` — no
scientific logic changed, only location (git mv) and internal path
references.

Stage numbers stay fixed across every state pipeline; Pernambuco has no
age/vaccine/counterfactual/impact work, so stages 06-09 do not exist here.

| Step | Folder | Main action | Refits Stan? |
|---|---|---|:---:|
| 01 | `01_input_preparation` | Pernambuco data audit | NO |
| 02 | `02_transmission_fitting` | Global-q model variants (modelA case-only, modelB U14 serology, modelB U14 full) + fixed-q sweep (q=0.05-0.30) | **YES** |
| 03 | `03_model_diagnostics` | Diagnostics/figures for the global-q and fixed-q-sweep fits, CE/BA/PE cross-state comparison | NO |
| 04 | `04_sensitivity_identification` | Global-q and curved ridge-reparameterisation identification diagnostics (computationally improved coordinate systems for the frozen v4.9 model, not new models) | **YES**† |
| 05 | `05_spatial_validation` | Municipality-level turnover diagnostic, 5-strata/Spatial-1 hierarchical-region attempt, Spatial-0 common-q baseline, Q2/regional-R feasibility stops | **YES**† |

† Stage 04/05's fits (modelC/D ridge/curved reparam, 5-strata pilot) are
already on disk except one documented halt (see below); `PERNAMBUCO_ALLOW_REFIT`
gates any refit exactly as in Stage 02, and the top-level orchestrator only
auto-enables that gate for Stage 02 — a Stage 04/05 refit requires the user
to `Sys.setenv(PERNAMBUCO_ALLOW_REFIT = "TRUE")` themselves first.

## Halted work, kept live here (not archived)

Per explicit decision, Pernambuco's stopped spatial-model attempts stay
inside the live `05_spatial_validation` stage, documented as halted/failed
in place, rather than being moved to `90_development_archive`:

- **5-strata / "Spatial-1" hierarchical-region model** (steps 09-16): the
  4-chain pilot's chain 1 hit a severe, unrecovered warmup pathology
  (78/1500 iterations vs. 1500/1500 for chains 2-4) and the run was stopped
  before a final fit was ever saved. See `PE_SPATIAL1_GEOMETRY_DIAGNOSTIC.md`.
  Steps 13-15 (HMC/structural/PPC diagnostics that need the missing bundled
  fit) are skipped by the stage runner; step 16 is the diagnostic that
  explains the failure from the raw (partial) chain CSVs, which do exist.
- **Q2 (differential-ascertainment) feasibility** (step 21): a
  deterministic, no-Stan-refit feasibility scan concluded the Recife
  case/serology conflict is not solvable by (q_Recife, q_rest) alone; no Q2
  Stan model was ever built. See `PE_SPATIAL0_Q2_FEASIBILITY.md`.
- **Regional transmission-intercept feasibility** (step 22): a
  deterministic, no-Stan-refit profile/optimisation test concluded absolute
  regional susceptibility inference should be stopped. See
  `PE_REGIONAL_R_FEASIBILITY.md`.

Spatial-0 (steps 17-20, the "common-q" model-criticism baseline) is
**not** halted — it is a completed track whose accepted artifact is its raw
per-chain CSVs (`outputs/spatial0/chains/`), never bundled into a single
`.rds` by design, since the branch-agreement diagnostic (step 19) needs the
separate LOW/INTERMEDIATE/HIGH-initialised chains directly.

```r
# safe default (Stage 3-5, never touches Stan)
source("02_Script/04_pernambuco_pipeline/99_run_pernambuco_pipeline.R")

# include fitting (only refits whatever is actually missing)
Sys.setenv(RUN_FROM_STAGE = "1", RUN_TO_STAGE = "5")
source("02_Script/04_pernambuco_pipeline/99_run_pernambuco_pipeline.R")
```

Outputs remain at their original locations under
`40_renewal_model/33_pernambuco_v4_9_replication/outputs*/` and
`03_Output/tables|figures/pernambuco_v4_9_replication/` — this refactor did
not move or rename any existing result (a separate, not-yet-completed step
will migrate these to `03_Output/model_fits/pernambuco/...` under the
project's extension-based output rule).

See `REFACTOR_MIGRATION_MAP.csv` (in `02_ceara_pipeline/`, extended with a
Pernambuco section) for the full old -> new script mapping.
