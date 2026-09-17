# MATO GROSSO PIPELINE

Consolidated from `40_renewal_model/35_mato_grosso_v4_9_replication/` -- no
scientific logic changed, only location (git mv) and internal path
references.

Stage numbers stay fixed across every state pipeline; Mato Grosso has no
sensitivity/spatial/age/vaccine/counterfactual/impact work, so stages 04-09
do not exist here.

| Step | Folder | Main action | Refits Stan? |
|---|---|---|:---:|
| 01 | `01_input_preparation` | Data/episode/input audits | NO |
| 02 | `02_transmission_fitting` | Global-q case-only fit (+ ad098 rescue), fixed-q profile (q=0.020-0.070) | **YES** |
| 03 | `03_model_diagnostics` | Diagnostics/figures for the global-q fit, HMC diagnostics, six-panel plot, fixed-q profile analysis | NO |

```r
# safe default (Stage 3 only, never touches Stan)
source("02_Script/06_mato_grosso_pipeline/99_run_mato_grosso_pipeline.R")

# include fitting (only refits whatever is actually missing)
Sys.setenv(RUN_FROM_STAGE = "1", RUN_TO_STAGE = "3")
source("02_Script/06_mato_grosso_pipeline/99_run_mato_grosso_pipeline.R")
```

Outputs remain at their original locations under
`40_renewal_model/35_mato_grosso_v4_9_replication/outputs/` and
`03_Output/tables|figures/mato_grosso_v4_9_replication/` -- this refactor
did not move or rename any existing result (a separate, not-yet-completed
step will migrate these to `03_Output/model_fits/mato_grosso/...` under the
project's extension-based output rule).

See `REFACTOR_MIGRATION_MAP.csv` (in `02_ceara_pipeline/`, extended with a
Mato Grosso section) for the full old -> new script mapping.
