# BAHIA PIPELINE

Consolidated from `40_renewal_model/30_bahia_v4_9_replication/`,
`31_bahia_spatial_turnover_diagnostic/`, and
`32_bahia_global_q_local_consistency_audit/` — no scientific logic changed,
only location (git mv) and internal path references.

Stage numbers stay fixed across every state pipeline; Bahia has no
age/vaccine/counterfactual/impact work, so stages 06-09 do not exist here.

| Step | Folder | Main action | Refits Stan? |
|---|---|---|:---:|
| 01 | `01_input_preparation` | Bahia data audit | NO |
| 02 | `02_transmission_fitting` | Case-only v4.9 (q=0.05 primary, q-sweep 0.10-0.30 available), multisite serology fit, global-q multisite fit | **YES** |
| 03 | `03_model_diagnostics` | Six-panel plots, wave PPC, Ceará-vs-Bahia comparison, multisite/global-q diagnostics, sero-ablation comparison | NO |
| 04 | `04_sensitivity_identification` | Global-q (estimated, not fixed) identification experiment: prior-predictive check, identification diagnostics, plots | NO |
| 05 | `05_spatial_validation` | Municipality-level turnover diagnostic (01-06) + global-q local consistency audit (07-13) | NO |

```r
# safe default (Stage 3-5, never touches Stan)
source("02_Script/03_bahia_pipeline/99_run_bahia_pipeline.R")

# include fitting (only refits whatever is actually missing)
Sys.setenv(RUN_FROM_STAGE = "1", RUN_TO_STAGE = "5")
source("02_Script/03_bahia_pipeline/99_run_bahia_pipeline.R")
```

Outputs remain at their original locations under
`40_renewal_model/30_bahia_v4_9_replication/outputs*/` and
`03_Output/tables|figures/{bahia_v4_9_replication,bahia_spatial_turnover_diagnostic,bahia_global_q_local_consistency_audit}/`
— this refactor did not move or rename any existing result (a separate,
not-yet-completed step will migrate these to `03_Output/model_fits/bahia/...`
under the project's extension-based output rule).

See `REFACTOR_MIGRATION_MAP.csv` (in `02_ceara_pipeline/`, extended with a
Bahia section) for the full old -> new script mapping.
