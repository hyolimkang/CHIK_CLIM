# RIO DE JANEIRO PIPELINE

Consolidated from `40_renewal_model/34_rio_de_janeiro_v4_9_replication/` --
no scientific logic changed, only location (git mv) and internal path
references. State and city are both inside this one pipeline (there is no
separate city folder in the source material).

Stage numbers stay fixed across every state pipeline; Rio de Janeiro has no
sensitivity/spatial/age/vaccine/counterfactual/impact work, so stages 04-09
do not exist here.

| Step | Folder | Main action | Refits Stan? |
|---|---|---|:---:|
| 01 | `01_input_preparation` | State data/episode/input audits + city data audit | NO |
| 02 | `02_transmission_fitting` | State global-q case-only fit (+ ad098 rescue), targeted fixed-q sweep (q=0.010-0.020), city global-q case-only fit, city U10-serology fit (+ ad099 rescue) | **YES** |
| 03 | `03_model_diagnostics` | Diagnostics/figures for state and city fits, state-vs-city-U10 external plausibility check, city A/B final comparison | NO |

```r
# safe default (Stage 3 only, never touches Stan)
source("02_Script/05_rio_de_janeiro_pipeline/99_run_rio_de_janeiro_pipeline.R")

# include fitting (only refits whatever is actually missing)
Sys.setenv(RUN_FROM_STAGE = "1", RUN_TO_STAGE = "3")
source("02_Script/05_rio_de_janeiro_pipeline/99_run_rio_de_janeiro_pipeline.R")
```

Outputs remain at their original locations under
`40_renewal_model/34_rio_de_janeiro_v4_9_replication/outputs/` and
`03_Output/tables|figures/rio_de_janeiro_v4_9_replication/{,city}` -- this
refactor did not move or rename any existing result (a separate,
not-yet-completed step will migrate these to
`03_Output/model_fits/rio_de_janeiro/...` under the project's
extension-based output rule).

See `REFACTOR_MIGRATION_MAP.csv` (in `02_ceara_pipeline/`, extended with a
Rio de Janeiro section) for the full old -> new script mapping.
