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

Outputs live under `03_Output/05_rio_de_janeiro_pipeline/`: fitted Stan
objects in `model_fits/v4_9_replication/` (state and city fits together),
tables/figures in `tables|figures/rio_de_janeiro_v4_9_replication/`, and
design/result reports in `reports/`.

See `REFACTOR_MIGRATION_MAP.csv` (in `02_ceara_pipeline/`, extended with a
Rio de Janeiro section) for the full old -> new script mapping.
