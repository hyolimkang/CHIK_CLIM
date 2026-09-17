# 02_Script

All human-written source code for the CHIK_CLIM project: `.R` scripts and
active `.stan` model sources. Generated artifacts (`.rds` fits, tables,
figures, reports) never live here -- see `03_Output/README.md`.

```
02_Script/
├── 00_shared/                    <- code shared by 2+ pipelines (functions, active Stan models)
├── 01_data_pipeline/             <- fetch/build/clean raw and processed data (writes to 01_Data)
├── 02_ceara_pipeline/            <- Ceara: accepted baseline, climate-forced extension, vaccine counterfactual
├── 03_bahia_pipeline/            <- Bahia: v4.9 replication + spatial/identification diagnostics
├── 04_pernambuco_pipeline/       <- Pernambuco: v4.9 replication + spatial validation
├── 05_rio_de_janeiro_pipeline/   <- Rio de Janeiro (state + city): v4.9 replication
├── 06_mato_grosso_pipeline/      <- Mato Grosso: v4.9 replication
├── 07_national_pipeline/         <- Brazil-wide: episode census, wave analysis, susceptibility, annual FOI
└── 90_development_archive/       <- frozen model-development history (see its own README.md)
```

## The numbered pipelines

Each of `02_ceara_pipeline` through `07_national_pipeline` is a sequence of
numbered stage folders (`01_input_preparation`, `02_transmission_fitting`,
...), each with a `00_run_stage_XX.R` orchestrator that checks its own
prerequisites before running, plus a top-level `99_run_<pipeline>.R`
orchestrator and a `00_README_EXECUTION_ORDER.md` with the exact execution
order, inputs, and outputs for that pipeline. `01_data_pipeline` follows
the same numbered-stage convention but has no `99_run_*` (its stages are
largely independent fetch/build scripts).

**A Stan refit is never triggered merely by sourcing an orchestrator.**
Every pipeline gates any refit behind an explicit, pipeline-specific
environment variable (`CEARA_ALLOW_REFIT`, `BAHIA_ALLOW_REFIT`,
`PERNAMBUCO_ALLOW_REFIT`, `RIO_DE_JANEIRO_ALLOW_REFIT`,
`MATO_GROSSO_ALLOW_REFIT`, `NATIONAL_ALLOW_REFIT`) that is never set
implicitly -- refitting is always a separate, explicit decision.

## 00_shared

Code depended on by more than one pipeline: reusable R functions, and
active Stan model sources (`00_shared/stan/`, the one canonical Stan
source location -- see `00_shared/stan/README.md` and
`00_shared/stan/MODEL_REGISTRY.md`). See `00_shared/README.md` for the
full breakdown, including why `legacy_model_functions/` is live rather
than archived.

## 90_development_archive

Records model development, superseded analyses, and legacy workflows.
**Current production pipelines must not depend on it.** See
`90_development_archive/README.md`.

## Outputs

Every pipeline's outputs live under the correspondingly-numbered
`03_Output/<pipeline>/` folder (e.g. `02_ceara_pipeline` outputs live
under `03_Output/02_ceara_pipeline/`) -- see `03_Output/README.md`.
