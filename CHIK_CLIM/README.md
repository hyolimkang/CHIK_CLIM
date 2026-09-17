# CHIK_CLIM

Chikungunya transmission dynamics in Brazil, via climate-sensitive renewal
modelling: susceptibility reconstruction, regional replication across five
states, and a vaccination counterfactual analysis.

- **Chikungunya transmission dynamics in Brazil** -- weekly renewal-model
  reconstruction of transmission (R0, susceptibility, recurrence) from
  reported case data.
- **Climate-sensitive renewal modelling** -- temperature/precipitation
  anomaly covariates extending the accepted Ceara baseline.
- **Susceptibility reconstruction** -- episode- and annual-resolution
  dynamic FOI/susceptibility models, nationally and per state.
- **Regional replication** -- the accepted v4.9 model structure replicated
  against Bahia, Pernambuco, Rio de Janeiro (state + city), and Mato
  Grosso.
- **Vaccination counterfactual analysis** -- historical age-12 routine
  vaccination impact and number-needed-to-vaccinate for Ceara.

## Quick start

1. Open `CHIK_CLIM.Rproj`.
2. `.Rprofile` automatically sources `00_project_setup.R`, which locates
   the project root, loads common packages, configures RStan, and defines
   `project_path()`/`output_path()`/`DIR_*`/`OUT_*` helpers -- no manual
   `setwd()`/`library()` calls needed.
3. For current analyses, use:
   - `02_Script/01_data_pipeline`
   - `02_Script/02_ceara_pipeline`
   - `02_Script/03_bahia_pipeline`
   - `02_Script/04_pernambuco_pipeline`
   - `02_Script/05_rio_de_janeiro_pipeline`
   - `02_Script/06_mato_grosso_pipeline`
   - `02_Script/07_national_pipeline`
4. Do **NOT** execute `02_Script/90_development_archive` as part of
   current production analysis -- it is frozen history (see below).

## Project structure

```
CHIK_CLIM/
├── 00_project_setup.R   <- sourced automatically on project open; ROOT + path helpers only
├── 01_Data/              <- data and model inputs
├── 02_Script/            <- human-written source code
└── 03_Output/            <- generated analytical artifacts
```

The rule: **`01_Data` = inputs, `02_Script` = code, `03_Output` =
everything code generates.** Nothing under `03_Output` is hand-written;
nothing under `02_Script` should need hand-editing to match a specific
saved result.

## Current production code vs. development/archive

`02_Script/00_shared` through `02_Script/07_national_pipeline` are the
current, executable production pipelines -- see `02_Script/README.md` for
the full architecture and `02_Script/90_development_archive/README.md`
for what "frozen history" means here and its few documented exceptions.
`02_Script/00_shared/stan/MODEL_REGISTRY.md` lists every Stan model's
current-vs-legacy status individually.

## Pipeline execution order and outputs

Each pipeline's exact execution order, Stan-fitting stages, expected
inputs/outputs, and sensitivity analyses are documented in its own
`00_README_EXECUTION_ORDER.md` (not duplicated here):

- [`02_Script/01_data_pipeline/00_README_EXECUTION_ORDER.md`](02_Script/01_data_pipeline/00_README_EXECUTION_ORDER.md)
- [`02_Script/02_ceara_pipeline/00_README_EXECUTION_ORDER.md`](02_Script/02_ceara_pipeline/00_README_EXECUTION_ORDER.md)
- [`02_Script/03_bahia_pipeline/00_README_EXECUTION_ORDER.md`](02_Script/03_bahia_pipeline/00_README_EXECUTION_ORDER.md)
- [`02_Script/04_pernambuco_pipeline/00_README_EXECUTION_ORDER.md`](02_Script/04_pernambuco_pipeline/00_README_EXECUTION_ORDER.md)
- [`02_Script/05_rio_de_janeiro_pipeline/00_README_EXECUTION_ORDER.md`](02_Script/05_rio_de_janeiro_pipeline/00_README_EXECUTION_ORDER.md)
- [`02_Script/06_mato_grosso_pipeline/00_README_EXECUTION_ORDER.md`](02_Script/06_mato_grosso_pipeline/00_README_EXECUTION_ORDER.md)
- [`02_Script/07_national_pipeline/00_README_EXECUTION_ORDER.md`](02_Script/07_national_pipeline/00_README_EXECUTION_ORDER.md)

See also [`REPOSITORY_ARCHITECTURE.md`](REPOSITORY_ARCHITECTURE.md) for
the full `02_Script` <-> `03_Output` correspondence and file-type
placement rules, and
[`REPOSITORY_FINAL_CLEANUP_MANIFEST.csv`](REPOSITORY_FINAL_CLEANUP_MANIFEST.csv)
for the complete record of this reorganisation.

Note: `analysis_plan.md`, referenced in earlier project notes, is not
currently present at the project root.
