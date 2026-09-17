# 03_Output

Every file here is a **generated artifact** -- regenerable by rerunning
the corresponding `02_Script/` code against `01_Data/`. Nothing here is
hand-written, and nothing under `02_Script/` should ever need to be
manually edited to match something saved here.

Outputs mirror `02_Script`'s pipeline numbering one-to-one:

```
03_Output/
├── 02_ceara_pipeline/
│   ├── model_fits/{initialization,baseline,climate_forced}/
│   ├── diagnostics/  results/  tables/  reports/
│   └── figures/{transmission,vaccine}/
├── 03_bahia_pipeline/            <- model_fits/ diagnostics/ results/ tables/ figures/ reports/
├── 04_pernambuco_pipeline/       <- (same shape)
├── 05_rio_de_janeiro_pipeline/   <- (same shape)
├── 06_mato_grosso_pipeline/      <- (same shape)
├── 07_national_pipeline/         <- model_fits/ tables/ figures/ (no per-state model_fits split)
└── 90_development_archive/       <- frozen historical outputs, mirrors 02_Script/90_development_archive/
```

`01_data_pipeline` has no folder here yet: all of its outputs (fetched and
cleaned data) are written to `01_Data/`, not `03_Output/`. A folder will
be created here if a data-pipeline-specific diagnostic or QA report is
ever added.

## Category meaning (per pipeline)

- **`model_fits/`** -- fitted Stan objects (`.rds` posterior draws) and
  their per-chain sampler-progress CSVs. The single largest category by
  disk size; almost entirely gitignored (regenerate by rerunning the
  corresponding fitting stage with its `*_ALLOW_REFIT` gate set).
- **`tables/`** -- `.csv`/tabular `.rds` results.
- **`figures/`** -- `.png`/`.pdf` plots.
- **`results/`** -- bundled analysis results that are more than a single
  table (e.g. a paired-draws counterfactual bundle).
- **`diagnostics/`** -- HMC/PPC/identifiability diagnostic output, kept
  separate from `results/` because it answers "is this fit trustworthy?"
  rather than "what does this fit imply?"
- **`reports/`** -- written `.md` assessment/design documents.

Only categories that a pipeline actually produces exist under it -- empty
category folders are not created purely for symmetry.

## The one exception: processed model inputs

A `.rds` that is a **processed input to a Stan fit** (a prepared data
list, a period/config definition) rather than a fitted result or a raw
data file belongs under `01_Data/processed/model_inputs/`, not here --
see `REPOSITORY_ARCHITECTURE.md` for the full placement rule.

## 90_development_archive

Outputs of the frozen model-development history recorded in
`02_Script/90_development_archive/`. Structured by output type
(`figures/`, `tables/`, `model_fits/`) rather than by pipeline number,
since archived versions predate the current pipeline numbering. Not read
by any current production script.
