# Repository Architecture

Final structure after the post-refactor architecture cleanup (see
`REPOSITORY_FINAL_CLEANUP_MANIFEST.csv` for the file-by-file record).

## The three top-level roles

```
01_Data     = data / model inputs        (never hand-edit; regenerate via 02_Script)
02_Script   = human-written source code  (.R, .stan)
03_Output   = generated analytical artifacts (never hand-edit; regenerate via 02_Script)
```

## File-type placement rule

| Extension / content | Belongs under | Exception |
|---|---|---|
| `.R` | `02_Script/` | -- |
| `.stan` (active) | `02_Script/00_shared/stan/current/<subgroup>/` | legacy -> `02_Script/90_development_archive/08_old_stan_models/` |
| `.rds` fitted Stan object (posterior draws) | `03_Output/<pipeline>/model_fits/...` | legacy -> `03_Output/90_development_archive/model_fits/<version>/` |
| `.rds` processed model **input** (a prepared Stan data list, a period/config definition -- not a fit, not raw data) | `01_Data/processed/model_inputs/` | -- |
| `.csv` / `.rds` table | `03_Output/<pipeline>/tables/` | -- |
| `.png` / `.pdf` figure | `03_Output/<pipeline>/figures/` | -- |
| `.md` report/assessment | `03_Output/<pipeline>/reports/` | -- |
| Raw/fetched data (`.rds`, `.csv`, cache directories) | `01_Data/` | -- |

The one file type that is genuinely ambiguous by extension alone is
`.rds`: it can be a fitted model, a processed input, or a raw-data cache,
and this repository has all three. It was never classified by extension
during this cleanup -- every `.rds` under a Stan source directory was
individually checked against its producing/consuming script (and, where
that was ambiguous, its file size relative to the ~1.7-1.9MB signature of
an `rstan_options(auto_write = TRUE)` compiled-model cache -- see
`02_Script/00_shared/stan/MODEL_REGISTRY.md`'s closing note).

## 02_Script <-> 03_Output correspondence

```
02_Script/02_ceara_pipeline            <->  03_Output/02_ceara_pipeline
02_Script/03_bahia_pipeline            <->  03_Output/03_bahia_pipeline
02_Script/04_pernambuco_pipeline       <->  03_Output/04_pernambuco_pipeline
02_Script/05_rio_de_janeiro_pipeline   <->  03_Output/05_rio_de_janeiro_pipeline
02_Script/06_mato_grosso_pipeline      <->  03_Output/06_mato_grosso_pipeline
02_Script/07_national_pipeline         <->  03_Output/07_national_pipeline
02_Script/90_development_archive       <->  03_Output/90_development_archive
```

`02_Script/00_shared` and `02_Script/01_data_pipeline` are the two
exceptions: shared code has no output bucket of its own (its live
model-fitting sub-projects' outputs attribute to whichever pipeline
originated them -- e.g. `legacy_model_functions/`'s Ceara-origin fits live
under `03_Output/02_ceara_pipeline/model_fits/`), and `01_data_pipeline`'s
outputs go to `01_Data/`, not `03_Output/`.

## Verification performed

- Every current pipeline script's `source()`, Stan-source, `readRDS()`,
  and write-path (`saveRDS()`/`write_csv()`/`ggsave()`) target was checked
  to exist (for reads) or have a valid parent directory (for writes)
  after every relocation in this cleanup.
- No `.rds` remains under any Stan source directory
  (`02_Script/00_shared/stan/`, and `02_Script/stan/` no longer exists).
- No unexplained loose analytical artifact remains directly under the old
  `03_Output/{figures,tables,results,diagnostics,reports}` roots -- those
  directories no longer exist; every file that lived there was routed to
  a specific pipeline's or the archive's output tree.
- No current pipeline was refit; no Stan model equation, prior, seed
  rule, q value, climate definition, demographic bookkeeping, or
  vaccination/age mechanic was changed. This was a code-location and
  output-location change only.
