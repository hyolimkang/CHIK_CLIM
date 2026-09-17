# 00_shared

Code that more than one pipeline (`01_data_pipeline` through
`07_national_pipeline`) depends on.

```
00_shared/
├── functions/               <- reusable R helpers (state-weekly aggregation, NA handling, thinning)
├── legacy_model_functions/   <- still-live shared model-fitting code (see its own README.md)
├── stan/                     <- the one canonical active Stan source location (see its own README.md)
├── plotting/                 <- reserved for shared plotting helpers (currently empty)
└── validation/               <- reserved for shared validation helpers (currently empty)
```

## What belongs here

- Reusable R functions sourced by two or more pipelines (e.g.
  `build_state_weekly()` in `functions/01_build_ceara_state_weekly.R`,
  used by every state's input-preparation stage).
- Active Stan model source code that more than one state fits verbatim, or
  that a state's pipeline still depends on for shared fitting machinery
  (`legacy_model_functions/`, `stan/current/`).
- Shared plotting or validation helpers, once any are extracted from a
  single pipeline into common use (`plotting/`, `validation/` are
  placeholders for this, not yet populated).

## What does NOT belong here

- Fitted model objects (`.rds` posterior draws) -- these belong under
  `03_Output/<pipeline>/model_fits/...`.
- Generated tables, figures, reports, or diagnostics -- these belong under
  the corresponding `03_Output/<pipeline>/` category.
- Region-specific outputs of any kind, even if the *code* that produced
  them lives here (e.g. `legacy_model_functions/17_v4_0.../`'s fitted
  objects live in `03_Output/02_ceara_pipeline/model_fits/initialization/`,
  not next to the script).

`legacy_model_functions/` is a special case: it looks like frozen version
history (numbered the same as when it lived under the now-deleted
`40_renewal_model/`), but every folder in it is a **live** dependency --
every state's Stage 02 transmission-fitting script sources one or more of
them. See `legacy_model_functions/README.md` before moving or renaming
anything there.
