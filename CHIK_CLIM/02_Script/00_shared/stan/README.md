# Shared Stan sources

This is the **one canonical location for active Stan source code** in the
project. `02_Script/stan/` (the old flat, mixed source+output directory) no
longer exists.

```
00_shared/stan/
├── README.md              <- this file
├── MODEL_REGISTRY.md       <- full file-by-file status table
└── current/
    ├── initialization/     <- minimal model used to initialise every state's fit
    ├── ceara/              <- Ceara-specific baseline / climate-forced / ablation models
    ├── multistate/         <- generic models shared verbatim by Bahia/PE/RJ/MT
    ├── national/           <- national-scope models (episode Re, annual FOI shape)
    ├── pernambuco/         <- Pernambuco-specific spatial/reparameterisation extensions
    └── rio_de_janeiro/     <- Rio de Janeiro city-level extension
```

**Only active Stan source code (`.stan`) lives here** -- see
`MODEL_REGISTRY.md` for exactly which file is used by which pipeline.
Compiled-model caches (`rstan_options(auto_write = TRUE)` saves a small
~1.7-1.9MB `.rds` next to its `.stan` source) travel alongside their source
file here so recompilation is never required by this reorganisation.

Fitted Stan objects (posterior draws) are **not** here -- per the
project's extension-based placement rule they live under
`03_Output/<pipeline>/model_fits/...`. Processed Stan *inputs* (data lists,
period definitions) that are not raw data live under
`01_Data/processed/model_inputs/`.

Legacy/superseded Stan sources (and their compiled-model caches) live under
`02_Script/90_development_archive/08_old_stan_models/`; their old fitted
results live under `03_Output/90_development_archive/model_fits/`. Nothing
here duplicates a generic model per state -- state-specific extensions only
exist where the model itself is genuinely different (see `multistate/` vs
`pernambuco/`/`rio_de_janeiro/` in the registry).

Do not move or rename any file here without first checking
`MODEL_REGISTRY.md`'s "Used by" column and re-running
`grep -rl "<filename>" 02_Script/` to confirm nothing else depends on it.
