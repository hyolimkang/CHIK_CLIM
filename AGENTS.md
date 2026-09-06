# Project instructions

This is a scientific modelling repository for chikungunya, primarily using R and Stan. Prioritize scientific clarity, reproducibility, and preservation of existing work.

## Repository layout and current workflow

- The Git repository root contains this file, `docs/`, and the inner `CHIK_CLIM/` R project. The R project root is `CHIK_CLIM/`, containing `CHIK_CLIM.Rproj`, `01_Data/`, `02_Script/`, and `03_Output/`. Paths below are relative to that R project root unless stated otherwise.
- `02_Script/00_data_prep/`: SINAN download and cleaning, IBGE municipal and age-specific population, annual births, and municipality boundaries/centroids. SINAN cleaning produces individual records and municipality-week/month panels, currently named for 2015-2024.
- `02_Script/10_sir_transmission_model/`: state-week panel construction, SIR fitting, diagnostics, and age-structured simulations. The fitter sources `01_build_uf_weekly_panel.R`; diagnostics and simulations depend on objects retained in the same R session.
- `02_Script/20_vaccine_impact/`: additional prevaccination and vaccination simulations with session/helper dependencies. `03_R_simulator_counterfact.R` checks R reproduction of Stan dynamics; it does not itself perform vaccination simulations.
- `02_Script/30_climate_covariates_dlnm/`: ERA5-Land temperature and precipitation downloads, municipality-week climate aggregation, and joining climate, cases, and population into the analysis panel. `05_dlnm_model.R` currently only inspects that panel. `00_current_covariates_baseline.R` belongs to a separate spatial covariate workflow and is not a prerequisite for the ERA5 pipeline.
- `02_Script/40_renewal_model/`: `02_fit_renewal_v1.R` sources `01_build_ceara_state_weekly.R`, aggregates the combined panel to Ceara, and fits `02_Script/stan/renewal_ceara_v1.stan`; `03_plot_renewal_v1_fit.R` reads the saved fit. Version 1 uses constant R, lognormal process noise, negative-binomial observations, and no susceptible depletion or climate effect.
- `02_Script/90_exploratory/`: descriptive, recurrence, synchrony, FOI, burden, and climate analyses. Most are independent; the asynchrony supplement consumes the epidemic-wave audit table. Climate plots consume the combined climate/case panel.
- `02_Script/Functions/`: `fixNAs_adj.R` and `thinning_chik.R` support spatial covariate processing. Shared transmission and renewal helpers also live inside their respective pipeline folders.
- `02_Script/stan/`: current and historical model sources plus serialized R objects. The SIR `.rds` files inspected here are compiled `stanmodel` caches, whereas `renewal_ceara_v1.rds` is a posterior `stanfit`. Inspect object classes before using a cache as a fit.
- `03_Output/figures/` and `03_Output/tables/`: generated results, potentially from earlier code versions. Some simulations produce only in-memory objects; SIR plot export is controlled by `SAVE_PLOTS`.
- `02_Script/archive/`, top-level duplicate preparation scripts, `01_Data/current_covariates_baseline.R`, and `WorldClim.R` include historical or separate workflows. Presence or numbering alone does not establish execution order or active status.
- `analysis_plan.md` and the repository-root `docs/proposal_climate_renewal_immunisation.md` describe proposed work; do not assume their full nowcasting or climate-renewal-immunisation pipelines are implemented.

Recheck these observations against the code before editing. In the inspected SIR fitter, the configured state is Minas Gerais (`31`, `mg`), the actual fit call uses `ceara_weekly_v12f.stan`, the result is named `fit_v12e`, and `fit_v12c` is a compatibility alias. The call hardcodes v12f despite a separate model-selection variable. The v12d pulse fit is conditional. Do not infer model identity from filenames, object names, output prefixes, or version numbers alone.

## Repository safety

- Never modify raw data, including downloaded SINAN records and ERA5-Land source files.
- Do not move, rename, delete, or consolidate apparently duplicated or legacy scripts unless the user explicitly requests repository cleanup.
- Treat `archive/` files as historical reference; do not normally modify them.
- Existing outputs may come from earlier script versions. Verify provenance before treating them as evidence about current code.
- Large inputs and intermediate files are intentionally excluded from Git. Do not delete or recreate them unless explicitly requested; being ignored does not make them disposable.
- Untracked renewal scripts, Stan sources, fitted objects, and figures are user work, not disposable artifacts.
- Do not overwrite existing fitted results or outputs as an incidental consequence of a test. Use separate test destinations when possible.

## Scientific modelling

- Do not change scientific assumptions without explicitly identifying the change.
- Always distinguish mathematical/model formulation issues, implementation/coding issues, statistical identifiability issues, and limitations of the available data. Separate definite bugs from modelling choices and unresolved hypotheses.
- Prefer incremental changes and the simplest model that addresses the scientific question over redesigning the entire model.
- Flag potential confounding and identifiability problems, particularly reporting fraction versus infection scale, process versus observation noise, and time-varying transmission versus climate, seasonality, or susceptible depletion.
- State what each quantity represents: observed cases, latent incidence, infectious prevalence, reporting fraction, or reproduction number. Do not interchange these interpretations.
- Explicitly identify changes to generation intervals, time aggregation, seeding, priors, demography, reporting, or observation delays. Check whether multiplicative noise is centered on a mean or a median before interpreting R.
- The renewal input currently comes from the climate-complete panel. Check municipality exclusions when comparing its state totals with the SIR branch's case-panel aggregation.

## R

- Reuse existing project functions where possible, including the state-panel helpers, `check_stan_panel()`, and renewal `build_state_weekly()`.
- Avoid new packages unless necessary. There is no central dependency lockfile in the inspected project; inspect each script's package requirements.
- Prefer relative/project-aware paths, and verify that `here::here()` resolves to the inner R project rather than the outer Git root. Avoid new machine-specific absolute paths.
- Be cautious with external fallback paths and implicit session dependencies in older scripts. Identify required objects and helper definitions before running them.
- Several climate and renewal scripts guard execution with `if (sys.nframe() == 0)`. Ordinary `source()` skips those main blocks; use the appropriate execution mode. SIR downstream scripts instead expect a shared fitting session.
- Do not assume a script is active simply because it remains in `02_Script/`.
- Preserve downstream compatibility when changing columns, dimensions, object names, or serialized structures. Trace consumers before changing `ce_fit`, `fit_v12e`, `fit_v12c`, or `stan_data_v12e`.
- Inspect scripts before sourcing: some install packages, download data, or write files immediately.

## Stan

- Check indexing, dimensions, bounds, and parameter constraints against the R data list. For renewal models, verify complete weekly ordering, valid generation weights, and `G <= seed_weeks < N` for the current implementation.
- Explicitly inspect priors and parameter identifiability; check scale and distribution parameterization.
- Distinguish the latent process, observation likelihood, and generated quantities. Distinguish latent credible intervals from posterior predictive intervals, and observed-period replication from future forecasting.
- Do not silently alter the statistical model to solve a computational problem. Explain whether an intervention changes parameterization, numerical evaluation, or the scientific formulation.
- Compile every modified Stan model. When feasible, run a small sampling test after compilation.
- Inspect divergences, R-hat, ESS, treedepth, and obvious posterior pathologies. A short sampling test is a computational check, not evidence of adequate convergence or scientific validity.
- If compilation or sampling cannot be completed, report the exact blocker and what remains unverified.
- When the scientific model changes, create a new model version rather than overwriting an existing Stan model; keep previous sources and results available for comparison.

## Validation

For substantial model changes, consider comparison with the previous version, posterior sanity checks, posterior predictive checks, sensitivity to important assumptions, and out-of-sample or held-out validation when scientifically appropriate.

Do not automatically add complex validation. First propose it when it would materially affect the scientific conclusion. Match checks to the change, and distinguish successful execution from statistical adequacy. Never assume existing output files validate newly edited code.

## Workflow

Before substantial edits:

1. Inspect Git status and all relevant files, including callers, helpers, model sources, and downstream consumers.
2. Reconstruct the current workflow or model from executable code.
3. Explain what the current implementation does.
4. Identify definite bugs separately from modelling choices, identifiability issues, and data limitations.
5. Propose a minimal implementation plan.
6. Wait for the user's approval before making major scientific changes. Existing approval for the specific change remains valid; do not request it again unnecessarily.

When implementing an approved change:

- Modify only the requested component; avoid unrelated refactoring or cleanup.
- Preserve existing working behavior and downstream compatibility where possible.
- Create a new Stan version for changes to the scientific model, and make corresponding R model selection and result destinations explicit.

After editing, report:

- Files changed.
- Exact implementation changes.
- Scientific assumptions changed, if any.
- Commands/tests executed.
- Validation results, including checks not completed.
- Remaining concerns.
- Any unexpected downstream consequences.

## Git

- Inspect `git status` before substantial edits and preserve pre-existing modifications, deletions, and untracked files.
- Do not modify or stage unrelated files.
- Do not commit unless the user explicitly asks.
- Before a commit, summarize the Git diff and identify accidental or unrelated changes without reverting the user's work.
