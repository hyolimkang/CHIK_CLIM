# Active scripts moved to 38_ceara_current_pipeline/

As part of the Céara production-pipeline refactor, every active script
that used to live in `37_ceara_historical_age12_vaccine_counterfactual/scripts/`
was moved (via `git mv`, history preserved) into:

    02_Script/40_renewal_model/38_ceara_current_pipeline/
    (mainly stages 06_vaccine_counterfactual, 07_vaccine_validation,
    08_impact_nnv_summaries, 09_publication_figures)

The old ambiguous `scripts/99_run_all.R` orchestrator was retired and
replaced by `38_ceara_current_pipeline/06_vaccine_counterfactual/00_run_stage_06.R`.

See `38_ceara_current_pipeline/REFACTOR_MIGRATION_MAP.csv` for the full
old -> new mapping, and `38_ceara_current_pipeline/00_README_EXECUTION_ORDER.md`
for the canonical execution order.

**Nothing else in this folder changed.** `results/`, `tables/`,
`diagnostics/`, and `reports/` (including
`HISTORICAL_AGE12_VACCINE_COUNTERFACTUAL_ASSESSMENT.md` and
`NNV_ASSESSMENT.md`) remain exactly where they were, and the moved scripts'
`01_config.R` (formerly `00_config.R`) still points `ANALYSIS_DIR` at this
same folder — this was a code-location-only refactor (no scientific
change, no output relocation).
