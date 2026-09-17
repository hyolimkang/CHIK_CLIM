# Active scripts moved to 38_ceara_current_pipeline/

As part of the Céara production-pipeline refactor, every active script
that used to live in `36_climate_forced_v4_9/scripts/` was moved (via
`git mv`, history preserved) into:

    02_Script/40_renewal_model/38_ceara_current_pipeline/

See `38_ceara_current_pipeline/REFACTOR_MIGRATION_MAP.csv` for the full
old -> new mapping, and `38_ceara_current_pipeline/00_README_EXECUTION_ORDER.md`
for the canonical execution order.

**Nothing else in this folder changed.** `outputs/` (including the
accepted fit `outputs/ce_canary/ce_climate_forced_canary_q0.05.rds`) and
the phase reports (`PHASE0_AUDIT.md`, `CLIMATE_V49_CANARY_ASSESSMENT.md`,
etc.) remain exactly where they were — this was a code-location-only
refactor (no scientific change, no output relocation).
