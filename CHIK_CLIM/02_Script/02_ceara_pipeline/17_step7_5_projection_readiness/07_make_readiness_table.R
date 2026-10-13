# PHASE 5/6 -- authored projection-readiness decision table (synthesis
# across Step 6/7 and this module's Phase 0-4 findings; not a new
# simulation output). Saved as its own CSV per the task's Phase 6 deliverable.

suppressPackageStartupMessages({ library(here); library(readr); library(tibble) })
source(file.path(ROOT, "02_Script/02_ceara_pipeline/17_step7_5_projection_readiness/00_config_step7_5.R"))

readiness <- tribble(
  ~objective, ~supported_conclusions, ~unsupported_conclusions, ~dependence_on_unvalidated_assumptions, ~evidence_required_before_stronger_claims,

  "A. Absolute incidence forecasting (2026-2045 baseline)",
  "The generator produces internally consistent, numerically correct weekly dynamics (R0/Reff/demography/observation all pass consistency checks); susceptibility and peak transmission potential are well-estimated and stable across independent holdouts (H1 vs. H2 nearly identical).",
  "That the predictive distribution's MEDIAN or central tendency is a reliable point forecast of absolute case counts: H1 underpredicts the one major historical resurgence by 369x at median even under the most favourable tested importation assumption (lambda=4). Coverage of extreme events (2022-scale) is not achieved at any tested constant-rate lambda.",
  "Entirely dependent on the as-yet-unvalidated introduction-process structure/magnitude (Step 7.5 Phase 1-3) and on the still-provisional tau_year HMC reliability (upper-tail width).",
  "A structurally revised, jointly-refit introduction process (Phase 2's finding: equal-mean episodic bursts at tested magnitudes do NOT close the gap -- larger, rarer events or genomic-anchored event-intensity data are needed) PLUS resolved HMC reliability before ANY absolute-incidence claim is defensible.",

  "B. Conditional climate-scenario analysis (reference vs. SSP, holding epidemic assumptions fixed and labelled)",
  "Climate enters the model through a mechanistically verified, numerically consistent log-R0 pathway (R0 component reconstruction matches to 4.4e-16); RELATIVE differences between two climate inputs under the SAME epidemic-generator assumptions are computed from the same internally consistent machinery throughout.",
  "That a reference-vs-SSP CONTRAST is free of the introduction-process limitation: if the generator structurally underestimates standing infection/persistence, a climate-driven amplification of an already-too-low baseline may itself be biased in ways not characterised here. Absolute magnitudes under any climate scenario inherit all of (A)'s limitations.",
  "Depends on treating the epidemic-generator assumptions (tau_year, constant-rate lambda=1) as an explicitly LABELLED, fixed conditioning choice, not a validated prediction -- and on the climate-standardisation leakage (Step 6/7 Phase 0 item 4, small but uncorrected) being immaterial to the RELATIVE contrast specifically (not separately tested here).",
  "A sensitivity check of the reference-vs-SSP contrast ACROSS the tau_year prior range and the importation structural alternatives (Step 7's lambda grid, this module's episodic arms) to confirm the relative climate signal is not an artefact of one specific, unvalidated epidemic-generator configuration.",

  "C. Vaccination strategy sensitivity (comparing strategies across a transparent range of epidemic assumptions)",
  "A TRANSPARENT range of plausible epidemic trajectories (e.g. the full Step 7/7.5 lambda and structural grid) can be used as alternative CONDITIONING scenarios; comparing vaccination strategies WITHIN each such scenario, never claiming one scenario is 'the' true future, is methodologically sound given everything else in this report.",
  "That relative vaccination effects are AUTOMATICALLY reliable simply because they are relative: if the underlying epidemic dynamics are structurally misspecified (e.g. wrong introduction-process magnitude), the TIMING and AGE/SPATIAL distribution of transmission -- which can matter materially for vaccination impact -- may also be misspecified in ways a purely aggregate case-count comparison would not reveal.",
  "Depends on vaccination-relevant quantities (e.g. age-stratified attack rates, timing of peak incidence relative to a campaign) not being disproportionately sensitive to the SAME introduction-process limitation documented here -- NOT tested in this module.",
  "Explicit testing of whether vaccination-strategy RANKINGS (not just aggregate case counts) are stable across the Step 7/7.5 structural grid, plus any age/spatial structure the eventual vaccination model requires, validated independently of the aggregate-incidence calibration problem documented here."
)
write_csv(readiness, file.path(DIR_TABLES_STEP75, "PHASE5_projection_readiness_decision.csv"))
message("Projection-readiness decision table saved -> ", file.path(DIR_TABLES_STEP75, "PHASE5_projection_readiness_decision.csv"))
print(as.data.frame(readiness[, c("objective")]), row.names = FALSE)
