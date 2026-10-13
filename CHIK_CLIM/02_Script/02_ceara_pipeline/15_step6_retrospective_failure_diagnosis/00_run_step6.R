# Run the complete additive Step 6 analysis. No Stan fitting is performed.
source(file.path(ROOT, "02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis/01_build_diagnostics.R"))
source(file.path(ROOT, "02_Script/02_ceara_pipeline/15_step6_retrospective_failure_diagnosis/02_make_figures_and_report.R"))
message("STEP 6 complete: ", file.path(DIR_REPORTS_STEP6, "STEP6_RETROSPECTIVE_FAILURE_DIAGNOSIS.md"))

