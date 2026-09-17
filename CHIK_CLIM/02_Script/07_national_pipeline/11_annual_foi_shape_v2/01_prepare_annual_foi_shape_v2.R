# Prepare a separate input/audit bundle for the AR(1)-free annual FOI revision.
project_root_v2 <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  file.path(root, "CHIK_CLIM")
}
source(file.path(project_root_v2(), "02_Script", "90_development_archive", "06_legacy_susceptibility", "prepare_dynamic_annual_foi_data.R"))
prepared_v2 <- build_dynamic_annual_foi_data(write_outputs = FALSE)
out_dir <- project_path("03_Output", "07_national_pipeline", "tables", "annual_foi_shape_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(prepared_v2, file.path(out_dir, "annual_foi_shape_v2_input.rds"))
readr::write_csv(prepared_v2$annual, file.path(out_dir, "annual_foi_shape_v2_input_audit.csv"))
readr::write_csv(prepared_v2$q_scenarios, file.path(out_dir, "annual_foi_shape_v2_q_prior_audit.csv"))
readr::write_csv(tibble::tibble(state = "Ceara", foi_definition = "foi_equiv", valid_population_coverage = prepared_v2$foi_anchor$population_coverage), file.path(out_dir, "annual_foi_shape_v2_foi_coverage_audit.csv"))
message("[prepare-v2] saved separate annual FOI shape input bundle")
