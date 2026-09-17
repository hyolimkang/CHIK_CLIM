# ============================================================
# NATIONAL PIPELINE -- TOP-LEVEL ORCHESTRATOR
# ============================================================
# Unlike the state pipelines, this pipeline is not one sequential stage
# chain -- it bundles several thematically distinct national-scope
# analyses (see 00_README_EXECUTION_ORDER.md for what each folder is).
# Folder numbers (01-11) reflect topic grouping inherited from the source
# material, NOT execution order; EXECUTION_ORDER below is the actual
# dependency-respecting run order (census tables first, then the analyses
# that read them). ALL outputs for every folder already exist on disk as
# of this refactor, so a default run just prints "reusing" for each step.
#
# NATIONAL_ALLOW_REFIT gates every fit/rebuild step across all 11 folders
# (a single project-wide gate, not per-folder, since none of these tracks
# is as individually expensive/risky as a full state v4.9 Stan refit) --
# it is NEVER set automatically by this orchestrator; set it yourself
# first if you deliberately want to allow rebuilding something missing.

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this orchestrator.")
  source(setup_path, local = .GlobalEnv)
}

PIPELINE_DIR <- project_path("02_Script/07_national_pipeline")

EXECUTION_ORDER <- c(
  "04_national_episode_census",
  "05_national_wave_census_v2",
  "01_episode_renewal",
  "02_major_episode_renewal",
  "03_episode_susceptibility_reconstruction",
  "06_national_wave_analysis",
  "07_provisional_climate_susceptibility_checkpoint",
  "08_first_epidemic_climate_transmission",
  "09_recurrent_climate_susceptibility_phase",
  "10_grid_susceptibility_reconstruction",
  "11_annual_foi_shape_v2"
)

RUN_ONLY <- Sys.getenv("NATIONAL_RUN_ONLY", "")
folders_to_run <- if (nzchar(RUN_ONLY)) {
  requested <- trimws(strsplit(RUN_ONLY, ",")[[1]])
  missing <- setdiff(requested, EXECUTION_ORDER)
  if (length(missing) > 0) stop("Unknown folder(s) in NATIONAL_RUN_ONLY: ", paste(missing, collapse = ", "), call. = FALSE)
  requested
} else {
  EXECUTION_ORDER
}

message("[pipeline] Running National pipeline folders: ", paste(folders_to_run, collapse = ", "), "\n")
for (folder in folders_to_run) {
  runner_path <- file.path(PIPELINE_DIR, folder, "00_run_folder.R")
  message(strrep("=", 60)); message(sprintf("FOLDER: %s", folder)); message(strrep("=", 60))
  source(runner_path, local = .GlobalEnv)
}
message("\n[pipeline] Done. Folders run: ", paste(folders_to_run, collapse = ", "))
