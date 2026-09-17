# ============================================================
# 10_grid_susceptibility_reconstruction (National pipeline, orchestrator only)
# ============================================================
# Bottom-up (grid-cell) reconstruction of Brazil-wide susceptibility from
# the existing UF-level FOI ensemble. Existence-gated; only rebuilds when
# missing AND NATIONAL_ALLOW_REFIT="TRUE".

if (!exists("ROOT")) {
  setup_path <- normalizePath(file.path(getwd(), "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) setup_path <- normalizePath(file.path(getwd(), "..", "..", "00_project_setup.R"), mustWork = FALSE)
  if (!file.exists(setup_path)) stop("Cannot find 00_project_setup.R -- open CHIK_CLIM.Rproj first, or source it manually before this runner.")
  source(setup_path, local = .GlobalEnv)
}
FOLDER_DIR <- project_path("02_Script/07_national_pipeline/10_grid_susceptibility_reconstruction")
ALLOW_REFIT <- identical(Sys.getenv("NATIONAL_ALLOW_REFIT", "FALSE"), "TRUE")
Sys.setenv(GRID_SUSC_OUTPUT_TAG = "allfoi_s1_spatial_scale_v2") # matches the actual saved output tag; both scripts' own internal default ("") does not correspond to where the accepted results live
TABLE_DIR <- project_path("03_Output/tables/national_grid_susceptibility/allfoi_s1_spatial_scale_v2")

message("\n[10_grid_susceptibility_reconstruction] 01/02: reconstruct_brazil_grid_susceptibility ...")
fp <- file.path(TABLE_DIR, "brazil_chik_uf_bottomup_susceptibility.csv")
if (file.exists(fp)) {
  message("  Found: ", fp, " -- reusing (not rebuilt).")
} else if (ALLOW_REFIT) {
  source(file.path(FOLDER_DIR, "01_reconstruct_brazil_grid_susceptibility.R"), local = .GlobalEnv) # defines run_brazil_grid_susceptibility(); does not auto-run under source()
  run_brazil_grid_susceptibility()
} else {
  stop("Missing prerequisite:\n  ", fp, "\n\nSet Sys.setenv(NATIONAL_ALLOW_REFIT = \"TRUE\") to allow rebuilding, then re-source this runner.", call. = FALSE)
}

message("\n[10_grid_susceptibility_reconstruction] 02/02: plot_brazil_grid_susceptibility_maps ...")
source(file.path(FOLDER_DIR, "02_plot_brazil_grid_susceptibility_maps.R"), local = .GlobalEnv) # defines run_grid_susceptibility_maps(); does not auto-run under source()
run_grid_susceptibility_maps()

message("\n[10_grid_susceptibility_reconstruction] Complete.")
