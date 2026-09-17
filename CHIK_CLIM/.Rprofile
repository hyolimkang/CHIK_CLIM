# Auto-source project initialisation (Step 00) when CHIK_CLIM.Rproj opens.
# All actual setup logic lives in 00_project_setup.R -- keep this file tiny.

setup_file <- file.path(getwd(), "00_project_setup.R")

if (file.exists(setup_file)) {
  source(setup_file, local = .GlobalEnv)
}

rm(setup_file)
