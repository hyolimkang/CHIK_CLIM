# ============================================================
# STEP 00 -- CHIK_CLIM PROJECT INITIALISATION
# ============================================================
#
# Run automatically when CHIK_CLIM.Rproj opens via .Rprofile.
#
# This script:
#   - identifies project root
#   - loads common libraries
#   - configures RStan
#   - defines project paths
#
# This script DOES NOT:
#   - prepare data
#   - fit models
#   - run simulations
#   - generate results
#
# All scientific workflows begin AFTER Step 00.
#
# Safe to source more than once in the same session (idempotent): no
# rm(list = ls()), no data loading, no package installation, no model
# fitting. For a clean session use RStudio's Session > Restart R.

# ----------------------------------------------------------------
# 1. Locate project root (no hard-coded, user-specific path)
# ----------------------------------------------------------------

#' Find the CHIK_CLIM project root by searching for CHIK_CLIM.Rproj.
#'
#' Starts at `start` (default: current working directory) and:
#'   1. checks whether start/CHIK_CLIM/CHIK_CLIM.Rproj exists (the case
#'      where start is the outer git repository and the RStudio project
#'      lives one level down, in ./CHIK_CLIM);
#'   2. otherwise walks upward from start, at each level checking both
#'      <dir>/CHIK_CLIM.Rproj (dir IS the project root) and
#'      <dir>/CHIK_CLIM/CHIK_CLIM.Rproj (dir contains the project root).
#' Returns the normalised absolute project root path, or stops with a
#' clear error if CHIK_CLIM.Rproj cannot be found within `max_up` levels.
find_chik_clim_root <- function(start = getwd(), max_up = 10L) {
  outer_candidate <- file.path(start, "CHIK_CLIM", "CHIK_CLIM.Rproj")
  if (file.exists(outer_candidate)) {
    return(normalizePath(dirname(outer_candidate), winslash = "/", mustWork = TRUE))
  }

  dir <- normalizePath(start, winslash = "/", mustWork = TRUE)
  for (i in seq_len(max_up)) {
    direct_candidate <- file.path(dir, "CHIK_CLIM.Rproj")
    if (file.exists(direct_candidate)) {
      return(normalizePath(dir, winslash = "/", mustWork = TRUE))
    }
    nested_candidate <- file.path(dir, "CHIK_CLIM", "CHIK_CLIM.Rproj")
    if (file.exists(nested_candidate)) {
      return(normalizePath(dirname(nested_candidate), winslash = "/", mustWork = TRUE))
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) break # reached filesystem root
    dir <- parent
  }

  stop(
    "find_chik_clim_root(): could not locate CHIK_CLIM.Rproj by searching ",
    "upward from:\n  ", start, "\n",
    "Expected either <dir>/CHIK_CLIM.Rproj or <dir>/CHIK_CLIM/CHIK_CLIM.Rproj ",
    "at some level above (or equal to) the starting directory.",
    call. = FALSE
  )
}

ROOT <- find_chik_clim_root()
setwd(ROOT)

# ----------------------------------------------------------------
# 2. Common package set (check availability; never silently install)
# ----------------------------------------------------------------

PROJECT_PACKAGES <- c(
  "here", "rstan", "dplyr", "readr", "tibble", "tidyr",
  "lubridate", "ggplot2", "posterior", "patchwork", "scales"
)

missing_packages <- PROJECT_PACKAGES[
  !vapply(PROJECT_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "CHIK_CLIM project setup: missing required package(s):\n",
    paste0("  ", missing_packages, collapse = "\n"), "\n\n",
    "Install with:\n",
    "  install.packages(c(", paste0('"', missing_packages, '"', collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(here)
  library(rstan)
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(posterior)
  library(patchwork)
  library(scales)
})

here::i_am("CHIK_CLIM.Rproj")

# ----------------------------------------------------------------
# 3. Common RStan configuration (configuration only -- no compiling,
#    no fitting)
# ----------------------------------------------------------------

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores()))

# ----------------------------------------------------------------
# 4. Common path helpers
# ----------------------------------------------------------------

project_path <- function(...) file.path(ROOT, ...)

DIR_DATA    <- project_path("01_Data")
DIR_SCRIPT  <- project_path("02_Script")
DIR_OUTPUT  <- project_path("03_Output")
DIR_STAN    <- project_path("02_Script", "00_shared", "stan", "current")

# Per-pipeline output roots, mirroring 02_Script's numbered pipeline folders.
# `output_path("02_ceara_pipeline", "model_fits", "baseline")` is equivalent
# to `project_path("03_Output", "02_ceara_pipeline", "model_fits", "baseline")`.
output_path <- function(pipeline, ...) file.path(DIR_OUTPUT, pipeline, ...)

OUT_CEARA        <- output_path("02_ceara_pipeline")
OUT_BAHIA        <- output_path("03_bahia_pipeline")
OUT_PERNAMBUCO   <- output_path("04_pernambuco_pipeline")
OUT_RIO          <- output_path("05_rio_de_janeiro_pipeline")
OUT_MATO_GROSSO  <- output_path("06_mato_grosso_pipeline")
OUT_NATIONAL     <- output_path("07_national_pipeline")

# ----------------------------------------------------------------
# 5. Lightweight Ceara core-input check (reports only; never halts startup)
# ----------------------------------------------------------------

#' Check whether the key upstream files used by the accepted Ceara
#' pipeline currently exist on disk. Does NOT stop project startup if any
#' are missing -- only reports. Call explicitly when you need to know
#' whether downstream analysis inputs are available.
check_ceara_core_inputs <- function() {
  files <- c(
    "01_Data/chik_dlnm_panel_muni_week_2015_2025.rds",
    "01_Data/ibge_population_projection_uf_2024revision.rds",
    "01_Data/ceara_weekly_demography_2015_2025.rds",
    "01_Data/ibge_pop_uf_single_age_expanded_2015_2024.rds",
    "03_Output/07_national_pipeline/tables/national_wave_analysis/brazil_chik_uf_weekly_climate.csv"
  )
  tibble::tibble(file = files, exists = file.exists(project_path(files)))
}

# ----------------------------------------------------------------
# 6. Detailed version info, available on request only (not printed
#    automatically)
# ----------------------------------------------------------------

project_session_info <- function() {
  cat("R version:      ", R.version.string, "\n")
  cat("RStan version:  ", as.character(utils::packageVersion("rstan")), "\n")
  cat("rstan_options(\"auto_write\"): ", rstan::rstan_options("auto_write"), "\n")
  cat("mc.cores:       ", getOption("mc.cores"), "\n")
  cat("Project root:   ", ROOT, "\n\n")
  print(sessionInfo())
}

# ----------------------------------------------------------------
# 7. Concise startup message (no automatic sessionInfo() dump)
# ----------------------------------------------------------------

cat(strrep("-", 52), "\n")
cat("CHIK_CLIM PROJECT INITIALISED\n")
cat("Root:     ", ROOT, "\n")
cat("R:        ", paste(R.version$major, R.version$minor, sep = "."), "\n")
cat("RStan:    ", as.character(utils::packageVersion("rstan")), "\n")
cat("Cores:    ", getOption("mc.cores"), "\n")
cat("Packages: OK (", length(PROJECT_PACKAGES), " loaded)\n", sep = "")
cat(strrep("-", 52), "\n")
