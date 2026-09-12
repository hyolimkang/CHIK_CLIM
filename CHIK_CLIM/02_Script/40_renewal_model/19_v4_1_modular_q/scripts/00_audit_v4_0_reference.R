# First Action 1-2 (task spec): locate and verify the exact HMC-passing
# v4.0 td14 source, and verify the new v4.1 modular-q Stan model is
# mathematically identical to it conditional on q (i.e. differs ONLY in
# renaming q_fixed -> q_external; no other line changes).
#
# v4.0 td14 is the reference because it is the run that PASSED the strict
# HMC gate (divergences=0, treedepth hits=0, max Rhat=1.00897, min bulk
# ESS=480.8, BFMI all >0.30) using
# 17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R with
# RENEWAL_V4_RUN_TAG=td14, RENEWAL_V4_MAX_TREEDEPTH=14 (all other settings
# default). This script does not refit anything -- it only verifies
# provenance and byte-level structural equivalence before any node is run.

audit_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_audit_v4_0_reference <- function() {
  root <- audit_root()
  v4_0_stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_0_minimal_no_vaccine.stan")
  v4_0_fit_script <- file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R")
  v4_0_td14_fit_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_0_minimal_no_vaccine_fit_td14.rds")
  v4_1_stan_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_1_modular_q.stan")

  for (p in c(v4_0_stan_path, v4_0_fit_script, v4_0_td14_fit_path)) {
    if (!file.exists(p)) stop("Required v4.0 reference file is missing: ", p)
  }
  message("[audit] v4.0 Stan source:      ", v4_0_stan_path)
  message("[audit] v4.0 fit script:       ", v4_0_fit_script)
  message("[audit] v4.0 td14 fit bundle:  ", v4_0_td14_fit_path)

  td14_bundle <- readRDS(v4_0_td14_fit_path)
  hmc <- td14_bundle$hmc
  message(sprintf(
    "[audit] td14 stored HMC gate: pass=%s | divergences=%d | max_treedepth_hits=%d | max_rhat=%.5f | min_bulk_ess=%.1f | bfmi_range=[%.3f, %.3f]",
    hmc$hmc_pass, hmc$divergences, hmc$max_treedepth_hits, hmc$maximum_rhat, hmc$minimum_bulk_ess, min(hmc$bfmi), max(hmc$bfmi)
  ))
  if (!isTRUE(hmc$hmc_pass)) stop("td14 bundle does not record hmc_pass=TRUE -- wrong reference file.")
  if (!identical(td14_bundle$config$max_treedepth, 14L)) stop("td14 bundle max_treedepth is not 14 -- wrong reference file.")
  message("[audit] confirmed: this IS the HMC-passing max_treedepth=14 reference run (td14).")
  message("[audit] td14 fixed inputs: q_fixed=", td14_bundle$config$q_fixed,
          ", imports_per_week=", td14_bundle$config$imports_per_week,
          ", year_effect_prior_sd=", td14_bundle$config$year_effect_prior_sd)

  if (!file.exists(v4_1_stan_path)) {
    message("[audit] v4.1 modular-q Stan file not created yet -- run this again after 00b creates it, or proceed to create it now.")
    return(invisible(list(v4_0_lines = readLines(v4_0_stan_path), v4_1_exists = FALSE)))
  }

  v4_0_lines <- readLines(v4_0_stan_path)
  v4_1_lines <- readLines(v4_1_stan_path)
  # Normalise ONLY the deliberate rename (q_fixed -> q_external) and the
  # top-of-file comment block, then require byte-for-byte equality of
  # everything else (data/transformed data/parameters/transformed
  # parameters/model/generated quantities).
  strip_comment_header <- function(lines) {
    first_code_line <- which(!startsWith(trimws(lines), "//") & nzchar(trimws(lines)))[1]
    lines[first_code_line:length(lines)]
  }
  v4_0_body <- strip_comment_header(v4_0_lines)
  v4_1_body <- strip_comment_header(v4_1_lines)
  v4_1_body_renamed_back <- gsub("q_external", "q_fixed", v4_1_body, fixed = TRUE)

  if (length(v4_0_body) != length(v4_1_body_renamed_back)) {
    stop("v4.1 modular-q Stan body has a different number of lines than v4.0 -- structural change detected, this violates 'preserve proven geometry'.")
  }
  diffs <- which(v4_0_body != v4_1_body_renamed_back)
  if (length(diffs)) {
    message("[audit] FAIL: ", length(diffs), " line(s) differ beyond the q_fixed -> q_external rename:")
    for (i in diffs) {
      message("  line ", i, ":\n    v4.0: ", v4_0_body[i], "\n    v4.1: ", v4_1_body_renamed_back[i])
    }
    stop("v4.1 modular-q Stan model is NOT structurally identical to the HMC-passing v4.0 td14 model.")
  }
  message("[audit] PASS: v4.1 modular-q Stan model is byte-for-byte identical to v4.0/td14, aside from the deliberate q_fixed -> q_external rename.")
  invisible(list(v4_0_lines = v4_0_lines, v4_1_lines = v4_1_lines, identical_conditional_on_q = TRUE))
}

if (sys.nframe() == 0L) run_audit_v4_0_reference()
