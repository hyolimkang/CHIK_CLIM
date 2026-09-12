# Verifies renewal_ceara_v4_1_q_only.stan is a clean single-factor ablation
# of the ARCHIVED, checksum-verified v4.0/td14 source: identical line-for-
# line except for (a) q_fixed removed from `data`, (b) `q` added to
# `parameters` with its prior in `model`, (c) `expected_reported_cases[t] =
# q_fixed * X[t]` -> `q * X[t]`. No other line may differ -- in particular
# the annual-effect block (`year_effect = year_effect_prior_sd * (z_year -
# mean(z_year))`) must be untouched, with NO sigma_year and NO Y-1 basis.

verify_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_verify_ablation_identity <- function() {
  root <- verify_root()
  archived_path <- file.path(root, "02_Script", "40_renewal_model", "20_v4_1_q_only", "archive", "renewal_ceara_v4_0_minimal_no_vaccine_ARCHIVED.stan")
  new_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_1_q_only.stan")
  if (!file.exists(archived_path)) stop("Archived v4.0 source missing: ", archived_path)
  if (!file.exists(new_path)) stop("v4.1 q-only source missing: ", new_path)

  archived_sha <- tools::md5sum(archived_path) # quick local integrity check (SHA-256 recorded in CHECKSUM_MANIFEST.md)
  message("[verify] archived v4.0 source md5: ", archived_sha)

  strip_comment_header <- function(lines) {
    first_code_line <- which(!startsWith(trimws(lines), "//") & nzchar(trimws(lines)))[1]
    lines[first_code_line:length(lines)]
  }
  v4_0_body <- strip_comment_header(readLines(archived_path))
  v4_1_body <- strip_comment_header(readLines(new_path))

  # Expected, and ONLY expected, differences:
  expected_removed <- 'real<lower=0, upper=1> q_fixed;'
  expected_added_param <- 'real<lower=0, upper=1> q;'
  expected_prior_line <- 'q ~ beta(16.2, 108.7);'
  expected_old_assignment <- 'expected_reported_cases[t] = q_fixed * X[t] + 1e-9;'
  expected_new_assignment <- 'expected_reported_cases[t] = q * X[t] + 1e-9;'

  v4_0_norm <- trimws(v4_0_body)
  v4_1_norm <- trimws(v4_1_body)

  # Remove the one data-line difference and align.
  v4_0_no_qfixed <- v4_0_norm[v4_0_norm != expected_removed]
  v4_1_no_qparam <- v4_1_norm[!(v4_1_norm %in% c(expected_added_param, expected_prior_line))]
  v4_1_no_qparam <- sub(expected_new_assignment, expected_old_assignment, v4_1_no_qparam, fixed = TRUE)

  if (length(v4_0_no_qfixed) != length(v4_1_no_qparam)) {
    stop("Line count mismatch after removing the expected q-related lines -- an unexpected structural change is present. v4.0 (minus q_fixed): ",
         length(v4_0_no_qfixed), " lines; v4.1 (minus q/prior, assignment reverted): ", length(v4_1_no_qparam), " lines.")
  }
  diffs <- which(v4_0_no_qfixed != v4_1_no_qparam)
  if (length(diffs)) {
    message("[verify] FAIL: unexpected differences beyond the q ablation:")
    for (i in diffs) message("  v4.0: ", v4_0_no_qfixed[i], "\n  v4.1: ", v4_1_no_qparam[i])
    stop("v4.1 q-only is NOT a clean single-factor ablation of the archived v4.0/td14 source.")
  }

  # Explicit negative checks: confirm sigma_year / B_year / Y-1 basis are NOT present.
  forbidden_patterns <- c("sigma_year", "B_year", "z_year_free")
  found_forbidden <- forbidden_patterns[vapply(forbidden_patterns, function(p) any(grepl(p, v4_1_body, fixed = TRUE)), logical(1))]
  if (length(found_forbidden)) stop("Forbidden pattern(s) found in v4.1 q-only Stan source (this must be a q-ONLY ablation): ", paste(found_forbidden, collapse = ", "))

  message("[verify] PASS: renewal_ceara_v4_1_q_only.stan differs from the archived, checksum-verified v4.0/td14 source ONLY by:")
  message("  (1) `", expected_removed, "` (data) removed")
  message("  (2) `", expected_added_param, "` (parameter) added, with prior `", expected_prior_line, "`")
  message("  (3) `q_fixed * X[t]` -> `q * X[t]` in expected_reported_cases")
  message("[verify] Confirmed absent: sigma_year, B_year, z_year_free (no simultaneous annual-effect changes).")
  invisible(TRUE)
}

if (sys.nframe() == 0L) run_verify_ablation_identity()
