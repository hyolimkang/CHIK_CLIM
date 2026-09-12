# Verifies renewal_ceara_v4_2_q_juazeiro_serology.stan differs from the
# archived, checksum-verified v4.0/td14 source ONLY by: (a) the same
# q_fixed -> q change already verified for the q-only ablation, and (b) a
# strictly ADDITIVE serology block (new data fields, new transformed-
# parameter lines computing p_state_sero_at_anchor/alpha_sero/beta_sero, one
# new model-block sampling statement, two new generated-quantities lines).
# No existing v4.0 line (the recursion, priors, NB likelihood) may be
# altered.

verify_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_verify_ablation_identity <- function() {
  root <- verify_root()
  archived_path <- file.path(root, "02_Script", "40_renewal_model", "21_v4_2_q_juazeiro_serology", "archive", "renewal_ceara_v4_0_minimal_no_vaccine_ARCHIVED.stan")
  new_path <- file.path(root, "02_Script", "stan", "renewal_ceara_v4_2_q_juazeiro_serology.stan")
  if (!file.exists(archived_path)) stop("Archived v4.0 source missing: ", archived_path)
  if (!file.exists(new_path)) stop("v4.2 source missing: ", new_path)

  strip_all_comments_and_blanks <- function(lines) {
    trimmed <- trimws(lines)
    trimmed[!startsWith(trimmed, "//") & nzchar(trimmed)]
  }
  v4_0_body <- strip_all_comments_and_blanks(readLines(archived_path))
  v4_2_body <- strip_all_comments_and_blanks(readLines(new_path))

  # Lines that are pure ADDITIONS in v4.2 (not present, in any form, in
  # v4.0) -- the serology block plus q's own new lines.
  addition_patterns <- c(
    "^int<lower=1, upper=N> t_sero;", "^int<lower=0> sero_pos;", "^int<lower=1> sero_n;", "^real<lower=0> kappa_sero;",
    "^real<lower=0, upper=1> q;$", "^real p_state_sero_at_anchor;", "^real p_sero_safe;", "^real alpha_sero;", "^real beta_sero;",
    "^p_state_sero_at_anchor = immune_prop\\[t_sero\\];", "^p_sero_safe = fmin", "^alpha_sero = p_sero_safe", "^beta_sero = \\(1 - p_sero_safe\\)",
    "^q ~ beta\\(16.2, 108.7\\);", "^sero_pos ~ beta_binomial", "^real log_lik_serology;", "^int sero_pred;",
    "^log_lik_serology = beta_binomial_lpmf", "^sero_pred = beta_binomial_rng",
    "^real log_prior_q;", "^log_prior_q = beta_lpdf",
    "^//", "^$" # comment/blank lines (already stripped, but guard)
  )
  is_addition <- function(line) any(vapply(addition_patterns, function(p) grepl(p, line), logical(1)))

  v4_2_non_addition <- v4_2_body[!vapply(v4_2_body, is_addition, logical(1))]
  # Revert the one expected substitution (q_fixed <-> q) for direct comparison.
  v4_2_reverted <- sub("^expected_reported_cases\\[t\\] = q \\* X\\[t\\] \\+ 1e-9;$", "expected_reported_cases[t] = q_fixed * X[t] + 1e-9;", v4_2_non_addition)
  v4_0_no_qfixed_decl <- v4_0_body[v4_0_body != "real<lower=0, upper=1> q_fixed;"]

  if (length(v4_0_no_qfixed_decl) != length(v4_2_reverted)) {
    stop("Line count mismatch after removing additions -- unexpected structural change. v4.0(-q_fixed decl): ",
         length(v4_0_no_qfixed_decl), " lines; v4.2(-additions, q reverted): ", length(v4_2_reverted), " lines.")
  }
  diffs <- which(v4_0_no_qfixed_decl != v4_2_reverted)
  if (length(diffs)) {
    message("[verify] FAIL: unexpected differences beyond the documented q + serology additions:")
    for (i in diffs) message("  v4.0: ", v4_0_no_qfixed_decl[i], "\n  v4.2: ", v4_2_reverted[i])
    stop("v4.2 is NOT a clean additive extension of the archived v4.0/td14 source.")
  }

  forbidden_patterns <- c("sigma_year", "B_year", "z_year_free", "quixad", "Quixad")
  found_forbidden <- forbidden_patterns[vapply(forbidden_patterns, function(p) any(grepl(p, readLines(new_path), ignore.case = FALSE)), logical(1))]
  if (length(found_forbidden)) stop("Forbidden pattern(s) found: ", paste(found_forbidden, collapse = ", "))

  message("[verify] PASS: renewal_ceara_v4_2_q_juazeiro_serology.stan = archived v4.0/td14 backbone")
  message("  + q (parameter, Beta(16.2,108.7) prior) replacing q_fixed (data) -- identical to the q-only ablation's change")
  message("  + ONE additive serology block (t_sero/sero_pos/sero_n/kappa_sero as data; p_state_sero_at_anchor/alpha_sero/beta_sero computed; beta_binomial likelihood; log_lik_serology/sero_pred generated quantities)")
  message("  Confirmed absent: sigma_year, B_year, z_year_free, any Quixada reference.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) run_verify_ablation_identity()
