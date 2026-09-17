# Standalone numerical check for the Y x (Y-1) orthonormal sum-to-zero
# annual-effect basis used by v4.1 (Part 4 of the task spec): confirms
# t(B) %*% B = I and colSums(B) = 0 to floating-point tolerance, for the
# actual Y used in the Ceará 2015-2025 fit, and additionally checks that
# every posterior draw's sum(year_effect) is ~0 once a fit exists.

v4_1_check_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

source(file.path(
  v4_1_check_root(), "02_Script", "90_development_archive", "02_historical_renewal_versions", "18_v4_1_minimal_no_vaccine", "01_fit_v4_1_minimal_no_vaccine.R"
))

run_check_orthonormal_year_basis <- function(Y = 11L, fit_path = NULL) {
  B <- build_orthonormal_year_basis(Y)
  orthonormality_error <- max(abs(crossprod(B) - diag(Y - 1L)))
  sum_to_zero_error <- max(abs(colSums(B)))
  message(sprintf(
    "[year-basis check] Y=%d | max|t(B)B - I| = %.3e | max|colSums(B)| = %.3e",
    Y, orthonormality_error, sum_to_zero_error
  ))
  if (orthonormality_error > 1e-8 || sum_to_zero_error > 1e-8) {
    stop("Orthonormal sum-to-zero basis failed numerical verification.")
  }

  if (is.null(fit_path)) fit_path <- v4_1_paths()$fit
  if (file.exists(fit_path)) {
    bundle <- readRDS(fit_path)
    sum_year_effect <- rstan::extract(bundle$fit, pars = "sum_year_effect")$sum_year_effect
    message(sprintf(
      "[year-basis check] sum(year_effect) across %d posterior draws: max|.| = %.3e",
      length(sum_year_effect), max(abs(sum_year_effect))
    ))
    if (max(abs(sum_year_effect)) > 1e-6) {
      stop("Posterior draws of sum(year_effect) are not ~0 -- basis or fit is inconsistent.")
    }
  } else {
    message("[year-basis check] no fit found at ", fit_path, "; basis-only check complete.")
  }
  invisible(list(B = B, orthonormality_error = orthonormality_error, sum_to_zero_error = sum_to_zero_error))
}

if (sys.nframe() == 0L) run_check_orthonormal_year_basis()
