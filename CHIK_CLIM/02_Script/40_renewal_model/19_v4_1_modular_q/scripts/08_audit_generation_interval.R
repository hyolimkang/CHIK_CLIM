# Generation-interval audit (task Section 15) -- documentation only. Does
# NOT change the generation interval in v4.1; w and G remain exactly as in
# v4.0/td14 (verified byte-identical by 00_audit_v4_0_reference.R).

required_packages <- c("here", "tibble", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(tibble); library(readr) })

gi_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

# Exact copy of the function in
# 17_v4_0_minimal_no_vaccine/01_fit_v4_0_minimal_no_vaccine.R (source file
# of record for G and w) -- reproduced here read-only for documentation;
# NOT a redefinition used by any fit.
generation_weights_source_copy <- function() {
  weights <- diff(pgamma(0:8, shape = 4, rate = 2))
  weights / sum(weights)
}

run_audit_generation_interval <- function() {
  root <- gi_root()
  out_dir <- file.path(root, "02_Script", "40_renewal_model", "19_v4_1_modular_q")
  source_file <- file.path(root, "02_Script", "40_renewal_model", "17_v4_0_minimal_no_vaccine", "01_fit_v4_0_minimal_no_vaccine.R")

  G <- 8L
  shape <- 4; rate <- 2 # Gamma(shape, rate), units: WEEKS (weights are applied to X[t-g] for weekly-indexed g)
  w_raw <- diff(pgamma(0:G, shape = shape, rate = rate))
  retained_probability_before_renorm <- pgamma(G, shape = shape, rate = rate)
  w <- w_raw / sum(w_raw)
  stopifnot(isTRUE(all.equal(w, generation_weights_source_copy())))

  mean_weeks <- shape / rate
  sd_weeks <- sqrt(shape) / rate

  summary_tbl <- tibble(
    field = c(
      "source_file", "distribution", "shape", "rate", "mean_weeks", "mean_days", "sd_weeks", "sd_days",
      "discretisation_method", "truncation_G_weeks", "total_probability_retained_before_renormalisation"
    ),
    value = c(
      source_file, "Gamma(shape, rate), discretised weekly", as.character(shape), as.character(rate),
      as.character(mean_weeks), as.character(mean_weeks * 7), as.character(sd_weeks), as.character(sd_weeks * 7),
      "w_g = [pgamma(g, shape, rate) - pgamma(g-1, shape, rate)] for g=1..G, then renormalised to sum to 1",
      as.character(G), as.character(retained_probability_before_renorm)
    )
  )
  weekly_weights_tbl <- tibble(g = seq_len(G), w_g_raw = w_raw, w_g_renormalised = w)

  write_csv(summary_tbl, file.path(out_dir, "generation_interval_audit_summary.csv"))
  write_csv(weekly_weights_tbl, file.path(out_dir, "generation_interval_audit_weights.csv"))

  message("[GI audit] source file: ", source_file)
  message(sprintf("[GI audit] Gamma(shape=%d, rate=%d) in WEEKS -> mean=%.2f weeks (%.1f days), sd=%.2f weeks (%.1f days)",
                   shape, rate, mean_weeks, mean_weeks * 7, sd_weeks, sd_weeks * 7))
  message(sprintf("[GI audit] truncated at G=%d weeks; total probability mass retained before renormalisation = %.5f (%.3f%% discarded)",
                   G, retained_probability_before_renorm, 100 * (1 - retained_probability_before_renorm)))
  message("[GI audit] weekly weights (renormalised):")
  print(weekly_weights_tbl, n = Inf)
  message("[GI audit] NOTE: this is an audit only -- the generation interval is UNCHANGED in v4.1 (verified byte-identical to v4.0/td14).")
  invisible(list(summary = summary_tbl, weights = weekly_weights_tbl))
}

if (sys.nframe() == 0L) run_audit_generation_interval()
