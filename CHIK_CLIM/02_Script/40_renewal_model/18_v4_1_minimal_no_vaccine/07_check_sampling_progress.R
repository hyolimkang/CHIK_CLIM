# Live iteration-progress check for a running v4.1 fit.
#
# rstan parallelises chains on Windows via PSOCK worker processes, which do
# NOT stream Stan's refresh() console messages back to the parent R session
# -- so stdout redirection alone cannot show "Iteration: k / n" while chains
# are running. Instead, 01_fit_v4_1_minimal_no_vaccine.R now passes
# `sample_file`, which makes each chain write its own CSV of draws to disk
# incrementally (flushed row-by-row as Stan's C++ writer completes each
# iteration). This script just counts how many data rows exist in each
# chain's CSV right now -- safe to run at any time while the fit is in
# progress, including from a different R session.
#
# Usage: RENEWAL_V4_1_RUN_TAG=<tag> Rscript 07_check_sampling_progress.R

check_progress_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

count_csv_progress <- function(csv_path, warmup, iter) {
  if (!file.exists(csv_path)) return(list(exists = FALSE, warmup_done = NA, sampling_done = NA, phase = "not started"))
  lines <- readLines(csv_path, warn = FALSE)
  data_lines <- lines[!startsWith(lines, "#")]
  # First non-comment line is the CSV header; every line after that is one
  # written iteration (warmup rows come before Stan's own
  # "# Adaptation terminated" comment marker, sampling rows after).
  n_written <- max(0L, length(data_lines) - 1L)
  adaptation_line <- grep("Adaptation terminated", lines)
  if (length(adaptation_line)) {
    # Rows written before that comment line, minus the header, = warmup rows.
    lines_before <- lines[seq_len(adaptation_line[1] - 1L)]
    warmup_done <- max(0L, length(lines_before[!startsWith(lines_before, "#")]) - 1L)
    sampling_done <- max(0L, n_written - warmup_done)
    phase <- if (sampling_done >= (iter - warmup)) "complete" else "sampling"
  } else {
    warmup_done <- min(n_written, warmup)
    sampling_done <- 0L
    phase <- if (n_written > 0) "warmup" else "not started"
  }
  list(exists = TRUE, warmup_done = warmup_done, sampling_done = sampling_done, phase = phase, n_written = n_written)
}

run_check_sampling_progress <- function() {
  root <- check_progress_root()
  run_tag <- Sys.getenv("RENEWAL_V4_1_RUN_TAG", unset = "base")
  suffix <- if (identical(run_tag, "base")) "" else paste0("_", run_tag)
  progress_dir <- file.path(root, "02_Script", "stan", "sampling_progress", paste0("v4_1", suffix))
  iter <- as.integer(Sys.getenv("RENEWAL_V4_1_ITER", "2000"))
  warmup <- as.integer(Sys.getenv("RENEWAL_V4_1_WARMUP", "1000"))
  chains <- as.integer(Sys.getenv("RENEWAL_V4_1_CHAINS", "4"))

  if (!dir.exists(progress_dir)) {
    message("[progress] no progress directory found at ", progress_dir, " -- fit may not have started yet, or predates sample_file tracking.")
    return(invisible(NULL))
  }

  message(sprintf("[progress] run_tag=%s | target: %d warmup + %d sampling per chain", run_tag, warmup, iter - warmup))
  for (c in seq_len(chains)) {
    csv_path <- file.path(progress_dir, paste0("chain_", c, ".csv"))
    p <- count_csv_progress(csv_path, warmup, iter)
    if (!p$exists) {
      message(sprintf("  chain %d: not started (no file yet)", c))
    } else {
      pct <- round(100 * p$n_written / iter, 1)
      message(sprintf("  chain %d: %s | warmup %d/%d | sampling %d/%d | overall %s%%",
                       c, p$phase, p$warmup_done, warmup, p$sampling_done, iter - warmup, pct))
    }
  }
}

if (sys.nframe() == 0L) run_check_sampling_progress()
