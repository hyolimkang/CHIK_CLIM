# Live, read-only monitor for a running v4.1 modular-q node fit.
#
# Run in a SEPARATE terminal/session while 02_run_q_node.R is sampling:
#   Q_NODE_ID=q50 Rscript 03_monitor_v4_1_fit.R
#
# Reads the growing per-chain CmdStan-style CSVs written by sample_file
# (outputs/<node_id>/chains/chain_<id>.csv) every 60 seconds, WITHOUT
# modifying them, and safely ignores an incomplete final line (the row
# currently being written by Stan's C++ writer). Does not wait for the fit
# to finish -- exits on its own once it observes the node's fit .rds has
# been saved, or can be interrupted at any time with no side effects since
# it never writes to the chain CSVs themselves.
#
# Only the first 7 diagnostic columns (lp__, accept_stat__, stepsize__,
# treedepth__, n_leapfrog__, divergent__, energy__) are parsed per row --
# never the full row, which for this model has ~7000 columns (log_R0, R0_t,
# R_eff_t, force_of_infection, X, S, U, S_prop, immune_prop,
# expected_reported_cases, C_pred, log_lik_cases are each length N=573) --
# parsing only the diagnostic prefix keeps this monitor cheap even as the
# chain files grow to hundreds of MB.

monitor_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

n_diagnostic_cols <- 7L
diagnostic_col_names <- c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")

# Extracts only the first `n` comma-separated fields of a line without
# splitting the (potentially huge) remainder -- much faster than
# strsplit() on a ~7000-column row.
first_n_fields <- function(line, n) {
  commas <- gregexpr(",", line, fixed = TRUE)[[1]]
  if (length(commas) < n) return(NULL) # malformed/incomplete row
  cutoff <- commas[n]
  as.numeric(strsplit(substr(line, 1, cutoff - 1), ",", fixed = TRUE)[[1]])
}

# Reads one chain's CSV as-is right now. Returns NULL if the file does not
# exist yet or has no data rows. The last row is dropped if it does not
# parse as `n_diagnostic_cols` numbers (Stan's writer was still mid-write).
read_chain_snapshot <- function(csv_path, warmup, iter) {
  info <- if (file.exists(csv_path)) file.info(csv_path) else NULL
  if (is.null(info)) return(list(exists = FALSE))
  lines <- tryCatch(readLines(csv_path, warn = FALSE), error = function(e) character(0))
  data_lines <- lines[!startsWith(lines, "#")]
  if (length(data_lines) < 2) {
    return(list(exists = TRUE, n_rows = 0L, mtime = info$mtime, size_bytes = info$size, warmup_done = FALSE, diag = NULL))
  }
  rows <- data_lines[-1] # drop header
  parsed <- lapply(rows, first_n_fields, n = n_diagnostic_cols)
  ok <- !vapply(parsed, is.null, logical(1))
  parsed <- parsed[ok] # silently drops one incomplete trailing row, if any
  if (!length(parsed)) {
    return(list(exists = TRUE, n_rows = 0L, mtime = info$mtime, size_bytes = info$size, warmup_done = FALSE, diag = NULL))
  }
  m <- do.call(rbind, parsed)
  colnames(m) <- diagnostic_col_names
  adaptation_idx <- grep("Adaptation terminated", lines)
  warmup_done <- length(adaptation_idx) > 0
  n_warmup_rows <- if (warmup_done) {
    lines_before <- lines[seq_len(adaptation_idx[1] - 1L)]
    max(0L, length(lines_before[!startsWith(lines_before, "#")]) - 1L)
  } else {
    min(nrow(m), warmup)
  }
  list(
    exists = TRUE, n_rows = nrow(m), mtime = info$mtime, size_bytes = info$size,
    warmup_done = warmup_done, n_warmup_rows = n_warmup_rows,
    n_sampling_rows = max(0L, nrow(m) - n_warmup_rows),
    diag = m, diag_post_warmup = if (warmup_done) m[seq(n_warmup_rows + 1L, nrow(m)), , drop = FALSE] else m[0, , drop = FALSE]
  )
}

fmt_pct <- function(x, total) if (total > 0) sprintf("%.1f%%", 100 * min(x, total) / total) else "NA"

run_monitor_v4_1_fit <- function(node_id = Sys.getenv("Q_NODE_ID", "q50"),
                                  poll_seconds = 60, max_polls = Inf) {
  root <- monitor_root()
  base_dir <- file.path(root, "02_Script", "90_development_archive", "02_historical_renewal_versions", "19_v4_1_modular_q")
  out_dir <- file.path(base_dir, "outputs", node_id)
  chains_dir <- file.path(out_dir, "chains")
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_1_modular_q_fit_", node_id, ".rds"))
  status_path <- file.path(out_dir, "live_sampler_status.tsv")

  benchmark_path <- file.path(root, "03_Output", "tables", "renewal_v4_0_minimal_no_vaccine_td14", "hmc_diagnostics.csv")
  benchmark_seconds <- if (file.exists(benchmark_path)) {
    b <- read.csv(benchmark_path)
    b$value[b$metric == "elapsed_seconds"]
  } else NA_real_

  warmup <- 1000L; iter <- 2000L; chains <- 4L
  message(sprintf("[monitor:%s] watching %s (poll every %ds; Ctrl+C to stop, read-only, safe to interrupt anytime)", node_id, chains_dir, poll_seconds))

  poll <- 0L
  repeat {
    poll <- poll + 1L
    started_at <- Sys.time()
    snapshots <- lapply(seq_len(chains), function(c) read_chain_snapshot(file.path(chains_dir, sprintf("chain_%d.csv", c)), warmup, iter))

    rows_out <- list()
    any_gate_violation <- FALSE
    for (c in seq_len(chains)) {
      s <- snapshots[[c]]
      if (!isTRUE(s$exists) || is.null(s$diag)) {
        message(sprintf("[%s] chain %d: not started yet", format(started_at, "%H:%M:%S"), c))
        next
      }
      d <- s$diag
      dpw <- s$diag_post_warmup
      post_warmup_divergent <- if (nrow(dpw)) sum(dpw[, "divergent__"]) else 0L
      post_warmup_td_max <- if (nrow(dpw)) sum(dpw[, "treedepth__"] == 14) else 0L
      if (post_warmup_divergent > 0 || post_warmup_td_max > 0) any_gate_violation <- TRUE

      row <- list(
        timestamp = format(started_at, "%Y-%m-%d %H:%M:%S"), node_id = node_id, chain = c,
        n_rows = s$n_rows, warmup_progress = fmt_pct(s$n_warmup_rows, warmup),
        sampling_progress = fmt_pct(s$n_sampling_rows, iter - warmup),
        mtime = format(s$mtime, "%Y-%m-%d %H:%M:%S"), size_bytes = s$size_bytes,
        divergent_total = sum(d[, "divergent__"]), divergent_post_warmup = post_warmup_divergent,
        treedepth_median = median(d[, "treedepth__"]), treedepth_p95 = quantile(d[, "treedepth__"], .95, names = FALSE),
        treedepth_eq14_total = sum(d[, "treedepth__"] == 14), treedepth_eq14_post_warmup = post_warmup_td_max,
        n_leapfrog_median = median(d[, "n_leapfrog__"]), n_leapfrog_p95 = quantile(d[, "n_leapfrog__"], .95, names = FALSE),
        n_leapfrog_max = max(d[, "n_leapfrog__"]), stepsize_last = tail(d[, "stepsize__"], 1),
        accept_stat_median = median(d[, "accept_stat__"])
      )
      rows_out[[length(rows_out) + 1L]] <- row
      message(sprintf(
        "[%s] chain %d: rows=%d | warmup=%s sampling=%s | TD[med=%.1f p95=%.1f eq14=%d(post-wu:%d)] | leapfrog[med=%.0f p95=%.0f max=%d] | div[total=%d post-wu:%d] | stepsize=%.5f | accept_med=%.2f",
        row$timestamp, c, row$n_rows, row$warmup_progress, row$sampling_progress,
        row$treedepth_median, row$treedepth_p95, row$treedepth_eq14_total, row$treedepth_eq14_post_warmup,
        row$n_leapfrog_median, row$n_leapfrog_p95, row$n_leapfrog_max,
        row$divergent_total, row$divergent_post_warmup, row$stepsize_last, row$accept_stat_median
      ))
      if (post_warmup_divergent > 0) message("  *** WARNING: chain ", c, " has ", post_warmup_divergent, " POST-WARMUP divergent transition(s) -- strict gate already violated. ***")
      if (post_warmup_td_max > 0) message("  *** WARNING: chain ", c, " has ", post_warmup_td_max, " POST-WARMUP treedepth==14 hit(s) -- strict gate already violated. ***")
      if (row$n_leapfrog_p95 > 8000) message("  *** WARNING: chain ", c, " n_leapfrog__ 95th percentile is very large (", round(row$n_leapfrog_p95), ") -- geometry may be difficult. ***")
    }

    if (length(rows_out)) {
      total_rows <- sum(vapply(rows_out, function(r) r$n_rows, integer(1)))
      elapsed_so_far <- as.numeric(difftime(Sys.time(), started_at, units = "secs")) # negligible; real elapsed tracked by caller
      message(sprintf("[%s] TOTAL rows across chains: %d / %d target", format(started_at, "%H:%M:%S"), total_rows, chains * iter))
      if (!is.na(benchmark_seconds)) {
        implied_seconds_per_row <- NA # cannot know true elapsed without the fit script's own timer; report benchmark for reference only
        message(sprintf("[%s] v4.0/td14 benchmark for reference: %.0f seconds total (4 chains, 2000 iter each)", format(started_at, "%H:%M:%S"), benchmark_seconds))
      }
    }

    if (any_gate_violation) {
      message(">>> STRICT GATE ALREADY FAILED -- SAFE TO STOP THIS NODE <<<")
    }

    if (length(rows_out)) {
      snapshot_df <- do.call(rbind.data.frame, rows_out)
      write.table(snapshot_df, status_path, sep = "\t", row.names = FALSE,
                  col.names = !file.exists(status_path), append = file.exists(status_path))
    }

    if (file.exists(fit_path)) {
      message(sprintf("[monitor:%s] fit .rds now exists (%s) -- node finished, stopping monitor.", node_id, fit_path))
      break
    }
    if (poll >= max_polls) {
      message(sprintf("[monitor:%s] reached max_polls=%d -- stopping (fit may still be running; rerun to keep watching).", node_id, max_polls))
      break
    }
    Sys.sleep(poll_seconds)
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) run_monitor_v4_1_fit()
