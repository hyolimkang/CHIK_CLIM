# Live, read-only monitor for a running v4.1 q-only stage fit. Run in a
# SEPARATE terminal while 02_fit_v4_1_q_only.R is sampling:
#   STAGE=B Rscript 03_monitor_v4_1_q_only.R
#
# Reads the growing per-chain CmdStan-style CSVs every 60 seconds WITHOUT
# modifying them, safely ignoring an incomplete final line. Only the first 7
# diagnostic columns are parsed per row (never the ~7000-column full row) to
# stay cheap as files grow. Does not wait for the fit to finish.

monitor_root <- function() {
  root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

n_diagnostic_cols <- 7L
diagnostic_col_names <- c("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__")

first_n_fields <- function(line, n) {
  commas <- gregexpr(",", line, fixed = TRUE)[[1]]
  if (length(commas) < n) return(NULL)
  as.numeric(strsplit(substr(line, 1, commas[n] - 1), ",", fixed = TRUE)[[1]])
}

read_chain_snapshot <- function(csv_path, warmup) {
  info <- if (file.exists(csv_path)) file.info(csv_path) else NULL
  if (is.null(info)) return(list(exists = FALSE))
  lines <- tryCatch(readLines(csv_path, warn = FALSE), error = function(e) character(0))
  data_lines <- lines[!startsWith(lines, "#")]
  if (length(data_lines) < 2) return(list(exists = TRUE, n_rows = 0L, mtime = info$mtime, size_bytes = info$size, diag = NULL))
  rows <- data_lines[-1]
  parsed <- lapply(rows, first_n_fields, n = n_diagnostic_cols)
  ok <- !vapply(parsed, is.null, logical(1))
  parsed <- parsed[ok]
  if (!length(parsed)) return(list(exists = TRUE, n_rows = 0L, mtime = info$mtime, size_bytes = info$size, diag = NULL))
  m <- do.call(rbind, parsed)
  colnames(m) <- diagnostic_col_names
  adaptation_idx <- grep("Adaptation terminated", lines)
  warmup_done <- length(adaptation_idx) > 0
  n_warmup_rows <- if (warmup_done) {
    lines_before <- lines[seq_len(adaptation_idx[1] - 1L)]
    max(0L, length(lines_before[!startsWith(lines_before, "#")]) - 1L)
  } else min(nrow(m), warmup)
  list(
    exists = TRUE, n_rows = nrow(m), mtime = info$mtime, size_bytes = info$size,
    warmup_done = warmup_done, n_warmup_rows = n_warmup_rows, n_sampling_rows = max(0L, nrow(m) - n_warmup_rows),
    diag = m, diag_post_warmup = if (warmup_done) m[seq(n_warmup_rows + 1L, nrow(m)), , drop = FALSE] else m[0, , drop = FALSE]
  )
}

fmt_pct <- function(x, total) if (total > 0) sprintf("%.1f%%", 100 * min(x, total) / total) else "NA"

run_monitor_v4_1_q_only <- function(stage = Sys.getenv("STAGE", "B"), poll_seconds = 60, max_polls = Inf) {
  root <- monitor_root()
  base_dir <- file.path(root, "02_Script", "40_renewal_model", "20_v4_1_q_only")
  out_dir <- file.path(base_dir, "outputs", paste0("stage", stage))
  chains_dir <- file.path(out_dir, "chains")
  fit_path <- file.path(out_dir, paste0("renewal_ceara_v4_1_q_only_fit_stage", stage, ".rds"))
  status_path <- file.path(out_dir, "live_sampler_status.tsv")

  stage_warmup <- switch(stage, A = 5L, B = 300L, C = 500L, D = 1000L)
  stage_sampling <- switch(stage, A = 5L, B = 150L, C = 300L, D = 1000L)
  chains <- 4L

  benchmark_path <- file.path(root, "03_Output", "tables", "renewal_v4_0_minimal_no_vaccine_td14", "hmc_diagnostics.csv")
  benchmark_seconds <- if (file.exists(benchmark_path)) { b <- read.csv(benchmark_path); b$value[b$metric == "elapsed_seconds"] } else NA_real_

  message(sprintf("[monitor:stage%s] watching %s (poll every %ds; read-only, safe to interrupt anytime)", stage, chains_dir, poll_seconds))

  poll <- 0L
  chain_q_history <- vector("list", chains)
  repeat {
    poll <- poll + 1L
    started_at <- Sys.time()
    snapshots <- lapply(seq_len(chains), function(c) read_chain_snapshot(file.path(chains_dir, sprintf("chain_%d.csv", c)), stage_warmup))

    rows_out <- list(); any_gate_violation <- FALSE
    for (c in seq_len(chains)) {
      s <- snapshots[[c]]
      if (!isTRUE(s$exists) || is.null(s$diag)) { message(sprintf("[%s] chain %d: not started yet", format(started_at, "%H:%M:%S"), c)); next }
      d <- s$diag; dpw <- s$diag_post_warmup
      post_warmup_divergent <- if (nrow(dpw)) sum(dpw[, "divergent__"]) else 0L
      post_warmup_td_max <- if (nrow(dpw)) sum(dpw[, "treedepth__"] == 14) else 0L
      if (post_warmup_divergent > 0 || post_warmup_td_max > 0) any_gate_violation <- TRUE

      row <- list(
        timestamp = format(started_at, "%Y-%m-%d %H:%M:%S"), stage = stage, chain = c,
        n_rows = s$n_rows, warmup_progress = fmt_pct(s$n_warmup_rows, stage_warmup), sampling_progress = fmt_pct(s$n_sampling_rows, stage_sampling),
        mtime = format(s$mtime, "%Y-%m-%d %H:%M:%S"), size_bytes = s$size_bytes,
        divergent_total = sum(d[, "divergent__"]), divergent_post_warmup = post_warmup_divergent,
        treedepth_median = median(d[, "treedepth__"]), treedepth_p90 = quantile(d[, "treedepth__"], .90, names = FALSE),
        treedepth_p95 = quantile(d[, "treedepth__"], .95, names = FALSE), treedepth_p99 = quantile(d[, "treedepth__"], .99, names = FALSE),
        treedepth_eq14_total = sum(d[, "treedepth__"] == 14), treedepth_eq14_post_warmup = post_warmup_td_max,
        n_leapfrog_median = median(d[, "n_leapfrog__"]), n_leapfrog_p95 = quantile(d[, "n_leapfrog__"], .95, names = FALSE),
        n_leapfrog_max = max(d[, "n_leapfrog__"]), stepsize_last = tail(d[, "stepsize__"], 1),
        accept_stat_median = median(d[, "accept_stat__"])
      )
      rows_out[[length(rows_out) + 1L]] <- row
      message(sprintf(
        "[%s] chain %d: rows=%d | warmup=%s sampling=%s | TD[med=%.1f p90=%.1f p95=%.1f p99=%.1f eq14=%d(post-wu:%d)] | leapfrog[med=%.0f p95=%.0f max=%d] | div[total=%d post-wu:%d] | stepsize=%.5f | accept_med=%.2f",
        row$timestamp, c, row$n_rows, row$warmup_progress, row$sampling_progress,
        row$treedepth_median, row$treedepth_p90, row$treedepth_p95, row$treedepth_p99, row$treedepth_eq14_total, row$treedepth_eq14_post_warmup,
        row$n_leapfrog_median, row$n_leapfrog_p95, row$n_leapfrog_max, row$divergent_total, row$divergent_post_warmup, row$stepsize_last, row$accept_stat_median
      ))
      if (post_warmup_divergent > 0) message("  *** WARNING: chain ", c, " has ", post_warmup_divergent, " POST-WARMUP divergent transition(s). ***")
      if (post_warmup_td_max > 0) message("  *** WARNING: chain ", c, " has ", post_warmup_td_max, " POST-WARMUP treedepth==14 hit(s). ***")
      if (row$n_leapfrog_p95 > 8000) message("  *** WARNING: chain ", c, " n_leapfrog__ 95th pct is very large (", round(row$n_leapfrog_p95), "). ***")

      # Track alpha_R proxy via lp__ drift is unreliable; instead track q if
      # present among the first columns is not guaranteed -- use per-chain
      # divergence/treedepth trend as the practical multi-chain-separation
      # proxy available cheaply here (full parameter separation is checked
      # post-hoc in 04_check_hmc_gate.R / 05_compare_vs_v4_0.R).
      chain_q_history[[c]] <- c(chain_q_history[[c]], row$treedepth_median)
    }

    if (length(rows_out)) {
      total_rows <- sum(vapply(rows_out, function(r) r$n_rows, integer(1)))
      message(sprintf("[%s] TOTAL rows across chains: %d / %d target", format(started_at, "%H:%M:%S"), total_rows, chains * (stage_warmup + stage_sampling)))
      if (!is.na(benchmark_seconds)) message(sprintf("[%s] v4.0/td14 benchmark for reference: %.0f seconds (full 1000+1000 run)", format(started_at, "%H:%M:%S"), benchmark_seconds))
      td_meds <- vapply(rows_out, function(r) r$treedepth_median, numeric(1))
      if (length(td_meds) == 4 && (max(td_meds) - min(td_meds)) > 3) {
        message("  *** WARNING: chains show markedly different median treedepth (range ", round(max(td_meds) - min(td_meds), 1), ") -- possible chain separation into different regions. ***")
      }
    }
    if (any_gate_violation) message(">>> STRICT GATE ALREADY FAILED (this stage) <<<")

    if (length(rows_out)) {
      snapshot_df <- do.call(rbind.data.frame, rows_out)
      write.table(snapshot_df, status_path, sep = "\t", row.names = FALSE, col.names = !file.exists(status_path), append = file.exists(status_path))
    }
    if (file.exists(fit_path)) { message(sprintf("[monitor:stage%s] fit .rds now exists -- stage finished, stopping monitor.", stage)); break }
    if (poll >= max_polls) { message(sprintf("[monitor:stage%s] reached max_polls=%d -- stopping.", stage, max_polls)); break }
    Sys.sleep(poll_seconds)
  }
  invisible(NULL)
}

if (sys.nframe() == 0L) run_monitor_v4_1_q_only()
