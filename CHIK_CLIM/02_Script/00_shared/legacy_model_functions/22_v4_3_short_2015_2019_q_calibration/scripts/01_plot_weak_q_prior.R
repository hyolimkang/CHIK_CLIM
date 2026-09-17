# Pre-registers the weak calibration prior on q (Section 2), BEFORE any
# fitting: eta_q ~ normal(logit(0.05), 1.25), q = inv_logit(eta_q). Report
# only -- this prior is not tuned after seeing posterior results.

required_packages <- c("here", "tibble", "readr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing package(s): ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({ library(tibble); library(readr); library(ggplot2) })

prior_root <- function() {
  root <- here::here()
  if (file.exists(file.path(root, "CHIK_CLIM.Rproj"))) return(root)
  candidate <- file.path(root, "CHIK_CLIM")
  if (file.exists(file.path(candidate, "CHIK_CLIM.Rproj"))) return(candidate)
  stop("Could not identify the inner CHIK_CLIM project root.")
}

run_plot_weak_q_prior <- function() {
  root <- prior_root()
  out_dir <- file.path(root, "03_Output", "model_fits", "ceara", "v4_3_short_2015_2019_q_calibration")

  n_mc <- 2e6
  set.seed(20260912L)
  eta_q_draws <- rnorm(n_mc, mean = qlogis(0.05), sd = 1.25)
  q_draws <- plogis(eta_q_draws)

  summary_tbl <- tibble(
    statistic = c("median", "q25", "q75", "q025", "q975", "P(q<0.03)", "P(q<0.05)", "P(q>0.10)", "P(q>0.20)"),
    value = c(
      median(q_draws), quantile(q_draws, .25, names = FALSE), quantile(q_draws, .75, names = FALSE),
      quantile(q_draws, .025, names = FALSE), quantile(q_draws, .975, names = FALSE),
      mean(q_draws < 0.03), mean(q_draws < 0.05), mean(q_draws > 0.10), mean(q_draws > 0.20)
    )
  )
  write_csv(summary_tbl, file.path(out_dir, "weak_q_prior_summary.csv"))
  message("[weak q prior] eta_q ~ normal(logit(0.05)=", round(qlogis(0.05), 4), ", 1.25); q = inv_logit(eta_q)")
  print(as.data.frame(summary_tbl))

  p <- ggplot(tibble(q = q_draws[q_draws < 0.6]), aes(q)) +
    geom_histogram(bins = 100, fill = "#0072B2", alpha = .7) +
    geom_vline(xintercept = median(q_draws), linetype = 2, colour = "#D55E00") +
    labs(title = "v4.3 weak calibration prior on q (pre-registered, not tuned post-hoc)",
         subtitle = "eta_q ~ normal(logit(0.05), 1.25), q = inv_logit(eta_q)", x = "q", y = "Monte Carlo density (count)") +
    theme_classic(base_size = 9)
  figure_dir <- file.path(root, "03_Output", "figures", "renewal_v4_3_short_2015_2019_q_calibration")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figure_dir, "weak_q_prior.png"), p, width = 180, height = 120, units = "mm", dpi = 300)
  message("[weak q prior] figure saved: ", file.path(figure_dir, "weak_q_prior.png"))
  invisible(summary_tbl)
}

if (sys.nframe() == 0L) run_plot_weak_q_prior()
